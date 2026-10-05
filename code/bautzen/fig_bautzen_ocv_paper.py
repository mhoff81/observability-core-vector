#!/usr/bin/env python3
"""OCV-Paper (Bautzen bridge site) — figures A–E, JSON result and English report.

Question: does the observability core vector
``OCV = [gamma^2, A, D, F, S, P]`` (deck-band coherence gamma^2 plus the mask
channels area A, density D, fragmentation F, shift S and persistence P of the
**deck-edge mask**) respond to the two documented interventions at the Bautzen
OpenLabs bridge (states RS -> DS1 -> DS2) — and does the *control* area 500 m
east of the bridge respond the same way (i.e. is the response seasonal rather
than structural)?

This is the plotting counterpart of ``code/lumo/fig_lumo_ocv_paper.py``: the same
structure (recompute everything from the committed channel table, pin it against
the committed reference JSONs, one figure PNG set, one result JSON, one Markdown
report, no timestamps, exit code != 0 on any pin failure). Only the mask layer
and the state model differ, and both differences are documented here and in the
report:

  * the mask layer is the **deck-edge mask** of ``bautzen_ocv_masks.py`` (the
    global echo rule of the LUMO site collapses to n ~ 0 at this clutter-
    dominated rect chip, see that module), so the OCV "observability" of this
    site is the visibility/contrast of the deck line;
  * the three states are **time intervals of one deployment** whose boundaries
    fall in the middle of a day (``core.state_of`` assigns them by overpass
    time). State and season are therefore confounded by construction — which is
    exactly why the four series (bridge/control x ASC/DESC) are carried through
    every figure as the counter-check.

Figures:
  A  OCV components per state, bridge and control of both orbit directions
  B  deck-edge visibility: mask share, band contrast, detectability gate,
     bridge-vs-control discrimination per state
  C  intervention response: median steps at the two cuts against the
     +/-2..12 week placebo cuts, for band contrast and for gamma^2
  D  bridge minus control: paired per-date differences per state and the
     step/placebo reading of that difference
  E  seasonality and environment controls: monthly medians, and gamma^2 against
     air temperature / relative humidity at the overpass instant

Usage:
  python3 code/bautzen/fig_bautzen_ocv_paper.py            # recompute everything
  python3 code/bautzen/fig_bautzen_ocv_paper.py --quick    # placebo sweeps copied
"""

from __future__ import annotations

import argparse
import copy
import datetime as dt
import json
import math
import os
import platform
import sys
import textwrap
import time
from collections import Counter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

try:                                    # scipy is optional in bautzen_ocv_stats
    import scipy
except ImportError:                     # pragma: no cover
    scipy = None

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import bautzen_ocv_core as core    # noqa: E402
import bautzen_ocv_masks as masks  # noqa: E402
import bautzen_ocv_stats as st     # noqa: E402

FIGDIR = core.FIGDIR
OUT_JSON = os.path.join(core.DATA, "fig_bautzen_ocv_paper.json")
OUT_MD = os.path.join(FIGDIR, "fig_bautzen_ocv_paper.md")

# The committed blocks whose placebo sweeps are the heavy part of a run;
# ``--quick`` copies them from ``bautzen_deck_edge_state.json`` (and marks them
# as copied) instead of recomputing them.
QUICK_COPIED = ("placebo", "gamma2_placebo")

# ---------------------------------------------------------------------------
# Style (identical to the LUMO figures: DejaVu Sans + STIX, 183 mm column width)
# ---------------------------------------------------------------------------
plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "mathtext.fontset": "stix",
    "font.size": 8.5,
    "axes.titlesize": 9.0,
    "axes.labelsize": 8.5,
    "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5,
    "legend.fontsize": 7.5,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.alpha": 0.25,
    "grid.linewidth": 0.4,
    "figure.dpi": 110,
    "savefig.dpi": 600,
    "savefig.bbox": "tight",
    "pdf.fonttype": 42,
})
MM = 1.0 / 25.4
COL_WIDTH = 183 * MM

# State colours/markers (three states, in the severity order of STATE_ORDER).
STATE_COLORS = {"RS": "#4C72B0", "DS1": "#DD8452", "DS2": "#C44E52"}
STATE_MARKERS = {"RS": "o", "DS1": "s", "DS2": "^"}
# The four rect series (committed German keys).
SERIES_COLORS = {
    "ASC Brücke": "#4C72B0",
    "DESC Brücke": "#55A868",
    "ASC Kontrolle 500 m östlich": "#DD8452",
    "DESC Kontrolle 500 m östlich": "#8172B3",
}
SERIES_MARKERS = {
    "ASC Brücke": "o",
    "DESC Brücke": "s",
    "ASC Kontrolle 500 m östlich": "^",
    "DESC Kontrolle 500 m östlich": "D",
}
CHANCE = 1.0 / 3.0

# The OCV vector of this site = the site's coherence channel + the five mask
# channels (``core.DIMS``); S is measured against the deck-local peak, see
# bautzen_ocv_masks.py.
CHANNELS = ["gamma2"] + core.DIMS
CHANNEL_AXIS = {"gamma2": "gamma^2", "A": "A  [px]", "D": "D  [1]",
                "F": "F  [components]", "S": "S  [px]", "P": "P  [1]"}
CHANNEL_TITLE = {
    "gamma2": "deck-band coherence gamma^2",
    "A": "mask area (deck-edge pixels)", "D": "mask density",
    "F": "mask fragmentation", "S": "centroid-to-peak shift",
    "P": "mask persistence",
}

# The gamma^2 variants of the committed gamma2_step_tests / gamma2_placebo
# blocks: (block key, column, row gate, gate label used by the reference).
G2_VARIANTS = (
    (core.GAMMA2_COLUMN, core.GAMMA2_COLUMN, None, None),
    (f"{core.GAMMA2_COLUMN}_n_ge_5", core.GAMMA2_COLUMN,
     lambda r: (r["n_band"] or 0) >= masks.LOW_N_GATE,
     f"n_band >= {masks.LOW_N_GATE}"),
    ("gamma2_house_floor_corrected", "gamma2_house_floor_corrected", None, None),
)
# Short names of the gamma^2 variants for the tables of the Markdown report.
G2_SHORT = {core.GAMMA2_COLUMN: "band",
            f"{core.GAMMA2_COLUMN}_n_ge_5": "band, n_band >= 5",
            "gamma2_house_floor_corrected": "house"}


def series_label(lab):
    """'ASC Brücke' -> 'ASC bridge' (the report is English, keys stay German)."""
    return core.SERIES_EN.get(lab, lab)


def save_fig(fig, name, out_dir=None):
    """Write one 600 dpi PNG (the device set keeps the binary out of git)."""
    out_dir = out_dir or FIGDIR
    png = os.path.join(out_dir, name + ".png")
    fig.savefig(png)
    plt.close(fig)
    print(f"  {os.path.relpath(png)}")
    return png


def panel(ax, letter, text, fs=7.4):
    """Panel caption ``(a) text``, wrapped to the width of its own axes.

    The captions are left-aligned and would otherwise grow into the caption of
    the neighbouring panel; the character budget follows from the axes width,
    so the three-column panels of figure A wrap tighter than the two-column
    panels of figures B-E.  ``textwrap`` breaks on whole words only.
    """
    fig = ax.get_figure()
    width_in = ax.get_position().width * fig.get_size_inches()[0] - 0.06
    n = max(16, int(width_in / (fs * 0.0093)))
    title = f"({letter}) {text}"
    lines = textwrap.wrap(title, n) or [title]
    ax.set_title("\n".join(lines), fontsize=fs, loc="left")


def jitter(n, width, seed):
    rng = np.random.default_rng(seed)
    return rng.uniform(-width, width, n)


def md_table(header, rows):
    head = "| " + " | ".join(header) + " |"
    rule = "|" + "|".join(["---"] * len(header)) + "|"
    body = ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]
    return "\n".join([head, rule] + body)


def fmt(x, nd=3):
    if x is None:
        return "n/a"
    if isinstance(x, bool):
        return "yes" if x else "no"
    if isinstance(x, float) and math.isnan(x):
        return "n/a"
    if isinstance(x, (int, np.integer)):
        return str(int(x))
    return f"{float(x):.{nd}f}"


def signed(x, nd=3):
    """Signed value for the prose: '+0.139', '-0.053', '0.000'."""
    if x is None:
        return "n/a"
    return f"{x:+.{nd}f}" if x else f"{x:.{nd}f}"




def mask_definition(key):
    """One committed mask-layer constant (``bautzen_ocv_channels_meta.json``)."""
    return core.load_csv_meta()["mask_definition"][key]


FLOAT_TOL = 1e-6
TOL_KEYS = {"perm_p", "p", "p_perm", "welch_p", "mannwhitney_p", "rho"}


class Pins:
    """Collects pin comparisons (recomputed vs. committed) and reports deviations.

    Comparison rule: exact equality for everything deterministic (delta, SDS,
    medians, counters, masks, gamma2); 1e-6 for permutation p-values and scipy
    results (library version / normalisation). Deviations above the tolerance
    and missing keys => exit code != 0.
    """

    def __init__(self):
        self.checks = []
        self.fails = []
        self.extra_keys = []
        self.label_translations = 0

    def add(self, path, got, want, tol):
        # Committed German doc strings (see core.LABEL_EN) are compared in their
        # English translation; the count is reported in the summary so that no
        # translation can silently hide a genuine deviation.
        if isinstance(got, str) and got in core.LABEL_EN:
            got = core.LABEL_EN[got]
            self.label_translations += 1
        if isinstance(want, str) and want in core.LABEL_EN:
            want = core.LABEL_EN[want]
            self.label_translations += 1
        # The mask unions are Python tuples here but lists in the committed JSON;
        # normalise before the leaf comparison (documented deviation from the
        # LUMO Pins, which never meets a tuple at a leaf).
        if isinstance(got, (list, tuple)) and isinstance(want, (list, tuple)):
            got, want = list(got), list(want)
        exact = bool(got == want)
        dev = 0.0
        if not exact:
            numeric = (isinstance(got, (int, float)) and isinstance(want, (int, float))
                       and not isinstance(got, bool) and not isinstance(want, bool))
            dev = abs(float(got) - float(want)) if numeric else None
        ok = exact or (dev is not None and dev <= tol)
        self.checks.append({"path": path, "kind": "exact" if exact else "tolerance",
                            "abs_deviation": dev, "tolerance": tol, "ok": ok})
        if not ok:
            self.fails.append({"path": path, "recomputed": repr(got),
                               "committed": repr(want), "abs_deviation": dev})

    def cmp(self, prefix, got, want, key=None):
        """Recursive comparison: `want` (committed) dictates structure and expectation."""
        tol = FLOAT_TOL if key in TOL_KEYS else 0.0
        if isinstance(want, dict):
            if not isinstance(got, dict):
                self.add(prefix, f"<{type(got).__name__}>", "<dict>", 0.0)
                return
            for k in want:
                if k not in got:
                    self.add(f"{prefix}.{k}", "<missing>", "<present>", 0.0)
                    continue
                self.cmp(f"{prefix}.{k}", got[k], want[k], k)
            for k in got:
                if k not in want:
                    self.extra_keys.append(f"{prefix}.{k}")
            return
        if isinstance(want, list):
            glen = len(got) if hasattr(got, "__len__") else -1
            if not isinstance(got, (list, tuple)) or glen != len(want):
                self.add(prefix, f"<list len {glen}>", f"<list len {len(want)}>", 0.0)
                return
            for i, (g, w) in enumerate(zip(got, want)):
                self.cmp(f"{prefix}[{i}]", g, w, key)
            return
        self.add(prefix, got, want, tol)

    def summary(self):
        devs = [c["abs_deviation"] for c in self.checks
                if c["abs_deviation"] not in (None, 0.0)]
        return {
            "n_checks": len(self.checks),
            "n_exact": sum(1 for c in self.checks if c["kind"] == "exact"),
            "n_within_tolerance": sum(1 for c in self.checks if c["kind"] == "tolerance"),
            "n_failed": len(self.fails),
            "max_abs_deviation": max(devs) if devs else 0.0,
            "label_translations": self.label_translations,
            "float_tolerance": FLOAT_TOL,
            "tolerance_keys": sorted(TOL_KEYS),
            "failures": self.fails,
            "recomputed_keys_not_in_reference": sorted(set(self.extra_keys)),
        }


# ---------------------------------------------------------------------------
# Recomputation layer — every pinned number, from the committed channel table
# ---------------------------------------------------------------------------
# The committed German doc strings of the state reference (``core.LABEL_EN``
# maps each of them to an English text; the pin block goes through that map, so
# a stale translation shows up as a failure instead of passing silently). The
# strings below are used verbatim; the metrics-level gamma^2 strings are not
# among them -- ``masks.state_metrics`` states those itself (in English, and
# ``core.LABEL_EN`` maps the committed German ones onto exactly that text).
DOC_DE = {
    "mask_source": "analyze_bautzen_deck_edge_mask.py (unverändert übernommen)",
    "formula": "|sum(z)|^2 / (N * sum(|z|^2)) ueber die Maskenpixel",
    # config.gamma2_detail.floor_corrected. Note that the metrics block carries a
    # *second* German floor-correction string (LABEL_EN maps both; the second
    # translation is masks.GAMMA2_FLOOR_CORRECTION, which state_metrics() emits).
    "floor_corrected": ("(N*gamma2 - 1)/(N - 1), inkohaerentes Nullmodell "
                        "(zirkular-Gauss, unabhaengige Looks); keine Klemmung, "
                        "negative Werte bleiben erhalten"),
    "channel_role": "komplementaerer Kanal, nicht unabhaengig",
    "mask_band": "Deckbandmaske = strukturell attribuiert (Deckkanten-Aussage)",
    "mask_house": ("Haus-/Peakmaske = signalbasiert (|z|^2 >= 0.3*Peak and "
                   ">= 5.0*Median, min 2 px), NICHT strukturell attribuiert -> "
                   "keine Deckkanten-Aussage"),
    "interpolation_generic": ("linear zwischen den Nachbarstunden des "
                              "Überflugzeitpunkts (keine Stundenrundung)"),
    "interpolation_geom": ("linear zwischen den Nachbarstunden der Überflugzeit "
                           "({hh:02d}:{mm:02d} lokal)"),
    "role": ("Gegenprobe zur Saisondeutung: Zustand und Jahreszeit sind "
             "konfundiert, T/RH sind der einfachste Saison-Ersatz"),
}
# The committed keys of one state-reference per-date row (exact key set, see
# bautzen_ocv_channels_csv.py; the mask JSON row is its 10/12-key subset).
STATE_ROW_KEYS = ("date", "state", "band_contrast", "peak_row", "peak_row_offset",
                  "n_masked", "has_mask", "centroid_row", "gamma2_band_raw",
                  "gamma2_band_floor_corrected", "n_band", "gamma2_house_raw",
                  "gamma2_house_floor_corrected", "n_house")

_METEO_META = {}


def meteo_meta():
    """Committed provenance of the meteo file (timezone / units)."""
    if not _METEO_META:
        with open(core.METEO_PATH) as fh:
            d = json.load(fh)
        _METEO_META.update({"file": os.path.basename(core.METEO_PATH),
                            "timezone": d.get("timezone"),
                            "units": d.get("hourly_units", {})})
    return dict(_METEO_META)


def interp_string(geom):
    """The committed German interpolation note of one orbit direction."""
    hh, mm = core.OVERPASS_LOCAL[geom]
    return DOC_DE["interpolation_geom"].format(hh=hh, mm=mm)


def config_state():
    """The committed ``config`` block of bautzen_deck_edge_state.json."""
    mm = meteo_meta()
    return {
        "cut_ds1": core.CUT_DS1.isoformat(),
        "cut_ds2": core.CUT_DS2.isoformat(),
        "overpass_local": {g: "%02d:%02d" % core.OVERPASS_LOCAL[g]
                           for g in ("ASC", "DESC")},
        "states": list(core.STATES),
        "series_order": list(core.SERIES_ORDER),
        "bridge_by_geom": dict(core.BRIDGE_BY_GEOM),
        "control_by_geom": dict(core.CONTROL_BY_GEOM),
        "mask_source": DOC_DE["mask_source"],
        "center_row": masks.CENTER_ROW,
        "center_col": masks.CENTER_COL,
        "col_half": masks.COL_HALF,
        "band_half": masks.BAND_HALF,
        "k_local": masks.K_LOCAL,
        "min_n_masked": masks.MIN_N_MASKED,
        "contrast_detect": masks.CONTRAST_DETECT,
        "n_perm": st.N_PERM,
        "n_boot": st.N_BOOT,
        "rng_seed": st.RNG_SEED,
        "boot_seed": st.BOOT_SEED,
        "placebo_weeks": list(st.PLACEBO_WEEKS),
        "gamma2_detail": {
            "formula": DOC_DE["formula"],
            "floor_corrected": DOC_DE["floor_corrected"],
            "channel_role": DOC_DE["channel_role"],
            "mask_band": DOC_DE["mask_band"],
            "mask_house": DOC_DE["mask_house"],
            "low_n_gate": masks.LOW_N_GATE,
        },
        "environment_source": {
            "file": mm["file"],
            "timezone": mm["timezone"],
            "variables": [{"name": n, "tag": t, "unit": u}
                          for n, t, u in core.ENV_VARS],
            "interpolation": DOC_DE["interpolation_generic"],
            "role": DOC_DE["role"],
        },
    }


def mask_row(row):
    """The committed mask-reference per-date row (its exact key set).

    Dates without a deck-edge pixel carry no centroid -- the committed rows are
    not equally wide on every date. The channel table holds more per-date
    columns than the committed mask block (``A``/``D``/``F``/``S``, band peak,
    bbox): they belong to the state reference and to the figures, so the mask
    block stays a byte-for-byte copy of the committed one.
    """
    out = {
        "anchor_row": masks.CENTER_ROW,
        "peak_row": int(row["peak_row"]),
        "peak_row_offset": int(row["peak_row_offset"]),
        "n_masked": int(row["n_masked"]),
        "band_contrast": float(row["band_contrast"]),
        "n_masked_house": int(row["n_masked_house"]),
        "n_components": int(row["n_components"]),
        "largest_component": int(row["largest_component"]),
        "has_mask": bool(row["has_mask"]),
    }
    if row["n_masked"] > 0:
        out["centroid_row"] = float(row["centroid_row"])
        out["centroid_col"] = float(row["centroid_col"])
    return out


def series_mask_block(rows):
    """The mask-reference block of one series, rebuilt from the channel table.

    The channel table stores every ingredient of ``date_mask`` verbatim (that is
    what ``bautzen_ocv_channels_csv.py --verify-only`` proves), so the series
    summary and the two persistence unions are rebuilt here without reading the
    rect files again.
    """
    per_date = [dict({"date": r["date"]}, **mask_row(r)) for r in rows]
    n = len(rows)
    contrast = np.asarray([r["band_contrast"] for r in rows], float)
    offsets = np.asarray([r["peak_row_offset"] for r in rows], float)
    nm = np.asarray([r["n_masked"] for r in rows], float)
    u50, u30 = masks.unions(masks.counts_of(rows), n)
    with_mask = int(sum(1 for r in rows if r["has_mask"]))
    summary = {
        "n_dates": n,
        "n_dates_with_mask": with_mask,
        "frac_dates_with_mask": with_mask / n if n else None,
        "band_contrast_median": float(np.median(contrast)) if n else None,
        "band_contrast_max": float(contrast.max()) if n else None,
        "frac_dates_contrast_gt_thr": (float((contrast >= masks.CONTRAST_DETECT).mean())
                                       if n else None),
        "peak_row_median": (float(np.median([r["peak_row"] for r in rows]))
                            if n else None),
        "peak_row_offset_median": (float(np.median(np.abs(offsets))) if n else None),
        "frac_dates_peak_within_1": (float((np.abs(offsets) <= 1).mean())
                                     if n else None),
        "n_masked_median": float(np.median(nm)) if n else None,
        "union50_n": len(u50),
        "union30_n": len(u30),
        "n_union_dates": n,
    }
    summary["deck_edge_detectable"] = bool(
        n > 0 and summary["band_contrast_median"] >= masks.CONTRAST_DETECT
        and summary["frac_dates_contrast_gt_thr"] >= 0.5
        and with_mask >= 0.5 * n)
    return {"file": rows[0]["rect_file"], "per_date": per_date, "summary": summary,
            "union50": u50, "union30": u30}


def jaccard_block(rows):
    """Jaccard overlap of the >= 50 % persistence unions of the three states.

    An empty union (both states without a single mask pixel) is undefined, not
    0 -- the committed reference stores null there.
    """
    sets = {}
    for s in core.STATES:
        sub = [r for r in rows if r["state"] == s]
        u50, _u30 = masks.unions(masks.counts_of(sub), len(sub))
        sets[s] = {(r, c) for r, c, _f in u50}
    out = {}
    for a, b in (("RS", "DS1"), ("DS1", "DS2"), ("RS", "DS2")):
        union = sets[a] | sets[b]
        out[f"{a}|{b}"] = (len(sets[a] & sets[b]) / len(union)) if union else None
    return out


def series_state_block(rows, lab):
    """The state-reference block of one series (committed key set)."""
    per_state = {s: [r for r in rows if r["state"] == s] for s in core.STATES}
    metrics = {}
    for s in list(core.STATES) + ["ALL"]:
        sub = rows if s == "ALL" else per_state[s]
        metrics[s] = masks.state_metrics(sub, masks.counts_of(sub))
    return {
        "file": rows[0]["rect_file"],
        "geom": core.geom_of(lab),
        "n_dates": len(rows),
        "n_by_state": {s: len(per_state[s]) for s in core.STATES},
        "dates_by_state": {s: [r["date"] for r in per_state[s]]
                           for s in core.STATES},
        "metrics": metrics,
        "jaccard": jaccard_block(rows),
        "per_date": [{k: r[k] for k in STATE_ROW_KEYS} for r in rows],
    }


def _keys_vals(rows, col, gate=None):
    """(overpass datetimes, column values, selected rows) of one series.

    ``gate`` filters rows by their own ingredients (e.g. ``n_band >= 5``); rows
    without a value in ``col`` are always dropped -- the ported estimators cannot
    read a null.
    """
    sel = [r for r in rows if r[col] is not None and (gate is None or gate(r))]
    geom = core.geom_of(rows[0]["series"])
    return ([core.acq_dt(r["date"], geom) for r in sel], [r[col] for r in sel], sel)


def step_block(rows, col, gate=None):
    """Median step at the two documented cuts (DS1/DS2) for one channel."""
    keys, vals, _sel = _keys_vals(rows, col, gate)
    return {"DS1": st.step_test(keys, vals, core.CUT_DS1),
            "DS2": st.step_test(keys, vals, core.CUT_DS2)}


def placebo_block(rows, col, gate=None):
    """The same two cuts against the +/-2..12 week placebo cuts."""
    keys, vals, _sel = _keys_vals(rows, col, gate)
    return {"DS1": st.placebo_sweep(keys, vals, core.CUT_DS1),
            "DS2": st.placebo_sweep(keys, vals, core.CUT_DS2)}


def season_block(rows, col, tag):
    """Per-month count + median of one channel (the season table of the site)."""
    monthly = {}
    for r in rows:
        if r[col] is None:
            continue
        monthly.setdefault(r["month"], []).append(r[col])
    return {m: {"n": len(v), f"{tag}_median": float(np.median(v))}
            for m, v in sorted(monthly.items())}


def gamma2_variant_blocks(rows, quick=False):
    """The three gamma^2 variants of gamma2_step_tests / gamma2_placebo."""
    steps, placebo = {}, {}
    for key, col, gate, label in G2_VARIANTS:
        keys, vals, sel = _keys_vals(rows, col, gate)
        steps[key] = {"min_n_filter": label, "n_values": len(sel),
                      "DS1": st.step_test(keys, vals, core.CUT_DS1),
                      "DS2": st.step_test(keys, vals, core.CUT_DS2)}
        placebo[key] = (None if quick else
                        {"DS1": st.placebo_sweep(keys, vals, core.CUT_DS1),
                         "DS2": st.placebo_sweep(keys, vals, core.CUT_DS2)})
    return steps, placebo


def paired_diff_series(by, col):
    """Per-date bridge-minus-control differences of one channel (both geoms)."""
    out = {}
    for geom in ("ASC", "DESC"):
        br = {r["date"]: r for r in by[core.BRIDGE_BY_GEOM[geom]]}
        ck = {r["date"]: r for r in by[core.CONTROL_BY_GEOM[geom]]}
        dates = sorted(set(br) & set(ck))
        if col != "band_contrast":
            dates = [d for d in dates
                     if br[d][col] is not None and ck[d][col] is not None]
        out[geom] = {"geom": geom, "channel": col, "dates": dates,
                     "diffs": [br[d][col] - ck[d][col] for d in dates],
                     "state": [br[d]["state"] for d in dates]}
    return out


def delta_vs_ctrl(by, col, with_feature=False, quick=False):
    """Bridge minus control, date-paired, in the committed block shape.

    ``col == 'band_contrast'`` keeps every date the two series have in common
    (band contrast is defined everywhere); a gamma^2 difference needs a mask on
    both sides, so those dates are dropped.
    """
    out = {}
    for geom, ser in paired_diff_series(by, col).items():
        dates, vals = ser["dates"], ser["diffs"]
        keys = [core.acq_dt(d, geom) for d in dates]
        states = {}
        for s in core.STATES:
            v = [val for val, stt in zip(vals, ser["state"]) if stt == s]
            states[s] = ({"n": len(v), "median": float(np.median(v)),
                          "iqr": [float(np.percentile(v, 25)),
                                  float(np.percentile(v, 75))]}
                         if v else None)
        blk = {"n_pairs": len(dates), "states": states,
               "step": {"DS1": st.step_test(keys, vals, core.CUT_DS1),
                        "DS2": st.step_test(keys, vals, core.CUT_DS2)},
               "placebo": (None if quick else
                           {"DS1": st.placebo_sweep(keys, vals, core.CUT_DS1),
                            "DS2": st.placebo_sweep(keys, vals, core.CUT_DS2)})}
        out[geom] = dict({"feature": col}, **blk) if with_feature else blk
    return out


def environment_block(rows, lab, env):
    """gamma^2 against the season proxy T/RH at the overpass instant (committed)."""
    geom = core.geom_of(lab)
    monthly = {}
    t_vals, rh_vals, g_vals = [], [], []
    n_env = 0
    for r in rows:
        e = st.env_at(env, core.acq_dt(r["date"], geom))
        if r[core.GAMMA2_COLUMN] is None:
            continue
        if e is not None:
            n_env += 1
        m = monthly.setdefault(r["month"], {"g": [], "T": [], "RH": []})
        m["g"].append(r[core.GAMMA2_COLUMN])
        if e is None:
            continue
        t_vals.append(e["T"])
        rh_vals.append(e["RH"])
        g_vals.append(r[core.GAMMA2_COLUMN])
        m["T"].append(e["T"])
        m["RH"].append(e["RH"])
    mo = {}
    for m, d in sorted(monthly.items()):
        mo[m] = {"n": len(d["g"]), "gamma2_median": float(np.median(d["g"])),
                 "n_T": len(d["T"]),
                 "T_mean": float(np.mean(d["T"])) if d["T"] else None,
                 "n_RH": len(d["RH"]),
                 "RH_mean": float(np.mean(d["RH"])) if d["RH"] else None}
    vars_ = {}
    gy = np.asarray(g_vals, float)
    for tag, x in (("T", np.asarray(t_vals, float)),
                   ("RH", np.asarray(rh_vals, float))):
        unit = [u for _n, t, u in core.ENV_VARS if t == tag][0]
        vars_[tag] = {"unit": unit, "n": int(x.size),
                      "x_min": float(x.min()) if x.size else None,
                      "x_max": float(x.max()) if x.size else None,
                      "points": [[float(a), float(b)] for a, b in zip(x, gy)],
                      "fit": (st.env_fit(x, gy) if x.size else None)}
    return {"feature": core.GAMMA2_COLUMN, "interpolation": interp_string(geom),
            "n_dates_with_gamma2": sum(1 for r in rows
                                       if r[core.GAMMA2_COLUMN] is not None),
            "n_dates_with_env": n_env, "vars": vars_, "monthly": mo}


def persistence_block(rows):
    """The ``persistence`` sub-block of the channel-table meta file.

    ``masks.add_persistence`` builds its frequency map in a majority chip shape
    and marks the deviating dates, so exactly these three numbers describe it.
    """
    shapes = Counter(tuple(r["mask"].shape) for r in rows)
    majority, n_majority = shapes.most_common(1)[0]
    return {"majority_shape": list(majority), "n_majority_shape": n_majority,
            "n_excluded_other_shape": len(rows) - n_majority}


# The committed ``mask_definition`` block, rebuilt from the mask module.
MASK_DEF_KEYS = {
    "anchor_half": masks.ANCHOR_HALF,
    "band_half": masks.BAND_HALF,
    "center_col": masks.CENTER_COL,
    "center_row": masks.CENTER_ROW,
    "col_half": masks.COL_HALF,
    "contrast_detect": masks.CONTRAST_DETECT,
    "house_median_mult": masks.GAMMA2_MEDIAN_MULT,
    "house_min_masked": masks.GAMMA2_MIN_MASKED,
    "house_peak_frac": masks.GAMMA2_PEAK_FRAC,
    "k_local": masks.K_LOCAL,
    "layer": "deck_edge",
    "low_n_gate": masks.LOW_N_GATE,
    "min_n_masked": masks.MIN_N_MASKED,
}


def compute(quick=False):
    """Recompute every plotted number and every pinned block from the CSV."""
    rows = core.load_csv()
    by = core.by_series(rows)
    ref = core.load_reference()
    meta = core.load_csv_meta()
    env = core.load_meteo()
    mask_ref = ref["deck_edge_mask"]
    state_ref = ref["deck_edge_state"]

    # ---- the two committed reference blocks, rebuilt from the table --------
    mask_series = {lab: series_mask_block(by[lab]) for lab in core.SERIES_ORDER}
    state_series = {lab: series_state_block(by[lab], lab) for lab in core.SERIES_ORDER}
    disc = {}
    for g in ("ASC", "DESC"):
        disc[g] = masks.discrimination(
            {s: [r for r in by[core.BRIDGE_BY_GEOM[g]] if r["state"] == s]
             for s in core.STATES},
            {s: [r for r in by[core.CONTROL_BY_GEOM[g]] if r["state"] == s]
             for s in core.STATES})
    steps, placebo, seasons = {}, {}, {}
    g2_steps, g2_placebo, g2_seasons, environ = {}, {}, {}, {}
    for lab in core.SERIES_ORDER:
        rs = by[lab]
        steps[lab] = step_block(rs, "band_contrast")
        seasons[lab] = season_block(rs, "band_contrast", "contrast")
        g2_seasons[lab] = season_block(rs, core.GAMMA2_COLUMN, core.GAMMA2_COLUMN)
        environ[lab] = environment_block(rs, lab, env)
        g2_steps[lab], g2_placebo[lab] = gamma2_variant_blocks(rs, quick=quick)
        if not quick:
            placebo[lab] = placebo_block(rs, "band_contrast")
    dvc = delta_vs_ctrl(by, "band_contrast", quick=quick)
    g2_dvc = delta_vs_ctrl(by, core.GAMMA2_COLUMN, with_feature=True, quick=quick)

    # ---- --quick: the placebo sweeps are the heavy part of a run ----------
    copied = []
    if quick:
        placebo = copy.deepcopy(state_ref["placebo"])
        g2_placebo = copy.deepcopy(state_ref["gamma2_placebo"])
        for g in ("ASC", "DESC"):
            dvc[g]["placebo"] = copy.deepcopy(state_ref["delta_vs_ctrl"][g]["placebo"])
            g2_dvc[g]["placebo"] = copy.deepcopy(
                state_ref["gamma2_delta_vs_ctrl"][g]["placebo"])
        copied = list(QUICK_COPIED) + ["delta_vs_ctrl[*].placebo",
                                       "gamma2_delta_vs_ctrl[*].placebo"]

    # ---- provenance of the inputs + cross-checks of the meta file ---------
    data = {
        "site": meta["site"],
        "csv": os.path.relpath(core.CSV_PATH),
        "csv_meta": os.path.relpath(core.CSV_META_PATH),
        "csv_sha256": core.sha256(core.CSV_PATH),
        "meteo": os.path.relpath(core.METEO_PATH),
        "n_rows": len(rows),
        "n_by_series": {lab: len(by[lab]) for lab in core.SERIES_ORDER},
        "n_by_state": {s: sum(1 for r in rows if r["state"] == s)
                       for s in core.STATES},
        "n_by_series_state": {lab: dict(state_series[lab]["n_by_state"])
                              for lab in core.SERIES_ORDER},
        "rect_files": meta["rect_files"],
        "rect_file_sha256": meta["rect_file_sha256"],
        "reference_source": meta["reference_source"],
        "mask_definition": meta["mask_definition"],
        "states": {s: meta["states"][s]["window"] for s in core.STATES},
        "cuts": {"DS1": core.CUT_DS1.isoformat(), "DS2": core.CUT_DS2.isoformat()},
        "environment_source": config_state()["environment_source"],
    }
    checks = []

    def check(name, got, want):
        entry = {"check": name, "ok": bool(got == want)}
        if not entry["ok"]:
            entry["recomputed"] = repr(got)[:300]
            entry["committed"] = repr(want)[:300]
        checks.append(entry)

    for lab in core.SERIES_ORDER:
        m = meta["series"][lab]
        check(f"meta.series.{lab}.file", mask_series[lab]["file"], m["file"])
        check(f"meta.series.{lab}.geom", state_series[lab]["geom"], m["geom"])
        check(f"meta.series.{lab}.n_dates", state_series[lab]["n_dates"], m["n_dates"])
        check(f"meta.series.{lab}.n_by_state", state_series[lab]["n_by_state"],
              m["n_by_state"])
        check(f"meta.series.{lab}.dates_by_state", state_series[lab]["dates_by_state"],
              m["dates_by_state"])
        check(f"meta.series.{lab}.summary", mask_series[lab]["summary"], m["summary"])
        for key in ("union50", "union30"):
            check(f"meta.series.{lab}.{key}",
                  [list(x) for x in mask_series[lab][key]], m[key])
        check(f"meta.series.{lab}.persistence", persistence_block(by[lab]),
              m["persistence"])
    check("meta.mask_definition", dict(MASK_DEF_KEYS), meta["mask_definition"])
    for fname, want in meta["rect_file_sha256"].items():
        check(f"rect_sha256.{fname}",
              core.sha256(os.path.join(core.DATA, fname)), want)
    for kind, mkey in (("deck_edge_mask", "mask"), ("deck_edge_state", "state")):
        check(f"reference_source.{kind}", core.REFERENCE_FILES[kind],
              os.path.basename(meta["reference_source"][mkey]))
    data["cross_checks"] = checks
    data["n_cross_checks_failed"] = sum(1 for c in checks if not c["ok"])

    # ---- figure payloads (every number the five figures draw) -------------
    plot_order = [core.BRIDGE_BY_GEOM["ASC"], core.BRIDGE_BY_GEOM["DESC"],
                  core.CONTROL_BY_GEOM["ASC"], core.CONTROL_BY_GEOM["DESC"]]
    fig_a = {"order": plot_order, "channels": {}, "descriptives": {}}
    for key in CHANNELS:
        per_series = {lab: {s: [r[key] for r in by[lab]
                                if r["state"] == s and r[key] is not None]
                            for s in core.STATES} for lab in core.SERIES_ORDER}
        fig_a["channels"][key] = per_series
        fig_a["descriptives"][key] = {
            lab: {s: st._stats(per_series[lab][s]) for s in core.STATES}
            for lab in core.SERIES_ORDER}

    fig_b = {"contrast_detect": masks.CONTRAST_DETECT, "order": plot_order,
             "discrimination": disc, "mask_summary": {}, "visibility": {},
             "band_contrast": {}, "n_masked": {}, "gamma2_effects": {}}
    for lab in core.SERIES_ORDER:
        s = mask_series[lab]["summary"]
        p, lo, hi = st.wilson(s["n_dates_with_mask"], s["n_dates"])
        fig_b["mask_summary"][lab] = s
        fig_b["visibility"][lab] = {"n": s["n_dates"],
                                    "n_with_mask": s["n_dates_with_mask"],
                                    "share": p, "wilson": [lo, hi]}
        fig_b["band_contrast"][lab] = {s_: [r["band_contrast"] for r in by[lab]
                                            if r["state"] == s_]
                                       for s_ in core.STATES}
        fig_b["n_masked"][lab] = {s_: [r["n_masked"] for r in by[lab]
                                       if r["state"] == s_]
                                  for s_ in core.STATES}
        grp = st.by_state(by[lab])
        fig_b["gamma2_effects"][lab] = {
            ds: st.delta_block(st.vals(grp, st.REFERENCE_STATE, key=core.GAMMA2_COLUMN),
                               st.vals(grp, ds, key=core.GAMMA2_COLUMN),
                               f"RS_vs_{ds}", n_boot=st.BOOT_STRATA,
                               seed=st.BOOT_SEED)
            for ds in st.DAMAGED}

    fig_c = {"order": plot_order, "placebo_weeks": list(st.PLACEBO_WEEKS),
             "step_tests": steps, "placebo": placebo,
             "gamma2_step_tests": g2_steps, "gamma2_placebo": g2_placebo}

    pairs = {col: paired_diff_series(by, col)
             for col in ("band_contrast", core.GAMMA2_COLUMN)}
    fig_d = {"order": [("ASC", core.BRIDGE_BY_GEOM["ASC"]),
                       ("DESC", core.BRIDGE_BY_GEOM["DESC"])],
             "delta_vs_ctrl": dvc, "gamma2_delta_vs_ctrl": g2_dvc,
             "pairs": pairs}

    fig_e = {"order": plot_order, "seasonality": seasons,
             "gamma2_seasonality": g2_seasons, "gamma2_environment": environ,
             "months": {lab: sorted(seasons[lab]) for lab in core.SERIES_ORDER}}

    # ---- one result dict: figures + the state block + provenance ----------
    state = {"config": config_state(), "series": state_series,
             "discrimination": disc, "step_tests": steps, "placebo": placebo,
             "delta_vs_ctrl": dvc, "seasonality": seasons,
             "gamma2_step_tests": g2_steps, "gamma2_placebo": g2_placebo,
             "gamma2_delta_vs_ctrl": g2_dvc, "gamma2_seasonality": g2_seasons,
             "gamma2_environment": environ}
    mask = {"config": dict(masks.CONFIG), "series": mask_series}
    result = {
        "figure": "fig_bautzen_ocv_paper",
        "site": meta["site"],
        "quick": bool(quick), "copied_blocks": copied,
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": getattr(scipy, "__version__", None),
        "channels": CHANNELS, "chance": CHANCE, "states": list(core.STATES),
        "state_windows": {s: list(core.STATE_WINDOWS[s]) for s in core.STATES},
        "data": data, "mask": mask, "state": state,
        "fig_a": fig_a, "fig_b": fig_b, "fig_c": fig_c, "fig_d": fig_d,
        "fig_e": fig_e,
    }
    return result, mask_ref, state_ref


# ---------------------------------------------------------------------------
# Pin block — recomputed against the two committed reference JSONs
# ---------------------------------------------------------------------------
# The channel table is the single input of this script; the two reference JSONs
# are the committed outputs of the independent Bautzen analysis. The pin block
# compares the whole committed structure (mask config + 4 series blocks; state
# config + 4 series blocks + the nine cross-series blocks) leaf by leaf.
def pin_all(result, mask_ref, state_ref):
    """Compare every recomputed block against the two committed references."""
    pins = Pins()
    pins.cmp("mask", result["mask"], mask_ref)
    pins.cmp("state", result["state"], state_ref)
    summary = pins.summary()
    # The channel-table meta file is a second, independent witness of the same
    # table (written by the generator, not by this script).
    failed = [c for c in result["data"]["cross_checks"] if not c["ok"]]
    summary["cross_checks"] = {"n": len(result["data"]["cross_checks"]),
                               "n_failed": len(failed), "failed": failed}
    summary["ok"] = bool(summary["n_failed"] == 0 and not failed)
    return summary


def json_default(o):
    """JSON fallback for the numpy scalars/arrays coming out of the statistics."""
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.floating):
        return float(o)
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    raise TypeError(f"not JSON serialisable: {type(o)}")


# ---------------------------------------------------------------------------
# Rendering — figures A–E, 183 mm column width, 600 dpi PNG
# ---------------------------------------------------------------------------
def series_tick(lab):
    """Short tick label of one of the four series: 'bridge ASC', 'control DESC'."""
    role = "bridge" if core.is_bridge(lab) else "control"
    return f"{role} {core.geom_of(lab)}"


def series_state_dots(ax, values, order, seed0=10, span=0.13, fs=6.0):
    """Dots per overpass (state-coloured) + median bar: three states per series.

    ``values`` is ``{series: {state: [value, ...]}}`` (the ``fig_a.channels``
    layout of the result JSON).
    """
    for gi, lab in enumerate(order):
        for si, s in enumerate(core.STATES):
            v = np.asarray(values[lab][s], float)
            v = v[np.isfinite(v)]
            x0 = gi + (si - 1) * 0.22
            if v.size:
                ax.plot(x0 + jitter(v.size, span, seed0 + 7 * gi + si), v,
                        STATE_MARKERS[s], ms=2.0, mfc="none", mec=STATE_COLORS[s],
                        mew=0.5, alpha=0.65, ls="none")
                ax.hlines(float(np.median(v)), x0 - 0.085, x0 + 0.085,
                          color=STATE_COLORS[s], lw=1.6)
    ax.axvline(1.5, color="0.85", lw=0.8, zorder=0)
    ax.set_xlim(-0.62, len(order) - 0.38)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([series_tick(l) for l in order], fontsize=fs, rotation=22,
                       ha="right", rotation_mode="anchor")


def state_legend(fig, y=-0.01, extra=()):
    """Figure-level legend: the three states (+ optionally other proxy handles)."""
    handles = [plt.Line2D([], [], marker=STATE_MARKERS[s], ls="none", mfc="none",
                          mec=STATE_COLORS[s], mew=0.8, ms=3.4, label=s)
               for s in core.STATES] + list(extra)
    fig.legend(handles=handles, loc="lower center", ncol=len(handles),
               frameon=False, bbox_to_anchor=(0.5, y))


def fig_a(res, out_dir=None):
    """Figure A — the six OCV components per series and structural state."""
    fa = res["fig_a"]
    order = fa["order"]
    fig, axes = plt.subplots(2, 3, figsize=(COL_WIDTH, 0.82 * COL_WIDTH))
    for ax, key in zip(axes.ravel(), CHANNELS):
        series_state_dots(ax, fa["channels"][key], order,
                          seed0=10 + 7 * CHANNELS.index(key))
        if key == "A":
            ax.set_yscale("log")
        elif key == "P":
            ax.set_ylim(0.0, 1.02)
        elif key == "gamma2":
            ax.axhline(0.0, color="0.6", lw=0.8, ls=":")
        ax.set_ylabel(CHANNEL_AXIS[key], fontsize=7.0)
        panel(ax, "abcdef"[CHANNELS.index(key)], CHANNEL_TITLE[key])
    n_by_state = res["data"]["n_by_state"]
    fig.suptitle(
        "Figure A — the six OCV components per overpass, per rect series and "
        "structural state  "
        f"({res['data']['n_rows']} overpasses: "
        + ", ".join(f"{s} n={n_by_state[s]}" for s in core.STATES)
        + ")\n"
        "dots = single overpass (state-coloured), bar = median of the series "
        "and state;  input: data/bautzen/bautzen_ocv_channels.csv",
        fontsize=7.6)
    state_legend(fig, y=0.002)
    fig.tight_layout(rect=(0, 0.035, 1, 0.90))
    return save_fig(fig, "fig_bautzen_ocv_paper_A", out_dir)


def fig_b(res, out_dir=None):
    """Figure B — the deck-edge mask layer: is the edge visible at all?"""
    fb = res["fig_b"]
    order = fb["order"]
    thr = fb["contrast_detect"]
    fig, axes = plt.subplots(2, 2, figsize=(COL_WIDTH, 0.78 * COL_WIDTH))
    ax = axes[0][0]
    series_state_dots(ax, fb["band_contrast"], order, seed0=41, fs=5.8)
    ax.axhline(thr, color="0.4", lw=0.9, ls="--")
    ax.text(0.02, thr, f" detect threshold {thr:g}", fontsize=5.6, va="bottom",
            color="0.35", transform=ax.get_yaxis_transform())
    ax.set_ylabel("band contrast", fontsize=7.0)
    panel(ax, "a", "deck-band contrast per series and state")

    ax = axes[0][1]
    vis = fb["visibility"]
    ks = np.arange(len(order))
    ax.bar(ks, [vis[l]["share"] for l in order], width=0.6,
           color=[SERIES_COLORS[l] for l in order], alpha=0.75)
    for i, l in enumerate(order):
        lo, hi = vis[l]["wilson"]
        ax.plot([i, i], [lo, hi], color="0.25", lw=1.0)
        ax.text(i, min(1.03, hi + 0.02), f"{vis[l]['n_with_mask']}/{vis[l]['n']}",
                fontsize=5.6, ha="center")
    ax.set_ylim(0.0, 1.05)
    ax.set_xticks(ks)
    ax.set_xticklabels([series_tick(l) for l in order], fontsize=5.8, rotation=22,
                       ha="right", rotation_mode="anchor")
    ax.set_ylabel("share of dates with a mask", fontsize=7.0)
    panel(ax, "b", "mask visibility (Wilson 95% CI)")

    ax = axes[1][0]
    series_state_dots(ax, {l: {s: [v + 1.0 for v in fb["n_masked"][l][s]]
                               for s in core.STATES} for l in order},
                      order, seed0=61, fs=5.8)
    ax.set_yscale("log")
    ax.set_ylabel("n_masked + 1", fontsize=7.0)
    panel(ax, "c", "deck-edge pixel count per series and state (+1, log)")

    ax = axes[1][1]
    eff = fb["gamma2_effects"]
    w = 0.34
    for si, ds in enumerate(st.DAMAGED):
        ax.barh(ks + (si - 0.5) * w, [eff[l][ds]["delta"] for l in order],
                height=w, color=STATE_COLORS[ds], alpha=0.85, label=f"RS vs {ds}")
    for i, l in enumerate(order):
        for si, ds in enumerate(st.DAMAGED):
            b = eff[l][ds]
            ax.text(b["delta"], i + (si - 0.5) * w,
                    f" {b['delta']:+.2f} (SDS {b['sds']:.1f})", fontsize=5.2,
                    va="center", ha="left" if b["delta"] >= 0 else "right")
    ax.axvline(0.0, color="0.4", lw=0.8)
    lim = max(abs(eff[l][ds]["delta"]) for l in order for ds in st.DAMAGED) or 1.0
    ax.set_xlim(-1.65 * lim, 1.65 * lim)
    ax.set_yticks(ks)
    ax.set_yticklabels([series_tick(l) for l in order], fontsize=5.8)
    ax.set_xlabel("Cliff's delta  (RS as reference)", fontsize=7.0)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.30), ncol=2,
              fontsize=5.6, frameon=False)
    panel(ax, "d", "gamma^2 effect size within each series")
    fig.suptitle(
        "Figure B — the deck-edge mask layer: the bridge edge is detectable in "
        "every state, the control edge is not\n"
        "band contrast = deck-band median / local median;  detect rule of the "
        f"origin analysis = contrast >= {thr:g} and >= "
        f"{mask_definition('min_n_masked')} deck-edge pixels", fontsize=7.6)
    state_legend(fig, y=0.002)
    fig.tight_layout(rect=(0, 0.04, 1, 0.89))
    return save_fig(fig, "fig_bautzen_ocv_paper_B", out_dir)


def placebo_delta_panel(ax, keys, blocks, cut, ylab, colours, ticks, fs=5.8):
    """Observed step delta (diamond) against the shifted-cut placebos of a series.

    ``blocks`` = ``{series: {cut: {"delta", "p", "placebo": [delta, ...]}}}``.
    The grey tokens are the placebo cuts of +/-2..12 weeks, the bar spans them;
    ``k/n`` counts the placebo cuts whose |delta| reaches the observed one.
    """
    for gi, lab in enumerate(keys):
        b = blocks[lab][cut]
        pl = b["placebo"]
        if pl:
            ax.plot([gi, gi], [min(pl), max(pl)], color="0.72", lw=4.0, zorder=1,
                    solid_capstyle="butt")
            ax.plot([gi] * len(pl), pl, "o", ms=2.2, mfc="none", mec="0.5",
                    mew=0.5, ls="none", zorder=2)
        ax.plot([gi], [b["delta"]], "D", ms=4.2, color=colours[lab], zorder=3)
        top = max([b["delta"]] + pl) if pl else b["delta"]
        k = sum(1 for x in pl if abs(x) >= abs(b["delta"]))
        ax.text(gi, top, f"p={b['p']:.3f}\n{k}/{len(pl)}" if pl else
                f"p={b['p']:.3f}", fontsize=5.0, ha="center", va="bottom")
    ax.axhline(0.0, color="0.5", lw=0.8, ls=":")
    ax.set_xlim(-0.6, len(keys) - 0.4)
    ax.set_xticks(range(len(keys)))
    ax.set_xticklabels(ticks, fontsize=fs, rotation=22, ha="right",
                       rotation_mode="anchor")
    ax.set_ylabel(ylab, fontsize=7.0)
    # headroom for the two-line ``p=.../k/n`` annotation above the placebos, so
    # it cannot grow into the (left-aligned, two-line) panel caption
    ax.autoscale_view(scalex=False)
    lo, hi = ax.get_ylim()
    ax.set_ylim(lo, hi + 0.15 * (hi - lo))


def fig_c(res, out_dir=None):
    """Figure C — the cut dates against cuts shifted by +/-2..12 weeks."""
    fc = res["fig_c"]
    order = fc["order"]
    g2 = core.GAMMA2_COLUMN

    def blocks(steps, placebos):
        return {l: {c: {"delta": steps[l][c]["delta"],
                        "p": steps[l][c]["p_perm"],
                        "placebo": [p["delta"] for p in
                                    (placebos[l][c]["placebo"] or [])]}
                    for c in ("DS1", "DS2")} for l in order}

    rows = (
        ("band contrast", "band contrast",
         blocks(fc["step_tests"], fc["placebo"])),
        ("gamma^2 (deck band, floor-corrected)", "gamma^2",
         blocks({l: fc["gamma2_step_tests"][l][g2] for l in order},
                {l: fc["gamma2_placebo"][l][g2] for l in order})),
    )
    fig, axes = plt.subplots(2, 2, figsize=(COL_WIDTH, 0.80 * COL_WIDTH))
    for ri, (ylab, short, blk) in enumerate(rows):
        for ci, (cut, cutname) in enumerate((("DS1", "DS1 (12 May 2025)"),
                                             ("DS2", "DS2 (29 Sep 2025)"))):
            ax = axes[ri][ci]
            placebo_delta_panel(ax, order, blk, cut, ylab, SERIES_COLORS,
                                [series_tick(l) for l in order])
            panel(ax, "abcd"[2 * ri + ci],
                  f"{cutname}: median jump (after - before), {short}")
    fig.legend(handles=[
        plt.Line2D([], [], marker="D", ls="none", color="0.3", ms=4.0,
                   label="observed cut"),
        plt.Line2D([], [], marker="o", ls="none", mfc="none", mec="0.5", ms=2.6,
                   label="placebo cut (+/-2..12 weeks)")],
        loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(0.5, -0.004))
    fig.suptitle(
        "Figure C — is the jump at the cut date special?  Each observed cut "
        "against cuts shifted by +/-2..12 weeks\n"
        "grey = placebo cuts of the same series (bar = their range), "
        "p = two-sided permutation p of the observed step", fontsize=7.6)
    fig.tight_layout(rect=(0, 0.04, 1, 0.88))
    return save_fig(fig, "fig_bautzen_ocv_paper_C", out_dir)


# Bridge-vs-control reading: the two differences are per orbit direction, so the
# geom (not the series) is the natural key of figure D.
GEOM_COLORS = {"ASC": "#4C72B0", "DESC": "#55A868"}
DAY0 = dt.date(2025, 1, 1)          # x=0 of the figure-D time axis


def _day(d):
    """Day number (since 2025-01-01) of an ISO date string or ``datetime.date``."""
    if isinstance(d, str):
        d = dt.date.fromisoformat(d)
    return (d - DAY0).days


def state_windows_panel(ax, x0, x1, label_spans=True, fs=5.4):
    """Shade the three structural states on a time axis of day numbers."""
    bounds = [x0, _day(core.CUT_DS1.date()), _day(core.CUT_DS2.date()), x1]
    for i, s in enumerate(core.STATES):
        lo, hi = bounds[i], bounds[i + 1]
        if hi <= lo:
            continue
        ax.axvspan(lo, hi, color=STATE_COLORS[s], alpha=0.07, zorder=0)
        if label_spans:
            ax.text((lo + hi) / 2.0, 0.97, s, transform=ax.get_xaxis_transform(),
                    fontsize=fs, color=STATE_COLORS[s], ha="center", va="top")


def state_effect_panel(ax, blocks, states, ylab):
    """Median + IQR of a paired effect per state, one marker per geom."""
    ks = np.arange(len(states))
    for gi, geom in enumerate(blocks):
        blk = blocks[geom]
        off = (gi - 0.5) * 0.26
        for i, s in enumerate(states):
            b = blk["states"][s]
            if not b:
                continue
            ax.vlines(i + off, b["iqr"][0], b["iqr"][1], color=GEOM_COLORS[geom],
                      lw=2.4, alpha=0.55)
            ax.plot([i + off], [b["median"]], "D", ms=4.0, color=GEOM_COLORS[geom],
                    label=geom if i == 0 else None)
            ax.text(i + off, b["iqr"][1], f"{b['median']:+.2f}\nn={b['n']}",
                    fontsize=4.6, ha="center", va="bottom")
    ax.axhline(0.0, color="0.5", lw=0.8, ls=":")
    ax.set_xticks(ks)
    ax.set_xticklabels(states, fontsize=6.4)
    ax.set_xlim(-0.55, len(states) - 0.45)
    ax.set_ylabel(ylab, fontsize=7.0)
    ax.legend(loc="best", fontsize=6.0, frameon=False)


def fig_d(res, out_dir=None):
    """Figure D — bridge minus control: is the bridge-specific jump visible?"""
    fd = res["fig_d"]
    fig, axes = plt.subplots(2, 2, figsize=(COL_WIDTH, 0.78 * COL_WIDTH))
    state_effect_panel(axes[0][0], fd["delta_vs_ctrl"], core.STATES,
                       "band contrast difference")
    panel(axes[0][0], "a", "bridge - control, band contrast, per state")
    state_effect_panel(axes[0][1], fd["gamma2_delta_vs_ctrl"], core.STATES,
                       "gamma^2 difference")
    panel(axes[0][1], "b", "bridge - control, gamma^2")

    ax = axes[1][0]
    pairs = fd["pairs"]["band_contrast"]
    all_x = [_day(d) for p in pairs.values() for d in p["dates"]]
    for geom in pairs:
        p = pairs[geom]
        xs = [_day(d) for d in p["dates"]]
        ax.plot(xs, p["diffs"], "-", color=GEOM_COLORS[geom], lw=0.7, alpha=0.4)
        for s in core.STATES:
            ax.plot([x for x, ss in zip(xs, p["state"]) if ss == s],
                    [y for y, ss in zip(p["diffs"], p["state"]) if ss == s],
                    STATE_MARKERS[s], ms=2.8, mfc="none", mec=GEOM_COLORS[geom],
                    mew=0.6, ls="none", alpha=0.85)
    ax.axhline(0.0, color="0.5", lw=0.8, ls=":")
    ax.set_xlim(min(all_x) - 4, max(all_x) + 4)
    state_windows_panel(ax, min(all_x) - 4, max(all_x) + 4)
    ax.set_ylabel("band contrast difference", fontsize=7.0)
    ax.set_xlabel("day of year 2025", fontsize=7.0)
    panel(ax, "c", "per-date difference (colour = geom, marker = state)")

    ax = axes[1][1]
    blocks = {g: {c: {"delta": fd["delta_vs_ctrl"][g]["step"][c]["delta"],
                      "p": fd["delta_vs_ctrl"][g]["step"][c]["p_perm"],
                      "placebo": [p["delta"] for p in
                                  (fd["delta_vs_ctrl"][g]["placebo"][c]["placebo"]
                                   or [])]}
                  for c in ("DS1", "DS2")} for g in fd["delta_vs_ctrl"]}
    placebo_delta_panel(ax, ["ASC", "DESC"], blocks, "DS1", "difference jump",
                        GEOM_COLORS, ["ASC", "DESC"])
    panel(ax, "d", "step/placebo of the difference at DS1")

    fig.suptitle(
        "Figure D — bridge minus control, date by date: the two orbit directions "
        "are the control and the counter-check\n"
        "dots = state median with the 25/75 percentile bar;  "
        "shaded bands of panel (c) = RS / DS1 / DS2, dashed line = DS1", fontsize=7.6)
    fig.tight_layout(rect=(0, 0.01, 1, 0.89))
    return save_fig(fig, "fig_bautzen_ocv_paper_D", out_dir)


def fig_e(res, out_dir=None):
    """Figure E — seasonality and environment: the confounders of the states."""
    fe = res["fig_e"]
    order = fe["order"]
    fig, axes = plt.subplots(2, 2, figsize=(COL_WIDTH, 0.78 * COL_WIDTH))

    rows = (("seasonality", "contrast_median", "monthly median band contrast"),
            ("gamma2_seasonality", f"{core.GAMMA2_COLUMN}_median",
             "monthly median gamma^2"))
    for ri, (key, field, ylab) in enumerate(rows):
        ax = axes[0][ri]
        months = sorted({m for lab in order for m in fe[key][lab]})
        xs = np.arange(len(months))
        for lab in order:
            blk = fe[key][lab]
            ax.plot(xs, [blk[m][field] if m in blk else np.nan for m in months],
                    "-", color=SERIES_COLORS[lab], lw=1.0,
                    marker=SERIES_MARKERS[lab], ms=2.8, alpha=0.9,
                    label=series_tick(lab))
        ax.set_xticks(xs)
        ax.set_xticklabels(months, fontsize=5.6, rotation=35, ha="right")
        ax.set_ylabel(ylab, fontsize=7.0)
        panel(ax, "ab"[ri], ylab)

    for ci, tag in enumerate(("T", "RH")):
        ax = axes[1][ci]
        lines = []
        for lab in order:
            var = fe["gamma2_environment"][lab]["vars"][tag]
            pts = np.asarray(var["points"], float)
            if pts.size:
                ax.plot(pts[:, 0], pts[:, 1], SERIES_MARKERS[lab], ms=2.6,
                        mfc="none", mec=SERIES_COLORS[lab], mew=0.5, ls="none",
                        alpha=0.85)
            fit = var["fit"]
            if fit:
                gx = np.array([var["x_min"], var["x_max"]])
                ax.plot(gx, fit["intercept"] + fit["slope"] * gx,
                        color=SERIES_COLORS[lab], lw=1.0)
                lines.append(f"{series_tick(lab)}: r={fit['r']:+.2f} "
                             f"(p={fit['p']:.2f}, n={fit['n']})")
        ax.text(0.02, 0.98, "\n".join(lines), transform=ax.transAxes, fontsize=4.8,
                va="top", ha="left")
        unit = fe["gamma2_environment"][order[0]]["vars"][tag]["unit"]
        name = "air temperature" if tag == "T" else "relative humidity"
        ax.set_xlabel(f"{name} at the overpass instant  [{unit}]", fontsize=7.0)
        ax.set_ylabel("gamma^2 (deck band, floor-corrected)", fontsize=7.0)
        panel(ax, "cd"[ci], f"gamma^2 vs. {name}")

    handles = [plt.Line2D([], [], color=SERIES_COLORS[l], marker=SERIES_MARKERS[l],
                          ls="-", lw=1.0, ms=2.8, label=series_tick(l))
               for l in order]
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False,
               bbox_to_anchor=(0.5, -0.004))
    fig.suptitle(
        "Figure E — seasonality and environment proxies: state and season are "
        "confounded (one deployment)\n"
        "top = median of the month, bottom = one dot per overpass "
        "(line = OLS fit); the control series is the counter-check", fontsize=7.6)
    fig.tight_layout(rect=(0, 0.04, 1, 0.89))
    return save_fig(fig, "fig_bautzen_ocv_paper_E", out_dir)


# ---------------------------------------------------------------------------
# Result JSON + English Markdown report
# ---------------------------------------------------------------------------
def build_json(res, pins, quick):
    """Self-contained result JSON: everything the figures plot, plus the pins."""
    out = dict(res)
    out["meta"] = {
        "generator": "code/bautzen/fig_bautzen_ocv_paper.py",
        "quick": bool(quick),
        "copied_blocks": res["copied_blocks"],
        "figures": [f"figures/bautzen/fig_bautzen_ocv_paper_{k}.png"
                    for k in "ABCDE"],
        "report": "figures/bautzen/fig_bautzen_ocv_paper.md",
        "versions": {"python": platform.python_version(), "numpy": np.__version__,
                     "matplotlib": matplotlib.__version__,
                     "scipy": getattr(scipy, "__version__", None) if scipy else None},
        "seeds": {"rng_seed": st.RNG_SEED, "boot_seed": st.BOOT_SEED,
                  "n_perm": st.N_PERM, "n_boot": st.N_BOOT},
        "note": "no timestamps on purpose: repeated runs must be identical",
    }
    out["pins"] = pins
    return out


def report(res, pins, quick):
    """The English Markdown report of the five figures, built from the result."""
    L = []
    A = L.append
    d = res["data"]
    mask, state = res["mask"], res["state"]
    mask_series, series = mask["series"], state["series"]
    disc, steps = state["discrimination"], state["step_tests"]
    placebo = state["placebo"]
    dvc, g2_dvc = state["delta_vs_ctrl"], state["gamma2_delta_vs_ctrl"]
    seasons, g2_seasons = state["seasonality"], state["gamma2_seasonality"]
    environ = state["gamma2_environment"]
    g2_steps, g2_placebo = state["gamma2_step_tests"], state["gamma2_placebo"]
    order = res["fig_a"]["order"]
    g2 = core.GAMMA2_COLUMN
    variants = [v[0] for v in G2_VARIANTS]

    A("# Bautzen deck-edge OCV — figures A–E")
    A("")
    A("Recomputation of the five figures of the Bautzen deck-edge observability "
      "package from the committed channel table. The table is the only input; "
      "every number below is recomputed here and pinned against the two "
      "committed reference JSONs (section 5).")
    A("")
    n_by_state = ", ".join(f"{s} n={d['n_by_state'][s]}" for s in core.STATES)
    A(md_table(["", ""], [
        ["input", f"`{d['csv']}` (sha256 `{d['csv_sha256'][:16]}…`)"],
        ["rows", f"{d['n_rows']} overpasses in 4 rect series ({n_by_state})"],
        ["mask layer", f"`{os.path.basename(core.REFERENCE_FILES['deck_edge_mask'])}` "
                       "and `bautzen_ocv_channels_meta.json` (rect-file sha256s verified)"],
        ["mode", ("`--quick` (placebo sweeps copied from the reference)"
                  if quick else "full recomputation")],
        ["figures", "`figures/bautzen/fig_bautzen_ocv_paper_A..E.png` (600 dpi)"],
        ["result JSON", "`data/bautzen/fig_bautzen_ocv_paper.json`"],
        ["pins", f"{pins['n_checks']} checks, {pins['n_failed']} failures, "
                 f"{pins['label_translations']} doc-string translations"],
        ["cross-checks", f"{len(d['cross_checks'])} meta checks, "
                         f"{d['n_cross_checks_failed']} failed"],
    ]))
    A("")

    # ---- 1. question and design ------------------------------------------
    A("## 1. Question and design")
    A("")
    A("The site is a pedestrian bridge (`bautzen_openbridge`): a bare deck edge "
      "crosses the SAR pixel grid of the InSAR OCV (observability core vector), "
      "so a bright, geometrically fixed line is expected in the coherence image. "
      "The deck edge was modified twice in the deployment "
      f"(`{core.CUT_DS1:%d %b %Y}` and `{core.CUT_DS2:%d %b %Y}`); the three "
      "resulting states are RS (before, "
      f"{d['states']['RS']}), DS1 ({d['states']['DS1']}) and "
      f"DS2 ({d['states']['DS2']}).")
    A("")
    A("Two design properties decide how the numbers below may be read:")
    A("")
    A("1. **Four series, two of them counter-checks.** The same geometry was "
      "analysed for both orbit directions and for a control rectangle 500 m east "
      "of the bridge. Only a difference that appears in the bridge series and "
      "*not* in the control is evidence about the deck edge; the control fixes "
      "orbit, season and acquisition effects.")
    A("2. **State and season are confounded.** RS, DS1 and DS2 are time "
      "intervals of one deployment, so a state effect and a seasonal effect are "
      "the same contrast unless the control series moves with the bridge. The "
      f"single cut dates are therefore tested against cuts shifted by "
      f"+/-{', '.join(str(w) for w in st.PLACEBO_WEEKS)} weeks (section C), and "
      "the season is probed directly with air temperature / humidity at the "
      "overpass instant (figure E).")
    A("")
    A("The vector examined here is the site's coherence channel `gamma^2` of the "
      "deck band *plus* the five mask-morphology channels "
      f"`{'`, `'.join(core.DIMS)}` of `bautzen_ocv_masks.py` (area, density, "
      "fragmentation, peak shift, persistence). The mask layer is a documented "
      "deviation from the LUMO package (section 6).")
    A("")

    # ---- 2. data ----------------------------------------------------------
    A("## 2. Data")
    A("")
    A("One row per overpass and rect series (the committed channel table). "
      "The mask share counts the dates on which the deck-edge mask rule of the "
      "origin analysis fires (band contrast and the deck-edge pixel count above "
      "the detection gate).")
    A("")
    rows = []
    for lab in order:
        s = mask_series[lab]["summary"]
        m = series[lab]
        rows.append([series_label(lab), m["n_dates"],
                     "/".join(str(m["n_by_state"][c]) for c in core.STATES),
                     f"{s['n_dates_with_mask']}/{s['n_dates']}",
                     fmt(s["band_contrast_median"], 2),
                     fmt(s["frac_dates_contrast_gt_thr"], 2),
                     fmt(s["peak_row_offset_median"], 1),
                     fmt(s["frac_dates_peak_within_1"], 2),
                     fmt(s["n_masked_median"], 1),
                     "yes" if s["deck_edge_detectable"] else "no"])
    A(md_table(["series", "n", "n RS/DS1/DS2", "dates with mask",
                "median contrast", "share >= thr", "abs peak offset med",
                "share within 1 px", "n_masked med", "detectable"], rows))
    A("")
    A(f"Detection gate of the origin analysis: band contrast >= "
      f"{mask_definition('contrast_detect'):g} and at least "
      f"{mask_definition('min_n_masked')} deck-edge pixels. The `gamma^2` "
      "channel is the floor-corrected deck-band coherence "
      f"(`{g2}`); the raw estimator and the house/peak variant stay in the "
      "table as audit columns (`bautzen_ocv_channels.py`).")
    A("")

    # ---- 3. result at a glance -------------------------------------------
    A("## 3. Result at a glance")
    A("")
    A("**Bridge/control discrimination of the band contrast** — ratio of the "
      "bridge and control medians, paired per date, percentile bootstrap 95% CI:")
    A("")
    A(md_table(["orbit direction", *core.STATES], [
        [geom] + [f"{fmt(disc[geom][c]['ratio'], 2)} "
                  f"[{fmt(disc[geom][c]['ratio_ci95'][0], 2)}, "
                  f"{fmt(disc[geom][c]['ratio_ci95'][1], 2)}] "
                  f"(n={disc[geom][c]['n_pairs']})" for c in core.STATES]
        for geom in ("ASC", "DESC")]))
    A("")
    A("**Step at the two cut dates** — median after − median before, two-sided "
      "permutation p, and the number of shifted-cut placebos (`k/n`) whose "
      "|delta| reaches the observed one:")
    A("")
    rows = []
    for lab in order:
        row = [series_label(lab)]
        for cut in ("DS1", "DS2"):
            s_ = steps[lab][cut]
            pl = placebo[lab][cut]
            k = sum(1 for p in pl["placebo"]
                    if abs(p["delta"]) >= abs(s_["delta"]))
            row.append(f"{signed(s_['delta'], 2)}  p={fmt(s_['p_perm'], 3)}  "
                       f"{k}/{pl['n_placebo']}")
        rows.append(row)
    A(md_table(["band contrast step", "DS1 (12 May 2025)", "DS2 (29 Sep 2025)"],
               rows))
    A("")
    rows = []
    for lab in order:
        row = [series_label(lab)]
        for cut in ("DS1", "DS2"):
            parts = []
            for var in variants:
                b = g2_steps[lab][var][cut]
                pl = g2_placebo[lab][var][cut]
                k = sum(1 for p in pl["placebo"]
                        if abs(p["delta"]) >= abs(b["delta"]))
                parts.append(f"{G2_SHORT[var]} {signed(b['delta'], 3)} "
                             f"(p={fmt(b['p_perm'], 3)}, {k}/{pl['n_placebo']})")
            row.append("; ".join(parts))
        rows.append(row)
    A(md_table(["gamma^2 step", "DS1", "DS2"], rows))
    A("")
    A("**Season proxy** — gamma^2 against the air temperature at the overpass "
      "instant (OLS, one row per series; `r` and the two-sided p of the slope):")
    A("")
    A(md_table(["series", "n", "r(T)", "p(T)", "r(RH)", "p(RH)"], [
        [series_label(lab),
         environ[lab]["n_dates_with_env"],
         signed(environ[lab]["vars"]["T"]["fit"]["r"], 2),
         fmt(environ[lab]["vars"]["T"]["fit"]["p"], 3),
         signed(environ[lab]["vars"]["RH"]["fit"]["r"], 2),
         fmt(environ[lab]["vars"]["RH"]["fit"]["p"], 3)]
        for lab in order]))
    A("")

    # ---- 4. figures -------------------------------------------------------
    A("## 4. Figures")
    A("")
    A("### Figure A — the six OCV components per series and state")
    A("")
    A("One dot per overpass (state-coloured), the bar is the median of the "
      "series and state. The bridge series carries a fixed geometric edge in "
      "every state; the control series is the same rectangle shifted 500 m east "
      "and only serves as the counter-check.")
    A("")
    A("Median `gamma^2` (deck band, floor-corrected) per series and state:")
    A("")
    A(md_table(["series", *[f"{c} (n={series[order[0]]['n_by_state'][c]})"
                            for c in core.STATES]], [
        [series_label(lab)] + [
            fmt((res["fig_a"]["descriptives"]["gamma2"][lab][c] or {}).get("median"), 3)
            for c in core.STATES] for lab in order]))
    A("")
    A("The descriptives of all six channels (n, median, mean, standard deviation, "
      "the 5/25/75/95 percentiles) are in `fig_a.descriptives` of the result "
      "JSON.")
    A("")
    A("### Figure B — the deck-edge mask layer")
    A("")
    A("Panels: (a) band contrast per series and state against the detection "
      "gate, (b) share of dates with a mask and its Wilson 95% CI, (c) the "
      "number of deck-edge pixels, (d) the in-series gamma^2 effect size "
      "(Cliff's delta with the RS state as reference).")
    A("")
    A(md_table(["series"] + [f"{c}: contrast med" for c in core.STATES]
               + ["share with mask", "detectable"], [
        [series_label(lab)]
        + [fmt(series[lab]["metrics"][c]["contrast_median"], 2)
           for c in core.STATES]
        + [f"{mask_series[lab]['summary']['n_dates_with_mask']}/"
           f"{mask_series[lab]['summary']['n_dates']}",
           "yes" if mask_series[lab]["summary"]["deck_edge_detectable"] else "no"]
        for lab in order]))
    A("")
    A("### Figure C — intervention response against shifted cuts")
    A("")
    A("The observed cut is special only if it is the most extreme cut of the "
      "series. `k/n` counts the shifted cuts whose |delta| reaches the observed "
      "one, so `k = n` means 'no placebo cut is as large' and `real is most "
      "extreme` says whether the observed cut leads the ranking.")
    A("")
    rows = []
    for lab in order:
        for cut in ("DS1", "DS2"):
            pl = placebo[lab][cut]
            s_ = steps[lab][cut]
            k = sum(1 for p in pl["placebo"]
                    if abs(p["delta"]) >= abs(s_["delta"]))
            rows.append([series_label(lab), cut, signed(s_["delta"], 3),
                         fmt(s_["p_perm"], 3), f"{k}/{pl['n_placebo']}",
                         fmt(pl["max_abs_placebo_delta"], 2),
                         "yes" if pl["real_is_most_extreme"] else "no"])
    A(md_table(["series", "cut", "band contrast delta", "p",
                "placebos >= observed", "max abs placebo delta",
                "real most extreme"], rows))
    A("")
    A("### Figure D — bridge minus control, date by date")
    A("")
    A("The paired difference removes everything the two rectangles have in "
      "common (orbit, date, weather). A bridge-specific intervention effect "
      "must show up here, in the difference of the two bridge rectangles, and "
      "*not* in the single-series step tests of the two control rectangles.")
    A("")
    dvc_pick = {"band_contrast": dvc, "gamma2": g2_dvc}
    ch_name = {"band_contrast": "band contrast", "gamma2": g2}
    A(md_table(["channel", "geom", *core.STATES], [
        [ch_name[ch], geom] + [
            (f"{signed(dvc_pick[ch][geom]['states'][c]['median'], 3)}"
             f" (n={dvc_pick[ch][geom]['states'][c]['n']})"
             if dvc_pick[ch][geom]["states"][c] else "n/a")
            for c in core.STATES]
        for ch in ("band_contrast", "gamma2") for geom in ("ASC", "DESC")]))
    A("")
    A(md_table(["channel", "geom", "n pairs", "step DS1 (p)", "step DS2 (p)",
                "placebo at DS1"], [
        [ch_name[ch], geom, dvc_pick[ch][geom]["n_pairs"],
         f"{signed(dvc_pick[ch][geom]['step']['DS1']['delta'], 3)} "
         f"(p={fmt(dvc_pick[ch][geom]['step']['DS1']['p_perm'], 3)})",
         f"{signed(dvc_pick[ch][geom]['step']['DS2']['delta'], 3)} "
         f"(p={fmt(dvc_pick[ch][geom]['step']['DS2']['p_perm'], 3)})",
         "most extreme"
         if dvc_pick[ch][geom]["placebo"]["DS1"]["real_is_most_extreme"]
         else "not extreme"]
        for ch in ("band_contrast", "gamma2") for geom in ("ASC", "DESC")]))
    A("")
    A("### Figure E — seasonality and environment")
    A("")
    A("Monthly medians (top) and gamma^2 against temperature / humidity at the "
      "overpass instant (bottom). Because the three states are consecutive time "
      "intervals, any state difference is *also* a seasonal difference; this "
      "figure makes that confound visible instead of hiding it.")
    A("")
    for label, blk, field in (("band contrast", seasons, "contrast_median"),
                              ("gamma^2", g2_seasons,
                               f"{core.GAMMA2_COLUMN}_median")):
        months = sorted({m for lab in order for m in blk[lab]})
        A(f"Monthly median {label} (empty where the series has no overpass):")
        A("")
        A(md_table(["series", *months], [
            [series_label(lab)] + [fmt(blk[lab][m][field], 2)
                                   if m in blk[lab] else "–" for m in months]
            for lab in order]))
        A("")

    # ---- 5. verification --------------------------------------------------
    A("## 5. Verification against the committed Bautzen reference JSONs")
    A("")
    A("Every block of the two committed reference files is recomputed from the "
      "channel table and compared leaf by leaf:")
    A("")
    A(f"- `{os.path.basename(core.REFERENCE_FILES['deck_edge_mask'])}` — the "
      "`config` block (10 constants of the mask module) and the four series "
      "blocks (`file`, `per_date` with its exact per-date key set, `summary`, "
      "`union50`, `union30`).")
    A(f"- `{os.path.basename(core.REFERENCE_FILES['deck_edge_state'])}` — the "
      "`config` block (22 keys incl. the two German doc strings), the four "
      "series blocks (`per_date`, `metrics` with the `gamma2` sub-block, "
      "`jaccard`, `n_by_state`, `dates_by_state`) and the ten cross-series "
      "blocks (`discrimination`, `step_tests`, `placebo`, `delta_vs_ctrl`, "
      "`seasonality`, `gamma2_step_tests`, `gamma2_placebo`, "
      "`gamma2_delta_vs_ctrl`, `gamma2_seasonality`, `gamma2_environment`).")
    A("")
    A(md_table(["pin result", "value"], [
        ["checks", pins["n_checks"]],
        ["exact (byte for byte)", pins["n_exact"]],
        ["within tolerance", f"{pins['n_within_tolerance']} "
                             f"(<= {pins['float_tolerance']:g})"],
        ["failures", pins["n_failed"]],
        ["max |deviation|", pins["max_abs_deviation"]],
        ["doc strings translated", pins["label_translations"]],
    ]))
    A("")
    A("**Comparison rules.** Deterministic numbers (medians, effect sizes, SDS, "
      "counters, mask sets, unions, gamma^2 values) must match exactly. "
      "Permutation and scipy p-values — "
      f"`{'`, `'.join(pins['tolerance_keys'])}` — are compared with a tolerance "
      f"of {pins['float_tolerance']:g}, because they depend on the library's RNG "
      "and normalisation internals. Committed German doc strings are compared "
      "through `core.LABEL_EN` (an unmapped German string mismatches and becomes "
      "a failure). Committed JSON lists and recomputed tuples are normalised "
      "before the comparison.")
    A("")
    A(f"**Meta cross-checks.** {len(d['cross_checks'])} independent checks "
      "against `bautzen_ocv_channels_meta.json` "
      f"({d['n_cross_checks_failed']} failed): the file names and geometries of "
      "the four series, the per-date row counts, the mask summary, both "
      "persistence unions, the mask definition constants, the sha256 of the four "
      "rect files and the reference-source names.")
    A("")
    extra = pins["recomputed_keys_not_in_reference"]
    A("**Keys beyond the reference.** The recomputed blocks carry "
      f"{len(extra)} additional keys (the figure payloads and the per-variant "
      "gamma^2 gates); they are listed in "
      "`pins.recomputed_keys_not_in_reference` and are not failures — the "
      "reference files do not contain these values, they are added here.")
    A("")
    if pins["failures"]:
        A("**Failures:**")
        A("")
        A(md_table(["path", "recomputed", "committed", "|dev|"], [
            [f["path"], f["recomputed"], f["committed"], fmt(f["abs_deviation"])]
            for f in pins["failures"]]))
    else:
        A("**Failures: none.** Every pinned number of the two committed "
          "reference files is reproduced from the channel table alone.")
    A("")

    # ---- 6. deviations ----------------------------------------------------
    A("## 6. Documented deviations from the LUMO figure package")
    A("")
    A("The LUMO package (`code/lumo/fig_lumo_ocv_paper.py`) is the template of "
      "this script: same layout, same pin philosophy, same figure style. Four "
      "things differ by design, all of them documented in the header of this "
      "script and in `bautzen_ocv_masks.py`:")
    A("")
    A("1. **The mask layer is site-specific.** The LUMO mask rules are echo "
      "based (the deck echo decides the mask). At this site that rule collapses "
      "— the deck-edge echo is not a separated peak — so the mask is "
      "*structurally* attributed: a deck band around the known deck-edge "
      "geometry (`bautzen_ocv_masks.py`). The house/peak variant stays in the "
      "table as a signal-based audit channel and is explicitly *not* a "
      "deck-edge statement.")
    A("2. **Four series instead of one state series.** Here the four rect series "
      "(bridge/control × ASC/DESC) are the counter-check, so every effect is "
      "read per series and the paired bridge-minus-control difference is the "
      "headline comparison (figure D).")
    A("3. **State and season are confounded.** RS, DS1 and DS2 are consecutive "
      "time intervals of a single deployment, so state and season cannot be "
      "separated by design. The script therefore reports the placebo sweep "
      "(figure C) and the environment proxy (figure E) instead of treating the "
      "state contrast as causal.")
    A("4. **The committed references are German, this report is English.** Doc "
      "strings are compared through `core.LABEL_EN`; the number of translated "
      f"strings is reported in the pin summary ({pins['label_translations']} in "
      "this run). Series keys stay German — they are the committed file names.")
    A("")
    A(f"`--quick` copies the placebo sweeps of the reference (`placebo`, "
      "`gamma2_placebo`, `delta_vs_ctrl[*].placebo`, "
      "`gamma2_delta_vs_ctrl[*].placebo`) instead of recomputing them; the "
      "copied blocks are recorded in `meta.copied_blocks`. This report was "
      f"produced in **{'--quick' if quick else 'full'}** mode.")
    A("")

    # ---- 7. reproduce -----------------------------------------------------
    A("## 7. Reproduce")
    A("")
    A("```sh")
    A("cd <repo root>")
    A("python3 code/bautzen/fig_bautzen_ocv_paper.py            # full run")
    A("python3 code/bautzen/fig_bautzen_ocv_paper.py --quick    # placebo copied")
    A("bash figures/bautzen/fig_bautzen_ocv_paper.sh          # the same, via the shell wrapper")
    A("```")
    A("")
    A("The figures need numpy, matplotlib (600 dpi PNG) and scipy. This report "
      f"was generated with python {platform.python_version()}, "
      f"numpy {np.__version__}, matplotlib {matplotlib.__version__}, scipy "
      f"{getattr(scipy, '__version__', None) if scipy else None}. All random "
      "draws use fixed seeds and no output carries a timestamp or a run time, "
      "so repeated runs produce identical JSON, Markdown and PNG files.")
    A("")

    # ---- 8. caveats -------------------------------------------------------
    A("## 8. Caveats")
    A("")
    A(f"- **Small samples.** {d['n_by_state']['RS']} overpasses are in RS in "
      f"total, {d['n_by_state']['DS2']} in DS2; the per-series, per-state cells "
      "go down to "
      f"{min(series[lab]['n_by_state'][c] for lab in order for c in core.STATES)} "
      "overpasses. Every interval in this report (Wilson, bootstrap) is wide for "
      "that reason, and the permutation p-values are quantised by the number of "
      "permutations and — in figure C — by the number of placebo cuts.")
    A("- **The vector is not independent.** `gamma^2` and the mask dimensions "
      "describe the same band (`gamma^2` is the coherence *inside* the mask); "
      "the DIMS are complementary, not orthogonal, evidence.")
    A("- **Thresholds are inherited.** The detectability gate (band contrast, "
      "deck-edge pixel count) and the placebo offsets (+/-2..12 weeks) are taken "
      "over from the origin analysis, not re-tuned here.")
    A("- **State assignment.** States are assigned by overpass time "
      "(`core.state_of`), so a date inside the DS1 window is labelled DS1 even "
      "if the intervention left no trace in the image. Figures B and C test "
      "exactly that, they do not assume it.")
    A("- **Season is not removed.** The temperature/humidity regression of "
      "figure E is a probe, not a correction: the states are months apart and "
      "the covariate is collinear with them by construction.")
    A("")
    A("---")
    A("")
    A("Generated by `code/bautzen/fig_bautzen_ocv_paper.py` "
      f"({'--quick' if quick else 'full'} mode, "
      f"{pins['n_checks']} pins, {pins['n_failed']} failures).")
    return "\n".join(L) + "\n"

# ---------------------------------------------------------------------------
# Command line
# ---------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Recompute the five Bautzen deck-edge OCV figures (A–E).")
    ap.add_argument("--quick", action="store_true",
                    help="copy the placebo sweeps from the committed reference "
                         "instead of recomputing them")
    ap.add_argument("--figdir", default=FIGDIR,
                    help="output directory of the PNG figures")
    ap.add_argument("--json", default=OUT_JSON, help="path of the result JSON")
    ap.add_argument("--md", default=OUT_MD, help="path of the Markdown report")
    ap.add_argument("--no-figures", action="store_true",
                    help="only recompute and pin (no matplotlib output)")
    args = ap.parse_args(argv)

    t0 = time.time()
    print("Bautzen deck-edge OCV figures A–E — "
          f"{'quick' if args.quick else 'full'} mode ...")
    res, mask_ref, state_ref = compute(quick=args.quick)
    print(f"  compute: {time.time() - t0:.1f}s "
          f"({res['data']['n_rows']} overpasses)")

    pins = pin_all(res, mask_ref, state_ref)
    print(f"  pins: {pins['n_checks']} checks, {pins['n_exact']} exact, "
          f"{pins['n_within_tolerance']} within tolerance, "
          f"{pins['n_failed']} failures, "
          f"{pins['label_translations']} doc strings translated, "
          f"{pins['cross_checks']['n']} meta cross-checks "
          f"({pins['cross_checks']['n_failed']} failed)")
    for f in pins["failures"]:
        print(f"    FAIL {f['path']}: {f['recomputed']} != {f['committed']}")
    for c in pins["cross_checks"]["failed"]:
        print(f"    FAIL {c['check']}: {c['recomputed']} != {c['committed']}")

    if not args.no_figures:
        os.makedirs(args.figdir, exist_ok=True)
        print("  writing figures:")
        for fn in (fig_a, fig_b, fig_c, fig_d, fig_e):
            fn(res, args.figdir)

    os.makedirs(os.path.dirname(os.path.abspath(args.json)), exist_ok=True)
    with open(args.json, "w", encoding="utf-8") as fh:
        json.dump(build_json(res, pins, args.quick), fh, indent=1,
                  sort_keys=True, default=json_default)
        fh.write("\n")
    print(f"  {os.path.relpath(args.json)}")

    os.makedirs(os.path.dirname(os.path.abspath(args.md)), exist_ok=True)
    with open(args.md, "w", encoding="utf-8") as fh:
        fh.write(report(res, pins, args.quick))
    print(f"  {os.path.relpath(args.md)}")

    print(f"done in {time.time() - t0:.1f}s, {pins['n_failed']} pin failures")
    return 0 if pins["ok"] else 1


if __name__ == "__main__":
    raise SystemExit(main())




