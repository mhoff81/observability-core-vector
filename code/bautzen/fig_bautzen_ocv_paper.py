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
import json
import math
import os
import platform
import sys
import time

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
    "gamma2": "deck-band coherence gamma^2 (floor-corrected)",
    "A": "mask area (deck-edge pixels)", "D": "mask density",
    "F": "mask fragmentation", "S": "centroid<->deck-local peak shift",
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


def panel(ax, letter, text):
    ax.set_title(f"({letter}) {text}", fontsize=8.0, loc="left")


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
# two gamma^2 strings are produced by ``masks.state_metrics`` itself, so their
# English constants are used here directly.
DOC_DE = {
    "mask_source": "analyze_bautzen_deck_edge_mask.py (unverändert übernommen)",
    "formula": "|sum(z)|^2 / (N * sum(|z|^2)) ueber die Maskenpixel",
    "floor_corrected": masks.GAMMA2_FLOOR_CORRECTION,
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

    Dates without a deck-edge pixel carry no centroid, and ``D``/``S`` are then
    null -- the committed rows are not equally wide on every date.
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
    out.update({
        "band_peak_r": int(row["band_peak_r"]),
        "band_peak_c": int(row["band_peak_c"]),
        "A": int(row["A"]), "D": row["D"], "F": int(row["F"]), "S": row["S"],
        "row_span": int(row["row_span"]), "col_span": int(row["col_span"]),
        "bbox_area": int(row["bbox_area"]),
    })
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
        "peak_row_offset_median": float(np.median(offsets)) if n else None,
        "frac_dates_peak_within_1": (float((offsets <= 1).mean()) if n else None),
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


def gamma2_variant_blocks(rows):
    """The three gamma^2 variants of gamma2_step_tests / gamma2_placebo."""
    steps, placebo = {}, {}
    for key, col, gate, label in G2_VARIANTS:
        keys, vals, sel = _keys_vals(rows, col, gate)
        steps[key] = {"min_n_filter": label, "n_values": len(sel),
                      "DS1": st.step_test(keys, vals, core.CUT_DS1),
                      "DS2": st.step_test(keys, vals, core.CUT_DS2)}
        placebo[key] = {"DS1": st.placebo_sweep(keys, vals, core.CUT_DS1),
                        "DS2": st.placebo_sweep(keys, vals, core.CUT_DS2)}
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


def delta_vs_ctrl(by, col, with_feature=False):
    """Bridge minus control, date-paired, in the committed block shape.

    ``col == 'band_contrast'`` keeps every date the two series have in common
    (band contrast is defined everywhere); a gamma^2 difference needs a mask on
    both sides, so those dates are dropped.
    """
    out = {}
    for geom, ser in paired_diff_series(by, col).items():
        dates, vals = ser["dates"], ser["diffs"]
        diffs = dict(zip(dates, vals))
        keys = [core.acq_dt(d, geom) for d in dates]
        states = {}
        for s in core.STATES:
            v = [diffs[d] for d in dates if br[d]["state"] == s]
            states[s] = ({"n": len(v), "median": float(np.median(v)),
                          "iqr": [float(np.percentile(v, 25)),
                                  float(np.percentile(v, 75))]}
                         if v else None)
        blk = {"n_pairs": len(dates), "states": states,
               "step": {"DS1": st.step_test(keys, vals, core.CUT_DS1),
                        "DS2": st.step_test(keys, vals, core.CUT_DS2)},
               "placebo": {"DS1": st.placebo_sweep(keys, vals, core.CUT_DS1),
                           "DS2": st.placebo_sweep(keys, vals, core.CUT_DS2)}}
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
        if e is not None:
            n_env += 1
        if r[core.GAMMA2_COLUMN] is None:
            continue
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
    t0 = time.time()
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
        placebo[lab] = placebo_block(rs, "band_contrast")
        seasons[lab] = season_block(rs, "band_contrast", "contrast")
        g2_seasons[lab] = season_block(rs, core.GAMMA2_COLUMN, core.GAMMA2_COLUMN)
        g2_steps[lab], g2_placebo[lab] = gamma2_variant_blocks(rs)
        environ[lab] = environment_block(rs, lab, env)
    dvc = delta_vs_ctrl(by, "band_contrast")
    g2_dvc = delta_vs_ctrl(by, core.GAMMA2_COLUMN, with_feature=True)

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
    for kind in ("deck_edge_mask", "deck_edge_state"):
        check(f"reference_source.{kind}", os.path.basename(core.REFERENCE_FILES[kind]),
              os.path.basename(meta["reference_source"][kind]))
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
        "runtime_seconds": None,      # filled in by main() (no timestamp)
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": getattr(scipy, "__version__", None),
        "channels": CHANNELS, "chance": CHANCE, "states": list(core.STATES),
        "state_windows": {s: list(core.STATE_WINDOWS[s]) for s in core.STATES},
        "data": data, "mask": mask, "state": state,
        "fig_a": fig_a, "fig_b": fig_b, "fig_c": fig_c, "fig_d": fig_d,
        "fig_e": fig_e,
    }
    result["runtime_seconds"] = round(time.time() - t0, 2)
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








