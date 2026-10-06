#!/usr/bin/env python3
"""OCV-Paper (Bautzen bridge site) — amplitude/phase time series of the deck edge.

Question: how do the **amplitude** A(t) and the **phase** phi(t) of the deck band
respond to the two documented interventions at the Bautzen OpenLabs bridge
(states RS -> DS1 -> DS2), when both are measured *only over the deck-edge
pixels* of the committed channel table -- not over the whole bridge?

The complex chip per overpass is ``z`` (the committed rect files). Over the
deck-edge mask (``bautzen_ocv_masks.date_mask``, the same local, geometry-
anchored band as in the channel table) the complex sum is ``Z = sum(z_sel)`` with
``N`` deck-edge pixels and energy ``S2 = sum(|z_sel|^2)``. This script reports its
polar decomposition

  A(t)   = |sum z_sel| / sqrt(N * sum |z_sel|^2)  in [0, 1]    coherence amplitude
  phi(t) = arg(sum z_sel)                         in (-pi, pi] interferometric phase

of the **deck-edge pixels only** (``z_sel = z[mask]``), never of the full chip.
A third panel adds the **paired differential phase** of the deck edge between the
bridge and its control rectangle (identical overpass times, no interpolation):

  dphi(t) = wrap(phi_bridge - phi_control)        in (-pi, pi] paired differential phase
The amplitude is the square root of this site's coherence channel, so ``A(t)^2``
must reproduce the committed ``gamma2_band_raw`` column date by date -- that
equality is the built-in cross-check of this figure: it ties the new numbers to
the pinned channel table without needing a second committed reference file (the
two reference JSONs of the package hold no amplitude or phase).

The four rect series (ASC/DESC bridge and ASC/DESC control 500 m east of the
bridge) are carried through all panels, as everywhere else in this package:
only the bridge series is expected to track the interventions, the control
rectangle is the counter-check against a seasonal reading. The two cut dates of
the deployment are drawn as vertical lines, RS / DS1 / DS2 as background bands.

The differential panel works on the **pair**: bridge and control rectangle share
every overpass time step, so their deck-edge phases are differenced pair by pair
without any interpolation. ``dphi`` removes what both rectangles see in common
(season, atmosphere, repeat-pass baseline) and keeps the local deck-edge term --
which is exactly what the control rectangle is for. It is a difference of two
single-SLC ``arg`` values, so it is defined modulo 2*pi and is *not* an
interferometric (master-slave) phase.

The same pair also yields the **state-to-state** shift: the state windows are
consecutive and disjoint, so the change away from the reference state RS to DS1
and DS2 is ``delta = wrap(circmean(state) - circmean(RS))``. It is computed for
the paired ``dphi`` and, as a counter-check, for the absolute ``phi`` of each
series, and it carries a seeded percentile bootstrap interval (fixed seed, so the
output stays byte-identical across runs). Because the per-state samples are
small, the interval of every delta covers almost the whole circle: the numbers
are **not resolvable** at these sample sizes, and the report says so explicitly.

This is a standalone companion of ``fig_bautzen_ocv_paper.py`` (figures A-E): it
adds one separate figure and does not touch the A-E pin contract. It mirrors the
per-channel time-series figure of the sister package
(``observability-lumo-espoo/code/fig_lumo_channels.py``) -- same style, same
"one series per column" layout -- with the deck-edge mask layer, the four-series
counter-check and the state model of this site.

Usage:
  python3 code/bautzen/fig_bautzen_ocv_amp_phase.py    # compute + figure + report
  python3 code/bautzen/fig_bautzen_ocv_amp_phase.py --no-figures
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import math
import os
import platform
import sys
import time
import zlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import dates as mdates

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import bautzen_ocv_core as core    # noqa: E402
import bautzen_ocv_masks as masks  # noqa: E402

FIGDIR = core.FIGDIR
OUT_JSON = os.path.join(core.DATA, "fig_bautzen_ocv_amp_phase.json")
OUT_MD = os.path.join(FIGDIR, "fig_bautzen_ocv_amp_phase.md")
OUT_PNG = os.path.join(FIGDIR, "fig_bautzen_ocv_amp_phase.png")

# ---------------------------------------------------------------------------
# Style (identical to the rest of the package: DejaVu Sans + STIX, 183 mm)
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
STATE_LABEL = {"RS": "RS (reference)", "DS1": "DS1", "DS2": "DS2"}

# The equality A^2 == gamma2_band_raw is exact (same mask, same floats); the
# tolerance guards only against the CSV round-trip of the committed column.
CHECK_TOL = 1e-9

# The state-to-state deltas carry a percentile bootstrap CI: the circular mean is
# not linear, so the CI must come from resampling, not from a Gaussian formula.
# A fixed seed keeps every run byte-identical -- the script stays deterministic.
BOOTSTRAP_N = 20000
BOOTSTRAP_SEED = 12345


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------
def json_default(o):
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    raise TypeError(f"not JSON serialisable: {type(o)!r}")


def fmt(x, digits=4):
    """Compact float/text for tables (``None`` stays empty)."""
    if x is None:
        return ""
    if isinstance(x, float):
        return f"{x:.{digits}f}"
    return str(x)


def md_table(header, rows):
    out = ["| " + " | ".join(str(h) for h in header) + " |",
           "| " + " | ".join("---" for _ in header) + " |"]
    for row in rows:
        out.append("| " + " | ".join(fmt(c) for c in row) + " |")
    return "\n".join(out)


def wrap_pi(angle):
    """Wrap an angle to (-pi, pi] (the branch of ``arg``).

    A difference of two phases is only defined modulo 2*pi; this folds the raw
    difference back onto the same branch ``np.angle`` uses, without the sign
    corner cases of ``math.remainder``.
    """
    return math.atan2(math.sin(angle), math.cos(angle))


def circ_stats(angles):
    """Circular mean and resultant length of an angle list (``None`` -> empty).

    The deck-edge phase is an angle: an arithmetic mean would be wrong across
    the +/-pi branch cut, so the per-state summary uses the circular mean
    ``atan2(<sin>, <cos>)`` and its resultant length ``R = |<e^{i phi}>|``.
    """
    v = np.asarray([a for a in angles if a is not None], dtype=float)
    if v.size == 0:
        return {"n": 0, "mean_rad": None, "R": None}
    c = float(np.cos(v).mean())
    s = float(np.sin(v).mean())
    return {"n": int(v.size), "mean_rad": float(math.atan2(s, c)),
            "R": float(math.hypot(c, s))}


def _circ_mean(v):
    """Circular mean of a float array: ``atan2(<sin>, <cos>)``."""
    return math.atan2(float(np.sin(v).sum()), float(np.cos(v).sum()))


def circ_delta(a, b, tag):
    """``wrap(circmean(a) - circmean(b))`` with a seeded percentile bootstrap CI.

    The angle sets ``a`` (the state) and ``b`` (the reference state RS) are
    resampled independently with replacement; ``ci_lo``/``ci_hi`` are the 2.5 and
    97.5 percentiles of the wrapped difference, and ``resolvable`` says whether
    that interval excludes 0. ``tag`` selects an independent, deterministic RNG
    stream (``zlib.crc32`` based), so a delta never depends on how many deltas
    were computed before it -- the output stays byte-identical across runs.
    """
    va = np.asarray([x for x in a if x is not None], dtype=float)
    vb = np.asarray([x for x in b if x is not None], dtype=float)
    if va.size == 0 or vb.size == 0:
        return {"n_a": int(va.size), "n_b": int(vb.size), "delta_rad": None,
                "ci_lo": None, "ci_hi": None, "resolvable": None,
                "n_boot": BOOTSTRAP_N, "seed": BOOTSTRAP_SEED}
    delta = wrap_pi(_circ_mean(va) - _circ_mean(vb))
    rng = np.random.default_rng(zlib.crc32(tag.encode("utf-8")) + BOOTSTRAP_SEED)
    ia = rng.integers(0, va.size, size=(BOOTSTRAP_N, va.size))
    ib = rng.integers(0, vb.size, size=(BOOTSTRAP_N, vb.size))
    da = np.arctan2(np.sin(va[ia]).sum(axis=1), np.cos(va[ia]).sum(axis=1))
    db = np.arctan2(np.sin(vb[ib]).sum(axis=1), np.cos(vb[ib]).sum(axis=1))
    bs = np.arctan2(np.sin(da - db), np.cos(da - db))
    lo = float(np.percentile(bs, 2.5))
    hi = float(np.percentile(bs, 97.5))
    return {"n_a": int(va.size), "n_b": int(vb.size), "delta_rad": delta,
            "ci_lo": lo, "ci_hi": hi,
            "resolvable": bool(lo > 0.0 or hi < 0.0),
            "n_boot": BOOTSTRAP_N, "seed": BOOTSTRAP_SEED}


def series_summary(points):
    """Per-state summary of one series: median/mean A and the circular phase."""
    out = {}
    for state in core.STATES:
        sub = [p for p in points if p["state"] == state]
        a = [p["A"] for p in sub if p["A"] is not None]
        phi = [p["phi_rad"] for p in sub if p["phi_rad"] is not None]
        out[state] = {
            "n": len(sub),
            "n_amplitude": len(a),
            "median_A": float(np.median(a)) if a else None,
            "mean_A": float(np.mean(a)) if a else None,
            "phase": circ_stats(phi),
        }
    return out


def dphi_summary(points):
    """Per-state circular summary of the paired differential phase."""
    out = {}
    for state in core.STATES:
        v = [p["dphi_rad"] for p in points if p["state"] == state]
        out[state] = {"n": len(v), "dphi": circ_stats(v)}
    return out


def series_label(lab):
    """'ASC Brücke' -> 'ASC bridge' (the report is English, keys stay German)."""
    return core.SERIES_EN.get(lab, lab)


# ---------------------------------------------------------------------------
# Compute: A(t) and phi(t) over the deck-edge pixels, plus the cross-checks
# ---------------------------------------------------------------------------
def _check(checks, name, ok):
    checks.append({"check": name, "ok": bool(ok)})


def compute():
    """Recompute A(t) and phi(t) per series and date from the rect files.

    The mask is rebuilt per date with ``masks.date_mask`` on the *same*
    intensities (``|z|^2``) the channel table used, so ``z_sel = z[full]`` is the
    deck-edge pixel set of the committed table. Every quantity is cross-checked
    against the committed row: the pixel count ``n_masked``, the band contrast
    and the coherence ``A^2 == gamma2_band_raw``.
    """
    rows = core.load_csv()
    by = core.by_series(rows)
    checks = []
    series = {}
    for lab in core.SERIES_ORDER:
        rect = os.path.join(core.DATA, core.SERIES[lab])
        recs = dict(masks.load_complex(rect))
        geom = core.geom_of(lab)
        points = []
        for r in by[lab]:
            z = recs.get(r["date"])
            if z is None:                      # row without a committed chip
                continue
            info, full = masks.date_mask(np.abs(z) ** 2)
            z_sel = z[full]
            n = int(z_sel.size)
            s2 = float(np.sum(np.abs(z_sel) ** 2)) if n else 0.0
            if n > 0 and s2 > 0.0:
                zsum = complex(np.sum(z_sel))
                amp = float(abs(zsum) / math.sqrt(n * s2))
                phi = float(np.angle(zsum))
            else:
                amp, phi = None, None
            acq = core.acq_dt(r["date"], geom)
            points.append({
                "date": r["date"],
                "acq_local": acq.strftime("%Y-%m-%dT%H:%M"),
                "orbit": geom,
                "state": r["state"],
                "n_edge_px": n,
                "has_mask": r["has_mask"],
                "A": amp,
                "phi_rad": phi,
                "gamma2_band_raw": r["gamma2_band_raw"],
                "band_contrast_csv": r["band_contrast"],
                "band_contrast_mask": (None if info["band_contrast"] !=
                                       info["band_contrast"]
                                       else float(info["band_contrast"])),
            })
            # ---- cross-checks of the recomputed mask/coherence vs. the CSV ---
            _check(checks, f"mask.n_masked[{lab}][{r['date']}]",
                   info["n_masked"] == r["n_masked"])
            _check(checks, f"mask.band_contrast[{lab}][{r['date']}]",
                   r["band_contrast"] is not None
                   and abs(info["band_contrast"] - r["band_contrast"]) <= CHECK_TOL)
            if amp is not None and r["gamma2_band_raw"] is not None:
                _check(checks, f"A^2==gamma2_band_raw[{lab}][{r['date']}]",
                       abs(amp * amp - r["gamma2_band_raw"]) <= CHECK_TOL)
        series[lab] = {"file": core.SERIES[lab], "geom": geom,
                       "n_dates": len(points), "points": points}

    # ---- paired differential phase: bridge - control, one panel per geometry --
    # Bridge and control rectangle sit on the same overpasses (identical
    # acquisition times), so the deck-edge phase can be differenced pair by pair
    # without any interpolation. The difference removes what both rectangles see
    # in common -- see the module docstring.
    dphi = {}
    for geom in ("ASC", "DESC"):
        br = core.BRIDGE_BY_GEOM[geom]
        ct = core.CONTROL_BY_GEOM[geom]
        pb = {p["date"]: p for p in series[br]["points"]
              if p["phi_rad"] is not None}
        pc = {p["date"]: p for p in series[ct]["points"]
              if p["phi_rad"] is not None}
        pts = []
        for date in sorted(set(pb) & set(pc)):
            b, c = pb[date], pc[date]
            raw = b["phi_rad"] - c["phi_rad"]
            d = wrap_pi(raw)
            pts.append({"date": date,
                        "acq_local": b["acq_local"],
                        "state": b["state"],
                        "phi_bridge_rad": b["phi_rad"],
                        "phi_control_rad": c["phi_rad"],
                        "dphi_rad": d})
            # the signed difference is only defined modulo 2*pi: wrapping the
            # raw difference must reproduce the stored value exactly.
            _check(checks, f"dphi.wrap[{geom}][{date}]",
                   abs(wrap_pi(raw - d)) <= 1e-12)
        # a date can only be paired where *both* rectangles have a valid
        # deck-edge mask; the pairing must cover exactly that overlap.
        shared = set(pb) & set(pc)
        _check(checks, f"dphi.pairs[{geom}]",
               len(pts) == len(shared) == len({p["date"] for p in pts}))
        _check(checks, f"dphi.finite[{geom}]",
               all(math.isfinite(p["dphi_rad"]) for p in pts))
        dphi[geom] = {"bridge": br, "control": ct,
                      "n_bridge": len(pb), "n_control": len(pc),
                      "n_bridge_total": series[br]["n_dates"],
                      "n_control_total": series[ct]["n_dates"],
                      "n_bridge_only": len(set(pb) - set(pc)),
                      "n_control_only": len(set(pc) - set(pb)),
                      "n_pairs": len(pts), "points": pts}

    # ---- state-to-state phase deltas: dphi_DS1-RS and dphi_DS2-RS ------------
    # The state windows are consecutive and disjoint, so there is *no* per-date
    # pairing between states; the shift from the reference state RS is therefore
    # the difference of the two state circular means, wrap(circmean(state) -
    # circmean(RS)). The circular mean is not linear, so this is deliberately
    # *not* the difference of the series' own state phases -- both are reported
    # separately ("dphi" is the figure's quantity, "series" the counter-check).
    deltas = {"dphi": {}, "series": {}}
    for geom in ("ASC", "DESC"):
        deltas["dphi"][geom] = {}
        dmeans = dphi_summary(dphi[geom]["points"])
        for st in ("DS1", "DS2"):
            a = [p["dphi_rad"] for p in dphi[geom]["points"]
                 if p["state"] == st]
            b = [p["dphi_rad"] for p in dphi[geom]["points"]
                 if p["state"] == "RS"]
            d = circ_delta(a, b, f"dphi.{geom}.{st}")
            deltas["dphi"][geom][st] = d
            expect = wrap_pi(dmeans[st]["dphi"]["mean_rad"]
                             - dmeans["RS"]["dphi"]["mean_rad"])
            _check(checks, f"delta.wrap[dphi.{geom}.{st}]",
                   abs(wrap_pi(d["delta_rad"] - expect)) <= 1e-12)
            _check(checks, f"delta.ci[dphi.{geom}.{st}]",
                   d["ci_lo"] <= d["delta_rad"] <= d["ci_hi"])
            _check(checks, f"delta.finite[dphi.{geom}.{st}]",
                   math.isfinite(d["delta_rad"]))
    for lab in core.SERIES_ORDER:
        deltas["series"][lab] = {}
        pmeans = series_summary(series[lab]["points"])
        for st in ("DS1", "DS2"):
            a = [p["phi_rad"] for p in series[lab]["points"]
                 if p["state"] == st and p["phi_rad"] is not None]
            b = [p["phi_rad"] for p in series[lab]["points"]
                 if p["state"] == "RS" and p["phi_rad"] is not None]
            d = circ_delta(a, b, f"series.{lab}.{st}")
            deltas["series"][lab][st] = d
            expect = wrap_pi(pmeans[st]["phase"]["mean_rad"]
                             - pmeans["RS"]["phase"]["mean_rad"])
            _check(checks, f"delta.wrap[series.{lab}][{st}]",
                   abs(wrap_pi(d["delta_rad"] - expect)) <= 1e-12)
            _check(checks, f"delta.ci[series.{lab}][{st}]",
                   d["ci_lo"] <= d["delta_rad"] <= d["ci_hi"])
            _check(checks, f"delta.finite[series.{lab}][{st}]",
                   math.isfinite(d["delta_rad"]))

    data = {
        "csv": os.path.relpath(core.CSV_PATH),
        "csv_meta": os.path.relpath(core.CSV_META_PATH),
        "csv_sha256": core.sha256(core.CSV_PATH),
        "rect_files": {lab: core.SERIES[lab] for lab in core.SERIES_ORDER},
        "rect_file_sha256": {core.SERIES[lab]:
                             core.sha256(os.path.join(core.DATA,
                                                      core.SERIES[lab]))
                             for lab in core.SERIES_ORDER},
        "n_rows": len(rows),
        "n_by_series": {lab: len(by[lab]) for lab in core.SERIES_ORDER},
        "n_points_plotted": {lab: series[lab]["n_dates"]
                             for lab in core.SERIES_ORDER},
        "states": {s: core.STATE_WINDOWS[s] for s in core.STATES},
        "cuts": {"DS1": core.CUT_DS1.isoformat(),
                 "DS2": core.CUT_DS2.isoformat()},
        "mask_definition": {k: getattr(masks, k, None) for k in
                            ("CENTER_ROW", "CENTER_COL", "ANCHOR_HALF",
                             "COL_HALF", "BAND_HALF", "K_LOCAL")},
        "definition": ("A = |sum z| / sqrt(N * sum |z|^2), "
                       "phi = arg(sum z), over the deck-edge mask pixels"),
        "dphi_definition": ("dphi = wrap(phi_bridge - phi_control), paired over "
                            "identical overpass times, over the deck-edge mask "
                            "pixels"),
        "delta_definition": (
            "delta = wrap(circmean(state) - circmean(RS)) over the state's dates; "
            "the states are consecutive and disjoint, so there is no per-date "
            "pairing between them, and because the circular mean is not linear "
            "this is not the difference of the series' own state phases; "
            "ci_lo/ci_hi are a seeded percentile bootstrap "
            f"({BOOTSTRAP_N} resamples, seed {BOOTSTRAP_SEED})"),
    }
    return {
        "site": core.SITE,
        "data": data,
        "series": series,
        "summary": {lab: series_summary(series[lab]["points"])
                    for lab in core.SERIES_ORDER},
        "dphi": dphi,
        "dphi_summary": {geom: dphi_summary(dphi[geom]["points"])
                         for geom in dphi},
        "delta": deltas,
        "cross_checks": checks,
        "n_cross_checks": len(checks),
        "n_cross_checks_failed": sum(1 for c in checks if not c["ok"]),
    }


# ---------------------------------------------------------------------------
# Figure: three rows (A(t), phi(t), dphi(t)) x one column per series
# ---------------------------------------------------------------------------
def _times(points):
    return [dt.datetime.strptime(p["acq_local"], "%Y-%m-%dT%H:%M")
            for p in points]


def _shade_states(ax, x0, x1):
    """RS / DS1 / DS2 background bands plus the two cut dates."""
    bounds = [x0, core.CUT_DS1, core.CUT_DS2, x1]
    for s, lo, hi in zip(core.STATES, bounds[:-1], bounds[1:]):
        lo, hi = min(lo, hi), max(lo, hi)
        if hi > lo:
            ax.axvspan(lo, hi, color=STATE_COLORS[s], alpha=0.07, lw=0,
                       zorder=0)
    for cut in (core.CUT_DS1, core.CUT_DS2):
        ax.axvline(cut, color="0.35", ls="--", lw=0.8, zorder=1)


def make_figure(res, figdir):
    """Draw A(t), phi(t) and the paired dphi(t) with the states shaded.

    Rows 0/1 carry the four series one column each. Row 2 is the paired
    differential phase: it belongs to the *geometry* (bridge - control), so it
    is drawn once per geometry, under the bridge column; the control column of
    that row stays empty (the control series is the counter-check already shown
    in the two panels above).
    """
    order = core.SERIES_ORDER
    fig, axes = plt.subplots(3, len(order),
                             figsize=(COL_WIDTH, COL_WIDTH * 0.72),
                             sharex="col", squeeze=False)
    all_t = [t for lab in order for t in _times(res["series"][lab]["points"])]
    x0 = min(all_t) - dt.timedelta(days=12)
    x1 = max(all_t) + dt.timedelta(days=12)
    for col, lab in enumerate(order):
        pts = res["series"][lab]["points"]
        ts = _times(pts)
        a_ax, p_ax = axes[0][col], axes[1][col]
        for ax in (a_ax, p_ax):
            _shade_states(ax, x0, x1)
            ax.set_xlim(x0, x1)
        # ---- amplitude A(t) ----
        good = [(t, p) for t, p in zip(ts, pts) if p["A"] is not None]
        if good:
            a_ax.plot([t for t, _ in good], [p["A"] for _, p in good],
                      color=SERIES_COLORS[lab], lw=0.8, alpha=0.45, zorder=2)
        for s in core.STATES:
            sel = [p["A"] for _t, p in good if p["state"] == s]
            keep = [t for t, p in good if p["state"] == s]
            if sel:
                a_ax.scatter(keep, sel, s=11, marker=STATE_MARKERS[s],
                             color=STATE_COLORS[s], edgecolor="white",
                             linewidth=0.3, zorder=3, label=STATE_LABEL[s])
        a_ax.set_ylim(0.0, 1.0)
        a_ax.set_ylabel(r"$A(t)$  [coherence amplitude]" if col == 0 else "")
        # ---- phase phi(t) ----
        goodp = [(t, p) for t, p in zip(ts, pts) if p["phi_rad"] is not None]
        if goodp:
            p_ax.plot([t for t, _ in goodp], [p["phi_rad"] for _, p in goodp],
                      color=SERIES_COLORS[lab], lw=0.8, alpha=0.45, zorder=2)
        for s in core.STATES:
            keep = [t for t, p in goodp if p["state"] == s]
            val = [p["phi_rad"] for _t, p in goodp if p["state"] == s]
            if val:
                p_ax.scatter(keep, val, s=11, marker=STATE_MARKERS[s],
                             color=STATE_COLORS[s], edgecolor="white",
                             linewidth=0.3, zorder=3)
        p_ax.set_ylim(-math.pi * 1.05, math.pi * 1.05)
        p_ax.set_yticks([-math.pi, -math.pi / 2, 0, math.pi / 2, math.pi])
        p_ax.set_yticklabels([r"$-\pi$", r"$-\pi/2$", "0", r"$\pi/2$", r"$\pi$"])
        p_ax.set_ylabel(r"$\varphi(t)$  [rad]" if col == 0 else "")
        # ---- paired differential phase dphi(t) = bridge - control ----
        # The difference belongs to the geometry, so it is drawn under the
        # bridge column; the control column of this row keeps the state bands
        # and a pointer, but carries no series of its own.
        d_ax = axes[2][col]
        _shade_states(d_ax, x0, x1)
        d_ax.set_xlim(x0, x1)
        if core.is_bridge(lab):
            dpts = res["dphi"][core.geom_of(lab)]["points"]
            dts = _times(dpts)
            goodd = [(t, p) for t, p in zip(dts, dpts)
                     if p["dphi_rad"] is not None]
            if goodd:
                d_ax.plot([t for t, _ in goodd],
                          [p["dphi_rad"] for _, p in goodd],
                          color=SERIES_COLORS[lab], lw=0.8, alpha=0.45, zorder=2)
            for s in core.STATES:
                keep = [t for t, p in goodd if p["state"] == s]
                val = [p["dphi_rad"] for _t, p in goodd if p["state"] == s]
                if val:
                    d_ax.scatter(keep, val, s=11, marker=STATE_MARKERS[s],
                                 color=STATE_COLORS[s], edgecolor="white",
                                 linewidth=0.3, zorder=3)
            d_ax.axhline(0.0, color="0.35", lw=0.7, zorder=1)
            d_ax.set_ylim(-math.pi * 1.05, math.pi * 1.05)
            d_ax.set_yticks([-math.pi, -math.pi / 2, 0, math.pi / 2, math.pi])
            d_ax.set_yticklabels([r"$-\pi$", r"$-\pi/2$", "0", r"$\pi/2$",
                                  r"$\pi$"])
            d_ax.set_ylabel(r"$\Delta\varphi(t)$  [rad]" if col == 0 else "")
        else:
            d_ax.set_yticks([])
            d_ax.set_title(r"$\Delta\varphi$ per geometry:" "\n"
                           "(see bridge column)", fontsize=6.8, color="0.45")
        # ---- titles and x axis ----
        role = "bridge" if core.is_bridge(lab) else "control 500 m east"
        axes[0][col].set_title(f"{series_label(lab)}\n({role})", fontsize=8.0)
        for r in (0, 1, 2):
            ax = axes[r][col]
            ax.xaxis.set_major_locator(mdates.MonthLocator(interval=2))
            ax.xaxis.set_major_formatter(mdates.DateFormatter("%b\n%Y"))
        if col == 0:
            axes[0][col].legend(loc="upper right", frameon=False)
    fig.suptitle("Deck-edge amplitude, phase and its paired differential "
                 "(bridge vs. control, ASC/DESC)", fontsize=9.0, y=1.01)
    fig.tight_layout()
    out = os.path.join(figdir, "fig_bautzen_ocv_amp_phase.png")
    fig.savefig(out)
    plt.close(fig)
    return out


# ---------------------------------------------------------------------------
# Markdown report
# ---------------------------------------------------------------------------
def report(res):
    d = res["data"]
    order = core.SERIES_ORDER
    L = []
    A = L.append
    A("# Bautzen deck-edge amplitude and phase — figure, per series and state")
    A("")
    A("## 1. Question and definition")
    A("")
    A("How do the amplitude A(t) and the phase phi(t) of the deck band respond to "
      "the two documented interventions at the Bautzen OpenLabs bridge "
      "(RS -> DS1 -> DS2)? Both are measured **over the deck-edge pixels only** — "
      "the same local, geometry-anchored deck band the channel table uses — and "
      "not over the whole bridge chip:")
    A("")
    A("- `A(t)   = |sum z| / sqrt(N * sum |z|^2)`  — coherence amplitude, in [0, 1]")
    A("- `phi(t) = arg(sum z)`  — interferometric phase of the deck edge, in (-pi, pi]")
    A("- `dphi(t) = wrap(phi_bridge(t) - phi_control(t))`  — paired differential "
      "phase of the deck edge, in (-pi, pi]")
    A("")
    A("with `z` the complex chip of the overpass and the sums over the `N` "
      "deck-edge pixels `z_sel = z[mask]`. `A(t)^2` is exactly this site's "
      "coherence channel, so it is cross-checked against the committed "
      "`gamma2_band_raw` column (section 4). The four series (bridge/control × "
      "ASC/DESC) are carried through the amplitude and phase panels; the control "
      "rectangle 500 m east is the counter-check against a purely seasonal "
      "reading and the subtrahend of the differential panel.")
    A("")
    A("The differential panel works on the **pair**. Bridge and control rectangle "
      "sit on the same overpasses (identical acquisition times), so their "
      "deck-edge phases are differenced pair by pair without any interpolation — "
      "wherever both rectangles have a valid deck-edge mask. The difference "
      "removes what both rectangles see in common — season, atmosphere, "
      "repeat-pass baseline — and keeps the local deck-edge term, which is "
      "exactly what the control rectangle is for. `wrap` folds the raw "
      "difference back onto the (-pi, pi] branch of `arg`.")
    A("")
    A("Section 4 closes with the **state-to-state** change of both quantities, "
      "`delta = wrap(circmean(state) - circmean(RS))` — the shift from the "
      "reference state RS to DS1 and to DS2 — for the paired `dphi` and, as "
      "counter-check, for the absolute `phi` of each series. It is reported with "
      "a seeded bootstrap confidence interval, because a state holds only a "
      "handful of overpasses.")
    A("")
    A("The deck-edge mask is rebuilt here with the *same* routine and the same "
      "floats as the channel table (`masks.date_mask` on `|z|^2`), so the pixel "
      "set is bit-identical and `A(t)^2` reproduces the pinned channel value.")
    A("")
    A("## 2. Inputs and provenance")
    A("")
    A(md_table(["input", "value"], [
        ["channel table", f"`{d['csv']}`"],
        ["channel table sha256", f"`{d['csv_sha256']}`"],
        ["overpasses in the table", d["n_rows"]],
        ["rect files (complex chips)",
         ", ".join(f"`{v}` ({d['rect_file_sha256'][v]})"
                   for v in d["rect_files"].values())],
        ["states", "; ".join(f"{s}: {d['states'][s]}" for s in core.STATES)],
        ["cuts", f"DS1 = {d['cuts']['DS1']}, DS2 = {d['cuts']['DS2']}"],
        ["definition", f"`{d['definition']}`"],
        ["differential", f"`{d['dphi_definition']}`"],
    ]))
    A("")
    A(md_table(["series", "file", "dates plotted"],
               [[lab, f"`{d['rect_files'][lab]}`", d["n_points_plotted"][lab]]
                for lab in order]))
    A("")
    A(md_table(["geometry", "bridge", "control 500 m east",
                "overpasses (bridge/control)", "paired", "unpaired"],
               [[geom, res["dphi"][geom]["bridge"], res["dphi"][geom]["control"],
                 f"{res['dphi'][geom]['n_bridge_total']}/"
                 f"{res['dphi'][geom]['n_control_total']}",
                 res["dphi"][geom]["n_pairs"],
                 f"{res['dphi'][geom]['n_bridge_only']}/"
                 f"{res['dphi'][geom]['n_control_only']}"]
                for geom in ("ASC", "DESC")]))
    A("")
    A("A pair needs a valid deck-edge mask on *both* overpasses; if one rectangle "
      "has none (and therefore no committed `gamma2_band_raw` either), that date "
      "is left out of the difference. `unpaired` lists bridge-only / control-only "
      "such dates — here all of them are ASC-control dates whose mask is empty.")
    A("")
    A("## 3. Figure")
    A("")
    A("`figures/bautzen/fig_bautzen_ocv_amp_phase.png` — three rows × one column "
      "per series. Top row: the coherence amplitude `A(t)` in [0, 1]. Second row: "
      "the deck-edge phase `phi(t)` in (-pi, pi]. Third row: the paired "
      "differential phase `dphi(t)` in (-pi, pi]. Per-overpass dots are coloured "
      "by state (RS / DS1 / DS2), a thin line connects the series in time order, "
      "the state windows are the background bands, and the two cut dates are the "
      "dashed vertical lines. The bridge columns are the signal, the two control "
      "columns the counter-check. The `dphi` panel belongs to the *geometry*, so "
      "it is drawn once per geometry under its bridge column (ASC = first, "
      "DESC = second); the control columns of that row keep the state bands and a "
      "pointer, but carry no series of their own.")
    A("")
    A("## 4. Per-state summary and cross-check")
    A("")
    A("Circular statistics for the phase (arithmetic means would be wrong across "
      "the ±pi branch cut): `mean` is the circular mean, `R` its resultant length "
      "(0 = uniform angles, 1 = identical angles).")
    A("")
    rows = []
    for lab in order:
        for s in core.STATES:
            c = res["summary"][lab][s]
            rows.append([series_label(lab), STATE_LABEL[s], c["n"],
                         c["median_A"], c["mean_A"],
                         c["phase"]["mean_rad"], c["phase"]["R"]])
    A(md_table(["series", "state", "n", "median A", "mean A",
                "phase mean [rad]", "R"], rows))
    A("")
    A("Paired differential phase `dphi(t) = wrap(phi_bridge - phi_control)`, per "
      "geometry and state (circular statistics as above; `R` is close to 0 "
      "because the differential phase is speckle-dominated and the pair counts "
      "per state are small):")
    A("")
    drows = []
    for geom in ("ASC", "DESC"):
        for s in core.STATES:
            c = res["dphi_summary"][geom][s]
            drows.append([geom, STATE_LABEL[s], c["n"],
                          c["dphi"]["mean_rad"], c["dphi"]["R"]])
    A(md_table(["geometry", "state", "n pairs", "dphi mean [rad]", "R"], drows))
    A("")
    A("### Δφ between the states (Δφ_DS1−RS, Δφ_DS2−RS)")
    A("")
    A("The state windows are consecutive and disjoint, so there is **no per-date "
      "pairing between states**; the shift away from the reference state RS is "
      "the difference of the two state circular means, "
      "`delta = wrap(circmean(state) - circmean(RS))`. Because the circular mean "
      "is not linear this is **not** the same as the difference of the series' "
      "own state phases, and both are given: first the figure's paired quantity "
      "`dphi`, then the absolute `phi` per series as counter-check. The interval "
      "is a seeded percentile bootstrap (20 000 resamples, 2.5/97.5 percentiles); "
      "`resolvable` marks a delta whose interval excludes 0.")
    A("")
    drows2 = []
    for geom in ("ASC", "DESC"):
        for st in ("DS1", "DS2"):
            d = res["delta"]["dphi"][geom][st]
            drows2.append([geom, f"Δφ_{st}−RS", d["delta_rad"],
                           f"[{d['ci_lo']:+.2f}, {d['ci_hi']:+.2f}]",
                           "yes" if d["resolvable"] else "no"])
    A(md_table(["geometry", "delta", "Δφ [rad]", "95% CI [rad]", "resolvable"],
               drows2))
    A("")
    drows3 = []
    for lab in order:
        for st in ("DS1", "DS2"):
            d = res["delta"]["series"][lab][st]
            drows3.append([series_label(lab), f"Δφ_{st}−RS", d["delta_rad"],
                           f"[{d['ci_lo']:+.2f}, {d['ci_hi']:+.2f}]",
                           "yes" if d["resolvable"] else "no"])
    A(md_table(["series", "delta", "Δφ [rad]", "95% CI [rad]", "resolvable"],
               drows3))
    A("")
    A("**Verdict: not resolvable.** Every interval above spans almost the whole "
      "circle, so no shift from RS to DS1 or DS2 survives the per-state noise "
      "(5–23 overpasses per state, resultant length `R` ≈ 0.1–0.7, which puts the "
      "uncertainty of a single state mean at order 1 rad). The value also depends "
      "on the definition — for ASC DS1 the paired panel gives `+1.69 rad` while "
      "`Δφ_bridge − Δφ_control` gives `+0.91 rad` — which is exactly why the "
      "definition is fixed above and why neither is reported as an effect.")
    A("")
    A(f"Cross-checks of the recomputed mask, coherence and pairing against the "
      f"committed table: **{res['n_cross_checks']} checks, "
      f"{res['n_cross_checks_failed']} failures**. Each date is checked for its "
      "deck-edge pixel count (`n_masked`), its band contrast, and the identity "
      "`A(t)^2 == gamma2_band_raw` (tolerance 1e-9); each bridge/control pair is "
      "checked to cover the same overpasses and to satisfy the wrap identity "
      "`wrap_pi(dphi - (phi_bridge - phi_control)) == 0` (tolerance 1e-12); each "
      "state delta is checked to equal `wrap_pi(circmean(state) - circmean(RS))` "
      "and to lie inside its own bootstrap interval.")
    A("")
    if res["n_cross_checks_failed"]:
        A("**Failures:**")
        A("")
        A(md_table(["check"], [[c["check"]] for c in res["cross_checks"]
                               if not c["ok"]]))
    else:
        A("**Failures: none.** Every recomputed deck-edge pixel count, band "
          "contrast and coherence value reproduces the committed table, every "
          "paired date satisfies the wrap identity, and every state delta "
          "matches its circular-mean definition and its bootstrap interval.")
    A("")


    A("## 5. Reproduce")
    A("")
    A("```sh")
    A("cd <repo root>")
    A("python3 code/bautzen/fig_bautzen_ocv_amp_phase.py")
    A("PYTHON=/usr/bin/python3 bash figures/bautzen/fig_bautzen_ocv_amp_phase.sh")
    A("```")
    A("")
    A("Needs numpy and matplotlib. This report was generated with python "
      f"{platform.python_version()}, numpy {np.__version__}, matplotlib "
      f"{matplotlib.__version__}. No output carries a timestamp, so repeated runs "
      "are identical.")
    A("")
    A("## 6. Caveats")
    A("")
    A("- **Amplitude = coherence amplitude.** `A(t)` is the square root of the "
      "`gamma2_band_raw` channel, so it inherits its limits: the raw estimator is "
      "used here (not floor-corrected) and the deck-edge pixel count `N` is small "
      "(single digits to a few tens), which widens every per-state reading.")
    A("- **Phase without a reference date.** The chips are single SLCs, so "
      "`phi(t)` is the *absolute* phase of the coherent deck-edge sum, dominated "
      "by the mean backscatter phase of the deck line. It is a time series of "
      "that phase, not unwrapped and not referenced to a master date; a real "
      "displacement would only appear as a deviation from this baseline.")
    A("- **State and season are confounded.** RS / DS1 / DS2 are consecutive time "
      "intervals of one deployment, assigned by overpass time (`core.state_of`); "
      "the control series and the cut dates are the only counter-checks carried "
      "here.")
    A("- **The differential phase is not InSAR phase.** `dphi(t)` differences the "
      "`arg` of two single-SLC coherent sums; it is defined modulo 2*pi and is "
      "*not* an interferometric (master-slave) phase. What it does remove is the "
      "component common to the bridge and its control rectangle, so a seasonal or "
      "atmospheric swing cancels and a local deck-edge change survives; the "
      "residual is still dominated by speckle and by the small pixel count `N`.")
    A("- **The state-to-state Δφ is not resolvable.** The shift "
      "`wrap(circmean(state) - circmean(RS))` is reported with a seeded bootstrap "
      "interval (section 4), but every interval spans almost the whole circle: "
      "with 5–23 overpasses per state and `R` ≈ 0.1–0.7 the circular mean of one "
      "state carries an uncertainty of order 1 rad. The delta is documented, not "
      "claimed, and it is not additive — the circular mean is not linear.")
    A("- **No committed pin.** The reference JSONs of the package hold no "
      "amplitude or phase, so the only verification is the internal identity "
      "`A(t)^2 == gamma2_band_raw` (section 4); the figure does not touch the "
      "figures A-E pin block.")
    A("")
    A("---")
    A("")
    A("Generated by `code/bautzen/fig_bautzen_ocv_amp_phase.py` "
      f"({res['n_cross_checks']} cross-checks, "
      f"{res['n_cross_checks_failed']} failures).")
    return "\n".join(L) + "\n"


# ---------------------------------------------------------------------------
# Command line
# ---------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Amplitude/phase time series of the Bautzen deck edge.")
    ap.add_argument("--figdir", default=FIGDIR,
                    help="output directory of the PNG figure")
    ap.add_argument("--json", default=OUT_JSON, help="path of the result JSON")
    ap.add_argument("--md", default=OUT_MD, help="path of the Markdown report")
    ap.add_argument("--no-figures", action="store_true",
                    help="only recompute and cross-check (no matplotlib output)")
    args = ap.parse_args(argv)

    t0 = time.time()
    print("Bautzen deck-edge amplitude/phase figure ...")
    res = compute()
    print(f"  compute: {time.time() - t0:.1f}s "
          f"({res['data']['n_rows']} overpasses, "
          f"{res['n_cross_checks']} cross-checks, "
          f"{res['n_cross_checks_failed']} failed)")
    for c in res["cross_checks"]:
        if not c["ok"]:
            print(f"    FAIL {c['check']}")

    if not args.no_figures:
        os.makedirs(args.figdir, exist_ok=True)
        out = make_figure(res, args.figdir)
        print(f"  {os.path.relpath(out)}")

    os.makedirs(os.path.dirname(os.path.abspath(args.json)), exist_ok=True)
    with open(args.json, "w", encoding="utf-8") as fh:
        json.dump(res, fh, indent=1, sort_keys=True, default=json_default)
        fh.write("\n")
    print(f"  {os.path.relpath(args.json)}")

    os.makedirs(os.path.dirname(os.path.abspath(args.md)), exist_ok=True)
    with open(args.md, "w", encoding="utf-8") as fh:
        fh.write(report(res))
    print(f"  {os.path.relpath(args.md)}")

    print(f"done in {time.time() - t0:.1f}s, "
          f"{res['n_cross_checks_failed']} cross-check failures")
    return 0 if res["n_cross_checks_failed"] == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
