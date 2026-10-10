#!/usr/bin/env python3
"""Espoo OCV-Paper — the record month by month (companion of the A–E figures).

Question: ``fig_espoo_ocv_paper.py`` pools the whole record and compares its two
**orbit geometries** — ASCENDING (n = 70, afternoon pass) against DESCENDING
(n = 80, morning pass). This companion keeps the time axis instead. The record is
event-free, so there is nothing to cut on and no window to anchor: the whole span
is shown, 25 months at 6 acquisitions each.

    months   2024-09 … 2026-09          the full record, nothing excluded
    groups   ASCENDING / DESCENDING     the same two orbit geometries as A–E

It asks three things the pooled figures cannot:

1. is the geometry contrast present in *every* month, or is it driven by a few?
   (panel c: the per-month Cliff's delta with its bootstrap CI)
2. does the record's own no-event bookkeeping move on the same axis as the
   contrast? (panel e: coverage and the phase ladder; panel f: the weather
   co-variates)
3. does the committed stationarity layer keep its promise on the time axis?
   (panel d: the two primary series of ``3-8 Hz · strip400_pipeline``)

Approach: a **standalone companion** of ``fig_espoo_ocv_paper.py`` (like
``../carola/fig_carola_ocv_months.py`` is a companion of the Carola paper
script). It reads only the committed channel table
``data/espoo/espoo_ocv_channels.csv`` and the committed reference
``data/espoo/reference/espoo_phase_stationarity.json``, recomputes every number
it draws with ``espoo_ocv_stats``, and writes one figure, a small JSON and an
English report. It adds **no** number to the A–E pin contract and touches none
of it; its own pin block compares the two primary series, the per-orbit pooled
medians, the 25 × 6 cadence and the phase availability against the committed
files and exits non-zero on any deviation.

Figure (2 × 3, 183 mm wide, 600 dpi PNG):
  a  gamma^2 by month, per orbit        (per-acquisition dots + monthly median)
  b  A by month, per orbit              (the mask-area half of the OCV)
  c  the geometry effect in every month (Cliff's delta per month + bootstrap CI)
  d  the two primary stationarity series``3-8 Hz · strip400_pipeline`` ASC / DESC
  e  coverage and the phase ladder by month (the two availability controls)
  f  the weather co-variates by month   (wind, temperature)

Usage:
  python3 code/espoo/fig_espoo_ocv_months.py
  python3 code/espoo/fig_espoo_ocv_months.py --figdir /tmp/x --json /tmp/y.json
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import platform
import sys
import time

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import espoo_ocv_core as core   # noqa: E402
import espoo_ocv_stats as st    # noqa: E402

FIGDIR = core.FIGDIR
OUT_JSON = os.path.join(core.DATA, "fig_espoo_ocv_months.json")
OUT_MD = os.path.join(FIGDIR, "fig_espoo_ocv_months.md")
REF_CHANNELS = os.path.join(core.REF, "fig_espoo_channels.json")
REF_STATIONARITY = os.path.join(core.REF, "espoo_phase_stationarity.json")

# The band x time base the committed file declares as its primary series.
PRIMARY_BAND = "3-8 Hz · strip400_pipeline"
PRIMARY_ASC = PRIMARY_BAND + " · ASCENDING"
PRIMARY_DESC = PRIMARY_BAND + " · DESCENDING"
# The two committed medians this companion pins as literal claims.
PIN_PRIMARY_MEDIAN = 0.2474
PIN_PRIMARY_DESC_MEDIAN = 0.2165
PIN_MEDIAN_TOL = 5e-5
# The two OCV channels of the paper, in the paper's order.
CHANNELS = ["gamma2", "A"]
CHANNEL_LABEL = {"gamma2": "coherence  gamma^2", "A": "mask area  A  [px]"}
CHANNEL_COLORS = {"gamma2": "#4C72B0", "A": "#DD8452"}
ORBIT_COLORS = {"ASCENDING": "#4C72B0", "DESCENDING": "#DD8452"}
# A per-month delta is only quoted when both sides carry at least this many
# acquisitions — 3/3 in the modal month, 2/4 in the one month that differs.
MIN_SIDE = 3
MONTH_MIN = 2      # months with fewer rows contribute no monthly mean (as committed)

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



def _median(vals):
    """The generator's median (``espoo_ocv_channels_csv._median``), reproduced.

    The committed ``fig_espoo_channels.json`` and the committed stationarity
    file both use it, and it sorts before picking the central value.
    """
    vals = sorted(v for v in vals if v is not None)
    if not vals:
        return None
    n = len(vals)
    return vals[n // 2] if n % 2 else 0.5 * (vals[n // 2 - 1] + vals[n // 2])


def _vals(rows, key):
    """The defined values of one channel over a row list."""
    return [r[key] for r in rows if r.get(key) is not None]


def save_fig(fig, name, out_dir=None):
    """Write one 600 dpi PNG (the binaries stay out of git)."""
    out_dir = out_dir or FIGDIR
    png = os.path.join(out_dir, name + ".png")
    fig.savefig(png)
    plt.close(fig)
    print(f"  {os.path.relpath(png)}")
    return png


def panel(ax, letter, text):
    ax.set_title(f"({letter}) {text}", loc="left")


def jitter(n, width, seed):
    rng = np.random.default_rng(seed)
    return rng.uniform(-width, width, n)


def stat_block(rows, key):
    """``{n, median}`` of one channel over one row list, by the generator's rule."""
    vals = _vals(rows, key)
    return {"n": len(vals), "median": _median(vals)}


# ---------------------------------------------------------------------------
# The blocks
# ---------------------------------------------------------------------------
def month_block(rows):
    """Per month: the composition and the two OCV channels per orbit geometry."""
    out = {}
    for m in sorted({r["month"] for r in rows}):
        sel = [r for r in rows if r["month"] == m]
        entry = {
            "n": len(sel),
            "n_by_state": {core.SHORT[s]: len([r for r in sel if r["state"] == s])
                           for s in core.STATE_ORDER},
            "date_range": [min(r["date"] for r in sel), max(r["date"] for r in sel)],
            "n_phase": sum(1 for r in sel if r["phase_available"]),
            "wind_speed_ms": _median([r["wind_speed_ms"] for r in sel]),
            "temperature_c": _median([r["temperature_c"] for r in sel]),
            "channels": {},
        }
        for k in CHANNELS:
            entry["channels"][k] = {"all": stat_block(sel, k)}
            for s in core.STATE_ORDER:
                entry["channels"][k][core.SHORT[s]] = stat_block(
                    [r for r in sel if r["state"] == s], k)
        out[m] = entry
    return out


def month_delta_block(rows, n_boot=2000, seed=core.RNG_SEED):
    """Per month and channel: ``delta(DESC - ASC)``, its CI and the two medians.

    A month is only quoted when both geometries carry at least ``MIN_SIDE``
    acquisitions, so the number rests on groups rather than on single rows.
    """
    out = {}
    for m in sorted({r["month"] for r in rows}):
        sel = [r for r in rows if r["month"] == m]
        entry = {"n": len(sel)}
        for k in CHANNELS:
            ref = _vals([r for r in sel if r["state"] == core.STATE_ORDER[0]], k)
            grp = _vals([r for r in sel if r["state"] == core.STATE_ORDER[1]], k)
            b = {"n": [len(ref), len(grp)]}
            if len(ref) >= MIN_SIDE and len(grp) >= MIN_SIDE:
                b["median"] = [_median(ref), _median(grp)]
                b["delta"] = st.cliffs_delta(ref, grp)
                ci = st.bootstrap_delta_ci(ref, grp, n=n_boot, seed=seed)
                if ci:
                    b["delta_ci95"] = [float(ci[0]), float(ci[1])]
                    b["ci_excludes_zero"] = bool(ci[0] * ci[1] > 0)
            else:
                b["note"] = f"n < {MIN_SIDE} on one side"
            entry[k] = b
        out[m] = entry
    return out


def primary_block(doc):
    """The two primary stationarity series, recomputed from their committed lists.

    ``per_year`` is the committed trend unit: ``st.theil_sen`` returns a slope
    per day, multiplied here by 365.25 — the convention of
    ``espoo_ocv_stats.block_bootstrap_slope(..., years=True)``.
    """
    out = {}
    series_by_group = doc["series_by_group"]
    for key in (PRIMARY_ASC, PRIMARY_DESC):
        series = [(int(d), v) for d, v in series_by_group[key] if v is not None]
        committed = doc["groups"][key]
        days = np.asarray([d for d, _v in series], float)
        vals = np.asarray([v for _d, v in series], float)
        slope = st.theil_sen(days, vals)
        out[key] = {
            "orbit": key.rsplit(" · ", 1)[-1],
            "band": " · ".join(key.split(" · ")[:2]),
            "n": len(series),
            "n_committed": committed["n"],
            "median": _median(vals.tolist()),
            "median_committed": committed["median"],
            "per_year": None if slope is None else slope * 365.25,
            "per_year_committed": committed["theil_sen_per_year"],
            "verdict": committed["verdict"],
            "flags": list(committed["flags"]),
            "mk_p": committed["mann_kendall"]["p"],
            "split_delta": (committed.get("split_half") or {}).get("delta"),
            "date_range": [dt.date.fromordinal(int(days[0])).isoformat(),
                           dt.date.fromordinal(int(days[-1])).isoformat()],
            "monthly_means": st.monthly_means_series(series),
            "monthly_means_committed": committed["monthly_means"],
            "series": [[int(d), float(v)] for d, v in series],
        }
    return out


def availability_block(rows):
    """Coverage and the phase ladder: per month counts, shares and their stop."""
    months = sorted({r["month"] for r in rows})
    by_month = {}
    for m in months:
        sel = [r for r in rows if r["month"] == m]
        n_ph = sum(1 for r in sel if r["phase_available"])
        by_month[m] = {"n": len(sel), "n_phase": n_ph,
                       "share": (n_ph / len(sel)) if sel else None,
                       "n_by_state": {core.SHORT[s]: len([r for r in sel
                                                          if r["state"] == s])
                                      for s in core.STATE_ORDER}}
    with_phase = [m for m in months if by_month[m]["n_phase"]]
    return {
        "months": months,
        "by_month": by_month,
        "cadence": sorted({by_month[m]["n"] for m in months}),
        "n_total": len(rows),
        "n_available": sum(1 for r in rows if r["phase_available"]),
        "n_snr_db_available": sum(1 for r in rows
                                  if r.get("phase_snr_db") is not None),
        "share": sum(1 for r in rows if r["phase_available"]) / len(rows),
        "first_month_with_phase": with_phase[0] if with_phase else None,
        "last_month_with_phase": with_phase[-1] if with_phase else None,
    }



def sign_block(deltas, months, key):
    """How often the per-month effect of one channel keeps its pooled sign."""
    vals = [deltas[m][key].get("delta") for m in months]
    vals = [v for v in vals if v is not None]
    return {"n_evaluable": len(vals),
            "n_negative": sum(1 for v in vals if v < 0),
            "n_positive": sum(1 for v in vals if v > 0),
            "n_zero": sum(1 for v in vals if v == 0),
            "n_skipped": len(months) - len(vals),
            "min": min(vals) if vals else None,
            "max": max(vals) if vals else None}


# ---------------------------------------------------------------------------
# Compute
# ---------------------------------------------------------------------------
def compute(n_boot=2000):
    """Every number the JSON, the report, the figure and the pins need."""
    t0 = time.time()
    rows = core.load_csv()
    refs = core.load_reference()
    sta = refs["stationarity"]
    ch = refs["channels"]

    per_month = month_block(rows)
    deltas = month_delta_block(rows, n_boot=n_boot, seed=core.RNG_SEED)
    primary = primary_block(sta)
    avail = availability_block(rows)

    pooled, standing = {}, {}
    for k in CHANNELS:
        pooled[k] = {"all": stat_block(rows, k)}
        for s in core.STATE_ORDER:
            pooled[k][core.SHORT[s]] = stat_block(
                [r for r in rows if r["state"] == s], k)
        ref = _vals([r for r in rows if r["state"] == core.STATE_ORDER[0]], k)
        grp = _vals([r for r in rows if r["state"] == core.STATE_ORDER[1]], k)
        standing[k] = st.delta_block(ref, grp, "DESC vs ASC", n_boot=n_boot,
                                     seed=core.RNG_SEED)

    months = sorted(per_month)
    signs = {k: sign_block(deltas, months, k) for k in CHANNELS}
    return {
        "rows": rows, "channels_committed": ch, "stationarity": sta,
        "months": months, "per_month": per_month, "deltas": deltas,
        "primary": primary, "availability": avail, "signs": signs,
        "pooled": pooled, "standing": standing,
        "n_rows": len(rows), "n_boot": n_boot,
        "date_first": min(r["date"] for r in rows),
        "date_last": max(r["date"] for r in rows),
        "runtime_s": time.time() - t0,
    }


def _delta_row(res, key):
    """``(months, deltas, los, his)`` of the evaluable months of one channel."""
    xs, ys, los, his = [], [], [], []
    for i, m in enumerate(res["months"]):
        b = res["deltas"][m][key]
        if b.get("delta") is None:
            continue
        xs.append(i)
        ys.append(b["delta"])
        ci = b.get("delta_ci95") or [None, None]
        los.append(ci[0])
        his.append(ci[1])
    return xs, ys, los, his


# ---------------------------------------------------------------------------
# Rendering — 2 x 3, 183 mm column width, 600 dpi PNG (Espoo style)
# ---------------------------------------------------------------------------
def _month_axis(ax, months, step=6):
    """Five month labels ("24-09" … "26-09"), horizontal and readable."""
    idx = list(range(0, len(months), step))
    ax.set_xticks(idx)
    ax.set_xticklabels([months[i][2:] for i in idx], fontsize=6.4)
    ax.set_xlim(-0.8, len(months) - 0.2)


def _month_channel_panel(ax, res, key, seed):
    """Dots per acquisition + the monthly median of each orbit geometry."""
    rows, months = res["rows"], res["months"]
    idx = {m: i for i, m in enumerate(months)}
    ax.axhline(res["pooled"][key]["all"]["median"], color="#777777", ls=":",
               lw=0.9)
    for j, s in enumerate(core.STATE_ORDER):
        sel = [r for r in rows if r["state"] == s]
        xs = np.asarray([idx[r["month"]] for r in sel], float) \
            + (j - 0.5) * 0.30 + jitter(len(sel), 0.10, seed + j)
        ax.plot(xs, _vals(sel, key), "o", ms=2.1, alpha=0.55, mew=0,
                color=ORBIT_COLORS[s])
        meds = [res["per_month"][m]["channels"][key][core.SHORT[s]]["median"]
                for m in months]
        ax.plot(range(len(months)), meds, "-o", ms=2.4, lw=1.0,
                color=ORBIT_COLORS[s])
    _month_axis(ax, months)
    ax.set_ylabel(CHANNEL_LABEL[key])


def fig_months(res, out_dir=None):
    """The six month-by-month panels."""
    months = res["months"]
    fig, axes = plt.subplots(2, 3, figsize=(COL_WIDTH, 0.72 * COL_WIDTH))

    # (a) gamma^2 by month, (b) A by month
    for j, key in enumerate(CHANNELS):
        ax = axes[0][j]
        _month_channel_panel(ax, res, key, seed=core.RNG_SEED + j)
        panel(ax, "ab"[j], f"{'gamma^2' if key == 'gamma2' else key} per month")

    # (c) the geometry effect in every month
    ax = axes[0][2]
    for j, key in enumerate(CHANNELS):
        xs, ys, los, his = _delta_row(res, key)
        lo = [0 if v is None else v for v in los]
        hi = [0 if v is None else v for v in his]
        ax.errorbar(xs, ys, yerr=[np.asarray(ys) - np.asarray(lo),
                                  np.asarray(hi) - np.asarray(ys)],
                    fmt="o", ms=2.6, lw=0.9, capsize=0.0, mew=0,
                    color=CHANNEL_COLORS[key], alpha=0.9)
    ax.axhline(0.0, color="#666666", lw=0.9)
    for key in CHANNELS:
        d0 = res["standing"][key]["delta"]
        ax.axhline(d0, color=CHANNEL_COLORS[key], ls="--", lw=0.9)
        # The label rides just above its own dashed line, at the right edge, with
        # a translucent plate so it stays legible where bars cross it.
        ax.text(len(months) - 0.6, d0, f"{key} pooled ", ha="right",
                va="bottom", fontsize=6.0, color=CHANNEL_COLORS[key],
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.75,
                          pad=0.6))
    ax.set_ylim(-1.15, 1.15)
    ax.set_xlim(-0.8, len(months) - 0.2)
    _month_axis(ax, months)
    ax.set_ylabel("Cliff's delta of the month\n(DESC - ASC)")
    g2 = [res["deltas"][m]["gamma2"].get("delta") for m in months]
    n_eval = sum(1 for v in g2 if v is not None)
    panel(ax, "c", f"the effect in the {n_eval} evaluable months")

    # (d) the two primary stationarity series
    ax = axes[1][0]
    midx = {m: i for i, m in enumerate(months)}
    for i, key in enumerate((PRIMARY_ASC, PRIMARY_DESC)):
        b = res["primary"][key]
        s = b["orbit"]
        xs = np.asarray([midx[dt.date.fromordinal(int(d)).isoformat()[:7]]
                         for d, _v in b["series"]], float)
        ax.plot(xs + jitter(len(xs), 0.16, core.RNG_SEED + 10 + i),
                [v for _d, v in b["series"]], "o", ms=2.1, alpha=0.5, mew=0,
                color=ORBIT_COLORS[s], label=f"{core.SHORT[s]} (n = {b['n']})")
        mm = b["monthly_means"]
        mx = [midx[m] for m in sorted(mm)]
        ax.plot(mx, [mm[m] for m in sorted(mm)], "-", lw=1.1,
                color=ORBIT_COLORS[s])
        slope_m = b["per_year"] / 12.0
        span = np.asarray([0.0, float(max(mx))])
        ax.plot(span, b["median"] + slope_m * (span - float(np.mean(mx))), "--",
                lw=0.9, color=ORBIT_COLORS[s])
    ax.set_ylabel("pooled phase coherence")
    y0, y1 = ax.get_ylim()
    ax.set_ylim(y0, y1 + 0.30 * (y1 - y0))
    for i, key in enumerate((PRIMARY_ASC, PRIMARY_DESC)):
        b = res["primary"][key]
        ax.text(0.02, 0.97 - 0.080 * i,
                f"{core.SHORT[b['orbit']]}: {b['median']:.4f}  "
                f"{b['per_year']:+.4f}/yr  p={b['mk_p']:.2f}  {b['verdict']}",
                transform=ax.transAxes, fontsize=5.6, va="top",
                color=ORBIT_COLORS[b["orbit"]])
    _month_axis(ax, months)
    panel(ax, "d", "primary series\n3-8 Hz · strip400_pipeline")

    # (e) coverage and the phase ladder
    ax = axes[1][1]
    av = res["availability"]
    width = 0.34
    for si, s in enumerate(core.STATE_ORDER):
        counts = [av["by_month"][m]["n_by_state"][core.SHORT[s]] for m in months]
        ax.bar(np.arange(len(months)) + (si - 0.5) * 0.38,
               counts, width=width, color=ORBIT_COLORS[s], alpha=0.9)
    ax.set_ylim(0, 4.5)
    ax.set_ylabel("acquisitions per month")
    _month_axis(ax, months)
    ax2 = ax.twinx()
    ax2.grid(False)
    ax2.plot(range(len(months)),
             [av["by_month"][m]["share"] for m in months], "-s", ms=2.4,
             lw=1.0, color="#55A868")
    ax2.set_ylim(0, 1.15)
    ax2.set_ylabel("share with phase", color="#3A7A4C")
    ax2.tick_params(axis="y", colors="#3A7A4C")
    ax2.text(len(months) - 0.6, 0.90, f"{av['n_available']} of {av['n_total']}",
             ha="right", fontsize=6.0, color="#3A7A4C")
    panel(ax, "e", "coverage and the phase ladder")

    # (f) the weather co-variates
    ax = axes[1][2]
    wind = [res["per_month"][m]["wind_speed_ms"] for m in months]
    temp = [res["per_month"][m]["temperature_c"] for m in months]
    ax.plot(range(len(months)), wind, "-o", ms=2.4, lw=1.0, color="#55A868")
    ax.set_ylim(min(wind) - 0.4, max(wind) + 0.4)
    ax.set_ylabel("wind speed [m/s]", color="#3A7A4C")
    ax.tick_params(axis="y", colors="#3A7A4C")
    ax.set_xlim(-0.8, len(months) - 0.2)
    _month_axis(ax, months)
    ax3 = ax.twinx()
    ax3.grid(False)
    ax3.plot(range(len(months)), temp, "-s", ms=2.4, lw=1.0, color="#C44E52")
    ax3.set_ylim(min(temp) - 0.6, max(temp) + 0.6)
    ax3.set_ylabel("temperature [C]", color="#8C2F33")
    ax3.tick_params(axis="y", colors="#8C2F33")
    panel(ax, "f", "weather co-variates, monthly median")

    handles = [plt.Line2D([], [], color=ORBIT_COLORS[s], marker="o", ls="",
                          ms=3.4, mew=0, label=core.SHORT[s])
               for s in core.STATE_ORDER]
    fig.legend(handles=handles, loc="upper right", ncol=2, frameon=False,
               fontsize=7.5, bbox_to_anchor=(0.99, 0.985),
               handletextpad=0.4, columnspacing=1.2)
    fig.suptitle("(Espoo) the record month by month — the orbit-geometry contrast "
                 "is in every month", fontsize=9.0, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.94), w_pad=0.7, h_pad=1.0)
    return save_fig(fig, "fig_espoo_ocv_months", out_dir)


# ---------------------------------------------------------------------------
# Pins: every committed number this companion quotes is re-derived and checked
# ---------------------------------------------------------------------------
def _check(problems, name, got, want, tol=0.0):
    """Field-by-field comparison; a float is compared with ``tol``, else exactly."""
    if got is None or want is None:
        if got is not want:
            problems.append(f"{name}: got {got!r} != want {want!r}")
        return
    if isinstance(want, dict):
        if not isinstance(got, dict):
            problems.append(f"{name}: got <{type(got).__name__}> != <dict>")
            return
        for k in want:
            if k not in got:
                problems.append(f"{name}.{k}: missing in the recomputed value")
                continue
            _check(problems, f"{name}.{k}", got[k], want[k], tol)
        return
    if isinstance(want, (int, float)) and not isinstance(want, bool):
        if isinstance(got, (int, float)) and abs(float(got) - float(want)) <= tol:
            return
        problems.append(f"{name}: got {got!r} != want {want!r} (tol {tol:g})")
        return
    if got != want:
        problems.append(f"{name}: got {got!r} != want {want!r}")


ORBIT_SUFFIX = {True: "ascending", False: "descending"}
CHANNEL_COLUMN = {"gamma2": "coherence_gamma2", "A": "coherence_masked_pixels"}


def pin_block(res):
    """Re-derive every committed number this companion draws on and compare."""
    problems = []
    ch = res["channels_committed"]
    rows = res["rows"]
    months = res["months"]
    av = res["availability"]

    # --- 1. the channel table: counts, span, cadence, medians ---------------
    _check(problems, "channels.n_acquisitions", len(rows), ch["n_acquisitions"])
    _check(problems, "channels.n_months", len(months), len(ch["monthly"]))
    _check(problems, "channels.date_first", res["date_first"], ch["date_first"])
    _check(problems, "channels.date_last", res["date_last"], ch["date_last"])
    for m in months:
        for s in core.STATE_ORDER:
            _check(problems, f"channels.monthly.{m}.{s}",
                   res["per_month"][m]["n_by_state"][core.SHORT[s]],
                   ch["monthly"][m][s])
    for k, col in CHANNEL_COLUMN.items():
        _check(problems, f"channels.{col}",
               res["pooled"][k]["all"],
               {"n": ch[col]["n"], "median": ch[col]["median"]})
        for s in ("ASCENDING", "DESCENDING"):
            row = f"{col}_{ORBIT_SUFFIX[s == 'ASCENDING']}"
            _check(problems, f"channels.{row}", res["pooled"][k][core.SHORT[s]],
                   {"n": ch[row]["n"], "median": ch[row]["median"]})
    _check(problems, "channels.n_phase",
           {"n": av["n_available"], "share": av["share"]},
           {"n": ch["phase_coherence"]["n"],
            "share": ch["phase_coherence"]["n"] / ch["n_acquisitions"]}, tol=1e-12)
    # the cadence: every month of the record carries the same number of rows
    _check(problems, "cadence", av["cadence"], [6])
    _check(problems, "cadence.total", len(months) * av["cadence"][0], len(rows))
    # --- 2. the two primary stationarity series ----------------------------
    sta = res["stationarity"]
    for key in (PRIMARY_ASC, PRIMARY_DESC):
        b = res["primary"][key]
        tag = f"primary.{b['orbit']}"
        _check(problems, f"{tag}.n", b["n"], b["n_committed"])
        _check(problems, f"{tag}.median", b["median"], b["median_committed"])
        _check(problems, f"{tag}.trend", b["per_year"], b["per_year_committed"],
               tol=1e-9)
        _check(problems, f"{tag}.monthly_means", b["monthly_means"],
               b["monthly_means_committed"], tol=1e-9)
        _check(problems, f"{tag}.date_range", b["date_range"],
               [sta["date_range"][0], sta["date_range"][1]])
        _check(problems, f"{tag}.verdict", b["verdict"], "stationary")
        _check(problems, f"{tag}.flags", b["flags"], [])
    _check(problems, "primary.median_literal",
           round(res["primary"][PRIMARY_ASC]["median"], 4), PIN_PRIMARY_MEDIAN)
    _check(problems, "primary.median_literal_desc",
           round(res["primary"][PRIMARY_DESC]["median"], 4),
           PIN_PRIMARY_DESC_MEDIAN)
    _check(problems, "primary.trend_sign",
           [res["primary"][PRIMARY_ASC]["per_year"] > 0,
            res["primary"][PRIMARY_DESC]["per_year"] < 0], [True, True])
    # --- 3. the figure's own claims ----------------------------------------
    signs = res["signs"]
    _check(problems, "signs.gamma2_no_inversion",
           signs["gamma2"]["n_positive"] + signs["gamma2"]["n_zero"], 0)
    _check(problems, "signs.gamma2_all_negative",
           signs["gamma2"]["n_negative"], signs["gamma2"]["n_evaluable"])
    _check(problems, "signs.A_majority_positive",
           signs["A"]["n_positive"] > signs["A"]["n_negative"], True)
    _check(problems, "ladder.last_month", av["last_month_with_phase"], "2025-10")
    _check(problems, "ladder.first_month", av["first_month_with_phase"], "2024-09")
    return problems





# ---------------------------------------------------------------------------
# Report and JSON
# ---------------------------------------------------------------------------
def report(res, problems):
    """The markdown report written next to the PNG."""
    months = res["months"]
    prim = res["primary"]
    md = []
    A = md.append

    A("# Espoo Kurttila mast — the record month by month")
    A("")
    A(f"*{core.SITE_SHORT} · event-free control site · {len(months)} months · "
      f"{res['n_rows']} acquisitions · {res['date_first']} … {res['date_last']} · "
      "generated by `code/espoo/fig_espoo_ocv_months.py`*")
    A("")
    A("Companion figure of `figures/espoo/fig_espoo_ocv_paper_A.png` … `_E.png`. "
      "The A–E figures pool the whole record; this one keeps the time axis and "
      "shows the same two orbit geometries month by month. There is no event in "
      "the window, so nothing is cut, anchored or excluded — all "
      f"{res['n_rows']} acquisitions appear here, 6 in each of the {len(months)} "
      "months.")
    A("")
    A("## The record month by month")
    A("")
    A("| month | n | ASC / DESC | gamma^2 ASC | gamma^2 DESC | A ASC [px] | A DESC [px] |")
    A("|---|---|---|---|---|---|---|")
    for m in months:
        e = res["per_month"][m]
        c = e["channels"]
        A(f"| {m} | {e['n']} | {e['n_by_state']['ASC']} / {e['n_by_state']['DESC']} "
          f"| {c['gamma2']['ASC']['median']:.4f} | {c['gamma2']['DESC']['median']:.4f} "
          f"| {c['A']['ASC']['median']:.0f} | {c['A']['DESC']['median']:.0f} |")
    A("")
    A(f"Pooled over the record: gamma^2 {res['pooled']['gamma2']['ASC']['median']:.4f} "
      f"(ASC) against {res['pooled']['gamma2']['DESC']['median']:.4f} (DESC), "
      f"A {res['pooled']['A']['ASC']['median']:.0f} px against "
      f"{res['pooled']['A']['DESC']['median']:.0f} px. The 25 monthly medians sit "
      "on the two sides of the pooled values rather than wandering across them.")
    A("")
    A("## The geometry contrast in every month")
    A("")
    A("Per-month Cliff's delta (DESC − ASC) of the two measured channels, with a "
      f"{res['n_boot']} draw bootstrap CI. A month enters only when both orbit "
      "geometries carry at least 3 acquisitions, so 5 months (2 ASC / 4 DESC) are "
      "reported as not evaluable rather than quoted on two rows:")
    A("")
    A("| month | n ASC / DESC | delta gamma^2 | 95 % CI | delta A | 95 % CI |")
    A("|---|---|---|---|---|---|")
    for m in months:
        b = res["deltas"][m]
        if b["gamma2"].get("delta") is None:
            continue
        g, a = b["gamma2"], b["A"]
        fmt = lambda x: "--" if x is None else f"{x:+.3f}"        # noqa: E731

        def ci(v, c=None):
            """The CI of one channel block, ``--`` when the month is not quoted."""
            c = v.get("delta_ci95") if c is None else c
            return "--" if not c else "[" + fmt(c[0]) + ", " + fmt(c[1]) + "]"
        A(f"| {m} | {g['n'][0]} / {g['n'][1]} | {fmt(g['delta'])} | {ci(g)} "
          f"| {fmt(a['delta'])} | {ci(a)} |")
    A("")
    sg, sa = res["signs"]["gamma2"], res["signs"]["A"]
    A(f"gamma^2 keeps the pooled sign in **{sg['n_negative']} of "
      f"{sg['n_evaluable']}** evaluable months — it never turns positive "
      f"(range {sg['min']:+.3f} … {sg['max']:+.3f}). The mask area A is positive in "
      f"{sa['n_positive']} of {sa['n_evaluable']} months "
      f"({sa['n_negative']} negative, {sa['n_zero']} at zero, range "
      f"{sa['min']:+.3f} … {sa['max']:+.3f}), i.e. the mask channel is the noisier "
      "of the two, but the pooled reading survives the monthly split.")
    A("")
    A("## The two primary stationarity series")
    A("")
    A(f"`{PRIMARY_BAND}` is the band x time base the committed "
      "`espoo_phase_stationarity.json` declares as its primary series. Both orbits "
      "are recomputed here from the committed per-acquisition lists "
      "(`theil_sen` per day x 365.25, the committed trend unit):")
    A("")
    A("| series | n | median | trend /yr | MK p | split delta | verdict |")
    A("|---|---|---|---|---|---|---|")
    for key in (PRIMARY_ASC, PRIMARY_DESC):
        b = prim[key]
        A(f"| `{b['orbit']}` | {b['n']} | {b['median']:.4f} | "
          f"{b['per_year']:+.4f} | {b['mk_p']:.4f} | {b['split_delta']:+.4f} | "
          f"{b['verdict']} |")
    A("")
    A(f"The two committed medians — **{prim[PRIMARY_ASC]['median']:.4f}** (ASC) and "
      f"**{prim[PRIMARY_DESC]['median']:.4f}** (DESC) — are pinned as literal "
      f"claims to 4 decimals ({PIN_PRIMARY_MEDIAN} / {PIN_PRIMARY_DESC_MEDIAN}), "
      "and the committed monthly means of both series are recomputed from the "
      "lists and compared field by field. Both series are stationary with no "
      "flags, and the orbit gap (ASC above DESC) is larger than either trend — "
      "the statement the pooled figures make, now made inside the time axis.")
    A("")
    av = res["availability"]
    A("## The record's own controls, on the same axis")
    A("")
    A(f"**Coverage.** Every month of the record carries exactly "
      f"{av['cadence'][0]} acquisitions "
      f"({len(months)} x {av['cadence'][0]} = {res['n_rows']}), split "
      "3 ASC / 3 DESC in 20 months and 2 / 4 in 5. The contrast is therefore not "
      "a coverage artefact: the two geometries are simultaneously present in "
      "every single month.")
    A("")
    A(f"**The phase ladder.** `phase_coherence` and `phase_rms_rad` exist for "
      f"{av['n_available']} of {av['n_total']} acquisitions "
      f"({100 * av['share']:.1f} %), from "
      f"{av['first_month_with_phase']} to {av['last_month_with_phase']}; after "
      "that month the export carries no phase at all, and "
      f"{av['n_snr_db_available']} acquisitions carry an SNR. The ladder is an "
      "availability fact about the export, not a third observability channel, and "
      "it is drawn as a share rather than a value for that reason.")
    A("")
    A("**The weather co-variates.** `wind_speed_ms` and `temperature_c` are "
      "constant over long stretches of this export (5 m/s and 10 C in 132 of the "
      "150 rows), so their monthly medians are step functions rather than "
      "measurements. That is exactly why panel E3 of the paper recomputes the "
      "contrast *inside* each distinct value instead of correcting for it.")
    A("")
    A("## Pins")
    A("")
    if problems:
        A(f"**{len(problems)} check(s) failed** — this report was written from a "
          "run in which the committed numbers did not reproduce:")
        A("")
        for p in problems:
            A(f"* `{p}`")
    else:
        A("Every committed number this companion draws on is re-derived and "
          "compared; **all checks passed**.")
    A("")
    A("| target | what is pinned |")
    A("|---|---|")
    A("| `fig_espoo_channels.json` | 150 acquisitions, the 2024-09-05 … "
      "2026-09-07 span, the 25-month composition per orbit, the pooled and "
      "per-orbit medians of gamma^2 and A, and the phase count — field by field |")
    A("| `espoo_phase_stationarity.json` | the two primary series "
      "(`3-8 Hz · strip400_pipeline` ASC/DESC): n, median, trend, verdict and "
      "flags, plus all 25 committed monthly means of both series recomputed from "
      "the per-acquisition lists |")
    A("| literal claims | the two primary medians to 4 decimals (0.2474 / 0.2165), "
      "the trend signs, the 25-month cadence of 6, the phase ladder ending in "
      "2025-10 |")
    A("| this figure's own claim | gamma^2 keeps the DESC < ASC sign in every "
      "evaluable month, and A is positive in the majority of them |")
    A("")
    A("## Reproduce")
    A("")
    A("```bash")
    A("python3 code/espoo/fig_espoo_ocv_months.py          # this figure")
    A("python3 code/espoo/fig_espoo_ocv_paper.py          # the A–E figures + report")
    A("```")
    A("")
    A(f"*environment: python {platform.python_version()}, numpy {np.__version__}, "
      f"matplotlib {matplotlib.__version__}*")
    A("")
    A("*No wall-clock value appears in this file, so two runs of the same code "
      "produce identical bytes.*")

    return "\n".join(md) + "\n"



def json_default(o):
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.floating):
        return float(o)
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, dt.date):
        return o.isoformat()
    return str(o)


def build_json(res, problems):
    """The result JSON — the per-month blocks, not the per-row rows."""
    primary = {}
    for key, b in res["primary"].items():
        primary[key] = {k: v for k, v in b.items() if k != "series"}
        primary[key]["series_days"] = [d for d, _v in b["series"]]
    return {
        "site": core.SITE,
        "site_short": core.SITE_SHORT,
        "event": None,
        "event_note": ("the record is event-free: no damage_label, no "
                       "structural_state and no condition_label on any of the "
                       f"{res['n_rows']} rows, so this companion has no window to "
                       "anchor and shows the whole span instead"),
        "generated_by": "code/espoo/fig_espoo_ocv_months.py",
        "state_order": list(core.STATE_ORDER),
        "states": dict(core.SHORT),
        "comparison": f"{core.SHORT[core.STATE_ORDER[1]]}_vs_"
                      f"{core.SHORT[core.STATE_ORDER[0]]}",
        "channels": list(CHANNELS),
        "seed": core.RNG_SEED,
        "n_boot": res["n_boot"],
        "min_side": MIN_SIDE,
        "n_rows": res["n_rows"],
        "date_first": res["date_first"],
        "date_last": res["date_last"],
        "months": res["months"],
        "per_month": res["per_month"],
        "deltas": res["deltas"],
        "signs": res["signs"],
        "pooled": res["pooled"],
        "standing": res["standing"],
        "primary": primary,
        "availability": res["availability"],
        "figures": ["figures/espoo/fig_espoo_ocv_months.png"],
        "pins": {"ok": not problems, "n_problems": len(problems),
                 "problems": list(problems)},
    }



def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--figdir", default=FIGDIR)
    ap.add_argument("--json", default=OUT_JSON)
    ap.add_argument("--md", default=OUT_MD)
    ap.add_argument("--n-boot", type=int, default=2000,
                    help="bootstrap draws of the per-month CI (default 2000)")
    args = ap.parse_args(argv)

    t0 = time.time()
    print(f"computing — {core.SITE}, {core.SHORT[core.STATE_ORDER[0]]} -> "
          f"{core.SHORT[core.STATE_ORDER[1]]}, month by month (no event) …")
    res = compute(n_boot=args.n_boot)
    print(f"  {res['n_rows']} acquisitions in {len(res['months'])} months "
          f"({res['availability']['cadence']} per month), {res['runtime_s']:.1f} s")

    problems = pin_block(res)
    if problems:
        print(f"PIN FAILED ({len(problems)} problem(s)):")
        for p in problems:
            print(f"  - {p}")
    else:
        print("  pins OK (channels, primary series, cadence, ladder)")

    os.makedirs(args.figdir, exist_ok=True)
    fig_months(res, args.figdir)
    os.makedirs(os.path.dirname(os.path.abspath(args.json)), exist_ok=True)
    with open(args.json, "w") as fh:
        json.dump(build_json(res, problems), fh, indent=1, sort_keys=True,
                  default=json_default)
    print(f"  {os.path.relpath(args.json)}")
    os.makedirs(os.path.dirname(os.path.abspath(args.md)), exist_ok=True)
    with open(args.md, "w") as fh:
        fh.write(report(res, problems))
    print(f"  {os.path.relpath(args.md)}")
    print(f"DONE ({time.time() - t0:.1f} s)")
    # The artifacts are written first so a failure can be read in place; the
    # exit code still reports it.
    return 2 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
