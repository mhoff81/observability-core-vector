#!/usr/bin/env python3
"""Bautzen OCV time series by dimension and orbit.

The paper figures A-E of this package summarise the six OCV channels per state
and per rect series.  This companion asks the same channels as a *time series*
and separates the two Sentinel-1 pass geometries: one row per OCV dimension
(``gamma2``, ``P``, ``D``, ``F``, ``A``, ``S``) and one column per orbit
(ASC / DESC), i.e. twelve panels in a 6 x 2 grid.  The row order is the formal
OCV definition ``[gamma2, P, D, F, A, S]``, not the column order of the CSV.

The figure is drawn from the committed channel table
``data/bautzen/bautzen_ocv_channels.csv`` only.  The x coordinate is the
**actual local overpass instant** (``core.acq_dt``), not midnight: the state
boundaries fall in the middle of a day (DESC flies at 07:08, ASC at 18:51
local), so a date-only axis would put the boundary-day points into the wrong
state band.  The three state windows RS / DS1 / DS2 are background bands with
different colours, separated by dashed lines at ``core.CUT_DS1`` and
``core.CUT_DS2``.

Only the bridge series is drawn.  The control rectangle 500 m east is the
counter-check of the paper figure, but it is deliberately left out here: its
ASC ``gamma2`` range (-0.94 .. +0.93) is much wider than the bridge range
(-0.14 .. +0.24), so a shared panel would compress exactly the bridge signal
this figure is meant to show.  The control remains in the channel table and in
figures A-E.

Missing channel values are never interpolated.  The bridge series happens to be
complete; any gap would be drawn as a gap.  ``gamma2`` is the committed
floor-corrected deck-band coherence (``gamma2_band_floor_corrected``);
``A`` is the deck-edge pixel count, ``D`` its density, ``F`` its 8-neighbourhood
fragmentation, ``S`` the centroid-to-peak shift and ``P`` its persistence -- the
same definitions as ``bautzen_ocv_core.DIMS``.

Everything is recomputed from the committed table and checked against
``bautzen_ocv_channels_meta.json``: the row count, the state rule, the per-series
date/state counts and the rect-file sha256s.  Any failed check makes the run
exit non-zero.  No output carries a wall-clock timestamp, so repeated runs are
byte-identical.

Usage:
  python3 code/bautzen/fig_bautzen_ocv_orbits.py
  python3 code/bautzen/fig_bautzen_ocv_orbits.py --no-figures
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

import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import bautzen_ocv_core as core    # noqa: E402

FIGDIR = core.FIGDIR
FIG_NAME = "fig_bautzen_ocv_orbits"
FIG_DPI = 630
OUT_JSON = os.path.join(core.DATA, FIG_NAME + ".json")
OUT_MD = os.path.join(FIGDIR, FIG_NAME + ".md")

# ---------------------------------------------------------------------------
# Style (the same block as the other Bautzen companion figures)
# ---------------------------------------------------------------------------
plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "mathtext.fontset": "stix",
    "font.size": 8.5,
    "axes.titlesize": 8.0,
    "axes.labelsize": 8.0,
    "xtick.labelsize": 7.0,
    "ytick.labelsize": 7.0,
    "legend.fontsize": 7.2,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.alpha": 0.25,
    "grid.linewidth": 0.4,
    "figure.dpi": 110,
    "savefig.dpi": FIG_DPI,
    "savefig.bbox": "tight",
    "pdf.fonttype": 42,
})
MM = 1.0 / 25.4
COL_WIDTH = 183 * MM

# State colours/markers (three states, in the severity order of STATE_ORDER).
STATE_COLORS = {"RS": "#4C72B0", "DS1": "#DD8452", "DS2": "#C44E52"}
STATE_MARKERS = {"RS": "o", "DS1": "s", "DS2": "^"}
STATE_LABEL = {"RS": "RS (reference)", "DS1": "DS1", "DS2": "DS2"}
# Pastel fills make the three windows distinguishable without hiding the data.
STATE_FILL = {"RS": "#E4EEF8", "DS1": "#FBE7D5", "DS2": "#F6D8D8"}

ORBITS = ("ASC", "DESC")
ORBIT_LABEL = {
    "ASC": "ASC — IW2, ~18:51 local",
    "DESC": "DESC — IW3, ~07:08 local",
}
ORBIT_COLOR = {"ASC": "#4C72B0", "DESC": "#55A868"}

# The formal OCV definition order: [gamma2, P, D, F, A, S].  This is the panel
# row order of the figure; it deliberately differs from the CSV column order.
CHANNELS = ["gamma2", "P", "D", "F", "A", "S"]
CHANNEL_AXIS = {
    "gamma2": "gamma^2  [1]",
    "P": "P  [1]",
    "D": "D  [1]",
    "F": "F  [components]",
    "A": "A  [px]",
    "S": "S  [px]",
}
FIG_CHANNEL_AXIS = {k: v.replace("gamma^2", r"$\gamma^2$") for k, v in CHANNEL_AXIS.items()}
CHANNEL_TITLE = {
    "gamma2": "deck-band coherence",
    "P": "mask persistence",
    "D": "mask density",
    "F": "mask fragmentation",
    "A": "mask area",
    "S": "centroid-to-peak shift",
}
FIG_CHANNEL_TITLE = {k: v.replace("gamma^2", r"$\gamma^2$") for k, v in CHANNEL_TITLE.items()}
PANEL_LETTERS = "abcdefghijkl"


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------
def json_default(o):
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    raise TypeError(f"not JSON serialisable: {type(o)!r}")


def fmt(x, digits=3):
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


def finite(value):
    return value is not None and math.isfinite(float(value))


def check(checks, name, got, want):
    """Append one pass/fail entry; failures carry both sides for the report."""
    entry = {"check": name, "ok": bool(got == want)}
    if not entry["ok"]:
        entry["recomputed"] = repr(got)[:300]
        entry["committed"] = repr(want)[:300]
    checks.append(entry)


# ---------------------------------------------------------------------------
# Recompute the plotted payload and the committed cross-checks
# ---------------------------------------------------------------------------
def compute(csv_path=None, meta_path=None):
    """Load the channel table, rebuild every plotted point and pin the counts."""
    rows = core.load_csv(csv_path) if csv_path else core.load_csv()
    meta = core.load_csv_meta(meta_path) if meta_path else core.load_csv_meta()
    csv_path = csv_path or core.CSV_PATH
    meta_path = meta_path or core.CSV_META_PATH
    # The rect files live beside the table that enumerates them, so the sha pins
    # follow the CSV: a ring table must check the ring chips, not the committed
    # ones it was derived from.
    data_dir = os.path.dirname(os.path.abspath(csv_path))
    by = {lab: [r for r in rows if r["series"] == lab]
          for lab in core.SERIES_ORDER}
    checks = []

    check(checks, "data.n_rows", len(rows), 142)
    check(checks, "meta.site", meta["site"], core.SITE)
    check(checks, "meta.states", {s: meta["states"][s]["window"]
                                  for s in core.STATES},
          {s: core.STATE_WINDOWS[s] for s in core.STATES})

    # The state must be reproducible from the overpass time, not from the date.
    bad_state = [(r["series"], r["date"], r["orbit"], r["state"],
                  core.state_of(r["date"], r["orbit"]))
                 for r in rows
                 if core.state_of(r["date"], r["orbit"]) != r["state"]]
    check(checks, "state_rule.n_bad", len(bad_state), 0)

    for lab in core.SERIES_ORDER:
        sub = by[lab]
        m = meta["series"][lab]
        check(checks, f"meta.series.{lab}.n_dates", len(sub), m["n_dates"])
        check(checks, f"meta.series.{lab}.n_by_state",
              {s: sum(1 for r in sub if r["state"] == s) for s in core.STATES},
              m["n_by_state"])
        check(checks, f"meta.series.{lab}.dates_by_state",
              {s: [r["date"] for r in sub if r["state"] == s]
               for s in core.STATES},
              m["dates_by_state"])
        check(checks, f"dates_unique.{lab}",
              len({r["date"] for r in sub}), len(sub))
        check(checks, f"A_equals_n_masked.{lab}",
              sum(1 for r in sub if r["A"] != r["n_masked"]), 0)
        fname = core.SERIES[lab]
        check(checks, f"rect_sha256.{fname}",
              core.sha256(os.path.join(data_dir, fname)),
              meta["rect_file_sha256"][fname])

    for dim in CHANNELS:
        for geom in ORBITS:
            lab = core.BRIDGE_BY_GEOM[geom]
            n = sum(1 for r in by[lab] if finite(r[dim]))
            check(checks, f"channel_present.{dim}.{geom}", n > 0, True)

    # The report draws only the bridge series and states that they are complete;
    # pin that statement too rather than letting it drift with the table.
    check(checks, "bridge_complete.n_missing",
          sum(1 for geom in ORBITS for r in by[core.BRIDGE_BY_GEOM[geom]]
              for dim in CHANNELS if not finite(r[dim])), 0)

    # The plotted payload: one point per overpass, in committed file order.
    series = {}
    for geom in ORBITS:
        lab = core.BRIDGE_BY_GEOM[geom]
        points = []
        for r in by[lab]:
            acq = core.acq_dt(r["date"], geom)
            points.append({
                "date": r["date"],
                "acq_local": acq.strftime("%Y-%m-%dT%H:%M"),
                "state": r["state"],
                "values": {dim: r[dim] for dim in CHANNELS},
            })
        series[geom] = {
            "label": lab,
            "n_dates": len(points),
            "n_by_state": {s: sum(1 for p in points if p["state"] == s)
                           for s in core.STATES},
            "points": points,
        }

    # Descriptive per-state summaries; these are reported, not tested.
    summary = {}
    for geom in ORBITS:
        summary[geom] = {}
        for dim in CHANNELS:
            summary[geom][dim] = {}
            for state in core.STATES:
                vals = [p["values"][dim] for p in series[geom]["points"]
                        if p["state"] == state and finite(p["values"][dim])]
                v = np.asarray(vals, dtype=float)
                summary[geom][dim][state] = {
                    "n": int(v.size),
                    "median": float(np.median(v)) if v.size else None,
                    "mean": float(np.mean(v)) if v.size else None,
                    "min": float(np.min(v)) if v.size else None,
                    "max": float(np.max(v)) if v.size else None,
                }

    data = {
        "site": core.SITE,
        "csv": os.path.relpath(csv_path),
        "csv_meta": os.path.relpath(meta_path),
        "csv_sha256": core.sha256(csv_path),
        "reference": meta.get("reference", "neighbour"),
        "n_rows": len(rows),
        "n_by_series": {lab: len(by[lab]) for lab in core.SERIES_ORDER},
        "n_by_state": {s: sum(1 for r in rows if r["state"] == s)
                       for s in core.STATES},
        "n_by_series_state": {
            lab: {s: sum(1 for r in by[lab] if r["state"] == s)
                  for s in core.STATES}
            for lab in core.SERIES_ORDER},
        "rect_files": {lab: core.SERIES[lab] for lab in core.SERIES_ORDER},
        "rect_file_sha256": {
            core.SERIES[lab]: core.sha256(os.path.join(data_dir,
                                                       core.SERIES[lab]))
            for lab in core.SERIES_ORDER},
        "states": {s: core.STATE_WINDOWS[s] for s in core.STATES},
        "cuts": {"DS1": core.CUT_DS1.isoformat(),
                 "DS2": core.CUT_DS2.isoformat()},
        "overpass_local": {g: list(core.OVERPASS_LOCAL[g]) for g in ORBITS},
        "channels": CHANNELS,
        "definition": (
            "gamma2 = floor-corrected deck-band coherence; P = mean pixel "
            "frequency across the dates of the same series; "
            "D = A / bbox_area; F = number of 8-connected mask components; "
            "A = deck-edge pixel count; "
            "S = distance mask centroid <-> deck-band peak"),
        "state_rule": (
            "RS / DS1 / DS2 assigned by overpass time (DESC 07:08, ASC 18:51 "
            "local), not by date; cuts at 2025-05-12 10:00 and 2025-09-29 "
            "10:00 local"),
    }
    return {
        "site": core.SITE,
        "data": data,
        "series": series,
        "summary": summary,
        "cross_checks": checks,
        "n_cross_checks": len(checks),
        "n_cross_checks_failed": sum(1 for c in checks if not c["ok"]),
    }


# ---------------------------------------------------------------------------
# Figure: six dimension rows x two orbit columns
# ---------------------------------------------------------------------------
def _times(points):
    return [dt.datetime.strptime(p["acq_local"], "%Y-%m-%dT%H:%M")
            for p in points]


def _state_bounds(x0, x1):
    bounds = [x0, core.CUT_DS1, core.CUT_DS2, x1]
    return list(zip(core.STATES, bounds[:-1], bounds[1:]))


def _shade_states(ax, x0, x1):
    """Draw the RS / DS1 / DS2 background bands and the two cut lines."""
    for state, lo, hi in _state_bounds(x0, x1):
        lo, hi = min(lo, hi), max(lo, hi)
        if hi > lo:
            ax.axvspan(lo, hi, color=STATE_FILL[state], alpha=0.9, lw=0,
                       zorder=0)
    for cut in (core.CUT_DS1, core.CUT_DS2):
        ax.axvline(cut, color="0.45", ls="--", lw=0.7, zorder=1)


def _label_states(ax, x0, x1):
    """Put the three state names near the top of the first panel."""
    for state, lo, hi in _state_bounds(x0, x1):
        if hi <= lo:
            continue
        ax.text(lo + (hi - lo) / 2.0, 0.965, state,
                transform=ax.get_xaxis_transform(), ha="center", va="top",
                fontsize=6.2, color=STATE_COLORS[state], fontweight="bold",
                zorder=4)


def _channel_ylim(dim, series):
    """Row y limits: common to ASC/DESC, dimensional, data-driven."""
    vals = [float(p["values"][dim])
            for geom in ORBITS for p in series[geom]["points"]
            if finite(p["values"][dim])]
    if dim == "D":
        return 0.0, 1.0
    if dim == "P":
        return 0.0, 1.02
    if dim == "F":
        return 0.0, max(1.0, max(vals) + 0.5)
    if dim in ("A", "S"):
        hi = max(vals) if vals else 0.0
        return 0.0, hi * 1.12 if hi > 0 else 1.0
    lo, hi = min(vals), max(vals)
    lo, hi = min(0.0, lo), max(0.0, hi)
    pad = (hi - lo) * 0.08 or 0.05
    return lo - pad, hi + pad


def _draw_panel(ax, points, dim, geom):
    """One dimension/orbit panel: line plus state-coloured overpass markers."""
    ts = _times(points)
    ys = [np.nan if not finite(p["values"][dim])
          else float(p["values"][dim]) for p in points]
    ax.plot(ts, ys, color=ORBIT_COLOR[geom], lw=0.85, alpha=0.55,
            zorder=2)
    for state in core.STATES:
        keep = [(t, float(p["values"][dim]))
                for t, p in zip(ts, points)
                if p["state"] == state and finite(p["values"][dim])]
        if keep:
            ax.plot([t for t, _v in keep], [v for _t, v in keep],
                    linestyle="none", marker=STATE_MARKERS[state], ms=4.0,
                    mfc=STATE_COLORS[state], mec="white", mew=0.45,
                    zorder=3)


def make_figure(res, figdir):
    """Draw the 6 x 2 time-series grid and return the JPEG path."""
    series = res["series"]
    all_t = [t for geom in ORBITS for t in _times(series[geom]["points"])]
    x0 = min(all_t) - dt.timedelta(days=10)
    x1 = max(all_t) + dt.timedelta(days=10)

    nrows, ncols = len(CHANNELS), len(ORBITS)
    fig, axes = plt.subplots(
        nrows, ncols, figsize=(COL_WIDTH, 245 * MM),
        sharex="all", sharey="row", squeeze=False)

    for row, dim in enumerate(CHANNELS):
        ylim = _channel_ylim(dim, series)
        for col, geom in enumerate(ORBITS):
            ax = axes[row][col]
            _shade_states(ax, x0, x1)
            ax.set_xlim(x0, x1)
            ax.set_ylim(*ylim)
            if dim == "gamma2":
                ax.axhline(0.0, color="0.55", lw=0.7, ls=":", zorder=1)
            _draw_panel(ax, series[geom]["points"], dim, geom)

            # One y label per row, one x label per column, panel letters in
            # the corner; with sharex/sharey the inner tick labels are hidden.
            if col == 0:
                ax.set_ylabel(FIG_CHANNEL_AXIS[dim])
            ax.tick_params(labelleft=(col == 0),
                           labelbottom=(row == nrows - 1))
            if dim == "F":
                ax.yaxis.set_major_locator(MaxNLocator(integer=True))
            ax.text(0.012, 0.955, f"({PANEL_LETTERS[row * ncols + col]})",
                    transform=ax.transAxes, ha="left", va="top",
                    fontsize=6.8, color="0.25", zorder=5)
            if row == 0:
                ax.set_title(ORBIT_LABEL[geom], fontsize=8.0, pad=5)
            if row == 0 and col == 0:
                _label_states(ax, x0, x1)

    for ax in axes[-1]:
        ax.xaxis.set_major_locator(mdates.MonthLocator(interval=1))
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%b"))
        ax.set_xlabel("2025 — local overpass time")

    handles = [Patch(facecolor=STATE_FILL[s], edgecolor=STATE_COLORS[s],
                     label=STATE_LABEL[s]) for s in core.STATES]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False,
               bbox_to_anchor=(0.5, 0.004))
    fig.suptitle(
        "Bautzen OpenLabs bridge — OCV time series by dimension and orbit\n"
        "bridge deck edge; markers = overpasses (shape and colour by state), "
        "line = acquisition order\n"
        "state bands RS / DS1 / DS2; dashed cuts 2025-05-12 and 2025-09-29, "
        "10:00 local",
        fontsize=8.3, y=0.985)
    fig.tight_layout(rect=(0, 0.035, 1, 0.925))
    os.makedirs(figdir, exist_ok=True)
    out = os.path.join(figdir, FIG_NAME + ".jpg")
    fig.savefig(out, dpi=FIG_DPI, facecolor="white")
    plt.close(fig)
    return out


# ---------------------------------------------------------------------------
# Markdown report
# ---------------------------------------------------------------------------
def report(res):
    d = res["data"]
    L = []
    A = L.append
    A("# Bautzen OCV time series by dimension and orbit")
    A("")
    A("## 1. Question and layout")
    A("")
    A("How do the six OCV components of the Bautzen OpenLabs bridge deck edge "
      "move through the three documented states RS, DS1 and DS2, and do the "
      "two Sentinel-1 pass geometries tell the same story? The figure draws "
      "one panel per **dimension x orbit**: six rows in the formal OCV order "
      "(`gamma2`, `P`, `D`, `F`, `A`, `S`) and two columns (ASC IW2, DESC IW3), "
      "i.e. twelve panels.")
    A("")
    A("The x coordinate is the **actual local overpass instant**, not "
      "midnight. The state boundaries fall in the middle of a day (DESC "
      "07:08, ASC 18:51 local), so a date-only axis would place the "
      "boundary-day points in the wrong state band. The three state windows "
      "are background bands with different colours and the two cut dates are "
      "dashed vertical lines.")
    A("")
    A("Only the bridge series is drawn. The control rectangle 500 m east is "
      "the counter-check of the paper figure, but its ASC `gamma2` range "
      "(-0.94 .. +0.93) is much wider than the bridge range (-0.14 .. +0.24); "
      "a shared panel would compress the bridge signal this figure is meant "
      "to show. The control remains in `data/bautzen/bautzen_ocv_channels.csv` "
      "and in figures A-E.")
    A("")
    A("## 2. Data and state model")
    A("")
    A(f"- source: `{d['csv']}`")
    A(f"- committed metadata: `{d['csv_meta']}`")
    A(f"- overpasses: **{d['n_rows']}** in four rect series "
      f"({', '.join(f'{lab}: {n}' for lab, n in d['n_by_series'].items())})")
    A(f"- states: **{d['n_by_state']['RS']} RS**, "
      f"**{d['n_by_state']['DS1']} DS1**, **{d['n_by_state']['DS2']} DS2**")
    A(f"- cuts: DS1 `{d['cuts']['DS1']}`, DS2 `{d['cuts']['DS2']}` local")
    A("")
    A("State assignment is by overpass time, not date:")
    A("")
    for s in core.STATES:
        A(f"- **{s}**: {d['states'][s]}")
    A("")
    A("The bridge series used by the figure:")
    A("")
    rows = []
    for geom in ORBITS:
        s = res["series"][geom]
        rows.append([geom, s["label"], s["n_dates"],
                     s["n_by_state"]["RS"], s["n_by_state"]["DS1"],
                     s["n_by_state"]["DS2"]])
    A(md_table(["orbit", "series", "dates", "RS", "DS1", "DS2"], rows))
    A("")
    A("## 3. Per-state medians")
    A("")
    A("Median of each channel inside each state, with the number of finite "
      "overpasses in parentheses. The bridge panels are complete (all 36 ASC "
      "and 35 DESC overpasses); the plotting path would leave any missing "
      "value as a gap, and nothing is interpolated.")
    A("")
    rows = []
    for geom in ORBITS:
        for dim in CHANNELS:
            row = [geom, dim]
            for state in core.STATES:
                c = res["summary"][geom][dim][state]
                row.append(f"{fmt(c['median'])} (n={c['n']})")
            rows.append(row)
    A(md_table(["orbit", "dimension", "RS median (n)", "DS1 median (n)",
                "DS2 median (n)"], rows))
    A("")
    A("The corresponding means, minima and maxima are in the result JSON "
      "under `summary.<orbit>.<dimension>.<state>`.")
    A("")
    A("## 4. Cross-checks")
    A("")
    A(f"The script recomputes every plotted point from the committed channel "
      f"table and checks it against the committed metadata: "
      f"**{res['n_cross_checks']} checks, "
      f"{res['n_cross_checks_failed']} failures**. The checks cover the row "
      "count, the state rule, the per-series date/state counts, uniqueness of "
      "dates, `A == n_masked`, the presence and completeness of every channel "
      "on both bridge orbits, and the sha256 of every rect file.")
    A("")
    if res["n_cross_checks_failed"]:
        A("**Failures:**")
        A("")
        A(md_table(["check", "recomputed", "committed"],
                   [[c["check"], c.get("recomputed", ""),
                     c.get("committed", "")]
                    for c in res["cross_checks"] if not c["ok"]]))
    else:
        A("**Failures: none.** Every state value reproduces `core.state_of`, "
          "every per-series count and date list reproduces the committed "
          "metadata, and every rect-file sha256 matches.")
    A("")
    A("## 5. Reproduce")
    A("")
    A("```sh")
    A("cd <repo root>")
    A("python3 code/bautzen/fig_bautzen_ocv_orbits.py")
    A("PYTHON=/usr/bin/python3 bash figures/bautzen/fig_bautzen_ocv_orbits.sh")
    A("```")
    A("")
    A("Needs numpy and matplotlib. This report was generated with python "
      f"{platform.python_version()}, numpy {np.__version__}, matplotlib "
      f"{matplotlib.__version__}. No output carries a timestamp, so repeated "
      "runs are identical.")
    A("")
    A("## 6. Caveats")
    A("")
    A("- **States are consecutive time intervals, not a randomized "
      "treatment.** RS, DS1 and DS2 are one deployment observed over six "
      "months; state and season are confounded. The background bands make "
      "that structure visible, they do not remove it.")
    A("- **No control series in this figure.** The control rectangle is the "
      "counter-check in figures A-E; it is omitted here only because its "
      "ASC `gamma2` range would flatten the bridge signal in a shared panel.")
    A("- **Missing values are gaps.** The bridge series are complete: all 36 "
      "ASC and 35 DESC overpasses carry all six channels. The 14 missing "
      "`gamma2` rows and the 8 missing `D`, `S`, `P` rows are in the control "
      "series, which this figure does not draw; no value is interpolated.")
    A("- **`gamma2` is the floor-corrected deck-band coherence** "
      "(`gamma2_band_floor_corrected`), the same channel as paper figure A, "
      "not the raw estimator.")
    A("- **The x axis is local overpass time.** ASC is drawn at 18:51 and DESC "
      "at 07:08 local; this is what makes the state assignment at the two cut "
      "days correct.")
    A("")
    A("---")
    A("")
    A("Generated by `code/bautzen/fig_bautzen_ocv_orbits.py` "
      f"({res['n_cross_checks']} cross-checks, "
      f"{res['n_cross_checks_failed']} failures).")
    return "\n".join(L) + "\n"


# ---------------------------------------------------------------------------
# Command line
# ---------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Bautzen OCV time series per dimension and orbit.")
    ap.add_argument("--figdir", default=FIGDIR,
                    help="output directory of the JPEG figure")
    ap.add_argument("--json", default=OUT_JSON,
                    help="path of the result JSON")
    ap.add_argument("--md", default=OUT_MD,
                    help="path of the Markdown report")
    ap.add_argument("--no-figures", action="store_true",
                    help="only recompute and cross-check (no matplotlib output)")
    ap.add_argument("--csv", default=None,
                    help="channel table to plot (default: the committed "
                         "bautzen_ocv_channels.csv)")
    ap.add_argument("--meta", default=None,
                    help="its metadata JSON (default: the committed one); the rect "
                         "sha pins are read from the CSV's own directory")
    args = ap.parse_args(argv)

    t0 = time.time()
    print("Bautzen OCV time series by dimension and orbit ...")
    res = compute(args.csv, args.meta)
    print(f"  compute: {time.time() - t0:.1f}s "
          f"({res['data']['n_rows']} overpasses, "
          f"{res['n_cross_checks']} cross-checks, "
          f"{res['n_cross_checks_failed']} failed)")
    for c in res["cross_checks"]:
        if not c["ok"]:
            print(f"    FAIL {c['check']}")

    if not args.no_figures:
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
