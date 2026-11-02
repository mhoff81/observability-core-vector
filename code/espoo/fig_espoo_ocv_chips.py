#!/usr/bin/env python3
"""Espoo on the **chip** — the OCV vector of the 75 support chips.

This is a second instrument beside the site's pipeline figure (retired; its
frozen lock survives as ``data/espoo/fig_espoo_ocv_orbits.json`` and its export
as ``data/espoo/espoo_channels.csv``), and the two are not interchangeable:

* the pipeline figure reads the committed export (``coherence_gamma2``,
  ``coherence_masked_pixels``, ...) over 150 overpasses;
* this one reads ``data/espoo/support/espoo_ocv_channels.csv`` -- 75 chips cut
  with the reference rule (support 12 m base + height layover, guard 30 m, ring
  30-200 m, ``ring_complete: true``) -- and carries ``A, D, F, S, P`` and
  ``gamma2_mask`` for the first time at this site.

They agree only loosely: r(log A, log pipeline mask size) = -0.24 and
r(gamma2_mask, coherence_gamma2) = 0.13. Neither is wrong; they measure
different things.

What the chip series is, structurally
-------------------------------------
* **75 of 150 overpasses.** The mast sits on the far range edge of IW1 in two of
the four relative orbits, so the cut is restricted to the two tracks that place
  it inside the swath: **R160 (ASC, 31 dates)** and **R51 (DESC, 44 dates)**.
  The two tracks have different dates -- only **3 days are shared** -- so the
  orbit contrast is **unpaired** here, unlike the 69 same-day pairs of the
  pipeline series.
* **`A` is orbit-dependent by a factor of seven** (ASC median 76, DESC 528).
* **Two dimensions are functions of `A`**: ``D = A / bbox_area`` exactly, and
  ``gamma2_mask`` carries ``n = A`` in its denominator. A drift in ``D`` or
  ``gamma2_mask`` is therefore partly a drift in the mask size, not an
  independent observation.

The tests are the same three as in the pipeline figure -- drift, year repeat
(two 365-day campaign windows), orbit contrast -- with the pairing dropped and
that stated in the report.

Usage
-----
  python3 code/espoo/fig_espoo_ocv_chips.py
  python3 code/espoo/fig_espoo_ocv_chips.py --selftest
  python3 code/espoo/fig_espoo_ocv_chips.py --verify

Outputs
-------
  figures/espoo/fig_espoo_ocv_chips.jpg   (630 dpi JPEG, the 170 mm canvas)
  figures/espoo/fig_espoo_ocv_chips.pdf   (vector, Type-42 text: the LaTeX copy)
  figures/espoo/fig_espoo_ocv_chips.eps   (vector, alpha pre-blended over white)
  figures/espoo/fig_espoo_ocv_chips.md    (the caption and the figure audit)
  data/espoo/fig_espoo_ocv_chips.json

The plate
---------
One authored canvas of 170 mm (:data:`figure_style.COL_WIDTH`) and three rows of
two panels, exported **whole**.  ``bbox_inches="tight"`` would make the file as
wide as the layout happens to be and shrink every glyph with it, so the bounds
are *asserted* instead of cropped: the canvas is the print width, and the drawn
content has to stay inside it with a visible margin
(:func:`figure_style.assert_fits`).  The audit that records this -- target and
native width, ``figsize``, dpi, the final point sizes, the EPS alpha policy and
one pinned row per panel and track -- is written into the ``.md`` beside the
caption, so the scientific JSON stays exactly what it was.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import math
import os
import shutil
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
for path in (HERE, os.path.join(REPO, "code"), os.path.join(REPO, "code", "lumo")):
    if path not in sys.path:
        sys.path.insert(0, path)

import espoo_ocv_stats as st  # noqa: E402
import figure_style as style  # noqa: E402
from lumo_ocv_stats import cliffs_delta  # noqa: E402

SCHEMA = "espoo-ocv-chips-v1"
SITE = "espoo"
LABEL = "Espoo mast (OSM 2233589804) — support chips"
TABLE = os.path.join(REPO, "data", SITE, "support", f"{SITE}_ocv_channels.csv")
#: 630 dpi JPEG, like the other companion plates of this repository
FIG_DPI = 630
OUT_FIG = os.path.join(REPO, "figures", SITE, "fig_espoo_ocv_chips.jpg")
OUT_MD = os.path.join(REPO, "figures", SITE, "fig_espoo_ocv_chips.md")
OUT_JSON = os.path.join(REPO, "data", SITE, "fig_espoo_ocv_chips.json")
OUT_PDF = os.path.splitext(OUT_FIG)[0] + ".pdf"
OUT_EPS = os.path.splitext(OUT_FIG)[0] + ".eps"

DIMS = {
    "A": ("A", "masked cells (echo area)", False),
    "D": ("D", "density (A / bbox)", False),
    "F": ("F", "fragmentation (components)", False),
    "S": ("S", "shift (centroid to peak)", False),
    "P": ("P", "persistence (pixel frequency)", False),
    "gamma2_mask": ("gamma2_mask", "echo coherence", False),
}
ORBITS = (("ASCENDING", "ASC R160", "#2b6cb0"),
          ("DESCENDING", "DESC R51", "#c53030"))
BLOCK = 5
BOOT = 2000
SEED = 13


# ---------------------------------------------------------------------------
# The plate: one 170 mm canvas, authored and asserted, never cropped
# ---------------------------------------------------------------------------
#: Canvas in inches.  ``style.COL_WIDTH`` is exactly 170 mm; the height is what
#: three rows of two panels need once the statistics band, the date labels and
#: the one-line header are paid for -- six panels of 74.6 x 45.7 mm.
PLATE_IN = (style.COL_WIDTH, 7.55)
PLATE_DPI = FIG_DPI
#: Grid fractions.  ``wspace`` is set by the widest y tick label of the right
#: column ("0.175"): it has to clear the left panel by >= 1.5 mm.  ``hspace``
#: is set by the two statistics rows plus the date labels of the row below.
PLATE_GRID = dict(left=0.055, right=0.994, top=0.915, bottom=0.048,
                  wspace=0.15, hspace=0.336)
#: The one-line header.  Counts, dates, orbit ids and the trend semantics live
#: in the caption: the plate states what it shows, nothing about the pipeline.
HEADER = "Espoo mast — OCV trajectories by acquisition track"
#: How each dimension is labelled on the plate: canonical term, unit in
#: parentheses, and no column name -- ``gamma2_mask`` is the coherence
#: ``gamma**2``, in math text, because that is what it measures.
PANEL_LABEL = {
    "A": "A — observable area (cells)",
    "D": "D — density (fraction of bbox)",
    "F": "F — fragmentation (components)",
    "S": "S — spatial shift (cells)",
    "P": "P — persistence (pixel frequency)",
    "gamma2_mask": r"$\gamma^2$ — coherence",
}
#: Reading order of the panels -- pinned by the audit and by the caption.
PANEL_ORDER = tuple(DIMS)
#: One fixed statistics band above every axes: DESC above, ASC second,
#: flush with the left spine, identical in all six panels.  Points from the
#: axes top; ``TITLE_PAD`` leaves ~4 mm between the DESC row and the title.
STAT_PT = 6.6
STAT_ROW = {"DESCENDING": 10.0, "ASCENDING": 1.0}
TITLE_PAD = 20.0
#: Marks, sized for the final 170 mm plate: a white edge keeps the outline of a
#: 3.6 pt dot visible, the trend line is heavier than the grid but lighter than
#: the markers (the markers are drawn on top of it).
MARKER_PT = 3.6
MARKER_EDGE_PT = 0.5
TREND_LW = 1.2
GRID_ALPHA = 0.18
GRID_LW = 0.5
#: The drawn content must stay this far inside the canvas: the raster is
#: written from the canvas, so a glyph at the edge would be clipped -- a tight
#: bbox would have hidden that by cropping, which is the defect it must not be.
MARGIN_MIN_MM = 0.3


# ---------------------------------------------------------------------------
# Data and the three tests
# ---------------------------------------------------------------------------
def _number(row, column):
    raw = (row.get(column) or "").strip()
    if raw in ("", "NA", "null", "NaN"):
        return None
    try:
        return float(raw)
    except ValueError:
        return None


def load(table=None):
    """One dict per chip: date, orbit, and the six chip dimensions."""
    table = table or TABLE
    rows = []
    with open(table, newline="") as handle:
        for raw in csv.DictReader(handle):
            stamp = (raw.get("acquisition_date") or raw.get("acquisition_ts") or "")[:10]
            when = dt.date.fromisoformat(stamp)
            rows.append({"date": when, "ordinal": when.toordinal(),
                         "year": when.year, "month": when.month,
                         "orbit": raw["orbit_direction"],
                         "values": {name: _number(raw, column)
                                    for name, (column, _l, _s) in DIMS.items()}})
    rows.sort(key=lambda row: (row["ordinal"], row["orbit"]))
    return rows


def series_for(rows, orbit, dim):
    out = [[row["ordinal"], row["values"][dim]] for row in rows
           if row["orbit"] == orbit and row["values"][dim] is not None]
    out.sort()
    return out


def drift(rows, boot=BOOT):
    """Test 1: the ported stationarity block plus a ties diagnostic."""
    out = {}
    for dim in DIMS:
        out[dim] = {}
        for orbit, name, _colour in ORBITS:
            series = series_for(rows, orbit, dim)
            if len(series) < 12:
                out[dim][orbit] = {"available": False, "n": len(series)}
                continue
            block = st.stationarity_of(series, block=BLOCK, n_boot=boot)
            block["available"] = True
            block["orbit_label"] = name
            rounded = [round(value, 9) for _day, value in series]
            counted = {}
            for value in rounded:
                counted[value] = counted.get(value, 0) + 1
            top = max(counted, key=counted.get)
            block["top_value"] = top
            block["top_value_share"] = float(counted[top] / len(rounded))
            block["degenerate"] = bool(
                block["top_value_share"] >= 0.5
                and not block["theil_sen_per_year"]
                and block["slope_ci95_per_year"] in ([0.0, 0.0], None))
            out[dim][orbit] = block
    return out


def bootstrap_median_ci(values, block=BLOCK, n=BOOT, seed=SEED):
    array = np.asarray(values, float)
    if array.size < 4:
        return None
    rng = np.random.default_rng(seed)
    block = max(1, min(block, array.size))
    draws = np.empty(n)
    for draw in range(n):
        starts = rng.integers(0, array.size - block + 1,
                              int(np.ceil(array.size / block)))
        index = np.concatenate([np.arange(s, s + block) for s in starts])[:array.size]
        draws[draw] = float(np.median(array[index]))
    return [float(np.percentile(draws, 2.5)), float(np.percentile(draws, 97.5))]


def year_repeat(rows, boot=BOOT):
    """Test 2: the same calendar month in the two 365-day campaign windows."""
    start = min(row["date"] for row in rows)
    out = {"by_dim": {}}
    if (max(row["date"] for row in rows) - start).days < 700:
        return out
    for row in rows:
        row["campaign_year"] = 0 if (row["date"] - start).days < 365 else 1
    for dim in DIMS:
        out["by_dim"][dim] = {}
        for orbit, name, _colour in ORBITS:
            pairs = []
            for month in range(1, 13):
                left = [row["values"][dim] for row in rows
                        if row["orbit"] == orbit and row["campaign_year"] == 0
                        and row["month"] == month and row["values"][dim] is not None]
                right = [row["values"][dim] for row in rows
                         if row["orbit"] == orbit and row["campaign_year"] == 1
                         and row["month"] == month and row["values"][dim] is not None]
                if left and right:
                    pairs.append(float(np.median(right) - np.median(left)))
            if len(pairs) < 6:
                out["by_dim"][dim][orbit] = {"available": False,
                                             "months": len(pairs)}
                continue
            positive = sum(1 for value in pairs if value > 0)
            out["by_dim"][dim][orbit] = {
                "available": True, "orbit_label": name, "months": len(pairs),
                "median_delta": float(np.median(pairs)),
                "ci": bootstrap_median_ci(pairs, n=boot),
                "months_positive": int(positive),
                "months_negative": int(len(pairs) - positive),
            }
    return out


def orbit_contrast(rows, boot=BOOT):
    """Test 3, **unpaired**: ASC against DESC as two groups, not by date.

    The track filter that makes the chip possible puts the two orbits on
    different dates (only three days are shared), so there is no same-day pair
    to difference. The contrast is therefore a two-sample comparison, and the
    report says so instead of pretending otherwise.
    """
    out = {}
    for dim in DIMS:
        left = [row["values"][dim] for row in rows
                if row["orbit"] == "ASCENDING" and row["values"][dim] is not None]
        right = [row["values"][dim] for row in rows
                 if row["orbit"] == "DESCENDING" and row["values"][dim] is not None]
        if len(left) < 8 or len(right) < 8:
            out[dim] = {"available": False, "n": [len(left), len(right)]}
            continue
        # ``cliffs_delta(a, b)`` is positive when b > a; pass DESC as the
        # reference so that a positive delta means "ASC above DESC".
        delta = cliffs_delta(right, left)
        rng = np.random.default_rng(SEED)
        draws = np.empty(boot)
        for draw in range(boot):
            a = np.asarray(left)[rng.integers(0, len(left), len(left))]
            b = np.asarray(right)[rng.integers(0, len(right), len(right))]
            draws[draw] = cliffs_delta(b, a)
        test = st.compare_orbit(left, right) or {}
        out[dim] = {
            "available": True, "n": [len(left), len(right)],
            "median": [float(np.median(left)), float(np.median(right))],
            "delta": delta,
            "delta_ci95": [float(np.percentile(draws, 2.5)),
                           float(np.percentile(draws, 97.5))],
            "mannwhitney_p": test.get("mannwhitney_p"),
            "welch_p": test.get("welch_p"),
            "paired": False,
            "shared_dates": None,
        }
    shared = ({row["date"] for row in rows if row["orbit"] == "ASCENDING"}
              & {row["date"] for row in rows if row["orbit"] == "DESCENDING"})
    for dim in out:
        out[dim]["shared_dates"] = len(shared)
    return out


def analyse(table=None, boot=BOOT):
    rows = load(table)
    table = table or TABLE
    import hashlib
    digest = hashlib.sha256(open(table, "rb").read()).hexdigest()
    by_orbit = {orbit: sum(1 for row in rows if row["orbit"] == orbit)
                for orbit, _name, _colour in ORBITS}
    shared = len({row["date"] for row in rows if row["orbit"] == "ASCENDING"}
                 & {row["date"] for row in rows if row["orbit"] == "DESCENDING"})
    return {
        "schema": SCHEMA, "site": SITE, "label": LABEL,
        "table": {"path": os.path.relpath(table, REPO), "sha256": digest,
                  "n_rows": len(rows)},
        "by_orbit": by_orbit, "shared_dates": shared,
        "dates": [min(row["date"] for row in rows).isoformat(),
                  max(row["date"] for row in rows).isoformat()],
        "dimensions": {name: {"column": column, "label": label}
                       for name, (column, label, _s) in DIMS.items()},
        "drift": drift(rows, boot=boot),
        "year_repeat": year_repeat(rows, boot=boot),
        "orbit_contrast": orbit_contrast(rows, boot=boot),
        "drawing_rule": "the Theil-Sen line is solid only where the ported "
                        "stationarity flag says trend, dashed otherwise",
        "rules": {"block": BLOCK, "boot": int(boot), "seed": SEED,
                  "states": "none: the record carries no damage label",
                  "tracks": "ASC = relative orbit 160, DESC = 51; the other "
                            "two tracks put the mast on the far range edge",
                  "pairing": "none: the tracks have 3 shared dates",
                  "tautology": "D = A / bbox_area exactly; gamma2_mask "
                               "contains n = A in its denominator",
                  "instrument": "chip values; not the pipeline columns "
                                "(r(log A, log coherence_masked_pixels) = -0.24)"},
    }


# ---------------------------------------------------------------------------
# The figure
# ---------------------------------------------------------------------------
def statistics_line(block):
    """One statistics row: slope, interval, MK p, verdict -- nothing else.

    The words that would not fit beside a 74 mm panel at 6.6 pt ("slope",
    "95% interval", "Mann-Kendall") are in the caption; the row itself is one
    line of numbers, formatted identically in all six panels and both tracks.
    ``trend`` is the ported flag rule firing, ``flat`` is it not firing -- the
    same word the console report uses.
    """
    slope = block.get("theil_sen_per_year")
    ci = block.get("slope_ci95_per_year") or [float("nan"), float("nan")]
    p = (block.get("mann_kendall") or {}).get("p")
    flagged = "trend" in (block.get("flags") or [])
    track = (block.get("orbit_label") or "").split()
    return (f"{(track[0] if track else '---')} {slope:+.3g}/a "
            f"[{ci[0]:+.3g}, {ci[1]:+.3g}] p {_p(p)} "
            f"{'trend' if flagged else 'flat'}")


def trend_endpoints(block, series):
    """The Theil-Sen line exactly as it is drawn: first date to last date."""
    x = np.array([day for day, _value in series], float)
    y = np.array([value for _day, value in series], float)
    per_day = block["theil_sen_per_year"] / 365.25
    intercept = float(np.median(y - per_day * x))
    xs = np.array([x.min(), x.max()])
    return xs, intercept + per_day * xs


def make_figure(report, rows, fig_path=None):
    """Draw the 170 mm plate, write jpg + pdf + eps; return ``(audit, drawn)``.

    Nothing here decides a scientific value: every number comes from *report*,
    every point from *rows*, and the two are cross-checked before the file is
    written (:func:`assert_figure_matches_report`).
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.dates as mdates
    plt = style.apply(dpi=PLATE_DPI)
    from matplotlib.transforms import offset_copy

    fig_path = fig_path or OUT_FIG
    os.makedirs(os.path.dirname(os.path.abspath(fig_path)), exist_ok=True)
    # The canvas is authored, not measured from the content: 170 mm is the
    # print width, and the six panels are 74.6 x 45.7 mm inside it.  A tight
    # bbox would re-cut this after the fact and shrink every glyph with it.
    fig, grid = plt.subplots(3, 2, figsize=PLATE_IN)
    fig.subplots_adjust(**PLATE_GRID)
    axes = list(grid.ravel())
    stats_rows = {}          # per axes, for the clearance guards below
    titles = {}              # the left-aligned panel title, for the same reason
    drawn = {"panel_order": list(PANEL_ORDER), "panels": {}}
    for index, dim in enumerate(PANEL_ORDER):
        ax = axes[index]
        ax.set_axisbelow(True)
        ax.grid(alpha=GRID_ALPHA, lw=GRID_LW)
        tracks = {}
        for orbit, name, colour in ORBITS:
            series = series_for(rows, orbit, dim)
            block = report["drift"][dim].get(orbit) or {}
            if not series or not block.get("available"):
                continue
            dates = [dt.date.fromordinal(int(day)) for day, _value in series]
            values = [value for _day, value in series]
            # Markers above the trend line: the observations are the data, the
            # line is the summary, so the dots must not be crossed by it.
            ax.plot(dates, values, linestyle="none", marker="o", ms=MARKER_PT,
                    mfc=colour, mec="white", mew=MARKER_EDGE_PT, alpha=0.9,
                    zorder=4, label=f"{name} (n={len(values)})")
            xs, ys = trend_endpoints(block, series)
            flagged = "trend" in (block.get("flags") or [])
            ax.plot([dt.date.fromordinal(int(v)) for v in xs], ys, color=colour,
                    lw=TREND_LW, ls="-" if flagged else "--", alpha=0.9, zorder=3)
            tracks[orbit] = {
                "n": len(values),
                "first": dates[0].isoformat(),
                "last": dates[-1].isoformat(),
                "median": float(np.median(values)),
                "values": [[int(day), float(value)] for day, value in series],
                "trend": {"x": [float(xs[0]), float(xs[1])],
                          "y": [float(ys[0]), float(ys[1])],
                          "flagged": flagged},
            }
        drawn["panels"][dim] = {"label": PANEL_LABEL[dim], "tracks": tracks}
        # One fixed band per panel: DESC above, ASC second, flush with the left
        # spine -- the same rule in all six, so the panels can be compared.
        stats_rows[ax] = [
            ax.text(0.0, 1.0, statistics_line(report["drift"][dim][orbit]),
                    transform=offset_copy(ax.transAxes, fig=fig,
                                          y=STAT_ROW[orbit], units="points"),
                    va="bottom", ha="left", fontsize=STAT_PT, color=colour,
                    clip_on=False)
            for orbit, _name, colour in ORBITS
            if (report["drift"][dim].get(orbit) or {}).get("available")]
        # ``set_title`` returns the artist, and a left title is not ``ax.title``:
        # the guards need the box that is actually drawn.
        titles[ax] = ax.set_title(PANEL_LABEL[dim], loc="left", pad=TITLE_PAD,
                                  fontsize=style.FS["title"])
        ax.tick_params(labelsize=style.FS["tick"], length=2.4, width=0.5, pad=2)
        ax.xaxis.set_major_locator(mdates.AutoDateLocator(minticks=4, maxticks=6))
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m"))
    for ax in axes[3:]:                 # one date label per column, at the foot
        ax.set_xlabel("acquisition date", fontsize=style.FS["label"])
    # One figure legend, frameless, at the right end of the header line: the
    # counts belong to the track, not to a panel, and six small legends were
    # six places to read the same two numbers.
    handles = [plt.Line2D([], [], marker="o", linestyle="none", ms=MARKER_PT,
                          mfc=colour, mec="white", mew=MARKER_EDGE_PT,
                          label=f"{name} (n={report['by_orbit'][orbit]})")
               for orbit, name, colour in ORBITS]
    legend = fig.legend(handles=handles, loc="upper right", ncol=2,
                        bbox_to_anchor=(PLATE_GRID["right"], 0.996),
                        frameon=False, fontsize=style.FS["legend"],
                        handletextpad=0.4, columnspacing=1.4,
                        borderaxespad=0.0, labelspacing=0.2)
    header = fig.suptitle(HEADER, x=0.006, y=0.996, ha="left", va="top",
                          fontsize=style.FS["suptitle"])
    # --- guards -------------------------------------------------------------
    # These are the checks a tight bbox hides: with the canvas as the export
    # box, a glyph at the edge is either inside it or cut off, and a band that
    # grew by half a point collides instead of pushing the frame wider.
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    to_mm = 25.4 / fig.dpi
    for ax in axes:
        box = ax.get_window_extent(renderer)
        rows_here = stats_rows.get(ax) or []
        for artist in rows_here:
            span = artist.get_window_extent(renderer)
            if span.x0 < box.x0 or span.x1 > box.x1:
                raise RuntimeError(
                    "%r: the statistics line runs %.1f mm past the axes edge"
                    % (titles[ax].get_text(),
                       max(box.x0 - span.x0, span.x1 - box.x1) * to_mm))
        if rows_here:
            top = max(a.get_window_extent(renderer).y1 for a in rows_here)
            if titles[ax].get_window_extent(renderer).y0 < top + 0.5:
                raise RuntimeError("%r: the title touches the DESC row"
                                   % titles[ax].get_text())
    header_box = header.get_window_extent(renderer)
    legend_box = legend.get_window_extent(renderer)
    gap = (legend_box.x0 - header_box.x1) * to_mm
    if gap < 2.0:
        raise RuntimeError("the header and the legend leave no gap (%.1f mm): "
                           "shorten the header" % gap)
    for name, span in (("header", header_box), ("legend", legend_box)):
        if (span.x0 < 0 or span.x1 > fig.bbox.width
                or span.y0 < 0 or span.y1 > fig.bbox.height):
            raise RuntimeError("the %s leaves the canvas" % name)
    title_top = max(titles[ax].get_window_extent(renderer).y1 for ax in axes)
    if (header_box.y0 - title_top) * to_mm < 1.5:
        raise RuntimeError("the header clears the panel titles by only %.1f mm"
                           % ((header_box.y0 - title_top) * to_mm))
    for dim in PANEL_ORDER:
        for orbit, name, _colour in ORBITS:
            track = drawn["panels"][dim]["tracks"].get(orbit)
            if track and track["n"] != report["by_orbit"][orbit]:
                raise RuntimeError(
                    "%s/%s: %d points drawn, the legend says n=%d"
                    % (dim, name, track["n"], report["by_orbit"][orbit]))
    tight = fig.get_tightbbox(renderer)                       # inches
    edges = {"left": tight.x0, "right": PLATE_IN[0] - tight.x1,
             "bottom": tight.y0, "top": PLATE_IN[1] - tight.y1}
    worst = min(edges.values()) * 25.4
    if worst < MARGIN_MIN_MM:
        raise RuntimeError(
            "the drawn content is %.2f mm from the canvas edge (min %.1f mm): "
            "re-cut the grid -- do not crop the export" % (worst, MARGIN_MIN_MM))
    assert_figure_matches_report(report, drawn)
    translucent = style.alpha_artists(fig)     # counted before they are blended
    limits = {dim: [axes[PANEL_ORDER.index(dim)].get_xlim(),
                    axes[PANEL_ORDER.index(dim)].get_ylim()]
              for dim in PANEL_ORDER}
    style.save(fig, fig_path, dpi=PLATE_DPI, eps="blend")
    plt.close(fig)
    return plate_audit(report, drawn, fig_path, tight, edges, translucent,
                       limits), drawn
def plate_audit(report, drawn, fig_path, tight, edges, translucent, limits):
    """The machine-readable record of the plate, written into the ``.md``.

    It answers what a reader of the figure cannot see by looking: how wide the
    file really is (measured from the file, all three containers), what
    happened to the alpha in the EPS, and which numbers each panel and track
    pins.  Nothing here is recomputed -- these are the values that were drawn.
    """
    stem = os.path.splitext(fig_path)[0]
    plate = {
        "target_width_mm": style.FINAL_WIDTH_MM,
        "target_height_mm": round(PLATE_IN[1] * 25.4, 2),
        "figsize_inches": [round(PLATE_IN[0], 6), round(PLATE_IN[1], 6)],
        "native_dpi": PLATE_DPI,
        #: Keyed by the ``figure_style`` role, so a check can read its band.
        "final_font_sizes_pt": {
            "suptitle": style.FS["suptitle"],
            "title": style.FS["title"],
            "label": style.FS["label"],
            "tick": style.FS["tick"],
            "legend": style.FS["legend"],
            "annotation": STAT_PT,
        },
        "native_width_mm": {},
        "native_height_mm": {},
        "latex_display_width_mm": style.FINAL_WIDTH_MM,
        "latex_scaling": "the file is 170 mm wide — included with "
                         "`width=\\linewidth`; no trim, clip or resizebox",
        "bounds": {
            "export": "canvas-authorised (figure_style.save); bbox_inches is "
                      "never 'tight'",
            "canvas_mm": [round(PLATE_IN[0] * 25.4, 2),
                          round(PLATE_IN[1] * 25.4, 2)],
            "content_mm": [round(tight.width * 25.4, 1),
                           round(tight.height * 25.4, 1)],
            "edge_mm": {key: round(value * 25.4, 2)
                        for key, value in edges.items()},
            "minimum_edge_mm": MARGIN_MIN_MM,
        },
        "eps": {"alpha_artists": len(translucent), "strategy": "blend",
                "background": "white", "rasterised_artists": 0,
                "vector_text": "Type 42 (pdf/ps fonttype 42)"},
        "panel_order": list(PANEL_ORDER),
        "statistics": {
            "pt": STAT_PT,
            "row_offset_pt": {orbit: STAT_ROW[orbit]
                              for orbit, _name, _colour in ORBITS},
            "anchor": "x=0.0 axes, ha=left, va=bottom -- one rule, six panels",
        },
        "panels": {},
    }
    for path in (fig_path, stem + ".pdf", stem + ".eps"):
        info = style.measure(path)
        plate["native_width_mm"][info["format"]] = round(info["mm"], 2)
        plate["native_height_mm"][info["format"]] = round(info["height_mm"], 2)
    for dim in PANEL_ORDER:
        x_lim, y_lim = limits[dim]
        panel = {"label": PANEL_LABEL[dim],
                 "x_limits": [round(x_lim[0], 4), round(x_lim[1], 4)],
                 "y_limits": [round(y_lim[0], 6), round(y_lim[1], 6)],
                 "tracks": {}}
        for orbit, _name, _colour in ORBITS:
            track = drawn["panels"][dim]["tracks"].get(orbit)
            block = report["drift"][dim].get(orbit) or {}
            if track is None or not block.get("available"):
                panel["tracks"][orbit] = {"n": 0, "available": False}
                continue
            panel["tracks"][orbit] = {
                "n": track["n"],
                "first": track["first"],
                "last": track["last"],
                "median": round(track["median"], 6),
                "iqr": block.get("iqr"),
                "theil_sen_per_year": block.get("theil_sen_per_year"),
                "slope_ci95_per_year": block.get("slope_ci95_per_year"),
                "mann_kendall_p": (block.get("mann_kendall") or {}).get("p"),
                "verdict": "trend" if track["trend"]["flagged"] else "no trend",
                "trend_line": {"first_ordinal": track["trend"]["x"][0],
                               "last_ordinal": track["trend"]["x"][1],
                               "dashed": not track["trend"]["flagged"]},
            }
        plate["panels"][dim] = panel
    return plate


def assert_figure_matches_report(report, drawn):
    """The plate is the report: same points, same dates, same line.

    Runs on every build, so a layout change can never quietly move a value: the
    drawn series has to be the series the JSON was computed from, the drawn
    dates the dates in ``date_range``, and the drawn trend line the line the
    recorded Theil-Sen slope describes.
    """
    problems = []
    for dim in drawn["panel_order"]:
        for orbit, name, _colour in ORBITS:
            block = report["drift"][dim][orbit]
            track = drawn["panels"][dim]["tracks"].get(orbit)
            where = f"{dim}/{name.split()[0]}"
            if not block.get("available"):
                if track is not None:
                    problems.append(f"{where}: drawn but not available")
                continue
            if track is None:
                problems.append(f"{where}: available but not drawn")
                continue
            if track["n"] != block["n"]:
                problems.append(f"{where}: n {track['n']} != {block['n']}")
            if [track["first"], track["last"]] != list(block["date_range"]):
                problems.append(
                    f"{where}: {track['first']}..{track['last']} != "
                    f"{block['date_range']}")
            if abs(track["median"] - block["median"]) > 1e-9:
                problems.append(f"{where}: median {track['median']!r} != "
                                f"{block['median']!r}")
            x = np.array([day for day, _value in track["values"]], float)
            y = np.array([value for _day, value in track["values"]], float)
            per_day = block["theil_sen_per_year"] / 365.25
            intercept = float(np.median(y - per_day * x))
            want = [intercept + per_day * edge for edge in track["trend"]["x"]]
            if not np.allclose(track["trend"]["y"], want, rtol=0, atol=1e-9):
                problems.append(f"{where}: the line is not the recorded slope")
            flagged = "trend" in (block.get("flags") or [])
            if flagged != track["trend"]["flagged"]:
                problems.append(
                    f"{where}: the line style does not match the verdict")
    if problems:
        raise RuntimeError("the figure does not match the report:\n  "
                           + "\n  ".join(problems))


def _f(value, digits=3, signed=True):
    if value is None:
        return "—"
    return f"{value:+.{digits}g}" if signed else f"{value:.{digits}g}"


def _ci(pair, digits=3):
    if not pair:
        return "—"
    return f"[{pair[0]:+.{digits}g}, {pair[1]:+.{digits}g}]"


def _p(value):
    if value is None:
        return "—"
    return "<0.001" if value < 0.001 else f"{value:.3f}"


def _g(value, digits=6):
    """A pinned float: enough digits to re-derive it, few enough to read."""
    if value is None:
        return "—"
    return f"{value:.{digits}g}"


def _pair(pair, digits=6):
    if not pair:
        return "—"
    return f"[{_g(pair[0], digits)}, {_g(pair[1], digits)}]"


def figure_section(figure):
    """The plate as data: size, type, bounds, alpha, one pinned row per cell.

    Written into the ``.md`` because these are properties of the *file*: the
    JSON stays the analysis, the markdown carries the audit of the drawing.
    """
    fonts = figure["final_font_sizes_pt"]
    bounds = figure["bounds"]
    edges = bounds["edge_mm"]
    widths, heights = figure["native_width_mm"], figure["native_height_mm"]
    order = " · ".join(f"`{figure['panels'][dim]['label']}`"
                       for dim in figure["panel_order"])
    lines = [
        "## The figure — one 170 mm plate", "",
        f"* **canvas** `{figure['figsize_inches'][0]:.4f} × "
        f"{figure['figsize_inches'][1]:.4f} in` = **{bounds['canvas_mm'][0]:.2f} × "
        f"{bounds['canvas_mm'][1]:.2f} mm**; measured from the written files: "
        + ", ".join(f"{fmt} {mm:.2f} × {heights[fmt]:.2f} mm"
                    for fmt, mm in sorted(widths.items()))
        + ". The plate is authored at the print width — it is not cropped to its "
        "content, so nothing here is rescaled.",
        f"* **typography** header {fonts['suptitle']} · panel title "
        f"{fonts['title']} · axis label {fonts['label']} · ticks {fonts['tick']} · "
        f"legend {fonts['legend']} · statistics {fonts['annotation']} pt — nothing "
        "below 6 pt, every role inside the `figure_style` band.",
        f"* **panels**, reading order: {order}. The column names of the table are "
        "not on the plate: `gamma2_mask` is drawn as the coherence it measures.",
        "* **statistics rows**: one fixed band above every panel, identical in all "
        f"six — DESC at +{figure['statistics']['row_offset_pt']['DESCENDING']:.1f} "
        f"pt, ASC at +{figure['statistics']['row_offset_pt']['ASCENDING']:.1f} pt, "
        f"{figure['statistics']['pt']} pt, flush with the left spine. Each row is "
        "the Theil–Sen slope per year, its bootstrap 95 % interval, the "
        "Mann–Kendall p and the verdict of the ported flag rule (`trend` when it "
        "fires, `flat` when it does not).",
        "* **marks** 3.6 pt markers with a 0.5 pt white edge above the 1.2 pt trend "
        "line; solid where that verdict is `trend`, dashed where it is `flat`; the "
        "colour is the track's. One figure legend carries the track counts at the "
        "right end of the header line.",
        f"* **bounds** content {bounds['content_mm'][0]:.1f} × "
        f"{bounds['content_mm'][1]:.1f} mm inside the {bounds['canvas_mm'][0]:.1f} "
        f"× {bounds['canvas_mm'][1]:.1f} mm canvas, "
        f"{min(edges.values()):.2f} mm to the nearest edge (floor "
        f"{bounds['minimum_edge_mm']} mm, asserted) — a guard, not a crop. "
        f"LaTeX: {figure['latex_scaling']}.",
        f"* **EPS** {figure['eps']['alpha_artists']} translucent artist(s) "
        "pre-blended over white (`eps=\"blend\"`), "
        f"{figure['eps']['rasterised_artists']} rasterised, text "
        f"{figure['eps']['vector_text']}: the EPS is vector, like the PDF.",
        "* **moved into this caption**: the orbit ids (ASC = relative orbit 160, "
        "DESC = 51), the three shared dates, the trend semantics and the "
        "definitions of the six dimensions.",
        "", "```json",
        json.dumps(figure, indent=1, ensure_ascii=False),
        "```", "",
        "| panel | track | n | first | last | median | IQR | slope / a | 95% CI | "
        "MK p | verdict | y-limits |",
        "|---|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    for dim in figure["panel_order"]:
        panel = figure["panels"][dim]
        for orbit, name, _colour in ORBITS:
            pin = panel["tracks"][orbit]
            if not pin.get("n"):
                continue
            lines.append(
                f"| {panel['label']} | {name} | {pin['n']} | {pin['first']} | "
                f"{pin['last']} | {_g(pin['median'])} | {_pair(pin['iqr'])} | "
                f"{_g(pin['theil_sen_per_year'])} | "
                f"{_pair(pin['slope_ci95_per_year'])} | "
                f"{_p(pin['mann_kendall_p'])} | {pin['verdict']} | "
                f"{_pair(panel['y_limits'], 8)} |")
    lines.append("")
    return lines


def to_markdown(report, figure=None):
    lines = []
    add = lines.append
    add(f"# {report['label']} — the OCV vector of the support chips")
    add("")
    add(f"* table `{report['table']['path']}` (sha256 "
        f"`{report['table']['sha256'][:16]}…`), {report['table']['n_rows']} chips")
    add(f"* {report['dates'][0]} … {report['dates'][1]}; orbits: "
        + ", ".join(f"{k} {v}" for k, v in sorted(report['by_orbit'].items()))
        + f"; **{report['shared_dates']} shared dates**")
    add("* **second instrument**: these are chip values, not the pipeline columns "
        "of the retired pipeline figure (frozen at "
        "`data/espoo/fig_espoo_ocv_orbits.json`: r(log A, log "
        "`coherence_masked_pixels`) = −0.24, r(`gamma2_mask`, "
        "`coherence_gamma2`) = 0.13)")
    add("* **two tracks**: ASC = relative orbit 160, DESC = 51; the other two "
        "tracks place the mast on the far range edge, where the 200 m ring runs "
        "past the swath")
    add("* **no pairing**: the tracks have three shared dates, so the orbit "
        "contrast below is a two-sample comparison")
    add("* **tautology**: `D = A / bbox_area` exactly, and `gamma2_mask` carries "
        "`n = A` in its denominator; `A` differs by a factor 7 between the "
        "orbits, so both dimensions are dominated by the mask size")
    add("")
    if figure:
        lines.extend(figure_section(figure))
        add("")
    add("## 1. Drift across the two years")
    add("")
    add("| dimension | orbit | n | slope / year | 95% CI | MK p | lag-1 | "
        "split-half delta | trend | serial flags | ties |")
    add("|---|---|---|---|---|---|---|---|---|---|---|")
    for dim in DIMS:
        for orbit, name, _colour in ORBITS:
            block = report["drift"][dim][orbit]
            if not block.get("available"):
                add(f"| `{dim}` | {name} | {block.get('n')} | — | — | — | — | — | "
                    "n too small | — | — |")
                continue
            flags = list(block.get("flags") or [])
            trend = "yes" if "trend" in flags else "no"
            serial = ", ".join(f for f in flags if f != "trend") or "none"
            split = (block.get("split_half") or {}).get("delta")
            add(f"| `{dim}` | {name} | {block['n']} | "
                f"{_f(block['theil_sen_per_year'])} | "
                f"{_ci(block['slope_ci95_per_year'])} | "
                f"{_p(block['mann_kendall'].get('p'))} | "
                f"{_f(block.get('lag1_autocorr'), 2)} | {_f(split, 2)} | "
                f"{trend} | {serial} | "
                f"{block.get('top_value_share', 0):.2f} at {block.get('top_value')}"
                + (" **(degenerate)**" if block.get("degenerate") else "") + " |")
    add("")
    add("## 2. Year repeat — the same calendar month in both campaign years")
    add("")
    add("| dimension | orbit | months | median Δ (y2 − y1) | 95% CI | months + / − |")
    add("|---|---|---|---|---|---|")
    for dim in DIMS:
        for orbit, name, _colour in ORBITS:
            block = report["year_repeat"]["by_dim"][dim][orbit]
            if not block.get("available"):
                add(f"| `{dim}` | {name} | {block.get('months')} | — | — | zu wenige "
                    "Monatspaare |")
                continue
            add(f"| `{dim}` | {name} | {block['months']} | {_f(block['median_delta'])} "
                f"| {_ci(block['ci'])} | {block['months_positive']} / "
                f"{block['months_negative']} |")
    add("")
    add("## 3. Orbit contrast, unpaired (ASC R160 vs DESC R51)")
    add("")
    add("| dimension | n ASC/DESC | median ASC | median DESC | Cliff's δ | 95% CI | "
        "Mann-Whitney p |")
    add("|---|---|---|---|---|---|---|")
    for dim in DIMS:
        block = report["orbit_contrast"][dim]
        if not block.get("available"):
            add(f"| `{dim}` | {block.get('n')} | — | — | — | — | — |")
            continue
        add(f"| `{dim}` | {block['n'][0]}/{block['n'][1]} | "
            f"{_f(block['median'][0], 3, False)} | {_f(block['median'][1], 3, False)} | "
            f"{_f(block['delta'], 2)} | {_ci(block['delta_ci95'], 2)} | "
            f"{_p(block['mannwhitney_p'])} |")
    add("")
    add("## What this can and cannot say")
    add("")
    add("* **It is a different instrument from the pipeline figure.** The two "
        "must not be mixed: they correlate at −0.24 (`A`) and 0.13 (`gamma2`).")
    add("* **The series is halved and selected by track.** 75 of 150 overpasses; "
        "the excluded tracks (ASC R58, DESC R153) are named in the site README.")
    add("* **No paired orbit contrast.** Three shared dates; the contrast is a "
        "two-sample comparison and is reported as such.")
    add("* **`D` and `gamma2_mask` are functions of `A`**, which itself differs "
        "by a factor 7 between the orbits.")
    add("")
    return "\n".join(lines)


def write_outputs(report, figure=None, fig_path=None, md_path=None,
                  json_path=None, quiet=False):
    """Write the analysis JSON unchanged, and the markdown with the plate.

    The JSON is the scientific record and does not know about the layout: the
    audit of the drawing goes into the markdown, next to the caption that
    belongs to it.
    """
    fig_path = fig_path or OUT_FIG
    md_path = md_path or OUT_MD
    json_path = json_path or OUT_JSON
    for path in (fig_path, md_path, json_path):
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(json_path, "w") as handle:
        json.dump(report, handle, indent=1)
        handle.write(chr(10))
    with open(md_path, "w") as handle:
        handle.write(to_markdown(report, figure))
        handle.write(chr(10))
    if not quiet:
        print(f"  wrote {os.path.relpath(json_path, REPO)}")
        print(f"  wrote {os.path.relpath(md_path, REPO)}")
    return fig_path, md_path, json_path


def _figure_block(text):
    """The audit block the markdown carries, or ``None``."""
    start = text.find("## The figure")
    if start < 0:
        return None
    fence = text.find("```json", start)
    if fence < 0:
        return None
    end = text.find("```", fence + len("```json"))
    if end < 0:
        return None
    return json.loads(text[fence + len("```json"):end])


def _close(left, right, tol=1e-6):
    return abs(left - right) <= tol * max(1.0, abs(left), abs(right))


def verify_figure(figure, report, table=None, paths=None):
    """Re-measure the three containers and re-derive every pinned number.

    The audit in the markdown is only worth something if it can fail: the three
    containers are measured again from disk, the point sizes are put back
    against the ``figure_style`` bands, and each pinned row is recomputed from
    the CSV -- n, dates, median, slope, verdict, and the y-limits the autoscaler
    must have produced from those points.
    """
    table = table or TABLE
    rows = load(table)
    failures, checked = [], 0

    def check(condition, message):
        nonlocal checked
        checked += 1
        if not condition:
            failures.append(message)

    check(_close(figure["target_width_mm"], style.FINAL_WIDTH_MM, 1e-9),
          "the target width is not the print width")
    check(_close(figure["figsize_inches"][0] * 25.4, style.FINAL_WIDTH_MM, 1e-6),
          "the figsize does not resolve to 170 mm")
    check(list(figure["panel_order"]) == list(PANEL_ORDER), "panel order")
    for path in (paths if paths is not None
                 else (OUT_FIG, OUT_PDF, OUT_EPS)):
        if not os.path.exists(path):
            failures.append(f"{os.path.basename(path)}: missing")
            continue
        info = style.measure(path)
        fmt = info["format"]
        check(abs(info["mm"] - style.FINAL_WIDTH_MM) <= style.WIDTH_TOL_MM,
              f"{fmt}: {info['mm']:.2f} mm, not 170 ± 0.5 mm")
        check(abs(info["height_mm"] - figure["target_height_mm"]) <= 0.5,
              f"{fmt}: {info['height_mm']:.2f} mm high")
        check(abs(info["mm"] - figure["native_width_mm"][fmt]) <= 0.01,
              f"{fmt}: the recorded width is not the measured one")
    for role, (low, high) in style.BANDS.items():
        value = figure["final_font_sizes_pt"].get(role)
        check(value is not None and low - 1e-9 <= value
              and (high is None or value <= high + 1e-9),
              f"{role}: {value} pt is outside the {low}..{high} band")
    for dim in figure["panel_order"]:
        panel = figure["panels"][dim]
        check(panel["label"] == PANEL_LABEL[dim], f"{dim}: label")
        check("gamma2_mask" not in panel["label"], f"{dim}: column name drawn")
        lows, highs = [], []
        for orbit, name, _colour in ORBITS:
            pin, block = panel["tracks"][orbit], report["drift"][dim][orbit]
            where = f"{dim}/{name.split()[0]}"
            series = series_for(rows, orbit, dim)
            check(pin["n"] == len(series),
                  f"{where}: n {pin['n']} != {len(series)}")
            if not series:
                continue
            first = dt.date.fromordinal(int(series[0][0]))
            last = dt.date.fromordinal(int(series[-1][0]))
            check(pin["first"] == first.isoformat()
                  and pin["last"] == last.isoformat(), f"{where}: dates")
            x = np.array([day for day, _value in series], float)
            y = np.array([value for _day, value in series], float)
            check(_close(pin["median"], float(np.median(y))), f"{where}: median")
            slope = st.theil_sen(x, y)
            check(slope is not None
                  and _close(pin["theil_sen_per_year"], slope * 365.25, 1e-9),
                  f"{where}: slope")
            check(_close(pin["theil_sen_per_year"],
                         block["theil_sen_per_year"], 1e-9),
                  f"{where}: the pin is not the report")
            flagged = "trend" in (block.get("flags") or [])
            check(pin["verdict"] == ("trend" if flagged else "no trend"),
                  f"{where}: verdict")
            check(pin["trend_line"]["dashed"] is not flagged,
                  f"{where}: line style")
            intercept = float(np.median(y - slope * x))
            ends = [intercept + slope * pin["trend_line"][key]
                    for key in ("first_ordinal", "last_ordinal")]
            lows += [float(y.min()), min(ends)]
            highs += [float(y.max()), max(ends)]
        if lows:
            low, high = min(lows), max(highs)
            pad = 0.05 * (high - low)
            check(_close(panel["y_limits"][0], low - pad, 1e-6)
                  and _close(panel["y_limits"][1], high + pad, 1e-6),
                  f"{dim}: y-limits {panel['y_limits']}, "
                  f"expected [{low - pad:.6g}, {high + pad:.6g}]")
    return failures, checked


def verify(report=None, table=None):
    if report is None:
        with open(OUT_JSON) as handle:
            report = json.load(handle)
    table = table or TABLE
    import hashlib
    if report["table"]["sha256"] != hashlib.sha256(open(table, "rb").read()).hexdigest():
        print("FAILED — the table on disk is not the one this report came from")
        return 1
    rows = load(table)
    failures, checked = [], 0
    for dim in DIMS:
        for orbit, _name, _colour in ORBITS:
            series = series_for(rows, orbit, dim)
            block = report["drift"][dim][orbit]
            if not block.get("available"):
                continue
            x = np.array([day for day, _v in series], float)
            y = np.array([value for _d, value in series], float)
            slope = st.theil_sen(x, y)
            checked += 3
            if int(y.size) != block["n"]:
                failures.append(f"{dim}/{orbit}: n")
            if abs(float(np.median(y)) - block["median"]) > 1e-9:
                failures.append(f"{dim}/{orbit}: median")
            if slope is None or abs(slope * 365.25 - (block["theil_sen_per_year"] or 0)) > 1e-9:
                failures.append(f"{dim}/{orbit}: slope")
        contrast = report["orbit_contrast"][dim]
        if contrast.get("available"):
            left = [r["values"][dim] for r in rows if r["orbit"] == "ASCENDING"
                    and r["values"][dim] is not None]
            right = [r["values"][dim] for r in rows if r["orbit"] == "DESCENDING"
                     and r["values"][dim] is not None]
            checked += 1
            if abs(cliffs_delta(right, left) - contrast["delta"]) > 1e-9:
                failures.append(f"{dim}: orbit delta")
    print(f"  {checked} field(s) recomputed from {os.path.relpath(table, REPO)}")
    if os.path.exists(OUT_MD):
        with open(OUT_MD) as handle:
            figure = _figure_block(handle.read())
        if figure is None:
            failures.append("figure: the markdown carries no audit block")
        else:
            fig_failures, fig_checked = verify_figure(figure, report, table)
            checked += fig_checked
            failures += fig_failures
            print(f"  {fig_checked} figure pin(s) re-measured and re-derived")
    else:
        print("  (no markdown: the figure audit is not checked)")
    if failures:
        print(f"FAILED — {len(failures)} mismatch(es): {failures[:5]}")
        return 1
    print("  all recomputed values match the JSON and the figure audit")
    return 0


def _synthetic(kind, n_days=80, seed=5):
    rng = np.random.default_rng(seed)
    start = dt.date(2024, 9, 5)
    rows = []
    for index in range(n_days):
        day = start + dt.timedelta(days=12 * index)
        drift = 0.01 * index if kind == "drift" else 0.0
        year = 0.6 if (kind == "year" and (day - start).days >= 365) else 0.0
        for orbit in ("ASCENDING", "DESCENDING"):
            offset = 3.0 if orbit == "ASCENDING" else 0.0
            base = 10.0 + drift + year + offset + rng.normal(0, 0.1)
            rows.append({"date": day, "ordinal": day.toordinal(), "year": day.year,
                         "month": day.month, "orbit": orbit,
                         "values": {dim: base for dim in DIMS}})
    rows.sort(key=lambda row: (row["ordinal"], row["orbit"]))
    return rows


def selftest():
    failures = []

    def check(label, condition, detail=""):
        print(f"  [{'ok  ' if condition else 'FAIL'}] {label}"
              + (f" — {detail}" if detail else ""))
        if not condition:
            failures.append(label)

    print("\n1. Flat, drifting, year shift")
    flat = _synthetic("flat")
    block = st.stationarity_of(series_for(flat, "ASCENDING", "A"), block=5, n_boot=300)
    check("a flat series has no trend flag", not block["slope_excludes_zero"],
          f"{block['theil_sen_per_year']:+.3g}")
    drifting = _synthetic("drift")
    block = st.stationarity_of(series_for(drifting, "ASCENDING", "A"),
                               block=5, n_boot=300)
    check("a planted drift is seen", block["slope_excludes_zero"],
          f"{block['theil_sen_per_year']:+.3g}")
    shifted = _synthetic("year")
    repeat = year_repeat(shifted, boot=300)
    cell = repeat["by_dim"]["A"]["ASCENDING"]
    check("the year shift is caught by the repeat test",
          cell.get("available") and cell["median_delta"] > 0.2,
          f"{cell.get('median_delta')}")

    print("\n2. The unpaired orbit contrast")
    flat = _synthetic("flat")
    contrast = orbit_contrast(flat, boot=200)["A"]
    check("the planted ASC offset comes out positive",
          contrast["delta"] > 0.5, f"delta {contrast['delta']:+.2f}")
    check("and it is marked as unpaired", contrast["paired"] is False)

    print("\n3. The real table")
    report = analyse(boot=100)
    check("75 chips are read", report["table"]["n_rows"] == 75,
          str(report["table"]["n_rows"]))
    check("the track split is reported",
          report["by_orbit"].get("ASCENDING") == 31
          and report["by_orbit"].get("DESCENDING") == 44,
          str(report["by_orbit"]))
    check("only three dates are shared", report["shared_dates"] == 3,
          str(report["shared_dates"]))
    text = to_markdown(report)
    check("the report names the second instrument", "second instrument" in text)
    check("and the tautology", "tautology" in text)

    print("\n4. The plate")
    check("the panels are labelled canonically, no column name is drawn",
          all("gamma2_mask" not in label for label in PANEL_LABEL.values())
          and [dim for dim in PANEL_LABEL] == list(PANEL_ORDER),
          " · ".join(PANEL_LABEL[dim] for dim in PANEL_ORDER))
    check("the coherence panel is math text",
          r"\gamma^2" in PANEL_LABEL["gamma2_mask"])
    check("the canvas is the print width",
          abs(PLATE_IN[0] * 25.4 - style.FINAL_WIDTH_MM) < 1e-9,
          f"{PLATE_IN[0]:.4f} in = {PLATE_IN[0] * 25.4:.2f} mm")
    check("the two statistics rows cannot touch",
          (STAT_ROW["DESCENDING"] - STAT_ROW["ASCENDING"]) > 1.3 * STAT_PT,
          f"{STAT_ROW} at {STAT_PT} pt")
    check("the title clears the DESC row",
          TITLE_PAD > STAT_ROW["DESCENDING"] + 1.3 * STAT_PT,
          f"pad {TITLE_PAD} pt")
    rows = load()
    lines = [statistics_line(report["drift"]["A"][orbit])
             for orbit, _name, _colour in ORBITS]
    check("a statistics row names the track and the verdict",
          lines[0].startswith("ASC ") and lines[0].endswith(("trend", "flat")),
          lines[0])
    scratch = tempfile.mkdtemp(prefix="espoo_plate_")
    try:
        plate, _drawn = make_figure(report, rows,
                                    os.path.join(scratch, "plate.jpg"))
    finally:
        shutil.rmtree(scratch, ignore_errors=True)
    check("all three containers measure the print width",
          sorted(plate["native_width_mm"]) == ["eps", "jpeg", "pdf"]
          and all(abs(mm - style.FINAL_WIDTH_MM) <= style.WIDTH_TOL_MM
                  for mm in plate["native_width_mm"].values()),
          str(plate["native_width_mm"]))
    check("every panel pins both tracks",
          all(sum(1 for pin in plate["panels"][dim]["tracks"].values()
                  if pin.get("n")) == 2 for dim in PANEL_ORDER))
    check("the alpha is blended, not rasterised",
          plate["eps"]["alpha_artists"] > 0
          and plate["eps"]["rasterised_artists"] == 0,
          f"{plate['eps']['alpha_artists']} artist(s)")
    with_audit = to_markdown(report, plate)
    check("the markdown carries the plate audit",
          "## The figure — one 170 mm plate" in with_audit
          and '"native_dpi": 630' in with_audit)
    block = _figure_block(with_audit) or {}
    check("and the block reads back in panel order",
          block.get("panel_order") == list(PANEL_ORDER))
    pin_failures, pin_checked = verify_figure(block, report, paths=())
    check("the pins re-derive from the table", not pin_failures,
          f"{pin_checked} pin(s), {pin_failures[:3]}")

    print()
    if failures:
        print(f"FAILED — {len(failures)} check(s): {failures}")
        return 1
    print("all checks passed")
    return 0


def build_parser():
    parser = argparse.ArgumentParser(
        description="Espoo on the support chips: drift, year repeat, unpaired "
                    "orbit contrast")
    parser.add_argument("--table", default=None)
    parser.add_argument("--fig", "--png", default=None,
                       help="output image (630 dpi JPEG by default)")
    parser.add_argument("--md", default=None)
    parser.add_argument("--json", default=None)
    parser.add_argument("--boot", type=int, default=BOOT)
    parser.add_argument("--selftest", action="store_true")
    parser.add_argument("--verify", action="store_true")
    parser.add_argument("--quiet", action="store_true")
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    if args.selftest:
        return selftest()
    if args.verify:
        return verify(table=args.table)
    print(f"{LABEL} — {LABEL and 'the OCV vector of the support chips'}")
    rows = load(args.table)
    report = analyse(args.table, boot=args.boot)
    figure, _drawn = make_figure(report, rows, args.fig)
    write_outputs(report, figure, args.fig, args.md, args.json, quiet=args.quiet)
    if not args.quiet:
        print(f"  {report['table']['n_rows']} chips, orbits {report['by_orbit']}, "
              f"{report['shared_dates']} shared dates")
        for dim in DIMS:
            line = f"  {dim:11s}"
            for orbit, _name, _colour in ORBITS:
                block = report["drift"][dim][orbit]
                if not block.get("available"):
                    line += "  n too small"
                    continue
                slope = block["theil_sen_per_year"]
                text = (f"{slope:+.2e}" if slope and abs(slope) < 1e-3
                        else f"{slope:+.3g}")
                flag = "trend" if "trend" in (block.get("flags") or []) else "flat "
                line += f"  {orbit[:4]} {text}/a {flag}"
            print(line)
        print(f"  wrote {os.path.relpath(args.fig or OUT_FIG, REPO)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
