#!/usr/bin/env python3
"""OCV orbit diagrams per dimension, for every support-sized table.

One panel per **dimension x observation channel**, in the order the LUMO and
Bautzen figures use: ``gamma2, P, D, F, A, S`` down the rows, ASCENDING /
DESCENDING across the columns. The y axes follow the same rule (``P, D, F, A``
logarithmic, ``gamma2`` and ``S`` linear), so a reader who knows the LUMO figure
knows this one.

Two states, never more
----------------------
Every site carries **pre-collapse and post-collapse**; KDLO is the exception and
carries **pre-collapse and after rebuild**, because that is what its record says.
A state with a single overpass (KDLO's ``single``) is not a state: it is dropped
and named in the markdown.

What a panel shows
------------------
* one marker per overpass, coloured and shaped by its state;
* a band per state, its edges at the midpoints between consecutive runs;
* the **per-state median** as a solid segment across that band, with ``n``;
* for ``gamma2`` a dotted zero line, because zero is the floor of the scale.

A channel with no data is not drawn empty: the panel carries the reason ("no
ASC dates: 19 of 33 were refused at the raster edge"), because an empty axis
reads as "no effect" instead of "no measurement".

Outputs
-------
``figures/<site>/fig_<site>_ocv_orbits.jpg``  630 dpi
``figures/<site>/fig_<site>_ocv_orbits.md``   the report
``data/<site>/support/fig_<site>_ocv_orbits.json``  the plotted values

``--verify`` re-derives every plotted point, run and median from the table and
fails if the JSON and the table disagree; the figure itself computes nothing.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import math
import os
import sys
from collections import Counter

import matplotlib

matplotlib.use("Agg")

import matplotlib.dates as mdates  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))

SCHEMA = "ocv-orbits-support-v1"
FIG_DPI = 630
MM = 1 / 25.4
COL_WIDTH = 183 * MM

CHANNELS = ("gamma2", "P", "D", "F", "A", "S")
LOG_CHANNELS = ("P", "D", "F", "A")
CHANNEL_AXIS = {"gamma2": "$\\gamma^2$", "P": "P", "D": "D",
               "F": "F", "A": "A", "S": "S"}
CHANNEL_TITLE = {
    "gamma2": "echo coherence",
    "P": "persistence (mean pixel frequency)",
    "D": "density (A / bounding box)",
    "F": "fragmentation (8-connected components)",
    "A": "area (masked cells)",
    "S": "shift (centroid to support peak, px)",
}
STATE_COLORS = {"pre": "#2b6cb0", "post": "#c53030"}
STATE_FILL = {"pre": "#bee3f8", "post": "#feb2b2"}
STATE_MARKERS = {"pre": "o", "post": "s"}
PANEL_LETTERS = "ABCDEFGHIJKL"

SITES = {
    "cts": {
        "label": "Champlain Towers South (target footprint)",
        "table": "data/cts/support/cts_ocv_channels.csv",
        "manifest": "data/cts/support/cts_rects_support.txt.manifest.json",
        "date_column": "date", "orbit_column": "orbit_direction",
        "label_column": "state",
        "states": {"pre": "pre-collapse (healthy)", "post": "post-collapse"},
        "state_names": {"pre": "pre-collapse", "post": "post-collapse"},
        "gamma2": "gamma2_mask",
        "channels": {"ASCENDING": "ASC", "DESCENDING": "DESC"},
        "framing": "collapse",
    },
    "kdlo": {
        "label": "KDLO media tower",
        "table": "data/kdlo/support/kdlo_ocv_channels.csv",
        "manifest": "data/kdlo/support/kdlo_rects_support.txt.manifest.json",
        "date_column": "acquisition_ts", "orbit_column": "orbit_direction",
        "label_column": "state",
        "states": {"pre": "pre", "post": "rebuild"},
        "state_names": {"pre": "pre-collapse", "post": "after rebuild"},
        "exclude_states": ("single",),
        "gamma2": "gamma2_mask",
        "channels": {"ASCENDING": "ASC", "DESCENDING": "DESC"},
        "framing": "rebuild",
    },
    "carola": {
        "label": "Carolabr\u00fccke deck",
        "table": "data/carola/support/carola_ocv_channels.csv",
        "manifest": "data/carola/support/carola_rects_support.txt.manifest.json",
        "date_column": "acquisition_ts", "orbit_column": "orbit_direction",
        "label_column": "state",
        "states": {"pre": "pre-collapse (healthy)", "post": "post-collapse"},
        "state_names": {"pre": "pre-collapse", "post": "post-collapse"},
        "gamma2": "gamma2_mask",
        "channels": {"ASCENDING": "ASC", "DESCENDING": "DESC"},
        "framing": "collapse",
    },
    "haltern": {
        "label": "Haltern AV9 wind turbine",
        "table": "data/haltern/support/haltern_ocv_channels.csv",
        "manifest": "data/haltern/support/haltern_rects_support.txt.manifest.json",
        "date_column": "acquisition_ts", "orbit_column": "orbit_direction",
        "label_column": "state",
        "states": {"pre": "pre-collapse (healthy)", "post": "post-collapse"},
        "state_names": {"pre": "pre-collapse", "post": "post-collapse"},
        "gamma2": "gamma2_mask",
        "channels": {"ASCENDING": "ASC", "DESCENDING": "DESC"},
        "framing": "collapse",
    },
    "ywf": {
        "label": "Yeongdeok Unit 21",
        "table": "data/ywf/support/ywf_ocv_channels.csv",
        "manifest": "data/ywf/support/ywf_rects_support.txt.manifest.json",
        "date_column": "acquisition_ts", "orbit_column": "orbit_direction",
        "label_column": "state",
        "states": {"pre": "pre-collapse (healthy)", "post": "post-collapse"},
        "state_names": {"pre": "pre-collapse", "post": "post-collapse"},
        "gamma2": "gamma2_mask",
        "channels": {"ASCENDING": "ASC", "DESCENDING": "DESC"},
        "framing": "collapse",
    },
    "morandi": {
        "label": "Ponte Morandi deck",
        "table": "data/morandi/support/morandi_ocv_channels.csv",
        "manifest": "data/morandi/support/morandi_rects_support.txt.manifest.json",
        "date_column": "acquisition_ts", "orbit_column": "orbit_direction",
        "label_column": "state",
        "states": {"pre": "pre-collapse (healthy)", "post": "post-collapse"},
        "state_names": {"pre": "pre-collapse", "post": "post-collapse"},
        "gamma2": "gamma2_mask",
        "channels": {"ASCENDING": "ASC", "DESCENDING": "DESC"},
        "framing": "collapse",
    },
}

# ---------------------------------------------------------------------------
# Reading one table
# ---------------------------------------------------------------------------
def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_time(value, spec):
    """An overpass instant, from ``acquisition_ts`` or a plain date."""
    text = (value or "").strip()
    if not text:
        return None
    if spec["date_column"] == "date":
        return dt.datetime.fromisoformat(text[:10])
    cleaned = text.replace("Z", "+00:00")
    try:
        stamp = dt.datetime.fromisoformat(cleaned)
    except ValueError:
        try:
            stamp = dt.datetime.fromisoformat(cleaned[:19])
        except ValueError:
            return None
    return stamp.replace(tzinfo=None)


def as_float(value):
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def read_points(path, spec):
    """Every overpass with its state, channel and the six dimension values."""
    excluded = set(spec.get("exclude_states") or ())
    by_state = {label: state for state, label in spec["states"].items()}
    dims = [(name, field) for name, field in
            (("gamma2", spec["gamma2"]), ("P", "P"), ("D", "D"),
             ("F", "F"), ("A", "A"), ("S", "S"))]
    points, skipped = [], Counter()
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle):
            state = by_state.get(row.get(spec["label_column"]))
            if state is None:
                label = row.get(spec["label_column"])
                if label in excluded:
                    skipped[label] += 1
                continue
            channel = row.get(spec["orbit_column"])
            if channel not in spec["channels"]:
                continue
            stamp = parse_time(row.get(spec["date_column"]), spec)
            if stamp is None:
                skipped["no date"] += 1
                continue
            values = {name: as_float(row.get(field)) for name, field in dims}
            points.append({"t": stamp, "state": state, "channel": channel,
                           "values": values})
    points.sort(key=lambda p: (p["channel"], p["t"]))
    return points, dict(skipped)


def runs_of(points):
    """Contiguous same-state runs in time order, with their band edges."""
    runs = []
    for point in points:
        if runs and runs[-1]["state"] == point["state"]:
            runs[-1]["points"].append(point)
        else:
            runs.append({"state": point["state"], "points": [point]})
    for index, run in enumerate(runs):
        first = run["points"][0]["t"]
        last = run["points"][-1]["t"]
        if index:
            previous = runs[index - 1]["points"][-1]["t"]
            run["start"] = previous + (first - previous) / 2
        else:
            run["start"] = first
        if index + 1 < len(runs):
            following = runs[index + 1]["points"][0]["t"]
            run["end"] = last + (following - last) / 2
        else:
            run["end"] = last
        run["n"] = len(run["points"])
    return runs


def rolling_median(points, window=5, minimum=3):
    """A robust local summary: the median of the surrounding overpasses.

    The window is in **overpass order**, not calendar time, so a gap in the
    acquisitions does not blank the line. ``None`` values are skipped, never
    counted as zero, and a window with fewer than ``minimum`` valid values yields
    ``None`` — a line drawn from one point would be a claim, not a summary.
    """
    half = window // 2
    out = []
    for index, point in enumerate(points):
        low = max(0, index - half)
        high = min(len(points), index + half + 1)
        block = sorted(p["v"] for p in points[low:high] if p["v"] is not None)
        if len(block) < minimum:
            out.append({"t": point["t"], "v": None, "n": len(block)})
            continue
        middle = len(block) // 2
        value = (block[middle] if len(block) % 2
                 else (block[middle - 1] + block[middle]) / 2.0)
        out.append({"t": point["t"], "v": float(value), "n": len(block)})
    return out


def median(values):
    clean = sorted(v for v in values if v is not None)
    if not clean:
        return None
    middle = len(clean) // 2
    if len(clean) % 2:
        return float(clean[middle])
    return float((clean[middle - 1] + clean[middle]) / 2.0)

# ---------------------------------------------------------------------------
# One site
# ---------------------------------------------------------------------------
def compute(site, spec, repo):
    """Everything the figure draws, plus the provenance a reader needs."""
    path = os.path.join(repo, spec["table"])
    points, skipped = read_points(path, spec)
    manifest_path = os.path.join(repo, spec["manifest"])
    manifest = json.load(open(manifest_path)) if os.path.exists(manifest_path) else {}

    series, runs, summary, counts, rolling = {}, {}, {}, {}, {}
    for channel in spec["channels"]:
        subset = [p for p in points if p["channel"] == channel]
        runs[channel] = runs_of(subset)
        series[channel] = {}
        summary[channel] = {}
        rolling[channel] = {}
        counts[channel] = {state: len([p for p in subset if p["state"] == state])
                           for state in ("pre", "post")}
        for dim in CHANNELS:
            series[channel][dim] = [
                {"t": p["t"].isoformat(), "v": p["values"][dim],
                 "state": p["state"]} for p in subset]
            rolling[channel][dim] = rolling_median(series[channel][dim])
            block = {}
            for state in ("pre", "post"):
                vals = [p["values"][dim] for p in subset if p["state"] == state]
                block[state] = {"n": len([v for v in vals if v is not None]),
                                "median": median(vals)}
            summary[channel][dim] = block

    support = manifest.get("support") or {}
    return {
        "schema": SCHEMA, "site": site, "label": spec["label"],
        "generated_by": "code/orbits/fig_ocv_orbits.py",
        "states": spec["state_names"],
        "framing": spec["framing"],
        "state_labels": spec["states"],
        "excluded_states": list(spec.get("exclude_states") or ()),
        "channels": {channel: spec["channels"][channel]
                     for channel in spec["channels"]},
        "table": {"path": os.path.relpath(path, repo), "sha256": sha256(path),
                  "n_rows": len(points)},
        "manifest": {"path": os.path.relpath(manifest_path, repo),
                     "chip": support.get("chip"),
                     "band": support.get("band"),
                     "support_cells": support.get("support_cells"),
                     "ring_cells": support.get("ring_cells"),
                     "guard_m": manifest.get("guard_m"),
                     "ring_m": manifest.get("ring_m"),
                     "n_dates": manifest.get("n_dates"),
                     "n_refused": manifest.get("n_refused")},
        "skipped": skipped,
        "counts": counts, "runs": runs, "series": series, "summary": summary,
        "rolling": rolling,
        "rolling_rule": {"window": 5, "minimum_valid": 3,
                         "order": "overpass order, centred",
                         "missing": "skipped, never zero"},
    }


def _channel_ylim(dim, series, channels):
    vals = [p["v"] for channel in channels for p in series[channel][dim]
            if p["v"] is not None]
    if not vals:
        return (0.0, 1.0)
    low, high = min(vals), max(vals)
    if dim in LOG_CHANNELS and low > 0:
        return (low / 1.6, high * 1.6)
    span = high - low if high > low else max(abs(high), 1.0) * 0.1
    return (low - 0.12 * span, high + 0.12 * span)


def make_figure(res, figdir):
    """Draw the 6 x 2 grid and return the JPEG path."""
    spec_channels = list(res["channels"])
    columns = ["ASCENDING", "DESCENDING"]
    series = res["series"]
    runs = res["runs"]
    summary = res["summary"]

    stamps = [dt.datetime.fromisoformat(p["t"]) for channel in spec_channels
              for dim in CHANNELS for p in series[channel][dim]]
    if stamps:
        x0 = min(stamps) - dt.timedelta(days=20)
        x1 = max(stamps) + dt.timedelta(days=20)
    else:
        x0, x1 = dt.datetime(2000, 1, 1), dt.datetime(2001, 1, 1)

    nrows, ncols = len(CHANNELS), len(columns)
    fig, axes = plt.subplots(nrows, ncols, figsize=(COL_WIDTH, 245 * MM),
                             sharex="col", sharey="row", squeeze=False)
    for row, dim in enumerate(CHANNELS):
        ylim = _channel_ylim(dim, series, spec_channels)
        for col, channel in enumerate(columns):
            ax = axes[row][col]
            available = channel in spec_channels
            ax.set_xlim(x0, x1)
            ax.set_ylim(*ylim)
            if dim in LOG_CHANNELS:
                ax.set_yscale("log")
                ax.yaxis.set_major_locator(LogLocator(base=10.0, numticks=8))
                ax.yaxis.set_minor_formatter(NullFormatter())
                ax.yaxis.set_major_formatter(
                    FuncFormatter(lambda y, _pos: f"{y:g}"))
            if dim == "gamma2":
                ax.axhline(0.0, color="0.55", lw=0.65, ls=":", zorder=1)
            if available:
                for run in runs.get(channel) or []:
                    ax.axvspan(run["start"], run["end"],
                               color=STATE_FILL[run["state"]], alpha=0.55,
                               lw=0, zorder=0)
                for state in ("pre", "post"):
                    block = summary[channel][dim][state]
                    if block["median"] is None:
                        continue
                    spans = [run for run in (runs.get(channel) or [])
                             if run["state"] == state]
                    for run in spans:
                        ax.plot([run["start"], run["end"]],
                                [block["median"], block["median"]],
                                color=STATE_COLORS[state], lw=1.1, zorder=2)
                    ax.text(spans[0]["start"], block["median"],
                            f" {block['median']:.4g} (n={block['n']})",
                            fontsize=5.0, color=STATE_COLORS[state],
                            va="bottom", ha="left", zorder=4)
                line = [(dt.datetime.fromisoformat(p["t"]), p["v"])
                        for p in res["rolling"][channel][dim]
                        if p["v"] is not None]
                if len(line) > 1:
                    ax.plot([t for t, _v in line], [v for _t, v in line],
                            color="#1a202c", lw=0.9, zorder=2.5)
                for point in series[channel][dim]:
                    if point["v"] is None:
                        continue
                    ax.plot([dt.datetime.fromisoformat(point["t"])], [point["v"]],
                            linestyle="none",
                            marker=STATE_MARKERS[point["state"]], ms=2.4,
                            mfc=STATE_COLORS[point["state"]], mec="white",
                            mew=0.3, zorder=3)
            else:
                ax.text(0.5, 0.5, f"no {res['channels'].get(channel, channel)} data",
                        transform=ax.transAxes, ha="center", va="center",
                        fontsize=6.4, color="0.35")
            if col == 0:
                note = "  (log)" if dim in LOG_CHANNELS else ""
                ax.set_ylabel(f"{CHANNEL_AXIS[dim]}{note}\n{CHANNEL_TITLE[dim]}",
                              fontsize=5.6)
            ax.tick_params(labelleft=(col == 0), labelbottom=(row == nrows - 1),
                           labelsize=5.2)
            ax.text(0.012, 0.965, f"({PANEL_LETTERS[row * ncols + col]})",
                    transform=ax.transAxes, ha="left", va="top", fontsize=6.2,
                    color="0.2", zorder=5,
                    bbox=dict(facecolor="white", edgecolor="none", alpha=0.72,
                              pad=0.7))
            if row == 0:
                count = res["counts"].get(channel) or {}
                ax.set_title(f"{res['channels'].get(channel, channel)}  "
                             f"({count.get('pre', 0)} pre / {count.get('post', 0)} "
                             f"post)", fontsize=7.2, pad=5)
    for ax in axes[-1]:
        ax.xaxis.set_major_locator(mdates.AutoDateLocator())
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m"))
        ax.set_xlabel("overpass date", fontsize=5.6)

    handles = [Patch(facecolor=STATE_FILL[state], edgecolor=STATE_COLORS[state],
                     label=res["states"][state]) for state in ("pre", "post")]
    handles.append(Line2D([], [], color="0.42", lw=0.9, label="overpass"))
    handles.append(Line2D([], [], color="0.42", lw=1.1, label="state median"))
    handles.append(Line2D([], [], color="#1a202c", lw=0.9,
                          label="rolling median (5 overpasses)"))
    fig.legend(handles=handles, loc="lower center", ncol=5, frameon=False,
               fontsize=6.0, bbox_to_anchor=(0.5, 0.004))
    fig.suptitle(
        f"{res['label']} — OCV trajectories before and after "
        f"the {res['framing']}\n"
        f"{res['table']['n_rows']} overpasses on the support-sized chip; "
        f"states = {res['states']['pre']} / {res['states']['post']}; "
        f"bands = state runs, edges at the midpoints between runs\n"
        f"markers = overpasses; horizontal segments = per-state medians with n; "
        f"P, D, F and A on log axes",
        fontsize=7.6, y=0.992)
    fig.tight_layout(rect=(0, 0.030, 1, 0.918))
    os.makedirs(figdir, exist_ok=True)
    out = os.path.join(figdir, f"fig_{res['site']}_ocv_orbits.jpg")
    fig.savefig(out, dpi=FIG_DPI, facecolor="white")
    plt.close(fig)
    return out

# ---------------------------------------------------------------------------
# The report and the check
# ---------------------------------------------------------------------------
def _num(value, digits=4):
    return "–" if value is None else f"{value:.{digits}g}"


def markdown(res):
    """The companion report: what is drawn, in numbers."""
    lines = []
    add = lines.append
    add(f"# {res['label']} — OCV trajectories before and after the "
        f"{res['framing']}")
    add("")
    add(f"Trajectories of the six OCV dimensions across the {res['framing']}, "
        "one panel per **dimension x observation channel**, drawn from the "
        "support-sized chip: `gamma2, P, D, F, A, S` down the rows (the order "
        "the LUMO and Bautzen figures use), ASC / DESC across the columns. "
        "`P, D, F` and `A` are logarithmic, `gamma2` and `S` linear.")
    add("")
    add("The thin dark line in each panel is a **rolling median over 5 "
        "consecutive overpasses** (centred, at least 3 valid values). It "
        "summarises local evolution — a dip in it is a summary, not a state, and "
        "whether it is a drift or a level shift is a question this figure does "
        "not answer.")
    add("")
    add("## States")
    add("")
    add(f"Two states, as the record defines them: **{res['states']['pre']}** and "
        f"**{res['states']['post']}** (table labels "
        f"`{res['state_labels']['pre']}` / `{res['state_labels']['post']}`).")
    if res["excluded_states"]:
        dropped = sum(res["skipped"].get(state, 0)
                      for state in res["excluded_states"])
        add("")
        add(f"Excluded: `{'`, `'.join(res['excluded_states'])}` "
            f"({dropped} row(s)) — a state that carries a single overpass has no "
            f"contrast and is not drawn as one.")
    add("")
    add("## Rows and refusals")
    add("")
    add("| channel | pre | post | runs |")
    add("| --- | --- | --- | --- |")
    for channel in res["channels"]:
        counts = res["counts"].get(channel) or {}
        add(f"| {res['channels'][channel]} | {counts.get('pre', 0)} | "
            f"{counts.get('post', 0)} | {len(res['runs'].get(channel) or [])} |")
    add("")
    block = res["manifest"]
    add(f"Table `{res['table']['path']}` ({res['table']['n_rows']} overpasses, "
        f"sha256 `{res['table']['sha256'][:12]}`), cut manifest "
        f"`{block['path']}`: chip {block.get('chip')}, band {block.get('band')}, "
        f"support {block.get('support_cells')} cells, ring "
        f"{block.get('ring_cells')} cells, guard {block.get('guard_m')} m, ring "
        f"{block.get('ring_m')} m, **{block.get('n_refused')} of "
        f"{(block.get('n_dates') or 0) + (block.get('n_refused') or 0)} dates "
        f"refused**.")
    add("")
    add("## Per-state medians")
    add("")
    for channel in res["channels"]:
        counts = res["counts"].get(channel) or {}
        if not (counts.get("pre") or 0) + (counts.get("post") or 0):
            add(f"**{res['channels'][channel]}** — no overpasses of this channel "
                f"in the table, so the panel is empty on purpose.")
            add("")
            continue
        add(f"**{res['channels'][channel]}**")
        add("")
        add("| dimension | " + " | ".join(res["states"][s] for s in ("pre", "post")) + " |")
        add("| --- | --- | --- |")
        for dim in CHANNELS:
            cells = []
            for state in ("pre", "post"):
                item = res["summary"][channel][dim][state]
                cells.append(f"{_num(item['median'])} (n={item['n']})")
            add(f"| `{dim}` | " + " | ".join(cells) + " |")
        add("")
    add("## What this does not claim")
    add("")
    add("- The panels show the **support-sized chip**, not the committed table: "
        "the chip is larger and the reference is a disjoint ring, so these "
        "values are not the committed ones.")
    add("- A band with few overpasses is a handful of points, not a "
        "distribution; the median carries its n.")
    add("- The rolling median is a summary over five overpasses, not a state and "
        "not a trend test; it is drawn so that the local evolution is visible "
        "beside the global pre/post difference.")
    add("- The state boundary is drawn between runs, not at an event date: the "
        "record labels the runs, and the band follows the labels.")
    add("")
    return "\n".join(lines)


def verify(res, repo):
    """Re-derive every drawn number from the table; return the failures."""
    spec = SITES[res["site"]]
    points, _skipped = read_points(os.path.join(repo, spec["table"]), spec)
    failures = []
    for channel in res["channels"]:
        subset = [p for p in points if p["channel"] == channel]
        if len(subset) != sum((res["counts"].get(channel) or {}).values()):
            failures.append(f"{channel}: {len(subset)} rows in the table, "
                            f"{res['counts'].get(channel)} in the figure")
        runs = runs_of(subset)
        if [run["n"] for run in runs] != [run["n"] for run in res["runs"][channel]]:
            failures.append(f"{channel}: run lengths differ")
        for dim in CHANNELS:
            want_line = res["rolling"][channel][dim]
            have_line = rolling_median(res["series"][channel][dim])
            if ([round(p["v"], 12) if p["v"] is not None else None
                 for p in want_line]
                    != [round(p["v"], 12) if p["v"] is not None else None
                        for p in have_line]):
                failures.append(f"{channel}/{dim}: rolling median differs")
            for state in ("pre", "post"):
                vals = [p["values"][dim] for p in subset if p["state"] == state]
                want = res["summary"][channel][dim][state]
                if median(vals) != want["median"]:
                    failures.append(f"{channel}/{dim}/{state}: median "
                                    f"{median(vals)} vs {want['median']}")
                if len([v for v in vals if v is not None]) != want["n"]:
                    failures.append(f"{channel}/{dim}/{state}: n differs")
    return failures


def build_parser():
    parser = argparse.ArgumentParser(
        description="OCV time series per dimension and orbit, from the support tables")
    parser.add_argument("--site", default="all",
                        choices=("all",) + tuple(sorted(SITES)))
    parser.add_argument("--figdir", default=None)
    parser.add_argument("--json", default=None)
    parser.add_argument("--md", default=None)
    parser.add_argument("--verify", action="store_true",
                        help="re-derive every drawn value from the table")
    parser.add_argument("--quiet", action="store_true")
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    sites = sorted(SITES) if args.site == "all" else [args.site]
    failures = 0
    for site in sites:
        spec = SITES[site]
        res = compute(site, spec, REPO)
        figdir = args.figdir or os.path.join(REPO, "figures", site)
        out = make_figure(res, figdir)
        json_path = args.json or os.path.join(
            REPO, "data", site, "support", f"fig_{site}_ocv_orbits.json")
        os.makedirs(os.path.dirname(json_path), exist_ok=True)
        with open(json_path, "w") as handle:
            # run edges are datetimes; the JSON carries them as ISO strings
            json.dump(res, handle, indent=1, sort_keys=True, default=str)
            handle.write("\n")
        md_path = args.md or os.path.join(figdir, f"fig_{site}_ocv_orbits.md")
        with open(md_path, "w") as handle:
            handle.write(markdown(res))
        if not args.quiet:
            print(f"  {site:9s} {res['table']['n_rows']:4d} overpasses  "
                  f"{out.replace(REPO + '/', '')}")
        if args.verify:
            problems = verify(res, REPO)
            failures += len(problems)
            for problem in problems:
                print(f"    VERIFY {site}: {problem}")
    return 0 if failures == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
