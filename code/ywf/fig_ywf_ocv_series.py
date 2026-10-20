#!/usr/bin/env python3
"""YWF OCV-Series — the echo-mask vector over the *whole* record, on one dated axis.

Question: what does the observability core vector ``OCV = [gamma2, P, D]`` — and
the whole six-dimensional echo-mask vector ``x = [gamma2, P, D, A, F, S]`` it is
read off — do over the **entire** committed record of the Yeongdeok Wind Farm
(Changpo Wind Power Complex, Samgye-ri), from the first committed acquisition
(``2025-08-06``) to the last one (``2026-04-27``), with the collapse of the 100 m
monopole of Unit 21 (``2026-02-02``) drawn as a cut inside that axis? The weeks
companion (``fig_ywf_ocv_weeks.py``) opens the 26 pre-event days; the pooled
figure (``fig_ywf_ocv_paper.py``) collapses them into one pre median. This is the
time series of the whole record: 33 chips on 32 acquisition days, one dated axis.

What this axis carries, and what it does not
-------------------------------------------
  * **the record, not the states.** 33 unique 7x7 chips on 32 acquisition days
    (57 committed windows, 24 of them byte-identical duplicates of another
    request for the same chip). The event cuts them 27/6 — and it falls
    **inside** the 12-day step ``2026-01-27 -> 2026-02-08``: no chip sits on the
    collapse date, the nearest pre chip is 6 days before it and the nearest post
    chip 6 days after it.
  * **the post side is a coverage fact, not a level.** Of the 6 post-event chips
    exactly **1** carries an echo (``2026-03-16``, ascending, ``A = 5`` px of 49,
    ``intermediate``); the other 5 carry no echo at all and are drawn on the floor
    of every panel. Nothing on the post side is a distribution — one chip is not a
    level — and the figure never draws a line through that region except across
    the single 6-day step its neighbours share.
  * **a step longer than the cadence is a hole, never an interpolation.** The
    record's step is 6 days (the ascending and the descending pass alternate);
    **9** steps are longer — 7 of 12 days and 2 of 24 days — and every panel draws
    them as bands. A line joins two neighbouring chips **only** when they are one
    cadence step apart, so the drawn line breaks at each hole and at the collapse
    instead of bridging them.
  * **the y axis is the pre-window median, so the ordinates are the paper's.** Each
    dimension is drawn as ``value / pre-window median`` of that dimension, the
    median taken over the **16** pre-event echo-bearing chips — the same 16 the
    paper figure's pre summary is computed on, which is why every drawn ordinate of
    the pre side is a pinned number and not a new one. The post dot is therefore
    *one pre-event-scale reading*, not a post median.
  * **``P`` is a window-wide reference, not a per-chip quantity.** It is the mean
    pixel frequency inside the own mask over the **17** echo-bearing chips of the
    whole committed record (the cache's own reference count), so the chips of this
    figure are not independent draws of it — the report says so.
  * **no trend claim.** The figure draws the record and reports the per-state
    summaries the paper already commits (``summary_by_state``); it runs no test of
    its own and adds no number to the A-E pin contract of the paper script.

Reading the panels
------------------
  x = the acquisition day itself (a dated axis, one dot per chip, no interpolated
      point; the two orbits of a day are separated by a sub-day offset)
  y = value / pre-window median of that dimension, 1.0 drawn
  a  gamma2 (log)      masked coherence        d  D (log)  mask density A / bbox
  b  A (log)           masked pixels of 49     e  F (log)  8-connected components
  c  P (linear)        mask persistence        f  S (log)  ||centroid - peak||, px
  the panels are the family's series convention (``fig_kdlo_ocv_series.py``): the
  bounded pixel frequency ``P`` is the only linear panel, the five multiplicative
  mask magnitudes are logarithmic — all six ratios stay inside 0.12..2.91 of 1.0,
  so no axis hides a level. Dots are coloured by pass geometry (ascending /
  descending) because that is the record's only cross-cut beside the state; the one
  post-event echo is drawn hollow, the 16 echo-free chips on the panel floor.

Cross-checks: every run re-derives the record from the committed channel table and
pins it, field by field, against five committed layers — the generated metadata
(guard counts, the dedup counts, the payload accounting, the coverage per state
*and* per pass geometry), the mask-cache manifest, the mask cache itself (the 6D
vector of all 17 echo rows, the mask rule re-verified row by row), the paper result
JSON (the record counts, the coverage table and its Fisher test, the 32-day
timeline, the per-state channel summaries) and the analysis reference (the same
timelines, the coverage blocks, the sections block). A deviation refuses the write.

Usage:
  python3 code/ywf/fig_ywf_ocv_series.py                 # figure + JSON + report
  python3 code/ywf/fig_ywf_ocv_series.py --no-figures     # cross-checks only
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import math
import os
import sys
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as mdates  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.ticker import NullLocator  # noqa: E402
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import ywf_ocv_core as core     # noqa: E402
import ywf_ocv_masks as masks   # noqa: E402
import ywf_ocv_stats as st      # noqa: E402

FIGDIR = core.FIGDIR
OUT_JSON = os.path.join(core.DATA, "fig_ywf_ocv_series.json")
OUT_MD = os.path.join(FIGDIR, "fig_ywf_ocv_series.md")
FIG_NAME = "fig_ywf_ocv_series"

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

EVENT = dt.date.fromisoformat(core.COLLAPSE_DATE)
OCV = list(core.OCV)
FEATURES = list(core.FEATURES)
# Panel order (a..f) and the y scale of each: the family's series convention
# (``PANEL_SCALE`` of ``fig_kdlo_ocv_series.py``) — the bounded pixel frequency
# ``P`` linear, the five multiplicative mask magnitudes logarithmic.
PANELS = ["gamma2", "A", "P", "D", "F", "S"]
PANEL_GRID = (3, 2)       # the panel grid: three rows, two columns (a..f)
PANEL_BAND = (0.12, 2.91)  # the ratio band every drawn ordinate stays inside
PANEL_SCALE = {"gamma2": "log", "A": "log", "P": "linear",
               "D": "log", "F": "log", "S": "log"}
SCALE_RULE = ("P (a bounded pixel frequency) is the only linear panel; gamma2, A, "
              "D, F and S are multiplicative mask magnitudes drawn on log axes, as "
              "in the KDLO series companion")

# The few y positions every log panel labels, as plain numbers (the ratio band is
# narrower than a decade, so matplotlib would otherwise label every minor tick).
LOG_TICKS = (0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0)

CADENCE_DAYS = 6          # the repeat step of the two alternating orbits
FLOAT_TOL = 1e-12
WINDOW_PX = 49            # 7x7: the mask window, the denominator of D and A

# The co-variates the paper's own control block carries (``core.CONTROLS``), and
# the two committed prose notes this figure re-states: both are pinned word for
# word, so a reworded commit fires a failure instead of silently diverging.
CONTROLS = list(core.CONTROLS)
POWER_NOTE = ("1 of 6 post-event chips carries an echo, so there is no post-state "
              "distribution: the collapse contrast of this site is an echo-coverage "
              "statement with this caveat attached, not a damage verdict.")
ORBIT_NOTE = ("both orbits straddle the collapse date and the only post-event echo "
              "chip is an ascending pass: this block describes where the echoes sit "
              "over the two pass geometries, it is not an orbit contrast of the states")
PIPELINE_NOTE = ("a different generation: the column was computed on the asset-point "
                 "chip of the whole-scene window, while this package recomputes the "
                 "echo mask of the committed 7x7 window")

# The collapse cut (the two states) and the pass geometry (the record's second
# cross-cut — the colour of the dots).
PRE, POST = core.PRE, core.POST
ORBITS = list(core.ORBITS)
ORBIT_COLOUR = {"ASCENDING": "#4C72B0", "DESCENDING": "#55A868"}
STATE_COLOUR = {PRE: "#4C72B0", POST: "#DD8452"}
POST_COLOUR = "#C44E52"
NO_ECHO_COLOR = "#B0B0B0"
LINE_COLOUR = "#9A9A9A"
CUT_COLOR = "#888888"
GAP_COLOR = "#E8E8E8"
MONTH_COLOR = "#DDDDDD"
PRE_TINT_ALPHA = 0.050
POST_TINT_ALPHA = 0.085
MONTHS = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct",
          "Nov", "Dec"]


def _fail(msg):
    """Abort without writing: the artifact must never outrun its evidence."""
    raise SystemExit(f"REFUSING TO WRITE — {msg}")


# ---------------------------------------------------------------------------
# Small helpers (the shape of the family's companion scripts)
# ---------------------------------------------------------------------------
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


def day_date(day):
    """``YYYY-MM-DD`` -> ``datetime.date``."""
    return dt.date.fromisoformat(str(day))


def day_span(first, last):
    """Inclusive day count of ``[first, last]``."""
    return (day_date(last) - day_date(first)).days + 1


def shift_day(day, days):
    """ISO day ``days`` after ``day``."""
    return (day_date(day) + dt.timedelta(days=days)).isoformat()


def orbit_dx(orbit, unit=0.18):
    """A deterministic sub-day offset for the two orbits of one acquisition day.

    Exactly one day of this record carries two chips (``2025-08-12``: one
    ascending, one descending) and both are echo-free, so the offset only keeps
    their two floor ticks from landing on the same pixel of a 265-day axis.
    """
    lab = str(orbit or "").upper()
    if lab.startswith("ASC"):
        return unit
    if lab.startswith("DESC"):
        return -unit
    return 0.0


def x_of(day, orbit=None):
    """The x coordinate of a chip: its acquisition day, plus the orbit offset."""
    off = orbit_dx(orbit) if orbit else 0.0
    return day_date(day) + dt.timedelta(days=off)


def row_view(r):
    """The columns of one chip this figure uses (the 6D vector + the co-variates)."""
    d = {
        "id": r.get("id"), "day": core.day_of(r), "orbit": core.orbit_of(r),
        "pass": r.get("pass_label"), "status": r.get("status"),
        "windows": int(core._num(r.get("payload_group_size")) or 1),
        "echo": bool(r["echo_mask_present"]),
        "echo_mode": r.get("echo_mode"),
        "gamma2": core._num(r.get("gamma2_mask")), "P": core._num(r.get("P")),
        "D": core._num(r.get("D")), "A": core._num(r.get("A")),
        "F": core._num(r.get("F")), "S": core._num(r.get("S")),
        "coherence_masked_pixels": core._num(r.get("coherence_masked_pixels")),
        "intensity": core._num(r.get("intensity")),
        "coherence": core._num(r.get("coherence")),
        "coherence_gamma2": core._num(r.get("coherence_gamma2")),
        "brightness_ratio": core._num(r.get("brightness_ratio")),
        "cumulative_mm": core._num(r.get("cumulative_mm")),
        "structural_frequency_hz": core._num(r.get("structural_frequency_hz")),
        "precipitation_mm": core._num(r.get("precipitation_mm")),
        "relative_humidity_pct": core._num(r.get("relative_humidity_pct")),
        "section": r.get("segment"),
        "section_index": core._num(r.get("section_index")),
        "temperature_c": core._num(r.get("temperature_c")),
        "wind_speed_ms": core._num(r.get("wind_speed_ms")),
        "wind_gust_ms": core._num(r.get("wind_gust_ms")),
    }
    d["state"] = core.STATE_SHORT[core.state_of(d["day"])]
    d["pre"] = bool(r["pre_collapse"])
    if d["pre"] != (d["state"] == core.STATE_SHORT[PRE]):
        raise ValueError(f"{d['day']}: pre_collapse {d['pre']} contradicts "
                         f"{d['state']}")
    return d


def values(rows, key):
    """The non-empty values of one column."""
    return [r[key] for r in rows if r.get(key) is not None]


def stat_block(rows, key):
    """``st._stats`` over one column of the given rows."""
    return st._stats(values(rows, key))


def median_of(rows, key):
    v = values(rows, key)
    return float(np.median(v)) if v else None


def sorted_chips(rows):
    """The chips in acquisition order (day, then orbit label)."""
    return sorted(rows, key=lambda d: (d["day"], str(d["orbit"] or "")))


def ratio_values(rows, key, base):
    """``value / base`` of one dimension, over the rows carrying it."""
    return [r[key] / base for r in rows if r.get(key) is not None and base]


# ---------------------------------------------------------------------------
# The record as blocks (one entry per step, day, stretch, orbit …)
# ---------------------------------------------------------------------------
def days_of(rows):
    """The acquisition days of a set of chips, in order."""
    return sorted({d["day"] for d in rows})


def step_block(rows):
    """Every acquisition step of a chip set, with its length in days."""
    days = days_of(rows)
    return [{"after": days[i - 1], "before": days[i],
             "days": (day_date(days[i]) - day_date(days[i - 1])).days}
            for i in range(1, len(days))]


def gap_block(rows):
    """The steps longer than the cadence: holes in the record, never bridged."""
    out = []
    for s in step_block(rows):
        if s["days"] > CADENCE_DAYS:
            out.append(dict(s, crosses_event=(day_date(s["after"]) < EVENT
                                              <= day_date(s["before"]))))
    return out


def cadence_bands(rows):
    """The chip pairs one cadence step apart — the only place a line is drawn.

    A line through neighbouring chips that are 12 or 24 days apart would be an
    interpolation over a hole, so the figure leaves those steps unjoined.
    """
    chips = sorted_chips(rows)
    out = []
    for a, b in zip(chips, chips[1:]):
        if (day_date(b["day"]) - day_date(a["day"])).days == CADENCE_DAYS:
            out.append({"after": a["day"], "before": b["day"]})
    return out


def stretch_block(rows):
    """The record's two stretches (the states), measured on the chips."""
    out = []
    for state in (PRE, POST):
        sub = sorted_chips([d for d in rows
                            if d["state"] == core.STATE_SHORT[state]])
        echo = [d for d in sub if d["echo"]]
        days = days_of(sub)
        steps = step_block(sub)
        out.append({
            "state": state, "short": core.STATE_SHORT[state],
            "first_day": days[0], "last_day": days[-1],
            "span_days": day_span(days[0], days[-1]),
            "n_chips": len(sub), "n_days": len(days),
            "n_windows": sum(d["windows"] for d in sub),
            "n_duplicated_windows": sum(d["windows"] - 1 for d in sub),
            "n_with_echo": len(echo), "n_without_echo": len(sub) - len(echo),
            "echo_rate": len(echo) / len(sub),
            "orbits": sorted({d["orbit"] for d in sub}),
            "n_steps": len(steps),
            "n_gap_steps": len([s for s in steps if s["days"] > CADENCE_DAYS]),
            "step_lengths": sorted({s["days"] for s in steps}),
        })
    return out


def day_block(rows):
    """One entry per acquisition day: what was requested and what echoed.

    The keys are the committed timeline's own keys (``n_windows``,
    ``n_unique_chips``, ``n_duplicated``, ``n_with_echo``, ``state``), so the
    same block is comparable, field by field, with the two committed timelines.
    """
    out = []
    for day in days_of(rows):
        sub = [d for d in rows if d["day"] == day]
        out.append({
            "day": day, "state": sub[0]["state"],
            "n_windows": sum(d["windows"] for d in sub),
            "n_unique_chips": len(sub),
            "n_duplicated": sum(d["windows"] - 1 for d in sub),
            "n_with_echo": sum(1 for d in sub if d["echo"]),
            "orbits": [d["orbit"] for d in sub],
            "echo_modes": sorted({d["echo_mode"] for d in sub if d["echo"]}),
        })
    return out


def orbit_block(rows):
    """The coverage record split by orbit: where the echoes sit, per geometry."""
    out = {}
    for orbit in ORBITS:
        sub = [d for d in rows if d["orbit"] == orbit]
        by_state = {}
        for state in (PRE, POST):
            s = [d for d in sub if d["state"] == core.STATE_SHORT[state]]
            k = sum(1 for d in s if d["echo"])
            by_state[core.STATE_SHORT[state]] = {
                "n_chips": len(s), "n_with_echo": k,
                "n_without_echo": len(s) - k,
                "echo_rate": (k / len(s)) if s else 0.0}
        k = sum(1 for d in sub if d["echo"])
        out[orbit] = {"n_chips": len(sub), "n_with_echo": k,
                      "n_without_echo": len(sub) - k,
                      "echo_rate": k / len(sub), "by_state": by_state}
    return out


def echo_mode_block(rows):
    """The mask-geometry classes of the echo chips, per state."""
    modes = sorted({d["echo_mode"] for d in rows if d["echo"]})
    out = {}
    for state in (PRE, POST):
        sub = [d for d in rows
               if d["state"] == core.STATE_SHORT[state] and d["echo"]]
        counts = {m: sum(1 for d in sub if d["echo_mode"] == m) for m in modes}
        out[core.STATE_SHORT[state]] = {
            "n_echo_chips": len(sub), "counts": counts,
            "fractions": {m: (counts[m] / len(sub)) if sub else 0.0
                          for m in modes}}
    return out


def composition_block(rows):
    """The record's composition: chips per month, per state, per section."""
    months = {}
    for d in rows:
        m = d["day"][5:7]
        months.setdefault(m, {})[d["state"]] = months.get(m, {}).get(d["state"], 0) + 1
    sections = {}
    for d in rows:
        s = sections.setdefault(d["section"], {"n": 0, "n_pre": 0, "n_post": 0,
                                              "n_with_echo": 0})
        s["n"] += 1
        s["n_pre" if d["pre"] else "n_post"] += 1
        if d["echo"]:
            s["n_with_echo"] += 1
    return {"months": months, "pre": {"n": len([d for d in rows if d["pre"]])},
            "post": {"n": len([d for d in rows if not d["pre"]])},
            "seasons": {}, "sections": sections}


def payload_group_block(unique):
    """The payload groups of the committed window file, from their own digest.

    One group per unique chip: the digest of its measurement window, the days it
    was requested on (one in this record), the committed windows behind it (the
    byte-identical duplicates) and whether the mask rule holds on it. The block is
    recomputed from ``core.chip_table`` — that is, from the raw committed window
    file, not from the channel table's ``payload_group_size`` column — so the two
    independent tallies of the same duplication meet here.
    """
    out = []
    for u in sorted(unique, key=lambda u: u["digest"]):
        segs = u.get("payload_shared_with_segments") or []
        out.append({"digest": u["digest"],
                    "n_windows": int(u["payload_group_size"]),
                    "days": [u["date"]],
                    "segments": sorted(int(s) for s in segs),
                    "with_echo": u["mask"] is not None})
    return out


def post_chips_needed(n_pre, k_pre, alpha=0.05, nmax=400):
    """Post chips needed for a significant Fisher result if *none* echoes.

    The paper's own power question (``fig_ywf_ocv_paper.py``): the pre-event echo
    rate is known from the record, and this is how many post-event chips it would
    take, all echo-free, to reject "the coverage is unchanged" at ``alpha``. The
    definition is repeated here so the committed number is recomputed rather than
    copied; ``None`` when even ``nmax`` chips would not do it.
    """
    for n in range(1, nmax + 1):
        p = core.fisher_p(k_pre, n_pre - k_pre, 0, n)
        if p is not None and p < alpha:
            return n
    return None


def power_block(cov, alpha=0.05):
    """The coverage table's power note, plus the alphabetic halo of it."""
    pre, post = cov["per_state"]["pre"], cov["per_state"]["post"]
    return {
        "n_pre_chips": pre["n_chips"], "n_pre_with_echo": pre["n_with_echo"],
        "n_post_chips": post["n_chips"], "n_post_with_echo": post["n_with_echo"],
        "echo_rate_pre": pre["echo_rate"], "echo_rate_post": post["echo_rate"],
        "fisher_p_two_sided": cov["fisher"]["p_two_sided"],
        "fisher_alpha": alpha,
        "post_chips_needed_at_zero_echo_for_fisher_p_lt_alpha":
            post_chips_needed(pre["n_chips"], pre["n_with_echo"], alpha),
        "note": POWER_NOTE,
    }


def coverage_by_orbit_numbers(rows):
    """Where the echoes sit over the two pass geometries, recomputed from the CSV.

    The Fisher test of the block compares the two orbits' *echo coverage* over the
    whole record; ``fisher_by_state`` repeats it *within* each state (pre vs post
    inside one orbit), which is the paper's decomposition. Neither is a state
    contrast of the orbits — both orbits straddle the collapse — and the note says
    so.
    """
    per = orbit_block(rows)
    table = [[per[o]["n_with_echo"], per[o]["n_without_echo"]] for o in ORBITS]
    fisher_state = {}
    for o in ORBITS:
        pre, post = per[o]["by_state"]["pre"], per[o]["by_state"]["post"]
        fisher_state[o] = core.fisher_p(pre["n_with_echo"], pre["n_without_echo"],
                                        post["n_with_echo"], post["n_without_echo"])
    return {
        "event": core.COLLAPSE_DATE, "n_orbits": len(ORBITS),
        "orbits": list(ORBITS), "n_unique_chips": len(rows),
        "window_px": WINDOW_PX, "per_orbit": per,
        "fisher": {"table_rows": list(ORBITS), "table_cols": ["echo", "no_echo"],
                   "table": table,
                   "p_two_sided": core.fisher_p(table[0][0], table[0][1],
                                                table[1][0], table[1][1])},
        "fisher_by_state": fisher_state, "note": ORBIT_NOTE,
    }


# ---------------------------------------------------------------------------
# The evidence
# ---------------------------------------------------------------------------
def compute(reference_paper=None):
    """Everything this figure draws, from the committed tables only.

    Inputs: the channel table ``data/ywf/ywf_ocv_channels.csv`` (33 unique chips),
    the committed window file ``ywf_windows_full.txt`` (57 windows, recompiled by
    ``core.chip_table``), the committed mask cache ``ywf_windows_mask_cache.txt``
    and its manifest, the generated ``ywf_ocv_channels_meta.json``, the committed
    paper result ``fig_ywf_ocv_paper.json`` and the committed analysis reference
    ``reference/ywf_ocv_findings.json``.
    """
    t0 = time.time()
    meta = core.load_csv_meta()
    reference = core.load_findings()
    paper = core.load_reference(reference_paper or
                                os.path.join(core.DATA, "fig_ywf_ocv_paper.json"))
    cache = core.load_cache()
    windows, unique = core.chip_table()

    rows = sorted_chips([row_view(r) for r in core.load_csv()])
    for d in rows:
        d["x"] = x_of(d["day"], d["orbit"])

    pre = [d for d in rows if d["pre"]]
    post = [d for d in rows if not d["pre"]]
    echo = [d for d in rows if d["echo"]]
    no_echo = [d for d in rows if not d["echo"]]
    pre_echo = [d for d in pre if d["echo"]]
    post_echo = [d for d in post if d["echo"]]
    days = days_of(rows)
    steps = step_block(rows)
    gaps = gap_block(rows)

    # The drawn ordinate: this dimension's value over its own pre-window median.
    pre_median = {k: median_of(pre_echo, k) for k in FEATURES}
    ratios = {d["id"]: {k: d[k] / pre_median[k] for k in FEATURES}
              for d in echo}
    normalised = {k: [{"day": d["day"], "orbit": d["orbit"], "x": d["x"],
                       "value": d[k], "ratio": ratios[d["id"]][k]}
                      for d in echo] for k in FEATURES}

    cov = coverage_numbers({
        "n_windows": sum(d["windows"] for d in rows),
        "n_pre_chips": len(pre), "n_pre_echo": len(pre_echo),
        "n_post_chips": len(post), "n_post_echo": len(post_echo)})
    blocks = {k: {core.STATE_SHORT[s]: stat_block([d for d in rows
                                                   if d["state"] == core.STATE_SHORT[s]],
                                                  k)
                  for s in (PRE, POST)} for k in FEATURES}
    return {
        "rows": rows, "pre": pre, "post": post, "echo": echo, "no_echo": no_echo,
        "pre_echo": pre_echo, "post_echo": post_echo,
        "chips": rows, "days": days, "steps": steps, "gaps": gaps,
        "bands": cadence_bands(rows),
        "first_day": days[0], "last_day": days[-1],
        "last_pre_day": days_of(pre)[-1], "first_post_day": days_of(post)[0],
        "event_weekday": EVENT.strftime("%A"),
        "span_days": day_span(days[0], days[-1]),
        "n_chips": len(rows), "n_days": len(days),
        "n_windows": sum(d["windows"] for d in rows),
        "n_duplicated_windows": sum(d["windows"] - 1 for d in rows),
        "n_echo": len(echo), "n_no_echo": len(no_echo),
        "n_pre_chips": len(pre), "n_pre_days": len(days_of(pre)),
        "n_pre_echo": len(pre_echo), "n_pre_no_echo": len(pre) - len(pre_echo),
        "n_post_chips": len(post), "n_post_days": len(days_of(post)),
        "n_post_echo": len(post_echo), "n_post_no_echo": len(post) - len(post_echo),
        "echo_rate_pre": len(pre_echo) / len(pre),
        "echo_rate_post": len(post_echo) / len(post),
        "n_gap_steps": len(gaps), "n_cadence_steps": len(steps) - len(gaps),
        "step_lengths": sorted({s["days"] for s in steps}),
        "collapse_step": [g for g in gaps if g["crosses_event"]],
        "stretches": stretch_block(rows), "days_block": day_block(rows),
        "coverage": cov, "power": power_block(cov),
        "coverage_by_orbit": coverage_by_orbit_numbers(rows),
        "echo_modes": echo_mode_block(rows), "composition": composition_block(rows),
        "payload_groups": payload_group_block(unique),
        "by_state": blocks,
        "pre_window": {k: {"median": pre_median[k],
                           "n": len(values(pre_echo, k)),
                           "min": min(ratio_values(pre_echo, k, pre_median[k])),
                           "p25": st._stats(ratio_values(pre_echo, k, pre_median[k]))["p25"],
                           "p75": st._stats(ratio_values(pre_echo, k, pre_median[k]))["p75"],
                           "max": max(ratio_values(pre_echo, k, pre_median[k]))}
                       for k in FEATURES},
        "pre_median": pre_median, "ratios": ratios, "normalised": normalised,
        "covariates": {core.STATE_SHORT[s]: {k: stat_block([d for d in rows
                                                            if d["state"] == core.STATE_SHORT[s]], k)
                                             for k in CONTROLS}
                       for s in (PRE, POST)},
        "n_P_set": meta["n_P_set"],
        "p_reference_chips": sorted({r["P_reference_chips"] for r in cache}),
        "n_windows_file": len(windows), "n_unique_file": len(unique),
        "coherence_unique": len({round(d["coherence"], 12) for d in rows
                                 if d["coherence"] is not None}),
        "meta": meta, "manifest": core.load_manifest(),
        "reference": reference, "paper": paper, "cache": cache,
        "runtime_s": time.time() - t0,
    }


# ---------------------------------------------------------------------------
# The pins: four committed layers fix the same record
# ---------------------------------------------------------------------------
TOL_KEYS = {"echo_rate", "rho", "p", "perm_p", "median", "mean", "std", "p5",
            "p25", "p75", "p95", "p_two_sided", "norm", "value"}


class Pins:
    """Collects pin comparisons (recomputed vs. committed) and reports them."""

    def __init__(self):
        self.checks = []
        self.fails = []

    def add(self, path, got, want, tol):
        try:
            exact = bool(got == want)
        except Exception:                       # pragma: no cover - numpy arrays
            exact = False
        dev = 0.0
        if not exact:
            numeric = (isinstance(got, (int, float))
                       and isinstance(want, (int, float))
                       and not isinstance(got, bool) and not isinstance(want, bool))
            dev = abs(float(got) - float(want)) if numeric else None
        ok = exact or (dev is not None and dev <= tol)
        self.checks.append({"path": path, "kind": "exact" if exact else "tolerance",
                            "abs_deviation": dev, "tolerance": tol, "ok": ok})
        if not ok:
            self.fails.append({"path": path, "recomputed": repr(got)[:200],
                               "committed": repr(want)[:200], "abs_deviation": dev})

    def cmp(self, prefix, got, want, key=None):
        """Recursive comparison: ``want`` dictates the structure and expectation."""
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
            return
        if isinstance(want, (list, tuple)):
            glen = len(got) if hasattr(got, "__len__") else -1
            if not isinstance(got, (list, tuple)) or glen != len(want):
                self.add(prefix, f"<len {glen}>", f"<len {len(want)}>", 0.0)
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
            "n_within_tolerance": sum(1 for c in self.checks
                                      if c["kind"] == "tolerance"),
            "n_failed": len(self.fails),
            "max_abs_deviation": max(devs) if devs else 0.0,
            "float_tolerance": FLOAT_TOL,
            "failures": self.fails[:40],
        }


def coverage_numbers(res):
    """The coverage block of the whole record, from the channel table alone."""
    n_pre, k_pre = res["n_pre_chips"], res["n_pre_echo"]
    n_post, k_post = res["n_post_chips"], res["n_post_echo"]
    p = core.fisher_p(k_pre, n_pre - k_pre, k_post, n_post - k_post)
    per_state = {
        "pre": {"n_chips": n_pre, "n_with_echo": k_pre,
                "n_without_echo": n_pre - k_pre, "echo_rate": k_pre / n_pre},
        "post": {"n_chips": n_post, "n_with_echo": k_post,
                 "n_without_echo": n_post - k_post, "echo_rate": k_post / n_post},
    }
    return {
        "event": core.COLLAPSE_DATE, "window_px": WINDOW_PX,
        # ``n_committed_windows`` is the 57 of the window file; the committed
        # blocks name the 33 *chips* of the block ``n_windows`` (the two are
        # pinned against each other in ``pin_all`` so the quirk stays visible).
        "n_committed_windows": res["n_windows"], "n_unique_chips": n_pre + n_post,
        "per_state": per_state,
        "fisher": {"table_rows": ["pre", "post"], "table_cols": ["echo", "no_echo"],
                   "table": [[k_pre, n_pre - k_pre], [k_post, n_post - k_post]],
                   "p_two_sided": p},
    }


def cache_cross_check(res):
    """The committed mask cache against the channel table, channel by channel.

    ``masks.vector_from_row`` re-verifies the mask *rule* on every cache row
    (``A`` against ``|px|``, the peak against ``max_masked``, the floor against
    ``max(0.3 * peak, 5 * median)``) and derives ``A, D, F, S, gamma2`` from the
    committed coordinates; ``P`` comes from the cache's own reference set. All six
    are compared with the CSV — on this site the cache holds *every* chip of the
    record, echo-bearing or not, so all 33 rows meet the channel table here.
    """
    by_id = {d["id"]: d for d in res["rows"]}
    rows, missing = [], 0
    for d in res["cache"]:
        r = by_id.get(d["mid"])
        if r is None:
            missing += 1
            continue
        v = masks.vector_from_row(d)
        rows.append({
            "mid": d["mid"], "date": d["date"], "orbit": d["orbit"],
            "echo": r["echo"],
            "A": (v["A"], r["A"]), "D": (v["D"], r["D"]), "F": (v["F"], r["F"]),
            "S": (v["S"], r["S"]), "gamma2": (v["gamma2"], r["gamma2"]),
            "P": (d["P"], r["P"]),
            "P_reference_chips": d["P_reference_chips"],
        })
    return {"n_cache_rows": len(res["cache"]), "n_compared": len(rows),
            "n_missing_in_csv": missing, "rows": rows,
            "P_reference_chips": sorted({r["P_reference_chips"] for r in rows})}


# ---------------------------------------------------------------------------
# The pins: four committed layers fix the same record
# ---------------------------------------------------------------------------
def committed_timeline(blob, state=None):
    """The entries of a committed ``timeline`` block, keyed by day."""
    out = {}
    for e in blob["timeline"]:
        if state is None or e.get("state") == state:
            out[e["day"]] = e
    return out


def pin_all(res):
    """Re-derive every number of this axis and check it against the commits."""
    meta, man = res["meta"], res["manifest"]
    paper, ref = res["paper"], res["reference"]
    cov, rows = res["coverage"], res["rows"]
    P = Pins()

    # 1) the record, from four committed artifacts ---------------------------
    P.cmp("record.n_rows", len(rows), meta["n_rows"])
    P.cmp("record.n_rows_paper", len(rows), paper["n_rows"])
    P.cmp("record.n_unique_chips", meta["n_unique_chips"], len(rows))
    P.cmp("record.n_unique_chips_paper", paper["n_unique_chips"], len(rows))
    P.cmp("record.n_unique_chips_file", res["n_unique_file"], len(rows))
    P.cmp("record.n_windows", res["n_windows"], meta["n_windows"])
    P.cmp("record.n_windows_pair", res["n_windows"], res["n_windows_file"])
    P.cmp("record.n_windows_paper", res["n_windows"], paper["n_windows"])
    P.cmp("record.n_duplicated_windows", res["n_duplicated_windows"],
          meta["n_duplicate_windows"])
    P.cmp("record.sum_windows", sum(d["windows"] for d in rows), res["n_windows"])
    P.cmp("record.n_echo", res["n_echo"], meta["n_echo_masks"])
    P.cmp("record.n_echo_paper", res["n_echo"], paper["n_echo_chips"])
    P.cmp("record.n_no_echo", res["n_no_echo"], meta["n_rows_without_echo_mask"])
    P.cmp("record.n_pre_chips", res["n_pre_chips"], meta["n_pre_collapse"])
    P.cmp("record.n_post_chips", res["n_post_chips"], meta["n_post_collapse"])
    P.cmp("record.n_pre_echo", res["n_pre_echo"], cov["per_state"]["pre"]["n_with_echo"])
    P.cmp("record.n_post_echo", res["n_post_echo"],
          cov["per_state"]["post"]["n_with_echo"])
    P.cmp("record.n_pre_plus_post", res["n_pre_chips"] + res["n_post_chips"],
          res["n_chips"])
    P.cmp("record.n_echo_plus_no_echo", res["n_echo"] + res["n_no_echo"],
          res["n_chips"])
    P.cmp("record.n_days", res["n_days"],
          meta["mask_cache_manifest"]["asset_date_groups"])
    P.cmp("record.n_days_ref", res["n_days"], ref["records"]["n_dates"])
    P.cmp("record.ref_dedup", ref["records"]["n_deduplicated_same_payload"],
          res["n_duplicated_windows"])
    P.cmp("record.ref_window_sizes", ref["records"]["window_sizes"],
          {meta["mask_window"]: res["n_chips"]})
    P.cmp("record.ref_segments_present", ref["records"]["sections_present"],
          sorted(res["composition"]["sections"]))
    P.cmp("record.ref_segment_resolved", ref["records"]["segment_resolved"], False)
    P.cmp("record.n_by_orbit_meta", meta["n_by_orbit"],
          {core.ORBIT_SHORT[o]: res["coverage_by_orbit"]["per_orbit"][o]["n_chips"]
           for o in ORBITS})
    P.cmp("record.date_range", meta["date_range"], [res["first_day"], res["last_day"]])
    P.cmp("record.first_day", res["first_day"], "2025-08-06")
    P.cmp("record.last_day", res["last_day"], "2026-04-27")
    P.cmp("record.span_days", res["span_days"],
          (day_date(res["last_day"]) - day_date(res["first_day"])).days + 1)
    P.cmp("record.sections_counts", meta["sections_counts"],
          {k: v["n"] for k, v in res["composition"]["sections"].items()})
    P.cmp("record.n_pre_days", res["n_pre_days"],
          len(committed_timeline(paper, "pre")))
    P.cmp("record.n_post_days", res["n_post_days"],
          len(committed_timeline(paper, "post")))
    P.cmp("record.n_pre_days_literal", res["n_pre_days"], 26)
    P.cmp("record.n_post_days_literal", res["n_post_days"], 6)
    P.cmp("record.p_reference_chips", res["p_reference_chips"], [meta["n_P_set"]])
    P.cmp("record.n_P_set", res["n_P_set"], len(res["echo"]))
    P.cmp("record.mask_window", meta["mask_window"], "7x7")
    P.cmp("record.window_px", ref["window_px"], WINDOW_PX)
    P.cmp("record.echo_mode_cuts", meta["echo_mode_cuts"], ref["echo_mode_cuts"])
    P.cmp("record.echo_mode_cuts_paper", meta["echo_mode_cuts"],
          paper["echo_mode_cuts"])
    for key in ("event", "site", "ocv", "features", "state_labels", "mask_rule"):
        P.cmp(f"ref.{key}", ref[key], paper[key])
    P.cmp("ref.state_short", ref["state_short"], dict(core.STATE_SHORT))
    P.cmp("record.states", meta["states"], ref["state_labels"])
    P.cmp("record.orbits", meta["orbits"], ORBITS)
    P.cmp("record.event_weekday", res["event_weekday"], "Monday")
    P.cmp("record.last_pre_day", res["last_pre_day"] < core.COLLAPSE_DATE, True)
    P.cmp("record.first_post_day", res["first_post_day"] > core.COLLAPSE_DATE, True)
    P.cmp("record.states_pre", sorted({d["state"] for d in res["pre"]}), ["pre"])
    P.cmp("record.states_post", sorted({d["state"] for d in res["post"]}), ["post"])
    P.cmp("record.collapse_inside_record",
          res["first_day"] < core.COLLAPSE_DATE < res["last_day"], True)


    # 2) the mask rule, the guard counts and the committed inputs -------------
    P.cmp("manifest.generator", man["generator"],
          meta["mask_cache_manifest"]["generator"])
    P.cmp("manifest.mask_rule", man["mask_rule"], meta["mask_rule"])
    P.cmp("manifest.n_rows", man["n_rows"], res["n_echo"])
    P.cmp("manifest.n_without_echo", man["n_without_echo"], res["n_no_echo"])
    P.cmp("manifest.n_without_index", man["n_without_index"], 0)
    P.cmp("manifest.n_deduplicated", man["n_deduplicated_same_payload"],
          res["n_duplicated_windows"])
    P.cmp("manifest.segment_resolved", man["segment_resolved"],
          ref["records"]["segment_resolved"])
    P.cmp("manifest.window_size", man["window_size"], 7)
    P.cmp("manifest.window_sizes", man["window_sizes"],
          {meta["mask_window"]: res["n_chips"]})
    P.cmp("manifest.source_windows", man["source_windows"],
          os.path.basename(core.WINDOWS_PATH))
    for k in sorted(man):
        if k not in ("source_windows_sha256",):
            P.cmp(f"meta.mask_cache_manifest.{k}", man[k],
                  meta["mask_cache_manifest"][k])
    P.cmp("csv.n_columns", meta["n_columns"], 81)
    P.cmp("csv.row_unit", meta["row_unit"][:32], "one row per unique chip (asset, ")
    P.cmp("csv.n_coherence_compared", meta["n_coherence_masked_pixels_compared"],
          res["n_echo"])
    # the pipeline column is a *different* generation (see ``PIPELINE_NOTE``):
    # the committed export matches the recomputed ``A`` on none of the 17 chips.
    P.cmp("csv.n_coherence_matched", meta["n_coherence_masked_pixels_matched"], 0)
    P.cmp("csv.n_coherence_differing", meta["n_coherence_masked_pixels_differing"],
          res["n_echo"])
    P.cmp("csv.ocv", meta["ocv"], OCV)
    P.cmp("csv.features", meta["features"], FEATURES)
    P.cmp("csv.state_definition", meta["state_definition"],
          "pre-collapse (healthy) = day < 2026-02-02 <= post-collapse")

    # the committed inputs, by content: the same five files this run read ------
    repo = os.path.abspath(os.path.join(core.DATA, os.pardir, os.pardir))
    srcs = ref["sources"]
    for key, path in (("channels_csv", core.CSV_PATH),
                      ("mask_cache", core.CACHE_PATH),
                      ("mask_cache_manifest", core.MANIFEST_PATH),
                      ("windows", core.WINDOWS_PATH),
                      ("segments", core.SEGMENTS_PATH)):
        P.cmp(f"inputs.{key}.path", os.path.relpath(path, repo), srcs[key]["path"])
        P.cmp(f"inputs.{key}.sha256", core.sha256(path), srcs[key]["sha256"])
    P.cmp("inputs.windows_sha_in_manifest", man["source_windows_sha256"],
          srcs["windows"]["sha256"])

    # 3) the mask vector: every committed cache row against the CSV -----------
    cc = cache_cross_check(res)
    P.cmp("cache.n_cache_rows", cc["n_cache_rows"],
          meta["cache_cross_check"]["n_cache_rows"])
    P.cmp("cache.n_compared", cc["n_compared"], res["n_echo"])
    P.cmp("cache.n_missing_in_csv", cc["n_missing_in_csv"], 0)
    P.cmp("cache.P_reference_chips", cc["P_reference_chips"], [meta["n_P_set"]])
    P.cmp("cache.rows_all_echo", sorted({r["echo"] for r in cc["rows"]}), [True])
    for row in cc["rows"]:
        tag = f"cache[{row['date']}/{row['orbit']}]"
        for k in ("A", "D", "F", "S", "gamma2", "P"):
            P.add(f"{tag}.{k}", row[k][0], row[k][1], FLOAT_TOL)

    # 4) the channel summaries: this figure's blocks against the commits ------
    for state in ("pre", "post"):
        for k in FEATURES:
            P.cmp(f"summary_by_state.{state}.{k}",
                  res["by_state"][k][state], meta["summary_by_state"][state][k])
            P.cmp(f"summary_by_state_paper.{state}.{k}",
                  res["by_state"][k][state], paper["summary_by_state"][state][k])
            P.cmp(f"summary_by_state_ref.{state}.{k}",
                  res["by_state"][k][state], ref["echo_mask_by_state"][state][k])
        P.cmp(f"controls.keys.{state}", sorted(res["covariates"][state]),
              sorted(paper["controls_by_state"][state]))
        for k in sorted(paper["controls_by_state"][state]):
            P.cmp(f"controls_by_state.{state}.{k}", res["covariates"][state][k],
                  paper["controls_by_state"][state][k])


    # 5) the coverage of the record, against the three committed coverage blocks
    for src, blob in (("meta", meta["coverage"]), ("paper", paper["coverage"])):
        P.cmp(f"coverage.{src}.event", cov["event"], blob["event"])
        P.cmp(f"coverage.{src}.window_px", cov["window_px"], blob["window_px"])
        P.cmp(f"coverage.{src}.per_state", cov["per_state"], blob["per_state"])
        P.cmp(f"coverage.{src}.fisher", cov["fisher"], blob["fisher"])
        # the committed block's ``n_windows`` counts the 33 *chips* of the block
        P.cmp(f"coverage.{src}.n_windows_is_chips", cov["n_unique_chips"],
              blob["n_windows"])
    ec = ref["echo_coverage"]
    P.cmp("coverage.ref.per_state", cov["per_state"], ec["per_state"])
    P.cmp("coverage.ref.fisher", cov["fisher"], ec["fisher"])
    P.cmp("coverage.ref.n_unique_chips", cov["n_unique_chips"], ec["n_unique_chips"])
    P.cmp("coverage.n_committed_windows", cov["n_committed_windows"], res["n_windows"])
    # the reference's ``records.n_windows`` (41) is the *cache's* own tally of the
    # same duplication — its 17 rows plus the 24 byte-identical duplicates of the
    # manifest (``fig_ywf_ocv_paper.py``) — not the CSV's payload grouping.
    P.cmp("coverage.ref_records.windows_of_cache_rows",
          res["n_echo"] + res["n_duplicated_windows"], ref["records"]["n_windows"])
    P.cmp("coverage.ref_records.n_unique_chips", res["n_echo"],
          ref["records"]["n_unique_chips"])
    # the CSV's other partition of the 57 committed windows: how many windows sit
    # behind one unique chip (the payload groups the digest of the window names).
    P.cmp("coverage.windows_per_chip",
          {str(k): len([d for d in rows if d["windows"] == k])
           for k in sorted({d["windows"] for d in rows})}, {"1": 9, "2": 24})
    # the 24 byte-identical duplicates, split by the echo flag of the chip they
    # belong to: the two partitions of the same 57 windows are pinned apart here.
    P.cmp("coverage.duplicated_windows_split",
          [sum(d["windows"] - 1 for d in res["echo"]),
           sum(d["windows"] - 1 for d in res["no_echo"])], [13, 11])
    P.cmp("coverage.duplicated_windows_sum",
          sum(d["windows"] - 1 for d in res["echo"])
          + sum(d["windows"] - 1 for d in res["no_echo"]),
          res["n_duplicated_windows"])
    for src, blob in (("paper", paper["power"]), ("ref", ec["power"])):
        P.cmp(f"power.{src}", res["power"], blob)
    P.cmp("power.note_is_note", res["power"]["note"], POWER_NOTE)
    P.cmp("power.post_chips_needed", res["power"][
        "post_chips_needed_at_zero_echo_for_fisher_p_lt_alpha"], 4)
    P.cmp("power.fisher_p_lt_alpha",
          res["power"]["fisher_p_two_sided"] < res["power"]["fisher_alpha"], False)

    cbo = res["coverage_by_orbit"]
    P.cmp("coverage_by_orbit.meta", cbo, meta["coverage_by_orbit"])
    P.cmp("coverage_by_orbit.paper", cbo, paper["coverage_by_orbit"])
    P.cmp("coverage_by_orbit.ref", cbo, ref["echo_coverage_by_orbit"])
    P.cmp("coverage_by_orbit.note_is_note", cbo["note"], ORBIT_NOTE)
    P.cmp("coverage_by_orbit.per_orbit_total",
          sum(cbo["per_orbit"][o]["n_chips"] for o in ORBITS), res["n_chips"])
    P.cmp("coverage_by_orbit.per_orbit_with_echo",
          sum(cbo["per_orbit"][o]["n_with_echo"] for o in ORBITS), res["n_echo"])

    # the channel table's own masked-pixel column: a *different* generation, and
    # the reference says so — the comparison is committed as 17 chips, 0 matches.
    pc = ref["pipeline_column"]
    P.cmp("pipeline_column.name", "coherence_masked_pixels", pc["name"])
    P.cmp("pipeline_column.n_compared", res["n_echo"], pc["n_compared"])
    P.cmp("pipeline_column.column_present_on_every_chip",
          len([d for d in rows if d["coherence_masked_pixels"] is not None]),
          res["n_chips"])
    matched = sum(1 for d in res["echo"]
                  if d["coherence_masked_pixels"] is not None
                  and abs(d["coherence_masked_pixels"] - d["A"]) <= FLOAT_TOL)
    P.cmp("pipeline_column.n_matched", matched, pc["n_matched"])
    for state in ("pre", "post"):
        P.cmp(f"pipeline_column.export_by_state.{state}",
              stat_block([d for d in rows if d["state"] == state],
                         "coherence_masked_pixels"),
              pc["export_by_state"][state])
        P.cmp(f"pipeline_column.recomputed_by_state.{state}",
              stat_block([d for d in res["echo"] if d["state"] == state], "A"),
              pc["recomputed_by_state"][state])
    P.cmp("pipeline_column.note", PIPELINE_NOTE, pc["note"])


    # 6) the committed timeline: one entry per acquisition day, field by field --
    mine_days = {d["day"]: d for d in res["days_block"]}
    for src, blob in (("paper", paper), ("ref", ref)):
        com = committed_timeline(blob)
        P.cmp(f"timeline.{src}.keys", sorted(mine_days), sorted(com))
        P.cmp(f"timeline.{src}.n_days", len(com), res["n_days"])
        P.cmp(f"timeline.{src}.sum_windows",
              sum(e["n_windows"] for e in com.values()), res["n_windows"])
        P.cmp(f"timeline.{src}.sum_unique_chips",
              sum(e["n_unique_chips"] for e in com.values()), res["n_chips"])
        P.cmp(f"timeline.{src}.sum_with_echo",
              sum(e["n_with_echo"] for e in com.values()), res["n_echo"])
        P.cmp(f"timeline.{src}.sum_duplicated",
              sum(e["n_duplicated"] for e in com.values()),
              res["n_duplicated_windows"])
        P.cmp(f"timeline.{src}.pre_days",
              sorted(d for d, e in com.items() if e["state"] == "pre"),
              days_of(res["pre"]))
        P.cmp(f"timeline.{src}.post_days",
              sorted(d for d, e in com.items() if e["state"] == "post"),
              days_of(res["post"]))
        for day in sorted(com):
            P.cmp(f"timeline.{src}[{day}]", mine_days.get(day, "<missing day>"),
                  com[day])
    P.cmp("timeline.two_chip_day",
          sorted(d for d, m in mine_days.items() if m["n_unique_chips"] == 2),
          ["2025-08-12"])
    P.cmp("timeline.two_chip_day_orbits", sorted(mine_days["2025-08-12"]["orbits"]),
          ORBITS)
    P.cmp("timeline.two_chip_day_echo_free",
          mine_days["2025-08-12"]["n_with_echo"], 0)

    # 7) the payload accounting and the mask-geometry classes ------------------
    mine_pg = {g["digest"]: g for g in res["payload_groups"]}
    for src, blob in (("paper", paper), ("ref", ref)):
        com_pg = {g["digest"]: g for g in blob["payload_groups"]}
        P.cmp(f"payload_groups.{src}.digests", sorted(mine_pg), sorted(com_pg))
        P.cmp(f"payload_groups.{src}.n_groups", len(com_pg), res["n_chips"])
        P.cmp(f"payload_groups.{src}.sum_windows",
              sum(g["n_windows"] for g in com_pg.values()), res["n_windows"])
        P.cmp(f"payload_groups.{src}.n_duplicated",
              sum(g["n_windows"] - 1 for g in com_pg.values()),
              res["n_duplicated_windows"])
        P.cmp(f"payload_groups.{src}.n_with_echo",
              sum(1 for g in com_pg.values() if g["with_echo"]), res["n_echo"])
        for digest in sorted(com_pg):
            P.cmp(f"payload_groups.{src}[{digest}]", mine_pg.get(digest, "<missing>"),
                  com_pg[digest])
    P.cmp("payload_groups.windows_of_groups",
          sorted(g["n_windows"] for g in res["payload_groups"]),
          sorted(d["windows"] for d in rows))
    P.cmp("payload_groups.days_are_the_chip_days",
          sorted(d for g in res["payload_groups"] for d in g["days"]),
          sorted(d["day"] for d in res["rows"]))
    P.cmp("payload_groups.every_group_one_day",
          sorted({len(g["days"]) for g in res["payload_groups"]}), [1])
    P.cmp("payload_groups.group_with_echo",
          sum(1 for g in res["payload_groups"] if g["with_echo"]), res["n_echo"])
    P.cmp("payload_groups.group_without_echo",
          sum(1 for g in res["payload_groups"] if not g["with_echo"]), res["n_no_echo"])
    P.cmp("echo_modes.paper", res["echo_modes"], paper["echo_modes"])
    P.cmp("echo_modes.ref", res["echo_modes"], ref["echo_modes_by_state"])
    P.cmp("echo_modes.n_echo_chips",
          sum(v["n_echo_chips"] for v in res["echo_modes"].values()), res["n_echo"])
    P.cmp("echo_modes.sum_counts",
          sum(sum(v["counts"].values()) for v in res["echo_modes"].values()),
          res["n_echo"])

    # composition: months, states, sections — the record's own book-keeping ----
    P.cmp("composition.paper", res["composition"], paper["composition"])
    P.cmp("composition.months_ref", res["composition"]["months"],
          ref["months_by_state"])
    P.cmp("composition.sections_ref", res["composition"]["sections"],
          ref["sections"])
    P.cmp("composition.n_pre", res["composition"]["pre"]["n"], res["n_pre_chips"])
    P.cmp("composition.n_post", res["composition"]["post"]["n"], res["n_post_chips"])
    P.cmp("composition.months_total",
          sum(v for m in res["composition"]["months"].values() for v in m.values()),
          res["n_chips"])
    P.cmp("composition.months_present", sorted(res["composition"]["months"]),
          ["01", "02", "03", "04", "08", "09", "10", "11", "12"])
    P.cmp("composition.months_are_the_chip_months",
          sorted(f"{d['day'][5:7]}" for d in res["rows"]),
          sorted(m for m in res["composition"]["months"]
                 for _ in range(sum(res["composition"]["months"][m].values()))))


    # 8) the steps of the record: a hole is a measured step, never a bridge ----
    P.cmp("steps.n", len(res["steps"]), res["n_days"] - 1)
    P.cmp("steps.after_before",
          [[s["after"], s["before"]] for s in res["steps"]],
          [[a, b] for a, b in zip(res["days"], res["days"][1:])])
    P.cmp("steps.lengths_are_days",
          [s["days"] for s in res["steps"]],
          [(day_date(b) - day_date(a)).days for a, b in zip(res["days"],
                                                            res["days"][1:])])
    P.cmp("steps.sum_days", sum(s["days"] for s in res["steps"]),
          (day_date(res["last_day"]) - day_date(res["first_day"])).days)
    P.cmp("steps.histogram", {str(k): len([s for s in res["steps"] if s["days"] == k])
                              for k in sorted({s["days"] for s in res["steps"]})},
          {"6": 22, "12": 7, "24": 2})
    P.cmp("steps.lengths", res["step_lengths"], [CADENCE_DAYS, 12, 24])
    P.cmp("steps.min", min(s["days"] for s in res["steps"]), CADENCE_DAYS)
    P.cmp("steps.max", max(s["days"] for s in res["steps"]), 24)
    P.cmp("steps.n_cadence", res["n_cadence_steps"], 22)
    P.cmp("steps.n_gap", res["n_gap_steps"], 9)
    P.cmp("steps.n_12", len([s for s in res["gaps"] if s["days"] == 12]), 7)
    P.cmp("steps.n_24", len([s for s in res["gaps"] if s["days"] == 24]), 2)
    P.cmp("steps.gaps", [[s["after"], s["before"], s["days"]] for s in res["gaps"]],
          [["2025-08-12", "2025-08-24", 12], ["2025-08-24", "2025-09-05", 12],
           ["2026-01-03", "2026-01-15", 12], ["2026-01-15", "2026-01-27", 12],
           ["2026-01-27", "2026-02-08", 12], ["2026-02-08", "2026-02-20", 12],
           ["2026-02-20", "2026-03-16", 24], ["2026-03-16", "2026-04-09", 24],
           ["2026-04-15", "2026-04-27", 12]])
    P.cmp("steps.gaps_crossing_event",
          [s["after"] for s in res["gaps"] if s["crosses_event"]], ["2026-01-27"])
    P.cmp("steps.collapse_step",
          [[s["after"], s["before"], s["days"]] for s in res["collapse_step"]],
          [["2026-01-27", "2026-02-08", 12]])
    # the collapse step is the committed timeline's own step: its two days are the
    # largest committed day before the event and the smallest one after it.
    com_days = sorted(committed_timeline(paper))
    P.cmp("steps.collapse_step_from_timeline",
          res["collapse_step"][0]["after"],
          max(d for d in com_days if d < core.COLLAPSE_DATE))
    P.cmp("steps.collapse_step_before_from_timeline",
          res["collapse_step"][0]["before"],
          min(d for d in com_days if d >= core.COLLAPSE_DATE))
    P.cmp("steps.n_cadence_bands", len(res["bands"]), 22)
    P.cmp("steps.band_days",
          sorted({(day_date(b["before"]) - day_date(b["after"])).days
                  for b in res["bands"]}), [CADENCE_DAYS])
    P.cmp("steps.band_days_are_chip_days",
          [[b["after"], b["before"]] for b in res["bands"]],
          [[a["day"], b["day"]] for a, b in zip(res["chips"], res["chips"][1:])
           if (day_date(b["day"]) - day_date(a["day"])).days == CADENCE_DAYS])
    P.cmp("record.no_chip_on_event", core.COLLAPSE_DATE in res["days"], False)
    P.cmp("record.last_pre_day_literal", res["last_pre_day"], "2026-01-27")
    P.cmp("record.first_post_day_literal", res["first_post_day"], "2026-02-08")
    P.add("record.collapse_distances",
          [(day_date(core.COLLAPSE_DATE) - day_date(res["last_pre_day"])).days,
           (day_date(res["first_post_day"]) - day_date(core.COLLAPSE_DATE)).days],
          [6, 6], 0.0)

    st_pre, st_post = res["stretches"][0], res["stretches"][1]
    for stretch, s in ((PRE, st_pre), (POST, st_post)):
        short = core.STATE_SHORT[stretch]
        P.cmp(f"stretches.{short}.short", s["short"], short)
        P.cmp(f"stretches.{short}.n_chips", s["n_chips"],
              len(res["pre"]) if stretch == PRE else len(res["post"]))
        P.cmp(f"stretches.{short}.n_days", s["n_days"],
              res["n_pre_days"] if stretch == PRE else res["n_post_days"])
        P.cmp(f"stretches.{short}.n_with_echo", s["n_with_echo"],
              cov["per_state"][short]["n_with_echo"])
        P.cmp(f"stretches.{short}.n_without_echo", s["n_without_echo"],
              cov["per_state"][short]["n_without_echo"])
        P.cmp(f"stretches.{short}.echo_rate", s["echo_rate"],
              cov["per_state"][short]["echo_rate"])
        P.cmp(f"stretches.{short}.orbits", s["orbits"], ORBITS)
        P.cmp(f"stretches.{short}.window_sum",
              sum(d["windows"] for d in (res["pre"] if stretch == PRE else res["post"])),
              s["n_windows"])
    P.cmp("stretches.windows_sum_to_record",
          st_pre["n_windows"] + st_post["n_windows"], res["n_windows"])
    P.cmp("stretches.pre_span", [st_pre["first_day"], st_pre["last_day"]],
          [res["first_day"], res["last_pre_day"]])
    P.cmp("stretches.post_span", [st_post["first_day"], st_post["last_day"]],
          [res["first_post_day"], res["last_day"]])


    # 9) the drawn ordinate: value / pre-window median, and the pin contract ----
    for k in FEATURES:
        base = res["pre_median"][k]
        pre_r = ratio_values(res["pre_echo"], k, base)
        all_r = [p["ratio"] for p in res["normalised"][k]]
        P.add(f"ordinate.{k}.base_is_the_pre_echo_median", base,
              median_of(res["pre_echo"], k), FLOAT_TOL)
        P.cmp(f"ordinate.{k}.base_n", res["pre_window"][k]["n"], res["n_pre_echo"])
        P.cmp(f"ordinate.{k}.n_points", len(all_r), res["n_echo"])
        P.cmp(f"ordinate.{k}.one_point_per_echo_chip",
              sorted(p["value"] for p in res["normalised"][k]),
              sorted(values(res["echo"], k)))
        P.cmp(f"ordinate.{k}.x_days",
              sorted(p["day"] for p in res["normalised"][k]),
              sorted(d["day"] for d in res["echo"]))
        P.cmp(f"ordinate.{k}.pre_side_n",
              len([p for p in res["normalised"][k]
                   if p["day"] < core.COLLAPSE_DATE]), res["n_pre_echo"])
        P.cmp(f"ordinate.{k}.post_side_n",
              len([p for p in res["normalised"][k]
                   if p["day"] >= core.COLLAPSE_DATE]), res["n_post_echo"])
        # the base is the median of the 16 pre-event echo chips, so the drawn
        # pre side has median 1.0 and its quartiles are the committed quartiles
        # divided by that median: every drawn ordinate is a pinned number.
        P.add(f"ordinate.{k}.median_of_pre_ratios", float(np.median(pre_r)), 1.0,
              FLOAT_TOL)
        P.add(f"ordinate.{k}.p25_is_the_committed_p25_over_base",
              res["pre_window"][k]["p25"],
              meta["summary_by_state"]["pre"][k]["p25"] / base, FLOAT_TOL)
        P.add(f"ordinate.{k}.p75_is_the_committed_p75_over_base",
              res["pre_window"][k]["p75"],
              meta["summary_by_state"]["pre"][k]["p75"] / base, FLOAT_TOL)
        P.cmp(f"ordinate.{k}.min_max_are_the_pre_ratios",
              [res["pre_window"][k]["min"], res["pre_window"][k]["max"]],
              [min(pre_r), max(pre_r)])
        P.add(f"ordinate.{k}.post_reading_is_the_committed_post_median_over_base",
              [p["ratio"] for p in res["normalised"][k]
               if p["day"] >= core.COLLAPSE_DATE],
              [meta["summary_by_state"]["post"][k]["median"] / base], FLOAT_TOL)
        P.cmp(f"ordinate.{k}.ratios_inside_the_panel_band",
              [min(all_r) >= PANEL_BAND[0], max(all_r) <= PANEL_BAND[1]], [True, True])
    P.cmp("floor.n_echo_free", res["n_no_echo"], 16)
    P.cmp("floor.days", days_of(res["no_echo"]),
          sorted(set(res["days"]) - set(days_of(res["echo"]))))
    P.cmp("floor.every_echo_free_chip_has_an_x",
          all(d["x"] is not None for d in res["no_echo"]), True)
    P.cmp("floor.non_echo_chips_carry_no_ordinate",
          [k for k in FEATURES if any(d[k] is None for d in res["no_echo"])],
          FEATURES)

    # the panel conventions and the frame of this figure ----------------------
    P.cmp("panels.order", PANELS, ["gamma2", "A", "P", "D", "F", "S"])
    P.cmp("panels.every_dimension_drawn", sorted(PANELS), sorted(FEATURES))
    P.cmp("panels.linear_dims", [k for k in PANELS if PANEL_SCALE[k] == "linear"],
          ["P"])
    P.cmp("panels.log_dims", [k for k in PANELS if PANEL_SCALE[k] == "log"],
          ["gamma2", "A", "D", "F", "S"])
    P.cmp("panels.ocv_in_panel_order", [k for k in PANELS if k in OCV], OCV)
    P.cmp("panels.grid", [PANEL_GRID[0], PANEL_GRID[1]], [3, 2])
    P.cmp("panels.ratio_band", PANEL_BAND, [0.12, 2.91])
    P.cmp("figure.name", FIG_NAME, "fig_ywf_ocv_series")
    P.cmp("figure.json", os.path.relpath(OUT_JSON, repo),
          "data/ywf/fig_ywf_ocv_series.json")
    P.cmp("figure.report", os.path.relpath(OUT_MD, repo),
          "figures/ywf/fig_ywf_ocv_series.md")
    P.cmp("figure.png", os.path.relpath(os.path.join(FIGDIR, FIG_NAME + ".png"),
                                        repo),
          "figures/ywf/fig_ywf_ocv_series.png")
    P.cmp("figure.dpi", float(plt.rcParams["savefig.dpi"]), 600.0)
    P.cmp("figure.width_mm", COL_WIDTH / MM, 183.0)

    # the paper's own A-E pin contract: this figure adds nothing to it ---------
    P.cmp("paper_contract.n_failed", paper["pins"]["n_failed"], 0)
    P.cmp("paper_contract.float_tolerance", paper["pins"]["float_tolerance"],
          FLOAT_TOL)
    P.cmp("paper_contract.n_checks", paper["pins"]["n_checks"], 1179)
    P.cmp("paper_contract.figures", paper["figures"],
          [f"figures/ywf/fig_ywf_ocv_paper_{letter}.png" for letter in "ABCDE"])
    P.cmp("paper_contract.this_figure_is_not_in_it",
          any(FIG_NAME in f for f in paper["figures"]), False)
    return P


# ---------------------------------------------------------------------------
# The figure: one dated axis, the whole record, six dimensions in a 3 x 2 grid
# ---------------------------------------------------------------------------
def dtm(day):
    """Midnight of an ISO day as a ``datetime`` (the x axis carries datetimes)."""
    return dt.datetime.combine(day_date(day), dt.time())


def shift_dt(day, days):
    """``days`` (fractional is fine) after an ISO day or a ``datetime``."""
    base = day if isinstance(day, dt.datetime) else dtm(day)
    return base + dt.timedelta(days=days)


def panel_ylim_span(values, scale):
    """The y window of one panel (the family's own padding rule)."""
    vals = [v for v in values if v is not None and v > 0]
    lo, hi = min(vals), max(vals)
    if scale == "log":
        return lo / 2.2, hi * 2.2
    pad = 0.10 * (hi - lo) if hi > lo else 0.10
    return lo - pad, hi + pad


def panels_block(res):
    """The JSON's description of the six panels (every drawing convention of them)."""
    dims = {}
    for i, dim in enumerate(PANELS):
        points = res["normalised"][dim]
        ratios = [p["ratio"] for p in points]
        vals = [p["value"] for p in points]
        lo, hi = panel_ylim_span(ratios + [1.0], PANEL_SCALE[dim])
        dims[dim] = {
            "letter": "abcdef"[i], "scale": PANEL_SCALE[dim],
            "title": core.CHANNEL_TITLE[dim], "base": res["pre_median"][dim],
            "ordinate": "value / pre-window median", "n_points": len(ratios),
            "ratio_min": min(ratios), "ratio_max": max(ratios),
            "value_min": min(vals), "value_max": max(vals), "ylim": [lo, hi],
            "reference_line": 1.0, "pre_window": res["pre_window"][dim],
        }
    return {
        "grid": list(PANEL_GRID), "band": list(PANEL_BAND),
        "scale_rule": SCALE_RULE, "ordinate": "value / pre-window median",
        "linear_dims": [k for k in PANELS if PANEL_SCALE[k] == "linear"],
        "log_dims": [k for k in PANELS if PANEL_SCALE[k] == "log"],
        "dot_colour": "pass geometry (ascending / descending)",
        "floor_ticks": "the echo-free chips, which carry no ordinate at all",
        "dimensions": dims,
    }


def add_months(day, months):
    """The first of the month ``months`` calendar months after ``day``'s month."""
    month = day.month - 1 + months
    return dt.date(day.year + month // 12, month % 12 + 1, 1)


def month_ticks(res, every=2):
    """The axis' month starts, every ``every`` months, each with its own label.

    The record's own months, not a calendar rule: the ladder starts at the first
    acquisition month and stops before the last one, and the label carries the
    year whenever the year changes — so a reader can place any dot without
    reading the subtitle, and five labels carry a 265-day axis.
    """
    first = dt.date.fromisoformat(res["first_day"])
    last = dt.date.fromisoformat(res["last_day"])
    out, month, year = [], first.replace(day=1), None
    while month <= last:
        if month.year != year:
            year = month.year
            out.append((month, f"{month:%b}\n{year}"))
        else:
            out.append((month, f"{month:%b}"))
        month = add_months(month, every)
    return out


def date_axis(ax, res):
    """A dated x axis: the record's own month starts, every second month.

    Not a calendar rule but the record's months, so the ladder starts at the
    first acquisition month and never runs past the last: the label carries the
    year whenever the year changes, and every second month in between.
    """
    ticks = month_ticks(res)
    ax.set_xticks([dtm(d) for d, _ in ticks], [lab for _, lab in ticks])
    ax.xaxis.set_minor_locator(mdates.MonthLocator())
    ax.tick_params(axis="x", which="minor", length=1.6)
    ax.tick_params(axis="x", which="major", length=2.6)


def draw_record_frame(ax, res):
    """The dated frame every panel shares: two stretches, nine holes, one cut.

    The pre-event stretch is tinted faintly and the post-event stretch more
    strongly — it is the region one chip speaks about — every step longer than the
    cadence is a grey band, because a hole is drawn and never interpolated over,
    and the collapse date is a dashed cut inside the 12-day step that holds it.
    """
    first, last = res["first_day"], res["last_day"]
    ax.set_xlim(shift_dt(first, -5), shift_dt(last, 5))
    ax.axvspan(dtm(first), shift_dt(res["last_pre_day"], 0.5),
               facecolor=STATE_COLOUR[PRE], alpha=PRE_TINT_ALPHA, lw=0.0, zorder=0.4)
    ax.axvspan(shift_dt(res["first_post_day"], -0.5), dtm(last),
               facecolor=STATE_COLOUR[POST], alpha=POST_TINT_ALPHA, lw=0.0,
               zorder=0.4)
    for gap in res["gaps"]:
        ax.axvspan(dtm(gap["after"]), dtm(gap["before"]), facecolor=GAP_COLOR,
                   lw=0.0, zorder=0.5)
    ax.axvline(EVENT, color=CUT_COLOR, lw=0.9, ls=(0, (2.5, 1.6)), zorder=1.6)


def draw_dimension(ax, res, dim):
    """One dimension of the record: cadence lines, orbit dots, the post dot.

    A line joins two chips only when they are one cadence step apart
    (``res["bands"]``), so it breaks at every hole and at the collapse; a band
    whose endpoint day carries no echo has no ordinate to join and stays broken.
    Dots are coloured by pass geometry, the single post-event echo is drawn hollow
    and the 16 echo-free chips are ticks on the panel floor.
    """
    points = {p["day"]: p for p in res["normalised"][dim]}
    for band in res["bands"]:
        a, b = points.get(band["after"]), points.get(band["before"])
        if a and b:
            ax.plot([a["x"], b["x"]], [a["ratio"], b["ratio"]], "-", lw=0.8,
                    color=LINE_COLOUR, alpha=0.85, zorder=2)
    for orbit in ORBITS:
        sel = [p for p in res["normalised"][dim]
               if p["day"] < core.COLLAPSE_DATE and str(p["orbit"]) == orbit]
        ax.plot([p["x"] for p in sel], [p["ratio"] for p in sel], "o", ms=3.0,
                color=ORBIT_COLOUR[orbit], markeredgewidth=0, zorder=4)
    post = [p for p in res["normalised"][dim]
            if p["day"] >= core.COLLAPSE_DATE]
    if post:
        ax.plot([p["x"] for p in post], [p["ratio"] for p in post], "o", ms=4.6,
                mfc="none", mec=POST_COLOUR, mew=1.1, zorder=5)


def panel_ylim(ax, values, scale):
    """The panel's own window. 1.0 is the pre side's median and is drawn."""
    lo, hi = panel_ylim_span(list(values) + [1.0], scale)
    if scale == "log":
        ax.set_yscale("log")
    ax.set_ylim(lo, hi)
    ax.axhline(1.0, color=CUT_COLOR, lw=0.8, zorder=1)
    ax.set_ylabel("value / pre-window median")


def set_y_ticks(ax, lo, hi, scale):
    """Few, plain y labels: this band is narrower than one decade.

    On a log axis matplotlib then labels *every* minor tick (``0.1``, ``0.2``,
    ``0.3`` …), which is unreadable at a ratio band of ``0.12 .. 2.91``. The
    ladder below is the family's own: the same few positions in every log panel,
    written as plain numbers, so nothing is read off a decade exponent.
    """
    if scale != "log":
        return
    ticks = [t for t in LOG_TICKS if lo < t < hi]
    ax.set_yticks(ticks, [f"{t:g}" for t in ticks])
    ax.yaxis.set_minor_locator(NullLocator())


def floor_ticks(ax, res, label=None):
    """The echo-free chips: ticks on the panel floor, no ordinate invented."""
    lo, hi = ax.get_ylim()
    for d in res["no_echo"]:
        ax.plot([d["x"], d["x"]], [lo, lo + 0.075 * (hi - lo)], "-", lw=0.7,
                color=NO_ECHO_COLOR, zorder=3,
                label=label if d is res["no_echo"][0] else None)


def fig_series(res, out_dir=None):
    """The 3 x 2 companion: one dated axis, the whole record, six dimensions.

    Every panel gets the identical treatment — the dated frame of the record then
    its own dimension — because all six dimensions exist on every chip of the
    record; the only difference between panels is the y scale and the label.
    """
    fig, axes = plt.subplots(PANEL_GRID[0], PANEL_GRID[1],
                             figsize=(COL_WIDTH, COL_WIDTH * 0.68))
    flat = list(np.ravel(axes))
    for letter, dim, ax in zip("abcdef", PANELS, flat):
        draw_record_frame(ax, res)
        draw_dimension(ax, res, dim)
        panel_ylim(ax, [p["ratio"] for p in res["normalised"][dim]],
                   PANEL_SCALE[dim])
        set_y_ticks(ax, *ax.get_ylim(), PANEL_SCALE[dim])
        floor_ticks(ax, res)
        date_axis(ax, res)
        panel(ax, letter, f"{core.CHANNEL_TITLE[dim]} ({PANEL_SCALE[dim]})")
    for ax in flat[:-2]:
        ax.tick_params(labelbottom=False)
    for ax in flat[-2:]:
        ax.set_xlabel("acquisition day")
    # the two readings of the cut, drawn once and not six times: a rotated label
    # inside the 12-day step of the bottom row, and the sentence that says what
    # that band is, in the panel with the room for it.
    step = res["collapse_step"][0]
    mid = dtm(step["after"]) + (dtm(step["before"]) - dtm(step["after"])) / 2
    for ax in flat[-2:]:
        ax.text(mid, 0.5, f"collapse {core.COLLAPSE_DATE}", rotation=90,
                transform=ax.get_xaxis_transform(), ha="center", va="center",
                fontsize=5.8, color="#5F5F5F")
    flat[-2].text(dtm(step["after"]), 0.985,
                  f"the event falls inside this {step['days']}-day step",
                  transform=flat[-2].get_xaxis_transform(), ha="right", va="top",
                  fontsize=6.0, color="#5F5F5F")
    handles = [Line2D([], [], marker="o", ls="", ms=3.4, color=ORBIT_COLOUR[o],
                      label=f"{o.title()} pass, pre-event") for o in ORBITS]
    handles.append(Line2D([], [], marker="o", ls="", ms=4.6, mfc="none",
                          mec=POST_COLOUR, mew=1.1,
                          label=f"post-event echo ({res['n_post_echo']} of "
                                f"{res['n_post_chips']} post chips)"))
    handles.append(Line2D([], [], marker="|", ls="", ms=7.0, color=NO_ECHO_COLOR,
                          label=f"echo-free chip, on the panel floor "
                                f"({res['n_no_echo']})"))
    # the key is a figure-level strip under the subtitle: inside a panel it either
    # covers the record or spills over the neighbouring axes, because this axis is
    # 265 days wide and the six panels are full of dots.
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.004, 0.972),
               ncol=2, frameon=False, fontsize=7.2, handletextpad=0.4,
               columnspacing=1.0, borderaxespad=0.0, labelspacing=0.28)
    fig.suptitle("YWF — the echo-mask vector over the whole record, collapse "
                 "inside the axis", x=0.004, y=1.040, ha="left", fontsize=10.0)
    # two short lines rather than one long one: a single 250-character line at
    # 7 pt is wider than the page and ``savefig.bbox = "tight"`` would widen the
    # canvas to fit it, stranding the grid in the left two thirds.
    fig.text(0.004, 1.012,
             "x = the acquisition day (one dot per chip, the two orbits of a day "
             "separated by a sub-day offset); y = the value over that dimension's\n"
             f"own pre-window median, taken over the {res['n_pre_echo']} pre-event "
             f"echo chips: {res['n_echo']} echo-bearing chips of {res['n_chips']} — "
             f"{res['n_pre_echo']} of {res['n_pre_chips']} pre, "
             f"{res['n_post_echo']} of {res['n_post_chips']} post",
             ha="left", va="top", fontsize=7.0, color="#444444", linespacing=1.4)
    fig.tight_layout(pad=0.5, h_pad=1.1, w_pad=1.3, rect=(0, 0, 1, 0.900))
    return save_fig(fig, FIG_NAME, out_dir)


# ---------------------------------------------------------------------------
# The artifacts: the JSON of every number and the report that reads them
# ---------------------------------------------------------------------------
def fmt(x, nd=3):
    """One cell of a Markdown table."""
    if x is None:
        return "-"
    if isinstance(x, float):
        return f"{x:.{nd}f}" if math.isfinite(x) else "n/a"
    return str(x)


def md_table(header, rows, nd=3):
    """A Markdown table with the family's own number format."""
    out = ["| " + " | ".join(str(h) for h in header) + " |",
           "| " + " | ".join("---" for _ in header) + " |"]
    for r in rows:
        out.append("| " + " | ".join(fmt(c, nd) for c in r) + " |")
    return out


def panel_rows(res):
    """One row per panel: what it draws, and what the committed tables hold."""
    out = []
    for i, dim in enumerate(PANELS):
        w = res["pre_window"][dim]
        post = res["by_state"][dim]["post"]
        reading = [p["ratio"] for p in res["normalised"][dim]
                   if p["day"] >= core.COLLAPSE_DATE]
        out.append(["abcdef"[i], dim, PANEL_SCALE[dim], res["pre_median"][dim],
                    w["n"], w["min"], w["p25"], w["p75"], w["max"], post["n"],
                    reading[0] if reading else None, post["median"]])
    return out


def step_rows(res):
    """The record's steps: each length pooled, the longer ones named."""
    out = []
    for length in res["step_lengths"]:
        sel = [s for s in res["steps"] if s["days"] == length]
        out.append([f"{length} d", len(sel),
                    "the cadence" if length == CADENCE_DAYS else "a hole, drawn as a "
                    "band" if len(sel) else "-",
                    ", ".join(f"{s['after']} -> {s['before']}" for s in sel)])
    return out


def timeline_rows(res):
    """The acquisition days of the record, in acquisition order."""
    return [[d["day"], d["state"], d["n_windows"], d["n_unique_chips"],
             d["n_duplicated"], d["n_with_echo"],
             ", ".join(str(o)[:4] for o in sorted(d["orbits"])),
             ", ".join(d["echo_modes"]) or "-"] for d in res["days_block"]]


def group_rows(res):
    """One row per unique chip: the payload groups of the window file."""
    return [[g["digest"][:8], ", ".join(g["days"]), g["n_windows"],
             "yes" if g["with_echo"] else "no",
             ", ".join(f"section-{s}" for s in g["segments"]) or "-"]
            for g in sorted(res["payload_groups"], key=lambda g: g["days"][0])]


def report(res, pins):
    """The report: what the figure draws, and every number behind it."""
    meta, ref, man = res["meta"], res["reference"], res["manifest"]
    cov, pw = res["coverage"], res["power"]
    pre_st, post_st = res["stretches"][0], res["stretches"][1]
    md = []
    md.append(f"# {core.SITE} — the echo-mask vector over the whole record")
    md.append("")
    md.append(f"**Site** {core.SITE} · **asset** {core.ASSET_NAME} · **event** "
              f"{core.COLLAPSE_DATE} ({res['event_weekday']}, inside the axis, "
              f"inside a step) · **record** {res['first_day']} -> "
              f"{res['last_day']} ({res['span_days']} days) · **chips** "
              f"{res['n_chips']} unique {meta['mask_window']} chips on "
              f"{res['n_days']} acquisition days, {res['n_echo']} carrying an echo")
    md.append("")
    md.append("```")
    md.append("OCV  = [gamma2, P, D]               observability core vector")
    md.append("x    = the acquisition day, one dot per chip (2025-08-06 -> 2026-04-27)")
    md.append("y    = the value over that dimension's own pre-window median, taken")
    md.append(f"       over the {res['n_pre_echo']} pre-event echo chips; 1.0 = that "
              "median")
    mr = meta["mask_rule"]
    md.append(f"mask = intensity >= {mr['peak_frac']:.2f}*peak AND >= "
              f"{mr['median_mult']}*{mr['median_kind']}(intensity), A >= "
              f"{mr['min_n_masked']} (window {meta['mask_window']})")
    md.append(f"step = {CADENCE_DAYS} d apart (the two orbits alternate); "
              f"{res['n_gap_steps']} of {len(res['steps'])} steps are longer and")
    md.append("       are drawn as bands — a hole is never interpolated over")
    md.append(f"post = {res['n_post_echo']} of {res['n_post_chips']} post-event "
              "chips carry an echo, so the post side carries")
    md.append("       one point and no line, no median and no trend")
    md.append("```")
    md.append("")
    md.append("## The record")
    md.append("")
    md.append("| fact | value |")
    md.append("| --- | --- |")
    md.append(f"| chips | **{res['n_chips']}** unique {meta['mask_window']} chips "
              f"behind **{res['n_windows']}** committed windows "
              f"({res['n_duplicated_windows']} byte-identical duplicates of the "
              f"windows file) |")
    md.append(f"| acquisition days | **{res['n_days']}**, {res['first_day']} -> "
              f"{res['last_day']} ({res['span_days']} days): a {CADENCE_DAYS}-day "
              f"cadence with {res['n_gap_steps']} longer steps and "
              f"{res['n_cadence_steps']} cadence steps |")
    md.append(f"| pre-event | **{res['n_pre_chips']} chips on {res['n_pre_days']} "
              f"days**, {pre_st['first_day']} -> {pre_st['last_day']} "
              f"({pre_st['span_days']} d), {res['n_pre_echo']} with an echo "
              f"({100 * res['echo_rate_pre']:.1f} %) |")
    md.append(f"| post-event | **{res['n_post_chips']} chips on {res['n_post_days']} "
              f"days**, {post_st['first_day']} -> {post_st['last_day']} "
              f"({post_st['span_days']} d), {res['n_post_echo']} with an echo "
              f"({100 * res['echo_rate_post']:.1f} %) |")
    md.append(f"| the event falls inside a step | {core.COLLAPSE_DATE} sits in the "
              f"**{res['collapse_step'][0]['days']}-day** step "
              f"{res['collapse_step'][0]['after']} -> "
              f"{res['collapse_step'][0]['before']}; no chip on the event day, so "
              f"the two stretches meet at a band and not at a point |")
    two_orbit = next(d for d in res["days_block"] if d["day"] == "2025-08-12")
    md.append(f"| the only two-orbit day | 2025-08-12 carries both passes "
              f"({', '.join(str(o).title() for o in two_orbit['orbits'])}), and "
              f"**neither** echoes |")
    pref = res["p_reference_chips"]
    md.append(f"| the `P` reference | `P` is a window-wide reference, not a per-chip "
              f"quantity: the committed cache carries "
              f"`P_reference_chips` = "
              f"{', '.join(str(v) for v in pref)} — the same "
              f"{res['n_P_set']} echo-bearing chips of the whole record of which the "
              f"mean pixel frequency inside the own mask is taken |")
    pc = ref["pipeline_column"]
    md.append(f"| the pipeline `coherence_masked_pixels` column | it runs "
              f"{fmt(pc['export_by_state']['pre']['p5'])} .. "
              f"{fmt(pc['export_by_state']['pre']['p95'])} masked pixels over the "
              f"{pc['export_by_state']['pre']['n']} pre-event rows (median "
              f"{fmt(pc['export_by_state']['pre']['median'])}) and agrees with the "
              f"recomputed `A` on **{pc['n_matched']} of {pc['n_compared']}** echo "
              f"chips — see the caveat below |")
    md.append(f"| the constant control `coherence` | "
              f"{fmt(res['covariates']['pre']['coherence']['median'], 12)} on every "
              f"one of the {res['n_chips']} chips ({res['coherence_unique']} distinct "
              f"value): reported as a control, never drawn |")
    md.append("")
    md.append("## The six panels")
    md.append("")
    md.append("Every panel draws the same axis: one dot per echo-bearing chip, y = "
              "that dimension's value over **its own** pre-window median. The "
              "committed quartiles of the pre window are divided by the same median "
              "to make the box, so every drawn ordinate is a pinned number.")
    md.append("")
    md += md_table(["panel", "dimension", "y scale", "base (pre median, raw)",
                    "pre n", "pre min (ratio)", "pre p25 (ratio)",
                    "pre p75 (ratio)", "pre max (ratio)", "post n",
                    "the one post chip (ratio)", "post median (raw)"],
                   panel_rows(res))
    md.append("")
    md.append(f"`P` is the only linear panel: it is the mean pixel frequency inside "
              f"the mask over the {res['n_P_set']} reference chips — a bounded "
              f"frequency, so its ratio is read as a difference. The other five are "
              f"multiplicative magnitudes ({WINDOW_PX} is the window, the denominator "
              f"of `A` and `D`); all six pre windows stay inside "
              f"{PANEL_BAND[0]}..{PANEL_BAND[1]} of 1.0, so no log axis hides a level. "
              f"The centre line of every box is 1.0 by construction: it *is* the "
              f"pre-window median the ordinate is divided by.")
    md.append("")
    md.append("## The steps of the record")
    md.append("")
    md.append(f"Acquisitions come in pairs: {res['n_gap_steps']} of the "
              f"{len(res['steps'])} steps between two consecutive acquisition days "
              f"are longer than the {CADENCE_DAYS}-day cadence. Each of them is "
              f"drawn as a grey band and **no line crosses it** — a line from one "
              f"chip to the next would be an interpolation over a hole this record "
              f"does not have data for.")
    md.append("")
    md += md_table(["step", "how many", "what it is", "the steps"], step_rows(res))
    md.append("")
    md.append("## The two stretches")
    md.append("")
    md += md_table(["stretch", "first day", "last day", "calendar days", "chips",
                    "acquisition days", "committed windows", "duplicates",
                    "with echo", "without echo", "echo rate", "orbits", "steps",
                    "longer steps", "step lengths"],
                   [[s["state"], s["first_day"], s["last_day"], s["span_days"],
                     s["n_chips"], s["n_days"], s["n_windows"],
                     s["n_duplicated_windows"], s["n_with_echo"], s["n_without_echo"],
                     s["echo_rate"], ", ".join(str(o).title() for o in s["orbits"]),
                     s["n_steps"], s["n_gap_steps"],
                     ", ".join(f"{d} d" for d in s["step_lengths"])]
                    for s in res["stretches"]])
    md.append("")
    md.append(f"The two stretches are unequal in every respect: "
              f"**{pre_st['span_days']} days and {pre_st['n_chips']} chips** before "
              f"the collapse against **{post_st['span_days']} days and "
              f"{post_st['n_chips']} chips** after it. The collapse step "
              f"({res['collapse_step'][0]['after']} -> "
              f"{res['collapse_step'][0]['before']}) is the only step of the record "
              f"that holds the event date, and the dashed cut of every panel falls "
              f"inside it: the two stretches of the figure meet at a grey band, "
              f"never at a point, and the nearest chips sit 6 days either side of "
              f"the event date.")
    md.append("")
    md.append("## The coverage: what the collapse left behind")
    md.append("")
    md += md_table(["state"] + [c.replace("_", " ")
                               for c in cov["fisher"]["table_cols"]],
                   [[cov["fisher"]["table_rows"][i]] + cov["fisher"]["table"][i]
                    for i in range(2)])
    md.append("")
    md.append(f"Fisher's exact test on that 2x2 table: **p = "
              f"{cov['fisher']['p_two_sided']:.4f}** "
              f"(two-sided, alpha = {pw['fisher_alpha']}). With "
              f"{pw['n_pre_with_echo']} of {pw['n_pre_chips']} pre-event chips "
              f"carrying an echo and {pw['n_post_with_echo']} of "
              f"{pw['n_post_chips']} post-event chips doing so, the record does not "
              f"reject \"the coverage is unchanged\" — and it cannot, with "
              f"{pw['n_post_chips']} post chips.")
    md.append("")
    md.append(f"**The power question of the record.** If all post-event chips came "
              f"back echo-free, it would take "
              f"**{pw['post_chips_needed_at_zero_echo_for_fisher_p_lt_alpha']} "
              f"post-event chips** for the same test to fall below alpha = "
              f"{pw['fisher_alpha']} against the pre-event rate of "
              f"{100 * pw['echo_rate_pre']:.1f} %. The record has "
              f"{pw['n_post_chips']} chips there, one of which echoes.")
    md.append("")
    md.append(pw["note"][0].upper() + pw["note"][1:])
    md.append("")
    md.append("## Where the echoes sit over the two pass geometries")
    md.append("")
    orb = res["coverage_by_orbit"]
    md += md_table(["orbit", "chips", "with echo", "without echo", "echo rate",
                    "pre chips", "pre echo rate", "post chips", "post echo rate",
                    "Fisher within the orbit"],
                   [[o, orb["per_orbit"][o]["n_chips"],
                     orb["per_orbit"][o]["n_with_echo"],
                     orb["per_orbit"][o]["n_without_echo"],
                     orb["per_orbit"][o]["echo_rate"],
                     orb["per_orbit"][o]["by_state"]["pre"]["n_chips"],
                     orb["per_orbit"][o]["by_state"]["pre"]["echo_rate"],
                     orb["per_orbit"][o]["by_state"]["post"]["n_chips"],
                     orb["per_orbit"][o]["by_state"]["post"]["echo_rate"],
                     orb["fisher_by_state"][o]] for o in ORBITS])
    md.append("")
    md.append(f"Across the whole record the two geometries are indistinguishable on "
              f"coverage (Fisher p = {orb['fisher']['p_two_sided']:.3f}). "
              f"{orb['note'][0].upper() + orb['note'][1:]}.")
    md.append("")
    md.append("## The composition of the record")
    md.append("")
    comp = res["composition"]
    md += md_table(["month", "pre chips", "post chips"],
                   [[MONTHS[int(m) - 1] if m.isdigit() else m,
                     comp["months"][m].get("pre", 0),
                     comp["months"][m].get("post", 0)]
                    for m in sorted(comp["months"])])
    md.append("")
    md.append(f"{comp['pre']['n']} of the {res['n_chips']} chips are pre-event and "
              f"{comp['post']['n']} post-event. The asset carries "
              f"{len(comp['sections'])} labelled section "
              f"({', '.join(f'`{s}`' for s in comp['sections'])}) — a single "
              f"turbine-mast section, so no panel of this figure is a comparison "
              f"across sections.")
    md.append("")
    md.append("## The echo geometries")
    md.append("")
    md += md_table(["state", "echo chips", "compact", "distributed", "intermediate"],
                   [[core.STATE_SHORT[s], res["echo_modes"][core.STATE_SHORT[s]]["n_echo_chips"]]
                    + [res["echo_modes"][core.STATE_SHORT[s]]["counts"][m]
                       for m in ("compact", "distributed", "intermediate")]
                    for s in (PRE, POST)])
    md.append("")
    md.append("## The record, day by day")
    md.append("")
    md.append("Every acquisition day of the record, with what was requested, what "
              "the mask rule found and which pass geometry carried it. This is the "
              "committed timeline (`n_windows`, `n_unique_chips`, `n_duplicated`, "
              "`n_with_echo`, `state`) as the channel table recomputes it.")
    md.append("")
    md += md_table(["day", "state", "committed windows", "unique chips", "duplicates",
                    "with echo", "orbits", "echo geometry"],
                   timeline_rows(res))
    md.append("")
    md.append("## The payload groups of the committed window file")
    md.append("")
    md.append(f"The committed window file holds **{res['n_windows_file']} windows** "
              f"for **{res['n_unique_file']} unique chips**: `payload_group_size` "
              f"counts the byte-identical duplicates of one measurement window. "
              f"{res['n_echo']} of the chips carry an echo and "
              f"{res['n_no_echo']} do not, so the {res['n_duplicated_windows']} "
              f"duplicates split "
              f"{sum(g['n_windows'] - 1 for g in res['payload_groups'] if g['with_echo'])}"
              f" onto echo-bearing chips and "
              f"{sum(g['n_windows'] - 1 for g in res['payload_groups'] if not g['with_echo'])}"
              f" onto echo-free ones. The cache's own tally of the same duplication "
              f"is {res['manifest']['n_deduplicated_same_payload']} rows against "
              f"{res['n_windows']} windows (17 + 24); the reference calls the same "
              f"quantity `n_windows` = {ref['records']['n_windows']}.")
    md.append("")
    md += md_table(["window digest", "day", "windows behind it", "echo",
                    "sections"], group_rows(res))
    md.append("")
    md.append("## The caveats of the record")
    md.append("")
    # the committed note carries its own label (`a different generation: …`), so
    # the report supplies the bold lead and uses the rest of the note verbatim.
    note = PIPELINE_NOTE.split(": ", 1)[1]
    md.append(f"**The pipeline column is a different generation.** "
              f"{note[0].upper()}{note[1:]}. "
              f"The two therefore disagree, and they do: the committed column matches "
              f"the recomputed `A` on **{pc['n_matched']} of {pc['n_compared']}** echo "
              f"chips. Its pre-event median is "
              f"{fmt(pc['export_by_state']['pre']['median'])} masked pixels against "
              f"{fmt(pc['recomputed_by_state']['pre']['median'])} recomputed; the "
              f"figure draws the recomputed `A`, because that is the quantity this "
              f"record's mask rule defines, and pins the committed column beside it "
              f"rather than plotting it.")
    md.append("")
    md.append(f"**A control with no variance.** `coherence` is "
              f"{fmt(res['covariates']['pre']['coherence']['median'], 12)} on all "
              f"{res['n_chips']} chips — {res['coherence_unique']} distinct value — so "
              f"it separates nothing: it is reported, never drawn. All "
              f"{len(CONTROLS)} co-variate controls of the channel table "
              f"({', '.join('`' + c + '`' for c in CONTROLS)}) are pinned in the "
              f"JSON, and one of them (`structural_frequency_hz`) is committed for "
              f"only {sum(1 for d in res['rows'] if d['structural_frequency_hz'] is not None)}"
              f" of the {res['n_chips']} chips, so it is not drawn either.")
    md.append("")
    md.append("## The pins and the provenance")
    md.append("")
    s = pins.summary()
    md.append(f"{s['n_checks']} pin comparisons against five committed layers (the "
              f"channel table's meta, the mask-cache manifest, the committed mask "
              f"cache, the paper result JSON and the analysis reference): "
              f"**{s['n_exact']} exact**, "
              f"**{s['n_within_tolerance']} within {s['float_tolerance']:g}**, "
              f"**{s['n_failed']} failed**. The recomputed number is in `recomputed`, "
              f"the committed one in `committed`, and no artifact is written before "
              f"they agree.")
    md.append("")
    md.append("| committed input | sha256 |")
    md.append("| --- | --- |")
    for key in sorted(ref["sources"]):
        md.append(f"| `{ref['sources'][key]['path']}` | "
                  f"`{ref['sources'][key]['sha256']}` |")
    md.append("")
    md.append(f"Figures: `figures/ywf/{FIG_NAME}.png` (600 dpi). JSON: "
              f"`data/ywf/{FIG_NAME}.json`.")
    md.append("")
    return "\n".join(md)


def json_default(o):
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.floating):
        return float(o)
    if isinstance(o, np.bool_):
        return bool(o)
    return str(o)


def json_safe(o):
    """Non-finite floats become ``None`` (JSON has no NaN literal)."""
    if isinstance(o, float):
        return o if math.isfinite(o) else None
    if isinstance(o, np.floating):
        return json_safe(float(o))
    if isinstance(o, dict):
        return {k: json_safe(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [json_safe(v) for v in o]
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.bool_):
        return bool(o)
    return o


def build_json(res, pins_summary):
    """The result JSON: the whole axis, every pinned number, no raw pixels."""
    meta = res["meta"]
    repo = os.path.abspath(os.path.join(core.DATA, os.pardir, os.pardir))
    return {
        "figure": FIG_NAME,
        "script": "code/ywf/fig_ywf_ocv_series.py",
        "site": core.SITE, "site_short": core.SITE_SHORT,
        "asset_id": core.ASSET_ID, "asset_name": core.ASSET_NAME,
        "generated_by": "code/ywf/fig_ywf_ocv_series.py",
        "title": "YWF — the echo-mask vector over the whole record, the collapse "
                 "inside the axis",
        "role": "companion: the 6D series of every echo-bearing chip of the record "
                "on one dated axis, with the echo-free chips on the panel floor",
        "ocv": list(OCV), "features": list(FEATURES),
        "mask_window": meta["mask_window"], "mask_rule": meta["mask_rule"],
        "axis": {
            "name": "acquisition_day",
            "start": res["first_day"], "end": res["last_day"],
            "event": core.COLLAPSE_DATE, "event_weekday": res["event_weekday"],
            "span_days": res["span_days"],
            "definition": "the acquisition day, one dot per chip; the two orbits of "
                          "one day are separated by a sub-day offset",
            "cadence_days": CADENCE_DAYS,
            "steps": res["steps"], "gap_steps": res["gaps"],
            "collapse_step": res["collapse_step"],
            "note": "the event falls inside a drawn step, so the two stretches meet "
                    "at a band and not at a point",
        },
        "normalisation": {
            "rule": "value / that dimension's own pre-window median",
            "base": "the median of the pre-event echo-bearing chips, per dimension",
            "pre_median": res["pre_median"],
            "reference_line": 1.0,
            "no_interpolation": "the nine steps longer than the cadence are drawn as "
                                "bands and no line crosses them",
        },
        "counts": {
            "n_chips": res["n_chips"], "n_days": res["n_days"],
            "n_windows": res["n_windows"],
            "n_duplicated_windows": res["n_duplicated_windows"],
            "n_echo": res["n_echo"], "n_no_echo": res["n_no_echo"],
            "n_pre_chips": res["n_pre_chips"], "n_pre_days": res["n_pre_days"],
            "n_pre_echo": res["n_pre_echo"], "n_pre_no_echo": res["n_pre_no_echo"],
            "n_post_chips": res["n_post_chips"], "n_post_days": res["n_post_days"],
            "n_post_echo": res["n_post_echo"],
            "n_post_no_echo": res["n_post_no_echo"],
            "echo_rate_pre": res["echo_rate_pre"],
            "echo_rate_post": res["echo_rate_post"],
            "n_cadence_steps": res["n_cadence_steps"],
            "n_gap_steps": res["n_gap_steps"],
            "step_lengths": res["step_lengths"],
        },
        "stretches": res["stretches"], "timeline": res["days_block"],
        "coverage": res["coverage"], "power": res["power"],
        "coverage_by_orbit": res["coverage_by_orbit"],
        "echo_modes": res["echo_modes"], "composition": res["composition"],
        "payload_groups": res["payload_groups"],
        "panels": panels_block(res), "pre_window": res["pre_window"],
        "by_state": res["by_state"], "normalised": res["normalised"],
        "covariates": res["covariates"],
        "cache": {"n_rows": len(res["cache"]), "n_P_set": res["n_P_set"],
                  "P_reference_chips": res["p_reference_chips"]},
        "provenance": {
            "note": "the five committed inputs, by content; the same digests the "
                    "analysis reference names, re-read here and pinned",
            "inputs": {key: {"path": os.path.relpath(path, repo),
                             "sha256": core.sha256(path)}
                       for key, path in
                       (("channels_csv", core.CSV_PATH),
                        ("mask_cache", core.CACHE_PATH),
                        ("mask_cache_manifest", core.MANIFEST_PATH),
                        ("windows", core.WINDOWS_PATH),
                        ("segments", core.SEGMENTS_PATH))},
        },
        "pins": pins_summary,
        # no wall-clock time: the figure is deterministic, so two runs of this
        # script must produce the identical JSON (the runtime stays on stdout).
        "figures": [f"figures/ywf/{FIG_NAME}.png"],
    }


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--figdir", default=FIGDIR)
    ap.add_argument("--json", default=OUT_JSON)
    ap.add_argument("--md", default=OUT_MD)
    ap.add_argument("--no-figures", action="store_true",
                    help="run the cross-checks and write no PNG")
    args = ap.parse_args(argv)

    t0 = time.time()
    print(f"computing — {core.SITE}, the whole record around "
          f"{core.COLLAPSE_DATE} …")
    res = compute()
    print(f"  {res['n_chips']} unique chips on {res['n_days']} acquisition days "
          f"({res['n_echo']} carry an echo: {res['n_pre_echo']} of "
          f"{res['n_pre_chips']} pre, {res['n_post_echo']} of "
          f"{res['n_post_chips']} post), {res['n_windows']} committed windows, "
          f"{res['n_gap_steps']} steps longer than the cadence, "
          f"{res['runtime_s']:.1f} s")

    # the pins are the gate: the five committed inputs are re-read and hashed,
    # and no artifact is written before every comparison agrees.
    pins = pin_all(res)
    summary = pins.summary()
    print(f"pins: {summary['n_checks']} checks, {summary['n_exact']} exact, "
          f"{summary['n_within_tolerance']} within "
          f"{summary['float_tolerance']:g}, {summary['n_failed']} failed")
    if summary["n_failed"]:
        for f in summary["failures"][:20]:
            print("PIN FAILED:", f["path"], "recomputed=", f["recomputed"],
                  "committed=", f["committed"])
        _fail(f"{summary['n_failed']} of {summary['n_checks']} pin comparisons "
              f"disagree with the committed layers")

    os.makedirs(args.figdir, exist_ok=True)
    if not args.no_figures:
        fig_series(res, args.figdir)

    os.makedirs(os.path.dirname(os.path.abspath(args.json)), exist_ok=True)
    with open(args.json, "w") as fh:
        json.dump(json_safe(build_json(res, summary)), fh, indent=1,
                  sort_keys=True, default=json_default)
        fh.write("\n")
    print(f"  {os.path.relpath(args.json)}")

    os.makedirs(os.path.dirname(os.path.abspath(args.md)), exist_ok=True)
    with open(args.md, "w") as fh:
        fh.write(report(res, pins))
    print(f"  {os.path.relpath(args.md)}")

    print(f"PIN VERIFICATION OK  ({time.time() - t0:.1f} s)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
