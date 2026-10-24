#!/usr/bin/env python3
"""WERDER HAVEL - the record's own time series: six channels, one dated axis, the repair drawn.

Figure: ``figures/werder_havel/fig_werder_havel_ocv_series.png`` (+ ``.md``)
Result: ``data/werder_havel/fig_werder_havel_ocv_series.json``

The Havelbrücke Werder (A 10) record is 491 target acquisitions from 2022-07-02 to
2025-12-04, and its one documented state change - the renewal of the joints, of the
guardrails, of the mastic asphalt and of the steel's corrosion protection, begun
2024-07-01 and finished 2025-09-05 (Autobahn GmbH) - divides it into the three strata
the export itself names (``W2_bestand`` / ``W3_sanierung`` / ``W4_post``). This figure
draws the six numeric channels the record carries on **every** acquisition, over the
three months before that repair and the month after it:

  axis        2024-04-01 .. 2025-10-05, 553 days, **241** of the record's 491 target
              acquisitions - one per day, because the record's chain is strictly
              interleaved (four covering pass geometries cycle, one acquisition each);
  pre         31 acquisitions 2024-04-01 .. 2024-06-29, all ``W2_bestand``, all S1A;
  repair      189 acquisitions 2024-07-03 .. 2025-09-04, all ``W3_sanierung``, drawn as
              the band 2024-07-01 .. 2025-09-05 (432 days);
  post        21 acquisitions 2025-09-05 .. 2025-10-05, ``W4_post`` truncated at the
              axis' end (that stratum continues to 2025-12-04).

Panels (3 x 2, 183 mm wide, one shared dated axis, every panel logarithmic):

  a  ``coherence_gamma2``          b  ``phase_scatterer_count``
  c  ``phase_coherence``           d  ``intensity``
  e  ``brightness_ratio``          f  ``coherence``

Every ordinate is that channel over the **pre window's own median** of it - the pinned
base, taken over 31 values for every one of the six (``coherence`` is empty on 8 of the
241 drawn rows, but on none of the pre window's). Four things are drawn and nothing
else: one dot per acquisition, coloured by pass geometry; each geometry's own
pre-window median as a solid rule inside the pre window; the post window's median as a
dashed rule; and the 8 empty rows as ticks on the panel floor. No line joins any two
dots - on a chain whose mean step is 2.3 days and which mixes four geometries, a line
would draw the interleave, not the record.

What the axis says, and what it does not.

  * The pooled post/pre ratio of ``coherence_gamma2`` is 1.182, but the four geometries
    disagree in sign (0.94 / 0.58 / 1.97 / 0.54 per track, and 0.67 pooled over the two
    IW1 tracks against their own reference chains) and the spread *between* them is
    11x inside the pre window alone (their own pre medians are 0.405 .. 4.411 of the
    pooled base). The pooled number is a mixture of geometries, not a level of the
    bridge. This figure therefore draws no trend line, quotes no test and takes no
    side; the per-geometry reading is in the report and in the JSON.
  * Two of the six channels are nearly geometry-stable - the four geometries' own pre
    medians span 1.07x on ``brightness_ratio`` (0.978 .. 1.045 of the pooled base) and
    1.22x on ``phase_scatterer_count`` (0.870 .. 1.059) - and the other four are not:
    1.86x (``phase_coherence``), 2.45x (``coherence``), 3.28x (``intensity``) and
    10.9x (``coherence_gamma2``). A reader who compares dots across the sawtooth is
    reading the geometry.
  * ``brightness_ratio`` is ``intensity`` divided by a *per-geometry constant*
    (20286.607 / 14510.520 / 6495.294 / 9935.879, constant to 1e-9 relative inside each
    geometry on all 241 rows), so panels d and e are one channel rescaled, not two
    measurements: the figure draws both because the export carries both and says so in
    the report instead of dropping one.
  * The record's *coverage* changes inside the band, and this is pinned rather than
    smoothed: S1C joins on 2025-03-30 (the only S1C acquisition of that month), and from
    the following month the axis' monthly count goes from 9-11 to 18-20, because all
    four geometries gain a 6-day revisit step where they had a 12-day one (three of the
    four drop the 12-day step entirely; ``desc_95_iw3`` had one 24-day miss before and
    ``asc_44_iw3`` keeps the 12-day step alongside the 6-day one). 57 of the 241 drawn
    acquisitions are S1C. It is drawn as a dotted vertical and read as a coverage fact,
    never as a level.
  * ``coherence`` is empty on exactly 8 of the 241 rows - the first four of the repair
    stratum and the first four of the post stratum; they are ticks on the panel floor,
    never zeros. ``phase_snr_db`` and ``displacement_los_m`` are empty on all 241 and
    are not drawn. ``coherence_masked_pixels`` is 640 = 10 % of the 6400-px window in
    every row of the record: a definition, not a measurement (never drawn).
  * The record carries 14 rows whose ``shift_y`` is about -174 m where every other row
    is within 3 m of 0 - all of them ``asc_44_iw3``, all S1C, 9 of them inside this axis.
    They are pinned in the JSON, never drawn and never removed: the four co-variate
    controls of the record are reported, not plotted, because they would need a second
    ordinate the plate does not have.

What this figure is not. It does not re-run the package's echo-mask rule: this record's
window is 80 x 80 px and its payload is the one external input of the package (~250 MB
of complex JSON), so ``coherence_gamma2`` drawn here is the **export's own column**,
not the rule's ``gamma2``, and the two are never conflated. It reads no SLC, no window
payload, no network: the committed measurement export, its segment list (digest only)
and the committed AOI reference (which commits the two documented event dates and the
export's own digest), all offline.

Usage:
  python3 code/werder_havel/fig_werder_havel_ocv_series.py
  python3 code/werder_havel/fig_werder_havel_ocv_series.py --figdir /tmp/x --json /tmp/y.json --md /tmp/z.md
  python3 code/werder_havel/fig_werder_havel_ocv_series.py --no-figures
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import math
import os
import statistics as st
import sys
import time
from collections import Counter

import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as mdates  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.ticker import NullLocator  # noqa: E402
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import werder_havel_ocv_masks as masks  # noqa: E402

ROOT = os.path.abspath(os.path.join(HERE, os.pardir, os.pardir))
FIGDIR = os.path.join(ROOT, "figures", "werder_havel")
DATA = os.path.join(ROOT, "data", "werder_havel")
NAME = "fig_werder_havel_ocv_series"
OUT_JSON = os.path.join(DATA, NAME + ".json")
OUT_MD = os.path.join(FIGDIR, NAME + ".md")

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "mathtext.fontset": "stix",
    "font.size": 8.5,
    "axes.titlesize": 9.0,
    "axes.labelsize": 8.5,
    "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5,
    "legend.fontsize": 7.0,
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

# ---------------------------------------------------------------------------
# The axis: the three months before the documented repair start, and the month
# after its documented finish. Both ends are the *record's* dates, not the
# calendar's: 2024-04-01 is the first acquisition of the pre window and
# 2025-10-05 the last acquisition this figure draws (W4_post runs on to
# 2025-12-04 and is truncated here, in the report and in the JSON).
# ---------------------------------------------------------------------------
AXIS_FIRST = "2024-04-01"
AXIS_LAST = "2025-10-05"
PRE_LAST = "2024-06-30"        # the pre window ends the day before the documented start
POST_FIRST = "2025-09-05"      # the documented reopening: the first post acquisition
REPAIR_START = masks.EVENT_START    # 2024-07-01, the documented start (Autobahn GmbH)
REPAIR_END = masks.EVENT_END        # 2025-09-05, the documented finish
STRETCHES = ("pre", "repair", "post")
# The six channels of the record this figure draws, in the panel order of the plate -
# the export's own column names, one dot per acquisition each, no derived channel and
# no composite. ``coherence_gamma2`` is the export's coherence column; it is *not* the
# package's echo-mask ``gamma2`` (see the module docstring).
PANELS = ["coherence_gamma2", "phase_scatterer_count", "phase_coherence",
          "intensity", "brightness_ratio", "coherence"]
PANEL_GRID = (3, 2)
CHANNEL_TITLE = {
    "coherence_gamma2": "gamma2 - coherence_gamma2",
    "phase_scatterer_count": "scatterers - phase_scatterer_count",
    "phase_coherence": "phase coherence - phase_coherence",
    "intensity": "intensity - the window's intensity",
    "brightness_ratio": "brightness ratio - brightness_ratio",
    "coherence": "coherence - the export's coherence",
}
FIG_CHANNEL_TITLE = {k: v.replace("gamma^2", r"$\gamma^2$") for k, v in CHANNEL_TITLE.items()}
CHANNEL_NOTE = {
    "coherence_gamma2": "the pipeline's coherence column of the window. The mask rule "
                        "of this package cannot be re-run offline on a record whose "
                        "payload is external, so this panel draws the export's number "
                        "and never claims it for the rule.",
    "phase_scatterer_count": "a count of the window's phase scatterers. It is the second "
                             "most geometry-stable channel of the six: the four "
                             "geometries' own pre medians span 1.22x (0.870 .. 1.059 of "
                             "the pooled base).",
    "phase_coherence": "the pipeline's phase coherence. Positive on all 241 drawn rows, "
                       "and geometry-dependent.",
    "intensity": "a window power (an unnormalised brightness), so it carries the "
                 "geometry and the sensor as well as the bridge.",
    "brightness_ratio": "a ratio of two powers - bounded, and the most geometry-stable "
                        "channel here (the four geometries' own pre medians span 1.07x). "
                        "It is intensity over a per-geometry constant, so as a time "
                        "series it says exactly what the intensity panel says.",
    "coherence": "a second, differently named coherence column of the export: empty on "
                 "8 of the 241 drawn rows (4 in the repair, 4 in the post stride).",
}
PANEL_SCALE = {k: "log" for k in PANELS}
SCALE_RULE = ("all six channels are drawn on a log ordinate: every one of them is "
              "strictly positive on every drawn row (a squared coherence, a count, a "
              "ratio, a power), so no ordinate of this plate can reach 0 and the "
              "gamma2 panel alone spans four decades.")
LOG_TICKS = (0.003, 0.01, 0.03, 0.1, 0.3, 0.5, 1.0, 3.0, 10.0, 30.0)
# The ratio band every drawn ordinate and every drawn level stays inside: it is the
# pinned extent of the 1438 drawn values (min 0.00207, max 21.72), not a rule of thumb.
PANEL_BAND = (0.002, 22.0)
# The floor ticks of the empty rows, as a fraction of the panel's own y span.
FLOOR_FRAC = 0.075
FLOAT_TOL = 1e-12
CADENCE_DAYS = 4               # steps longer than this are holes in the chain
# How many date pairs of one step length the .md names before it says "and N more";
# the JSON's ``steps`` list carries every pair.
STEP_LIST_MAX = 6
TRACK_ORDER = ["asc_146_iw1", "desc_168_iw1", "asc_44_iw3", "desc_95_iw3"]
# The two sub-swaths carry the geometry figure's hues; the light/dark pair inside a
# hue is the orbit (ascending / descending), so a reader keeps the sub-swath a colour.
TRACK_COLOR = {"asc_146_iw1": "#1f4e79", "desc_168_iw1": "#7FA8CC",
               "asc_44_iw3": "#b5651d", "desc_95_iw3": "#E4A96B"}
REPAIR_COLOR = "#EFE3CB"       # the documented repair span
REPAIR_EDGE = "#8C6D2E"
PRE_TINT = "#DCE6F1"
POST_TINT = "#FBE2CE"
TINT_ALPHA = 0.30
GAP_COLOR = "#DDDDDD"
CUT_COLOR = "#888888"
LINE_COLOUR = "#9A9A9A"
POST_COLOUR = "#C44E52"
EMPTY_COLOR = "#B0B0B0"
S1C_COLOR = "#6A3D9A"
CONTROLS = ["temperature_c", "precipitation_mm", "shift_x", "shift_y"]
CONTROL_TITLE = {"temperature_c": "temperature [degC]",
                 "precipitation_mm": "precipitation [mm]",
                 "shift_x": "shift_x [m]", "shift_y": "shift_y [m]"}
# ---------------------------------------------------------------------------
# Small helpers (the family's own, kept dumb on purpose)
# ---------------------------------------------------------------------------
def _fail(msg):
    """Abort without writing: the artifact must never outrun its evidence."""
    raise SystemExit(msg)


def day_date(day):
    """``YYYY-MM-DD`` (or already a ``date``/``datetime``) -> ``datetime.date``."""
    if isinstance(day, dt.datetime):
        return day.date()
    if isinstance(day, dt.date):
        return day
    return dt.date.fromisoformat(day)


def dtm(day):
    """Midnight of an ISO day (or of a ``date``) as a ``datetime`` (the x axis
    carries datetimes)."""
    return dt.datetime.combine(day_date(day), dt.time())


def shift_dt(day, days):
    """``days`` (fractional is fine) after an ISO day or a ``datetime``."""
    base = day if isinstance(day, dt.datetime) else dtm(day)
    return base + dt.timedelta(days=days)


def day_span(first, last):
    """Inclusive day count of ``[first, last]``."""
    return (day_date(last) - day_date(first)).days + 1


def values(rows, key):
    """The non-empty values of one column, over the given rows."""
    return [r["values"][key] for r in rows if r["values"][key] is not None]


def median_of(rows, key):
    """The median of one column over the given rows (``None`` if nothing is there)."""
    v = values(rows, key)
    return st.median(v) if v else None


def ratios(rows, key, base):
    """``value / base`` of one column, over the rows carrying a value."""
    return [r["values"][key] / base for r in rows if r["values"][key] is not None]


def stat_block(rows, key):
    """``min / p25 / median / p75 / max / n`` of one column, or ``None``."""
    v = sorted(values(rows, key))
    if not v:
        return None
    def q(p):
        i = p * (len(v) - 1)
        lo, hi = math.floor(i), math.ceil(i)
        return v[lo] if lo == hi else v[lo] + (v[hi] - v[lo]) * (i - lo)
    return {"n": len(v), "min": v[0], "p25": q(0.25), "median": st.median(v),
            "p75": q(0.75), "max": v[-1]}


def ratio_stat(rows, key, base):
    """``stat_block`` of one column over its own base, plus the base itself."""
    out = stat_block([{"values": {key: r["values"][key]}} for r in rows
                      if r["values"][key] is not None], key)
    if out is None:
        return None
    for k in ("min", "p25", "median", "p75", "max"):
        out[k] = out[k] / base
    out["base"] = base
    return out


def row_view(r):
    """The columns of one record row this figure draws."""
    out = {"id": r["id"], "date": r["date"], "track": r["track"],
           "orbit": r["orbit_direction"], "subswath": r["subswath"],
           "platform": r["platform"], "state": r["state_window"], "role": r["role"],
           "values": {}, "controls": {}}
    for k in PANELS:
        out["values"][k] = float(r[k]) if r.get(k, "") != "" else None
    for k in CONTROLS:
        out["controls"][k] = float(r[k]) if r.get(k, "") != "" else None
    return out


def sorted_rows(rows):
    """The rows in acquisition order (the chain is one acquisition per day)."""
    return sorted(rows, key=lambda d: (d["date"], d["track"]))
# ---------------------------------------------------------------------------
# The pins: the committed layers fix the same record, and this figure says so
# ---------------------------------------------------------------------------
TOL_KEYS = {"echo_rate", "rho", "p", "median", "mean", "std", "p5", "p25", "p75",
            "p95", "p_two_sided", "norm", "value", "masked_pixels_frac_of_window"}


class Pins:
    """Collects pin comparisons (recomputed vs. committed) and reports them."""

    def __init__(self):
        self.checks = []
        self.fails = []

    def add(self, path, got, want, tol):
        try:
            ok = abs(float(got) - float(want)) <= tol
        except (TypeError, ValueError):
            ok = got == want
        kind = "exact" if tol == 0.0 else "tolerance"
        dev = None
        try:
            dev = abs(float(got) - float(want))
        except (TypeError, ValueError):
            dev = None
        self.checks.append({"path": path, "expected": want, "found": got,
                            "kind": kind, "abs_deviation": dev, "ok": bool(ok)})
        if not ok:
            self.fails.append({"path": path, "expected": want, "found": got})

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
# ---------------------------------------------------------------------------
# The blocks: the axis, the chain, the geometries, the holes
# ---------------------------------------------------------------------------
def stretch_block(rows):
    """One stretch of the axis: its span, its size and what is drawn on it."""
    days = [d["date"] for d in rows]
    return {
        "n": len(rows), "n_days": len(set(days)),
        "first": days[0], "last": days[-1],
        "span_days": day_span(days[0], days[-1]),
        "strata": sorted({d["state"] for d in rows}),
        "tracks": dict(sorted(Counter(d["track"] for d in rows).items())),
        "subswaths": dict(sorted(Counter(d["subswath"] for d in rows).items())),
        "platforms": dict(sorted(Counter(d["platform"] for d in rows).items())),
    }


def step_block(rows):
    """Every step of the chain, marked where it straddles a documented boundary.

    The chain is what it is: the four geometries interleave, so the steps are 1..4
    days, and the one 7-day step is a hole drawn as a grey band and never
    interpolated over. A step is also the only thing a line could join, which is
    why no line is drawn at all.
    """
    out = []
    for a, b in zip(rows, rows[1:]):
        n = (day_date(b["date"]) - day_date(a["date"])).days
        out.append({
            "after": a["date"], "before": b["date"], "days": n,
            "after_track": a["track"], "before_track": b["track"],
            "after_platform": a["platform"], "before_platform": b["platform"],
            "crosses_repair_start": a["date"] < REPAIR_START <= b["date"],
            "crosses_repair_end": a["date"] < REPAIR_END <= b["date"],
        })
    return out


def track_block(rows, base):
    """One pass geometry: its size, its own pre level and its own post level.

    The two levels are what makes this plate readable: the four geometries do not
    share a level (on ``coherence_gamma2`` their own pre medians are 0.405 .. 4.411
    of the pooled base), so a dot is read against its own geometry's rule and never
    across the sawtooth. ``post_over_pre`` is the honest per-geometry change of the
    channel, and it is the number the pooled ratio mixes away.
    """
    out = {}
    for t in TRACK_ORDER:
        sel = [d for d in rows if d["track"] == t]
        pre = [d for d in sel if d["date"] <= PRE_LAST]
        rep = [d for d in sel if PRE_LAST < d["date"] < POST_FIRST]
        post = [d for d in sel if d["date"] >= POST_FIRST]
        rec = {"n": len(sel), "n_days": len({d["date"] for d in sel}),
               "n_pre": len(pre), "n_repair": len(rep), "n_post": len(post),
               "first": sel[0]["date"], "last": sel[-1]["date"],
               "orbit": sel[0]["orbit"], "subswath": sel[0]["subswath"],
               "platforms": dict(sorted(Counter(d["platform"] for d in sel).items()))}
        for k in PANELS:
            b = base[k]
            pre_med, post_med = median_of(pre, k), median_of(post, k)
            rec[k] = {
                "pre_median": pre_med, "post_median": post_med,
                "pre_ratio": (pre_med / b) if pre_med else None,
                "post_ratio": (post_med / b) if post_med else None,
                "post_over_pre": ((post_med / pre_med) if (pre_med and post_med)
                                  else None),
            }
        out[t] = rec
    return out


def holes_block(axis, facts):
    """The record's own holes on the drawn rows: empty cells, empty columns.

    The export writes empty fields and this figure never turns one into a 0: an
    empty ``coherence`` is a missing value and is drawn as a tick on the panel
    floor, and a column that is empty on all 991 rows of the record is not drawn at
    all. ``coherence_masked_pixels`` is carried here because it is the one column
    that is a *definition* (the brightest decile of the window), not a measurement.
    """
    per_channel = {}
    for k in PANELS:
        rows = [{"date": d["date"], "track": d["track"], "state": d["state"],
                 "platform": d["platform"]} for d in axis if d["values"][k] is None]
        per_channel[k] = {"n_empty": len(rows), "rows": rows}
    return {
        "n_rows": len(axis),
        "per_channel": per_channel,
        "empty_on_all_record_rows": sorted(
            k for k in facts["column_names"]
            if facts["empty_fields"].get(k) == facts["rows"]),
        "coherence_empty_in_record": facts["empty_fields"]["coherence"],
        "masked_pixels": facts["masked_pixels"],
        "masked_pixels_frac_of_window": facts["masked_pixels_frac_of_window"],
        "window_px": facts["window_px"],
        "record_rows": facts["rows"],
        "record_columns": facts["columns"],
    }
def control_stat(rows, key):
    """``stat_block`` of one co-variate column over the given rows."""
    sel = [{"values": {key: r["controls"][key]}} for r in rows
           if r["controls"][key] is not None]
    return stat_block(sel, key)


def coverage_block(axis, s1c_date=None):
    """The coverage of the axis: the month's acquisitions, and who flew them.

    The chain is not uniform. S1C joins the constellation on 2025-03-30, and from
    that month the axis' acquisitions per month go from 10 to 18-20 because each
    geometry's own revisit step halves (12 days .. 6 days). It is drawn - a dotted
    vertical - and never read as a level: a level step at a coverage change is the
    one artefact this plate could invent, so the report carries the two step sets
    side by side instead.
    """
    months = {}
    for d in axis:
        m = d["date"][:7]
        e = months.setdefault(m, {"n": 0, "platforms": Counter(), "tracks": Counter(),
                                  "first": d["date"], "last": d["date"]})
        e["n"] += 1
        e["platforms"][d["platform"]] += 1
        e["tracks"][d["track"]] += 1
        e["first"] = min(e["first"], d["date"])
        e["last"] = max(e["last"], d["date"])
    per_month = {m: {"n": e["n"], "first": e["first"], "last": e["last"],
                     "platforms": dict(sorted(e["platforms"].items())),
                     "tracks": dict(sorted(e["tracks"].items()))}
                 for m, e in sorted(months.items())}

    def step_set(sel):
        ds = sorted(d["date"] for d in sel)
        return sorted({(day_date(b) - day_date(a)).days for a, b in zip(ds, ds[1:])})

    per_track = {}
    for t in TRACK_ORDER:
        sel = [d for d in axis if d["track"] == t]
        if s1c_date:
            before = [d for d in sel if d["date"] < s1c_date]
            after = [d for d in sel if d["date"] >= s1c_date]
        else:
            before, after = sel, []
        per_track[t] = {"n": len(sel), "steps_before": step_set(before),
                        "steps_after": step_set(after)}
    s1c = [d["date"] for d in axis if d["platform"] == "S1C"]
    return {
        "per_month": per_month,
        "platforms": dict(sorted(Counter(d["platform"] for d in axis).items())),
        "first_s1c": (s1c[0] if s1c else None),
        "n_s1c": len(s1c),
        "first_s1c_month_count": (per_month[s1c[0][:7]]["n"] if s1c else None),
        "per_track": per_track,
        "days_with_more_than_one": len(axis) - len({d["date"] for d in axis}),
    }


REF_CHAINS = {"asc_146_iw1": ("ref_asc_146_iw1", "refout_asc_146_iw1"),
              "desc_168_iw1": ("ref_desc_168_iw1", "refout_desc_168_iw1")}


def reference_block(meas, axis, base):
    """The record's second layer: the IW1 dates carry their own reference chains.

    Two of the four geometries carry a target window *and* a reference window (and an
    outside-reference one) on the very same dates. On the 125 IW1 days of this axis
    that lets the target be read against its own reference instead of against its own
    past - the one comparison of this record that does not mix pass geometries (both
    windows are measured on the same pass, at different places). Only
    ``coherence_gamma2`` and ``phase_scatterer_count`` are carried (they are positive
    on every row); nothing here is drawn (a second ordinate would be needed) and
    nothing here is a test.
    """
    rows = []
    for d in axis:
        chains = REF_CHAINS.get(d["track"])
        if not chains:
            continue
        rec = {"date": d["date"], "track": d["track"], "stretch": d["stretch"],
               "platform": d["platform"]}
        for name, tid in (("target", d["id"]),
                          ("reference", f"{chains[0]}_{d['date']}"),
                          ("reference_outside", f"{chains[1]}_{d['date']}")):
            r = meas.get(tid)
            for k in ("coherence_gamma2", "phase_scatterer_count"):
                rec[f"{name}_{k}"] = (float(r[k]) if r and r.get(k, "") != "" else None)
        for name in ("reference", "reference_outside"):
            rec[f"ratio_{name}_gamma2"] = (
                rec["target_coherence_gamma2"] / rec[f"{name}_coherence_gamma2"]
                if rec[f"{name}_coherence_gamma2"] else None)
            rec[f"ratio_{name}_scatterers"] = (
                rec["target_phase_scatterer_count"] / rec[f"{name}_phase_scatterer_count"]
                if rec[f"{name}_phase_scatterer_count"] else None)
        rows.append(rec)
    per_stretch = {}
    for s in STRETCHES:
        sel = [r for r in rows if r["stretch"] == s]
        block = {"n": len(sel)}
        for name in ("target", "reference", "reference_outside"):
            for key in ("coherence_gamma2", "phase_scatterer_count"):
                block[f"{name}_{key}"] = median_of_block(sel, f"{name}_{key}")
        for name in ("reference", "reference_outside"):
            for key in ("gamma2", "scatterers"):
                block[f"ratio_{name}_{key}"] = median_of_block(
                    sel, f"ratio_{name}_{key}")
        per_stretch[s] = block
    return {"rows": rows, "per_stretch": per_stretch,
            "n_dates": len(rows),
            "n_dates_pre": sum(1 for r in rows if r["stretch"] == "pre"),
            "n_dates_repair": sum(1 for r in rows if r["stretch"] == "repair"),
            "n_dates_post": sum(1 for r in rows if r["stretch"] == "post"),
            "chains": {t: list(v) for t, v in REF_CHAINS.items()},
            "per_track": dict(sorted(Counter(r["track"] for r in rows).items())),
            "carried_columns": ["coherence_gamma2", "phase_scatterer_count"]}


def median_of_block(rows, key):
    """The median of a plain key of flat dicts (the reference block's own reader)."""
    v = sorted(r[key] for r in rows if r.get(key) is not None)
    return (st.median(v) if v else None)
# ---------------------------------------------------------------------------
# Everything the figure draws, from the committed export only
# ---------------------------------------------------------------------------
def stretch_of(date):
    """Which of the three drawn stretches an acquisition day belongs to."""
    if date <= PRE_LAST:
        return "pre"
    if date < POST_FIRST:
        return "repair"
    return "post"


def compute():
    """The axis, its three stretches, the six channels and the record's own holes.

    Reads the committed measurement export (and, for its digest only, the committed
    segment list) and re-verifies both with ``masks.check``: if the record no longer
    re-derives, this figure refuses to draw. It never opens the ~250 MB window
    payload, never re-runs the echo-mask rule and never reads the network.
    """
    t0 = time.time()
    meas = masks.read_measurements()
    problems, facts = masks.check(meas)
    if problems:
        _fail("the committed record no longer re-verifies, so this figure refuses "
              "to draw it:\n  - " + "\n  - ".join(problems))

    all_rows = [row_view(r) for r in meas.values()]
    rows = sorted_rows([d for d in all_rows if d["role"] == "target"])
    for d in rows:
        d["stretch"] = stretch_of(d["date"])
    axis = [d for d in rows if AXIS_FIRST <= d["date"] <= AXIS_LAST]
    pre = [d for d in axis if d["stretch"] == "pre"]
    repair = [d for d in axis if d["stretch"] == "repair"]
    post = [d for d in axis if d["stretch"] == "post"]
    base = {k: median_of(pre, k) for k in PANELS}

    channels = {}
    for k in PANELS:
        b = base[k]
        pr, re_, po = (ratio_stat(s, k, b) for s in (pre, repair, post))
        channels[k] = {
            "column": k, "title": CHANNEL_TITLE[k], "note": CHANNEL_NOTE[k],
            "panel": "abcdef"[PANELS.index(k)], "scale": PANEL_SCALE[k],
            "base": b, "base_n": len(values(pre, k)),
            "n_drawn": len(values(axis, k)), "n_empty": len(axis) - len(values(axis, k)),
            "pre": pr, "repair": re_, "post": po,
            "post_over_pre_pooled": (po["median"] / pr["median"]),
            "band": list(PANEL_BAND),
        }

    geometries = track_block(axis, base)
    steps = step_block(axis)
    holes = holes_block(axis, facts)
    s1c_record = min((d["date"] for d in all_rows if d["platform"] == "S1C"),
                     default=None)
    coverage = coverage_block(axis, s1c_date=s1c_record)
    reference = reference_block(meas, axis, base)

    controls = {s: {k: control_stat(sel, k) for k in CONTROLS}
                for s, sel in (("pre", pre), ("repair", repair), ("post", post))}
    artefacts = [{"date": d["date"], "track": d["track"], "platform": d["platform"],
                  "shift_x": d["controls"]["shift_x"],
                  "shift_y": d["controls"]["shift_y"]}
                 for d in rows if d["controls"]["shift_y"] is not None
                 and abs(d["controls"]["shift_y"]) > 10.0]
    strata = facts["state_window_ranges"]
    stretches = {}
    for s, sel, name in (("pre", pre, "W2_bestand"), ("repair", repair, "W3_sanierung"),
                         ("post", post, "W4_post")):
        blk = stretch_block(sel)
        blk["stratum"] = {"name": name, **strata[name]}
        blk["stratum_truncated"] = (blk["first"] != strata[name]["first_date"]
                                    or blk["last"] != strata[name]["last_date"])
        stretches[s] = blk

    return {
        "axis": axis, "rows": rows, "pre": pre, "repair": repair, "post": post,
        "channels": channels, "base": base, "geometries": geometries, "steps": steps,
        "holes": holes, "coverage": coverage, "reference": reference,
        "controls": controls, "artefacts": artefacts, "stretches": stretches,
        "axis_facts": {
            "first": axis[0]["date"], "last": axis[-1]["date"],
            "span_days": day_span(axis[0]["date"], axis[-1]["date"]),
            "n": len(axis), "n_days": len({d["date"] for d in axis}),
            "pre_window": [AXIS_FIRST, PRE_LAST],
            "repair_window": [REPAIR_START, REPAIR_END],
            "post_window": [POST_FIRST, AXIS_LAST],
            "repair_span_days": day_span(REPAIR_START, REPAIR_END),
            "n_outside_axis": len(rows) - len(axis),
            "n_before_axis": sum(1 for d in rows if d["date"] < AXIS_FIRST),
            "n_after_axis": sum(1 for d in rows if d["date"] > AXIS_LAST),
            "record_rows": facts["rows"], "record_columns": facts["columns"],
            "record_n_dates": facts["n_dates"], "record_n_tracks": facts["n_tracks"],
            "record_roles": facts["roles"],
            "record_state_windows": facts["state_windows"],
            "record_platforms": facts["platforms"], "record_tracks": facts["tracks"],
            "first_date": facts["first_date"], "last_date": facts["last_date"],
        },
        "s1c_record": s1c_record,
        "record": facts,
        "plate": {
            "grid": list(PANEL_GRID), "band": list(PANEL_BAND),
            "scale_rule": SCALE_RULE, "ordinate": "value / pre-window median",
            "log_ticks": list(LOG_TICKS), "floor_frac": FLOOR_FRAC,
            "dot_colour": "pass geometry (two sub-swaths x two orbits)",
            "levels_drawn": ["each geometry's own pre-window median (solid rule, "
                             "inside the pre window)",
                             "the post window's median (dashed rule)"],
            "lines_between_dots": "none - the four geometries alternate on a 2.3-day "
                                  "mean step, so a line would draw the interleave",
            "drawn_bands": ["the pre window (faint tint)",
                            "the documented repair (the band, edged)",
                            "the post window (faint tint)",
                            "every step longer than the cadence (grey)"],
            "S1C_rule": "a dotted vertical at the platform change, a coverage fact",
            "channels": PANELS, "n_rows_drawn": len(axis) * len(PANELS),
        },
        "runtime_s": time.time() - t0,
    }
# ---------------------------------------------------------------------------
# The plate: one dated axis, six channels, the documented repair as a band
# ---------------------------------------------------------------------------
def save_fig(fig, name, out_dir=None):
    """Write one 600 dpi PNG (the binaries stay out of git)."""
    out_dir = out_dir or FIGDIR
    os.makedirs(out_dir, exist_ok=True)
    png = os.path.join(out_dir, name + ".png")
    fig.savefig(png)
    plt.close(fig)
    print(f"  {os.path.relpath(png, ROOT)}")
    return png


def panel(ax, letter, text):
    ax.set_title(f"({letter}) {text}", loc="left")


def add_months(day, months):
    """The first day of the month ``months`` after the month of ``day``."""
    m = day.month - 1 + months
    return day.replace(year=day.year + m // 12, month=m % 12 + 1, day=1)


def month_ticks(res, every=2):
    """The record's own month starts, every second month, each with its label.

    Not a calendar rule but the axis' months: the ladder starts at the first
    acquisition month and stops before the last, and the label carries the year
    whenever the year changes - so a reader places a dot without the subtitle.
    """
    first, last = day_date(res["axis_facts"]["first"]), day_date(res["axis_facts"]["last"])
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
    """The dated x axis every panel shares: the axis' own month starts."""
    ticks = month_ticks(res)
    ax.set_xticks([dtm(d) for d, _ in ticks], [lab for _, lab in ticks])
    ax.xaxis.set_minor_locator(mdates.MonthLocator())
    ax.tick_params(axis="x", which="minor", length=1.6)
    ax.tick_params(axis="x", which="major", length=2.6)


def draw_frame(ax, res):
    """The dated frame every panel shares: two windows, one band, the one hole.

    The pre window is tinted faintly and the post window more strongly (it is the
    month the reopened bridge is read in), the documented repair is the band
    between them, and the one 7-day step of the chain is a grey band because a
    hole is drawn and never interpolated over. The single dotted vertical is the
    coverage change - S1C joining - and is never read as a level.
    """
    first, last = res["axis_facts"]["first"], res["axis_facts"]["last"]
    ax.set_xlim(shift_dt(first, -6), shift_dt(last, 6))
    ax.axvspan(dtm(first), shift_dt(PRE_LAST, 0.5), facecolor=PRE_TINT,
               alpha=TINT_ALPHA, lw=0.0, zorder=0.4)
    ax.axvspan(shift_dt(POST_FIRST, -0.5), shift_dt(last, 6), facecolor=POST_TINT,
               alpha=TINT_ALPHA, lw=0.0, zorder=0.4)
    ax.axvspan(dtm(REPAIR_START), shift_dt(REPAIR_END, -0.5), facecolor=REPAIR_COLOR,
               edgecolor=REPAIR_EDGE, lw=0.6, zorder=0.5)
    for s in res["steps"]:
        if s["days"] > CADENCE_DAYS:
            ax.axvspan(dtm(s["after"]), dtm(s["before"]), facecolor=GAP_COLOR,
                       lw=0.0, zorder=0.6)
    if res["s1c_record"]:
        ax.axvline(dtm(res["s1c_record"]), color=S1C_COLOR, lw=0.8,
                   ls=(0, (1.2, 1.4)), zorder=1.5)


def set_y_ticks(ax):
    """Few, plain y labels: the band is 22x wide and a gamma2 panel spans four decades.

    On a log axis matplotlib would otherwise label every minor tick (``0.01``,
    ``0.02``, ``0.03`` ...), which is unreadable here. The ladder is the family's
    own constant set of positions, written as plain numbers, so nothing is read off
    a decade exponent.
    """
    lo, hi = ax.get_ylim()
    ticks = [t for t in LOG_TICKS if lo < t < hi]
    ax.set_yticks(ticks, [f"{t:g}" for t in ticks])
    ax.yaxis.set_minor_locator(NullLocator())


def draw_channel(ax, res, key):
    """One channel of the axis: the three levels, the dots, nothing joining them.

    Four things and nothing else: one dot per acquisition, coloured by pass
    geometry; each geometry's *own* pre-window median as a solid rule inside the pre
    window (the four disagree, so the rule is drawn four times, not once); the post
    window's median as a dashed rule; and the rows where the column is empty as
    ticks on the panel floor. No line joins two dots: the geometries alternate on a
    2.3-day mean step, so a line would draw the interleave.
    """
    base = res["base"][key]
    ax.set_yscale("log")
    ax.set_ylim(*PANEL_BAND)
    ax.axhline(1.0, color=CUT_COLOR, lw=0.7, ls=(0, (3.0, 2.0)), zorder=1.2)
    for t in TRACK_ORDER:
        rec = res["geometries"][t]
        med = rec[key]["pre_median"]
        if not med or not rec["n_pre"]:
            continue
        ax.plot([dtm(AXIS_FIRST), shift_dt(PRE_LAST, 0.5)], [med / base] * 2, "-",
                lw=1.1, color=TRACK_COLOR[t], zorder=3)
    post_med = res["channels"][key]["post"]["median"]
    if post_med:
        ax.plot([shift_dt(POST_FIRST, -0.5), shift_dt(AXIS_LAST, 6)], [post_med] * 2,
                ls=(0, (3.2, 1.8)), lw=0.9, color=POST_COLOUR, zorder=3)
    for t in TRACK_ORDER:
        sel = [d for d in res["axis"] if d["track"] == t and d["values"][key] is not None]
        if not sel:
            continue
        ax.plot([dtm(d["date"]) for d in sel], [d["values"][key] / base for d in sel],
                "o", ms=2.8, color=TRACK_COLOR[t], markeredgewidth=0, zorder=4)
    ax.set_ylabel("value / its pre-window median")


def floor_ticks(ax, res, key):
    """The rows where this column is empty: ticks on the panel floor, never a zero.

    They stand on the lower spine in axes coordinates, so the tick means "this
    acquisition has no value here" in every panel at the same height and no
    ordinate is invented for it.
    """
    for d in res["axis"]:
        if d["values"][key] is None:
            x = dtm(d["date"])
            ax.plot([x, x], [0.0, FLOOR_FRAC], "-", lw=0.7, color=EMPTY_COLOR,
                    transform=ax.get_xaxis_transform(), zorder=3)
def legend_handles(res):
    """The one legend of the plate: what a dot, a rule, a band and a tick mean."""
    h = [(Line2D([], [], marker="o", ls="none", ms=3.4, color=TRACK_COLOR[t],
                 markeredgewidth=0),
          f"{t} ({res['geometries'][t]['orbit']}, {res['geometries'][t]['subswath']})")
         for t in TRACK_ORDER]
    h.append((Line2D([], [], color=LINE_COLOUR, lw=1.1),
              "each geometry's own pre-window median"))
    h.append((Line2D([], [], color=POST_COLOUR, lw=0.9, ls=(0, (3.2, 1.8))),
              f"the post window's median ({res['stretches']['post']['n']} acquisitions)"))
    h.append((Line2D([], [], color=EMPTY_COLOR, lw=0.9),
              "the column is empty on that acquisition"))
    h.append((Line2D([], [], color=REPAIR_EDGE, lw=0.7, marker="s", ms=4.4,
                     mfc=REPAIR_COLOR, ls="none"),
              f"the documented repair {REPAIR_START} .. {REPAIR_END}"))
    if res["s1c_record"]:
        h.append((Line2D([], [], color=S1C_COLOR, lw=0.9, ls=(0, (1.2, 1.4))),
                  f"S1C joins the constellation ({res['s1c_record']})"))
    return h


def fig_series(res, out_dir=None):
    """The 3 x 2 plate: one dated axis, six channels, one documented repair band.

    Every panel gets the identical treatment - the frame of the axis, then its own
    channel - because all six channels exist on every one of the 241 drawn
    acquisitions. The only differences between panels are the y label and the
    number of empty rows on the floor.
    """
    fig, axes = plt.subplots(PANEL_GRID[0], PANEL_GRID[1],
                             figsize=(COL_WIDTH, COL_WIDTH * 0.70))
    flat = list(np.ravel(axes))
    for letter, key, ax in zip("abcdef", PANELS, flat):
        draw_frame(ax, res)
        draw_channel(ax, res, key)
        set_y_ticks(ax)
        floor_ticks(ax, res, key)
        date_axis(ax, res)
        panel(ax, letter, f"{FIG_CHANNEL_TITLE[key]} ({PANEL_SCALE[key]})")
    for ax in flat[:-2]:
        ax.tick_params(labelbottom=False)
    for ax in flat[-2:]:
        ax.set_xlabel("acquisition day")
    a = res["axis_facts"]
    mid_band = shift_dt(REPAIR_START, a["repair_span_days"] / 2.0)
    flat[0].text(mid_band, 0.5, f"W3_sanierung: the documented repair, "
                                f"{a['repair_span_days']} days",
                 transform=flat[0].get_xaxis_transform(), ha="center", va="center",
                 fontsize=6.0, color=REPAIR_EDGE, zorder=6)
    handles = legend_handles(res)
    fig.legend([x[0] for x in handles], [x[1] for x in handles], loc="lower center",
               ncol=3, frameon=False, fontsize=6.2, handletextpad=0.45,
               columnspacing=0.9, labelspacing=0.22, bbox_to_anchor=(0.5, -0.006))
    fig.suptitle(f"{res['record']['asset_name']} - the record's own time series\n"
                 f"{a['n']} acquisitions on {a['n_days']} days, {a['span_days']} days, "
                 f"six channels, and the one documented repair as a band", fontsize=9.0)
    fig.tight_layout(rect=(0.0, 0.055, 1.0, 0.975))
    return save_fig(fig, NAME, out_dir)
# ---------------------------------------------------------------------------
# The teeth: what the figure asserts, and never draws
# ---------------------------------------------------------------------------
AOI_JSON = os.path.join(DATA, "reference", "werder_havel_aoi.json")
# The record's own committed fields this figure re-derives and pins, field by field,
# against the committed AOI reference - which is the layer that commits them.
REF_FIELDS = ("rows", "columns", "n_dates", "first_date", "last_date", "n_tracks",
              "tracks", "roles", "state_windows", "state_window_ranges", "window_px",
              "window_size", "masked_pixels", "masked_pixels_frac_of_window")


def self_check(res):
    """Every claim this figure makes, as a boolean - the drawing is downstream of this.

    Two layers: the plate's own contract numbers (241 / 553 / 31 / 189 / 21, the six
    bases, the per-geometry ratios) pinned in-code, and the record's committed fields
    re-derived from the export and pinned against the committed AOI reference. A
    failure is a stop: raise SystemExit before anything is written.
    """
    P = Pins()
    facts, a = res["record"], res["axis_facts"]
    checks = []

    def add(name, ok, detail=""):
        checks.append({"check": name, "ok": bool(ok), "detail": detail})

    def pin(path, got, want, rel=1e-12, floor=1e-9):
        P.add(path, got, want, max(floor, abs(want) * rel))

    # ---- the committed layers, re-derived (never copied) ----
    problems, _ = masks.check()
    add("the committed measurement export re-verifies field by field",
        not problems, f"ocv_masks.check(): {len(problems)} problem(s)" if problems
        else "ocv_masks.check(): no problem")
    with open(AOI_JSON, encoding="utf-8") as fh:
        aoi = json.load(fh)
    P.cmp("aoi.reference.record", {k: facts[k] for k in REF_FIELDS},
          {k: aoi["record"][k] for k in REF_FIELDS})
    P.add("aoi.reference.event.start", masks.EVENT_START, aoi["event"]["start"], 0.0)
    P.add("aoi.reference.event.end", masks.EVENT_END, aoi["event"]["end"], 0.0)
    P.add("record.measurements_sha256", masks._sha256(masks.MEAS_PATH),
          aoi["record"]["measurements_sha256"], 0.0)
    P.add("record.segments_sha256", masks._sha256(masks.SEGMENTS_PATH),
          facts["segments_sha256"], 0.0)

    # ---- the axis and its three stretches ----
    pin("axis.n", a["n"], 241)
    pin("axis.span_days", a["span_days"], 553)
    pin("axis.n_days", a["n_days"], 241)
    P.add("axis.first", a["first"], AXIS_FIRST, 0.0)
    P.add("axis.last", a["last"], AXIS_LAST, 0.0)
    add("the axis is 241 acquisitions over 553 days, one acquisition per day",
        a["n"] == 241 and a["span_days"] == 553 and a["n_days"] == a["n"],
        f"{a['n']} rows on {a['n_days']} distinct days, {a['first']} .. {a['last']} "
        f"({a['span_days']} d, mean step {a['span_days'] / (a['n'] - 1):.2f} d)")
    ns = {s: res["stretches"][s]["n"] for s in STRETCHES}
    P.cmp("stretches.n", ns, {"pre": 31, "repair": 189, "post": 21})
    add("the three stretches partition the axis, in order and without overlap",
        sum(ns.values()) == a["n"] and ns["pre"] + ns["repair"] + ns["post"] == 241,
        " + ".join(f"{s} {ns[s]}" for s in STRETCHES) + f" = {sum(ns.values())}")

    inside, strata_ok = [], True
    for s, name in zip(STRETCHES, ("W2_bestand", "W3_sanierung", "W4_post")):
        rng = facts["state_window_ranges"][name]
        blk = res["stretches"][s]
        strata_ok &= blk["strata"] == [name]
        inside.append(all(rng["first_date"] <= d["date"] <= rng["last_date"]
                          for d in res[s]))
    add("every drawn acquisition is inside its own committed stratum",
        strata_ok and all(inside),
        "; ".join(f"{s} = {res[s][0]['date']} .. {res[s][-1]['date']}" for s in STRETCHES)
        + " (only the repair covers its stratum end to end, so the other two are "
          "marked truncated in the JSON)")
    add("the repair band is drawn on the documented span and the other two on the axis",
        REPAIR_START == masks.EVENT_START and REPAIR_END == masks.EVENT_END
        and res["stretches"]["repair"]["stratum"] == {**{"name": "W3_sanierung"},
                                                      **facts["state_window_ranges"]["W3_sanierung"]},
        f"band {REPAIR_START} .. {REPAIR_END} = {day_span(REPAIR_START, REPAIR_END)} d")
    # ---- the six channels: the base, the drawn extent, the geometry ----
    drawn = {k: values(res["axis"], k) for k in PANELS}
    n_drawn = sum(len(v) for v in drawn.values())
    ratios_all = [v / res["base"][k] for k in PANELS for v in drawn[k]]
    for k in PANELS:
        ch = res["channels"][k]
        pin(f"channels.{k}.base", ch["base"], median_of(res["pre"], k))
        pin(f"channels.{k}.post_over_pre_pooled", ch["post_over_pre_pooled"],
            {"coherence_gamma2": 1.1821100586736395,
             "phase_scatterer_count": 0.9259259259259259,
             "phase_coherence": 1.045145145460116,
             "intensity": 0.9642363436706458,
             "brightness_ratio": 0.9931399305784077,
             "coherence": 0.8376848263344595}[k])
    add("the base of every channel is the pre window's own median, over 31 values",
        all(ch["base_n"] == 31 and ch["pre"]["median"] == 1.0
            for ch in res["channels"].values()),
        "six bases: " + ", ".join(f"{k} {res['base'][k]:g}" for k in PANELS))
    add("every drawn ordinate is positive and inside the pinned band",
        n_drawn == 1438 and min(ratios_all) > PANEL_BAND[0]
        and max(ratios_all) < PANEL_BAND[1] and min(ratios_all) > 0.0,
        f"{n_drawn} drawn values, {min(ratios_all):.4g} .. {max(ratios_all):.4g} of "
        f"their own base, band {PANEL_BAND[0]:g} .. {PANEL_BAND[1]:g}")

    g2 = {t: res["geometries"][t]["coherence_gamma2"] for t in TRACK_ORDER}
    signs = {("up" if g2[t]["post_over_pre"] > 1 else "down") for t in TRACK_ORDER}
    def spread(k):
        return (max(res["geometries"][t][k]["pre_ratio"] for t in TRACK_ORDER)
                / min(res["geometries"][t][k]["pre_ratio"] for t in TRACK_ORDER))
    add("the four geometries disagree on gamma2's post/pre sign and on its pre level",
        signs == {"up", "down"} and spread("coherence_gamma2") > 10.0,
        "post/pre " + ", ".join(f"{g2[t]['post_over_pre']:.3f}" for t in TRACK_ORDER)
        + f"; pre medians span {spread('coherence_gamma2'):.1f}x")
    add("brightness_ratio is the most geometry-stable channel and "
        "phase_scatterer_count the second",
        spread("brightness_ratio") < 1.10 < spread("phase_scatterer_count") < 1.30,
        f"pre-median spread: brightness_ratio {spread('brightness_ratio'):.3f}x, "
        f"phase_scatterer_count {spread('phase_scatterer_count'):.3f}x, "
        f"coherence_gamma2 {spread('coherence_gamma2'):.1f}x")

    ks = {}
    for t in TRACK_ORDER:
        v = [d["values"]["intensity"] / d["values"]["brightness_ratio"]
             for d in res["axis"] if d["track"] == t]
        ks[t] = (min(v), max(v))
    worst = max((b - a) / a for a, b in ks.values())
    add("brightness_ratio is intensity over a per-geometry constant", worst < 1e-6,
        "constants " + " / ".join(f"{(a + b) / 2:.3f}" for a, b in ks.values())
        + f", worst relative spread {worst:.1e}: panels d and e are one channel "
          "rescaled, not two measurements")

    # ---- the holes, the coverage ----
    empty = res["holes"]["per_channel"]
    coh = [d["date"] for d in res["axis"] if d["values"]["coherence"] is None]
    first4 = [res[s][i]["date"] for s in ("repair", "post") for i in range(4)]
    add("coherence is empty on the first four rows of the repair and of the post "
        "stretch, and on nothing else",
        coh == first4 and len(coh) == 8
        and empty["coherence"]["n_empty"] == 8
        and all(empty[k]["n_empty"] == 0 for k in PANELS if k != "coherence"),
        f"{len(coh)} empty coherence rows: {', '.join(coh)}; the other five channels "
        "are empty on 0 of 241")
    add("the two columns empty on every row of the record are not drawn",
        res["holes"]["empty_on_all_record_rows"] == ["displacement_los_m", "phase_snr_db"]
        and res["holes"]["coherence_empty_in_record"] == 24
        and res["holes"]["masked_pixels"] == 640
        and abs(res["holes"]["masked_pixels_frac_of_window"] - 0.1) < 1e-9,
        "empty on all 991 rows: displacement_los_m, phase_snr_db; coherence is empty "
        "on 24 of the record's 991 rows (8 of them inside this axis); "
        "coherence_masked_pixels = 640 = 10 % of the 6400 px window, one value")

    long_steps = [s for s in res["steps"] if s["days"] > CADENCE_DAYS]
    add("the chain has exactly one step longer than the cadence, and it is drawn",
        len(long_steps) == 1 and long_steps[0]["days"] == 7
        and long_steps[0]["after"] == "2024-08-28"
        and {s["days"] for s in res["steps"]} == {1, 2, 3, 4, 7},
        "steps " + ", ".join(f"{s['after']} -> {s['before']} ({s['days']} d)"
                             for s in long_steps) + "; every other step is 1..4 days")
    cov = res["coverage"]
    steps_before = {t: cov["per_track"][t]["steps_before"] for t in TRACK_ORDER}
    steps_after = {t: cov["per_track"][t]["steps_after"] for t in TRACK_ORDER}
    add("the coverage change (S1C joining) is pinned, drawn and never read as a level",
        cov["n_s1c"] == 57 and cov["first_s1c"] == "2025-03-30"
        and all(steps_before[t] == [12] for t in TRACK_ORDER if t != "desc_95_iw3")
        and all(6 in v for v in steps_after.values())
        and sum(1 for v in steps_after.values() if 12 in v) == 1
        and cov["days_with_more_than_one"] == 0,
        f"S1C from {cov['first_s1c']} ({cov['first_s1c_month_count']} acquisitions that "
        f"month), {cov['n_s1c']} of the 241 drawn rows; per-geometry revisit steps "
        + ", ".join(f"{t} {steps_before[t]} -> {steps_after[t]}" for t in TRACK_ORDER)
        + "; no day carries two acquisitions")

    ctrl = {s: {k: res["controls"][s][k]["n"] for k in CONTROLS} for s in STRETCHES}
    n_cells = sum(ctrl[s][k] for s in STRETCHES for k in CONTROLS)
    add("all four co-variate controls are carried on every drawn row and none is drawn",
        n_cells == 4 * 241 and all(ctrl[s][k] == res["stretches"][s]["n"]
                                   for s in STRETCHES for k in CONTROLS),
        f"temperature_c, precipitation_mm, shift_x, shift_y: {n_cells} of 964 cells "
        "filled, reported in the JSON, not plotted")
    art = res["artefacts"]
    add("the shift_y artefacts are reported, not removed",
        len(art) == 14 and {a0["track"] for a0 in art} == {"asc_44_iw3"}
        and {a0["platform"] for a0 in art} == {"S1C"}
        and sum(1 for a0 in art if a0["date"] <= AXIS_LAST) == 9,
        f"{len(art)} of the record's 991 rows carry shift_y ~ -174 m "
        f"({min(a0['shift_y'] for a0 in art):.1f} .. {max(a0['shift_y'] for a0 in art):.1f}), "
        "all asc_44_iw3, all S1C, 9 of them inside this axis; the other 477 rows stay "
        "within 2.21 m")

    # ---- the reference chains of the two IW1 geometries ----
    ref = res["reference"]
    add("the two IW1 geometries carry their own reference chain on the same dates",
        ref["n_dates"] == 125
        and ref["per_track"] == {"asc_146_iw1": 63, "desc_168_iw1": 62}
        and ref["n_dates_pre"] + ref["n_dates_repair"] + ref["n_dates_post"] == 125,
        f"{ref['n_dates']} of the axis' 241 acquisitions carry a target and a reference "
        f"window ({ref['per_track']}; pre {ref['n_dates_pre']}, repair "
        f"{ref['n_dates_repair']}, post {ref['n_dates_post']}), the other 143 do not")
    for s, want in (("pre", (1.3229529389735704, 1.3968218260024237, 16)),
                    ("repair", (1.2260732043806444, 1.2540255328010441, 98)),
                    ("post", (0.6663937765405988, 4.308582342652187, 11))):
        pin(f"reference.{s}.n", ref["per_stretch"][s]["n"], want[2])
        pin(f"reference.{s}.ratio_reference_gamma2",
            ref["per_stretch"][s]["ratio_reference_gamma2"], want[0])
        pin(f"reference.{s}.ratio_reference_outside_gamma2",
            ref["per_stretch"][s]["ratio_reference_outside_gamma2"], want[1])
    post_blk = ref["per_stretch"]["post"]
    pooled = res["channels"]["coherence_gamma2"]["post_over_pre_pooled"]
    add("the same-pass reference reading disagrees with the pooled post/pre ratio",
        post_blk["ratio_reference_gamma2"] < 1.0 < pooled,
        f"target/reference gamma2 in the post window "
        f"{post_blk['ratio_reference_gamma2']:.3f} (pre "
        f"{ref['per_stretch']['pre']['ratio_reference_gamma2']:.3f}, "
        f"outside-reference in the post window "
        f"{post_blk['ratio_reference_outside_gamma2']:.2f}), against the pooled post/pre "
        f"{pooled:.3f} of panel a - neither is drawn and neither is a test, both are in "
        "the JSON")

    summary = P.summary()
    add("every pinned number agrees with the layer it is pinned against",
        summary["n_failed"] == 0,
        f"{summary['n_checks']} comparisons ({summary['n_exact']} exact, "
        f"{summary['n_within_tolerance']} within tolerance): "
        f"{summary['n_failed']} failed, max |deviation| {summary['max_abs_deviation']:.3g}")
    return checks, P
# ---------------------------------------------------------------------------
# The result JSON and the report
# ---------------------------------------------------------------------------
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


def build_json(res, checks, pins_summary):
    """The whole drawn axis, every pinned number, and nothing that is not the export.

    The 241 drawn acquisitions are carried row by row, with the six channel values and
    the four co-variates exactly as the export wrote them (an empty field is ``null``,
    never a 0). No wall-clock time is written: two runs of this script over the same
    committed export must produce the identical JSON (the runtime is printed instead).
    """
    a, ref = res["axis_facts"], res["reference"]
    passed = sum(1 for c in checks if c["ok"])
    return {
        "figure": NAME,
        "script": "code/werder_havel/fig_werder_havel_ocv_series.py",
        "generated_by": NAME,
        "asset": masks.ASSET, "asset_id": masks.ASSET_ID,
        "asset_name": masks.ASSET_NAME, "asset_type": masks.ASSET_TYPE,
        "title": "Werder/Havel - the record's six vector channels over the three months "
                 "before the documented repair and the month after it",
        "role": "companion: every acquisition of the axis, one dot per channel per "
                "acquisition, each channel over the pre window's own median of it",
        "reads": [
            "data/werder_havel/werder_havel_measurements_full.txt - the committed "
            "export, re-verified by werder_havel_ocv_masks.check() before anything is "
            "drawn; its sha256 is pinned",
            "data/werder_havel/werder_havel_segments.txt - the segment list, digest only",
            "data/werder_havel/reference/werder_havel_aoi.json - the committed AOI "
            "reference: the two documented event dates, the export's digest and the "
            "record's own census, re-derived and pinned field by field",
        ],
        "not_read": [
            "the ~250 MB echo-window payload: never opened, so coherence_gamma2 here is "
            "the export's own column and never the mask rule's gamma2",
            "no SLC, no network, no cache",
        ],
        "channels_drawn": list(PANELS),
        "axis": {
            "name": "acquisition_day",
            "first": a["first"], "last": a["last"], "span_days": a["span_days"],
            "n_acquisitions": a["n"], "n_days": a["n_days"],
            "definition": "one dot per acquisition, on the acquisition day; the chain is "
                          "strictly interleaved (four pass geometries, one acquisition "
                          "each) and no day carries two acquisitions",
            "cadence_days": CADENCE_DAYS,
            "record_rows": a["record_rows"], "record_columns": a["record_columns"],
            "n_record_acquisitions": a["record_n_dates"],
            "n_outside_axis": a["n_outside_axis"],
            "n_before_axis": a["n_before_axis"], "n_after_axis": a["n_after_axis"],
            "record_first": a["first_date"], "record_last": a["last_date"],
        },
        "stretches": res["stretches"],
        "channels": res["channels"],
        "geometries": res["geometries"],
        "coverage": res["coverage"],
        "holes": res["holes"],
        "steps": res["steps"],
        "reference": ref,
        "controls": res["controls"],
        "artefacts": res["artefacts"],
        "plate": res["plate"],
        "drawn_rows": res["axis"],
        "self_check": {"n_checks": len(checks), "n_passed": passed, "checks": checks},
        "pins": pins_summary,
        "figures": [f"figures/werder_havel/{NAME}.png"],
    }
def fmt(x, nd=3):
    """One cell of a Markdown table: ``nd`` decimals, never a bare ``0.000``.

    ``coherence_gamma2`` runs at 1e-4 (its pre median is 0.000128534) and a fixed
    three decimals would print the target window, its reference window and the
    outside-reference window of the reference table all as ``0.000`` - the very
    numbers that section exists to report. Magnitudes below 0.01 therefore keep
    six significant digits instead of ``nd`` decimals.
    """
    if x is None:
        return "-"
    if isinstance(x, float):
        if not math.isfinite(x):
            return "n/a"
        if x != 0.0 and abs(x) < 0.01:
            return f"{x:.6g}"
        return f"{x:.{nd}f}"
    return str(x)


def md_table(header, rows, nd=3):
    """A Markdown table with the family's own number format."""
    out = ["| " + " | ".join(str(h) for h in header) + " |",
           "| " + " | ".join("---" for _ in header) + " |"]
    for r in rows:
        out.append("| " + " | ".join(fmt(c, nd) for c in r) + " |")
    return out


def stretch_rows(res):
    """One row per drawn stretch, with the committed stratum it sits in."""
    labels = {"pre": "the three months before the documented start",
              "repair": "the band (the documented repair), 2024-07-01 .. 2025-09-05",
              "post": "the month after the documented finish"}
    out = []
    for s in STRETCHES:
        b, st_ = res["stretches"][s], res["stretches"][s]["stratum"]
        out.append([labels[s], b["first"], b["last"], b["n"], b["span_days"],
                    st_["name"], st_["first_date"] + " .. " + st_["last_date"],
                    "yes" if b["stratum_truncated"] else "no",
                    ", ".join(f"{k} {v}" for k, v in b["platforms"].items()),
                    ", ".join(f"{k} {v}" for k, v in b["tracks"].items())])
    return out


def channel_rows(res):
    """One row per drawn channel: its base, its own pre window and its reading."""
    geo, out = res["geometries"], []
    for i, k in enumerate(PANELS):
        ch = res["channels"][k]
        pr = [geo[t][k]["pre_ratio"] for t in TRACK_ORDER]
        out.append(["abcdef"[i], k, ch["base"], ch["n_drawn"], ch["n_empty"],
                    ch["pre"]["min"], ch["pre"]["median"], ch["pre"]["max"],
                    ch["repair"]["median"], ch["post"]["median"],
                    ch["post_over_pre_pooled"], max(pr) / min(pr)])
    return out


def geometry_rows(res):
    """One row per pass geometry: its size, its own gamma2 rules."""
    out = []
    for t in TRACK_ORDER:
        g = res["geometries"][t]
        out.append([t, g["orbit"], g["subswath"],
                    ", ".join(f"{k} {v}" for k, v in g["platforms"].items()),
                    g["n_pre"], g["n_repair"], g["n_post"], g["first"], g["last"],
                    g["coherence_gamma2"]["pre_ratio"],
                    g["coherence_gamma2"]["post_over_pre"]])
    return out


def geometry_channel_rows(res):
    """The per-geometry post/pre ratio of all six channels, side by side."""
    return [[k] + [res["geometries"][t][k]["post_over_pre"] for t in TRACK_ORDER]
            for k in PANELS]


def coverage_rows(res):
    """One row per month of the axis: the acquisitions and who flew them."""
    out = []
    for m, e in res["coverage"]["per_month"].items():
        out.append([m, e["n"], ", ".join(f"{k} {v}" for k, v in e["platforms"].items()),
                    ", ".join(f"{k} {v}" for k, v in e["tracks"].items()),
                    e["first"], e["last"]])
    return out


def step_rows(res):
    """The chain's steps, pooled by length, with the long pools elided.

    The 3- and 4-day pools hold 66 and 29 steps: spelling every pair out turns one
    cell into a wall of 66 dates that hides the two numbers this table is for (the
    step length and its count). The first ``STEP_LIST_MAX`` pairs are named, the
    count and the pointer to the JSON's own ``steps`` list say what is left out.
    """
    by = {}
    for s in res["steps"]:
        by.setdefault(s["days"], []).append(f"{s['after']} -> {s['before']}")
    out = []
    for d, v in sorted(by.items()):
        shown = ", ".join(v[:STEP_LIST_MAX])
        if len(v) > STEP_LIST_MAX:
            shown += (f", ... and {len(v) - STEP_LIST_MAX} more - all {len(v)} are "
                      f"in the JSON (`steps`)")
        out.append([d, len(v),
                    "the cadence" if d <= CADENCE_DAYS else "a hole, drawn as a band",
                    shown])
    return out


def reference_rows(res):
    """The reference reading per stretch (target / same-pass reference window)."""
    out = []
    for s in STRETCHES:
        b = res["reference"]["per_stretch"][s]
        out.append([s, b["n"], b["target_coherence_gamma2"],
                    b["reference_coherence_gamma2"], b["reference_outside_coherence_gamma2"],
                    b["ratio_reference_gamma2"], b["ratio_reference_outside_gamma2"],
                    b["ratio_reference_scatterers"]])
    return out


def hole_rows(res):
    """One row per drawn channel: how many of the 241 rows are empty, and which."""
    out = []
    for k in PANELS:
        e = res["holes"]["per_channel"][k]
        out.append([k, e["n_empty"],
                    ", ".join(r["date"] for r in e["rows"]) if e["rows"] else "-"])
    return out


def control_rows(res):
    """The four co-variates of the record: medians per stretch, extreme on the axis."""
    out = []
    for k in CONTROLS:
        blocks = [res["controls"][s][k] for s in STRETCHES]
        out.append([CONTROL_TITLE[k]] + [b["median"] for b in blocks]
                   + [min(b["min"] for b in blocks), max(b["max"] for b in blocks)])
    return out


def artefact_rows(res):
    """The rows whose shift_y is two orders of magnitude off the rest, named."""
    return [[a["date"], a["track"], a["platform"], a["shift_x"], a["shift_y"]]
            for a in res["artefacts"]]
def report(res, checks, pins_summary):
    """The ``.md`` beside the plate: the same numbers in text, with their caveats."""
    a, rec, cov = res["axis_facts"], res["record"], res["coverage"]
    ref, plate = res["reference"], res["plate"]
    ns = {s: res["stretches"][s]["n"] for s in STRETCHES}
    n_drawn = sum(res["channels"][k]["n_drawn"] for k in PANELS)
    ratio = {k: res["channels"][k]["post_over_pre_pooled"] for k in PANELS}
    # the per-geometry constant that turns intensity into brightness_ratio
    bconst = {t: min(d["values"]["intensity"] / d["values"]["brightness_ratio"]
                     for d in res["axis"] if d["track"] == t) for t in TRACK_ORDER}
    def gspread(k):
        """The spread of one channel's four geometry-level pre medians."""
        return (max(res["geometries"][t][k]["pre_ratio"] for t in TRACK_ORDER)
                / min(res["geometries"][t][k]["pre_ratio"] for t in TRACK_ORDER))
    L = [f"# {NAME}\n"]
    L.append(f"**Site** {masks.ASSET_NAME} (`{masks.ASSET_ID}`) - Werder/Havel, the A 10 "
             f"crossing of the Havel. **Axis** {a['first']} -> {a['last']} "
             f"({a['span_days']} d, {a['n']} acquisitions, one per day). "
             f"**Stretches** pre {ns['pre']} (2024-04-01 .. 2024-06-29, all "
             f"`W2_bestand`, all S1A), the documented repair {masks.EVENT_START} .. "
             f"{masks.EVENT_END} ({a['repair_span_days']} d drawn, {ns['repair']} "
             f"acquisitions, all `W3_sanierung`), post {ns['post']} (2025-09-05 .. "
             f"2025-10-05, `W4_post` truncated). **Plate** six channels x {a['n']} "
             f"acquisitions, {n_drawn} drawn values, 183 mm, 600 dpi.\n")
    L.append("```")
    L.append("channels = [coherence_gamma2, phase_scatterer_count, phase_coherence,")
    L.append("            intensity, brightness_ratio, coherence]      the export's own columns")
    L.append(f"x        = the acquisition day ({a['first']} -> {a['last']}), one dot per")
    L.append("           acquisition; the four pass geometries alternate, no day carries two")
    L.append(f"y        = value / that channel's own pre-window median (1.0 = the median of")
    L.append(f"           the {ns['pre']} pre acquisitions); all six panels logarithmic")
    L.append(f"base     = " + " / ".join(f"{k} {res['base'][k]:g}" for k in PANELS[:3]))
    L.append(f"           " + " / ".join(f"{k} {res['base'][k]:g}" for k in PANELS[3:]))
    L.append(f"band     = {plate['band'][0]:g} .. {plate['band'][1]:g}: every drawn ordinate and every")
    L.append("           drawn level stays inside it (pinned, not a rule of thumb)")
    L.append("no line  = the mean step is "
             f"{a['span_days'] / (a['n'] - 1):.2f} d over four interleaved geometries, so a")
    L.append("           line would draw the interleave and not the record")
    L.append("event    = the documented repair, drawn as a band and never as a level")
    L.append("```\n")
    L.append("## The axis\n")
    L.extend(md_table(["stretch", "first", "last", "n", "days", "committed stratum",
                       "that stratum's committed span", "stratum truncated?",
                       "platforms", "geometries"],
                      stretch_rows(res)))
    L.append("")
    L.append(f"The record itself is {rec['rows']} rows x {rec['columns']} columns over "
             f"{rec['n_dates']} acquisition dates ({rec['first_date']} .. "
             f"{rec['last_date']}), of which this plate draws the {a['n']} between "
             f"2024-04-01 and 2025-10-05 ({a['n_before_axis']} earlier, "
             f"{a['n_after_axis']} later are not drawn). The chain inside the axis is "
             f"strictly interleaved - the four pass geometries cycle, "
             f"{cov['days_with_more_than_one']} days carry two acquisitions - and the "
             f"three stretches cut on the record's own strata, not on the calendar: "
             f"the documented start {masks.EVENT_START} falls between "
             f"`W2_bestand` and `W3_sanierung` and the documented finish "
             f"{masks.EVENT_END} on the first `W4_post` acquisition "
             f"({masks.EVENT_SOURCE}). Only the band covers its stratum end to end; the "
             f"pre and post stretches are marked truncated in the JSON because both "
             f"strata run past the axis this figure documents.\n")
    L.append("## The six panels\n")
    L.append(f"Each panel draws one column of the export, and every ordinate is that "
             f"column over the **pre window's own median** of it - re-derived here from "
             f"the {ns['pre']} pre acquisitions, never copied. Six bases: "
             + ", ".join(f"`{k}` {res['base'][k]:g}" for k in PANELS) + ". The pre "
             f"column of the table below is therefore a distribution of ratios about "
             f"1.0, and it is the only window in which 1.0 means anything.\n")
    L.extend(md_table(["panel", "column", "pre median (raw, the base)", "n drawn",
                       "n empty", "pre min", "pre median", "pre max", "repair median",
                       "post median", "post/pre (pooled)", "pre-median spread over the "
                       "four geometries"],
                      channel_rows(res)))
    L.append("")
    L.append(f"All six channels are strictly positive on every drawn row (a squared "
             f"coherence, a count, a ratio, a power), so every panel is logarithmic and "
             f"no ordinate of this plate can reach 0; the drawn extent of the "
             f"{n_drawn} values is pinned as {plate['band'][0]:g} .. {plate['band'][1]:g} "
             f"of their own base and the ordinate ticks are "
             + ", ".join(f"{t:g}" for t in plate["log_ticks"]) + ".\n")
    L.append(f"The **pooled post/pre ratio** in the last-but-one column is the one "
             f"number this plate refuses to read as a level. On `coherence_gamma2` it is "
             f"{ratio['coherence_gamma2']:.3f}, and the four geometries disagree about "
             f"it in **sign**: "
             + ", ".join(f"{res['geometries'][t]['coherence_gamma2']['post_over_pre']:.2f} "
                         f"({t})" for t in TRACK_ORDER)
             + f". Their own pre medians already span "
             f"{gspread('coherence_gamma2'):.1f}x "
             f"inside the pre window alone, so the pooled number is a **mixture of "
             f"geometries**, not a level of the bridge: a ratio of a median of four "
             f"populations is not the median of their four ratios, and with sign "
             f"disagreement it is not even the direction. `phase_scatterer_count` and "
             f"`brightness_ratio` are the two channels whose four geometries agree "
             f"(their pre medians span "
             f"{gspread('brightness_ratio'):.3f}x and "
             f"{gspread('phase_scatterer_count'):.3f}x "
             f"respectively), and even there the direction is per geometry: see the "
             f"table in the next section.\n")
    L.append(f"No trend line, no fit and no test is drawn: on a sawtooth that changes "
             f"its level by {gspread('coherence_gamma2'):.1f}x "
             f"with the pass geometry, any line through the dots would measure the "
             f"orbits. The honest per-geometry numbers are in the JSON "
             f"(`geometries`), and the reference-chain reading is its own section "
             f"below.\n")
    L.append(f"One caveat on the panel pair. `brightness_ratio` is `intensity` divided by "
             f"a **per-geometry constant** - "
             + ", ".join(f"{bconst[t]:.3f} ({t})" for t in TRACK_ORDER)
             + f" - so panels (d) and (e) are one channel rescaled and not two "
             f"measurements (the constant holds to 1e-9 relative on all 241 rows). Both "
             f"are drawn because the export carries both columns and this figure does not "
             f"choose for the reader; the JSON says so, and so does the self-check.\n")
    L.append("## The four pass geometries\n")
    L.extend(md_table(["geometry", "orbit", "sub-swath", "platforms",
                       "n pre", "n repair", "n post", "first", "last",
                       "own pre median of gamma2 (ratio to the pooled base)",
                       "own post/pre of gamma2"],
                      geometry_rows(res)))
    L.append("")
    L.append("The plate draws, inside the pre window only, **each geometry's own pre "
             "median** of that channel as a solid rule in that geometry's own colour, "
             "and the post window's pooled median as one dashed rule; the two are not "
             "the same quantity and the dashed rule is the pooled one, because a "
             "per-geometry post median over 5-6 acquisitions per geometry would be a "
             "number this figure cannot defend. A dot is read against its own "
             "geometry's rule, never across the sawtooth.\n")
    L.extend(md_table(["channel"] + list(TRACK_ORDER),
                      geometry_channel_rows(res)))
    L.append("")
    L.append("That table is the honest reading of the plate, and it is why the pooled "
             "ratio in the previous section is quarantined: the same four post weeks "
             "move `coherence_gamma2` **up** on "
             + ", ".join(t for t in TRACK_ORDER
                         if res["geometries"][t]["coherence_gamma2"]["post_over_pre"] > 1)
             + " and **down** on "
             + ", ".join(t for t in TRACK_ORDER
                         if res["geometries"][t]["coherence_gamma2"]["post_over_pre"] <= 1)
             + ", and the same is true of `phase_coherence`, `intensity` and "
             "`coherence`. Only `phase_scatterer_count` agrees in sign across all four "
             "geometries, and there the drop is not uniform either ("
             + ", ".join(f"{res['geometries'][t]['phase_scatterer_count']['post_over_pre']:.2f} "
                         f"({t})" for t in TRACK_ORDER) + ").\n")
    L.append("## The coverage change inside the band\n")
    L.append(f"The axis is not uniform, and this is drawn as a **dotted vertical and "
             f"never as a level**: S1C joins the constellation on {cov['first_s1c']} "
             f"({cov['first_s1c_month_count']} acquisitions that month, then "
             f"{cov['per_month']['2025-04']['n']} and "
             f"{cov['per_month']['2025-05']['n']} in the two months after), so "
             f"{cov['n_s1c']} of the 241 drawn acquisitions are S1C and the axis' "
             f"monthly count goes from 9-11 to 18-20. The mechanism is the revisit: "
             + ", ".join(f"{t} {cov['per_track'][t]['steps_before']} -> "
                         f"{cov['per_track'][t]['steps_after']}" for t in TRACK_ORDER)
             + f" days. Every geometry gains a 6-day step; three of the four drop the "
             f"12-day one, `asc_44_iw3` keeps it alongside and `desc_95_iw3` had one "
             f"24-day miss before. Two acquisitions more per month is not two more "
             f"observations of the bridge - it is the same bridge sampled by two "
             f"satellites, and the plate marks the change instead of smoothing it.\n")
    L.extend(md_table(["month", "acquisitions", "platforms", "geometries", "first", "last"],
                      coverage_rows(res)))
    L.append("")
    L.append("## The steps of the chain\n")
    L.extend(md_table(["step", "how many", "what it is", "the steps"], step_rows(res)))
    L.append("")
    L.append(f"No line crosses the plate and no line crosses a hole: the mean step is "
             f"{a['span_days'] / (a['n'] - 1):.2f} days over four interleaved "
             f"geometries, and the one step longer than the {CADENCE_DAYS}-day cadence "
             f"(2024-08-28 -> 2024-09-04, 7 days) is drawn as a grey band and never "
             f"interpolated over. The table names the first {STEP_LIST_MAX} pairs of "
             f"each step length; all {len(res['steps'])} steps of the chain are in the "
             f"JSON (`steps`).\n")
    L.append("## The reference chains of the two IW1 geometries\n")
    L.append(f"Two of the four geometries carry a target window **and** its own reference "
             f"window on the very same dates - {ref['n_dates']} of the axis' 241 "
             f"acquisitions ({ref['per_track']}; pre {ref['n_dates_pre']}, repair "
             f"{ref['n_dates_repair']}, post {ref['n_dates_post']}). That is the one "
             f"comparison of this record that does **not** mix pass geometries, because "
             f"both windows are measured on the same pass at different places. It is not "
             f"drawn - it would need a second ordinate - and it is not a test; it is "
             f"reported:\n")
    L.extend(md_table(["stretch", "n dates", "target gamma2 (median)", "reference gamma2",
                       "outside-reference gamma2", "target / reference",
                       "target / outside-reference", "target / reference (scatterers)"],
                      reference_rows(res)))
    L.append("")
    post_ref = ref["per_stretch"]["post"]
    pre_ref = ref["per_stretch"]["pre"]
    L.append(f"In the post window that reading is "
             f"**{post_ref['ratio_reference_gamma2']:.3f}** - the target window's gamma2 "
             f"over its own same-pass reference window, i.e. *below* the reference - "
             f"while the pooled post/pre ratio of the same column is "
             f"{ratio['coherence_gamma2']:.3f}. In the pre window the same reading is "
             f"{pre_ref['ratio_reference_gamma2']:.3f}. Both numbers are in the JSON and "
             f"neither is drawn: the figure's plate carries the six channels, and a "
             f"reader who wants the reference reading has it here rather than "
             f"reconstructed from a plot that does not show it.\n")
    L.append("## The holes and the artefacts of the record\n")
    L.append(f"An empty field in the export is a **missing value** and this figure never "
             f"turns one into a 0. Inside the axis:\n")
    L.extend(md_table(["channel", "empty rows of 241", "which dates"], hole_rows(res)))
    L.append("")
    L.append(f"The 8 empty `coherence` rows are the first four acquisitions of the repair "
             f"stratum and the first four of the post stratum - the two places where the "
             f"chain changes state - and they are drawn as ticks on the panel floor, "
             f"never as zeros. Across the whole record ({rec['rows']} rows) two columns "
             f"are empty on **every** row: "
             + ", ".join(f"`{k}`" for k in res["holes"]["empty_on_all_record_rows"])
             + f" - they are not drawn at all. `coherence` is empty on "
             f"{res['holes']['coherence_empty_in_record']} of the record's "
             f"{rec['rows']} rows. `coherence_masked_pixels` is "
             f"{res['holes']['masked_pixels']} = "
             f"{100 * res['holes']['masked_pixels_frac_of_window']:.0f} % of the "
             f"{res['holes']['window_px']}-px window in every row of the record: a "
             f"definition (the brightest decile), not a measurement, and never drawn. "
             f"The record's window note: {masks.WINDOW_NOTE}.\n")
    L.append(f"The four co-variate controls of the record are carried on every drawn "
             f"row and none of them is plotted (a second ordinate the plate does not "
             f"have):\n")
    L.extend(md_table(["control", "median pre", "median repair", "median post",
                       "axis min", "axis max"], control_rows(res)))
    L.append("")
    L.append(f"And the record carries {len(res['artefacts'])} rows whose `shift_y` is "
             f"about -174 m where every other row of the record stays within 2.21 m - "
             f"all of them `asc_44_iw3`, all S1C, "
             f"{sum(1 for a0 in res['artefacts'] if a0['date'] <= AXIS_LAST)} of them "
             f"inside this axis. They are pinned in the JSON, never drawn and never "
             f"removed:\n")
    L.extend(md_table(["date", "geometry", "platform", "shift_x [m]", "shift_y [m]"],
                      artefact_rows(res)))
    L.append("## What this figure asserts\n")
    L.append(f"{sum(1 for c in checks if c['ok'])} of {len(checks)} self-checks pass, "
             f"each pinned against the layer it belongs to, plus "
             f"{pins_summary['n_checks']} numeric pins "
             f"({pins_summary['n_exact']} exact, "
             f"{pins_summary['n_within_tolerance']} within tolerance, "
             f"{pins_summary['n_failed']} failed). The list is in the JSON under "
             f"`self_check`; a failed check aborts the run before anything is "
             f"written.\n")
    for c in checks:
        L.append(f"- **{'ok' if c['ok'] else 'FAIL'}** - {c['check']}: {c['detail']}")
    L.append("")
    L.append("## Provenance\n")
    L.append(f"- the committed export: `data/werder_havel/werder_havel_measurements_full.txt`, "
             f"sha256 `{rec['measurements_sha256']}`, {rec['rows']} rows x "
             f"{rec['columns']} columns, {rec['n_dates']} acquisition dates "
             f"({rec['first_date']} .. {rec['last_date']}), roles "
             + ", ".join(f"{k} {v}" for k, v in rec["roles"].items())
             + "; platforms " + ", ".join(rec["platforms"])
             + "; state windows "
             + ", ".join(f"{k} {v}" for k, v in rec["state_windows"].items()) + ".")
    L.append(f"- the committed segment list: "
             f"`data/werder_havel/werder_havel_segments.txt`, sha256 "
             f"`{rec['segments_sha256']}`.")
    L.append(f"- the committed AOI reference: "
             f"`data/werder_havel/reference/werder_havel_aoi.json`, re-derived field by "
             f"field and pinned; the documented change is "
             f"`{masks.EVENT_START}` .. `{masks.EVENT_END}` "
             f"({masks.EVENT_SOURCE}).")
    L.append(f"- the record's own window note: {masks.WINDOW_NOTE}")
    L.append(f"- the plate: `figures/werder_havel/{NAME}.png` (183 mm, 600 dpi, "
             f"{PANEL_GRID[0]} x {PANEL_GRID[1]} panels); the numbers: "
             f"`data/werder_havel/{NAME}.json`.")
    L.append(f"- this run: {res['runtime_s']:.1f} s, stdlib + numpy + matplotlib only, "
             f"offline, no cache read and none written.\n")
    return "\n".join(L)
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--figdir", default=FIGDIR, help="directory for the PNG")
    ap.add_argument("--json", default=OUT_JSON, help="path of the result JSON")
    ap.add_argument("--md", default=OUT_MD, help="path of the report")
    ap.add_argument("--no-figures", action="store_true",
                    help="write the JSON and the .md only (no matplotlib output)")
    ap.add_argument("--check", action="store_true",
                    help="run the self-checks and the pins against the committed "
                         "layers and write nothing")
    args = ap.parse_args(argv)

    t0 = time.time()
    res = compute()
    a = res["axis_facts"]
    print(f"{NAME}: {a['n']} acquisitions on {a['n_days']} days, {a['first']} -> "
          f"{a['last']} ({a['span_days']} d); stretches "
          + ", ".join(f"{s} {res['stretches'][s]['n']}" for s in STRETCHES)
          + f"; six channels, "
          f"{sum(res['channels'][k]['n_drawn'] for k in PANELS)} drawn values")

    checks, P = self_check(res)
    summary = P.summary()
    print(f"  pins: {summary['n_checks']} comparisons ({summary['n_exact']} exact, "
          f"{summary['n_within_tolerance']} within {summary['float_tolerance']:g}), "
          f"{summary['n_failed']} failed")
    for c in checks:
        print(f"  [{'ok' if c['ok'] else 'FAIL'}] {c['check']}: {c['detail']}")
    for f in summary["failures"][:20]:
        print("  PIN FAILED:", f["path"], "recomputed =", f["recomputed"],
              "committed =", f["committed"])
    if summary["n_failed"]:
        _fail(f"{summary['n_failed']} of {summary['n_checks']} pin comparisons disagree "
              f"with the committed layers")
    bad = [c for c in checks if not c["ok"]]
    if bad:
        _fail(f"{len(bad)} of {len(checks)} self-checks failed - nothing written")
    print(f"  self-checks: {len(checks)} passed")

    if args.check:
        print(f"CHECK OK  ({time.time() - t0:.1f} s, nothing written)")
        return 0

    if not args.no_figures:
        fig_series(res, args.figdir)

    for path in (args.json, args.md):
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(args.json, "w", encoding="utf-8") as fh:
        json.dump(json_safe(build_json(res, checks, summary)), fh, indent=1,
                  sort_keys=True, default=json_default)
        fh.write("\n")
    with open(args.md, "w", encoding="utf-8") as fh:
        fh.write(report(res, checks, summary))
    print(f"  {os.path.relpath(args.json, ROOT)}")
    print(f"  {os.path.relpath(args.md, ROOT)}")
    print(f"DONE  ({time.time() - t0:.1f} s)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
