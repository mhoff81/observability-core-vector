#!/usr/bin/env python3
"""KDLO OCV-Paper — the acquisition series, from the first scene to the collapse.

Question: what do the channels of the echo-mask vector do over the **whole**
catalogue-covered life of the site — from the first ascending IW2 scene over the
KDLO tower (``2015-03-20``) to the collapse (``2022-12-11``)? The A-E package's
record is 60 acquisitions on one descending track (from ``2022-06-02``); this
companion puts that record back into the seven years that precede it and draws
both on **one dated axis**.

What that axis can and cannot carry
-----------------------------------
The six dimensions are ``[gamma2, P, D, A, F, S]``, and only **two** of them
exist on both sides of the geometry break:

  * ``gamma2`` — the era's backend rows carry ``coherence_gamma2`` (162/162) and
    the record's carry ``coherence_gamma2``: the same column by name;
  * the mask size ``A`` — the era carries ``coherence_masked_pixels``, the
    pipeline's own masked-pixel count, which is the row the record's ``A`` is
    recomputed from and verified against (``kdlo_ocv_channels_csv.py --verify``);
  * the controls ``wind_speed_ms`` / ``temperature_c`` — also on both sides.

``P``, ``D``, ``F`` and ``S`` are **recomputed from the committed full-dwell
strips**, and ``data/kdlo/strips/`` holds only the 60 acquisitions of the record:
the era has no mask geometry here, and its backend rows carry no phase or
frequency layer either. Those four dimensions are therefore drawn **only where
they exist** — the 16 pre-collapse acquisitions — and the era is hatched in their
panels as *not measurable*: never interpolated, never extrapolated, never pooled.

The axis has four stretches, and the figure draws all four instead of blurring
them into one series (every date below is read from the committed availability
reference, not asserted from this prose):

  era             2015-03-20 .. 2021-12-19   162 days, ASCENDING IW2, relative
                  orbit 136, catalogue hour 00 — ``gamma2`` and the mask size only
  blind gap       2021-12-20 .. 2022-02-25    69 days, no scene of either track
  catalogue only  2022-02-26 .. 2022-06-01   the descending stretch exists in the
                  CDSE catalogue but has no committed row in this record
  record (pre)    2022-06-02 .. 2022-12-11    16 days, DESCENDING IW1, relative
                  orbit 12 — all six dimensions

Normalisation: **each stretch is drawn in units of its own median** (``value /
stretch median``, with 1.0 drawn), log where the dimension spans decades. The two
stretches sit on different tracks and swaths and are imaged at different times of
day, so their levels are not comparable by construction; dividing each by its own
median shows what each stretch *does* without asserting that the era "was" a level
the record "changed from". No statistic of this figure is computed across the
break — that is the point of the empty era region in the four record-only panels.

The era's own clock is re-derived here rather than asserted: the reference's
per-year table labels the era's hour block ``00``, while the export times its rows
on the pipeline's 05:00Z ``morning`` slot, and the one metadata-anomaly day
(``2020-08-13``) sits on the ``afternoon`` slot at 16:00Z because its catalogue
burst name carries the wrong hour. Both are self-checked and neither enters a
drawn value.

The era export is committed **evidence**, never rows of the record: no number of
the A-E pin contract depends on it and it is not pooled into the channel table.
This figure is its only reader, and it re-verifies the sha256 that the
availability reference pins for it before using a single value (aborting on any
drift), and hard-checks that the era really is ``P``/``D``/``F``/``S``-free. The
three committed inputs are:

  data/kdlo/fig_kdlo_0_raw.json                  the era's decoded backend rows
  data/kdlo/kdlo_ocv_channels.csv                the record's 60 acquisitions
  data/kdlo/reference/kdlo_s1_availability.json  the era's window/days/geometry
                                                 and the export's pinned sha256

Figure (3 x 2, 183 mm wide), one dated axis 2015-03-20 .. 2022-12-11:
  a  ``gamma2``       log y — era dots (blue) vs. record dots + line (orange)
  b  ``A`` (mask size) log y — same treatment
  c  ``P``            record only — the era region is hatched "not measurable"
  d  ``D``            record only
  e  ``F``            record only
  f  ``S``            record only
Every panel carries the blind gap, the catalogue-only stretch, the collapse cut
(2022-12-11) and, hollow at the right edge, the first collapsed acquisition
(2022-12-23), which lies outside this axis (labelled in panel a).

Usage:
  python3 code/kdlo/fig_kdlo_ocv_series.py
  python3 code/kdlo/fig_kdlo_ocv_series.py --figdir /tmp/x --json /tmp/y.json --md /tmp/z.md
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import sys
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as mdates  # noqa: E402
import matplotlib.pyplot as plt   # noqa: E402
import numpy as np                # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import kdlo_ocv_core as core   # noqa: E402
import kdlo_ocv_stats as st    # noqa: E402

FIGDIR = core.FIGDIR
# The three committed inputs: the era's decoded backend rows (evidence, pinned by
# the availability reference), the record's channel table, and that reference.
ERA_EXPORT = os.path.join(core.DATA, "fig_kdlo_0_raw.json")
CHANNEL_CSV = core.CSV_PATH
AVAILABILITY_REF = os.path.join(core.DATA, "reference",
                                core.REFERENCE_FILES["s1_availability"])
OUT_JSON = os.path.join(core.DATA, "fig_kdlo_ocv_series.json")
OUT_MD = os.path.join(FIGDIR, "fig_kdlo_ocv_series.md")

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
    "hatch.linewidth": 0.4,
    "figure.dpi": 110,
    "savefig.dpi": 600,
    "savefig.bbox": "tight",
    "pdf.fonttype": 42,
})
MM = 1.0 / 25.4
COL_WIDTH = 183 * MM

# The two epoch labels of the record, in chronological order. A KDLO "state" is a
# time window, not a damage label.
PRE = st.STATE_ORDER[0]                          # "pre" — standing mast
# The axis of this figure: the first scene the catalogue has over the site, to the
# last intact acquisition of the record.
AXIS_START = core.ERA_FIRST_DAY                  # 2015-03-20
COLLAPSE_DAY = core.STATE_WINDOWS[PRE][1]        # 2022-12-11, last intact day
AXIS_END = COLLAPSE_DAY
FIRST_COLLAPSED_DAY = core.SINGLETON_DAY         # 2022-12-23, outside the axis
# The stretches of the axis, from the reference-pinned geometry break.
ERA_END = core.GEOMETRY_BREAK["last_ascending_day_before"]           # 2021-12-19
DESCENT_START = core.GEOMETRY_BREAK["first_descending_day_of_stretch"]  # 2022-02-26
BLIND_GAP_DAYS = core.GEOMETRY_BREAK["blind_gap_days"]               # 69
RECORD_START = core.STATE_WINDOWS[PRE][0]        # 2022-06-02
# The two dimensions whose era column is the same column by name: the mask size
# the era's backend reports is the ``coherence_masked_pixels`` the record's ``A``
# is recomputed from.
ERA_CHANNEL = {"gamma2": "coherence_gamma2", "A": "coherence_masked_pixels"}
LONG_DIMS = ["gamma2", "A"]
RECORD_ONLY_DIMS = [k for k in core.FEATURES if k not in LONG_DIMS]   # P, D, F, S
PANELS = LONG_DIMS + RECORD_ONLY_DIMS
# gamma2 and A span decades on both sides; P is a fraction; D, F and S span one to
# two decades inside the record only.
PANEL_SCALE = {"gamma2": "log", "A": "log", "P": "linear",
               "D": "log", "F": "log", "S": "log"}
# The era export's ten non-null columns, verbatim. Asserting this set is what makes
# "the era has no P/D/F/S, no phase and no displacement" a measurement.
ERA_NONNULL_COLUMNS = ["acquisition_ts", "coherence_gamma2",
                       "coherence_masked_pixels", "id", "orbit_direction",
                       "pass_label", "status", "temperature_c",
                       "traffic_load_label", "wind_speed_ms"]
COLOUR = {"era": "#4C72B0", "era_med": "#1F4E79", "record": "#DD8452",
          "record_med": "#B2601F", "cut": "#888888", "singleton": "#C44E52",
          "band": "#B0B0B0"}
ERA_TINT = 0.10


def _fail(msg):
    """Abort without writing: the artifact must never outrun its evidence."""
    raise SystemExit(f"REFUSING TO WRITE — {msg}")


def sha(path):
    """The sha256 of a committed input, as the artifacts record it."""
    return core.sha256(path)


def D(day):
    """``datetime.date`` of an ISO day string."""
    return dt.date(int(day[:4]), int(day[5:7]), int(day[8:10]))


def shift_day(day, days):
    """ISO day ``days`` after ``day`` (a day, not a calendar interval)."""
    return (D(day) + dt.timedelta(days=days)).isoformat()


def span_days(start, end):
    """Inclusive day count of ``[start, end]``."""
    return (D(end) - D(start)).days + 1


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


def date_axis(ax):
    """Year ticks, quarter minors — the whole axis is seven and a half years."""
    ax.xaxis.set_major_locator(mdates.YearLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    ax.xaxis.set_minor_locator(mdates.MonthLocator(bymonth=(1, 4, 7, 10)))
    ax.tick_params(axis="x", which="minor", length=2)


# ---------------------------------------------------------------------------
# The evidence, re-verified before a value of it is used
# ---------------------------------------------------------------------------
def load_availability():
    """The committed availability reference — the authority of the axis."""
    with open(AVAILABILITY_REF) as fh:
        return json.load(fh)


def geometry_break(av):
    """The break that separates the era from the record, from the reference.

    It is re-read from the committed reference and checked key by key against
    ``core.GEOMETRY_BREAK``, so no figure can draw a break the packaged constants
    disagree with.
    """
    gb = (av.get("epochs") or {}).get("geometry_break") or {}
    for k, v in core.GEOMETRY_BREAK.items():
        if gb.get(k) != v:
            _fail(f"the availability reference's geometry break disagrees with "
                  f"core.GEOMETRY_BREAK on {k}: {gb.get(k)!r} != {v!r}")
    if gb.get("ascending_days_inside_gap"):
        _fail("the reference reports ascending days inside the blind gap: "
              f"{gb['ascending_days_inside_gap']!r} — the gap is not blind")
    return gb


def era_rows(av):
    """The era's decoded backend rows, verified against the pinned reference.

    The checks run in the order in which they can fail: the file is the one the
    reference pins (by sha256), it holds the measured days of ``figure_0``, its
    days are the era's own window, and its **non-null columns are exactly the ten
    of a backend row without a phase or frequency layer** — which is what makes
    "the era carries ``gamma2`` and the mask size and no P/D/F/S" a measurement
    rather than a claim. The metadata-anomaly day is checked on its slot.
    """
    pin = (av.get("pipeline") or {}).get("figure_0_raw_export") or {}
    if pin.get("status") != "ok":
        _fail("the availability reference does not report the era export as 'ok' "
              f"(status={pin.get('status')!r})")
    got = sha(ERA_EXPORT)
    if got != pin.get("sha256"):
        _fail("data/kdlo/fig_kdlo_0_raw.json is not the file the availability "
              f"reference pins: {got} != {pin.get('sha256')}")
    with open(ERA_EXPORT) as fh:
        blob = json.load(fh)
    rows = blob.get("rows") or []
    if len(rows) != core.ERA_DAYS_EXPECTED:
        _fail(f"the era export holds {len(rows)} rows, expected "
              f"{core.ERA_DAYS_EXPECTED}")
    days = sorted({str(r.get("acquisition_ts"))[:10] for r in rows})
    if len(days) != len(rows):
        _fail("the era export carries more than one row for a day")
    if days != sorted(av["figure_0"]["acquisition_days"]):
        _fail("the era export's days differ from figure_0.acquisition_days")
    if days[0] != core.ERA_FIRST_DAY or days[-1] != core.ERA_WINDOW[1]:
        _fail(f"the era runs {days[0]}..{days[-1]}, expected "
              f"{core.ERA_FIRST_DAY}..{core.ERA_WINDOW[1]}")
    nonnull = sorted({k for r in rows for k, v in r.items() if v is not None})
    if nonnull != ERA_NONNULL_COLUMNS:
        _fail("the era export's non-null columns are " + ", ".join(nonnull)
              + " — expected " + ", ".join(ERA_NONNULL_COLUMNS))
    dirs = {str(r.get("orbit_direction")) for r in rows}
    if dirs != {"ASCENDING"}:
        _fail(f"the era is not purely ascending: {sorted(dirs)}")
    anom = set(core.ERA_METADATA_ANOMALY_DAYS)
    for d in sorted(anom):
        slots = {r["pass_label"] for r in rows
                 if str(r["acquisition_ts"])[:10] == d}
        if slots != {"afternoon"}:
            _fail(f"the era's {d} is on the {sorted(slots)} slot, expected the "
                  "alternate 'afternoon' slot")
    off = {r["pass_label"] for r in rows
           if str(r["acquisition_ts"])[:10] not in anom}
    if off != {"morning"}:
        _fail(f"the era carries non-nominal pass slots: {sorted(off)}")
    out = []
    for r in rows:
        day = str(r["acquisition_ts"])[:10]
        g, a = r.get(ERA_CHANNEL["gamma2"]), r.get(ERA_CHANNEL["A"])
        if g is None or a is None:
            _fail(f"the era row {day} has no {ERA_CHANNEL['gamma2']} or no "
                  "mask size")
        out.append({"day": day, "stretch": "era", "gamma2": float(g),
                    "A": float(a), "wind_speed_ms": r.get("wind_speed_ms"),
                    "temperature_c": r.get("temperature_c"),
                    "pass_label": r.get("pass_label"),
                    "hour": str(r["acquisition_ts"])[11:13]})
    out.sort(key=lambda r: r["day"])
    return out


def record_rows(csv_rows):
    """The record's acquisitions inside the axis: the 16 pre-collapse days.

    The four record-only dimensions must be complete on every one of them: the
    panels (c)-(f) are empty on the era side by construction and may not be patchy
    on the record side either.
    """
    rows = []
    for r in csv_rows:
        day = core.day_of(r)
        if day is None or day > AXIS_END:
            continue
        if core.state_of(day) != PRE:
            _fail(f"{day} is inside the axis but outside the '{PRE}' window")
        if r.get("orbit_direction") != "DESCENDING":
            _fail(f"the record row {day} is not descending")
        for k in core.FEATURES:
            if r.get(k) is None:
                _fail(f"the record row {day} has no {k}")
        rows.append(r)
    rows.sort(key=lambda r: core.day_of(r))
    if not rows:
        _fail("no record acquisition is inside the axis")
    if core.day_of(rows[0]) != RECORD_START or core.day_of(rows[-1]) != AXIS_END:
        _fail(f"the record stretch inside the axis runs {core.day_of(rows[0])}.."
              f"{core.day_of(rows[-1])}, expected {RECORD_START}..{AXIS_END}")
    return rows


def singleton_row(csv_rows):
    """The one acquisition between the epochs (2022-12-23), outside the axis.

    It is drawn hollow and labelled: it is the collapse's own first acquisition,
    and leaving it out would make the cut look like the end of the record.
    """
    rows = [r for r in csv_rows if core.day_of(r) == FIRST_COLLAPSED_DAY]
    if len(rows) != 1:
        _fail(f"expected exactly one acquisition on {FIRST_COLLAPSED_DAY}, "
              f"found {len(rows)}")
    for k in core.FEATURES:
        if rows[0].get(k) is None:
            _fail(f"the {FIRST_COLLAPSED_DAY} acquisition has no {k}")
    return rows[0]


def segments(av, era, record):
    """The four stretches of the axis, with the dates they actually cover.

    The two geometries are read off the availability reference rather than typed
    in: the era's orbit, swath and catalogue hour block come from its own part of
    the reference (the hour block from the reference's per-year table, over the
    years the drawn era rows actually cover), and the record's orbit and swath
    from the record's part. The two stretches sit on different tracks, so their
    levels are not comparable by construction — one of the reasons no statistic of
    this figure is ever taken over both.
    """
    gb = geometry_break(av)
    years = sorted({r["day"][:4] for r in era})
    hours = sorted({h for y, blk in av["epochs"]["per_year"].items()
                    if y in years for h in blk["utc_hours"]})
    era_geom = (f"{list(av['figure_0']['orbits'])[0]} · "
                f"{list(av['figure_0']['swaths'])[0]} · catalogue hour {hours[0]}")
    rec_geom = (f"{list(av['figure_1']['orbits'])[0]} · "
                f"{list(av['figure_1']['swaths'])[0]}")
    return [
        {"key": "era", "label": "the ascending IW2 era",
         "start": AXIS_START, "end": ERA_END,
         "span_days": span_days(AXIS_START, ERA_END), "filled_days": len(era),
         "geometry": era_geom,
         "dims": list(LONG_DIMS),
         "note": f"{len(era)} measured days, gamma2 and the mask size only (the "
                 "four mask-geometry dimensions are not recomputable: no "
                 "full-dwell strips in data/kdlo/strips/ for this geometry)"},
        {"key": "blind_gap", "label": "the blind gap",
         "start": shift_day(ERA_END, 1), "end": shift_day(DESCENT_START, -1),
         "span_days": span_days(shift_day(ERA_END, 1), shift_day(DESCENT_START, -1)),
         "filled_days": 0, "catalogue_gap_days": gb.get("blind_gap_days"),
         "geometry": "—", "dims": [],
         "note": "no scene of either track; the reference's "
                 f"{gb.get('blind_gap_days')}-day blind gap counts through "
                 f"{DESCENT_START}, the first descending day, so the empty stretch "
                 "drawn here is the day before it"},
        {"key": "catalogue_only", "label": "in the catalogue, not in the record",
         "start": DESCENT_START, "end": shift_day(RECORD_START, -1),
         "span_days": span_days(DESCENT_START, shift_day(RECORD_START, -1)),
         "filled_days": 0, "geometry": f"{rec_geom} (catalogue)",
         "dims": [],
         "note": "the first descending stretch of the track that carries the "
                 "collapse; the record's committed rows start 2022-06-02"},
        {"key": "record", "label": "the record's pre-collapse stretch",
         "start": RECORD_START, "end": AXIS_END, "filled_days": len(record),
         "span_days": span_days(RECORD_START, AXIS_END),
         "geometry": f"{rec_geom} · one per day",
         "dims": list(core.FEATURES),
         "note": f"all six dimensions, one acquisition per day, cut at the "
                 f"collapse {COLLAPSE_DAY}"},
    ]


# ---------------------------------------------------------------------------
# The series
# ---------------------------------------------------------------------------
def _rel(path):
    """A path as the artifacts name it: relative to the repository root."""
    root = os.path.abspath(os.path.join(HERE, os.pardir, os.pardir))
    return os.path.relpath(path, root)


def _median(rows, dim):
    vals = [r[dim] for r in rows if r.get(dim) is not None]
    return st._stats(vals) if vals else None


def compute():
    """The era rows, the record rows, the per-stretch medians, the segments.

    The medians are the 1.0 of the figure: **each stretch is divided by its own
    median** (the era's 162-day median for the era, the 16-day median of
    2022-06-02..2022-12-11 for the record), so no level is compared across the
    geometry break. The record's first month (2022-06) is kept separately: it is
    the baseline the months companion uses, and the cross-read of the two
    artifacts hangs on it.
    """
    t0 = time.time()
    av = load_availability()
    era = era_rows(av)
    csv_rows = core.load_csv()
    record = record_rows(csv_rows)
    single = singleton_row(csv_rows)
    later = sorted(core.day_of(r) for r in csv_rows if core.day_of(r) > AXIS_END)

    groups = {"era": era,
              "record": record,
              "record_first_month": [r for r in record
                                     if core.day_of(r)[:7] == RECORD_START[:7]]}
    dims_of = {"era": LONG_DIMS, "record": list(core.FEATURES),
               "record_first_month": list(core.FEATURES)}
    med = {}
    for g, rows in groups.items():
        med[g] = {}
        for k in dims_of[g]:
            vals = [r[k] for r in rows]
            med[g][k] = st._stats(vals)["median"] if vals else None
    stats = {"era": {k: _median(era, k) for k in LONG_DIMS},
             "record": {k: _median(record, k) for k in core.FEATURES}}

    era_out = [{"day": r["day"], "stretch": "era",
                "pass_label": r["pass_label"], "hour": r["hour"],
                "wind_speed_ms": r["wind_speed_ms"],
                "temperature_c": r["temperature_c"],
                "channels": {k: r[k] for k in LONG_DIMS},
                "ratio": {k: r[k] / med["era"][k] for k in LONG_DIMS}}
               for r in era]
    record_out = [{"day": core.day_of(r), "stretch": PRE,
                   "pass_label": r.get("pass_label"),
                   "wind_speed_ms": r.get("wind_speed_ms"),
                   "temperature_c": r.get("temperature_c"),
                   "channels": {k: r[k] for k in core.FEATURES},
                   "ratio": {k: r[k] / med["record"][k] for k in core.FEATURES}}
                  for r in record]
    single_out = {"day": FIRST_COLLAPSED_DAY,
                  "wind_speed_ms": single.get("wind_speed_ms"),
                  "temperature_c": single.get("temperature_c"),
                  "channels": {k: single[k] for k in core.FEATURES},
                  "ratio": {k: single[k] / med["record"][k]
                            for k in core.FEATURES}}

    # The era's own year-by-year medians: the only aggregation over the era, and
    # it stays inside the era (a year is never compared with the record).
    per_year = {}
    for r in era:
        per_year.setdefault(r["day"][:4], []).append(r)
    per_year = {y: {"n": len(rs), "first_day": rs[0]["day"],
                    "last_day": rs[-1]["day"],
                    "channels": {k: st._stats([x[k] for x in rs])
                                 for k in LONG_DIMS}}
                for y, rs in sorted(per_year.items())}

    sources = {
        "era_export": {"path": _rel(ERA_EXPORT), "sha256": sha(ERA_EXPORT),
                       "rows": len(era),
                       "pinned_by": "availability.pipeline."
                                    "figure_0_raw_export.sha256"},
        "channel_table": {"path": _rel(CHANNEL_CSV), "sha256": sha(CHANNEL_CSV)},
        "availability_reference": {"path": _rel(AVAILABILITY_REF),
                                   "sha256": sha(AVAILABILITY_REF)},
    }
    avail = {
        "era": {k: av["figure_0"][k] for k in
                ("window", "n_products", "n_bursts", "n_acquisition_days",
                 "orbits", "swaths", "polarisations", "platforms")},
        "record": {k: av["figure_1"][k] for k in
                   ("window", "n_products", "n_acquisition_days", "orbits",
                    "swaths", "polarisations", "platforms")},
        "per_year": av["epochs"]["per_year"],
    }
    return {
        "era_rows": era_out,
        "record_rows": record_out,
        "singleton": single_out,
        "later_days": later,
        "per_year": per_year,
        "medians": med,
        "stats": stats,
        "n_era": len(era),
        "n_record": len(record),
        "segments": segments(av, era, record),
        "geometry_break": av["epochs"]["geometry_break"],
        "availability": avail,
        "sources": sources,
        "n_rows": len(era) + len(record),
        "runtime_s": time.time() - t0,
    }


# ---------------------------------------------------------------------------
# The figure
# ---------------------------------------------------------------------------
def draw_stretches(ax, hatch_era):
    """The four stretches of the axis, drawn — never interpolated across.

    ``hatch_era`` is set for the four record-only panels. There the era is not a
    series with holes but a region where the dimension **does not exist**, so it
    is hatched and labelled instead of being left blank or, worse, filled.
    """
    ax.set_xlim(D(AXIS_START), D(FIRST_COLLAPSED_DAY))
    # the era: tinted where the two shared dimensions live, hatched where the four
    # record-only dimensions cannot be measured
    if hatch_era:
        ax.axvspan(D(AXIS_START), D(ERA_END), facecolor=COLOUR["era"],
                   alpha=ERA_TINT, hatch="//////", edgecolor="#8FA8C8", lw=0.0,
                   zorder=0.5)
    else:
        ax.axvspan(D(AXIS_START), D(ERA_END), facecolor=COLOUR["era"],
                   alpha=0.055, lw=0.0, zorder=0.5)
    # the two stretches that fill the era and the record apart: no scene of either
    # track, and the descending stretch that the record does not carry
    for lo, hi, hatch, col in (
            (shift_day(ERA_END, 1), shift_day(DESCENT_START, -1), "\\\\\\",
             "#9A9A9A"),
            (DESCENT_START, shift_day(RECORD_START, -1), "...", "#B5B5B5")):
        ax.axvspan(D(lo), D(hi), facecolor="#EDEDED", hatch=hatch,
                   edgecolor=col, lw=0.0, zorder=0.4)
        mid = D(lo) + (D(hi) - D(lo)) / 2
        ax.text(mid, 0.5,
                "blind gap" if hatch.startswith("\\") else "catalogue only",
                rotation=90, ha="center", va="center", fontsize=5.0,
                color="#5F5F5F", transform=ax.get_xaxis_transform(), zorder=1,
                clip_on=True)
    # the record
    ax.axvspan(D(RECORD_START), D(AXIS_END), facecolor=COLOUR["record"],
               alpha=0.075, lw=0.0, zorder=0.6)
    # the collapse: the cut the whole package is built on
    ax.axvline(D(AXIS_END), color=COLOUR["cut"], lw=0.8, ls=(0, (2.5, 1.6)),
               zorder=1.5)


def panel_ylim(ax, values, scale):
    """The panel's own window. 1.0 is each stretch's median and is drawn."""
    vals = [v for v in values if v is not None and v > 0]
    lo, hi = min(vals + [1.0]), max(vals + [1.0])
    if scale == "log":
        ax.set_yscale("log")
        ax.set_ylim(lo / 2.2, hi * 2.2)
    else:
        pad = 0.10 * (hi - lo) if hi > lo else 0.10
        ax.set_ylim(lo - pad, hi + pad)
    ax.axhline(1.0, color=COLOUR["cut"], lw=0.8, zorder=1)
    ax.set_ylabel("value / own-stretch median")


def draw_long(ax, res, dim):
    """A dimension that exists on **both** stretches: era dots and record line.

    Both stretches are divided by their **own** median, so the two colours are two
    levels, not one series. The era's year medians are drawn as bars over each
    year's own day range: that is an aggregation *inside* the era (and the only one
    here), never a bridge towards the record.
    """
    era, rec = res["era_rows"], res["record_rows"]
    ex = [D(r["day"]) for r in era]
    ey = [r["ratio"][dim] for r in era]
    rx = [D(r["day"]) for r in rec]
    ry = [r["ratio"][dim] for r in rec]
    for blk in res["per_year"].values():
        m = blk["channels"][dim]["median"] / res["medians"]["era"][dim]
        ax.plot([D(blk["first_day"]), D(blk["last_day"])], [m, m], "-", lw=2.2,
                color=COLOUR["era_med"], alpha=0.5, zorder=2,
                solid_capstyle="butt")
    ax.plot(ex, ey, "o", ms=2.2, color=COLOUR["era"], alpha=0.7,
            markeredgewidth=0, zorder=3,
            label=f"era day ({res['n_era']})")
    ax.plot(rx, ry, "-", lw=0.9, color=COLOUR["record"], alpha=0.9, zorder=3,
            label=f"record day ({res['n_record']})")
    ax.plot(rx, ry, "o", ms=3.4, color=COLOUR["record"], zorder=4)
    ax.plot([D(FIRST_COLLAPSED_DAY)], [res["singleton"]["ratio"][dim]], "o",
            ms=4.4, mfc="none", mec=COLOUR["singleton"], mew=1.0, zorder=5)
    panel_ylim(ax, ey + ry + [res["singleton"]["ratio"][dim]],
               PANEL_SCALE[dim])


def draw_record_only(ax, res, dim):
    """One of the four dimensions the era cannot carry: the record, and a hatch.

    No era value is drawn, none is derived: the hatched region says *not
    measurable here*, and the reason is written into it.
    """
    rec = res["record_rows"]
    rx = [D(r["day"]) for r in rec]
    ry = [r["ratio"][dim] for r in rec]
    ax.plot(rx, ry, "-", lw=0.9, color=COLOUR["record"], alpha=0.9, zorder=3,
            label=f"record day ({res['n_record']})")
    ax.plot(rx, ry, "o", ms=3.4, color=COLOUR["record"], zorder=4)
    ax.plot([D(FIRST_COLLAPSED_DAY)], [res["singleton"]["ratio"][dim]], "o",
            ms=4.4, mfc="none", mec=COLOUR["singleton"], mew=1.0, zorder=5)
    panel_ylim(ax, ry + [res["singleton"]["ratio"][dim]], PANEL_SCALE[dim])
    ax.text(0.30, 0.5, "not measurable in the era —\nIW2/ASC has no full-dwell "
            "strips\n(P, D, F, S are never interpolated)",
            transform=ax.transAxes, ha="center", va="center", fontsize=6.2,
            color="#4F4F4F", zorder=6,
            bbox=dict(fc="white", ec="none", alpha=0.85, pad=2.2))


def fig_series(res, out_dir=None):
    """The 3 x 2 companion: one dated axis, four stretches, six dimensions."""
    fig, axes = plt.subplots(3, 2, figsize=(COL_WIDTH, COL_WIDTH * 0.64))
    flat = list(np.ravel(axes))
    for letter, dim, ax in zip("abcdef", PANELS, flat):
        draw_stretches(ax, hatch_era=(dim in RECORD_ONLY_DIMS))
        if dim in LONG_DIMS:
            draw_long(ax, res, dim)
        else:
            draw_record_only(ax, res, dim)
        date_axis(ax)
        panel(ax, letter, core.CHANNEL_TITLE[dim])
    for ax in flat[:4]:
        ax.tick_params(labelbottom=False)
    for ax in flat[4:]:
        ax.set_xlabel("acquisition day")
    a = flat[0]
    a.legend(loc="lower left", frameon=False, ncol=3, handletextpad=0.3,
             columnspacing=0.8, borderaxespad=0.2)
    a.annotate("collapse " + COLLAPSE_DAY, xy=(D(AXIS_END), a.get_ylim()[0]),
               xytext=(-3, 3), textcoords="offset points", ha="right",
               va="bottom", fontsize=6.5, color="#5F5F5F")
    a.annotate("first collapsed\n" + FIRST_COLLAPSED_DAY,
               xy=(D(FIRST_COLLAPSED_DAY),
                   res["singleton"]["ratio"]["gamma2"]),
               xytext=(-7, 14), textcoords="offset points", ha="right",
               va="bottom", fontsize=6.2, color=COLOUR["singleton"],
               arrowprops=dict(arrowstyle="-", lw=0.6,
                               color=COLOUR["singleton"]))
    flat[2].legend(loc="lower left", frameon=False, handletextpad=0.3,
                   borderaxespad=0.2)
    fig.suptitle("KDLO — the echo-mask vector from the first IW2 scene to the "
                 "collapse", x=0.004, y=1.015, ha="left", fontsize=10.0)
    fig.text(0.004, 0.995, "each stretch in units of its own median (log where a "
             "dimension spans decades): 162 ascending-IW2 days, gamma^2 and the "
             "mask size only, then the 16 descending-IW1 days of the record — "
             "pooled across the break nowhere", ha="left", va="top",
             fontsize=7.0, color="#444444")
    fig.tight_layout(pad=0.5, h_pad=1.1, w_pad=1.3)
    return save_fig(fig, "fig_kdlo_ocv_series", out_dir)


# ---------------------------------------------------------------------------
# The artifacts
# ---------------------------------------------------------------------------
def md_table(header, rows):
    out = ["| " + " | ".join(header) + " |",
           "| " + " | ".join("---" for _ in header) + " |"]
    for r in rows:
        out.append("| " + " | ".join(str(c) for c in r) + " |")
    return out


def _fmt(x, nd=4):
    """House formatting: a fixed significant-digit count, None as an en dash."""
    if x is None:
        return "–"
    if isinstance(x, float):
        return f"{x:.{nd}g}"
    return str(x)


def build_json(res):
    """The committed artifact: every row, every median, every pinned date."""
    return {
        "figure": "fig_kdlo_ocv_series",
        "script": "code/kdlo/fig_kdlo_ocv_series.py",
        "site": "KDLO media tower (Garden City SD)",
        "title": "KDLO OCV-Paper: the echo-mask vector from the first IW2 scene "
                 "to the collapse",
        "role": "companion — the A-E record (60 descending acquisitions, "
                "2022-06-02..2023-05-…, released 2025-11-24) drawn inside the "
                "seven and a half years of catalogue coverage that precede it",
        "axis": {
            "start": AXIS_START, "end": AXIS_END,
            "first_collapsed_day": FIRST_COLLAPSED_DAY,
            "note": "one dated axis; the collapse cut is 2022-12-11 and the first "
                    "collapsed acquisition (2022-12-23) is drawn hollow at the "
                    "right edge",
        },
        "dims": {"long": LONG_DIMS, "record_only": RECORD_ONLY_DIMS,
                 "order": list(core.FEATURES)},
        "era_channel": dict(ERA_CHANNEL),
        "normalisation": {
            "rule": "each stretch divided by its own median",
            "era": f"median of the {res['n_era']} era days "
                   f"({AXIS_START}..{ERA_END})",
            "record": f"median of the {res['n_record']} record days inside the "
                      f"axis ({RECORD_START}..{AXIS_END})",
            "no_pooling": "the era (ASCENDING, IW2, relative orbit 136, catalogue "
                          "hour 00) and the record (DESCENDING, IW1, relative orbit "
                          "12) are two geometries read off the availability "
                          "reference; no number of this figure is computed over both, "
                          "and the 69-day blind gap is never bridged",
        },
        "segments": res["segments"],
        "geometry_break": res["geometry_break"],
        "availability": res["availability"],
        "sources": res["sources"],
        "counts": {"era_rows": res["n_era"], "record_rows": res["n_record"],
                   "rows": res["n_rows"],
                   "record_days_after_cut": len(res["later_days"])},
        "medians": res["medians"],
        "stats": res["stats"],
        "era_rows": res["era_rows"],
        "record_rows": res["record_rows"],
        "singleton": res["singleton"],
        "later_days": res["later_days"],
        "per_year": res["per_year"],
        "checks": res["checks"],
        "n_checks": len(res["checks"]),
    }


def report(res):
    """The companion's own reading — descriptive, with the empty region named."""
    med, segs = res["medians"], res["segments"]
    md = []
    md.append("# KDLO media tower (Garden City SD) — the echo-mask vector from the "
              "first IW2 scene to the collapse")
    md.append("")
    md.append(f"**Site** KDLO media tower (Garden City SD) · **axis** {AXIS_START} "
              f"… {AXIS_END} (the collapse cut) · **rows** {res['n_rows']} = era "
              f"{res['n_era']} + record {res['n_record']}")
    md.append("")
    md.append("```")
    md.append("x = [gamma2, P, D, A, F, S]   echo-mask vector")
    for key, lab in (("era", "era"), ("blind_gap", "gap"), ("catalogue_only", "cat"),
                     ("record", "record")):
        s = next(x for x in segs if x["key"] == key)
        md.append(f"{lab:<7} {s['start']}..{s['end']}  {s['span_days']:>5} d  "
                  f"{s['filled_days']:>4} acq  {s['geometry']}")
    md.append("```")
    md.append("")
    md.append(f"The availability reference calls this blind gap **{BLIND_GAP_DAYS} "
              f"days** because it counts *through* {DESCENT_START}, the first "
              f"descending day of the track that carries the collapse; the empty "
              f"stretch drawn here is the "
              f"{span_days(shift_day(ERA_END, 1), shift_day(DESCENT_START, -1))} "
              f"days before it. The two agree on where the gap ends.")
    md.append("")
    md.append("The six dimensions exist together on only **one** of those "
              "stretches. `gamma2` and the mask size `A` exist on both — the era's "
              "backend rows carry `coherence_gamma2` and "
              "`coherence_masked_pixels`, the same columns the record carries. "
              "`P`, `D`, `F` and `S` are recomputed from the committed full-dwell "
              "strips, which exist only for the record: the era has no mask "
              "geometry, and its rows have no phase or frequency layer either.")
    md.append("")
    md.append("## The four stretches of the axis")
    md.append("")
    body = [[s["label"], s["start"], s["end"], s["span_days"], s["filled_days"],
             ", ".join(s["dims"]) or "–"] for s in segs]
    md += md_table(["stretch", "first day", "last day", "calendar days",
                    "days with data", "dimensions that exist here"], body)
    md.append("")
    md.append("## How the axis is read")
    md.append("")
    md.append(f"Each stretch is divided by **its own** median: era `gamma2` "
              f"{_fmt(med['era']['gamma2'])}, era `A` {_fmt(med['era']['A'])} over "
              f"the {res['n_era']} era days; record `gamma2` "
              f"{_fmt(med['record']['gamma2'])}, record `A` "
              f"{_fmt(med['record']['A'])} over the {res['n_record']} record days "
              f"of the axis. The 1.0 line is therefore **two different levels in "
              f"one panel**, and the two colours of panels (a) and (b) are two "
              f"series of self-ratios, not one series with a step in it. Levels "
              f"are not comparable across the break by construction (different "
              f"track, different swath, different pass); what each stretch *does* "
              f"is.")
    md.append("")
    md.append(f"The record-only ratio for the four mask-geometry dimensions uses "
              f"the same {res['n_record']}-day median of the record "
              f"({RECORD_START}..{AXIS_END}). Log axes where a dimension spans "
              f"decades, so a value below 1 stays visible.")
    md.append("")
    md.append("## The record inside the axis, day by day")
    md.append("")
    body = []
    for r in res["record_rows"]:
        body.append([r["day"], _fmt(r["wind_speed_ms"], 3),
                     _fmt(r["temperature_c"], 4)]
                    + [_fmt(r["channels"][k]) for k in core.FEATURES]
                    + [_fmt(r["ratio"][k], 3) for k in LONG_DIMS])
    md += md_table(["day", "wind m/s", "T degC"] + core.FEATURES
                   + [k + " /med" for k in LONG_DIMS], body)
    md.append("")
    md.append(f"The record's own window starts {RECORD_START}: the "
              f"{span_days(DESCENT_START, shift_day(RECORD_START, -1))} descending "
              f"days between the blind gap and it ({DESCENT_START}.."
              f"{shift_day(RECORD_START, -1)}) exist in the catalogue but carry no "
              f"committed row in this record, so they are hatched, not filled.")
    md.append("")
    md.append(f"## The era inside the axis, year by year ({len(res['per_year'])} "
              f"years)")
    md.append("")
    body = []
    for y, blk in res["per_year"].items():
        body.append([y, blk["n"], blk["first_day"] + " … " + blk["last_day"]]
                    + [_fmt(blk["channels"][k]["median"]) for k in LONG_DIMS]
                    + [_fmt(blk["channels"][k]["median"] / med["era"][k], 3)
                       for k in LONG_DIMS])
    md += md_table(["year", "days", "first … last"] + [k + " median" for k in LONG_DIMS]
                   + [k + " /era med" for k in LONG_DIMS], body)
    md.append("")
    md.append("These year medians are the only aggregation drawn over the era "
              "(bars in panels a and b), and they stay **inside** the era: a year "
              "is never compared with the record. They are era-track levels.")
    md.append("")
    md.append("## Why panels (c)-(f) are empty on the era side")
    md.append("")
    md.append("A line, a dashed interpolation or a greyed continuation over "
              "2015-03-20..2021-12-19 in panels (c)-(f) would be an invention. "
              "`P`, `D`, `F` and `S` are recomputed here from the committed "
              "full-dwell strips, and `data/kdlo/strips/` holds the 60 "
              "acquisitions of the **descending** record only; the era's backend "
              "rows carry no phase or frequency layer at all (their non-null "
              "columns are exactly the ten of the export, checked against the "
              "pinned reference before anything is drawn). The era is therefore "
              "drawn in those four panels as a hatched region read **not "
              "measurable here** — with the reason written inside it — and the "
              "figure carries no value of those four dimensions before "
              f"{RECORD_START}.")
    md.append("")
    md.append("## The collapse cut and the first collapsed acquisition")
    md.append("")
    md.append(f"The axis stops at the collapse, **{AXIS_END}** — the last intact "
              f"acquisition of the record — and the record's last "
              f"{res['n_record']} days run up to it. The next acquisition, "
              f"**{FIRST_COLLAPSED_DAY}**, belongs to *neither* epoch: it is the "
              f"first collapsed scene and lies outside this axis, so it is drawn "
              f"hollow at the right edge in every panel (labelled in panel a).")
    md.append("")
    single = res["singleton"]
    body = []
    for k in core.FEATURES:
        body.append([k, _fmt(single["channels"][k]), _fmt(single["ratio"][k], 3),
                     _fmt(res["record_rows"][-1]["channels"][k]),
                     _fmt(res["record_rows"][-1]["ratio"][k], 3)])
    md += md_table(["dimension", f"{FIRST_COLLAPSED_DAY}", f"/record median",
                    f"{AXIS_END} (last intact)", "/record median"], body)
    md.append("")
    n_below = sum(1 for k in core.FEATURES if single["ratio"][k] < 1.0)
    n_above = len(core.FEATURES) - n_below
    md.append(f"On {FIRST_COLLAPSED_DAY} {n_below} of the six dimensions sit "
              f"**below** their own record median and {n_above} above; `gamma2` is "
              f"{_fmt(single['ratio']['gamma2'], 3)}× the record median, `A` is "
              f"{_fmt(single['ratio']['A'], 3)}×. That one scene is the only "
              f"post-collapse acquisition this figure draws, and it is drawn as a "
              f"single point: n = 1 is a point, not a level.")
    md.append("")
    md.append("## What this figure does not say")
    md.append("")
    md.append("It does not say that the era \"was\" quiet and the record "
              "\"declined\": the two stretches are different geometries "
              f"({segs[0]['geometry']} vs. {segs[3]['geometry']}), the season "
              "coverage differs, and no statistic here is computed over both. It "
              "does not extend any of the six dimensions across the blind gap. It "
              "does not carry the collapse date's meaning into the era: the cut is "
              "drawn where the record's last intact acquisition is, not where a "
              "cause is.")
    md.append("")
    md.append("## Inputs (all committed; all re-read before they are drawn)")
    md.append("")
    body = []
    for key, s in res["sources"].items():
        body.append([key, "`" + s["path"] + "`", "`" + s["sha256"][:16] + "…`"])
    md += md_table(["input", "path", "sha256 (first 16)"], body)
    md.append("")
    md.append(f"The era export's sha256 is the one the availability reference "
              f"pins (`pipeline.figure_0_raw_export`), and this figure re-verifies "
              f"it before using a single era value; the geometry break "
              f"({res['geometry_break']['last_ascending_day_before']} → "
              f"{res['geometry_break']['first_descending_day_of_stretch']}, "
              f"{res['geometry_break']['blind_gap_days']} days) is re-read from the "
              f"same reference and cross-checked against the packaged constants. "
              f"The figure carries **{len(res['checks'])}** self-checks, all of "
              f"which abort before any file is written.")
    md.append("")
    md.append("```bash")
    md.append("python3 code/kdlo/kdlo_ocv_channels_csv.py --verify-only   # the table")
    md.append("python3 code/kdlo/fig_kdlo_ocv_series.py                   # this figure")
    md.append("```")
    return "\n".join(md) + "\n"


# ---------------------------------------------------------------------------
# Self-checks — recomputed from the raw inputs, never from the results above
# ---------------------------------------------------------------------------
def self_checks(res):
    """Re-derive what the figure draws from the three inputs, and abort on drift.

    Nothing here reads the pipeline's intermediate values: the counts, the day
    ranges, both medians, the per-year medians and the single collapsed
    acquisition are recomputed from the committed files and compared with what is
    about to be drawn. The last three checks tie this artifact to the A-E package
    by its own pinned numbers (the collapse pair's `gamma2` and mask size), so a
    drift in either direction is caught. Any failure exits before a file is
    written.
    """
    checks = {}

    def eq(name, got, want):
        ok = got == want
        checks[name] = {"ok": bool(ok), "got": got, "expected": want}
        if not ok:
            _fail(f"self-check '{name}': got {got!r}, expected {want!r}")

    with open(ERA_EXPORT) as fh:
        raw = json.load(fh)["rows"]
    av = load_availability()
    csv_rows = core.load_csv()

    # --- the era: what it carries, and what it cannot carry -----------------
    eq("era.rows", len(raw), core.ERA_DAYS_EXPECTED)
    days = sorted({str(r["acquisition_ts"])[:10] for r in raw})
    eq("era.distinct_days", len(days), core.ERA_DAYS_EXPECTED)
    eq("era.first_last", [days[0], days[-1]],
       [core.ERA_FIRST_DAY, core.ERA_WINDOW[1]])
    eq("era.window_inside_reference", days == sorted(av["figure_0"]["acquisition_days"]),
       True)
    eq("era.nonnull_columns",
       sorted({k for r in raw for k, v in r.items() if v is not None}),
       ERA_NONNULL_COLUMNS)
    for k in ("phase_coherence", "phase_rms_rad", "phase_detected_frequency_hz",
              "measured_frequency_hz", "baseline_frequency_hz",
              "displacement_los_m", "ndvi", "sub_aperture_modulation"):
        eq(f"era.{k}_never_measured", sorted({r.get(k) for r in raw}), [None])
    eq("era.directions", sorted({str(r["orbit_direction"]) for r in raw}),
       ["ASCENDING"])
    eq("era.orbit", sorted(av["figure_0"]["orbits"]), [core.ERA_ORBIT])
    eq("era.swath", list(av["figure_0"]["swaths"]), [core.ERA_SWATH])
    eq("era.slot_hours_nominal",
       sorted({str(r["acquisition_ts"])[11:13] for r in raw
               if r["pass_label"] == "morning"}), ["05"])
    eq("era.hour_anomaly_off_nominal",
       "05" not in {str(r["acquisition_ts"])[11:13] for r in raw
                    if r["pass_label"] == "afternoon"}, True)
    eq("era.slot_hours_drawn", sorted({r["hour"] for r in res["era_rows"]}),
       ["05", "16"])
    eq("era.reference_hour_blocks",
       sorted({h for y, blk in av["epochs"]["per_year"].items()
               if y in res["per_year"] for h in blk["utc_hours"]}), ["00"])
    eq("era.reference_years", sorted(res["per_year"]),
       sorted(str(y) for y in av["epochs"]["per_year"] if str(y) <= ERA_END[:4]))
    eq("era.pass_slots",
       sorted({r["pass_label"] for r in raw}), ["afternoon", "morning"])
    eq("era.anomaly_days",
       sorted({str(r["acquisition_ts"])[:10] for r in raw
               if r["pass_label"] == "afternoon"}),
       core.ERA_METADATA_ANOMALY_DAYS)
    eq("era.break.crosses_era", res["geometry_break"]["last_ascending_day_before"],
       ERA_END)
    eq("era.break.blind_gap_days", res["geometry_break"]["blind_gap_days"],
       BLIND_GAP_DAYS)
    # The reference counts the gap *through* the first descending day
    # (2022-02-26, a catalogue day the record does not carry); the empty stretch
    # this figure draws stops the day before it, so its inclusive count is one day
    # smaller — both agree on the day the gap ends.
    gap = next(s for s in res["segments"] if s["key"] == "blind_gap")
    eq("era.break.blind_gap_drawn_span", gap["span_days"],
       span_days(shift_day(ERA_END, 1), shift_day(DESCENT_START, -1)))
    eq("era.break.blind_gap_reference_counts", gap["span_days"] + 1,
       BLIND_GAP_DAYS)
    eq("era.break.blind_gap_catalogue_days", gap["catalogue_gap_days"],
       BLIND_GAP_DAYS)
    eq("era.break.no_ascending_inside_gap",
       res["geometry_break"]["ascending_days_inside_gap"], [])
    eq("era.per_year_days",
       {y: blk["n"] for y, blk in res["per_year"].items()},
       {y: len([r for r in raw if str(r["acquisition_ts"])[:4] == y])
        for y in res["per_year"]})
    eq("era.per_year_total",
       sum(blk["n"] for blk in res["per_year"].values()), core.ERA_DAYS_EXPECTED)
    eq("era.reference_products_total",
       sum(v["n"] for y, v in av["epochs"]["per_year"].items()
           if y in res["per_year"]), av["figure_0"]["n_products"])
    eq("era.reference_days_total", av["figure_0"]["n_acquisition_days"],
       core.ERA_DAYS_EXPECTED)
    eq("era.rows_in_axis", [r["day"] for r in res["era_rows"]], days)

    # --- the record: only complete days, and only inside the axis -----------
    inside = sorted([r for r in csv_rows if core.day_of(r) <= AXIS_END],
                    key=core.day_of)
    eq("record.days_in_axis", len(inside), res["n_record"])
    eq("record.first_last", [core.day_of(inside[0]), core.day_of(inside[-1])],
       [RECORD_START, AXIS_END])
    eq("record.directions", sorted({r["orbit_direction"] for r in inside}),
       ["DESCENDING"])
    eq("record.orbit_one", len(av["figure_1"]["orbits"]), 1)
    eq("record.orbit_is_descending",
       sorted(str(o).split("|")[0] for o in av["figure_1"]["orbits"]),
       ["DESCENDING"])
    eq("record.swath_one", len(av["figure_1"]["swaths"]), 1)
    eq("record.all_six_dims_complete",
       sorted({k for r in inside for k in core.FEATURES if r.get(k) is None}), [])
    eq("record.rows_drawn", [r["day"] for r in res["record_rows"]],
       [core.day_of(r) for r in inside])
    # Everything after the cut: first the collapse's own singleton acquisition,
    # then the rebuild window — and nothing else, and nothing inside the axis.
    rb = core.STATE_WINDOWS[st.STATE_ORDER[1]]
    eq("record.after_cut_first", res["later_days"][0], FIRST_COLLAPSED_DAY)
    eq("record.after_cut_outside_singleton",
       [d for d in res["later_days"] if not rb[0] <= d <= rb[1]],
       [FIRST_COLLAPSED_DAY])
    eq("record.after_cut_last", res["later_days"][-1], rb[1])
    eq("record.after_cut_inside_axis",
       [d for d in res["later_days"] if d <= AXIS_END], [])
    eq("record.after_cut_count", len(res["later_days"]),
       1 + len([r for r in csv_rows
                if core.state_of(core.day_of(r)) == st.STATE_ORDER[1]]))

    # --- the two medians and the per-year medians, recomputed --------------
    g_era = [float(r["coherence_gamma2"]) for r in raw]
    a_era = [float(r["coherence_masked_pixels"]) for r in raw]
    eq("median.era.gamma2", res["medians"]["era"]["gamma2"],
       float(np.median(g_era)))
    eq("median.era.A", res["medians"]["era"]["A"], float(np.median(a_era)))
    for k in core.FEATURES:
        eq(f"median.record.{k}", res["medians"]["record"][k],
           float(np.median([r[k] for r in inside])))
    first_month = [r for r in inside if core.day_of(r)[:7] == RECORD_START[:7]]
    eq("median.first_month.n", len(first_month) > 0, True)
    eq("median.first_month.A", res["medians"]["record_first_month"]["A"],
       float(np.median([r["A"] for r in first_month])))
    for y, blk in res["per_year"].items():
        eq(f"per_year.{y}.gamma2_median", blk["channels"]["gamma2"]["median"],
           float(np.median([float(r["coherence_gamma2"]) for r in raw
                            if str(r["acquisition_ts"])[:4] == y])))
        eq(f"per_year.{y}.A_median", blk["channels"]["A"]["median"],
           float(np.median([float(r["coherence_masked_pixels"]) for r in raw
                            if str(r["acquisition_ts"])[:4] == y])))

    # --- the ratios are self-ratios: value / own median, and 1.0 is that median
    eq("ratio.era.identity", all(
        abs(r["ratio"][k] * res["medians"]["era"][k] - r["channels"][k])
        <= 1e-9 * abs(r["channels"][k])
        for r in res["era_rows"] for k in LONG_DIMS), True)
    eq("ratio.record.identity", all(
        abs(r["ratio"][k] * res["medians"]["record"][k] - r["channels"][k])
        <= 1e-9 * abs(r["channels"][k])
        for r in res["record_rows"] for k in core.FEATURES), True)
    eq("ratio.era.median_is_one", abs(float(np.median(
        [r["ratio"]["gamma2"] for r in res["era_rows"]])) - 1.0) < 1e-12, True)
    eq("ratio.record.median_is_one", abs(float(np.median(
        [r["ratio"]["P"] for r in res["record_rows"]])) - 1.0) < 1e-12, True)

    # --- the four stretches tile the axis exactly ---------------------------
    seg = {s["key"]: s for s in res["segments"]}
    order = ("era", "blind_gap", "catalogue_only", "record")
    eq("segments.keys", sorted(seg), sorted(order))
    eq("segments.tile",
       [[seg[k]["start"], seg[k]["end"]] for k in order],
       [[AXIS_START, ERA_END], [shift_day(ERA_END, 1), shift_day(DESCENT_START, -1)],
        [DESCENT_START, shift_day(RECORD_START, -1)], [RECORD_START, AXIS_END]])
    eq("segments.no_gap_no_overlap",
       sum(s["span_days"] for s in res["segments"]),
       span_days(AXIS_START, AXIS_END))
    eq("segments.dims_per_stretch", [len(seg[k]["dims"]) for k in order],
       [len(LONG_DIMS), 0, 0, len(core.FEATURES)])
    eq("segments.filled_days", [seg[k]["filled_days"] for k in order],
       [res["n_era"], 0, 0, res["n_record"]])

    # --- the A-E package's own pinned numbers, at the cut -------------------
    eq("collapse_pair.gamma2",
       [res["record_rows"][-1]["channels"]["gamma2"],
        res["singleton"]["channels"]["gamma2"]],
       [0.019928279168188158, 0.0013393056912106268])
    eq("collapse_pair.A",
       [res["record_rows"][-1]["channels"]["A"],
        res["singleton"]["channels"]["A"]], [60.0, 63.0])
    eq("collapse_pair.days",
       [res["record_rows"][-1]["day"], res["singleton"]["day"]],
       [AXIS_END, FIRST_COLLAPSED_DAY])
    eq("singleton.below_record_median_gamma2",
       res["singleton"]["ratio"]["gamma2"] < 1.0, True)
    return checks


def json_default(o):
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.floating):
        return float(o)
    if isinstance(o, np.bool_):
        return bool(o)
    return str(o)


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--figdir", default=FIGDIR)
    ap.add_argument("--json", default=OUT_JSON)
    ap.add_argument("--md", default=OUT_MD)
    args = ap.parse_args(argv)

    print(f"computing — KDLO, the axis {AXIS_START}..{AXIS_END} in four stretches "
          f"(era {AXIS_START}..{ERA_END}, blind gap {BLIND_GAP_DAYS} d, catalogue "
          f"only, record {RECORD_START}..{AXIS_END}) …")
    res = compute()
    res["checks"] = self_checks(res)
    print(f"  era {res['n_era']} rows ({', '.join(LONG_DIMS)}) · record "
          f"{res['n_record']} rows ({', '.join(core.FEATURES)}) · "
          f"{len(res['checks'])} self-checks, all passed · "
          f"{res['runtime_s']:.1f} s")

    os.makedirs(args.figdir, exist_ok=True)
    fig_series(res, args.figdir)
    os.makedirs(os.path.dirname(os.path.abspath(args.json)), exist_ok=True)
    with open(args.json, "w") as fh:
        json.dump(build_json(res), fh, indent=1, sort_keys=True,
                  default=json_default)
        fh.write("\n")
    print(f"  {os.path.relpath(args.json)}")
    os.makedirs(os.path.dirname(os.path.abspath(args.md)), exist_ok=True)
    with open(args.md, "w") as fh:
        fh.write(report(res))
    print(f"  {os.path.relpath(args.md)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
