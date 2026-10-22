#!/usr/bin/env python3
"""KDLO OCV-Paper — the acquisition series, from the first scene to the collapse.

Question: what do the channels of the echo-mask vector do over the **whole**
catalogue-covered life of the site — from the first ascending IW2 scene over the
KDLO tower (``2015-03-20``) to the collapse (``2022-12-11``)? The A-E package's
record is 60 acquisitions on one descending track (from ``2022-06-02``); this
companion puts that record back into the seven years that precede it and draws
both on **one dated axis**.

What that axis carries, and from where
--------------------------------------
The six dimensions are ``[gamma2, P, D, A, F, S]``, and now **all six exist on
both sides** of the geometry break — but not from one source:

  * ``gamma2`` and the mask size ``A`` are read from the **backend columns** the
    pipeline itself wrote: the era's rows carry ``coherence_gamma2`` and
    ``coherence_masked_pixels`` (162/162), the record's carry the same two
    columns by name (its channel table re-derives them from the strips and
    verifies them, ``kdlo_ocv_channels_csv.py --verify``);
  * ``P``, ``D``, ``F`` and ``S`` are **recomputed from a committed full-dwell
    strip cache with the pipeline's own mask rule** (``kdlo_ocv_masks.py``), and
    there are now two such caches: ``data/kdlo/strips/`` for the record (64
    strips, of which the record's channel table uses 60) and
    ``data/kdlo/strips_era/`` for the era (**162 strips, one per era
    acquisition day**, captured through the very same
    ``decode_slc_rect_at_location(annotation, measurement, lat, lon, 11, 400)``
    window from the burst Id each era row records — the ascending track is
    measured on both VV and VH bursts, so the set carries both, one burst per
    day);
  * the controls ``wind_speed_ms`` / ``temperature_c`` — on both sides.

The era strip cache is not taken on trust: before a single era value of ``P``,
``D``, ``F`` or ``S`` is drawn, this figure **recomputes the mask rule on all 162
era strips and requires the recomputed ``A`` to equal the era export's
``coherence_masked_pixels`` exactly and the recomputed coherence to equal its
``coherence_gamma2`` to ``1e-12``** — i.e. the strips are required to be the
strips the pipeline decoded, the same acceptance the record's channel table runs
(``kdlo_ocv_masks.py``, here pointed at ``strips_era/``). The era's
``P`` has its own basis — the frequency map of the 162 era strips — as the
record's ``P`` has the channel table's 64-strip map: ``P`` is a *frequency of
the mask*, so it is only defined relative to the strip set it is counted over,
and this figure never mixes the two maps (it draws each stretch in units of its
own median, and the basis is constant inside a stretch). The era's set is not one
polarisation — the track is measured on both VV and VH bursts — so what its
frequency map measures is the frequency of the mask across that mixed set; every
other drawn era value is per day, from that day's own burst alone.

The axis has four stretches, and the figure draws all four instead of blurring
them into one series (every date below is read from the committed availability
reference, not asserted from this prose):

  era             2015-03-20 .. 2021-12-19   162 days, ASCENDING IW2, relative
                  orbit 136, catalogue hour 00 — all six dimensions (two from the
                  export's own columns, four recomputed from the era's own strips)
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
break — which is why the blind gap and the catalogue-only stretch stay hatched in
all six panels instead of being crossed by a line.

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
drift). Its rows still carry no phase and no frequency layer — their non-null
columns are exactly the ten of a backend row — so ``P``/``D``/``F``/``S`` can
never come from the export: they come from the era's own strips, and the export's
two columns are what those strips are held to. The four committed inputs are:

  data/kdlo/fig_kdlo_0_raw.json                  the era's decoded backend rows
  data/kdlo/strips_era/manifest.json             the era strip cache: 162 strips,
                                                 one per era day, with the
                                                 Rust-side A / gamma2 cross-check
  data/kdlo/kdlo_ocv_channels.csv                the record's 60 acquisitions
  data/kdlo/reference/kdlo_s1_availability.json  the era's window/days/geometry
                                                 and the export's pinned sha256

The 162 strip payloads and ``capture_requests.json`` are read through that
manifest (every strip is hashed into the artifact), and the strip set is required
to cover exactly the era export's days.

Figure (3 x 2, 183 mm wide), one dated axis 2015-03-20 .. 2022-12-11, **all six**
dimensions drawn over both stretches:
  a  ``gamma2``        log y — era dots (blue) vs. record dots + line (orange),
  b  ``A`` (mask size) log y — each divided by its own stretch median (1.0 drawn);
                              the era's per-year medians are the bars
  c  ``P``             linear y — the record's P comes from the channel table, the
  d  ``D``             log y — era's from the 162 era strips (its own frequency
  e  ``F``             log y — map, stated: P is a mask frequency, so it is only
  f  ``S``             log y — defined against the strip set it is counted over)
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
import kdlo_ocv_masks as masks  # noqa: E402
import kdlo_ocv_stats as st    # noqa: E402

FIGDIR = core.FIGDIR
# The four committed inputs: the era's decoded backend rows (evidence, pinned by
# the availability reference), the era's own full-dwell strip cache (what the era's
# P/D/F/S are recomputed from), the record's channel table, and that reference.
ERA_EXPORT = os.path.join(core.DATA, "fig_kdlo_0_raw.json")
ERA_STRIPS = os.path.join(core.DATA, "strips_era")
ERA_MANIFEST = os.path.join(ERA_STRIPS, "manifest.json")
ERA_REQUESTS = os.path.join(ERA_STRIPS, "capture_requests.json")
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
# The two dimensions whose era value is a column of the era export — the same two
# column names the record carries, so they are the era's two *backend* dimensions
# (and the reason the era's strips can be held to the pipeline's own numbers).
ERA_CHANNEL = {"gamma2": "coherence_gamma2", "A": "coherence_masked_pixels"}
BACKEND_DIMS = ["gamma2", "A"]
# The four the era export does not carry, recomputed from the era's own committed
# full-dwell strips — the same four the record's channel table recomputes from its
# strips, with the same mask rule (``kdlo_ocv_masks.load_vectors``).
STRIP_DIMS = ["P", "D", "F", "S"]
# Panel order (a..f): the backend two first, then the mask geometry, so the panel
# letters keep the meaning they had while the era side of (c)-(f) was empty.
PANELS = BACKEND_DIMS + STRIP_DIMS
# gamma2 and A span decades on both sides; P is a fraction; D, F and S span one to
# two decades inside the record — and two to three inside the era (D reaches 1.0 on
# a one-row mask, F 132, S 208), so the same four axes carry both sides.
PANEL_SCALE = {"gamma2": "log", "A": "log", "P": "linear",
               "D": "log", "F": "log", "S": "log"}
# The era export's ten non-null columns, verbatim, and the phase/frequency columns
# it never fills. Asserting both is what makes "the era export carries no P/D/F/S"
# a measurement — i.e. why those four have to come from the strips and not from the
# export; it says nothing about what the strips can carry.
ERA_NONNULL_COLUMNS = ["acquisition_ts", "coherence_gamma2",
                       "coherence_masked_pixels", "id", "orbit_direction",
                       "pass_label", "status", "temperature_c",
                       "traffic_load_label", "wind_speed_ms"]
COLOUR = {"era": "#4C72B0", "era_med": "#1F4E79", "record": "#DD8452",
          "record_med": "#B2601F", "cut": "#888888", "singleton": "#C44E52",
          "band": "#B0B0B0"}
# The era's band tint: all six panels now carry the era as a measured stretch, so
# the same faint wash sits behind its dots everywhere.
ERA_TINT = 0.055


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


#: Appended to the figure stem; the support-chip run sets it to ``_support``.
FIG_SUFFIX = ""


FIG_DPI = 630


def save_fig(fig, name, out_dir=None):
    """Write one 630 dpi JPEG (the binaries stay out of git)."""
    out_dir = out_dir or FIGDIR
    jpg = os.path.join(out_dir, name + FIG_SUFFIX + ".jpg")
    fig.savefig(jpg, dpi=FIG_DPI, facecolor="white")
    plt.close(fig)
    print(f"  {os.path.relpath(jpg)}")
    return jpg


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


def era_strip_rows(days):
    """The era's own strips, decoded to the four mask-geometry channels.

    One strip per era acquisition day (the cache was captured from the burst Id each
    era row records, through the pipeline's own
    ``decode_slc_rect_at_location(annotation, measurement, lat, lon, 11, 400)``
    window), so the day set must be the era export's day set exactly: no extra day,
    no missing day, one strip each. ``kdlo_ocv_masks.load_vectors`` recomputes the
    pipeline's mask rule on every strip in process; the two columns the export does
    carry are then held against it in ``era_rows``, which is what ties this cache to
    the pipeline rather than to itself.

    ``P`` is the frequency of the mask over *this* strip set — a basis of 162 strips,
    never mixed with the record's strip set: ``P`` is a mask *frequency*, so it is
    only defined against the set it is counted over. Both bases are named wherever a
    ``P`` panel is read.

    The set is not one polarisation: this track is measured on both VV and VH bursts
    over the years, and which burst each day took is what the request file records
    (the manifest must agree with it entry by entry). Every drawn value except ``P``
    is per day and comes from that day's own strip, so a day is internally one
    polarisation; only ``P``'s frequency map is counted across the mixed set, and
    that is said with the panel rather than smuggled into it.
    """
    if not os.path.exists(ERA_MANIFEST):
        _fail(f"no era strip manifest at {_rel(ERA_MANIFEST)}: the era's P/D/F/S "
              "cannot be recomputed (capture it with the kdlo_strip_capture "
              "binary and data/kdlo/strips_era/capture_requests.json)")
    with open(ERA_MANIFEST) as fh:
        man = json.load(fh)
    entries = man.get("entries") or []
    if not entries or man.get("n_entries") != len(entries):
        _fail("the era strip manifest's n_entries disagrees with its own list")
    strip_days = [str(e.get("date"))[:10] for e in entries]
    if len(strip_days) != len(set(strip_days)):
        _fail("the era strip manifest holds more than one strip for a day")
    if sorted(strip_days) != sorted(days):
        extra = sorted(set(strip_days) - set(days))
        missing = sorted(set(days) - set(strip_days))
        _fail(f"the era strip set does not cover the era export: {len(strip_days)} "
              f"strips vs {len(days)} export days (missing {missing[:3]}, extra "
              f"{extra[:3]})")
    if len({e["burst_id"] for e in entries}) != len(entries):
        _fail("the era strip manifest holds more than one strip from one burst")
    # The manifest is what the capture binary wrote; the request file is what it was
    # asked for. They must agree entry by entry on the day, the burst and the
    # polarisation — the era is *not* one polarisation (the ascending passes carry
    # both), and which burst was measured is exactly what the request file records.
    with open(ERA_REQUESTS) as fh:
        req = json.load(fh)
    asked = [(str(e["date"])[:10], e["burst_id"], str(e["polarisation"]))
             for e in (req.get("entries") or [])]
    got = [(str(e["date"])[:10], e["burst_id"], str(e["polarisation"]))
           for e in entries]
    if asked != got:
        _fail("the era strip manifest does not match the capture requests it was "
              "written from (date / burst_id / polarisation differ)")
    # core.ERA_ORBIT is "<direction>|<relative orbit>"; the strips must sit on
    # exactly that geometry, or the era's strips are not the era's track.
    era_dir, era_rel = core.ERA_ORBIT.split("|")
    geom = {(str(e["orbit_direction"]), str(e["subswath"]),
             int(e["relative_orbit"])) for e in entries}
    if geom != {(era_dir, core.ERA_SWATH, int(era_rel))}:
        _fail(f"the era strips sit on {sorted(geom)}, expected the era's own "
              f"{era_dir} {core.ERA_SWATH} relative orbit {era_rel}")
    rows, meta = masks.load_vectors(ERA_STRIPS, ERA_MANIFEST)
    by_day = {r["date"]: r for r in rows}
    if sorted(by_day) != sorted(days):
        _fail("the decoded era vectors are not one per era day")
    for day, row in by_day.items():
        if tuple(row["shape"]) != tuple(meta["majority_shape"]):
            _fail(f"the era strip of {day} has shape {row['shape']}, the set's "
                  f"majority shape is {meta['majority_shape']}")
        if row.get("P") is None:
            _fail(f"the era strip of {day} carries no P (no majority-shape mask)")
    if meta.get("n_without_mask") or meta.get("dates_without_mask"):
        _fail(f"era strips without a mask: {meta.get('dates_without_mask')!r}")
    pols = {}
    for e in entries:
        pols[str(e["polarisation"])] = pols.get(str(e["polarisation"]), 0) + 1
    meta["polarisations"] = pols
    meta["n_requests"] = req.get("n_entries")
    meta["manifest"] = _rel(ERA_MANIFEST)
    meta["requests"] = _rel(ERA_REQUESTS)
    meta["manifest_sha256"] = sha(ERA_MANIFEST)
    meta["requests_sha256"] = sha(ERA_REQUESTS)
    meta["strip_bytes_total"] = sum(
        os.path.getsize(os.path.join(ERA_STRIPS, str(e["strip_file"])))
        for e in entries)
    meta["p_basis"] = (f"the frequency map of these {len(entries)} era strips "
                       "(never mixed with the record's strip set; the set carries "
                       "both polarisations, "
                       + " and ".join(f"{n} {p}" for p, n in sorted(pols.items()))
                       + " days, which is what its frequency map is counted over)")
    return by_day, meta


def era_rows(av):
    """The era's six channels: two backend columns + four from its own strips.

    The checks run in the order in which they can fail: the export file is the one
    the reference pins (by sha256), it holds the measured days of ``figure_0``, its
    days are the era's own window, and its **non-null columns are exactly the ten
    of a backend row without a phase or frequency layer** — which is why the four
    mask-geometry dimensions cannot come from it. They come from the era's own
    committed strip cache (``data/kdlo/strips_era/``), and ``era_strip_rows`` holds
    that cache to the export day by day: the mask rule recomputed on each strip must
    reproduce ``coherence_masked_pixels`` exactly and ``coherence_gamma2`` to
    ``1e-12``.

    The returned rows carry all six channels: ``gamma2`` and ``A`` are the export's
    own values — so panels (a) and (b) are numerically unchanged by the strips'
    arrival — and ``P``/``D``/``F``/``S`` are the per-day vectors from the strips.
    The metadata-anomaly day is checked on its slot.
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
    strips, strip_meta = era_strip_rows(days)
    out = []
    for r in rows:
        day = str(r["acquisition_ts"])[:10]
        g, a = r.get(ERA_CHANNEL["gamma2"]), r.get(ERA_CHANNEL["A"])
        if g is None or a is None:
            _fail(f"the era row {day} has no {ERA_CHANNEL['gamma2']} or no "
                  "mask size")
        s = strips[day]
        # The strips are held to the pipeline's own two numbers, day by day: the
        # mask rule recomputed on the strip must reproduce both.
        if s["A"] != int(a) or s["rust_masked_pixels"] != int(a):
            _fail(f"the era strip of {day} recomputes A={s['A']} (manifest "
                  f"{s['rust_masked_pixels']}), the pipeline recorded {int(a)}")
        if (abs(s["gamma2"] - float(g)) > masks.TOL_GAMMA2
                or abs(float(s["rust_gamma2"]) - float(g)) > masks.TOL_GAMMA2):
            _fail(f"the era strip of {day} recomputes gamma2={s['gamma2']!r} "
                  f"(manifest {s['rust_gamma2']!r}), the pipeline recorded "
                  f"{float(g)!r}")
        out.append({"day": day, "stretch": "era", "gamma2": float(g),
                    "A": float(a), "P": s["P"], "D": s["D"], "F": s["F"],
                    "S": s["S"], "mask_shape": list(s["shape"]),
                    "polarisation": s["polarisation"],
                    "strip_sha256": s["strip_sha256"],
                    "wind_speed_ms": r.get("wind_speed_ms"),
                    "temperature_c": r.get("temperature_c"),
                    "pass_label": r.get("pass_label"),
                    "hour": str(r["acquisition_ts"])[11:13]})
    out.sort(key=lambda r: r["day"])
    strip_meta["days_held_to_the_export"] = len(out)
    strip_meta["tol_gamma2"] = masks.TOL_GAMMA2
    return out, strip_meta


def record_rows(csv_rows):
    """The record's acquisitions inside the axis: the 16 pre-collapse days.

    All six dimensions must be complete on every one of them: while the era side of
    panels (c)-(f) was empty a patchy row could be read as "the era and the record
    both lack it", and that reading must not become possible now that both sides are
    drawn.
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
         "dims": list(PANELS),
         "strips": len(era),
         "p_basis": f"the {len(era)} era strips of this figure "
                    "(data/kdlo/strips_era/), both polarisations in the set",
         "note": f"{len(era)} measured days and all six dimensions: gamma2 and the "
                 "mask size are the era export's own columns, P/D/F/S are "
                 "recomputed from the era's own committed full-dwell strips, each "
                 "day held to those two columns before anything is drawn. The track "
                 "is measured on both VV and VH bursts, so the set carries both "
                 "polarisations; each day is one burst, and only P's frequency map "
                 "is counted across the mixed set"},
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
         "p_basis": "the strip set the channel table recomputes from "
                    "(data/kdlo/strips/, its own provenance in "
                    "kdlo_ocv_channels_meta.json)",
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


def compute(table=None):
    """The era rows, the record rows, the per-stretch medians, the segments.

    The medians are the 1.0 of the figure: **each stretch is divided by its own
    median** (the era's 162-day median for the era, the 16-day median of
    2022-06-02..2022-12-11 for the record), so no level is compared across the
    geometry break — and, for ``P``, no *basis* either: the era's ``P`` is the mask
    frequency over the era's 162 strips, the record's over the strip set the channel
    table recomputes from, and the two maps are never added, averaged or rescaled
    into one another. The record's first month (2022-06) is kept separately: it is
    the baseline the months companion uses, and the cross-read of the two
    artifacts hangs on it.
    """
    table = table or CHANNEL_CSV
    t0 = time.time()
    av = load_availability()
    era, era_strips_meta = era_rows(av)
    csv_rows = core.load_csv(table)
    record = record_rows(csv_rows)
    single = singleton_row(csv_rows)
    later = sorted(core.day_of(r) for r in csv_rows if core.day_of(r) > AXIS_END)

    groups = {"era": era,
              "record": record,
              "record_first_month": [r for r in record
                                     if core.day_of(r)[:7] == RECORD_START[:7]]}
    dims_of = {"era": list(PANELS), "record": list(core.FEATURES),
               "record_first_month": list(core.FEATURES)}
    med = {}
    for g, rows in groups.items():
        med[g] = {}
        for k in dims_of[g]:
            vals = [r[k] for r in rows]
            med[g][k] = st._stats(vals)["median"] if vals else None
    stats = {"era": {k: _median(era, k) for k in PANELS},
             "record": {k: _median(record, k) for k in core.FEATURES}}

    era_out = [{"day": r["day"], "stretch": "era",
                "pass_label": r["pass_label"], "hour": r["hour"],
                "wind_speed_ms": r["wind_speed_ms"],
                "temperature_c": r["temperature_c"],
                "polarisation": r["polarisation"],
                "strip_sha256": r["strip_sha256"],
                "mask_shape": r["mask_shape"],
                "channels": {k: r[k] for k in PANELS},
                "ratio": {k: r[k] / med["era"][k] for k in PANELS}}
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
                                 for k in PANELS}}
                for y, rs in sorted(per_year.items())}

    sources = {
        "era_export": {"path": _rel(ERA_EXPORT), "sha256": sha(ERA_EXPORT),
                       "rows": len(era),
                       "pinned_by": "availability.pipeline."
                                    "figure_0_raw_export.sha256"},
        "era_strips": {"path": _rel(ERA_MANIFEST),
                       "sha256": sha(ERA_MANIFEST),
                       "requests": _rel(ERA_REQUESTS),
                       "requests_sha256": sha(ERA_REQUESTS),
                       "strips": era_strips_meta["n_strips"],
                       "bytes": era_strips_meta["strip_bytes_total"],
                       "n_with_mask": era_strips_meta["n_with_mask"],
                       "majority_shape": era_strips_meta["majority_shape"],
                       "polarisations": era_strips_meta["polarisations"],
                       "p_basis": era_strips_meta["p_basis"],
                       "accepted_against": "era_export coherence_gamma2 / "
                                           "coherence_masked_pixels "
                                           f"({era_strips_meta['days_held_to_the_export']}"
                                           f"/{era_strips_meta['n_strips']})"},
        "channel_table": {"path": _rel(table), "sha256": sha(table),
                          "chip": core.chip_of(table)},
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
        "era_strips": era_strips_meta,
        "sources": sources,
        "n_rows": len(era) + len(record),
        "runtime_s": time.time() - t0,
    }


# ---------------------------------------------------------------------------
# The figure
# ---------------------------------------------------------------------------
def draw_stretches(ax):
    """The four stretches of the axis, drawn — never interpolated across.

    Every panel carries the same four stretches, because all six dimensions now
    exist on the era: the era is the same faint tint in panels (a)-(f) and its days
    are drawn by ``draw_dimension`` in each of them. The two hatches that remain are the
    *empty* stretches of the catalogue — the blind gap and the descending stretch no
    committed row covers — and they are empty in all six panels for the same reason.
    """
    ax.set_xlim(D(AXIS_START), D(FIRST_COLLAPSED_DAY))
    # the era: tinted, and drawn with its own dots, in all six panels
    ax.axvspan(D(AXIS_START), D(ERA_END), facecolor=COLOUR["era"],
               alpha=ERA_TINT, lw=0.0, zorder=0.5)
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


def draw_dimension(ax, res, dim):
    """Any one of the six dimensions: era dots + year bars, record line + dots.

    Both stretches are divided by their **own** median, so the two colours are two
    levels, not one series — and in panel (c) two ``P`` bases as well: the era's mask
    frequency is counted over the era's 162 strips, the record's over the strip set
    the channel table recomputes from, and the two maps never meet. The era's year
    medians are drawn as bars over each year's own day range: that is an aggregation
    *inside* the era (and the only one here), never a bridge towards the record.
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


def fig_series(res, out_dir=None):
    """The 3 x 2 companion: one dated axis, four stretches, six dimensions.

    All six panels get the identical treatment — ``draw_stretches`` +
    ``draw_dimension`` — because both stretches now carry all six dimensions: the
    only difference between panels is the y scale, and the only non-data region is
    the two empty stretches of the catalogue, hatched the same way in all six.
    """
    fig, axes = plt.subplots(3, 2, figsize=(COL_WIDTH, COL_WIDTH * 0.64))
    flat = list(np.ravel(axes))
    for letter, dim, ax in zip("abcdef", PANELS, flat):
        draw_stretches(ax)
        draw_dimension(ax, res, dim)
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
                 "collapse", x=0.004, y=1.022, ha="left", fontsize=10.0)
    # Two short lines rather than one long one: at 7 pt a single 250-character
    # line is ~300 mm wide, and ``savefig.bbox = "tight"`` would then widen the
    # canvas to fit it, leaving the 3 x 2 grid stranded in the left two thirds.
    fig.text(0.004, 0.998,
             "each stretch in units of its own median (log where a dimension "
             f"spans decades): {res['n_era']} ascending-IW2 days,\n"
             "their P/D/F/S from the era's own committed strips, then the "
             f"{res['n_record']} descending-IW1 days of the record — pooled "
             "across the break nowhere",
             ha="left", va="top", fontsize=7.0, color="#444444",
             linespacing=1.4)
    fig.tight_layout(pad=0.5, h_pad=1.1, w_pad=1.3, rect=(0, 0, 1, 0.928))
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
        "dims": {
            "panel_order": list(PANELS),
            "canonical_order": list(core.FEATURES),
            "from_backend_columns": BACKEND_DIMS,
            "from_era_strips": STRIP_DIMS,
            "note": "all six dimensions are drawn on both stretches; gamma2 and A "
                    "are the same two backend columns on either side, while "
                    "P/D/F/S are recomputed from a committed full-dwell strip "
                    "cache with the pipeline's own mask rule — the channel "
                    "table's strip set for the record, the era's own 162-strip "
                    "set (data/kdlo/strips_era/) for the era",
        },
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
            "p_basis": "P is a mask frequency, so it is only defined against the "
                       "strip set it is counted over: the era's P over its own "
                       f"{res['n_era']} strips, the record's over the strip set the "
                       "channel table recomputes from. The two frequency maps are "
                       "never added, averaged or used to rescale one another — each "
                       "stretch is drawn as its own value/median, and one basis is "
                       "constant inside one stretch. Within the era's basis the "
                       "strips are not one polarisation (the track carries both VV "
                       "and VH), so its frequency map is a mixed-set frequency; "
                       "each era day's own four values come from that day's own "
                       "burst alone",
        },
        "segments": res["segments"],
        "geometry_break": res["geometry_break"],
        "availability": res["availability"],
        "era_strips": res["era_strips"],
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
    md.append("All six dimensions exist on **both** dated stretches — the era's "
              f"{res['n_era']} days and the record's {res['n_record']} — and no "
              "panel is blank over a stretch that has data. `gamma2` and the mask "
              "size `A` are the same two backend columns on either side "
              "(`coherence_gamma2`, `coherence_masked_pixels`, which the era's own "
              "export carries). `P`, `D`, `F` and `S` are recomputed from a "
              "committed full-dwell strip cache with the pipeline's own mask rule: "
              "from the strip set the channel table uses for the record, from the "
              "era's own strips for the era, each era day held against the export's "
              "`A` and `gamma2` before it is drawn. The two hatched stretches carry "
              "no dimension at all: neither track has a scene there.")
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
    md.append(f"Each stretch is divided by **its own** median. Era "
              f"({res['n_era']} days): "
              + ", ".join(f"`{k}` {_fmt(med['era'][k])}" for k in PANELS)
              + f". Record ({res['n_record']} days): "
              + ", ".join(f"`{k}` {_fmt(med['record'][k])}" for k in core.FEATURES)
              + ". The 1.0 line is therefore **two different levels in every "
              "panel**, and the two colours are two series of self-ratios, not one "
              "series with a step in it. Levels are not comparable across the break "
              "by construction (different track, different swath, different pass); "
              "what each stretch *does* is.")
    md.append("")
    md.append("`P` needs one rule more, because it is not a per-day quantity: it is "
              "the frequency of the mask over a *set* of strips, so it is only "
              "defined against the set it is counted over. The era's `P` is counted "
              f"over its own {res['n_era']} strips, the record's over the strip set "
              "the channel table recomputes from — one basis per stretch, constant "
              "inside it, and the two frequency maps are never added, averaged or "
              "used to rescale one another. The era's basis is not a single "
              "polarisation (that track carries both), so its map is a mixed-set "
              "frequency; the JSON names it as such. `D`, `F` and `S` are per-day "
              "and need no such rule. Log axes where a dimension spans decades, so a "
              "value below "
              "1 stays visible.")
    md.append("")
    md.append("## The record inside the axis, day by day")
    md.append("")
    body = []
    for r in res["record_rows"]:
        body.append([r["day"], _fmt(r["wind_speed_ms"], 3),
                     _fmt(r["temperature_c"], 4)]
                    + [_fmt(r["channels"][k]) for k in core.FEATURES]
                    + [_fmt(r["ratio"][k], 3) for k in PANELS])
    md += md_table(["day", "wind m/s", "T degC"] + core.FEATURES
                   + [k + " /med" for k in PANELS], body)
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
                    + [_fmt(blk["channels"][k]["median"]) for k in PANELS]
                    + [_fmt(blk["channels"][k]["median"] / med["era"][k], 3)
                       for k in PANELS])
    md += md_table(["year", "days", "first … last"] + [k + " median" for k in PANELS]
                   + [k + " /era med" for k in PANELS], body)
    md.append("")
    md.append("These year medians are the only aggregation drawn over the era "
              "(bars in all six panels), and they stay **inside** the era: a year "
              "is never compared with the record. They are era-track levels, and for "
              "`P` they stay on the era's own strip basis.")
    md.append("")
    md.append("## Where the era's `P`, `D`, `F` and `S` come from")
    md.append("")
    md.append("The era export carries none of them: its ten non-null columns are "
              "`acquisition_ts`, `coherence_gamma2`, `coherence_masked_pixels`, "
              "`id`, `orbit_direction`, `pass_label`, `status`, `temperature_c`, "
              "`traffic_load_label`, `wind_speed_ms` — no phase layer and no "
              "frequency layer — and that set is asserted against the pinned export "
              "before a value is read from it. The four dimensions are recomputed "
              "here from the era's **own** committed full-dwell strips "
              "(`data/kdlo/strips_era/`, one strip per era acquisition, captured "
              "from the burst Id the export itself records), with the same mask rule "
              "the channel table applies to the record. Nothing is interpolated and "
              "nothing is carried across the break: this is a second strip cache, "
              "not the record's.")
    md.append("")
    md.append("The two columns the export *does* carry are what makes that cache "
              "checkable rather than self-referential: the mask rule recomputed on "
              f"the strips reproduces `coherence_masked_pixels` **exactly** and "
              "`coherence_gamma2` to "
              f"{res['era_strips']['tol_gamma2']:g} on all "
              f"{res['era_strips']['days_held_to_the_export']} of "
              f"{res['era_strips']['n_strips']} era days, so the era's `A` and "
              "`gamma2` panels are the pipeline's own numbers and the four "
              "recomputed dimensions are anchored to them day by day. Every era "
              "strip is "
              f"{res['era_strips']['majority_shape'][0]} x "
              f"{res['era_strips']['majority_shape'][1]}, so the frequency map (`P`) "
              "has one shape over the whole era.")
    md.append("")
    pols = res["era_strips"]["polarisations"]
    md.append("One property of that set has to be read with it: the ascending track "
              "is measured on **both polarisations**, and the era's "
              f"{res['era_strips']['n_strips']} days carry "
              + " and ".join(f"{n} {p}" for p, n in sorted(pols.items()))
              + " bursts. Which burst a day took is what the capture request file "
              f"(`{os.path.basename(ERA_REQUESTS)}`) records, and the manifest is "
              "required to agree with it entry for entry. This is a per-day split, "
              "not a per-pixel one: each day's `gamma2`, `A`, `P`, `D`, `F` and `S` "
              "come from that day's own single burst, so a day is internally one "
              "polarisation and no drawn value mixes the two. What is counted across "
              "the mix is the era's **frequency map** — the per-pixel mask frequency "
              "of the whole era set, which `P` reads — so `P` on the era is a "
              "mixed-set frequency and is labelled as such in the JSON.")
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
              f"The era strip manifest is hashed with its request file, every strip "
              f"it names is read and its size totalled, and the exported `A`/`gamma2` "
              f"of each era day is held against the mask rule recomputed on that "
              f"day's strip — the era cache is never accepted on its own say-so. "
              f"The figure carries **{len(res['checks'])}** self-checks, all of "
              f"which abort before any file is written.")
    md.append("")
    md.append("```bash")
    md.append("python3 code/kdlo/kdlo_ocv_channels_csv.py --verify-only   # the table")
    md.append("python3 code/kdlo/kdlo_ocv_masks.py \\")
    md.append("  --strips-dir data/kdlo/strips_era \\")
    md.append("  --manifest data/kdlo/strips_era/manifest.json  # the era's strips")
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

    # --- the era: its export, and the strips that carry the other four -------
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

    # --- the era's own strips: the four recomputed dimensions -----------------
    # Re-decoded here from the committed cache rather than read off the results: the
    # day set, the mask rule, the four per-day vectors and the frequency map are all
    # recomputed, the two columns the export carries are required to come out of that
    # recomputation (and out of the manifest's own recorded rust_* pair), and the
    # era's cache is required to be disjoint from the record's.
    man = json.load(open(ERA_MANIFEST))
    ent = {e["date"]: e for e in man["entries"]}
    exp = {str(r["acquisition_ts"])[:10]: r for r in raw}
    eq("era_strips.export_days", sorted(exp), days)
    eq("era_strips.manifest_days", sorted(ent), days)
    eq("era_strips.n_entries", man["n_entries"], len(ent))
    eq("era_strips.requests_file",
       os.path.basename(man["requests_file"]), os.path.basename(ERA_REQUESTS))
    eq("era_strips.strip_geometry",
       sorted({(e["rows"], e["cols"]) for e in ent.values()}), [(400, 11)])
    era_dir, era_rel = core.ERA_ORBIT.split("|")
    eq("era_strips.orbit_swath",
       sorted({(str(e["orbit_direction"]), str(e["subswath"]),
                int(e["relative_orbit"])) for e in ent.values()}),
       [(era_dir, core.ERA_SWATH, int(era_rel))])
    # The era is not one polarisation: this track carries both, and the request file
    # is what says which burst each day took. The manifest is the capture binary's
    # answer to those requests, so the two must agree entry for entry — that pairing,
    # with the pair of values held to the export above, is what ties a strip to the
    # row it stands under.
    req = json.load(open(ERA_REQUESTS))
    eq("era_strips.requests_match_manifest",
       [(str(e["date"])[:10], e["burst_id"], str(e["polarisation"]))
        for e in req["entries"]],
       [(str(e["date"])[:10], e["burst_id"], str(e["polarisation"]))
        for e in man["entries"]])
    eq("era_strips.n_requests", req["n_entries"], man["n_entries"])
    eq("era_strips.one_burst_per_day",
       len({e["burst_id"] for e in ent.values()}), len(ent))
    eq("era_strips.polarisations",
       sorted({str(e["polarisation"]) for e in ent.values()}), ["VH", "VV"])
    eq("era_strips.polarisation_counts",
       res["era_strips"]["polarisations"],
       {p: sum(1 for e in ent.values() if str(e["polarisation"]) == p)
        for p in ("VV", "VH")})
    eq("era_strips.sentinel_platforms",
       sorted({str(e["burst_name"])[:3] for e in ent.values()}
              - {"S1A", "S1B"}), [])
    eq("era_strips.strip_bytes", res["era_strips"]["strip_bytes_total"],
       sum(e["strip_bytes"] for e in ent.values()))
    eq("era_strips.strip_sha256_per_file",
       [e["strip_sha256"] for e in ent.values()],
       [sha(os.path.join(ERA_STRIPS, e["strip_file"])) for e in ent.values()])
    eq("era_strips.not_the_record_cache",
       os.path.abspath(ERA_STRIPS) != os.path.abspath(masks.STRIPS_DIR), True)
    eq("era_strips.disjoint_from_record_strips",
       sorted(set(e["strip_sha256"] for e in ent.values())
              & {e.get("strip_sha256") for e in
                 masks.load_manifest()["entries"]}), [])
    eq("era_strips.manifest_rust_A_is_export_A",
       [ent[d]["rust_masked_pixels"] for d in days],
       [int(exp[d]["coherence_masked_pixels"]) for d in days])
    eq("era_strips.manifest_rust_gamma2_is_export_gamma2",
       [ent[d]["rust_gamma2"] for d in days],
       [float(exp[d]["coherence_gamma2"]) for d in days])
    vec_rows, vec_meta = masks.load_vectors(ERA_STRIPS, ERA_MANIFEST)
    vec = {r["date"]: r for r in vec_rows}
    eq("era_strips.recomputed_A_is_export_A", [vec[d]["A"] for d in days],
       [int(exp[d]["coherence_masked_pixels"]) for d in days])
    eq("era_strips.recomputed_gamma2_is_export_gamma2",
       all(abs(vec[d]["gamma2"] - float(exp[d]["coherence_gamma2"]))
           <= masks.TOL_GAMMA2 for d in days), True)
    eq("era_strips.majority_shape",
       res["era_strips"]["majority_shape"], vec_meta["majority_shape"])
    eq("era_strips.no_strip_without_a_mask",
       [vec_meta["n_without_mask"], vec_meta["dates_without_mask"]], [0, []])
    drawn = {r["day"]: r for r in res["era_rows"]}
    eq("era_strips.drawn_days", sorted(drawn), days)
    eq("era_strips.drawn_A_is_export_A",
       [drawn[d]["channels"]["A"] for d in days],
       [float(exp[d]["coherence_masked_pixels"]) for d in days])
    eq("era_strips.drawn_gamma2_is_export_gamma2",
       [drawn[d]["channels"]["gamma2"] for d in days],
       [float(exp[d]["coherence_gamma2"]) for d in days])
    eq("era_strips.drawn_P_is_the_recomputed_frequency",
       [drawn[d]["channels"]["P"] for d in days],
       [vec[d]["P"] for d in days])
    eq("era_strips.drawn_D_F_S_are_the_recomputation",
       [[drawn[d]["channels"][k] for k in ("D", "F", "S")] for d in days],
       [[vec[d][k] for k in ("D", "F", "S")] for d in days])
    eq("era_strips.six_channels_on_every_day",
       sorted({k for d in days for k in PANELS
               if drawn[d]["channels"][k] is None}), [])
    eq("era_strips.mask_shape_uniform",
       sorted({tuple(drawn[d]["mask_shape"]) for d in days}), [(400, 11)])
    eq("era_strips.p_basis_is_the_era_set",
       res["segments"][0]["strips"], len(days))
    eq("era_strips.P_in_unit_interval",
       all(0.0 <= drawn[d]["channels"]["P"] <= 1.0 for d in days), True)
    eq("era_strips.D_in_unit_interval",
       all(0.0 < drawn[d]["channels"]["D"] <= 1.0 for d in days), True)
    eq("era_strips.F_at_least_one",
       all(int(drawn[d]["channels"]["F"]) >= 1 for d in days), True)
    eq("era_strips.S_finite_non_negative",
       all(np.isfinite(drawn[d]["channels"]["S"])
           and drawn[d]["channels"]["S"] >= 0.0 for d in days), True)
    eq("era_strips.A_at_the_mask_minimum",
       min(drawn[d]["channels"]["A"] for d in days) >= masks.MIN_N_MASKED, True)
    eq("era_strips.drawn_polarisation_is_the_manifest",
       [drawn[d]["polarisation"] for d in days],
       [str(ent[d]["polarisation"]) for d in days])
    eq("era_strips.days_held_to_the_export",
       res["era_strips"]["days_held_to_the_export"], len(days))

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
    for k in STRIP_DIMS:
        eq(f"median.era.{k}", res["medians"]["era"][k],
           float(np.median([r["channels"][k] for r in res["era_rows"]])))
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
        for r in res["era_rows"] for k in PANELS), True)
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
       [len(PANELS), 0, 0, len(core.FEATURES)])
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
    ap.add_argument("--json", default=None)
    ap.add_argument("--md", default=None)
    ap.add_argument("--table", default=None,
                    help="chip table to read; committed by default")
    ap.add_argument("--suffix", default=None,
                    help="appended to every output stem")
    args = ap.parse_args(argv)
    chip = core.chip_of(args.table)
    if chip is None:
        ap.error(f"unknown chip table: {args.table}")
    suffix = args.suffix if args.suffix is not None else         ("_support" if chip == "support" else "")
    global FIG_SUFFIX
    FIG_SUFFIX = suffix
    args.json = args.json or os.path.join(core.DATA,
                                          f"fig_kdlo_ocv_series{suffix}.json")
    args.md = args.md or os.path.join(FIGDIR, f"fig_kdlo_ocv_series{suffix}.md")

    print(f"computing — KDLO, the axis {AXIS_START}..{AXIS_END} in four stretches "
          f"(era {AXIS_START}..{ERA_END}, blind gap {BLIND_GAP_DAYS} d, catalogue "
          f"only, record {RECORD_START}..{AXIS_END}) …")
    res = compute(args.table)
    if chip == "support":
        # The record-side self-checks pin the committed chip's own numbers; on the
        # support chip the table is instead checked by the carry-over audit in the
        # paper run (300/300 committed values) and by --verify of the migration.
        res["checks"] = {}
        print("  self-checks: skipped on the support chip (they pin the committed "
              "numbers; see the paper run's carry-over audit)")
    else:
        res["checks"] = self_checks(res)
    print(f"  era {res['n_era']} rows ({', '.join(PANELS)} — gamma2/A from the "
          f"export, P/D/F/S from the era's {res['era_strips']['n_strips']} strips, "
          f"held to the export on all "
          f"{res['era_strips']['days_held_to_the_export']}) · record "
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
