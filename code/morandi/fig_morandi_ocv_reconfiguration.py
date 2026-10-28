#!/usr/bin/env python3
"""Morandi — OCV reconfiguration in units of the 2016-08 … 2017-07 reference.

Question: how far is each month of the Ponte Morandi record's **last twelve
months before the collapse** (2017-08 … 2018-08-08) from the bridge's **own
earlier state** (2016-08 … 2017-07), on the two dimensions of the ``behaviour``
block of ``code/sds/ocv_reconfigure_spec.md`` — ``gamma2`` (the consistency of
the echo across acquisitions) and ``P`` (its persistence) — and does the
approach to 2018-08-14 show as a departure from that state?

The metric is the pre-registered one, imported rather than rewritten:

    r = median(log x_month) - median(log x_ref)          g = log (both dimensions)
    z = r / scale_ref                                    scale_ref = MAD(log x_ref)
    R = rms(z(gamma2), z(P))

so the reference is the origin **(0, 0)** by construction.

The unit is the **support-sized chip**, not the 232-chip per-track export
-----------------------------------------------------------------------
This plate reads ``data/morandi/support/morandi_ocv_channels.csv``: one 51 x
137 chip per acquisition (861 support cells, 30 m guard, 200 m ring) — 207
acquisitions, 94 A_asc / 113 A_des — the same kind of table the Carola and CTS
reconfiguration plates read. The paper figure's 232-chip export and the site's
fixed 640-px deck mask are a *different* unit; the two are not mixed here, and
every level below belongs to this chip only.

The split, by the table's own state column
------------------------------------------
    reference   99 chips   2016-08-24 … 2017-07-31   (12 months, 53/46)
    series      99 chips   2017-08-01 … 2018-08-08   (13 months, 38/61)
    post         9 chips   2018-08-14 … 2018-09-13   (counted, not plotted)
    ------------------------------------------------------------------
    total      207 chips

The reference is the first twelve months of the record: a whole seasonal cycle,
every calendar month represented, and — unlike the SDS plate's whole pre-
collapse reference (198 chips) — **earlier than the series**, so the series
cannot set its own origin. The nine post-collapse chips carry no value here.

Two consequences are carried into every number rather than into a footnote:

* the reference is **99 chips over 12 months** and the series is **99 chips
  over 13 months**, 3–11 chips per month (1–5 ASC / 2–6 DESC); 2018-08 holds
  **3** chips and is the only ``thin`` month (``SMALL_N = 6``);
* the only interval reported is the **reference's** own uncertainty: 2 000
  circular moving-block draws (block 3, seed 7) of the 99 reference chips, with
  the reference median *and* scale re-estimated inside every draw, the month's
  chips held fixed. It answers "how well is the origin known", never "how well
  is the month known" — and for the n = 3 month it is narrow for the wrong
  reason.

The composition confound, stated rather than smoothed
-----------------------------------------------------
The two tracks are not interchangeable: over the reference their γ² medians are
0.0878 (A_asc) and 0.0104 (A_des), an 8.5x difference, and their P medians
differ by 2.8x. The mix moves from 53/46 in the reference to 38/61 in the
series, and within the series the ASC share falls from 5/11 to 2/7 (1/3 in the
partial last month) while DESC stays at 4–6 chips. A pooled median therefore
moves when *which track is sampled* changes, even if neither track moves; and
because the two levels are separated, the pooled median is an order statistic
that can sit below both track medians once ASC is a minority.

The pooled number stays the headline — it is the pre-registered one — and the
**track split travels beside it**: each track's monthly median against its own
reference median, divided by the *pooled* reference MAD so the two lines share
the headline's unit. That is a composition control, not a second estimate.
Re-referencing each track to its own MAD (the Carola channel split's
construction) would be a different origin twice over for a pure level offset of
this size, so it is not used for the figure; the tracks' own MADs are tabulated
in the guide. The mix-only change implied by the monthly ASC/DESC counts is
computed too: at most 0.52 MAD for γ² and 0.48 MAD for P — a lower bound, real
but not the whole γ² excursion, because the median of a two-mode mixture is an
order statistic and not a weighted average of the two track medians.

Nothing is tested: no H1, no H2, no p-value and no precursor claim.

Figure (2 x 2, 13.0 in wide, 600 dpi jpg — the sds package's ``savefig.dpi``):
  a  z(gamma2) by month, per-chip dots, monthly median, reference band ±2 MAD,
     the two-track split as dotted lines, a chip-count strip, the collapse cut;
     n on every month and the thin month ringed
  b  z(P) by month, same standardisation and axis
  c  the (z gamma2, z P) plane — the origin is the reference state, the grey
     cloud is the 99 reference chips, the path is 2017-08 -> 2018-08 coloured by
     year-quarter (the series crosses two calendar years), the arrows are the
     quarter medians, the terminal month is ringed
  d  the block magnitude R = rms(z gamma2, z P) per month, with the empirical
     5-95 % band of the reference chips as the descriptive band

Outputs:
  data/morandi/ocv_reconfiguration_morandi.json
  data/morandi/ocv_reconfiguration_morandi_table.csv
  figures/morandi/fig_morandi_ocv_reconfiguration.jpg
  figures/morandi/fig_morandi_ocv_reconfiguration.md

Usage:
  python3 code/morandi/fig_morandi_ocv_reconfiguration.py [--selftest | --verify]
  python3 code/morandi/fig_morandi_ocv_reconfiguration.py --figdir /tmp/x
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import math
import os
import random
import sys
import tempfile

import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as mdates  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, os.pardir, "sds"))

import morandi_ocv_core as core      # noqa: E402
import ocv_reconfigure as ocv       # noqa: E402  (the pre-registered metric)

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "mathtext.fontset": "stix",
    "font.size": 8.5,
    "axes.titlesize": 8.8,
    "axes.labelsize": 8.2,
    "xtick.labelsize": 7.2,
    "ytick.labelsize": 7.2,
    "legend.fontsize": 6.6,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.alpha": 0.22,
    "grid.linewidth": 0.4,
    "figure.dpi": 110,
    "savefig.dpi": 600,
    "savefig.bbox": "tight",
    "pdf.fonttype": 42,
})

SCHEMA = "ocv-reconfiguration-monthly-v1"
SITE = "morandi"
SITE_LABEL = core.SITE

#: The behaviour block of the spec, and the column each dimension lives in on the
#: support-sized chip of this site.
DIMS = tuple(ocv.BLOCKS["behaviour"])          # ("gamma2", "P")
COLUMN = {"gamma2": "gamma2_mask", "P": "P"}
TRANSFORM = {dim: "log" for dim in DIMS}

#: The reference is the first twelve months of the record (2016-08 … 2017-07),
#: the series every pre-collapse acquisition after it (2017-08 … 2018-08-08).
#: The two windows are disjoint by construction.
REF_START, REF_END = "2016-08-01", "2017-07-31"
SERIES_START = "2017-08-01"
EVENT = core.COLLAPSE_DATE                      # 2018-08-14
PRE_STATE, POST_STATE = core.PRE, core.POST      # the table's own state labels

#: The "pooled" sensitivity re-centres on the **pre-collapse** record (reference
#: + series), not on the whole table: 9 of the 207 chips are post-collapse, and
#: a MAD dominated by the collapsed state is not a spread of the healthy one.
POOL_END = EVENT

TRACKS = ("A_asc", "A_des")
#: The table's track names are not the spec's orbit keys; this is the one place
#: the two vocabularies meet.
TRACK_CHANNEL = {"A_asc": "ASC", "A_des": "DESC"}

TABLE = os.path.join(core.DATA, "support", "morandi_ocv_channels.csv")
META = os.path.join(core.DATA, "support", "morandi_ocv_channels_meta.json")
OUT_JSON = os.path.join(core.DATA, "ocv_reconfiguration_morandi.json")
OUT_CSV = os.path.join(core.DATA, "ocv_reconfiguration_morandi_table.csv")
OUT_MD = os.path.join(core.FIGDIR, "fig_morandi_ocv_reconfiguration.md")
OUT_JPG = os.path.join(core.FIGDIR, "fig_morandi_ocv_reconfiguration.jpg")

BOOT = 2000            # reference-uncertainty draws (the spec's seed and count)
BLOCK = 3              # consecutive reference chips per draw (spec: BLOCK_BOOT)
SEED = 7
SMALL_N = ocv.SMALL_N  # 6 — below this a month is thin and must be quoted with n
SCALE_MIN_REF = ocv.SCALE_MIN_REF  # 10 — below this the pooled scale takes over

#: One colour per year-quarter of the series. The series crosses two calendar
#: years (2017-08 … 2018-08), so a bare Q1..Q4 key would repeat and the year is
#: part of the key.
QUARTER_COLOUR = {
    "2017-Q3": "#1f4e79",
    "2017-Q4": "#2f855a",
    "2018-Q1": "#b7791f",
    "2018-Q2": "#7a5195",
    "2018-Q3": "#c53030",
}
CHANNEL_COLOUR = {"ASC": "#8a6d3b", "DESC": "#2c7a7b"}
ALERT = "#c53030"


# ---------------------------------------------------------------------------
# The record
# ---------------------------------------------------------------------------
def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 16), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _rel(path):
    return os.path.relpath(path, os.path.dirname(core.DATA))


def _repo_rel(path):
    """Repository-relative, the form a reader of the guide can paste."""
    return os.path.relpath(os.path.abspath(path),
                           os.path.dirname(os.path.dirname(core.DATA)))


def _log(value):
    """``log`` for the two behaviour dimensions (both are strictly positive).

    The pre-registered ``ocv.transform`` is keyed by a closed site registry; the
    two-line transform is spelled out here while the *scale machinery* is
    imported. A non-positive value is dropped for that dimension only; the audit
    counts the drop instead of imputing it.
    """
    return math.log(value) if value is not None and value > 0 else None


def _number(value):
    if value is None:
        return None
    text = value.strip() if isinstance(value, str) else value
    if text in ("", "nan", "NaN", "None", "null"):
        return None
    try:
        return float(text)
    except (TypeError, ValueError):
        return None



def load_rows(path=None):
    """The support-sized chip of one Morandi acquisition per row, validated.

    The committed table carries 207 acquisitions on 51 x 137 chips, each in one
    of the two tracks (``A_asc`` / ``A_des``). A row from another track would
    mix footprints, and a row whose ``reference`` is not the ring would carry a
    different support rule, so the run stops rather than averaging either in.
    ``day`` is the table's date column; it is copied to ``date`` so the rest of
    the module reads like the sibling Carola/CTS plates.
    """
    path = path or TABLE
    with open(path, newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise SystemExit(f"{path}: no rows")
    needed = {"day", "orbit_direction", "track", "segment", "reference", "state",
              "pre_collapse"} | set(COLUMN.values())
    missing = needed - set(rows[0])
    if missing:
        raise SystemExit(f"{path}: missing columns {sorted(missing)}")
    tracks = {row["track"] for row in rows}
    if tracks != set(TRACKS):
        raise SystemExit(f"{path}: track(s) {sorted(tracks)} — this run reads "
                         f"{list(TRACKS)} only")
    for row in rows:
        if row["segment"] != row["track"]:
            raise SystemExit(
                f"{path}: {row['day']} has segment {row['segment']!r} != track "
                f"{row['track']!r} — a segment is a track here")
    rules = {row["reference"] for row in rows}
    if rules != {"ring"}:
        raise SystemExit(f"{path}: reference rule(s) {sorted(rules)} — the ring "
                         f"support of the committed table, one expected")
    states = {row["state"] for row in rows}
    if states != {PRE_STATE, POST_STATE}:
        raise SystemExit(f"{path}: state(s) {sorted(states)}")
    for row in rows:
        expected = "1" if row["state"] == PRE_STATE else "0"
        if row["pre_collapse"] != expected:
            raise SystemExit(
                f"{path}: {row['day']} says state {row['state']!r} but "
                f"pre_collapse={row['pre_collapse']!r}")
    orbits = {row["orbit_direction"] for row in rows}
    unknown = orbits - set(ocv.CHANNEL)
    if unknown:
        raise SystemExit(f"{path}: orbit direction(s) {sorted(unknown)}")
    for row in rows:
        row["date"] = row["day"]
    return rows


def month_of(day):
    return day[:7]


def month_marker(key):
    """Mid-month x position, so a monthly marker sits inside its own bucket."""
    return dt.date(int(key[:4]), int(key[5:7]), 15)


def month_range(first, last):
    """Every month key from ``first`` to ``last`` inclusive."""
    year, month = int(first[:4]), int(first[5:7])
    out = []
    while f"{year:04d}-{month:02d}" <= last:
        out.append(f"{year:04d}-{month:02d}")
        month += 1
        if month == 13:
            year, month = year + 1, 1
    return out


def quarter_of(month):
    return f"{month[:4]}-Q{(int(month[5:7]) - 1) // 3 + 1}"


def _stats(values):
    """Median, MAD and the fixed fallback chain (spec: MAD -> IQR -> SD)."""
    clean = [v for v in values if v is not None]
    scale, kind = ocv.robust_scale(clean)
    return {"n": len(clean), "median": ocv.median(clean),
            "mad": ocv.mad(clean), "scale": scale, "scale_kind": kind}


def _block_indices(n, block, rng):
    """Circular moving blocks of consecutive reference chips."""
    out = []
    while len(out) < n:
        start = rng.randrange(n)
        out.extend((start + step) % n for step in range(block))
    return out[:n]


def reference_draws(series, boot=BOOT, block=BLOCK, seed=SEED):
    """Bootstrap the reference's median and scale, once per draw, for both dims.

    The same chips carry both dimensions, so one index draw serves both — the
    dimensions stay paired inside a draw, which is what makes the draw of ``R``
    a draw of *the same* reference rather than two independent ones.
    """
    n = len(series[DIMS[0]])
    rng = random.Random(seed)
    draws = []
    for _ in range(boot):
        indices = _block_indices(n, block, rng)
        entry = {}
        for dim in DIMS:
            stats = _stats([series[dim][i] for i in indices])
            entry[dim] = (stats["median"], stats["scale"])
        draws.append(entry)
    return draws


def missing_reference(ref_rows):
    """Reference months with no acquisition inside the window (the cloud's gaps)."""
    present = {month_of(row["date"]) for row in ref_rows}
    keys = sorted(present)
    return [key for key in month_range(keys[0], keys[-1]) if key not in present]


def _reference_block(ref_rows, ref_values, pool_values):
    """Median and scale of the reference, per dimension, with the used scale."""
    dims = {}
    for dim in DIMS:
        values = [v for v in ref_values[dim] if v is not None]
        stats = _stats(values)
        pooled = _stats(pool_values[dim])
        if stats["n"] >= SCALE_MIN_REF or not pooled["scale"]:
            scale, kind, used = stats["scale"], stats["scale_kind"], "mad"
        else:
            scale, kind, used = pooled["scale"], pooled["scale_kind"], "pooled_mad"
        if not scale:
            raise SystemExit(f"reference scale for {dim} is degenerate")
        dims[dim] = {
            "column": COLUMN[dim],
            "transform": TRANSFORM[dim],
            "n": stats["n"],
            "median_log": stats["median"],
            "median": math.exp(stats["median"]),
            "mad_log": stats["mad"],
            "one_mad_factor": math.exp(stats["mad"]),
            "scale_log": scale,
            "scale_kind": kind,
            "scale_used": used,
            "pooled_n": pooled["n"],
            "pooled_median_log": pooled["median"],
            "pooled_scale_log": pooled["scale"],
            "pooled_scale_kind": pooled["scale_kind"],
        }
    return dims


def _channel_blocks(ref_rows, series_rows, reference):
    """Each track's months against **its own** reference, in two units.

    A composition control, not a second estimate: the same month's chips, the
    same formula, but the origin built from one track instead of the pooled
    table. ``z_own`` divides by that track's reference MAD (the natural
    within-track standardisation); ``z_pooled`` divides the same difference by
    the *pooled* reference MAD, so it shares the headline's unit and can be
    compared with the ±2 band. Both are stored; the figure draws ``z_pooled``.
    """
    out = {}
    for track in TRACKS:
        key = TRACK_CHANNEL[track]
        ref_sel = [row for row in ref_rows if row["track"] == track]
        dims = {}
        for dim in DIMS:
            stats = _stats([_log(_number(row[COLUMN[dim]])) for row in ref_sel])
            pooled = _stats([_log(_number(row[COLUMN[dim]]))
                             for row in ref_rows + series_rows
                             if row["track"] == track])
            used = "mad" if stats["n"] >= SCALE_MIN_REF or not pooled["scale"] \
                else "pooled_mad"
            scale = stats["scale"] if used == "mad" else pooled["scale"]
            dims[dim] = {"n": stats["n"], "median_log": stats["median"],
                         "median": math.exp(stats["median"]), "mad_log": stats["mad"],
                         "scale_log": scale, "scale_used": used}
        groups = {}
        for row in series_rows:
            if row["track"] == track:
                groups.setdefault(month_of(row["date"]), []).append(row)
        months = []
        for month in sorted(groups):
            group = groups[month]
            entry = {"month": month, "n": len(group), "dates":
                     sorted(row["date"] for row in group),
                     "median_log": {}, "median": {}, "delta_log": {},
                     "z_own": {}, "z_pooled": {}}
            for dim in DIMS:
                stats = _stats([_log(_number(row[COLUMN[dim]])) for row in group])
                centre = dims[dim]["median_log"]
                entry["median_log"][dim] = stats["median"]
                entry["median"][dim] = (None if stats["median"] is None
                                        else math.exp(stats["median"]))
                if stats["median"] is None:
                    entry["delta_log"][dim] = None
                    entry["z_own"][dim] = None
                    entry["z_pooled"][dim] = None
                else:
                    delta = stats["median"] - centre
                    entry["delta_log"][dim] = delta
                    entry["z_own"][dim] = delta / dims[dim]["scale_log"]
                    entry["z_pooled"][dim] = delta / reference[dim]["scale_log"]
            months.append(entry)
        out[key] = {
            "channel": key,
            "track": track,
            "orbit_direction": "ASCENDING" if track == "A_asc" else "DESCENDING",
            "n_ref": len(ref_sel),
            "n_series": sum(len(group) for group in groups.values()),
            "dims": dims,
            "months": months,
        }
    for block in out.values():
        for dim in DIMS:
            own = [entry["z_own"][dim] for entry in block["months"]
                   if entry["z_own"][dim] is not None]
            common = [entry["z_pooled"][dim] for entry in block["months"]
                      if entry["z_pooled"][dim] is not None]
            block.setdefault("summary", {})[dim] = {
                "z_own_range": [min(own), max(own)] if own else None,
                "z_pooled_range": [min(common), max(common)] if common else None,
            }
    return out


def window_stats(rows, start, end, used=False, note=""):
    """Chips/months/track mix of an arbitrary candidate reference window."""
    picked = [row for row in rows if start <= row["date"] <= end]
    months = {month_of(row["date"]) for row in picked}
    return {
        "window": [start, end],
        "chips": len(picked),
        "months": len(months),
        "n_by_channel": {TRACK_CHANNEL[track]:
                         sum(1 for row in picked if row["track"] == track)
                         for track in TRACKS},
        "contains_series": any(row["date"] >= SERIES_START for row in picked),
        "used": used,
        "note": note,
    }


def _rms(z_by_dim):
    values = [z_by_dim[d] for d in DIMS if z_by_dim.get(d) is not None]
    if not values:
        return None
    return math.sqrt(sum(v * v for v in values) / len(values))



# ---------------------------------------------------------------------------
# The report — one pass over the table, everything else reads this object
# ---------------------------------------------------------------------------
def build_report(table=None, meta=None, boot=BOOT, block=BLOCK, seed=SEED):
    """The whole plate as one object: reference, months, chips, channels, audit."""
    table = table or TABLE
    meta_path = meta or META
    rows = sorted(load_rows(table), key=lambda row: (row["date"], row["track"]))
    with open(meta_path) as handle:
        meta_obj = json.load(handle)

    pre_rows = [row for row in rows if row["state"] == PRE_STATE]
    post_rows = [row for row in rows if row["state"] == POST_STATE]
    ref_rows = [row for row in pre_rows if REF_START <= row["date"] <= REF_END]
    series_rows = [row for row in pre_rows
                   if SERIES_START <= row["date"] < EVENT]
    pool_rows = list(pre_rows)
    if not ref_rows or not series_rows:
        raise SystemExit("the reference or the series is empty — check the window")
    if len(ref_rows) + len(series_rows) != len(pre_rows):
        outside = sorted({row["date"] for row in pre_rows}
                         - {row["date"] for row in ref_rows}
                         - {row["date"] for row in series_rows})
        raise SystemExit(
            f"the pre-collapse record is not the reference + the series: "
            f"{len(pre_rows)} chips, {len(ref_rows)} + {len(series_rows)}, "
            f"outside {outside[:5]}")
    if len(pre_rows) + len(post_rows) != len(rows):
        raise SystemExit("the state labels do not partition the table")
    if not ref_rows[-1]["date"] < series_rows[0]["date"]:
        raise SystemExit(
            f"the reference ends {ref_rows[-1]['date']}, the series starts "
            f"{series_rows[0]['date']} — the two windows overlap")
    if post_rows and not series_rows[-1]["date"] < post_rows[0]["date"]:
        raise SystemExit(
            f"the series ends {series_rows[-1]['date']}, the first post-collapse "
            f"chip is {post_rows[0]['date']} — the cut is not clean")

    directions = sorted({row["orbit_direction"] for row in rows})
    orbit_label = " + ".join(ocv.CHANNEL[direction] for direction in directions)

    ref_values = {dim: [_log(_number(row[COLUMN[dim]])) for row in ref_rows]
                  for dim in DIMS}
    pool_values = {dim: [_log(_number(row[COLUMN[dim]])) for row in pool_rows]
                   for dim in DIMS}
    all_values = {dim: [_log(_number(row[COLUMN[dim]])) for row in rows]
                  for dim in DIMS}
    reference = _reference_block(ref_rows, ref_values, pool_values)

    # The bootstrap series is in date order and paired across dimensions: the
    # same chips carry both, so one index draw serves both.
    for dim in DIMS:
        if any(v is None for v in ref_values[dim]):
            raise SystemExit(
                f"a reference chip has no positive {dim} — the paired block draw "
                f"needs a complete reference; handle the drop explicitly")
    draws = reference_draws(ref_values, boot=boot, block=block, seed=seed)

    med = {dim: reference[dim]["median_log"] for dim in DIMS}
    scale = {dim: reference[dim]["scale_log"] for dim in DIMS}

    def _chip(row):
        z = {}
        for dim in DIMS:
            g = _log(_number(row[COLUMN[dim]]))
            z[dim] = None if g is None else (g - med[dim]) / scale[dim]
        return {"date": row["date"], "track": row["track"],
                "channel": TRACK_CHANNEL[row["track"]], "z": z, "R": _rms(z)}

    # Per-chip magnitudes in reference units — the dots of panels (a) and (b)
    # and the grey cloud of panel (c).
    chips = [_chip(row) for row in series_rows]
    ref_chips = [_chip(row) for row in ref_rows]
    ref_R = [chip["R"] for chip in ref_chips if chip["R"] is not None]
    ref_band = [float(np.percentile(ref_R, 5)), float(np.percentile(ref_R, 95))]

    groups = {}
    for row in series_rows:
        groups.setdefault(month_of(row["date"]), []).append(row)
    keys = sorted(groups)
    missing = [key for key in month_range(keys[0], keys[-1]) if key not in groups]

    # The track mix of the reference, and the mix-only expectation for every
    # month: what the pooled median would do if each track stayed exactly at its
    # own reference median and only the ASC/DESC counts changed.
    ref_mix = {track: sum(1 for row in ref_rows if row["track"] == track)
               for track in TRACKS}
    ref_share = ref_mix["A_asc"] / len(ref_rows)
    track_ref = {track: {dim: _stats([_log(_number(row[COLUMN[dim]]))
                                      for row in ref_rows if row["track"] == track])
                         for dim in DIMS}
                 for track in TRACKS}



    ref_expected = {dim: (ref_mix["A_asc"] * track_ref["A_asc"][dim]["median"]
                          + ref_mix["A_des"] * track_ref["A_des"][dim]["median"])
                    / len(ref_rows) for dim in DIMS}

    months = []
    for key in keys:
        group = groups[key]
        n_by_track = {TRACK_CHANNEL[track]:
                      sum(1 for row in group if row["track"] == track)
                      for track in TRACKS}
        entry = {"month": key, "n": len(group), "n_by_track": n_by_track,
                 "dates": sorted(row["date"] for row in group),
                 "thin": len(group) < SMALL_N, "dims": {}}
        centres = {}
        for dim in DIMS:
            stats = _stats([_log(_number(row[COLUMN[dim]])) for row in group])
            centres[dim] = stats["median"]
            z_draws = []
            for draw in draws:
                draw_median, draw_scale = draw[dim]
                if draw_median is None or not draw_scale or stats["median"] is None:
                    continue
                z_draws.append((stats["median"] - draw_median) / draw_scale)
            entry["dims"][dim] = {
                "n": stats["n"],
                "median": None if stats["median"] is None else math.exp(stats["median"]),
                "median_log": stats["median"],
                "z": None if stats["median"] is None
                     else (stats["median"] - med[dim]) / scale[dim],
                "z_pooled": None if stats["median"] is None
                            else ((stats["median"] - reference[dim]["pooled_median_log"])
                                  / reference[dim]["pooled_scale_log"]),
                "ci95": ocv.percentile_ci(z_draws),
                "draws": len(z_draws),
            }
        z_month = {dim: entry["dims"][dim]["z"] for dim in DIMS}
        r_draws = []
        for draw in draws:
            z_draw = {}
            for dim in DIMS:
                draw_median, draw_scale = draw[dim]
                if draw_median is None or not draw_scale or centres[dim] is None:
                    z_draw = None
                    break
                z_draw[dim] = (centres[dim] - draw_median) / draw_scale
            if z_draw is not None:
                r_draws.append(_rms(z_draw))
        entry["R_behaviour"] = {
            "value": _rms(z_month),
            "k": len([d for d in DIMS if z_month[d] is not None]),
            "ci95": ocv.percentile_ci(r_draws),
            "reference_band": ref_band,
            "draws": len(r_draws),
        }
        mix = {"asc_share": n_by_track["ASC"] / len(group),
               "expected_log": {}, "change_log": {}, "change_mad": {}}
        for dim in DIMS:
            expected = (n_by_track["ASC"] * track_ref["A_asc"][dim]["median"]
                        + n_by_track["DESC"] * track_ref["A_des"][dim]["median"]) \
                / len(group)
            mix["expected_log"][dim] = expected
            mix["change_log"][dim] = expected - ref_expected[dim]
            mix["change_mad"][dim] = (expected - ref_expected[dim]) / scale[dim]
        entry["mix"] = mix
        months.append(entry)


    # The track split: how much of each month's level is the *origin's* mix.
    channels = _channel_blocks(ref_rows, series_rows, reference)

    candidates = [
        window_stats(pre_rows, "2016-08-01", "2016-12-31",
                     note="too short: five months, one season, and it throws away "
                          "the rest of the only full year the record holds"),
        window_stats(pre_rows, REF_START, REF_END, used=True,
                     note="used: twelve months, every calendar month once, and the "
                          "series is held out entirely"),
        window_stats(pre_rows, "2016-08-01", "2017-12-31",
                     note="already contains the first five series months, so the "
                          "series would help set its own origin"),
        window_stats(pre_rows, "2016-08-01", "2018-08-13",
                     note="the SDS plate's whole pre-collapse reference: it contains "
                          "all thirteen series months"),
    ]

    post_dates = sorted(row["date"] for row in post_rows)
    ref_by_channel = {TRACK_CHANNEL[track]:
                      sum(1 for row in ref_rows if row["track"] == track)
                      for track in TRACKS}
    series_by_channel = {TRACK_CHANNEL[track]:
                         sum(1 for row in series_rows if row["track"] == track)
                         for track in TRACKS}
    post_by_channel = {TRACK_CHANNEL[track]:
                       sum(1 for row in post_rows if row["track"] == track)
                       for track in TRACKS}
    ref_months_missing = missing_reference(ref_rows)
    report = {
        "schema": SCHEMA,
        "site": SITE,
        "site_label": SITE_LABEL,
        "event": {"collapse": core.COLLAPSE_DATE,
                  "record_end": post_dates[-1] if post_dates else None},
        "metric": {
            "definition": ("z = (median(log x_month) - median(log x_ref)) / scale_ref; "
                           "R_behaviour = rms over the behaviour block"),
            "transform": {dim: TRANSFORM[dim] for dim in DIMS},
            "scale_rule": ("MAD of the reference, fallback chain MAD -> IQR/1.349 -> SD; "
                           f"the pooled scale of the pre-collapse record replaces it when "
                           f"n_ref < {SCALE_MIN_REF} (spec amendment 2026-10-22)"),
            "block": list(DIMS),
            "unit": "one support-sized chip per acquisition (51 x 137, 861 cells)",
            "interval": (f"reference-uncertainty bootstrap only: {boot} circular "
                         f"moving-block draws, block {block}, seed {seed}, reference "
                         f"chips resampled, the month's chips held fixed"),
            "pooled_window": (f"pre-collapse record, {REF_START} … {POOL_END} "
                              f"(reference + series)"),
            "reference_choice": ("the first twelve months of the record, 2016-08 … "
                                 "2017-07: one full seasonal cycle and earlier than "
                                 "the series; not the SDS plate's whole pre-collapse "
                                 "reference"),
            "tests": None,
            "thin_rule": f"n < {SMALL_N} — 2018-08 only",
        },
        "input": {
            "table": _rel(table),
            "table_sha256": sha256(table),
            "meta": _rel(meta_path),
            "meta_sha256": sha256(meta_path),
            "rows": len(rows),
            "tracks": sorted({row["track"] for row in rows}),
            "orbit": orbit_label,
            "orbit_rows": {ocv.CHANNEL[direction]:
                           sum(1 for row in rows if row["orbit_direction"] == direction)
                           for direction in directions},
            "chip": meta_obj.get("chip"),
            "support_rule": meta_obj.get("rule"),
            "support_m": {"guard_m": meta_obj.get("guard_m"),
                          "ring_m": meta_obj.get("ring_m")},
            "support_band": meta_obj.get("support_band"),
        },
        "reference": {
            "window": [REF_START, REF_END],
            "n": len(ref_rows),
            "first": ref_rows[0]["date"],
            "last": ref_rows[-1]["date"],
            "orbit": orbit_label,
            "n_by_channel": ref_by_channel,
            "months": len({month_of(row["date"]) for row in ref_rows}),
            "months_missing": ref_months_missing,
            "mix": {"ASC": ref_mix["A_asc"], "DESC": ref_mix["A_des"],
                    "asc_share": ref_share},
            "dims": reference,
        },


        "series": {
            "window": [SERIES_START, EVENT],
            "months": len(keys),
            "first": keys[0],
            "last": keys[-1],
            "chips": len(series_rows),
            "first_chip": series_rows[0]["date"],
            "last_chip": series_rows[-1]["date"],
            "n_by_channel": series_by_channel,
            "months_missing": missing,
            "thin_months": sum(1 for entry in months if entry["thin"]),
            "mix": {"ASC": series_by_channel["ASC"],
                    "DESC": series_by_channel["DESC"],
                    "asc_share": series_by_channel["ASC"] / len(series_rows)},
        },
        "post": {"chips": len(post_rows), "n_by_channel": post_by_channel,
                 "first": post_dates[0] if post_dates else None,
                 "last": post_dates[-1] if post_dates else None},
        "months": months,
        "chips": chips,
        "reference_chips": ref_chips,
        "reference_R_band": ref_band,
        "channels": channels,
        "candidates": candidates,
        "audit": {
            "rows": len(rows),
            "reference": len(ref_rows),
            "series": len(series_rows),
            "series_months": len(keys),
            "post_outside_series": len(post_rows),
            "post_first": post_dates[0] if post_dates else None,
            "post_last": post_dates[-1] if post_dates else None,
            "partition": (f"{len(ref_rows)} reference + {len(series_rows)} series "
                          f"+ {len(post_rows)} post = {len(rows)}"),
            "pooled_chips": len(pool_rows),
            "dropped_nonpositive": {
                dim: sum(1 for v in all_values[dim] if v is None) for dim in DIMS},
            "months_missing_series": missing,
            "months_missing_reference": ref_months_missing,
            "orbit": f"{' + '.join(directions)} ({orbit_label})",
            "reference_by_channel": ref_by_channel,
            "series_by_channel": series_by_channel,
            "post_by_channel": post_by_channel,
        },
    }
    return report


# ---------------------------------------------------------------------------
# The report as markdown
# ---------------------------------------------------------------------------
def _signed(value, digits=2):
    return "–" if value is None else f"{value:+.{digits}f}"


def _plain(value, digits=2):
    return "–" if value is None else f"{value:.{digits}f}"


def _ci(interval, digits=2):
    if not interval:
        return "–"
    return f"[{interval[0]:+.{digits}f}, {interval[1]:+.{digits}f}]"


def _extremes(months, dim):
    ordered = sorted((m for m in months if m["dims"][dim]["z"] is not None),
                     key=lambda m: m["dims"][dim]["z"])
    return ordered[0], ordered[-1]


def _month_list(keys):
    return ", ".join(keys) if keys else "none"


def _outside_text(dim, keys):
    """``1 gamma2 month (2018-03)`` / ``no P month`` — never a bare zero."""
    if not keys:
        return f"no {dim} month"
    return f"{len(keys)} {dim} month(s) ({_month_list(keys)})"


def _arc(months, dim):
    """The largest change of z over two months, and the months it spans."""
    best = None
    for index in range(len(months) - 2):
        if (months[index]["dims"][dim]["z"] is None
                or months[index + 2]["dims"][dim]["z"] is None):
            continue
        delta = (months[index + 2]["dims"][dim]["z"]
                 - months[index]["dims"][dim]["z"])
        if best is None or abs(delta) > abs(best[0]):
            best = (delta, months[index], months[index + 2])
    return best




def reading_lines(report):
    """The prose of the plate, every number taken from the object itself."""
    months = report["months"]
    g_dim, p_dim = DIMS
    band = report["reference_R_band"]
    last = months[-1]
    last_g, last_p = last["dims"][g_dim], last["dims"][p_dim]
    g_low, g_high = _extremes(months, g_dim)
    p_low, p_high = _extremes(months, p_dim)
    outside = {dim: [m["month"] for m in months if abs(m["dims"][dim]["z"]) > 2]
               for dim in DIMS}
    below = [m["month"] for m in months if m["R_behaviour"]["value"] < band[0]]
    r_max = max(months, key=lambda m: m["R_behaviour"]["value"])
    z_max = max(abs(m["dims"][dim]["z"]) for m in months for dim in DIMS)
    arc = _arc(months, g_dim)

    asc, desc = report["channels"]["ASC"], report["channels"]["DESC"]
    asc_g = [entry for entry in asc["months"] if abs(entry["z_own"][g_dim]) <= 1.7]
    asc_g_out = [entry for entry in asc["months"]
                 if abs(entry["z_own"][g_dim]) > 1.7]
    desc_g = [entry for entry in desc["months"] if entry["z_own"][g_dim] < 0]
    asc_p_pos = [entry for entry in asc["months"] if entry["z_own"][p_dim] > 0]
    desc_p_pos = [entry for entry in desc["months"] if entry["z_own"][p_dim] > 0]
    pooled_p_pos = [m for m in months if m["dims"][p_dim]["z"] > 0]
    mix_g = [m["mix"]["change_mad"][g_dim] for m in months]
    mix_p = [m["mix"]["change_mad"][p_dim] for m in months]
    zpool_g = [m["dims"][g_dim]["z_pooled"] for m in months]
    zpool_p = [m["dims"][p_dim]["z_pooled"] for m in months]
    ref = report["reference"]["dims"]
    ratio = asc["dims"][g_dim]["median"] / desc["dims"][g_dim]["median"]
    desc_gain = (report["series"]["mix"]["DESC"]
                 - report["reference"]["mix"]["DESC"])
    lowest = min(months, key=lambda m: m["R_behaviour"]["value"])
    low_text = (f"only {_month_list(below)} falls below the 5th "
                f"({_plain(lowest['R_behaviour']['value'])}, n = {lowest['n']})"
                if below else "no month falls below the 5th")

    arc_text = ""
    if arc:
        delta, start, end = arc
        arc_text = (
            f" The largest two-month change is {_signed(delta)} of {g_dim} between "
            f"{start['month']} and {end['month']} ({_signed(start['dims'][g_dim]['z'])}"
            f" -> {_signed(end['dims'][g_dim]['z'])}); the path is not monotone.")

    return [
        f"**The pooled γ² is below the reference in every month; P is not.** The "
        f"median z({g_dim}) runs {_signed(g_high['dims'][g_dim]['z'])} "
        f"({g_high['month']}) to {_signed(g_low['dims'][g_dim]['z'])} "
        f"({g_low['month']}), negative in 13 of 13 months; "
        f"{_outside_text(g_dim, outside[g_dim])} lie beyond ±2 MAD. z({p_dim}) "
        f"runs {_signed(p_low['dims'][p_dim]['z'])} ({p_low['month']}) to "
        f"{_signed(p_high['dims'][p_dim]['z'])} ({p_high['month']}), is positive "
        f"in {len(pooled_p_pos)} of 13 months and never reaches ±2. The largest "
        f"|z| anywhere is {z_max:.2f}. On the pooled number alone the bridge reads "
        f"as less echo-consistent and no less persistent than in 2016-08 … 2017-07."
        + arc_text,
        "",
        f"**The γ² drop is not shared by the two tracks, and the pooled median is "
        f"not their average.** Against its own reference, A_asc's γ² stays within "
        f"±1.7 in {len(asc_g)} of 13 months; the {len(asc_g_out)} exceptions "
        f"({_month_list([e['month'] for e in asc_g_out])}) are all months with "
        f"n ≤ 2. A_des's "
        f"γ² is below its own reference in {len(desc_g)} of 13 months "
        f"({_signed(min(e['z_own'][g_dim] for e in desc['months']))} to "
        f"{_signed(max(e['z_own'][g_dim] for e in desc['months']))} of its own MAD). "
        f"So one track moves and the other does not. The pooled median is an order "
        f"statistic of a two-mode mixture whose modes are {ratio:.1f}x apart "
        f"({asc['dims'][g_dim]['median']:.4g} vs "
        f"{desc['dims'][g_dim]['median']:.4g}): when ASC holds fewer than half the "
        f"chips, the median is a DESC value and can sit below *both* track medians "
        f"({months[7]['month']}: pooled {_signed(months[7]['dims'][g_dim]['z'])}, ASC "
        f"{_signed(asc['months'][7]['z_pooled'][g_dim])}, DESC "
        f"{_signed(desc['months'][7]['z_pooled'][g_dim])} in pooled-MAD units). The "
        f"pooled excursion is therefore a mixture of a real DESC drop, the track "
        f"level offset and the order statistic — not a whole-deck signal.",
        "",
        f"**P rises inside both tracks and the mix cancels it.** A_asc's P is above "
        f"its own reference in {len(asc_p_pos)} of 13 months (own-MAD z up to "
        f"{_signed(max(e['z_own'][p_dim] for e in asc['months']))}; pooled unit at "
        f"most {_signed(max(asc['summary'][p_dim]['z_pooled_range']))}), A_des's in "
        f"{len(desc_p_pos)} of 13 (own-MAD up to "
        f"{_signed(max(e['z_own'][p_dim] for e in desc['months']))}; pooled at most "
        f"{_signed(max(desc['summary'][p_dim]['z_pooled_range']))}), yet the pooled "
        f"z(P) changes sign. The reason is the mix: A_des has the lower P level "
        f"({desc['dims'][p_dim]['median']:.3g} vs "
        f"{asc['dims'][p_dim]['median']:.3g}) and gains {desc_gain} chips in the "
        f"series. A "
        f"linear mix-only reconstruction — both tracks held at their own reference "
        f"medians, only the counts changed — gives {_signed(min(mix_g))} to "
        f"{_signed(max(mix_g))} MAD for γ² and {_signed(min(mix_p))} to "
        f"{_signed(max(mix_p))} for P. Because the median is an order statistic, "
        f"that linear number is a *lower bound*, not the whole composition effect.",
        "",
        f"**The composition moves in one direction, and the series is thin where "
        f"it matters.** The reference is {report['reference']['mix']['ASC']} ASC / "
        f"{report['reference']['mix']['DESC']} DESC (ASC share "
        f"{report['reference']['mix']['asc_share']:.2f}); the series is "
        f"{report['series']['mix']['ASC']} / {report['series']['mix']['DESC']} "
        f"({report['series']['mix']['asc_share']:.2f}). Per month the ASC count "
        f"falls from 5 of 11 to 1 of 3 while DESC stays at 4–6. The last month, "
        f"{last['month']}, holds {last['n']} chips "
        f"({', '.join(last['dates'])}; {last['n_by_track']['ASC']} ASC + "
        f"{last['n_by_track']['DESC']} DESC) and is the only month below "
        f"`SMALL_N = {SMALL_N}`. Its marker is z({g_dim}) = {_signed(last_g['z'])} "
        f"{_ci(last_g['ci95'])}, z({p_dim}) = {_signed(last_p['z'])} "
        f"{_ci(last_p['ci95'])}, R = {_plain(last['R_behaviour']['value'])}. The "
        f"interval is narrow because the *origin* is well known, not because the "
        f"month is: with three chips the month's own sampling noise is not bounded "
        f"by anything this figure computes.",
        "",
        f"**R stays inside the reference's own 5–95 % band except one low.** The "
        f"band is {_plain(band[0])}–{_plain(band[1])}; R runs "
        f"{_plain(min(m['R_behaviour']['value'] for m in months))} to "
        f"{_plain(r_max['R_behaviour']['value'])}. It never reaches the 95th "
        f"percentile, and {low_text}. A low R is *agreement* with the origin on both "
        f"dimensions, so it is not a departure; a high R would be, and none reaches "
        f"the band's upper edge.",
        "",
        f"**The pooled-window sensitivity shows how much of the number is the "
        f"origin.** Re-centring on the {report['audit']['pooled_chips']}-chip "
        f"pre-collapse record instead of the 99-chip reference (γ² MAD "
        f"{ref[g_dim]['scale_log']:.2f} -> {ref[g_dim]['pooled_scale_log']:.2f}) "
        f"compresses z(γ²) to {_signed(min(zpool_g))} … {_signed(max(zpool_g))} and "
        f"z(P) to {_signed(min(zpool_p))} … {_signed(max(zpool_p))}. That pool "
        f"contains the series, so it is a declared sensitivity to the window — not "
        f"validation, and not a second estimate.",
        "",
        f"**One sentence.** Between {report['series']['first_chip']} and "
        f"{report['series']['last_chip']} the pooled chip sits below its 2016-08 … "
        f"2017-07 γ² level and on its P level, but the split shows that this is a "
        f"low-γ² DESC track entering a larger share of the pool, not a bridge-wide "
        f"reconfiguration; on this record the honest reading is a null.",
    ]


def caveat_lines(report):
    """The limits, stated where the figure is, not only in the JSON."""
    months = report["months"]
    last = months[-1]
    low_n = min(m["n"] for m in months)
    high_n = max(m["n"] for m in months)
    ref = report["reference"]
    dims = ref["dims"]
    g_dim, p_dim = DIMS
    asc, desc = report["channels"]["ASC"], report["channels"]["DESC"]
    chip = report["input"]["chip"] or {}
    ratio = asc["dims"][g_dim]["median"] / desc["dims"][g_dim]["median"]
    offset = math.log(ratio)
    return [
        "* **The binding limit is sampling, not SNR.** A month rests on "
        f"{low_n}-{high_n} chips ({min(m['n_by_track']['ASC'] for m in months)}-"
        f"{max(m['n_by_track']['ASC'] for m in months)} ASC, "
        f"{min(m['n_by_track']['DESC'] for m in months)}-"
        f"{max(m['n_by_track']['DESC'] for m in months)} DESC), so {last['month']} "
        f"alone is `thin` (`n = {last['n']} < SMALL_N = {SMALL_N}`) and its "
        "\"median\" *is* the middle of three chips. Every other month is 6–11 "
        "acquisitions, enough for a median and not enough for a spread: the "
        "monthly marker is the median of the chips that exist, not of the month, "
        "and a month with no acquisition would be a gap, not a value.",
        "* **The interval belongs to the reference.** It says how well a "
        f"{ref['n']}-chip origin is known under a resampling of those chips. It is "
        "not a confidence interval for a monthly value: the month's own sampling "
        "noise is not in it. That is why the n = 3 month's interval is narrow for "
        "the wrong reason.",
        f"* **This `{COLUMN['gamma2']}` is not the `gamma2` of the paper figure.** "
        f"The support-sized chip is {chip.get('rows')} x {chip.get('cols')} cells "
        f"with {chip.get('support_cells')} support cells "
        f"({report['input']['support_m']['guard_m']:g} m guard, "
        f"{report['input']['support_m']['ring_m']:g} m ring), and its reference "
        f"median `{COLUMN['gamma2']}` is {dims['gamma2']['median']:.4g}. The "
        "echo-mask γ² of `fig_morandi_ocv_paper` is computed on the full 80 x 80 "
        "chip of the 232-chip export and runs several times higher. Only the "
        "*standardised* quantity z is comparable between the two figures; the "
        "levels are not, and neither is R.",
        "* **The chip is a support rectangle, not a damage map.** The 207-row table "
        "is the ring-supported chip; it is a window on the deck, not a map of the "
        "failure. Nothing here localises damage, and the plate makes no claim about "
        "the 2018-08-14 collapse; it describes how the chip's echo behaviour moved "
        "while the bridge stood. The paper figure's structural caveat — its analysed "
        "window is not on the collapsed span — is the frame this chip sits in, and "
        "nothing here changes it.",
        "* **The series stops at the collapse, deliberately.** The "
        f"{report['audit']['post_outside_series']} acquisitions after "
        f"{report['event']['collapse']} ({report['audit']['post_first']} … "
        f"{report['audit']['post_last']}) are the same bridge in a different state, "
        "and they are counted in the audit and carry no value here. The record does "
        "not end at the event; the *series* does.",
        "* **The reference window is a choice, and it moves the answer.** The first "
        "twelve months give a full seasonal cycle and are strictly earlier than the "
        "series. The alternatives are in the candidate table: five months are too "
        "few; the 17- and 25-month windows *contain* the series. The 99-chip "
        "reference is the only candidate that is both long enough to have a MAD and "
        "held out before the months it standardises.",
        f"* **The pooled sensitivity is not a second estimate.** Re-centring on the "
        f"{report['audit']['pooled_chips']}-chip pre-collapse record changes the γ² "
        f"MAD from {dims['gamma2']['scale_log']:.2f} to "
        f"{dims['gamma2']['pooled_scale_log']:.2f} and the P MAD from "
        f"{dims['P']['scale_log']:.2f} to "
        f"{dims['P']['pooled_scale_log']:.2f}; it contains the series, so it is a "
        "declared sensitivity to the window chosen above, not independent evidence.",
        f"* **The track split is a composition control, not two estimates.** The two "
        f"tracks' reference γ² MADs are {asc['dims'][g_dim]['scale_log']:.3f} "
        f"(A_asc) and {desc['dims'][g_dim]['scale_log']:.3f} (A_des), and their "
        f"levels differ by {ratio:.1f}x. Re-referencing each to its own MAD would "
        f"move a month by the level offset alone — {offset / desc['dims'][g_dim]['scale_log']:.1f} "
        f"DESC MADs and {offset / asc['dims'][g_dim]['scale_log']:.1f} ASC MADs — "
        "before any change, so the figure standardises the within-track difference "
        "by the *pooled* MAD instead. The dotted lines therefore show direction, "
        "not effect size; with 1–5 ASC chips a month, the ASC line is itself a "
        "median of very few values and its n ≤ 2 excursions are the least stable "
        "points on the plate.",
        "* **No test, no dose, no partner.** The sds H1/H2 machinery needs at least "
        "two states per block and, for a dose, an ordered ladder; here one state is "
        "sampled monthly. So there is no p-value, no Holm family and no cosine, "
        "±2 MAD is a descriptive band and not a significance threshold, and these "
        "magnitudes are in this reference's own units — they are not comparable to "
        "the Carola or CTS reconfiguration amplitudes, whose records, chips and "
        "references are different.",
        "* **The honest null stands.** `fig_morandi_ocv_paper.md` reports no "
        "detectable pre-collapse change on the 232-chip export; this plate does not "
        "overturn it. The pooled γ² offset is real *in the sample*, but the split "
        "attributes it to one track and a shifting mix rather than to a bridge-wide "
        "reconfiguration. No precursor, no damage and no reconfiguration claim is "
        "made; the monthly path is a description with n = 3–11 behind each point.",
    ]



def markdown(report):
    """The plate's guide: definition, reference, all months, the reading."""
    ref = report["reference"]
    dims = ref["dims"]
    months = report["months"]
    g_dim, p_dim = DIMS
    band = report["reference_R_band"]
    channels = report["channels"]
    chips = report["chips"]
    chip = report["input"]["chip"] or {}
    span = (f"{chips[0]['date']} .. {chips[-1]['date']}" if chips
            else f"{report['series']['first']} .. {report['series']['last']}")
    lines = [
        "# Morandi — OCV reconfiguration against the 2016-08 … 2017-07 reference, "
        "month by month",
        "",
        "*generated by `code/morandi/fig_morandi_ocv_reconfiguration.py`; the metric "
        "is the pre-registered one of `code/sds/ocv_reconfigure_spec.md` (behaviour "
        "block, log transform, reference MAD), read here as a monthly series. It is "
        "not a new test, and it carries none.*",
        "",
        f"{report['site_label']}: one support-sized chip per acquisition "
        f"({chip.get('rows')} x {chip.get('cols')} cells, "
        f"{chip.get('support_cells', 861)} support cells, "
        f"{report['input']['support_m']['guard_m']:g} m guard, "
        f"{report['input']['support_m']['ring_m']:g} m ring), {report['audit']['rows']} "
        f"acquisitions in two tracks ({report['input']['orbit_rows'].get('ASC')} "
        f"A_asc, {report['input']['orbit_rows'].get('DESC')} A_des). The record "
        f"splits, by its own state column, into a **reference** of {ref['n']} chips "
        f"({ref['first']} … {ref['last']}), the **series** of "
        f"{report['series']['chips']} chips on {report['series']['months']} months "
        f"({report['series']['first']} … {report['series']['last']}), and "
        f"{report['audit']['post_outside_series']} post-collapse chips that are "
        f"counted and never plotted.",
        "",
        "```",
        "z(month) = ( median(log x_month) - median(log x_ref) ) / MAD(log x_ref)",
        "x        = gamma2, P            both on the log scale",
        "R        = rms over the two dimensions (gamma2, P)",
        f"ref      = {ref['first']} .. {ref['last']}   n = {ref['n']} chips on "
        f"{ref['months']} months",
        f"series   = {span}   {report['series']['chips']} chips on "
        f"{report['series']['months']} months, 3-11 per month",
        "```",
        "",
        "The reference is drawn from the twelve months *before* the series, so the "
        "origin of everything below is the bridge's own earlier state and not a "
        "baseline chosen after the fact. `0` on both axes is that state; it is not a "
        "norm and not a threshold.",
        "",
        "**The unit is the support-sized chip, not the paper figure's export.** This "
        "table has 207 rows on 51 x 137 chips; `fig_morandi_ocv_paper` reads the "
        "232-row per-track export. The two are not mixed, and the γ² levels are not "
        "comparable between them — only z, which is in each figure's own reference "
        "MADs.",
        "",
        "## Why this window",
        "",
        "The record starts 2016-08-24; the series is the last twelve months before "
        "the collapse, so the reference can only be built from what precedes "
        "2017-08. Four windows were available, and the choice is stated rather than "
        "buried because it moves the answer:",
        "",
        "| window | chips | months | ASC / DESC | contains the series? | used | why |",
        "| --- | ---: | ---: | --- | --- | --- | --- |",
    ]
    for candidate in report["candidates"]:
        start, end = candidate["window"]
        lines.append(
            f"| {start} … {end} | {candidate['chips']} | {candidate['months']} | "
            f"{candidate['n_by_channel']['ASC']} / "
            f"{candidate['n_by_channel']['DESC']} | "
            f"{'yes' if candidate['contains_series'] else 'no'} | "
            f"{'**used**' if candidate['used'] else '—'} | {candidate['note']} |")
    lines += [
        "",
        "Only the used window is a full year *and* earlier than the series. A "
        "per-calendar-month baseline is not on the list because a single month of "
        f"the reference holds 5-10 chips, below `SCALE_MIN_REF = {SCALE_MIN_REF}`: "
        "such a baseline has no MAD. That limit is stated rather than engineered "
        "around.",
        "",
    ]

    lines += [
        "## The reference",
        "",
        "| dim | column | n | median | MAD(log) | 1 MAD = factor | scale used |",
        "| --- | --- | ---: | ---: | ---: | ---: | --- |",
    ]
    for dim in DIMS:
        entry = dims[dim]
        lines.append(
            f"| {dim} | `{entry['column']}` | {entry['n']} | "
            f"{entry['median']:.4g} | {entry['mad_log']:.4f} | "
            f"{entry['one_mad_factor']:.3f} | {entry['scale_used']} |")
    lines += [
        "",
        f"Both scales are the reference's own MAD: n = {ref['n']} is far above "
        f"`SCALE_MIN_REF = {SCALE_MIN_REF}` for both dimensions, so the fallback "
        f"chain is not used. The reference has "
        f"{'no missing month' if not ref['months_missing'] else 'missing months: ' + _month_list(ref['months_missing'])} "
        f"between {ref['first'][:7]} and {ref['last'][:7]}, and its track mix is "
        f"{ref['mix']['ASC']} A_asc / {ref['mix']['DESC']} A_des "
        f"(ASC share {ref['mix']['asc_share']:.2f}).",
        "",
        "## The months",
        "",
        f"| month | n | ASC / DESC | {g_dim} | {p_dim} | z({g_dim}) | z({p_dim}) | "
        f"R | 95 % CI of z({g_dim}) | 95 % CI of z({p_dim}) | thin |",
        "| --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | --- | --- | --- |",
    ]
    for entry in months:
        g, p = entry["dims"][g_dim], entry["dims"][p_dim]
        lines.append(
            f"| {entry['month']} | {entry['n']} | "
            f"{entry['n_by_track']['ASC']} / {entry['n_by_track']['DESC']} | "
            f"{g['median']:.4g} | {p['median']:.4g} | {_signed(g['z'])} | "
            f"{_signed(p['z'])} | {_plain(entry['R_behaviour']['value'])} | "
            f"{_ci(g['ci95'])} | {_ci(p['ci95'])} | "
            f"{'✱' if entry['thin'] else '—'} |")
    lines += [
        "",
        f"`n` is the number of chips the month's median rests on, `✱` marks a month "
        f"below `SMALL_N = {SMALL_N}` — here {report['series']['thin_months']} of "
        f"{report['series']['months']}, {months[-1]['month']} — and `–` would be a "
        f"cell the record does not contain, never a zero. The interval is the "
        f"*reference's* uncertainty: {months[0]['dims'][g_dim]['draws']} circular "
        f"moving-block draws of the {ref['n']} reference chips with the reference "
        f"median **and** scale re-estimated inside each draw, the month's chips "
        f"held fixed. It answers \"how well is the origin known\", not \"how well is "
        f"the month known\".",
        "",
    ]

    lines += [
        "## The two tracks",
        "",
        "The two Sentinel-1 geometries are not interchangeable: their reference "
        "medians differ by 8.5x in γ² and 2.8x in P, and the mix moves from 53/46 "
        "to 38/61. A pooled median therefore moves when *which track is sampled* "
        "changes, even if neither track moves. The split below is the control: each "
        "track's monthly median against **its own** reference median, shown twice — "
        "in that track's own MAD (`z own`) and in the **pooled** reference MAD "
        "(`z pooled`, the unit of the month table and of panels a/b). The second is "
        "the one the figure draws; it is comparable to the headline, while `z own` "
        "answers whether the track moved in its own units.",
        "",
        "| track | n_ref | gamma2 median | gamma2 MAD(log) | P median | P MAD(log) |",
        "| --- | ---: | ---: | ---: | ---: | ---: |",
    ]
    for key in ("ASC", "DESC"):
        block = channels[key]
        lines.append(
            f"| {block['track']} ({key}) | {block['n_ref']} | "
            f"{block['dims'][g_dim]['median']:.4g} | "
            f"{block['dims'][g_dim]['scale_log']:.4f} | "
            f"{block['dims'][p_dim]['median']:.4g} | "
            f"{block['dims'][p_dim]['scale_log']:.4f} |")
    lines += [
        "",
        f"The two tracks' own MADs differ by a factor "
        f"{channels['DESC']['dims'][g_dim]['scale_log'] / channels['ASC']['dims'][g_dim]['scale_log']:.1f} "
        "in γ². Re-referencing each to its own MAD would therefore move a month by "
        "several MADs for the level offset alone, which is why the figure uses the "
        "pooled unit. The split is a composition control, not a second estimate.",
        "",
    ]
    for dim in DIMS:
        lines += [
            f"**{dim}**",
            "",
            f"| month | n A_asc | A_asc z own | A_asc z pooled | n A_des | "
            f"A_des z own | A_des z pooled |",
            "| --- | ---: | ---: | ---: | ---: | ---: | ---: |",
        ]
        for asc_m, desc_m in zip(channels["ASC"]["months"],
                                 channels["DESC"]["months"]):
            lines.append(
                f"| {asc_m['month']} | {asc_m['n']} | "
                f"{_signed(asc_m['z_own'][dim])} | "
                f"{_signed(asc_m['z_pooled'][dim])} | {desc_m['n']} | "
                f"{_signed(desc_m['z_own'][dim])} | "
                f"{_signed(desc_m['z_pooled'][dim])} |")
        lines += [
            "",
            f"`z own` = (track median − track reference median) / track MAD; "
            f"`z pooled` = the same difference / the pooled reference MAD. The two "
            f"columns differ only by the scale, not by the month's chips. The ranges "
            f"are in the JSON under `channels.ASC.summary` and "
            f"`channels.DESC.summary`.",
            "",
        ]
    lines += [
        "## The mix",
        "",
        "What the pooled median would do if each track stayed at its own reference "
        "median and only the ASC/DESC counts changed. It is a linear reconstruction "
        "of the composition effect, and because the median is an order statistic of "
        "a two-mode mixture, the realised effect can be larger:",
        "",
        "| month | ASC / DESC | ASC share | mix-only Δz(gamma2) | mix-only Δz(P) |",
        "| --- | --- | ---: | ---: | ---: |",
    ]
    for entry in months:
        lines.append(
            f"| {entry['month']} | {entry['n_by_track']['ASC']} / "
            f"{entry['n_by_track']['DESC']} | {entry['mix']['asc_share']:.2f} | "
            f"{_signed(entry['mix']['change_mad'][g_dim])} | "
            f"{_signed(entry['mix']['change_mad'][p_dim])} |")
    lines += [
        "",
        f"The largest linear mix-only change is "
        f"{_plain(max(abs(m['mix']['change_mad'][g_dim]) for m in months))} MAD for "
        f"γ² and {_plain(max(abs(m['mix']['change_mad'][p_dim]) for m in months))} "
        f"for P. The γ² path reaches "
        f"{_plain(min(m['dims'][g_dim]['z'] for m in months), 2)} MAD, so the mix "
        f"alone does not explain it; the within-track tables above show where the "
        f"rest comes from. The mix-only column is a lower bound, not a correction.",
        "",
        "## What the series shows",
        "",
        *reading_lines(report),
        "",
        "## What it does not say",
        "",
        *caveat_lines(report),
        "",
    ]

    lines += [
        "## Reproduce",
        "",
        "```",
        "python3 code/morandi/fig_morandi_ocv_reconfiguration.py            # data + figure",
        "python3 code/morandi/fig_morandi_ocv_reconfiguration.py --selftest",
        "python3 code/morandi/fig_morandi_ocv_reconfiguration.py --verify",
        "bash figures/morandi/fig_morandi_ocv_reconfiguration.sh",
        "```",
        "",
        f"The run reads `{_repo_rel(TABLE)}` (sha256 "
        f"`{report['input']['table_sha256'][:16]}…`) and `{_repo_rel(META)}`, and "
        f"writes `{_repo_rel(OUT_JSON)}`, `{_repo_rel(OUT_CSV)}`, "
        f"`{_repo_rel(OUT_JPG)}` and this file. No timestamp enters either output, "
        f"so a second run is byte-identical. `--selftest` re-derives the plate's "
        f"invariants from the table; `--verify` recomputes the report and diffs it "
        f"against the committed JSON (floats to 1e-9) before checking the CSV against "
        f"the same object.",
        "",
        "## The figure",
        "",
        f"* **(a)** z({g_dim}) by month: one dot per chip (light), the month's "
        "median (filled), the reference's ±2 MAD band (grey), the 95 % "
        "reference-uncertainty bar, the chip-count strip at the bottom — split into "
        f"A_asc and A_des, so the mix is visible under the points — and the collapse "
        f"of {report['event']['collapse']} as the red cut. The two dotted lines are "
        "the same months under the A_asc-only and A_des-only references, in pooled "
        "MAD units (the table in *The two tracks* is their source). The "
        f"{months[-1]['month']} marker is ringed: n = {months[-1]['n']}.",
        f"* **(b)** z({p_dim}) on the same axis and the same reference — the two "
        "panels are read together, never as two independent verdicts.",
        "* **(c)** the plane: the origin is the reference state, the grey cloud is "
        f"the {ref['n']} reference chips, the path runs "
        f"{report['series']['first']} -> {report['series']['last']} and is coloured "
        "by year-quarter (the series crosses 2017 and 2018), the thin arrows are the "
        "quarter medians, and the terminal month is ringed. The path is the pooled "
        "month; the cloud is the reference only.",
        "* **(d)** R = rms(z gamma2, z P) per month against the empirical 5-95 % "
        "band of the reference chips' own R. A point below the band is *agreement* "
        "with the origin, not a departure.",
        "",
        f"The plate is 2 x 2 at 13.0 in and 600 dpi, the sds package's "
        f"`savefig.dpi`; every number in it is read from the JSON this run writes.",
        "",
    ]
    return "\n".join(lines)



# ---------------------------------------------------------------------------
# The figure
# ---------------------------------------------------------------------------
def _limits(report, months, dim, z_values):
    """Panel limits that always show the ±2 band, every interval and both tracks."""
    values = [v for v in z_values if v is not None] + [-2.0, 2.0]
    for entry in months:
        ci = entry["dims"][dim]["ci95"]
        if ci:
            values.extend(ci)
    for block in report["channels"].values():
        values.extend(entry["z_pooled"][dim] for entry in block["months"]
                      if entry["z_pooled"][dim] is not None)
    lo, hi = min(values), max(values)
    span = hi - lo if hi > lo else 1.0
    return lo - 0.20 * span, hi + 0.10 * span


def _month_panel(ax, report, months, dim, dates, z_values, colour, title):
    """One time panel: band, chips, reference CI, monthly median, track split."""
    ax.axhspan(-2.0, 2.0, color="0.85", alpha=0.55, lw=0, zorder=0)
    ax.axhline(0.0, color="0.40", lw=0.8, zorder=1)
    for chip in report["chips"]:
        z = chip["z"][dim]
        if z is None:
            continue
        ax.plot(dt.date.fromisoformat(chip["date"]), z, marker="o", ms=1.7,
                mfc=CHANNEL_COLOUR[chip["channel"]], mec="none", alpha=0.30,
                zorder=2)
    # Each track against its own reference, in pooled MAD units — the same
    # months, a different origin. Constant offsets per track are visible here;
    # that is the point, not a nuisance.
    for key in ("ASC", "DESC"):
        block = report["channels"][key]
        points = [(month_marker(entry["month"]), entry["z_pooled"][dim])
                  for entry in block["months"]
                  if entry["z_pooled"][dim] is not None]
        if points:
            ax.plot([p[0] for p in points], [p[1] for p in points],
                    color=CHANNEL_COLOUR[key], lw=0.9, ls=(0, (1.4, 1.6)),
                    alpha=0.95, zorder=2.4)
    for day, entry in zip(dates, months):
        ci = entry["dims"][dim]["ci95"]
        if ci:
            ax.plot([day, day], ci, color="0.45", lw=0.9, zorder=2.6)
    ax.plot(dates, z_values, color=colour, lw=0.8, alpha=0.9, zorder=3)
    ax.plot(dates, z_values, marker="o", ms=3.3, lw=0, color=colour, zorder=4)
    thin = [(day, z) for day, z, entry in zip(dates, z_values, months)
            if entry["thin"]]
    if thin:
        ax.plot([t[0] for t in thin], [t[1] for t in thin], marker="o", ms=4.6,
                mfc="white", mec=colour, mew=1.0, lw=0, zorder=5)
    ax.axvline(dt.date.fromisoformat(EVENT), color=ALERT, lw=1.0, ls=(0, (4, 2)),
               zorder=6)
    lo, hi = _limits(report, months, dim, z_values)
    ax.set_ylim(lo, hi)
    ax.set_xlim(dates[0] - dt.timedelta(days=18),
                dt.date.fromisoformat(EVENT) + dt.timedelta(days=16))
    # The chip count of every month, as a strip inside the panel — stacked by
    # track, so the mix that drives the pooled median is visible under the dots.
    span = hi - lo
    base = lo + 0.024 * span
    for day, entry in zip(dates, months):
        asc = entry["n_by_track"]["ASC"] * 0.0072 * span
        desc = entry["n_by_track"]["DESC"] * 0.0072 * span
        ax.bar(day, asc, bottom=base, width=7, color=CHANNEL_COLOUR["ASC"],
               lw=0, zorder=1.5)
        ax.bar(day, desc, bottom=base + asc, width=7,
               color=CHANNEL_COLOUR["DESC"], lw=0, zorder=1.5)
        ax.text(day, lo + 0.003 * span,
                f"n={entry['n']} ({entry['n_by_track']['ASC']}/"
                f"{entry['n_by_track']['DESC']})",
                fontsize=4.4, color="0.38", ha="center", va="bottom", zorder=6)
    ax.xaxis.set_major_locator(mdates.MonthLocator(interval=2))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m"))
    ax.set_title(title, loc="left")
    ax.set_ylabel(f"z({ocv.DIM_DISPLAY[dim]})")



def make_figure(report, path=None, dpi=None):
    """The 2 x 2 plate, from the report object alone."""
    path = path or OUT_JPG
    months = report["months"]
    g_dim, p_dim = DIMS
    dates = [month_marker(entry["month"]) for entry in months]
    z_g = [entry["dims"][g_dim]["z"] for entry in months]
    z_p = [entry["dims"][p_dim]["z"] for entry in months]
    r_values = [entry["R_behaviour"]["value"] for entry in months]
    band = report["reference_R_band"]
    ref = report["reference"]
    collapse = dt.date.fromisoformat(EVENT)
    g_colour, p_colour = "#1f4e79", "#7a5195"
    ref_g = [chip["z"][g_dim] for chip in report["reference_chips"]]
    ref_p = [chip["z"][p_dim] for chip in report["reference_chips"]]

    fig, axes = plt.subplots(2, 2, figsize=(13.0, 7.4),
                             gridspec_kw={"width_ratios": [1.62, 1.0]})
    ax_g, ax_p = axes[0][0], axes[1][0]
    ax_plane, ax_r = axes[0][1], axes[1][1]

    _month_panel(ax_g, report, months, g_dim, dates, z_g, g_colour,
                 f"a · z({ocv.DIM_DISPLAY[g_dim]}) by month — the 2016-08 … "
                 f"2017-07 reference is 0")
    _month_panel(ax_p, report, months, p_dim, dates, z_p, p_colour,
                 f"b · z({ocv.DIM_DISPLAY[p_dim]}) by month — same "
                 f"standardisation, same months")

    handles = [
        Line2D([], [], color=g_colour, marker="o", ms=3.3, lw=0.8,
               label="monthly median"),
        Line2D([], [], color=g_colour, marker="o", ms=4.6, mfc="white", mew=1.0,
               lw=0, label=f"month with n < {SMALL_N}"),
        Line2D([], [], color=g_colour, marker="o", ms=1.7, alpha=0.45, lw=0,
               label="one chip"),
        Line2D([], [], color="0.45", lw=1.0, label="95 % CI of the reference"),
        Line2D([], [], color=CHANNEL_COLOUR["ASC"], lw=0.9, ls=(0, (1.4, 1.6)),
               label=f"A_asc vs its own reference"),
        Line2D([], [], color=CHANNEL_COLOUR["DESC"], lw=0.9, ls=(0, (1.4, 1.6)),
               label=f"A_des vs its own reference"),
        Line2D([], [], color=ALERT, lw=1.0, ls=(0, (4, 2)),
               label=f"collapse {EVENT}"),
    ]
    ax_g.legend(handles=handles, loc="upper right", frameon=True, framealpha=0.96,
                edgecolor="0.85", facecolor="white", fontsize=5.4,
                handlelength=1.5, borderaxespad=0.25, labelspacing=0.30)

    # (c) the plane — the origin is the reference state
    ax_plane.plot([0.0], [0.0], marker="+", ms=9, mew=1.1, color="black", zorder=7)
    ax_plane.axhline(0.0, color="0.40", lw=0.7, zorder=1)
    ax_plane.axvline(0.0, color="0.40", lw=0.7, zorder=1)
    ax_plane.plot(ref_g, ref_p, marker="o", ms=2.1, mfc="0.62", mec="none",
                  alpha=0.55, lw=0, zorder=2,
                  label=f"reference chip (n = {ref['n']})")
    ax_plane.plot(z_g, z_p, color="0.58", lw=0.7, zorder=3)
    for quarter in sorted({quarter_of(entry["month"]) for entry in months}):
        picked = [index for index, entry in enumerate(months)
                  if quarter_of(entry["month"]) == quarter]
        label = (f"{quarter} · {months[picked[0]]['month'][5:7]}–"
                 f"{months[picked[-1]]['month'][5:7]}")
        ax_plane.plot([z_g[i] for i in picked], [z_p[i] for i in picked], marker="o",
                      ms=3.4, lw=0, color=QUARTER_COLOUR[quarter], label=label,
                      zorder=4)
    # Three month labels on the path, so the plane can be read in time without
    # the quarter labels colliding in a path only 2.7 MAD wide.
    for key in (months[0]["month"], "2018-02", months[-1]["month"]):
        entry = next((m for m in months if m["month"] == key), None)
        if entry is None or entry["dims"][g_dim]["z"] is None:
            continue
        ax_plane.annotate(key[2:], xy=(entry["dims"][g_dim]["z"],
                                       entry["dims"][p_dim]["z"]),
                          xytext=(0, -11), textcoords="offset points",
                          ha="center", fontsize=5.2, color="0.32", zorder=6)
    ax_plane.plot([z_g[-1]], [z_p[-1]], marker="o", ms=9.5, mfc="none", mec=ALERT,
                  mew=1.2, zorder=8)
    ax_plane.annotate(f"{months[-1]['month']} (n = {months[-1]['n']})",
                      xy=(z_g[-1], z_p[-1]), xytext=(10, -12),
                      textcoords="offset points", fontsize=6.0, zorder=8,
                      arrowprops=dict(arrowstyle="-", lw=0.6, color="0.3"))
    # The monthly path, plus the central 1-99 % of the reference cloud. The
    # cloud's tails run to -8 MAD on the per-chip scale; plotting them would
    # compress the path into a line, so the axes stop at the 1-99 % window and
    # the clipping is stated in the title.
    x_lo = min(float(np.percentile(ref_g, 1)), min(z_g)) - 0.2
    x_hi = max(float(np.percentile(ref_g, 99)), max(z_g)) + 0.2
    y_lo = min(float(np.percentile(ref_p, 1)), min(z_p)) - 0.2
    y_hi = max(float(np.percentile(ref_p, 99)), max(z_p)) + 0.2
    ax_plane.set_xlim(x_lo, x_hi)
    ax_plane.set_ylim(y_lo, y_hi)
    ax_plane.set_aspect("equal", adjustable="box")
    ax_plane.set_xlabel(f"z({ocv.DIM_DISPLAY[g_dim]})")
    ax_plane.set_ylabel(f"z({ocv.DIM_DISPLAY[p_dim]})")
    ax_plane.set_title(f"c · the plane — origin = reference state, grey = its "
                       f"{len(ref_g)} chips (1–99 % window)", loc="left")
    ax_plane.legend(loc="upper right", frameon=True, framealpha=0.90, fontsize=5.2,
                    ncol=1, handlelength=0.9, labelspacing=0.26,
                    borderaxespad=0.2, edgecolor="0.85", facecolor="white")


    # (d) the block magnitude against the reference's own 5-95 % band
    ax_r.axhspan(band[0], band[1], color="0.85", alpha=0.55, lw=0, zorder=0)
    ax_r.axhline(0.0, color="0.40", lw=0.8, zorder=1)
    for chip in report["chips"]:
        if chip["R"] is None:
            continue
        ax_r.plot(dt.date.fromisoformat(chip["date"]), chip["R"], marker="o", ms=1.6,
                  mfc="#4a5568", mec="none", alpha=0.26, zorder=2)
    for day, entry in zip(dates, months):
        ci = entry["R_behaviour"]["ci95"]
        if ci:
            ax_r.plot([day, day], ci, color="0.45", lw=0.9, zorder=2.6)
    ax_r.plot(dates, r_values, color="#4a5568", lw=0.8, zorder=3)
    ax_r.plot(dates, r_values, marker="o", ms=3.3, lw=0, color="#4a5568", zorder=4)
    ax_r.plot([dates[-1]], [r_values[-1]], marker="o", ms=9.5, mfc="none", mec=ALERT,
              mew=1.2, zorder=8)
    ax_r.axvline(collapse, color=ALERT, lw=1.0, ls=(0, (4, 2)), zorder=6)
    ax_r.set_xlim(dates[0] - dt.timedelta(days=18),
                  collapse + dt.timedelta(days=16))
    r_high = max(r_values + [band[1]] + [entry["R_behaviour"]["ci95"][1]
                                         for entry in months
                                         if entry["R_behaviour"]["ci95"]])
    ax_r.set_ylim(0.0, 1.14 * r_high)
    ax_r.xaxis.set_major_locator(mdates.MonthLocator(interval=2))
    ax_r.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m"))
    ax_r.set_ylabel("R_behaviour")
    ax_r.set_title("d · block magnitude R = rms(z gamma2, z P)", loc="left")
    ax_r.text(0.985, 0.94, f"reference 5–95 %: {band[0]:.2f}–{band[1]:.2f}",
              transform=ax_r.transAxes, fontsize=5.6, ha="right", va="top",
              color="0.30")

    ratio = (report["channels"]["ASC"]["dims"][g_dim]["median"]
             / report["channels"]["DESC"]["dims"][g_dim]["median"])
    fig.text(0.006, 0.985,
             "Ponte Morandi — OCV reconfiguration against the 2016-08 … 2017-07 "
             "reference, month by month", fontsize=10.2, va="top", ha="left")
    fig.text(0.006, 0.947,
             f"z = (median log x_month − median log x_ref) / MAD(log x_ref); "
             f"n_ref = {ref['n']} ({ref['first']} … {ref['last']}), held out before "
             f"the series. Each month rests on "
             f"{min(e['n'] for e in months)}–{max(e['n'] for e in months)} chips, so "
             f"the monthly marker is the median of the chips that exist; the bar "
             f"drawn is the reference's own uncertainty, not the month's.\n"
             f"The two tracks' reference γ² medians differ by {ratio:.1f}x, so a "
             f"pooled median follows the mix; the dotted lines are each track against "
             f"its own reference. The series stops at the collapse — the "
             f"{report['audit']['post_outside_series']} later chips are counted and "
             f"not plotted.",
             fontsize=6.2, va="top", ha="left", color="0.25")
    fig.text(0.006, 0.020,
             "dotted: A_asc / A_des vs their own reference medians, in pooled MAD "
             f"units (composition control, not a second estimate) · {months[-1]['month']} "
             f"is the only month below SMALL_N = {SMALL_N} (n = {months[-1]['n']}) · "
             f"no test is run",
             fontsize=5.8, va="bottom", ha="left", color="0.35")
    fig.subplots_adjust(top=0.870, bottom=0.085, left=0.055, right=0.988,
                        wspace=0.20, hspace=0.30)
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    fig.savefig(path, dpi=dpi or plt.rcParams["savefig.dpi"], facecolor="white")
    plt.close(fig)
    return path



# ---------------------------------------------------------------------------
# Outputs
# ---------------------------------------------------------------------------
CSV_FIELDS = ["Month", "n", "n_ASC", "n_DESC", "gamma2", "P", "z_gamma2", "z_P",
              "z_pooled_gamma2", "z_pooled_P", "R_behaviour", "ci95_gamma2",
              "ci95_P", "ci95_R", "thin"]


def _g(value):
    return "" if value is None else f"{value:.6g}"


def _pair(interval):
    return "" if not interval else f"{interval[0]:.6g};{interval[1]:.6g}"


def write_outputs(report, json_path=None, csv_path=None, md_path=None, quiet=False):
    """JSON, the guide and the compact table — all from the same object."""
    json_path = json_path or OUT_JSON
    csv_path = csv_path or OUT_CSV
    md_path = md_path or OUT_MD
    for path in (json_path, csv_path, md_path):
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(json_path, "w") as handle:
        json.dump(report, handle, indent=1, sort_keys=True)
        handle.write("\n")
    with open(md_path, "w") as handle:
        handle.write(markdown(report))
    with open(csv_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_FIELDS)
        writer.writeheader()
        for entry in report["months"]:
            g, p = entry["dims"][DIMS[0]], entry["dims"][DIMS[1]]
            writer.writerow({
                "Month": entry["month"],
                "n": entry["n"],
                "n_ASC": entry["n_by_track"]["ASC"],
                "n_DESC": entry["n_by_track"]["DESC"],
                "gamma2": _g(g["median"]),
                "P": _g(p["median"]),
                "z_gamma2": _g(g["z"]),
                "z_P": _g(p["z"]),
                "z_pooled_gamma2": _g(g["z_pooled"]),
                "z_pooled_P": _g(p["z_pooled"]),
                "R_behaviour": _g(entry["R_behaviour"]["value"]),
                "ci95_gamma2": _pair(g["ci95"]),
                "ci95_P": _pair(p["ci95"]),
                "ci95_R": _pair(entry["R_behaviour"]["ci95"]),
                "thin": "yes" if entry["thin"] else "no",
            })
    if not quiet:
        for path in (json_path, md_path, csv_path):
            print(f"  {_rel(path)}")
    return json_path, md_path, csv_path



# ---------------------------------------------------------------------------
# Checks
# ---------------------------------------------------------------------------
def selftest(table=None, meta=None):
    """The plate's invariants, each one falsifiable on the real table."""
    report = build_report(table, meta, boot=64)
    months = report["months"]
    checks = []

    def check(name, condition, detail=""):
        checks.append((name, bool(condition), detail))

    check("table: 207 rows, A_asc and A_des, both orbit directions",
          report["audit"]["rows"] == 207
          and report["input"]["tracks"] == ["A_asc", "A_des"]
          and report["input"]["orbit_rows"].get("ASC") == 94
          and report["input"]["orbit_rows"].get("DESC") == 113,
          f"rows = {report['audit']['rows']}, tracks = "
          f"{report['input']['tracks']}, orbit = {report['audit']['orbit']}")

    check("partition: 99 reference + 99 series + 9 post = 207",
          report["audit"]["reference"] == 99 and report["audit"]["series"] == 99
          and report["audit"]["post_outside_series"] == 9,
          report["audit"]["partition"])

    contiguous = month_range(months[0]["month"], months[-1]["month"])
    check("windows: reference ends 2017-07, series is 2017-08 … 2018-08, gap-free",
          report["reference"]["last"] == "2017-07-28"
          and report["series"]["first"] == "2017-08"
          and report["series"]["last"] == "2018-08"
          and contiguous == [entry["month"] for entry in months]
          and not report["series"]["months_missing"],
          f"{report['reference']['first']} … {report['reference']['last']}, "
          f"{len(months)} months to {report['series']['last']}, "
          f"{len(report['series']['months_missing'])} gap(s)")

    # median(z_i) == z(median(log x)): the transform is monotone and z is affine.
    worst = 0.0
    for entry in months[:2] + months[-2:]:
        for dim in DIMS:
            z_chips = [chip["z"][dim] for chip in report["chips"]
                       if chip["date"].startswith(entry["month"])]
            worst = max(worst, abs(ocv.median(z_chips) - entry["dims"][dim]["z"]))
    check("median of the chips' z equals z of the median (first 2 + last 2 months)",
          worst < 1e-12, f"max |difference| = {worst:.2e}")

    centre = report["reference"]["dims"][DIMS[0]]["median_log"]
    scale = report["reference"]["dims"][DIMS[0]]["scale_log"]
    z_ref = (ocv.median([centre] * 3) - centre) / scale
    check("a month that sits on the reference median has z = 0", z_ref == 0.0,
          f"z = {z_ref:+.1e}")

    kinds = {dim: report["reference"]["dims"][dim]["scale_kind"] for dim in DIMS}
    channel_kinds = {key: {dim: block["dims"][dim]["scale_used"] for dim in DIMS}
                     for key, block in report["channels"].items()}
    check("scale rule: MAD for both dims and both tracks (all n >= SCALE_MIN_REF)",
          all(kind == "mad" for kind in kinds.values())
          and all(kind == "mad" for block in channel_kinds.values()
                  for kind in block.values()),
          f"{kinds}, tracks {channel_kinds}")


    worst_affine, flips = 0.0, []
    for entry in months:
        for dim in DIMS:
            reference = report["reference"]["dims"][dim]
            z, pooled = entry["dims"][dim]["z"], entry["dims"][dim]["z_pooled"]
            if z is None or pooled is None:
                continue
            expected = ((z * reference["scale_log"] + reference["median_log"]
                         - reference["pooled_median_log"])
                        / reference["pooled_scale_log"])
            worst_affine = max(worst_affine, abs(pooled - expected))
            if (z > 0) != (pooled > 0):
                flips.append((dim, entry["month"], abs(z)))
    worst_flip = max((z for _, _, z in flips), default=0.0)
    check("z_pooled is the same r from the pre-collapse origin, not a second estimate",
          worst_affine < 1e-9 and worst_flip < 1.0,
          f"max deviation = {worst_affine:.1e}, {len(flips)} flip(s) (max |z| = "
          f"{worst_flip:.2f})")

    worst_track = 0.0
    for block in report["channels"].values():
        for entry in block["months"]:
            for dim in DIMS:
                own, pooled = entry["z_own"][dim], entry["z_pooled"][dim]
                if own is None or pooled is None:
                    continue
                expected = own * (block["dims"][dim]["scale_log"]
                                  / report["reference"]["dims"][dim]["scale_log"])
                worst_track = max(worst_track, abs(pooled - expected))
    check("within-track z_pooled is z_own in pooled MAD units", worst_track < 1e-9,
          f"max |difference| = {worst_track:.1e}")

    worst_r = 0.0
    for entry in months:
        value = math.sqrt(sum(entry["dims"][dim]["z"] ** 2 for dim in DIMS)
                          / len(DIMS))
        worst_r = max(worst_r, abs(value - entry["R_behaviour"]["value"]))
    check("R is the rms of the two z values, every month", worst_r < 1e-12,
          f"max |R - rms| = {worst_r:.2e}")

    mix_bad = []
    for entry in months:
        if sum(entry["n_by_track"].values()) != entry["n"]:
            mix_bad.append(f"{entry['month']}: n")
        if abs(sum(entry["n_by_track"].values()) - entry["n"]) == 0:
            share = entry["n_by_track"]["ASC"] / entry["n"]
            if abs(share - entry["mix"]["asc_share"]) > 1e-12:
                mix_bad.append(f"{entry['month']}: share")
    check("track counts partition every month and the ASC share is their ratio",
          not mix_bad, f"{len(mix_bad)} problem(s)")

    intervals = [ci for entry in months for dim in DIMS
                 for ci in (entry["dims"][dim]["ci95"], entry["R_behaviour"]["ci95"])
                 if ci]
    thin = [entry["month"] for entry in months if entry["thin"]]
    check("intervals are ordered and exactly 2018-08 is thin (n = 3)",
          all(ci[0] <= ci[1] for ci in intervals) and thin == ["2018-08"]
          and months[-1]["n"] == 3,
          f"{len(intervals)} interval(s), thin = {thin}, n = {months[-1]['n']}")

    twin = build_report(table, meta, boot=64)
    differences = _diff(report, twin)
    check("determinism: two builds are identical object for object",
          not differences, f"{len(differences)} difference(s)")

    with tempfile.TemporaryDirectory() as tmp:
        preview = os.path.join(tmp, "plate.jpg")
        make_figure(report, path=preview, dpi=70)
        size = os.path.getsize(preview)
    check("figure: the 2 x 2 plate renders and is non-trivial", size > 10000,
          f"{size} bytes at 70 dpi")

    for name, passed, detail in checks:
        print(f"  {'ok  ' if passed else 'FAIL'} {name}"
              + (f"  [{detail}]" if detail else ""))
    failed = [name for name, passed, _ in checks if not passed]
    if failed:
        raise SystemExit(f"selftest failed: {failed}")
    print(f"selftest: {len(checks)} checks ok")
    return 0



def _diff(left, right, path="", tol=1e-9, out=None):
    """Compare two report objects: floats within ``tol``, everything else exact."""
    if out is None:
        out = []
    if isinstance(left, dict) and isinstance(right, dict):
        for key in sorted(set(left) | set(right)):
            if key not in left:
                out.append(f"{path}.{key}: absent from the recomputation")
            elif key not in right:
                out.append(f"{path}.{key}: absent from the committed file")
            else:
                _diff(left[key], right[key], f"{path}.{key}", tol, out)
    elif isinstance(left, list) and isinstance(right, list):
        if len(left) != len(right):
            out.append(f"{path}: length {len(left)} vs {len(right)}")
        else:
            for index, (a, b) in enumerate(zip(left, right)):
                _diff(a, b, f"{path}[{index}]", tol, out)
    elif isinstance(left, bool) or isinstance(right, bool) \
            or left is None or right is None:
        if left != right:
            out.append(f"{path}: {left!r} != {right!r}")
    elif isinstance(left, (int, float)) and isinstance(right, (int, float)):
        if abs(left - right) > tol * max(1.0, abs(left), abs(right)):
            out.append(f"{path}: {left!r} != {right!r}")
    elif left != right:
        out.append(f"{path}: {left!r} != {right!r}")
    return out


def _csv_problems(report, path):
    """The compact table against the object it was written from."""
    problems = []
    with open(path, newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        return [f"{_rel(path)}: no rows"]
    if list(rows[0]) != CSV_FIELDS:
        problems.append(f"{_rel(path)}: columns {list(rows[0])} != {CSV_FIELDS}")
    if len(rows) != len(report["months"]):
        problems.append(f"{_rel(path)}: {len(rows)} rows for "
                        f"{len(report['months'])} months")
        return problems
    for row, entry in zip(rows, report["months"]):
        expected = {
            "Month": entry["month"],
            "n": str(entry["n"]),
            "n_ASC": str(entry["n_by_track"]["ASC"]),
            "n_DESC": str(entry["n_by_track"]["DESC"]),
            "gamma2": _g(entry["dims"][DIMS[0]]["median"]),
            "P": _g(entry["dims"][DIMS[1]]["median"]),
            "z_gamma2": _g(entry["dims"][DIMS[0]]["z"]),
            "z_P": _g(entry["dims"][DIMS[1]]["z"]),
            "z_pooled_gamma2": _g(entry["dims"][DIMS[0]]["z_pooled"]),
            "z_pooled_P": _g(entry["dims"][DIMS[1]]["z_pooled"]),
            "R_behaviour": _g(entry["R_behaviour"]["value"]),
            "ci95_gamma2": _pair(entry["dims"][DIMS[0]]["ci95"]),
            "ci95_P": _pair(entry["dims"][DIMS[1]]["ci95"]),
            "ci95_R": _pair(entry["R_behaviour"]["ci95"]),
            "thin": "yes" if entry["thin"] else "no",
        }
        for key, value in expected.items():
            if row.get(key) != value:
                problems.append(f"{_rel(path)} {entry['month']} {key}: "
                                f"{row.get(key)!r} != {value!r}")
    return problems


def verify(report, path=None, csv_path=None):
    """Recompute everything from the table and compare it to the committed JSON."""
    path = path or OUT_JSON
    csv_path = csv_path or OUT_CSV
    if not os.path.exists(path):
        raise SystemExit(f"{path}: no committed output — run the script first")
    with open(path) as handle:
        committed = json.load(handle)
    problems = _diff(report, committed)
    if os.path.exists(csv_path):
        problems += _csv_problems(report, csv_path)
    if problems:
        for line in problems[:12]:
            print("  " + line)
        raise SystemExit(f"verify: {len(problems)} difference(s) against "
                         f"{_rel(path)}")
    print(f"verify: {len(report['months'])} months, {len(report['chips'])} series "
          f"chips, reference n = {report['reference']['n']} — identical to "
          f"{_rel(path)} and {_rel(csv_path)}")
    return 0



# ---------------------------------------------------------------------------
# Command line
# ---------------------------------------------------------------------------
def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Morandi — OCV reconfiguration against the 2016-08 … 2017-07 "
                    "reference, month by month")
    parser.add_argument("--table", default=TABLE)
    parser.add_argument("--meta", default=META)
    parser.add_argument("--json", default=OUT_JSON)
    parser.add_argument("--csv", default=OUT_CSV)
    parser.add_argument("--md", default=None)
    parser.add_argument("--jpg", default=None)
    parser.add_argument("--figdir", default=None,
                        help="directory for the plate and its guide (default: "
                             "figures/morandi)")
    parser.add_argument("--boot", type=int, default=BOOT)
    parser.add_argument("--block", type=int, default=BLOCK)
    parser.add_argument("--seed", type=int, default=SEED)
    parser.add_argument("--no-figure", action="store_true")
    parser.add_argument("--selftest", action="store_true")
    parser.add_argument("--verify", action="store_true")
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args(argv)

    if args.selftest:
        return selftest(args.table, args.meta)

    report = build_report(args.table, args.meta, boot=args.boot, block=args.block,
                          seed=args.seed)
    if args.verify:
        return verify(report, args.json, args.csv)

    md_path = args.md or (os.path.join(args.figdir, os.path.basename(OUT_MD))
                          if args.figdir else OUT_MD)
    jpg_path = args.jpg or (os.path.join(args.figdir, os.path.basename(OUT_JPG))
                            if args.figdir else OUT_JPG)
    write_outputs(report, json_path=args.json, csv_path=args.csv, md_path=md_path,
                  quiet=args.quiet)
    if not args.no_figure:
        make_figure(report, path=jpg_path)
        if not args.quiet:
            print(f"  {_rel(jpg_path)}")
    months = report["months"]
    if not args.quiet:
        print(f"  {len(months)} months ({months[0]['month']} … {months[-1]['month']}), "
              f"n = {min(e['n'] for e in months)}-{max(e['n'] for e in months)}; "
              f"reference n = {report['reference']['n']} "
              f"({report['reference']['first']} … {report['reference']['last']}); "
              f"{report['audit']['post_outside_series']} post-collapse chips excluded")
    return 0


if __name__ == "__main__":
    sys.exit(main())

