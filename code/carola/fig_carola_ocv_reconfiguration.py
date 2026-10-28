#!/usr/bin/env python3
"""Carola — OCV reconfiguration in units of the 2022-09 … 2023-12 reference.

Question: how far is each month of the Carolabrücke record's **last nine months
before the collapse** (2024-01 … 2024-09) from the bridge's **own earlier
state** (2022-09 … 2023-12), on the two dimensions of the ``behaviour`` block of
``code/sds/ocv_reconfigure_spec.md`` — ``gamma2`` (the consistency of the echo
across acquisitions) and ``P`` (its persistence) — and does the approach to
2024-09-11 show as a departure from that state?

The metric is the pre-registered one, imported rather than rewritten:

    r = median(log x_month) - median(log x_ref)          g = log (both dimensions)
    z = r / scale_ref                                    scale_ref = MAD(log x_ref)
    R = rms(z(gamma2), z(P))

so the reference is the origin **(0, 0)** by construction, exactly as on the
``fig_sds_ocv_reconfigure*`` plates.

The reference period is a decision, not a default — Option B
------------------------------------------------------------
The record starts 2022-09-03. The candidates were a short recent baseline
(2023-01 … 2023-12), a seasonally balanced one (2022-09 … 2023-12), or the whole
pre-collapse span (2022-09 … 2024-09-09, the SDS plate's reference). Option B is
used here because it is the only one that is both **long enough** to be a
population (121 chips over 16 months, 81 ASCENDING / 40 DESCENDING, every
calendar month of a full year plus four) and **earlier than the series**: 2024-01
… 2024-09 is held out entirely, so the series cannot set its own origin. The
same calendar months appear in the reference (2023-01 … 2023-09), so season is
represented; it is *pooled*, not resolved, because a single calendar month holds
~7 chips and 7 < ``SCALE_MIN_REF = 10`` — a per-month baseline has no MAD. That
limit is stated in the guide rather than engineered around.

Two consequences are carried into every number rather than into a footnote:

* the reference is **121 chips of 2022-09-03 … 2023-12** (16 months, no empty
  month), and the 63 chips of 2024-01-01 … 2024-09-11 are the series — 7 or 8
  chips per month, except 2024-09, which has **2** and is the only ``thin`` month
  (``SMALL_N = 6``);
* the only interval reported is the **reference's** own uncertainty: 2 000
  circular moving-block draws (block 3, seed 7) of the 121 reference chips, with
  the reference median *and* scale re-estimated inside every draw, the month's
  chips held fixed. It answers "how well is the origin known", never "how well is
  the month known" — and for the n = 2 month it is narrow for the wrong reason.

Because one reference window is a choice, the **channel split** travels beside
the pooled number: the same months re-standardised by the ASCENDING-only
reference (n = 81) and the DESCENDING-only one (n = 40), drawn as dashed lines in
panels (a) and (b) and tabulated in the guide. They are the same r read from a
different origin — a sensitivity, not a second estimate.

Nothing is tested. A month axis is not a dose and one bridge has no partner, so
this figure carries no H1, no H2, no p-value and no precursor claim; the
252 post-collapse chips (this record continues to 2026-08-23) are counted and
never plotted.

Figure (2 x 2, 13.0 in wide, 600 dpi jpg — the sds package's ``savefig.dpi``):
  a  z(gamma2) by month, per-chip dots, monthly median, reference band ±2 MAD,
     the channel-split redraws as dashed lines, a chip-count strip, the collapse
     cut; n on every month and the thin month ringed
  b  z(P) by month, same standardisation and axis
  c  the (z gamma2, z P) plane — the origin is the reference state, the grey
     cloud is the 121 reference chips, the path is 2024-01 -> 2024-09 coloured by
     quarter, the arrows are the quarter medians, the terminal month is ringed
  d  the block magnitude R = rms(z gamma2, z P) per month, with the empirical
     5-95 % band of the reference chips as the null band

Outputs (all inside the Carola package, beside its own data and figures):
  data/carola/ocv_reconfiguration_carola.json
  data/carola/ocv_reconfiguration_carola_table.csv
  figures/carola/fig_carola_ocv_reconfiguration_carola.jpg
  figures/carola/fig_carola_ocv_reconfiguration_carola.md

Usage:
  python3 code/carola/fig_carola_ocv_reconfiguration.py
  python3 code/carola/fig_carola_ocv_reconfiguration.py --selftest
  python3 code/carola/fig_carola_ocv_reconfiguration.py --verify
  python3 code/carola/fig_carola_ocv_reconfiguration.py --figdir /tmp/x --json /tmp/y.json
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

import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as mdates  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, os.pardir, "sds"))

import carola_ocv_core as core      # noqa: E402
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
SITE = "carola"
SITE_LABEL = core.SITE

#: The behaviour block of the spec, and the column each dimension lives in on the
#: support-sized chip of this site.
DIMS = tuple(ocv.BLOCKS["behaviour"])          # ("gamma2", "P")
COLUMN = {"gamma2": "gamma2_mask", "P": "P"}
TRANSFORM = {dim: "log" for dim in DIMS}

#: Option B: the reference is the first 16 months of the record (2022-09 …
#: 2023-12, seasonally balanced), the series is the last nine months before the
#: collapse (2024-01 … 2024-09-11). The two windows are disjoint by construction.
REF_START, REF_END = "2022-09-01", "2023-12-31"
SERIES_START = "2024-01-01"
EVENT = core.COLLAPSE_DATE                      # 2024-09-11

#: The "pooled" sensitivity re-centres on the **pre-collapse** record (reference
#: + series), not on the whole table: 252 of the 436 chips are post-collapse, and
#: a MAD dominated by the collapsed state is not a spread of the healthy one.
POOL_END = EVENT

TABLE = os.path.join(core.DATA, "support", "carola_ocv_channels.csv")
META = os.path.join(core.DATA, "support", "carola_ocv_channels_meta.json")
OUT_JSON = os.path.join(core.DATA, "ocv_reconfiguration_carola.json")
OUT_CSV = os.path.join(core.DATA, "ocv_reconfiguration_carola_table.csv")
OUT_MD = os.path.join(core.FIGDIR, "fig_carola_ocv_reconfiguration_carola.md")
OUT_JPG = os.path.join(core.FIGDIR, "fig_carola_ocv_reconfiguration_carola.jpg")

BOOT = 2000            # reference-uncertainty draws (the spec's seed and count)
BLOCK = 3              # consecutive reference chips per draw (spec: BLOCK_BOOT)
SEED = 7
SMALL_N = ocv.SMALL_N  # 6 — below this a month is thin and must be quoted with n
SCALE_MIN_REF = ocv.SCALE_MIN_REF  # 10 — below this the pooled scale takes over

#: One colour per quarter of the series; the series is a single calendar year, so
#: the path is coloured by quarter instead of by year.
QUARTER_COLOUR = {"Q1": "#1f4e79", "Q2": "#2f855a", "Q3": "#b7791f"}
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

    The pre-registered ``ocv.transform`` is keyed by a closed site registry that
    does not contain ``carola`` — adding it would touch the pinned sds run — so
    the two-line transform is spelled out here while the *scale machinery* is
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
    """The support-sized chip of one Carola acquisition per row, validated.

    The committed table carries **Girder A** only (436 rows, ``segment =
    Girder A``, support band 39/15 on a 69 x 125 chip), acquired in both orbit
    directions. A row from another girder would mix footprints, and a row whose
    ``reference`` is not the ring would carry a different support rule, so the
    run stops rather than averaging either in. ``day`` is the table's date
    column; it is copied to ``date`` so the rest of the module reads like the
    sibling CTS plate.
    """
    path = path or TABLE
    with open(path, newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise SystemExit(f"{path}: no rows")
    needed = {"day", "orbit_direction", "segment", "reference"} | set(COLUMN.values())
    missing = needed - set(rows[0])
    if missing:
        raise SystemExit(f"{path}: missing columns {sorted(missing)}")
    segments = {row["segment"] for row in rows}
    if segments != {"Girder A"}:
        raise SystemExit(
            f"{path}: segment(s) {sorted(segments)} — this run reads Girder A's "
            f"own support-sized chip only, never a second girder")
    rules = {row["reference"] for row in rows}
    if rules != {"ring"}:
        raise SystemExit(f"{path}: reference rule(s) {sorted(rules)} — the ring "
                         f"support of the committed table, one expected")
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
    return f"Q{(int(month[5:7]) - 1) // 3 + 1}"


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


def _rms(z_by_dim):
    values = [z_by_dim[d] for d in DIMS if z_by_dim.get(d) is not None]
    if not values:
        return None
    return math.sqrt(sum(v * v for v in values) / len(values))


# ---------------------------------------------------------------------------
# The report — one pass over the table, everything else reads this object
# ---------------------------------------------------------------------------
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
            "pooled_median_log": pooled["median"],
            "pooled_scale_log": pooled["scale"],
            "pooled_scale_kind": pooled["scale_kind"],
        }
    return dims


def _channel_blocks(ref_rows, series_rows):
    """The same months re-standardised by the ASCENDING- and DESCENDING-only refs.

    A sensitivity to *which chips* the origin is built from, not a second
    estimate: the month's chips are the same, only the median and MAD of the
    origin change. Both channels stay above ``SCALE_MIN_REF`` here (81 and 40),
    so each keeps its own MAD; the fallback is spelled out anyway so a future
    table with a thin channel degrades instead of dividing by zero.
    """
    pool_rows = ref_rows + series_rows
    out = {}
    for direction in sorted({row["orbit_direction"] for row in ref_rows}):
        key = ocv.CHANNEL[direction]
        ref_sel = [row for row in ref_rows if row["orbit_direction"] == direction]
        dims = {}
        for dim in DIMS:
            stats = _stats([_log(_number(row[COLUMN[dim]])) for row in ref_sel])
            pooled = _stats([_log(_number(row[COLUMN[dim]]))
                             for row in pool_rows
                             if row["orbit_direction"] == direction])
            used = "mad" if stats["n"] >= SCALE_MIN_REF or not pooled["scale"] \
                else "pooled_mad"
            scale = stats["scale"] if used == "mad" else pooled["scale"]
            dims[dim] = {"n": stats["n"], "median_log": stats["median"],
                         "median": math.exp(stats["median"]), "mad_log": stats["mad"],
                         "scale_log": scale, "scale_used": used}
        groups = {}
        for row in series_rows:
            if row["orbit_direction"] == direction:
                groups.setdefault(month_of(row["date"]), []).append(row)
        months = []
        for month in sorted(groups):
            group = groups[month]
            z = {}
            for dim in DIMS:
                stats = _stats([_log(_number(row[COLUMN[dim]])) for row in group])
                z[dim] = None if stats["median"] is None else (
                    (stats["median"] - dims[dim]["median_log"]) / dims[dim]["scale_log"])
            months.append({"month": month, "n": len(group), "z": z,
                           "R": _rms(z), "dates": sorted(row["date"] for row in group)})
        out[key] = {
            "channel": key,
            "orbit_direction": direction,
            "n_ref": len(ref_sel),
            "n_series": sum(len(group) for group in groups.values()),
            "dims": dims,
            "months": months,
        }
    return out


def build_report(table=None, meta=None, boot=BOOT, block=BLOCK, seed=SEED):
    """The whole plate as one object: reference, months, chips, channels, audit."""
    table = table or TABLE
    meta_path = meta or META
    rows = sorted(load_rows(table), key=lambda row: (row["date"], row["orbit_direction"]))
    with open(meta_path) as handle:
        meta_obj = json.load(handle)
    directions = sorted({row["orbit_direction"] for row in rows})
    orbit_label = " + ".join(ocv.CHANNEL[direction] for direction in directions)

    ref_rows = [row for row in rows if REF_START <= row["date"] <= REF_END]
    series_rows = [row for row in rows if SERIES_START <= row["date"] <= EVENT]
    post_rows = [row for row in rows if row["date"] > EVENT]
    if not ref_rows or not series_rows:
        raise SystemExit("the reference or the series is empty — check the window")
    if len(ref_rows) + len(series_rows) + len(post_rows) != len(rows):
        raise SystemExit("the three windows do not partition the table")
    if ref_rows[-1]["date"] >= series_rows[0]["date"]:
        raise SystemExit(
            f"the reference ends {ref_rows[-1]['date']}, the series starts "
            f"{series_rows[0]['date']} — the two windows overlap")

    ref_values = {dim: [_log(_number(row[COLUMN[dim]])) for row in ref_rows]
                  for dim in DIMS}
    pool_rows = ref_rows + series_rows          # the pre-collapse record
    pool_values = {dim: [_log(_number(row[COLUMN[dim]])) for row in pool_rows]
                   for dim in DIMS}
    all_values = {dim: [_log(_number(row[COLUMN[dim]])) for row in rows] for dim in DIMS}
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
        return {"date": row["date"], "z": z, "R": _rms(z)}

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

    months = []
    for key in keys:
        group = groups[key]
        entry = {"month": key, "n": len(group),
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
        months.append(entry)

    # The channel split: how much of each month's z is the *origin's* choice.
    channels = _channel_blocks(ref_rows, series_rows)
    lookup = {entry["month"]: entry for entry in months}
    for block in channels.values():
        shifts, flips, outside = [], [], []
        for entry in block["months"]:
            base = lookup[entry["month"]]["dims"]
            for dim in DIMS:
                z_channel, z_ref = entry["z"][dim], base[dim]["z"]
                if z_channel is None or z_ref is None:
                    continue
                shifts.append(abs(z_channel - z_ref))
                if (z_channel > 0) != (z_ref > 0):
                    flips.append(f"{dim} {entry['month']}")
                if abs(z_channel) > 2.0:
                    outside.append(f"{dim} {entry['month']}")
        block["summary"] = {
            "max_abs_shift": max(shifts, default=0.0),
            "sign_flips": flips,
            "outside_2mad": outside,
            "R_max": max((entry["R"] for entry in block["months"]
                          if entry["R"] is not None), default=None),
        }

    post_dates = sorted(row["date"] for row in post_rows)
    ref_months_missing = missing_reference(ref_rows)
    ref_by_channel = {ocv.CHANNEL[direction]:
                      sum(1 for row in ref_rows if row["orbit_direction"] == direction)
                      for direction in directions}
    series_by_channel = {ocv.CHANNEL[direction]:
                         sum(1 for row in series_rows
                             if row["orbit_direction"] == direction)
                         for direction in directions}
    report = {
        "schema": SCHEMA,
        "site": SITE,
        "site_label": SITE_LABEL,
        "event": {"collapse": core.COLLAPSE_DATE, "record_end": post_dates[-1]},
        "metric": {
            "definition": ("z = (median(log x_month) - median(log x_ref)) / scale_ref; "
                           "R_behaviour = rms over the behaviour block"),
            "transform": {dim: TRANSFORM[dim] for dim in DIMS},
            "scale_rule": ("MAD of the reference, fallback chain MAD -> IQR/1.349 -> SD; "
                           f"the pooled scale of the pre-collapse record replaces it when "
                           f"n_ref < {SCALE_MIN_REF} (spec amendment 2026-10-22)"),
            "block": list(DIMS),
            "unit": "one support-sized chip per acquisition (Girder A)",
            "interval": (f"reference-uncertainty bootstrap only: {boot} circular "
                         f"moving-block draws, block {block}, seed {seed}, reference "
                         f"chips resampled, the month's chips held fixed"),
            "pooled_window": (f"pre-collapse record, {REF_START} … {POOL_END} "
                              f"(reference + series)"),
            "reference_choice": ("Option B: the 16 months 2022-09 … 2023-12, seasonally "
                                 "balanced and earlier than the series; not the SDS "
                                 "plate's whole pre-collapse reference"),
            "tests": None,
            "thin_rule": f"n < {SMALL_N} — 2024-09 only",
        },
        "input": {
            "table": _rel(table),
            "table_sha256": sha256(table),
            "meta": _rel(meta_path),
            "meta_sha256": sha256(meta_path),
            "rows": len(rows),
            "segment": sorted({row["segment"] for row in rows}),
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
        },
        "months": months,
        "chips": chips,
        "reference_chips": ref_chips,
        "reference_R_band": ref_band,
        "channels": channels,
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
    """``1 gamma2 month (2024-02)`` / ``no P month`` — never a bare zero."""
    if not keys:
        return f"no {dim} month"
    return f"{len(keys)} {dim} month(s) ({_month_list(keys)})"


def _arc(months, dim):
    """The largest change of z over two months, and the months it spans.

    With nine months this is the only *shape* statement the series supports: a
    three-point window is the shortest span in which a monotone move could show
    up at all, and even then it cannot be separated from the month-to-month
    spread.
    """
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
    above = [m["month"] for m in months if m["R_behaviour"]["value"] > band[1]]
    below = [m["month"] for m in months if m["R_behaviour"]["value"] < band[0]]
    r_max = max(months, key=lambda m: m["R_behaviour"]["value"])
    z_max = max(abs(m["dims"][dim]["z"]) for m in months for dim in DIMS)
    r_values = [m["R_behaviour"]["value"] for m in months]
    arc = _arc(months, g_dim)
    arc_text = ""
    if arc:
        delta, start, end = arc
        arc_text = (
            f" The largest two-month change is {_signed(delta)} of {g_dim} between "
            f"{start['month']} and {end['month']} ({_signed(start['dims'][g_dim]['z'])}"
            f" -> {_signed(end['dims'][g_dim]['z'])}); over nine months that is a "
            f"stretch, not a trend.")
    channels = report["channels"]
    channel_bits = []
    flip_z = []
    for key in sorted(channels):
        block = channels[key]
        summary = block["summary"]
        counts = [entry["n"] for entry in block["months"]]
        channel_bits.append(
            f"{key} (n_ref = {block['n_ref']}; {min(counts)}-{max(counts)} chips "
            f"per month) shifts a month by at most "
            f"{summary['max_abs_shift']:.2f} MAD against the pooled origin, flips "
            f"{len(summary['sign_flips'])} sign(s) and puts "
            f"{len(summary['outside_2mad'])} month-dimension cell(s) beyond ±2")
        for name in summary["sign_flips"]:
            dim, month = name.split()
            entry = next(m for m in months if m["month"] == month)
            if entry["dims"][dim]["z"] is not None:
                flip_z.append(abs(entry["dims"][dim]["z"]))
    return [
        f"**No month of the series is a departure from the reference.** The largest "
        f"|z| anywhere in {months[0]['month']} … {last['month']} is {z_max:.2f}, and "
        f"beyond ±2 MAD lie {_outside_text(g_dim, outside[g_dim])} and "
        f"{_outside_text(p_dim, outside[p_dim])}. The reference band on R is "
        f"{_plain(band[0])}–{_plain(band[1])}; R peaks at "
        f"{_plain(r_max['R_behaviour']['value'])} in {r_max['month']} "
        f"(n = {r_max['n']}) and "
        + (f"exceeds the band in {_month_list(above)}." if above
           else "never reaches its 95th percentile.")
        + " On this record, a bridge two days before its collapse reads as the "
        "bridge of sixteen months earlier.",
        "",
        f"**The largest single excursion is a low, not a high.** {r_max['month']} "
        f"carries the largest R ({_plain(r_max['R_behaviour']['value'])}) with "
        f"z({g_dim}) = {_signed(r_max['dims'][g_dim]['z'])} "
        f"{_ci(r_max['dims'][g_dim]['ci95'])} and z({p_dim}) = "
        f"{_signed(r_max['dims'][p_dim]['z'])} "
        f"{_ci(r_max['dims'][p_dim]['ci95'])} — both *below* the origin, i.e. less "
        f"echo consistency and less persistence than 2022-09 … 2023-12, not more. "
        f"Its two intervals exclude zero because the month's median is far from the "
        f"origin; that is a statement about where the origin sits, not about the "
        f"month's typicality, since a single reference chip sits ~2 MAD out by "
        f"construction.",
        "",
        f"**The series oscillates on a two-to-three-month scale.** {g_dim} swings "
        f"from {_signed(g_low['dims'][g_dim]['z'])} ({g_low['month']}) to "
        f"{_signed(g_high['dims'][g_dim]['z'])} ({g_high['month']}) and {p_dim} from "
        f"{_signed(p_low['dims'][p_dim]['z'])} ({p_low['month']}) to "
        f"{_signed(p_high['dims'][p_dim]['z'])} ({p_high['month']}).{arc_text}",
        "",
        f"**R stays bounded, and its smallest value is an artefact.** R runs "
        f"{_plain(min(r_values))} … {_plain(max(r_values))} against a reference band "
        f"of {_plain(band[0])}–{_plain(band[1])}. "
        + (f"{_month_list(below)} falls *below* the band, but it holds "
           f"{last['n']} chips ({', '.join(last['dates'])}): two chips can agree with "
           f"the origin on both dimensions by luck, so a low R there is a small "
           f"sample, not a stabilised bridge. " if below else "")
        + "R is an RMS over two dimensions, so one coordinate can carry it — which "
        "is why panel (c), not R alone, is the reading.",
        "",
        f"**The last month before the collapse is not measurable to better than its "
        f"n.** {last['month']} holds {last['n']} chips ({', '.join(last['dates'])}), "
        f"the last of them two days before the event, and is the only month below "
        f"`SMALL_N = {SMALL_N}`. Its marker is z({g_dim}) = {_signed(last_g['z'])} "
        f"{_ci(last_g['ci95'])}, z({p_dim}) = {_signed(last_p['z'])} "
        f"{_ci(last_p['ci95'])}, R = {_plain(last['R_behaviour']['value'])}. The "
        f"interval is narrow because the *origin* is well known, not because the "
        f"month is: with two chips, the month's own sampling noise is not bounded by "
        f"anything this figure computes.",
        "",
        "**The channel split is the widest uncertainty in the plate.** "
        + "; ".join(channel_bits) + ". "
        + (f"Every sign flip sits within {max(flip_z):.2f} MAD of zero, so no "
           f"excursion moves, but the *sign* of a month close to the origin does "
           f"depend on which orbit builds it. The DESC redraw rests on 2-3 chips a "
           f"month, so part of that spread is the month's own n rather than the "
           f"origin's composition — this figure cannot separate the two, and does "
           f"not try. The pooled ASC+DESC reference is the pre-registered one; the "
           f"split is declared, not preferred."
           if flip_z else
           "No month changes sign under either channel reference."),
        "",
        f"**One sentence.** Between {report['series']['first_chip']} and "
        f"{report['series']['last_chip']} the bridge oscillates inside the state it "
        f"held from {report['reference']['first']} to "
        f"{report['reference']['last']}, and nothing in this figure separates the "
        f"last month from the first.",
    ]


def caveat_lines(report):
    """The limits, stated where the figure is, not only in the JSON."""
    months = report["months"]
    last = months[-1]
    low_n = min(m["n"] for m in months)
    high_n = max(m["n"] for m in months)
    dims = report["reference"]["dims"]
    g_dim, p_dim = DIMS
    channels = report["channels"]
    channel_text = ", ".join(
        f"{key}: {block['n_ref']} reference chips, {block['n_series']} in the series"
        for key, block in sorted(channels.items()))
    return [
        "* **The binding limit is sampling, not SNR.** A month rests on "
        f"{low_n}-{high_n} chips, so 2024-09 alone is `thin` "
        f"(`n = {last['n']} < SMALL_N = {SMALL_N}`) and its \"median\" *is* the "
        "midpoint of two chips. Every other month is 7 or 8 acquisitions, which is "
        "enough for a median and not enough for a spread: the monthly marker is the "
        "median of the chips that exist, not of the month, and a month with no "
        "acquisition would be a gap, not a value.",
        "* **The interval belongs to the reference.** It says how well a "
        f"{report['reference']['n']}-chip origin is known under a resampling of "
        "those chips. It is not a confidence interval for a monthly value: the "
        "month's own sampling noise is not in it. That is why 2024-09's interval is "
        "*narrow* — the origin is well known, the month is not — and why an interval "
        "that excludes zero (2024-02) still describes a month inside the spread of "
        "the reference chips themselves.",
        "* **The reference is seasonally pooled, not seasonal.** Its 16 months "
        "contain every calendar month of 2023 once and 2022's September to December "
        "once, so the ±2 band mixes seasonal spread with state spread. A baseline "
        "per calendar month would need ~7 chips per month, below "
        f"`SCALE_MIN_REF = {SCALE_MIN_REF}`: such a reference has no MAD. The band "
        "is therefore *wide for a reason*, and a null inside it is weaker evidence "
        "than a null inside a tight band would be.",
        "* **The origin is a choice, and the channel split prices it.** "
        f"{channel_text}. Rebuilding the origin from one orbit "
        f"moves a month by up to "
        f"{max(b['summary']['max_abs_shift'] for b in channels.values()):.2f} MAD and "
        f"puts "
        f"{sum(len(b['summary']['outside_2mad']) for b in channels.values())} "
        "month-dimension cell(s) beyond ±2 — the difference between \"no signal\" and "
        "\"two cells outside the band\" is which chips the reference is made of. The "
        "pooled ASC+DESC reference is the one declared in advance; the split is a "
        "declared sensitivity, not a menu.",
        f"* **This `{COLUMN['gamma2']}` is not the `gamma2` of the months figure.** "
        "The support-sized chip here is 69 x 125 cells with 585 support cells "
        f"({report['input']['support_m']['guard_m']:g} m guard, "
        f"{report['input']['support_m']['ring_m']:g} m ring), and its median "
        f"`{COLUMN['gamma2']}` is {dims['gamma2']['median']:.4g} over the reference. "
        "The echo-mask γ² of `fig_carola_ocv_months` is computed on the full 80 x 80 "
        "girder chip and runs an order of magnitude higher. Only the *standardised* "
        "quantity z is comparable between the two figures; the levels are not, and "
        "neither is R.",
        "* **The series stops at the collapse, deliberately.** The "
        f"{report['audit']['post_outside_series']} acquisitions after "
        f"{report['event']['collapse']} ({report['audit']['post_first']} … "
        f"{report['audit']['post_last']}) are the same bridge in a different state, "
        "and they are counted in the audit and carry no value here. The record does "
        "not end at the event; the *series* does.",
        "* **No precursor claim.** The metric is descriptive: it says where the last "
        f"available month ({last['month']}) sits, and that month is on the "
        "reference. The last chip is 2024-09-09 and the next acquisition is after "
        "the event, so nothing between them was measured; a series of 7-8 chips a "
        "month cannot time, detect or exclude a precursor, and no test is run that "
        "would pretend otherwise.",
        "* **No test, no dose, no partner.** The sds H1/H2 machinery needs at least "
        "two states per block and, for a dose, an ordered ladder; here one state is "
        "sampled monthly. So there is no p-value, no Holm family and no cosine, "
        "±2 MAD is a descriptive band and not a significance threshold, and these "
        "magnitudes are in this reference's own units — they are not comparable to "
        "the `fig_sds_ocv_reconfigure_carola` amplitudes, whose reference is the "
        "whole pre-collapse record (2022-09 … 2024-09-09).",
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
        "# Carola — OCV reconfiguration against the 2022-09 … 2023-12 reference, "
        "month by month",
        "",
        "*generated by `code/carola/fig_carola_ocv_reconfiguration.py`; the metric is "
        "the pre-registered one of `code/sds/ocv_reconfigure_spec.md` (behaviour "
        "block, log transform, reference MAD), read here as a monthly series. It is "
        "not a new test, and it carries none.*",
        "",
        f"{report['site_label']}: one support-sized chip per acquisition "
        f"({chip.get('rows')} x {chip.get('cols')} cells, "
        f"{chip.get('support_cells', 585)} support cells, "
        f"{report['input']['support_m']['guard_m']:g} m guard, "
        f"{report['input']['support_m']['ring_m']:g} m ring) of **Girder A**, "
        f"{report['audit']['rows']} acquisitions in both orbit directions "
        f"({report['input']['orbit_rows'].get('ASC')} ASC, "
        f"{report['input']['orbit_rows'].get('DESC')} DESC). The record splits into "
        f"a **reference** of {ref['n']} chips ({ref['first']} … {ref['last']}), the "
        f"**series** of {report['series']['chips']} chips on "
        f"{report['series']['months']} months ({report['series']['first']} … "
        f"{report['series']['last']}), and "
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
        f"{report['series']['months']} months, 2-8 per month",
        "```",
        "",
        "The reference is drawn from the 16 months *before* the series, so the origin "
        "of everything below is the bridge's own earlier state and not a baseline "
        "chosen after the fact. `0` on both axes is that state; it is not a norm and "
        "not a threshold.",
        "",
        "## Why this window",
        "",
        "The record holds 436 chips from 2022-09-03; the series is the last nine "
        "months before the collapse, so the reference can only be built from what "
        "precedes 2024-01. Three windows were available, and the choice is stated "
        "rather than buried because it moves the answer:",
        "",
        "| window | chips | months | ASC / DESC | why it is, or is not, used |",
        "| --- | ---: | ---: | --- | --- |",
        "| 2023-01 … 2023-12 | 91 | 12 | 61 / 30 | one calendar year, but the shortest "
        "reference and it would discard the 2022 chips for no gain in balance |",
        f"| **2022-09 … 2023-12** | **{ref['n']}** | **{ref['months']}** | "
        f"**{ref['n_by_channel']['ASC']} / {ref['n_by_channel']['DESC']}** | "
        "**used** — 16 months, every calendar month of a year represented, and the "
        "series held out entirely |",
        "| 2022-09 … 2024-09-09 | 184 | 25 | 124 / 60 | the `fig_sds_ocv_reconfigure_"
        "carola` reference: it *contains* the series, so it cannot be the origin the "
        "series is measured against |",
        "",
        "The third window is not wrong — it is the right choice for the SDS plate, "
        "which asks how the pre-collapse block differs from the post-collapse one. "
        "Here the question is a month against an earlier state, so the origin has to "
        "predate the series. The cost of that decision is stated with it: a 2022-09 "
        "… 2023-12 origin knows nothing about the 2024 acquisition conditions except "
        "through the seasons it shares with 2023.",
        "",
        "## The reference is the origin",
        "",
        "| dim | column | transform | n_ref | median (log) | median | MAD (log) | "
        "1 MAD = | pooled MAD (log) | scale used |",
        "| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |",
    ]
    for dim in DIMS:
        entry = dims[dim]
        lines.append(
            f"| {dim} | `{entry['column']}` | {entry['transform']} | {entry['n']} | "
            f"{_plain(entry['median_log'], 4)} | {entry['median']:.4g} | "
            f"{_plain(entry['mad_log'], 4)} | x{entry['one_mad_factor']:.3f} | "
            f"{_plain(entry['pooled_scale_log'], 4)} | `{entry['scale_used']}` |")
    lines += [
        "",
        f"Both dimensions keep the reference MAD: n_ref = {ref['n']} is above the "
        f"`SCALE_MIN_REF = {SCALE_MIN_REF}` of the 2026-10-22 amendment, so the pooled "
        f"MAD of the pre-collapse record ({report['audit']['pooled_chips']} chips, "
        f"{report['metric']['pooled_window']}) is **not** substituted. It travels "
        f"beside the reference as the declared sensitivity — `z_pooled` in the JSON "
        f"reads the same r from that pooled median and MAD, an affine image of z "
        f"rather than a second estimate. One MAD of {g_dim} is a factor "
        f"x{dims[g_dim]['one_mad_factor']:.3f} on that dimension and "
        f"x{dims[p_dim]['one_mad_factor']:.3f} on {p_dim}.",
        "",
        f"The reference is *pooled*, not a run of adjacent images: {ref['n']} chips "
        f"over {ref['months']} months, with {len(ref['months_missing'])} months empty"
        + (f" ({', '.join(ref['months_missing'])})" if ref["months_missing"] else "")
        + ". That is the population the origin stands for, and the moving-block draw "
        "resamples it in date order so a draw keeps neighbouring acquisitions "
        "together.",
        "",
        "## The months",
        "",
        f"| month | n | {g_dim} | {p_dim} | z({g_dim}) | z({p_dim}) | R | 95 % CI of "
        f"z({g_dim}) | 95 % CI of z({p_dim}) | thin |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- | --- |",
    ]
    for entry in months:
        g, p = entry["dims"][g_dim], entry["dims"][p_dim]
        lines.append(
            f"| {entry['month']} | {entry['n']} | {g['median']:.4g} | "
            f"{p['median']:.4g} | {_signed(g['z'])} | {_signed(p['z'])} | "
            f"{_plain(entry['R_behaviour']['value'])} | {_ci(g['ci95'])} | "
            f"{_ci(p['ci95'])} | {'✱' if entry['thin'] else '—'} |")
    lines += [
        "",
        f"`n` is the number of chips the month's median rests on, `✱` marks a month "
        f"below `SMALL_N = {SMALL_N}` — here {report['series']['thin_months']} of "
        f"{report['series']['months']}, 2024-09 — and `–` would be a cell the record "
        f"does not contain, never a zero. The interval is the *reference's* "
        f"uncertainty: {months[0]['dims'][g_dim]['draws']} circular moving-block "
        f"draws of the {ref['n']} reference chips with the reference median **and** "
        f"scale re-estimated inside each draw, the month's chips held fixed. It "
        f"answers \"how well is the origin known\", not \"how well is the month "
        f"known\".",
        "",
        "## The channel split",
        "",
        "The same months, re-standardised by an ASCENDING-only and a DESCENDING-only "
        "reference. The month's chips do not change — only the median and MAD of the "
        "origin — so any difference is the *reference composition*, entangled with "
        "the smaller per-month n of a single channel. It is the plate's declared "
        "sensitivity to the one window this figure had to choose:",
        "",
        "| channel | n_ref | n_series | ref median (log) | ref MAD (log) | max "
        "&#124;Δz&#124; | sign flips | cells beyond ±2 | R max |",
        "| --- | ---: | ---: | --- | --- | ---: | ---: | ---: | ---: |",
    ]
    for key in sorted(channels):
        block = channels[key]
        summary = block["summary"]
        dim_entry = block["dims"]
        lines.append(
            f"| {key} ({block['orbit_direction'].lower()}) | {block['n_ref']} | "
            f"{block['n_series']} | "
            f"{_plain(dim_entry[g_dim]['median_log'], 4)} / "
            f"{_plain(dim_entry[p_dim]['median_log'], 4)} | "
            f"{_plain(dim_entry[g_dim]['mad_log'], 4)} / "
            f"{_plain(dim_entry[p_dim]['mad_log'], 4)} | "
            f"{summary['max_abs_shift']:.2f} | {len(summary['sign_flips'])} | "
            f"{len(summary['outside_2mad'])} | {_plain(summary['R_max'])} |")
    lines += [
        "",
        "Δz is measured against the pooled-reference z of the same month; the two "
        "numbers per column are gamma2 / P. The flips are listed in the JSON under "
        "`channels.*.summary.sign_flips`. No month moves by more than "
        f"{max(b['summary']['max_abs_shift'] for b in channels.values()):.2f} MAD, and "
        "the cells that cross ±2 sit in months of 2-3 chips. The pooled ASC+DESC "
        "reference stays the headline; this table exists so a reader can see how much "
        "of \"no signal\" rests on that one choice.",
        "",
        "## What the series shows",
        "",
        *reading_lines(report),
        "",
        "## What it does not say",
        "",
        *caveat_lines(report),
        "",
        "## Reproduce",
        "",
        "```",
        "python3 code/carola/fig_carola_ocv_reconfiguration.py            # data + figure",
        "python3 code/carola/fig_carola_ocv_reconfiguration.py --selftest",
        "python3 code/carola/fig_carola_ocv_reconfiguration.py --verify",
        "bash figures/carola/fig_carola_ocv_reconfiguration_carola.sh",
        "```",
        "",
        f"The run reads `{_repo_rel(TABLE)}` (sha256 "
        f"`{report['input']['table_sha256'][:16]}…`) and `{_repo_rel(META)}`, and "
        f"writes `{_repo_rel(OUT_JSON)}`, `{_repo_rel(OUT_CSV)}`, "
        f"`{_repo_rel(OUT_JPG)}` and this file. No timestamp enters either output, so "
        f"a second run is byte-identical. `--selftest` re-derives the plate's "
        f"invariants from the table; `--verify` recomputes the report and diffs it "
        f"against the committed JSON (floats to 1e-9) before checking the CSV against "
        f"the same object.",
        "",
        "## The figure",
        "",
        "* **(a)** z(gamma2) by month: one dot per chip (light), the month's median "
        "(filled), the reference's ±2 MAD band (grey), the 95 % reference-uncertainty "
        "bar, the chip count `n` in the strip at the bottom, and the collapse of "
        f"{report['event']['collapse']} as the red cut. The dashed lines are the same "
        "months under the ASCENDING- and DESCENDING-only references (panel c's table "
        "is their source).",
        "* **(b)** z(P) on the same axis and the same reference — the two panels are "
        "read together, never as two independent verdicts.",
        "* **(c)** the plane: the origin is the reference state, the grey cloud is "
        f"the {ref['n']} reference chips, the path runs "
        f"{report['series']['first']} -> {report['series']['last']} and is coloured "
        "by quarter, the thin arrows are the quarter medians, and the terminal month "
        "is ringed.",
        "* **(d)** R = rms(z gamma2, z P) per month against the empirical 5-95 % band "
        "of the reference chips' own R.",
        "",
    ]
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# The figure
# ---------------------------------------------------------------------------
def _limits(months, dim, z_values):
    """Panel limits that always show the ±2 band and every interval."""
    values = [v for v in z_values if v is not None] + [-2.0, 2.0]
    for entry in months:
        ci = entry["dims"][dim]["ci95"]
        if ci:
            values.extend(ci)
    lo, hi = min(values), max(values)
    span = hi - lo if hi > lo else 1.0
    return lo - 0.17 * span, hi + 0.09 * span


def _month_panel(ax, report, months, dim, dates, z_values, colour, title):
    """One time panel: band, chips, reference CI, monthly median, channel split."""
    ax.axhspan(-2.0, 2.0, color="0.85", alpha=0.55, lw=0, zorder=0)
    ax.axhline(0.0, color="0.40", lw=0.8, zorder=1)
    for chip in report["chips"]:
        z = chip["z"][dim]
        if z is None:
            continue
        ax.plot(dt.date.fromisoformat(chip["date"]), z, marker="o", ms=1.7,
                mfc=colour, mec="none", alpha=0.32, zorder=2)
    # The same months under the single-orbit references — a sensitivity, drawn
    # thinner and dashed so it can never be mistaken for the headline series.
    for key, block in sorted(report["channels"].items()):
        points = [(month_marker(entry["month"]), entry["z"][dim])
                  for entry in block["months"] if entry["z"][dim] is not None]
        if points:
            ax.plot([p[0] for p in points], [p[1] for p in points],
                    color=CHANNEL_COLOUR.get(key, "0.45"), lw=0.65,
                    ls=(0, (2.4, 1.8)), alpha=0.9, zorder=2.4)
    for day, entry in zip(dates, months):
        ci = entry["dims"][dim]["ci95"]
        if ci:
            ax.plot([day, day], ci, color="0.45", lw=0.9, zorder=2.6)
    ax.plot(dates, z_values, color=colour, lw=0.8, alpha=0.9, zorder=3)
    ax.plot(dates, z_values, marker="o", ms=3.3, lw=0, color=colour, zorder=4)
    thin = [(day, z, entry) for day, z, entry in zip(dates, z_values, months)
            if entry["thin"]]
    if thin:
        ax.plot([t[0] for t in thin], [t[1] for t in thin], marker="o", ms=4.6,
                mfc="white", mec=colour, mew=1.0, lw=0, zorder=5)
    ax.axvline(dt.date.fromisoformat(EVENT), color=ALERT, lw=1.0, ls=(0, (4, 2)),
               zorder=6)
    lo, hi = _limits(months, dim, z_values)
    ax.set_ylim(lo, hi)
    ax.set_xlim(dt.date(2023, 12, 18), dt.date(2024, 10, 2))
    # The chip count of every month, as a strip inside the panel — the thin flag
    # never travels without its n.
    span = hi - lo
    for day, entry in zip(dates, months):
        ax.bar(day, entry["n"] * 0.008 * span, bottom=lo + 0.022 * span, width=7,
               color="0.72", lw=0, zorder=1.5)
        ax.text(day, lo + 0.003 * span, f"n={entry['n']}", fontsize=5.0,
                color="0.38", ha="center", va="bottom", zorder=6)
    ax.xaxis.set_major_locator(mdates.MonthLocator(interval=1))
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

    fig, axes = plt.subplots(2, 2, figsize=(13.0, 7.4),
                             gridspec_kw={"width_ratios": [1.62, 1.0]})
    ax_g, ax_p = axes[0][0], axes[1][0]
    ax_plane, ax_r = axes[0][1], axes[1][1]

    _month_panel(ax_g, report, months, g_dim, dates, z_g, g_colour,
                 f"a · z({ocv.DIM_DISPLAY[g_dim]}) by month — the 2022-09 … "
                 f"2023-12 reference is 0")
    _month_panel(ax_p, report, months, p_dim, dates, z_p, p_colour,
                 f"b · z({ocv.DIM_DISPLAY[p_dim]}) by month — same "
                 f"standardisation, same months")

    ax_g.text(0.985, 0.035,
              f"reference {ref['first']} – {ref['last']}\nn = {ref['n']} → 0 on both "
              f"axes (±2 MAD shaded)",
              transform=ax_g.transAxes, fontsize=5.6, va="bottom", ha="right",
              bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="0.80", lw=0.5,
                        alpha=0.92))
    handles = [
        Line2D([], [], color=g_colour, marker="o", ms=3.3, lw=0.8,
               label="monthly median"),
        Line2D([], [], color=g_colour, marker="o", ms=4.6, mfc="white", mew=1.0,
               lw=0, label=f"month with n < {SMALL_N}"),
        Line2D([], [], color=g_colour, marker="o", ms=1.7, alpha=0.45, lw=0,
               label="one chip"),
        Line2D([], [], color="0.45", lw=1.0, label="95 % CI of the reference"),
    ]
    for key in sorted(report["channels"]):
        handles.append(Line2D([], [], color=CHANNEL_COLOUR.get(key, "0.45"), lw=0.7,
                              ls=(0, (2.4, 1.8)),
                              label=f"origin = {key} ref only"))
    handles.append(Line2D([], [], color=ALERT, lw=1.0, ls=(0, (4, 2)),
                          label=f"collapse {EVENT}"))
    ax_g.legend(handles=handles, loc="upper right", frameon=True, framealpha=0.96,
                edgecolor="0.85", facecolor="white", fontsize=5.5, handlelength=1.5,
                borderaxespad=0.25, labelspacing=0.32)

    # (c) the plane — the origin is the reference state
    ref_g = [chip["z"][g_dim] for chip in report["reference_chips"]]
    ref_p = [chip["z"][p_dim] for chip in report["reference_chips"]]
    ax_plane.plot([0.0], [0.0], marker="+", ms=9, mew=1.1, color="black", zorder=7)
    ax_plane.axhline(0.0, color="0.40", lw=0.7, zorder=1)
    ax_plane.axvline(0.0, color="0.40", lw=0.7, zorder=1)
    ax_plane.plot(ref_g, ref_p, marker="o", ms=2.1, mfc="0.62", mec="none",
                  alpha=0.60, lw=0, zorder=2,
                  label=f"reference chip (n = {ref['n']})")
    ax_plane.plot(z_g, z_p, color="0.58", lw=0.7, zorder=3)
    offsets = {"Q1": (8, 4), "Q2": (8, -11), "Q3": (-11, -12)}
    for quarter in sorted({quarter_of(entry["month"]) for entry in months}):
        picked = [index for index, entry in enumerate(months)
                  if quarter_of(entry["month"]) == quarter]
        label = (f"{quarter} · {months[picked[0]]['month'][5:7]}–"
                 f"{months[picked[-1]]['month'][5:7]}")
        ax_plane.plot([z_g[i] for i in picked], [z_p[i] for i in picked], marker="o",
                      ms=3.4, lw=0, color=QUARTER_COLOUR[quarter], label=label,
                      zorder=4)
        mid = (float(np.median([z_g[i] for i in picked])),
               float(np.median([z_p[i] for i in picked])))
        ax_plane.annotate("", xy=mid, xytext=(0.0, 0.0),
                          arrowprops=dict(arrowstyle="-|>", lw=1.0,
                                          color=QUARTER_COLOUR[quarter], shrinkA=0,
                                          shrinkB=0, alpha=0.95), zorder=5)
        ax_plane.annotate(quarter, xy=mid, xytext=offsets.get(quarter, (8, 4)),
                          textcoords="offset points", fontsize=6.0,
                          color=QUARTER_COLOUR[quarter], zorder=6)
    # Quarter boundaries labelled on the path, so the plane can be read in time.
    for key in (months[0]["month"], "2024-04", "2024-07", months[-1]["month"]):
        entry = next((m for m in months if m["month"] == key), None)
        if entry is None or entry["dims"][g_dim]["z"] is None:
            continue
        ax_plane.annotate(key[5:7], xy=(entry["dims"][g_dim]["z"],
                                        entry["dims"][p_dim]["z"]),
                          xytext=(0, -11), textcoords="offset points",
                          ha="center", fontsize=5.4, color="0.30", zorder=6)
    ax_plane.plot([z_g[-1]], [z_p[-1]], marker="o", ms=9.5, mfc="none", mec=ALERT,
                  mew=1.2, zorder=8)
    ax_plane.annotate(f"{months[-1]['month']} (n = {months[-1]['n']})",
                      xy=(z_g[-1], z_p[-1]), xytext=(10, -10),
                      textcoords="offset points", fontsize=6.0, zorder=8,
                      arrowprops=dict(arrowstyle="-", lw=0.6, color="0.3"))
    radius = 1.06 * max(abs(v) for v in z_g + z_p + ref_g + ref_p)
    ax_plane.set_aspect("equal", adjustable="box")
    ax_plane.set_xlim(-radius, radius)
    ax_plane.set_ylim(-radius, radius)
    ax_plane.set_xlabel(f"z({ocv.DIM_DISPLAY[g_dim]})")
    ax_plane.set_ylabel(f"z({ocv.DIM_DISPLAY[p_dim]})")
    ax_plane.set_title(f"c · the plane — origin = reference state, grey = its "
                       f"{len(ref_g)} chips", loc="left")
    ax_plane.legend(loc="upper left", frameon=False, fontsize=5.6, ncol=1,
                    handlelength=0.9, labelspacing=0.3, borderaxespad=0.2)

    # (d) the block magnitude against the reference's own 5-95 % band
    ax_r.axhspan(band[0], band[1], color="0.85", alpha=0.55, lw=0, zorder=0)
    ax_r.axhline(0.0, color="0.40", lw=0.8, zorder=1)
    for chip in report["chips"]:
        if chip["R"] is None:
            continue
        ax_r.plot(dt.date.fromisoformat(chip["date"]), chip["R"], marker="o", ms=1.6,
                  mfc="#4a5568", mec="none", alpha=0.28, zorder=2)
    for day, entry in zip(dates, months):
        ci = entry["R_behaviour"]["ci95"]
        if ci:
            ax_r.plot([day, day], ci, color="0.45", lw=0.9, zorder=2.6)
    ax_r.plot(dates, r_values, color="#4a5568", lw=0.8, zorder=3)
    ax_r.plot(dates, r_values, marker="o", ms=3.3, lw=0, color="#4a5568", zorder=4)
    ax_r.plot([dates[-1]], [r_values[-1]], marker="o", ms=9.5, mfc="none", mec=ALERT,
              mew=1.2, zorder=8)
    ax_r.axvline(collapse, color=ALERT, lw=1.0, ls=(0, (4, 2)), zorder=6)
    ax_r.set_xlim(dt.date(2023, 12, 18), dt.date(2024, 10, 2))
    ax_r.set_ylim(0.0, 1.20 * max(max(r_values), band[1]))
    ax_r.xaxis.set_major_locator(mdates.MonthLocator(interval=2))
    ax_r.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m"))
    ax_r.set_ylabel("R_behaviour")
    ax_r.set_title("d · block magnitude R = rms(z gamma2, z P)", loc="left")
    ax_r.text(0.985, 0.94, f"reference 5–95 %: {band[0]:.2f}–{band[1]:.2f}",
              transform=ax_r.transAxes, fontsize=5.6, ha="right", va="top",
              color="0.30")

    fig.text(0.006, 0.985,
             "Carola Brücke — OCV reconfiguration against the 2022-09 … 2023-12 "
             "reference, month by month", fontsize=10.2, va="top", ha="left")
    fig.text(0.006, 0.947,
             f"z = (median log x_month − median log x_ref) / MAD(log x_ref); "
             f"n_ref = {ref['n']} ({ref['first']} … {ref['last']}), held out before "
             f"the series. Each month rests on "
             f"{min(e['n'] for e in months)}–{max(e['n'] for e in months)} chips, so "
             f"the monthly marker is the median of the chips that exist; the bar "
             f"drawn is the reference's own uncertainty, not the month's. The series "
             f"stops at the collapse — the "
             f"{report['audit']['post_outside_series']} later chips are the same "
             f"bridge after the event, counted and not plotted.",
             fontsize=6.4, va="top", ha="left", color="0.25")
    fig.text(0.006, 0.020,
             "dashed: origin rebuilt from one orbit only (sensitivity, not a second "
             "estimate) · 2024-09 is the only month below SMALL_N = "
             f"{SMALL_N} (n = {months[-1]['n']}) · no test is run",
             fontsize=5.8, va="bottom", ha="left", color="0.35")
    fig.subplots_adjust(top=0.878, bottom=0.085, left=0.055, right=0.988,
                        wspace=0.20, hspace=0.30)
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    fig.savefig(path, dpi=dpi or plt.rcParams["savefig.dpi"], facecolor="white")
    plt.close(fig)
    return path


# ---------------------------------------------------------------------------
# Outputs
# ---------------------------------------------------------------------------
CSV_FIELDS = ["Month", "n", "gamma2", "P", "z_gamma2", "z_P", "z_pooled_gamma2",
              "z_pooled_P", "R_behaviour", "ci95_gamma2", "ci95_P", "ci95_R", "thin"]


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

    check("table: 436 rows, Girder A only, both orbit directions",
          report["audit"]["rows"] == 436
          and report["input"]["segment"] == ["Girder A"]
          and report["input"]["orbit_rows"].get("ASC", 0) > 0
          and report["input"]["orbit_rows"].get("DESC", 0) > 0,
          f"rows = {report['audit']['rows']}, segments = "
          f"{report['input']['segment']}, orbit = {report['audit']['orbit']}")

    check("partition: 121 reference + 63 series + 252 post = 436",
          report["audit"]["reference"] == 121 and report["audit"]["series"] == 63
          and report["audit"]["post_outside_series"] == 252,
          report["audit"]["partition"])

    contiguous = month_range(months[0]["month"], months[-1]["month"])
    check("windows: reference ends 2023-12, series is 2024-01 … 2024-09, gap-free",
          report["reference"]["last"].startswith("2023-12")
          and report["series"]["first"] == "2024-01"
          and report["series"]["last"] == "2024-09"
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
    check("scale rule: MAD for both dims and both channels (all n >= SCALE_MIN_REF)",
          all(kind == "mad" for kind in kinds.values())
          and all(kind == "mad" for block in channel_kinds.values()
                  for kind in block.values()),
          f"{kinds}, channels {channel_kinds}")

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

    worst_r = 0.0
    for entry in months:
        value = math.sqrt(sum(entry["dims"][dim]["z"] ** 2 for dim in DIMS) / len(DIMS))
        worst_r = max(worst_r, abs(value - entry["R_behaviour"]["value"]))
    check("R is the rms of the two z values, every month", worst_r < 1e-12,
          f"max |R - rms| = {worst_r:.2e}")

    intervals = [ci for entry in months for dim in DIMS
                 for ci in (entry["dims"][dim]["ci95"], entry["R_behaviour"]["ci95"])
                 if ci]
    thin = [entry["month"] for entry in months if entry["thin"]]
    check("intervals are ordered and exactly 2024-09 is thin (n = 2)",
          all(ci[0] <= ci[1] for ci in intervals) and thin == ["2024-09"]
          and months[-1]["n"] == 2,
          f"{len(intervals)} interval(s), thin = {thin}, n = {months[-1]['n']}")

    twin = build_report(table, meta, boot=64)
    differences = _diff(report, twin)
    check("determinism: two builds are identical object for object",
          not differences, f"{len(differences)} difference(s)")

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
        description="Carola — OCV reconfiguration against the 2022-09 … 2023-12 "
                    "reference, month by month")
    parser.add_argument("--table", default=TABLE)
    parser.add_argument("--meta", default=META)
    parser.add_argument("--json", default=OUT_JSON)
    parser.add_argument("--csv", default=OUT_CSV)
    parser.add_argument("--md", default=None)
    parser.add_argument("--jpg", default=None)
    parser.add_argument("--figdir", default=None,
                        help="directory for the plate and its guide (default: "
                             "figures/carola)")
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
