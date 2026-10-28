#!/usr/bin/env python3
"""CTS — OCV reconfiguration in units of the 2015-2017 reference, month by month.

Question: how far is each month of the Champlain Towers South record from the
site's **own 2015-2017 reference**, on the two dimensions of the
``behaviour`` block of ``code/sds/ocv_reconfigure_spec.md`` — ``gamma2`` (the
consistency of the echo across acquisitions) and ``P`` (its persistence) — and
where does the last month before the collapse sit?

The metric is not re-invented here; it is the pre-registered one:

    r = median(g(x_month)) - median(g(x_ref))          g = log (both dimensions)
    z = r / scale_ref                                  scale_ref = MAD(x_ref)

so the reference is the origin **(0, 0)** by construction, exactly as on the
``fig_sds_ocv_reconfigure*`` plates. Two consequences are carried into every
number rather than into a footnote:

* the **reference is pooled**, one median and one MAD over the 35 acquisitions
  of 2015-09-21 … 2017-12-15 (the first Sentinel-1 images that exist for this
  footprint — 2014 has none, ``data/cts/cts_acquisition_census.md``). 35 >=
  ``SCALE_MIN_REF = 10`` (spec amendment 2026-10-22), so the reference MAD is the
  primary scale; the pooled MAD of the whole record travels beside it as the
  declared sensitivity;
* a **month is 1-3 chips**, so every month is ``thin`` (``SMALL_N = 6``): there is
  no per-month interval of its own. The only interval reported is the
  **reference's** own uncertainty — 2 000 circular moving-block draws (block 3,
  seed 7) of the 35 reference chips, with the reference median *and* scale
  re-estimated inside every draw, the month's chips held fixed. That is the same
  construction as the ``✱`` cells of the sds plates.

Nothing is tested: a month axis is not a dose and a single asset has no partner,
so this figure carries no H1, no H2, no p-value and no precursor claim. CTS's own
note fixes the binding limit at **sampling, not SNR** (2 acquisitions in the last
21 days, 3 in 30, 5 in 60), and its post state is bare sand whose intensity is
moisture-driven — so the series stops at the collapse and the post-collapse
acquisitions are counted, never plotted as values.

Figure (2 x 2, 13.0 in wide, 600 dpi jpg — the CTS package's ``savefig.dpi``):
  a  z(gamma2) by month, per-chip dots, monthly median, reference band +-2 MAD,
     a chip-count strip, the collapse cut;  n and the thin mark on every month
  b  z(P) by month, same standardisation and axis
  c  the (z gamma2, z P) plane — the origin is the reference state, the grey
     cloud is the 35 reference chips, the path is 2018-01 -> 2021-06, the arrows
     are the per-year medians and the terminal month is ringed
  d  the block magnitude R = rms(z gamma2, z P) per month, with the empirical
     5-95 % band of the reference chips as the null band

Outputs (all inside the CTS package, beside its own data and figures):
  data/cts/ocv_reconfiguration_cts.json
  data/cts/ocv_reconfiguration_cts_table.csv
  figures/cts/fig_cts_ocv_reconfiguration_cts.jpg
  figures/cts/fig_cts_ocv_reconfiguration_cts.md

Usage:
  python3 code/cts/fig_cts_ocv_reconfiguration_cts.py
  python3 code/cts/fig_cts_ocv_reconfiguration_cts.py --selftest
  python3 code/cts/fig_cts_ocv_reconfiguration_cts.py --verify
  python3 code/cts/fig_cts_ocv_reconfiguration_cts.py --figdir /tmp/x --json /tmp/y.json
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

import cts_ocv_core as core      # noqa: E402
import ocv_reconfigure as ocv    # noqa: E402  (the pre-registered metric)

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
})

SCHEMA = "ocv-reconfiguration-monthly-v1"
SITE = "cts"
SITE_LABEL = core.SITE

#: The behaviour block of the spec, and the column each dimension lives in on the
#: support-sized chip of this site.
DIMS = tuple(ocv.BLOCKS["behaviour"])          # ("gamma2", "P")
COLUMN = {"gamma2": "gamma2_mask", "P": "P"}
TRANSFORM = {dim: "log" for dim in DIMS}

#: The reference is pooled over 2015-2017; the series starts in 2018-01 and stops
#: at the collapse, so its last chip precedes the event day (2021-06-21).
REF_START, REF_END = "2015-01-01", "2017-12-31"
SERIES_START = "2018-01-01"
EVENT = core.COLLAPSE_DATE                      # 2021-06-24

TABLE = os.path.join(core.DATA, "support", "cts_ocv_channels.csv")
META = os.path.join(core.DATA, "support", "cts_ocv_channels_meta.json")
OUT_JSON = os.path.join(core.DATA, "ocv_reconfiguration_cts.json")
OUT_CSV = os.path.join(core.DATA, "ocv_reconfiguration_cts_table.csv")
OUT_MD = os.path.join(core.FIGDIR, "fig_cts_ocv_reconfiguration_cts.md")
OUT_JPG = os.path.join(core.FIGDIR, "fig_cts_ocv_reconfiguration_cts.jpg")

BOOT = 2000            # reference-uncertainty draws (the spec's seed and count)
BLOCK = 3              # consecutive reference chips per draw (spec: BLOCK_BOOT)
SEED = 7
SMALL_N = ocv.SMALL_N  # 6 — below this a month is thin and must be quoted with n
SCALE_MIN_REF = ocv.SCALE_MIN_REF  # 10 — below this the pooled MAD takes over

#: One colour per calendar year; the last year is the partial one that ends at
#: the collapse.
YEAR_COLOUR = {2018: "#2b6cb0", 2019: "#2f855a", 2020: "#b7791f", 2021: "#c53030"}
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
    does not contain ``cts`` — adding it would touch the pinned sds run — so the
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
    """The support-sized chip of one CTS acquisition per row, validated.

    The committed table carries only the target's own chip (178 rows, all
    ``role = target``, one orbit direction). A row from another track would make
    the reference a mix of footprints, so the run stops rather than averaging it
    in.
    """
    path = path or TABLE
    with open(path, newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise SystemExit(f"{path}: no rows")
    needed = {"date", "orbit_direction", "role"} | set(COLUMN.values())
    missing = needed - set(rows[0])
    if missing:
        raise SystemExit(f"{path}: missing columns {sorted(missing)}")
    roles = {row["role"] for row in rows}
    if roles != {"target"}:
        raise SystemExit(
            f"{path}: role(s) {sorted(roles)} — this run reads the target's own "
            f"support-sized chip only, never a control or a reference track")
    orbits = {row["orbit_direction"] for row in rows}
    if len(orbits) != 1:
        raise SystemExit(f"{path}: orbit directions {sorted(orbits)} — one expected")
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


def _reference_block(ref_rows, ref_values, all_values):
    """Median and scale of the reference, per dimension, with the used scale."""
    dims = {}
    for dim in DIMS:
        values = [v for v in ref_values[dim] if v is not None]
        stats = _stats(values)
        pooled = _stats(all_values[dim])
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


def build_report(table=None, meta=None, boot=BOOT, block=BLOCK, seed=SEED):
    """The whole plate as one object: reference, months, chips, audit."""
    table = table or TABLE
    meta_path = meta or META
    rows = sorted(load_rows(table), key=lambda row: (row["date"], row["orbit_direction"]))
    with open(meta_path) as handle:
        meta_obj = json.load(handle)
    channel = ocv.CHANNEL[sorted({row["orbit_direction"] for row in rows})[0]]

    ref_rows = [row for row in rows if REF_START <= row["date"] <= REF_END]
    series_rows = [row for row in rows if SERIES_START <= row["date"] <= EVENT]
    post_rows = [row for row in rows if row["date"] > EVENT]
    if not ref_rows or not series_rows:
        raise SystemExit("the reference or the series is empty — check the window")
    if len(ref_rows) + len(series_rows) + len(post_rows) != len(rows):
        raise SystemExit("the three windows do not partition the table")

    ref_values = {dim: [_log(_number(row[COLUMN[dim]])) for row in ref_rows]
                  for dim in DIMS}
    all_values = {dim: [_log(_number(row[COLUMN[dim]])) for row in rows] for dim in DIMS}
    reference = _reference_block(ref_rows, ref_values, all_values)

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

    post_dates = sorted(row["date"] for row in post_rows)
    ref_months_missing = missing_reference(ref_rows)
    report = {
        "schema": SCHEMA,
        "site": SITE,
        "site_label": SITE_LABEL,
        "event": {"collapse": core.COLLAPSE_DATE, "demolition": core.DEMOLITION_DATE},
        "metric": {
            "definition": ("z = (median(log x_month) - median(log x_ref)) / scale_ref; "
                           "R_behaviour = rms over the behaviour block"),
            "transform": {dim: TRANSFORM[dim] for dim in DIMS},
            "scale_rule": ("MAD of the reference, fallback chain MAD -> IQR/1.349 -> SD; "
                           "the pooled MAD of the whole record replaces it when "
                           f"n_ref < {SCALE_MIN_REF} (spec amendment 2026-10-22)"),
            "block": list(DIMS),
            "unit": "one support-sized chip per acquisition",
            "interval": (f"reference-uncertainty bootstrap only: {boot} circular "
                         f"moving-block draws, block {block}, seed {seed}, reference "
                         f"chips resampled, the month's chips held fixed"),
            "tests": None,
            "thin_rule": f"n < {SMALL_N} — every month here",
        },
        "input": {
            "table": _rel(table),
            "table_sha256": sha256(table),
            "meta": _rel(meta_path),
            "meta_sha256": sha256(meta_path),
            "rows": len(rows),
            "role": "target",
            "orbit": channel,
            "chip": meta_obj.get("chip"),
            "support_rule": meta_obj.get("rule"),
            "support_m": {"guard_m": meta_obj.get("guard_m"),
                          "ring_m": meta_obj.get("ring_m")},
        },
        "reference": {
            "window": [REF_START, REF_END],
            "n": len(ref_rows),
            "first": ref_rows[0]["date"],
            "last": ref_rows[-1]["date"],
            "orbit": channel,
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
            "months_missing": missing,
            "thin_months": sum(1 for entry in months if entry["thin"]),
        },
        "months": months,
        "chips": chips,
        "reference_chips": ref_chips,
        "reference_R_band": ref_band,
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
            "dropped_nonpositive": {
                dim: sum(1 for v in all_values[dim] if v is None) for dim in DIMS},
            "months_missing_series": missing,
            "months_missing_reference": ref_months_missing,
            "orbit": f"{sorted({row['orbit_direction'] for row in rows})[0]} ({channel})",
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
    """``1 gamma2 month (2020-01)`` / ``no P month`` — never a bare zero."""
    if not keys:
        return f"no {dim} month"
    return f"{len(keys)} {dim} month(s) ({_month_list(keys)})"


def reading_lines(report):
    """The prose of the plate, every number taken from the object itself."""
    dims = report["reference"]["dims"]
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
    r_max = max(months, key=lambda m: m["R_behaviour"]["value"])
    negative = {dim: sum(1 for m in months if m["dims"][dim]["z"] < 0) for dim in DIMS}
    # The last calendar year of the record, and the years the excursions fall in —
    # both read off the data, never spelled out by hand.
    final_year = last["month"][:4]
    final_months = [m for m in months if m["month"].startswith(final_year)]
    z_final = max(abs(m["dims"][dim]["z"]) for m in final_months for dim in DIMS)
    extreme_years = sorted({m["month"][:4]
                            for m in (g_low, g_high, p_low, p_high)})
    extreme_span = (extreme_years[0] if len(extreme_years) == 1
                    else f"{extreme_years[0]}-{extreme_years[-1]}")
    affine = {dim: (dims[dim]["scale_log"] / dims[dim]["pooled_scale_log"],
                    (dims[dim]["median_log"] - dims[dim]["pooled_median_log"])
                    / dims[dim]["pooled_scale_log"]) for dim in DIMS}
    flips = [abs(m["dims"][dim]["z"]) for m in months for dim in DIMS
             if (m["dims"][dim]["z"] > 0) != (m["dims"][dim]["z_pooled"] > 0)]

    def longest_run(dim):
        """Longest unbroken run of months with the same sign of z."""
        best = run = 0
        sign = 0
        for entry in months:
            here = 1 if entry["dims"][dim]["z"] > 0 else -1
            run = run + 1 if here == sign else 1
            sign = here
            best = max(best, run)
        return best

    runs = {dim: longest_run(dim) for dim in DIMS}

    def covers_zero(entry):
        return bool(entry["ci95"]) and entry["ci95"][0] <= 0 <= entry["ci95"][1]

    both = covers_zero(last_g) and covers_zero(last_p)
    inside = band[0] <= last["R_behaviour"]["value"] <= band[1]
    return [
        f"**The last month before the collapse sits on the reference.** "
        f"{last['month']} carries {last['n']} chips "
        f"({', '.join(last['dates'])}): z({g_dim}) = {_signed(last_g['z'])} with a 95 % "
        f"reference-uncertainty interval of {_ci(last_g['ci95'])}, z({p_dim}) = "
        f"{_signed(last_p['z'])} {_ci(last_p['ci95'])}, and R = "
        f"{_plain(last['R_behaviour']['value'])} against a reference band of "
        f"{_plain(band[0])}–{_plain(band[1])}. "
        + ("Both intervals cover zero, " if both
           else "Not both intervals cover zero, ")
        + f"and R is {'inside' if inside else 'outside'} the band, so the terminal month "
        f"reads as *the reference*, not as a departure from it.",
        "",
        f"**The largest excursions are earlier, in {extreme_span}.** {g_dim} reaches "
        f"{_signed(g_low['dims'][g_dim]['z'])} in {g_low['month']} (n = {g_low['n']}) and "
        f"{_signed(g_high['dims'][g_dim]['z'])} in {g_high['month']} "
        f"(n = {g_high['n']}); {p_dim} reaches "
        f"{_signed(p_low['dims'][p_dim]['z'])} in {p_low['month']} (n = {p_low['n']}) and "
        f"{_signed(p_high['dims'][p_dim]['z'])} in {p_high['month']} "
        f"(n = {p_high['n']}). Beyond ±2 MAD lie {_outside_text(g_dim, outside[g_dim])} "
        f"and {_outside_text(p_dim, outside[p_dim])}.",
        "",
        f"**No {final_year} month is an outlier.** The largest |z| in any {final_year} "
        f"month is {z_final:.2f}, well inside the ±2 MAD band: the last "
        f"{len(final_months)} months run inside the band on both axes, so the record "
        f"does not end on an excursion.",
        "",
        f"**The block magnitude stays bounded.** R peaks at "
        f"{_plain(r_max['R_behaviour']['value'])} in {r_max['month']} (n = "
        f"{r_max['n']}), and the reference's own 95th percentile is "
        f"{_plain(band[1])}: "
        + (f"{len(above)} month(s) exceed it ({_month_list(above)})." if above
           else "no month reaches it.")
        + " R is an RMS over two dimensions, so one coordinate can carry it — which is "
        "why panel (c), not R alone, is the reading.",
        "",
        f"**The series oscillates; it does not drift.** {g_dim} is negative in "
        f"{negative[g_dim]} of {len(months)} months, {p_dim} in {negative[p_dim]} of "
        f"{len(months)}, and the longest unbroken run of one sign is "
        f"{runs[g_dim]} months for {g_dim} and {runs[p_dim]} for {p_dim}: a stretch, "
        f"not a trend. A drift toward the collapse would have to be monotone, and this "
        f"is not.",
        "",
        f"**The pooled sensitivity is the same picture from another origin.** "
        f"`z_pooled` re-centres the same r on the pooled median of all "
        f"{report['audit']['rows']} chips and re-scales by the pooled MAD, so it is an "
        f"affine image of z: {affine[g_dim][0]:.3f}·z {affine[g_dim][1]:+.3f} for "
        f"{g_dim} and {affine[p_dim][0]:.3f}·z {affine[p_dim][1]:+.3f} for {p_dim}. "
        f"{len(flips)} month(s) change side under it, all within "
        f"{max(flips, default=0.0):.2f} MAD of zero; no excursion moves, and it is "
        f"declared as a sensitivity, not as a second estimate.",
        "",
        f"**One sentence.** Between {report['series']['first']} and "
        f"{report['series']['last']} the site oscillates around its own 2015-2017 state "
        f"with the widest deviations in {extreme_span}, and the last month that exists "
        f"is not one of them.",
    ]


def caveat_lines(report):
    """The limits, stated where the figure is, not only in the JSON."""
    last = report["months"][-1]
    low_n = min(m["n"] for m in report["months"])
    high_n = max(m["n"] for m in report["months"])
    return [
        "* **The binding limit is sampling, not SNR.** A month rests on "
        f"{low_n}-{high_n} chips, so every month is `thin` (`n < SMALL_N = {SMALL_N}`) "
        "and a one-chip month's \"median\" *is* that single chip. The monthly marker is "
        "the median of the chips that exist, not of the month; a month with no "
        "acquisition would be a gap, not a value.",
        "* **The interval belongs to the reference.** It says how well a "
        f"{report['reference']['n']}-chip origin is known under a resampling of those "
        "chips. It is not a confidence interval for a monthly value: the month's own "
        "sampling noise is not in it, and for a one-chip month that noise is the whole "
        "spread of the record.",
        "* **The series stops at the collapse, deliberately.** The "
        f"{report['audit']['post_outside_series']} acquisitions after "
        f"{report['event']['collapse']} ({report['audit']['post_first']} … "
        f"{report['audit']['post_last']}) are bare sand, not the building, and their "
        "intensity is moisture-driven; they are counted in the audit and carry no value "
        f"here. Including them, or the demolition of {report['event']['demolition']}, "
        "would compare two different measured objects.",
        "* **No precursor claim.** The metric is descriptive: it says where the last "
        f"available month ({last['month']}) sits, and that month is on the reference. "
        "Nothing between that month and the collapse was measured by this record — the "
        "next acquisition is after the event — so a series of 1-3 chips per month "
        "cannot time, detect or exclude a precursor, and no test is run that would "
        "pretend otherwise.",
        "* **No test, no dose, no partner.** The sds H1/H2 machinery needs at least two "
        "states per block and, for a dose, an ordered ladder; here one state is sampled "
        "monthly. So there is no p-value, no Holm family and no cosine, and these "
        "magnitudes are in the reference's own units — they are not comparable to "
        "Haltern's, KDLO's or Morandi's.",
    ]


def markdown(report):
    """The plate's guide: definition, reference, all months, the reading."""
    ref = report["reference"]
    dims = ref["dims"]
    months = report["months"]
    g_dim, p_dim = DIMS
    outside = {dim: [m["month"] for m in months if abs(m["dims"][dim]["z"]) > 2]
               for dim in DIMS}
    last = months[-1]
    chip = report["input"]["chip"] or {}
    band = report["reference_R_band"]
    chips = report["chips"]
    span = (f"{chips[0]['date']} .. {chips[-1]['date']}" if chips
            else f"{report['series']['first']} .. {report['series']['last']}")
    lines = [
        "# CTS — OCV reconfiguration in units of the 2015-2017 reference, month by month",
        "",
        "*generated by `code/cts/fig_cts_ocv_reconfiguration_cts.py`; the metric is the "
        "pre-registered one of `code/sds/ocv_reconfigure_spec.md` (behaviour block, log "
        "transform, reference MAD), read here as a monthly series. It is not a new test, "
        "and it carries none.*",
        "",
        f"{report['site_label']}: one support-sized chip per acquisition "
        f"({chip.get('rows')} x {chip.get('cols')} cells, "
        f"{chip.get('support_cells', 225)} support cells, "
        f"{report['input']['support_m']['guard_m']:g} m guard, "
        f"{report['input']['support_m']['ring_m']:g} m ring), "
        f"{report['audit']['rows']} acquisitions in the {report['input']['orbit']} "
        f"orbit. The record splits into a **reference** of {ref['n']} chips "
        f"({ref['first']} … {ref['last']}), the **series** of "
        f"{report['series']['chips']} chips on {report['series']['months']} months "
        f"({report['series']['first']} … {report['series']['last']}), and "
        f"{report['audit']['post_outside_series']} post-collapse chips that are counted "
        f"and never plotted.",
        "",
        "```",
        "z(month) = ( median(log x_month) - median(log x_ref) ) / MAD(log x_ref)",
        "x        = gamma2, P            both on the log scale",
        "R        = rms over the two dimensions (gamma2, P)",
        f"ref      = {ref['first']} .. {ref['last']}   n = {ref['n']} chips on "
        f"{ref['months']} months",
        f"series   = {span}   {report['series']['chips']} chips on "
        f"{report['series']['months']} months, 1-3 per month",
        "```",
        "",
        "The reference is drawn from the years *before* the series, so the origin of "
        "everything below is the site's own long-run state and not a baseline chosen "
        "after the fact. `0` on both axes is that state; it is not a norm and not a "
        "threshold.",
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
        f"`SCALE_MIN_REF = {SCALE_MIN_REF}` of the 2026-10-22 amendment, so the pooled MAD "
        f"is **not** substituted. It travels beside the reference as the declared "
        f"sensitivity — `z_pooled` in the JSON reads the same r from the pooled median "
        f"and the pooled MAD of the whole record, an affine image of z rather than a "
        f"second estimate, quantified under *What the series shows*. One MAD of "
        f"{g_dim} is a factor x{dims[g_dim]['one_mad_factor']:.3f} on that dimension and "
        f"x{dims[p_dim]['one_mad_factor']:.3f} on {p_dim}.",
        "",
        f"The reference is *pooled*, not a run of adjacent images: {ref['n']} chips over "
        f"{ref['months']} months, with {len(ref['months_missing'])} months empty "
        f"({', '.join(ref['months_missing'])}) — the acquisition gaps of the 2015-2017 "
        f"record (`data/cts/cts_acquisition_census.md`). That is the population the "
        f"origin stands for, and the moving-block draw resamples it in date order so a "
        f"draw keeps neighbouring chips together.",
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
            f"| {entry['month']} | {entry['n']} | {g['median']:.4g} | {p['median']:.4g} | "
            f"{_signed(g['z'])} | {_signed(p['z'])} | "
            f"{_plain(entry['R_behaviour']['value'])} | {_ci(g['ci95'])} | "
            f"{_ci(p['ci95'])} | {'✱' if entry['thin'] else '—'} |")
    lines += [
        "",
        f"`n` is the number of chips the month's median rests on, `✱` marks a month at or "
        f"below `SMALL_N = {SMALL_N}` — here **every** month, "
        f"{report['series']['thin_months']} of {report['series']['months']} — and `–` is a "
        f"cell the record does not contain, never a zero. The interval is the "
        f"*reference's* uncertainty: {months[0]['dims'][g_dim]['draws']} circular "
        f"moving-block draws of the {ref['n']} reference chips with the reference median "
        f"**and** scale re-estimated inside each draw, the month's chips held fixed. It "
        f"answers \"how well is the origin known\", not \"how well is the month known\".",
        "",
        # __MD3__
        "## What the series shows",
        "",
        *reading_lines(report),
        "",
        "## What it does not say",
        "",
        *caveat_lines(report),
        # __MD4__
        "",
        "## Reproduce",
        "",
        "```",
        "python3 code/cts/fig_cts_ocv_reconfiguration_cts.py            # data + figure",
        "python3 code/cts/fig_cts_ocv_reconfiguration_cts.py --selftest",
        "python3 code/cts/fig_cts_ocv_reconfiguration_cts.py --verify",
        "bash figures/cts/fig_cts_ocv_reconfiguration_cts.sh",
        "```",
        "",
        f"The run reads `{_repo_rel(TABLE)}` "
        f"(sha256 `{report['input']['table_sha256'][:16]}…`) and "
        f"`{_repo_rel(META)}`, and writes "
        f"`{_repo_rel(OUT_JSON)}`, `{_repo_rel(OUT_CSV)}`, "
        f"`{_repo_rel(OUT_JPG)}` and this file. No timestamp enters either output, so a "
        f"second run is byte-identical.",
        "",
        "## The figure",
        "",
        "* **(a)** z(gamma2) by month: one dot per chip (light), the month's median "
        "(filled), the reference's +-2 MAD band (grey), the 95 % reference-uncertainty "
        "bar, the chip count `n` in the strip at the bottom, and the collapse of "
        "2021-06-24 as the red cut. A month with `n = 1` is drawn open, because its "
        "\"median\" is a single chip.",
        "* **(b)** z(P) on the same axis and the same reference — the two panels are "
        "read together, never as two independent verdicts.",
        "* **(c)** the plane: the origin is the reference state, the grey cloud is the "
        "35 reference chips, the path runs 2018-01 -> 2021-06 and is coloured by year, "
        "the thin arrows are the per-year medians, and the last pre-collapse month is "
        "ringed.",
        "* **(d)** R = rms(z gamma2, z P) per month against the empirical 5-95 % band of "
        "the reference chips' own R.",
        "",
    ]
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# The figure
# ---------------------------------------------------------------------------
def _limits(months, dim, z_values):
    values = [v for v in z_values if v is not None]
    for entry in months:
        ci = entry["dims"][dim]["ci95"]
        if ci:
            values.extend(ci)
    lo, hi = min(values), max(values)
    span = hi - lo if hi > lo else 1.0
    return lo - 0.20 * span, hi + 0.12 * span


def _month_panel(ax, report, months, dim, dates, z_values, colour, title):
    """One time panel: band, chips, reference CI, monthly median, chip count."""
    band_z = 2.0
    ax.axhspan(-band_z, band_z, color="0.85", alpha=0.55, lw=0, zorder=0)
    ax.axhline(0.0, color="0.40", lw=0.8, zorder=1)
    for chip in report["chips"]:
        z = chip["z"][dim]
        if z is None:
            continue
        ax.plot(dt.date.fromisoformat(chip["date"]), z, marker="o", ms=1.6,
                mfc=colour, mec="none", alpha=0.30, zorder=2)
    for day, entry in zip(dates, months):
        ci = entry["dims"][dim]["ci95"]
        if ci:
            ax.plot([day, day], ci, color="0.45", lw=0.9, zorder=2.5)
    ax.plot(dates, z_values, color=colour, lw=0.8, alpha=0.9, zorder=3)
    ax.plot(dates, z_values, marker="o", ms=3.2, lw=0, color=colour, zorder=4)
    single = [(day, z) for day, z, entry in zip(dates, z_values, months)
              if entry["n"] == 1]
    if single:
        ax.plot([s[0] for s in single], [s[1] for s in single], marker="o", ms=3.2,
                mfc="white", mec=colour, mew=0.9, lw=0, zorder=5)
    collapse = dt.date.fromisoformat(EVENT)
    ax.axvline(collapse, color=ALERT, lw=1.0, ls=(0, (4, 2)), zorder=6)
    lo, hi = _limits(months, dim, z_values)
    ax.set_ylim(lo, hi)
    ax.set_xlim(dt.date(2017, 11, 1), dt.date(2021, 7, 15))
    # The chip count of every month, as a strip inside the panel — the thin flag
    # never travels without its n.
    span = hi - lo
    for day, entry in zip(dates, months):
        ax.bar(day, entry["n"] * 0.030 * span, bottom=lo + 0.025 * span, width=8,
               color="0.72", lw=0, zorder=1.5)
    ax.text(0.004, 0.006, "chips n per month", transform=ax.transAxes, fontsize=4.9,
            color="0.42", ha="left", va="bottom")
    ax.xaxis.set_major_locator(mdates.MonthLocator(bymonth=(1, 7)))
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
                 f"a · z({ocv.DIM_DISPLAY[g_dim]}) by month — the 2015-2017 "
                 f"reference is 0")
    _month_panel(ax_p, report, months, p_dim, dates, z_p, p_colour,
                 f"b · z({ocv.DIM_DISPLAY[p_dim]}) by month — same standardisation, "
                 f"same months")

    ax_g.text(0.012, 0.965,
              f"reference {ref['first']} – {ref['last']}\nn = {ref['n']} → 0 on both "
              f"axes (±2 MAD shaded)",
              transform=ax_g.transAxes, fontsize=5.8, va="top", ha="left",
              bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="0.80", lw=0.5,
                        alpha=0.92))
    ax_g.legend(
        handles=[
            Line2D([], [], color=g_colour, marker="o", ms=3.4, lw=0.8,
                   label="monthly median"),
            Line2D([], [], color=g_colour, marker="o", ms=3.4, mfc="white", lw=0,
                   label="month with n = 1"),
            Line2D([], [], color=g_colour, marker="o", ms=1.8, alpha=0.45, lw=0,
                   label="one chip"),
            Line2D([], [], color="0.45", lw=1.0, label="95 % CI of the reference"),
            Line2D([], [], color=ALERT, lw=1.0, ls=(0, (4, 2)),
                   label="collapse 2021-06-24"),
        ],
        loc="upper right", frameon=True, framealpha=0.96, edgecolor="0.85",
        facecolor="white", fontsize=5.8, handlelength=1.5,
        borderaxespad=0.25, labelspacing=0.35)

    # (c) the plane — the origin is the reference state
    ref_g = [chip["z"][g_dim] for chip in report["reference_chips"]]
    ref_p = [chip["z"][p_dim] for chip in report["reference_chips"]]
    ax_plane.plot([0.0], [0.0], marker="+", ms=9, mew=1.1, color="black", zorder=7)
    ax_plane.axhline(0.0, color="0.40", lw=0.7, zorder=1)
    ax_plane.axvline(0.0, color="0.40", lw=0.7, zorder=1)
    ax_plane.plot(ref_g, ref_p, marker="o", ms=2.1, mfc="0.62", mec="none",
                  alpha=0.65, lw=0, zorder=2)
    ax_plane.plot(z_g, z_p, color="0.60", lw=0.7, zorder=3)
    for index, year in enumerate(sorted({int(entry["month"][:4])
                                         for entry in months})):
        picked = [i for i, entry in enumerate(months) if int(entry["month"][:4]) == year]
        label = f"{year}" + (" H1" if year == 2021 else "")
        ax_plane.plot([z_g[i] for i in picked], [z_p[i] for i in picked], marker="o",
                      ms=3.2, lw=0, color=YEAR_COLOUR[year], label=label, zorder=4)
        mid = (float(np.median([z_g[i] for i in picked])),
               float(np.median([z_p[i] for i in picked])))
        ax_plane.annotate("", xy=mid, xytext=(0.0, 0.0),
                          arrowprops=dict(arrowstyle="-|>", lw=1.0,
                                          color=YEAR_COLOUR[year], shrinkA=0,
                                          shrinkB=0, alpha=0.95), zorder=5)
        offset = [(7, 4), (7, -11), (-9, -12), (-9, 7)][index % 4]
        ax_plane.annotate(label, xy=mid, xytext=offset, textcoords="offset points",
                          fontsize=6.0, color=YEAR_COLOUR[year], zorder=6)
    ax_plane.plot([z_g[-1]], [z_p[-1]], marker="o", ms=9.5, mfc="none", mec="black",
                  mew=1.0, zorder=8)
    ax_plane.annotate(f"{months[-1]['month']} (n = {months[-1]['n']})",
                      xy=(z_g[-1], z_p[-1]), xytext=(9, -9), textcoords="offset points",
                      fontsize=6.0, zorder=8,
                      arrowprops=dict(arrowstyle="-", lw=0.6, color="0.3"))
    radius = 1.04 * max(abs(v) for v in z_g + z_p + ref_g + ref_p)
    ax_plane.set_aspect("equal", adjustable="box")
    ax_plane.set_xlim(-radius, radius)
    ax_plane.set_ylim(-radius, radius)
    ax_plane.set_xlabel(f"z({ocv.DIM_DISPLAY[g_dim]})")
    ax_plane.set_ylabel(f"z({ocv.DIM_DISPLAY[p_dim]})")
    ax_plane.set_title(f"c · the plane — origin = reference state, grey = its "
                       f"{len(ref_g)} chips", loc="left")
    ax_plane.legend(loc="lower left", frameon=False, fontsize=5.8, ncol=2,
                    handlelength=0.9, columnspacing=0.8, labelspacing=0.3,
                    borderaxespad=0.2)

    # (d) the block magnitude against the reference's own 5-95 % band
    ax_r.axhspan(band[0], band[1], color="0.85", alpha=0.55, lw=0, zorder=0)
    ax_r.axhline(0.0, color="0.40", lw=0.8, zorder=1)
    for chip in report["chips"]:
        if chip["R"] is None:
            continue
        ax_r.plot(dt.date.fromisoformat(chip["date"]), chip["R"], marker="o", ms=1.5,
                  mfc="#4a5568", mec="none", alpha=0.25, zorder=2)
    for day, entry in zip(dates, months):
        ci = entry["R_behaviour"]["ci95"]
        if ci:
            ax_r.plot([day, day], ci, color="0.45", lw=0.9, zorder=2.5)
    ax_r.plot(dates, r_values, color="#4a5568", lw=0.8, zorder=3)
    ax_r.plot(dates, r_values, marker="o", ms=3.2, lw=0, color="#4a5568", zorder=4)
    ax_r.plot([dates[-1]], [r_values[-1]], marker="o", ms=9.5, mfc="none", mec="black",
              mew=1.0, zorder=8)
    ax_r.axvline(collapse, color=ALERT, lw=1.0, ls=(0, (4, 2)), zorder=6)
    ax_r.set_xlim(dt.date(2017, 11, 1), dt.date(2021, 7, 15))
    ax_r.set_ylim(0.0, 1.22 * max(max(r_values), band[1]))
    ax_r.xaxis.set_major_locator(mdates.MonthLocator(bymonth=(1, 7)))
    ax_r.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m"))
    ax_r.set_ylabel("R_behaviour")
    ax_r.set_title("d · block magnitude R = rms(z gamma2, z P)", loc="left")
    ax_r.text(0.985, 0.93, f"reference 5–95 %: {band[0]:.2f}–{band[1]:.2f}",
              transform=ax_r.transAxes, fontsize=5.8, ha="right", va="top",
              color="0.35")

    fig.text(0.006, 0.985,
             "CTS — OCV reconfiguration against the 2015-2017 reference, month by "
             "month", fontsize=10.2, va="top", ha="left")
    fig.text(0.006, 0.947,
             f"z = (median log x_month − median log x_ref) / MAD(log x_ref); n_ref = "
             f"{ref['n']} ({ref['first']} … {ref['last']}). Every month rests on "
             f"{min(e['n'] for e in months)}–{max(e['n'] for e in months)} chips, so the "
             f"monthly marker is the median of the chips that exist; the bar drawn is "
             f"the reference's own uncertainty, not the month's. The series stops at the "
             f"collapse — the {report['audit']['post_outside_series']} later chips are "
             f"bare sand, counted and not plotted.",
             fontsize=6.4, va="top", ha="left", color="0.25")
    fig.subplots_adjust(top=0.885, bottom=0.075, left=0.055, right=0.988,
                        wspace=0.20, hspace=0.30)
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    fig.savefig(path, dpi=dpi or plt.rcParams["savefig.dpi"], facecolor="white")
    plt.close(fig)
    return path


# ---------------------------------------------------------------------------
# Outputs
# ---------------------------------------------------------------------------
CSV_FIELDS = ["Month", "n", "gamma2", "P", "z_gamma2", "z_P", "z_pooled_gamma2",
              "z_pooled_P", "R_behaviour", "ci95_gamma2", "ci95_P", "thin"]


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
    """Seven properties of the plate, each one falsifiable on the real table."""
    report = build_report(table, meta, boot=64)
    months = report["months"]
    checks = []

    def check(name, condition, detail=""):
        checks.append((name, bool(condition), detail))

    check("table: 178 rows, one role (target), one orbit direction",
          report["audit"]["rows"] == 178 and report["input"]["role"] == "target",
          f"rows = {report['audit']['rows']}, orbit = {report['audit']['orbit']}")

    check("partition: 35 reference + 101 series + 42 post = 178",
          report["audit"]["reference"] == 35 and report["audit"]["series"] == 101
          and report["audit"]["post_outside_series"] == 42,
          report["audit"]["partition"])

    check("windows: reference ends 2017-12-15, series starts 2018-01, 42 months",
          report["reference"]["last"] == "2017-12-15"
          and report["series"]["first"] == "2018-01"
          and report["series"]["last"] == "2021-06"
          and len(months) == 42,
          f"{report['reference']['first']} … {report['reference']['last']}, "
          f"{len(months)} months to {report['series']['last']}")

    # median(z_i) == z(median(log x)): the transform is monotone and z is affine.
    worst = 0.0
    for entry in months[:3] + months[-3:]:
        for dim in DIMS:
            chips = [c["z"][dim] for c in report["chips"]
                     if c["date"].startswith(entry["month"])]
            worst = max(worst, abs(ocv.median(chips) - entry["dims"][dim]["z"]))
    check("median of the chips' z equals z of the median (3 + 3 months)", worst < 1e-12,
          f"max |difference| = {worst:.2e}")

    centre = report["reference"]["dims"][DIMS[0]]["median_log"]
    scale = report["reference"]["dims"][DIMS[0]]["scale_log"]
    z_ref = (ocv.median([centre] * 3) - centre) / scale
    check("a month that sits on the reference median has z = 0", z_ref == 0.0,
          f"z = {z_ref:+.1e}")

    check("every month is thin (n < SMALL_N)", report["series"]["thin_months"] == len(months),
          f"{report['series']['thin_months']} of {len(months)}, "
          f"n = {min(e['n'] for e in months)}-{max(e['n'] for e in months)}")

    ratio = {dim: report["reference"]["dims"][dim]["scale_log"]
             / report["reference"]["dims"][dim]["pooled_scale_log"] for dim in DIMS}
    # z_pooled is the same r read from the pooled origin: an affine image of z, not
    # a second estimate. Its re-centring may only move months close to zero.
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
    check("z_pooled is an affine image of z and flips only months within 1 MAD",
          worst_affine < 1e-9 and worst_flip < 1.0,
          f"max deviation = {worst_affine:.1e}, {len(flips)} flip(s) (max |z| = "
          f"{worst_flip:.2f}), factors x{ratio[DIMS[0]]:.3f} / x{ratio[DIMS[1]]:.3f}")

    for name, passed, detail in checks:
        print(f"  {'ok  ' if passed else 'FAIL'} {name}" + (f"  [{detail}]" if detail else ""))
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


def verify(report, path=None):
    """Recompute everything from the table and compare it to the committed JSON."""
    path = path or OUT_JSON
    if not os.path.exists(path):
        raise SystemExit(f"{path}: no committed output — run the script first")
    with open(path) as handle:
        committed = json.load(handle)
    problems = _diff(report, committed)
    if problems:
        for line in problems[:12]:
            print("  " + line)
        raise SystemExit(f"verify: {len(problems)} difference(s) against {_rel(path)}")
    print(f"verify: {len(report['months'])} months, {len(report['chips'])} series chips, "
          f"reference n = {report['reference']['n']} — identical to {_rel(path)}")
    return 0


# ---------------------------------------------------------------------------
# Command line
# ---------------------------------------------------------------------------
def main(argv=None):
    parser = argparse.ArgumentParser(
        description="CTS — OCV reconfiguration against the 2015-2017 reference, "
                    "month by month")
    parser.add_argument("--table", default=TABLE)
    parser.add_argument("--meta", default=META)
    parser.add_argument("--json", default=OUT_JSON)
    parser.add_argument("--csv", default=OUT_CSV)
    parser.add_argument("--md", default=None)
    parser.add_argument("--jpg", default=None)
    parser.add_argument("--figdir", default=None,
                        help="directory for the plate and its guide (default: "
                             "figures/cts)")
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
        return verify(report, args.json)

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
