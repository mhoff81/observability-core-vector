#!/usr/bin/env python3
"""Pre-collapse diagnostics for the Carolabrücke deck — trend, and late level.

The orbit figures draw a rolling median, and on this deck that line is **not**
flat before the collapse: some dimensions move. A line is not a test, so this
script asks the two questions the line raises, and only those:

  1. **trend** — inside the pre-collapse observations, does ``D``, ``F`` or
     ``A`` change systematically with time?
  2. **late level** — is the last stretch before the collapse (three months)
     different from the older pre-collapse baseline?

Both are computed **per observation channel and per dimension**. ASC and DESC
are never pooled: they see different geometries, and this project treats the
channel as the unit of the claim (see ``code/sds/README.md``).

Stage 1 — trend, pre-collapse only
----------------------------------
  slope  Theil-Sen (median of pairwise slopes — an outlier cannot move it),
         reported per year and as percent of the baseline median per year
  CI     moving-block bootstrap, block 5 overpasses (about one month ASC),
         ``espoo_ocv_stats.block_bootstrap_slope``; the blocks keep the serial
         dependence that an independent bootstrap would destroy
  test   Mann-Kendall, tie-corrected (tau_b), two-sided

**The block-bootstrap interval is the primary evidence; the p-value is not.**
Mann-Kendall assumes independent observations, and consecutive overpasses of the
same deck are not independent, so its p is anti-conservative here. It is
reported as a direction indicator beside the interval, never instead of it.
``lag1_autocorr`` and ``runs_p`` say how far the independence assumption is from
holding at all.

Stage 2 — late level, last three months versus the older baseline
-----------------------------------------------------------------
  baseline  pre-collapse overpasses with ``t < event - 3 months``
  late-pre  pre-collapse overpasses with ``t >= event - 3 months``
  test      Mann-Whitney U, two-sided (ties force scipy's asymptotic branch;
            the number of ties in the pooled sample is reported)
  effect    Cliff's delta, **positive = late-pre larger than baseline**, so a
            sinking dimension gives a negative delta
  CIs       moving-block bootstrap (block 3, 5000 resamples) for delta and for
            the median difference — the house ``bootstrap_delta_ci`` resamples
            independently and would understate the dependence
  correction Holm over the three dimensions, per channel and per stage

Both boundaries are conventions, so the report carries a sensitivity ladder
(2 / 3 / 4 / 6 months) and the event-date shift (-12 / 0 / +12 days). A finding
that exists at one window only is a finding about the window.

What the numbers do not say
---------------------------
* A negative delta is a statement about three months before a date, with
  ``n_late`` of 6 to 16 — not a prediction, not a causal link, and **not a new
  state**. This script writes no state label at all; the states stay the two the
  committed table already carries (``selftest`` asserts that no state word
  reaches the report).
* DESC has a gap directly before the collapse (its last pre-collapse overpass is
  2024-08-20), so its "last three months" ends three weeks before the event.

Usage
-----
  python3 code/carola/carola_ocv_precollapse.py
  python3 code/carola/carola_ocv_precollapse.py --boot 500 --json /tmp/x.json
  python3 code/carola/carola_ocv_precollapse.py --selftest
  python3 code/carola/carola_ocv_precollapse.py --verify

Outputs
-------
  data/carola/support/carola_ocv_precollapse.json
  data/carola/support/carola_ocv_precollapse.md
"""
from __future__ import annotations

import argparse
import calendar
import csv
import datetime as dt
import hashlib
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
for path in (os.path.join(REPO, "code", "espoo"),
             os.path.join(REPO, "code", "lumo")):
    if path not in sys.path:
        sys.path.insert(0, path)

import espoo_ocv_stats as espoo  # noqa: E402
from lumo_ocv_stats import cliffs_delta, welch_mw_test  # noqa: E402

SCHEMA = "ocv-precollapse-v1"
SITE = "carola"
LABEL = "Carolabrücke deck"
TABLE = os.path.join(REPO, "data", SITE, "support", f"{SITE}_ocv_channels.csv")
OUT_JSON = os.path.join(REPO, "data", SITE, "support",
                        f"{SITE}_ocv_precollapse.json")
OUT_MD = os.path.join(REPO, "data", SITE, "support",
                      f"{SITE}_ocv_precollapse.md")

#: Observation channels, in the order the report prints them.
CHANNELS = ("ASCENDING", "DESCENDING")
#: The three dimensions this test was specified for.
DIMS = ("D", "F", "A")
UNITS = {"D": "m", "F": "Hz", "A": "count"}
#: The documented collapse of the Carolabrücke.
EVENT = dt.datetime(2024, 9, 11)
WINDOW_MONTHS = 3
LADDER_MONTHS = (2, 3, 4, 6)
SHIFT_DAYS = (-12, 0, 12)
#: Moving-block lengths, in overpasses.
SLOPE_BLOCK = 5
DELTA_BLOCK = 3
BOOT_SLOPE = 2000
BOOT_DELTA = 5000
SEED_SLOPE = 13
SEED_DELTA = 7
SMALL_N = 8


# ---------------------------------------------------------------------------
# Input
# ---------------------------------------------------------------------------
def _naive(ts):
    """The table stores UTC; compare naive datetimes so the cut is unambiguous."""
    when = dt.datetime.fromisoformat(ts)
    if when.tzinfo is not None:
        when = when.astimezone(dt.timezone.utc).replace(tzinfo=None)
    return when


def load_pre(table=TABLE, channel="ASCENDING", state="pre"):
    """Pre-collapse rows of one channel, in acquisition order."""
    rows = []
    with open(table, newline="") as handle:
        for row in csv.DictReader(handle):
            if row.get("orbit_direction") != channel:
                continue
            if not (row.get("state") or "").startswith(state):
                continue
            rows.append(row)
    rows.sort(key=lambda row: _naive(row["acquisition_ts"]))
    return rows


def series(rows, dim):
    """``[(t, v), ...]`` — a missing value is dropped, never read as zero."""
    out = []
    for row in rows:
        try:
            value = float(row.get(dim, ""))
        except (TypeError, ValueError):
            continue
        out.append((_naive(row["acquisition_ts"]), value))
    return out


def months_before(when, months):
    """Calendar subtraction: shift the month, clamp the day to its length."""
    month = when.month - months
    year = when.year + (month - 1) // 12
    month = (month - 1) % 12 + 1
    day = min(when.day, calendar.monthrange(year, month)[1])
    return when.replace(year=year, month=month, day=day)


def sha256_of(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _iqr(values):
    if not values:
        return None
    low, high = np.percentile(np.asarray(values, float), [25, 75])
    return [float(low), float(high)]


# ---------------------------------------------------------------------------
# Correction and resampling — the two pieces the house did not have
# ---------------------------------------------------------------------------
def holm(pvals):
    """Holm's step-down adjustment.

    ``None`` is skipped and does **not** count towards the family size: an
    unavailable test must not make the available ones look better. Returns
    ``(adjusted, family_size)``; the family size is reported with the numbers.
    """
    known = [(index, p) for index, p in enumerate(pvals) if p is not None]
    out = [None] * len(pvals)
    family = len(known)
    if not family:
        return out, 0
    running = 0.0
    for rank, (index, p) in enumerate(sorted(known, key=lambda kv: kv[1]),
                                      start=1):
        running = max(running, min(1.0, (family - rank + 1) * float(p)))
        out[index] = running
    return out, family


def _block_resample(values, block, rng):
    """One moving-block resample of a **time-ordered** group."""
    array = np.asarray(values, float)
    m = array.size
    if m <= 1:
        return array.copy()
    block = max(1, min(block, m))
    starts = rng.integers(0, m - block + 1, int(np.ceil(m / block)))
    index = np.concatenate([np.arange(s, s + block) for s in starts])[:m]
    return array[index]


def bootstrap_delta_ci_block(ref, grp, block=DELTA_BLOCK, n=BOOT_DELTA,
                             seed=SEED_DELTA):
    """Percentile 95% CI of Cliff's delta with moving-block resampling.

    ``lumo_ocv_stats.bootstrap_delta_ci`` resamples every observation on its
    own; consecutive overpasses are not independent, so here the blocks of each
    group move together and the interval stops pretending otherwise.
    """
    a = np.asarray(ref, float)
    b = np.asarray(grp, float)
    if a.size < 2 or b.size < 2:
        return None
    rng = np.random.default_rng(seed)
    out = np.empty(n)
    for index in range(n):
        out[index] = cliffs_delta(_block_resample(a, block, rng),
                                  _block_resample(b, block, rng))
    return [float(np.percentile(out, 2.5)), float(np.percentile(out, 97.5))]


def bootstrap_median_diff_block(ref, grp, block=DELTA_BLOCK, n=BOOT_DELTA,
                                seed=SEED_DELTA):
    """Percentile 95% CI of ``median(late) - median(baseline)``, blocks kept."""
    a = np.asarray(ref, float)
    b = np.asarray(grp, float)
    if a.size < 2 or b.size < 2:
        return None
    rng = np.random.default_rng(seed)
    out = np.empty(n)
    for index in range(n):
        out[index] = (np.median(_block_resample(b, block, rng))
                      - np.median(_block_resample(a, block, rng)))
    return [float(np.percentile(out, 2.5)), float(np.percentile(out, 97.5))]


# ---------------------------------------------------------------------------
# Stage 1 — trend inside the pre-collapse observations
# ---------------------------------------------------------------------------
def trend(points, dim):
    """Theil-Sen slope + moving-block CI + Mann-Kendall, for one dimension."""
    out = {"dimension": dim, "unit": UNITS.get(dim), "n": len(points)}
    if len(points) < 2 * SLOPE_BLOCK or len(points) < SMALL_N:
        out.update({"available": False,
                    "reason": f"needs {max(2 * SLOPE_BLOCK, SMALL_N)} "
                              f"pre-collapse values"})
        return out
    t0 = points[0][0]
    x = np.array([(when - t0).total_seconds() / 86400.0 for when, _v in points])
    y = np.array([value for _when, value in points])
    slope_day = espoo.theil_sen(x, y)
    interval = espoo.block_bootstrap_slope(x, y, block=SLOPE_BLOCK,
                                           n=BOOT_SLOPE, seed=SEED_SLOPE,
                                           years=True)
    mk = espoo.mann_kendall(y)
    base = float(np.median(y))
    slope_year = None if slope_day is None else slope_day * 365.25
    out.update({
        "available": True,
        "t0": t0.isoformat(),
        "span_days": float(x[-1] - x[0]),
        "baseline_median": base,
        "slope_per_year": slope_year,
        "slope_ci95_per_year": interval,
        "slope_rel_pct_per_year": (None if slope_year is None or base == 0
                                   else slope_year / base * 100.0),
        "kendall_tau": mk.get("tau"),
        "kendall_p": mk.get("p"),
        "kendall_s": mk.get("S"),
        "lag1_autocorr": espoo.lag1_autocorr(y),
        "runs_p": (espoo.runs_test(y) or {}).get("p"),
        "block": SLOPE_BLOCK,
        "boot": BOOT_SLOPE,
        "seed": SEED_SLOPE,
        "small_n": len(points) < 2 * SLOPE_BLOCK,
    })
    return out


# ---------------------------------------------------------------------------
# Stage 2 — late level against the older baseline
# ---------------------------------------------------------------------------
def level(points, cut, dim, boot=BOOT_DELTA, with_ci=True):
    """Mann-Whitney + Cliff's delta + moving-block CIs, for one dimension."""
    base = [value for when, value in points if when < cut]
    late = [value for when, value in points if when >= cut]
    out = {"dimension": dim, "unit": UNITS.get(dim), "cut": cut.isoformat(),
           "n": {"baseline": len(base), "late_pre": len(late)}}
    if len(base) < 3 or len(late) < 3:
        out.update({"available": False,
                    "reason": "fewer than 3 values in a group"})
        return out
    test = welch_mw_test(base, late)
    delta = cliffs_delta(base, late)
    pooled = np.concatenate([np.asarray(base, float), np.asarray(late, float)])
    ties = int(pooled.size - np.unique(pooled).size)
    out.update({
        "available": True,
        "baseline_median": float(np.median(base)),
        "late_median": float(np.median(late)),
        "baseline_iqr": _iqr(base),
        "late_iqr": _iqr(late),
        "delta": delta,
        "median_diff": float(np.median(late) - np.median(base)),
        "mannwhitney_p": test.get("mannwhitney_p"),
        "welch_p": test.get("welch_p"),
        "ties_in_pooled": ties,
        "test": ("mann-whitney U, two-sided; scipy method=auto "
                 + ("(asymptotic forced by ties)" if ties else "(exact possible)")),
        "small_n": len(late) < SMALL_N,
    })
    if with_ci:
        out["delta_ci95"] = bootstrap_delta_ci_block(base, late, n=boot)
        out["median_diff_ci95"] = bootstrap_median_diff_block(base, late, n=boot)
    return out


# ---------------------------------------------------------------------------
# Assembling one channel: Holm within (channel, stage)
# ---------------------------------------------------------------------------
def analyse_channel(channel, table=TABLE, cut=None, boot=BOOT_DELTA):
    """Trend + late level for one channel, with the Holm family applied."""
    rows = load_pre(table, channel)
    if not rows:
        return {"available": False, "n_pre": 0}
    cut = cut or months_before(EVENT, WINDOW_MONTHS)
    times = [_naive(row["acquisition_ts"]) for row in rows]
    gaps = sorted((times[i + 1] - times[i]).days for i in range(len(times) - 1))
    block = {
        "available": True,
        "n_pre": len(rows),
        "first": times[0].isoformat(),
        "last": times[-1].isoformat(),
        "median_gap_days": gaps[len(gaps) // 2] if gaps else None,
        "gap_before_event_days": (EVENT - times[-1]).days,
        "cut": cut.isoformat(),
        "trend": {}, "level": {},
    }
    for dim in DIMS:
        points = series(rows, dim)
        block["trend"][dim] = trend(points, dim)
        block["level"][dim] = level(points, cut, dim, boot=boot)
    adjusted_level, size_level = holm([block["level"][dim].get("mannwhitney_p")
                                       for dim in DIMS])
    adjusted_trend, size_trend = holm([block["trend"][dim].get("kendall_p")
                                       for dim in DIMS])
    for dim, p in zip(DIMS, adjusted_level):
        block["level"][dim]["holm_p"] = p
    for dim, p in zip(DIMS, adjusted_trend):
        block["trend"][dim]["holm_p"] = p
    block["holm"] = {
        "rule": "Holm over the three dimensions, within this channel and stage",
        "sizes": {"trend": size_trend, "level": size_level},
        "alpha": 0.05,
    }
    return block


# ---------------------------------------------------------------------------
# The window is a convention — show what moving it does
# ---------------------------------------------------------------------------
def sensitivity(table=TABLE):
    """The same test at 2/3/4/6 months and at event -12/+0/+12 days.

    No bootstraps here on purpose: the ladder answers *does the verdict depend
    on the boundary*, so it reports delta, p and the Holm decision only. The
    intervals belong to the headline window and are reported once, with it.
    """
    out = {}
    for channel in CHANNELS:
        rows = load_pre(table, channel)
        if not rows:
            continue
        grid = {}
        for months in LADDER_MONTHS:
            for shift in SHIFT_DAYS:
                cut = months_before(EVENT + dt.timedelta(days=shift), months)
                deltas, pvals, counts = {}, {}, {}
                for dim in DIMS:
                    cell = level(series(rows, dim), cut, dim, with_ci=False)
                    deltas[dim] = cell.get("delta")
                    pvals[dim] = cell.get("mannwhitney_p")
                    counts[dim] = cell.get("n", {}).get("late_pre")
                adjusted, family = holm([pvals[dim] for dim in DIMS])
                grid[f"{months}m{shift:+d}d"] = {
                    "months": months, "shift_days": shift,
                    "cut": cut.isoformat(), "n_late": counts,
                    "delta": deltas, "p": pvals,
                    "holm_p": dict(zip(DIMS, adjusted)), "family": family,
                    "significant": [dim for dim, p in zip(DIMS, adjusted)
                                    if p is not None and p < 0.05],
                }
        out[channel] = grid
    return out


def analysis(table=TABLE, boot=BOOT_DELTA, ladder=True):
    """The whole report, as it is written to JSON."""
    cut = months_before(EVENT, WINDOW_MONTHS)
    report = {
        "schema": SCHEMA,
        "site": SITE,
        "label": LABEL,
        "table": {"path": os.path.relpath(table, REPO),
                  "sha256": sha256_of(table)},
        "event": EVENT.isoformat(),
        "window_months": WINDOW_MONTHS,
        "cut": cut.isoformat(),
        "dims": list(DIMS),
        "units": UNITS,
        "rules": {
            "primary": "the block-bootstrap interval; the p-values are "
                       "indicators, not the claim",
            "channels": "ASCENDING and DESCENDING are analysed separately and "
                        "never pooled",
            "trend_block": SLOPE_BLOCK,
            "delta_block": DELTA_BLOCK,
            "boot_slope": BOOT_SLOPE,
            "boot_delta": BOOT_DELTA,
            "seed_slope": SEED_SLOPE,
            "seed_delta": SEED_DELTA,
            "delta_sign": "positive = late-pre larger than baseline",
            "holm": "per channel and per stage, over the three dimensions",
            "small_n": SMALL_N,
            "state_labels_written": False,
        },
        "channels": {},
    }
    for channel in CHANNELS:
        report["channels"][channel] = analyse_channel(channel, table, cut=cut,
                                                      boot=boot)
    if ladder:
        report["sensitivity"] = sensitivity(table)
    return report


# ---------------------------------------------------------------------------
# Markdown
# ---------------------------------------------------------------------------
def _num(value, digits=3, signed=False):
    if value is None:
        return "—"
    if isinstance(value, str):
        return value
    return f"{value:+.{digits}f}" if signed else f"{value:.{digits}f}"


def _interval(pair, digits=3, signed=True):
    if not pair:
        return "—"
    return f"[{_num(pair[0], digits, signed)}, {_num(pair[1], digits, signed)}]"


def _p(value):
    if value is None:
        return "—"
    return "<0.001" if value < 0.001 else f"{value:.3f}"


def to_markdown(report):
    """The reviewer table, primary evidence first."""
    lines = []
    add = lines.append
    add(f"# {report['label']} — pre-collapse trend and late level")
    add("")
    add(f"* table: `{report['table']['path']}` "
        f"(sha256 `{report['table']['sha256'][:16]}…`)")
    add(f"* event (documented collapse): **{report['event']}**")
    add(f"* window: last **{report['window_months']} months** before the event — "
        f"baseline `t < {report['cut']}`, late-pre "
        f"`{report['cut']} <= t < {report['event']}`")
    add("* dimensions: "
        + ", ".join(f"`{dim}`" for dim in report["dims"])
        + " — the unit is the one the table carries")
    add("* the two channels are analysed **separately and never pooled**")
    add("")
    add("**How to read this.** The numbers to trust are the intervals: the "
        "moving-block bootstrap CI of the slope, and the CI of Cliff's delta. "
        "The p-values are printed beside them because they were asked for, but "
        "Mann-Kendall assumes independent observations and consecutive "
        "overpasses of one deck are not independent, so its p is "
        "anti-conservative. `lag-1` and `runs p` say how far that assumption is "
        "from holding.")
    add("")
    add("**Sign of the effect.** Cliff's delta is positive when the late-pre "
        "values are *larger* than the baseline, so a sinking dimension gives a "
        "**negative** delta. Nothing here is a new state: this script writes no "
        "state label, and the states stay the two the committed table carries.")
    add("")

    for channel, block in report["channels"].items():
        if not block.get("available"):
            add(f"## {channel} — no pre-collapse overpasses")
            add("")
            continue
        add(f"## {channel}")
        add("")
        add(f"*n* = {block['n_pre']} pre-collapse overpasses, "
            f"{block['first'][:10]} … {block['last'][:10]}, median revisit "
            f"{block['median_gap_days']} d. Last pre-collapse overpass "
            f"**{block['gap_before_event_days']} d before the event**.")
        add("")
        add("### Trend inside the pre-collapse series")
        add("")
        add("| dim | u | n | Theil-Sen slope | CI95 (block 5) | rel. %/a | "
            "Kendall τ | MK p | Holm p | lag-1 | runs p |")
        add("|---|---|---|---|---|---|---|---|---|---|---|")
        for dim in report["dims"]:
            cell = block["trend"][dim]
            if not cell.get("available"):
                add(f"| `{dim}` | {cell.get('unit', '—')} | {cell['n']} | "
                    f"*{cell.get('reason', 'n/a')}* | — | — | — | — | — | — | — |")
                continue
            mark = " †" if cell.get("small_n") else ""
            add(f"| `{dim}`{mark} | {cell['unit']} | {cell['n']} | "
                f"{_num(cell['slope_per_year'], 4, signed=True)} | "
                f"{_interval(cell['slope_ci95_per_year'], 4)} | "
                f"{_num(cell['slope_rel_pct_per_year'], 1, signed=True)} % | "
                f"{_num(cell['kendall_tau'], 3, signed=True)} | "
                f"{_p(cell['kendall_p'])} | {_p(cell['holm_p'])} | "
                f"{_num(cell['lag1_autocorr'], 2, signed=True)} | "
                f"{_p(cell['runs_p'])} |")
        add("")
        add("### Late level — last "
            f"{report['window_months']} months versus the older baseline")
        add("")
        add("| dim | n base / late | baseline median | late median | "
            "Δ median [CI95] | Cliff's δ [CI95] | MW p | Holm p | ties |")
        add("|---|---|---|---|---|---|---|---|---|")
        for dim in report["dims"]:
            cell = block["level"][dim]
            if not cell.get("available"):
                add(f"| `{dim}` | — | — | — | *{cell.get('reason', 'n/a')}* | "
                    "— | — | — | — |")
                continue
            mark = " †" if cell.get("small_n") else ""
            add(f"| `{dim}`{mark} | {cell['n']['baseline']} / "
                f"{cell['n']['late_pre']} | "
                f"{_num(cell['baseline_median'], 4)} | "
                f"{_num(cell['late_median'], 4)} | "
                f"{_num(cell['median_diff'], 4, signed=True)} "
                f"{_interval(cell['median_diff_ci95'], 4)} | "
                f"{_num(cell['delta'], 3, signed=True)} "
                f"{_interval(cell['delta_ci95'], 3)} | "
                f"{_p(cell['mannwhitney_p'])} | {_p(cell['holm_p'])} | "
                f"{cell['ties_in_pooled']} |")
        add("")
        add(f"Holm family: {block['holm']['sizes']['level']} tests (level), "
            f"{block['holm']['sizes']['trend']} tests (trend), alpha = "
            f"{block['holm']['alpha']}. † = fewer than "
            f"{report['rules']['small_n']} values in the late window.")
        add("")
    add("## Sensitivity — does the verdict depend on the boundary?")
    add("")
    add("Delta per dimension and the Holm decision per window; the event date is "
        "also shifted by -12 / +0 / +12 days. No intervals here on purpose: the "
        "ladder answers *does the window drive the answer*, not how large the "
        "effect is.")
    add("")
    add("| channel | window | cut | n late (D/F/A) | δ D | δ F | δ A | "
        "Holm-significant |")
    add("|---|---|---|---|---|---|---|---|")
    for channel, grid in (report.get("sensitivity") or {}).items():
        for key, cell in grid.items():
            counts = cell["n_late"]
            add(f"| {channel[:4]} | {key} | {cell['cut'][:10]} | "
                + "/".join(str(counts[dim]) for dim in report["dims"])
                + " | "
                + " | ".join(_num(cell["delta"][dim], 2, signed=True)
                             for dim in report["dims"])
                + f" | {', '.join(cell['significant']) or '—'} |")
    add("")
    add("## What these numbers do not say")
    add("")
    add("* a negative delta is a statement about the last three months before a "
        "date, with the `n` printed beside it — not a prediction and not a "
        "cause;")
    add("* **no state is created here.** A rolling median and these tests can "
        "both show a move without the deck entering a third state; this script "
        "writes no state label at all, and the guard in `--selftest` checks "
        "that none reaches the report;")
    add("* a result that exists at one window width only is a finding about the "
        "window — the ladder is the check, not a search for the significant "
        "cell;")
    add("* `DESCENDING` has a gap directly before the collapse (see "
        "`gap_before_event_days` above), so its late window ends before the "
        "event does.")
    add("")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Output, verification, self-test
# ---------------------------------------------------------------------------
def write_outputs(report, json_path=OUT_JSON, md_path=OUT_MD, quiet=False):
    for path in (json_path, md_path):
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(json_path, "w") as handle:
        json.dump(report, handle, indent=1)
        handle.write("\n")
    with open(md_path, "w") as handle:
        handle.write(to_markdown(report))
        handle.write("\n")
    if not quiet:
        print(f"  wrote {os.path.relpath(json_path, REPO)}")
        print(f"  wrote {os.path.relpath(md_path, REPO)}")
    return json_path, md_path


def _strings(node):
    """Every string value in a nested report — keys are not values."""
    if isinstance(node, str):
        yield node
    elif isinstance(node, dict):
        for value in node.values():
            yield from _strings(value)
    elif isinstance(node, (list, tuple)):
        for item in node:
            yield from _strings(item)


def _same(left, right, tol=1e-9):
    if left is None or right is None:
        return left is None and right is None
    if isinstance(left, (list, tuple)) or isinstance(right, (list, tuple)):
        if not (isinstance(left, (list, tuple)) and isinstance(right, (list, tuple))):
            return False
        return len(left) == len(right) and all(
            _same(a, b, tol) for a, b in zip(left, right))
    try:
        return abs(float(left) - float(right)) <= tol
    except (TypeError, ValueError):
        return left == right


def verify(report=None, table=TABLE):
    """Recompute every headline number from the table and compare with the JSON."""
    if report is None:
        with open(OUT_JSON) as handle:
            report = json.load(handle)
    failures = []
    if report.get("schema") != SCHEMA:
        failures.append(f"schema {report.get('schema')!r} != {SCHEMA!r}")
    if report["table"]["sha256"] != sha256_of(table):
        failures.append("the table on disk is not the table this report was "
                        "computed from (sha256 differs)")
        print(f"FAILED — {failures}")
        return 1
    fresh = analysis(table, boot=report["rules"]["boot_delta"], ladder=False)
    checked = 0
    for channel in CHANNELS:
        new = fresh["channels"][channel]
        old = report["channels"][channel]
        for key in ("n_pre", "first", "last", "median_gap_days",
                    "gap_before_event_days", "cut"):
            checked += 1
            if not _same(new.get(key), old.get(key)):
                failures.append(f"{channel}: {key} {new.get(key)} != {old.get(key)}")
        for dim in DIMS:
            for stage in ("trend", "level"):
                left, right = new[stage][dim], old[stage][dim]
                for key in ("available", "n", "slope_per_year", "kendall_tau",
                            "kendall_p", "slope_ci95_per_year",
                            "baseline_median", "late_median", "delta",
                            "delta_ci95", "median_diff", "median_diff_ci95",
                            "mannwhitney_p", "holm_p", "ties_in_pooled"):
                    if key not in right:
                        continue
                    checked += 1
                    if not _same(left.get(key), right.get(key)):
                        failures.append(f"{channel}/{dim}/{stage}: {key} "
                                        f"{left.get(key)} != {right.get(key)}")
    print(f"  {checked} field(s) recomputed from "
          f"{os.path.relpath(table, REPO)}")
    if failures:
        print(f"FAILED — {len(failures)} mismatch(es): {failures[:5]}")
        return 1
    print("  all recomputed values match the JSON")
    return 0


def selftest():
    """Offline checks: Holm, the block bootstrap, both stages, the guard."""
    failures = []

    def check(label, condition, detail=""):
        print(f"  [{'ok  ' if condition else 'FAIL'}] {label}"
              + (f" — {detail}" if detail else ""))
        if not condition:
            failures.append(label)

    print("\n1. Holm, by hand")
    adjusted, family = holm([0.01, 0.02, 0.03])
    check("0.01 / 0.02 / 0.03 adjust to 0.03 / 0.04 / 0.04",
          [round(v, 10) for v in adjusted] == [0.03, 0.04, 0.04], str(adjusted))
    check("the family is the three available tests", family == 3, str(family))
    adjusted, family = holm([0.01, None, 0.03])
    check("a missing test does not enlarge the family", family == 2, str(family))
    check("and it stays missing", adjusted[1] is None, str(adjusted))
    check("the adjustment never lowers a p-value",
          adjusted[0] >= 0.01 and adjusted[2] >= 0.03, str(adjusted))

    print("\n2. The moving-block bootstrap")
    rng = np.random.default_rng(11)
    a = rng.normal(10.0, 0.5, 60)
    b = rng.normal(9.0, 0.5, 12)
    first = bootstrap_delta_ci_block(a, b, n=400, seed=7)
    again = bootstrap_delta_ci_block(a, b, n=400, seed=7)
    check("same seed, same interval (deterministic)", first == again, str(first))
    other = bootstrap_delta_ci_block(a, b, n=400, seed=8)
    check("a different seed moves it", first != other, f"{first} vs {other}")
    drawn = _block_resample(a, 3, np.random.default_rng(5))
    check("a resample only ever repeats input values",
          set(np.round(drawn, 12)) <= set(np.round(a, 12)))
    check("and keeps the group size", drawn.size == a.size,
          f"{drawn.size} vs {a.size}")

    print("\n3. Stage 1 finds a planted trend, and not one in noise")
    days = np.arange(0, 40) * 6.0
    start = dt.datetime(2022, 9, 3)
    planted = (2.0 + 0.001 * days
               + np.random.default_rng(2).normal(0, 0.005, days.size))
    cells = [(start + dt.timedelta(days=float(d)), v)
             for d, v in zip(days, planted)]
    cell = trend(cells, "D")
    interval = cell["slope_ci95_per_year"]
    check("the CI contains the planted 0.365/a",
          interval is not None and interval[0] <= 0.36525 <= interval[1],
          f"CI {interval}")
    check("and Mann-Kendall sees the direction",
          cell["kendall_p"] < 0.05 and cell["kendall_tau"] > 0,
          f"tau {cell['kendall_tau']}, p {cell['kendall_p']}")
    noise = np.random.default_rng(3).normal(5.0, 1.0, days.size)
    quiet = [(start + dt.timedelta(days=float(d)), v) for d, v in zip(days, noise)]
    cell = trend(quiet, "D")
    interval = cell["slope_ci95_per_year"]
    check("noise gives a CI that contains zero",
          interval is not None and interval[0] <= 0 <= interval[1],
          f"CI {interval}")

    print("\n4. Stage 2 finds a planted drop, and not one in noise")
    cut = dt.datetime(2024, 6, 11)
    base = np.random.default_rng(5).normal(10.0, 0.4, 60)
    late = np.random.default_rng(6).normal(9.0, 0.4, 12)
    shifted = ([(cut - dt.timedelta(days=6 * (60 - i)), v)
                for i, v in enumerate(base)]
               + [(cut + dt.timedelta(days=6 * i), v)
                  for i, v in enumerate(late)])
    cell = level(shifted, cut, "D", boot=400)
    check("the drop gives a negative delta", cell["delta"] < 0,
          f"delta {cell['delta']:.3f}")
    check("its CI excludes zero", cell["delta_ci95"][1] < 0,
          str(cell["delta_ci95"]))
    check("the late window holds exactly the 12 planted points",
          cell["n"] == {"baseline": 60, "late_pre": 12}, str(cell["n"]))
    check("and the median difference points the same way",
          cell["median_diff"] < 0, f"{cell['median_diff']:.3f}")
    flat_late = np.random.default_rng(7).normal(10.0, 0.4, 12)
    flat = ([(cut - dt.timedelta(days=6 * (60 - i)), v)
             for i, v in enumerate(base)]
            + [(cut + dt.timedelta(days=6 * i), v)
               for i, v in enumerate(flat_late)])
    cell = level(flat, cut, "D", boot=400)
    check("no shift: the delta CI contains zero",
          cell["delta_ci95"][0] <= 0 <= cell["delta_ci95"][1],
          str(cell["delta_ci95"]))

    print("\n5. Window arithmetic")
    check("three months before 2024-09-11 is 2024-06-11",
          months_before(EVENT, 3).date() == dt.date(2024, 6, 11),
          months_before(EVENT, 3).date().isoformat())
    check("month-end clamps to a real day",
          months_before(dt.datetime(2024, 5, 31), 3).date() == dt.date(2024, 2, 29),
          months_before(dt.datetime(2024, 5, 31), 3).date().isoformat())
    check("and the clamp does not wrap into the next month",
          months_before(dt.datetime(2023, 1, 31), 1).date() == dt.date(2022, 12, 31),
          months_before(dt.datetime(2023, 1, 31), 1).date().isoformat())

    print("\n6. No state label is invented anywhere")
    report = analysis(boot=200, ladder=False)
    forbidden = ("state", "pre-collapse", "post-collapse", "damage", "healthy")
    found = sorted({text for text in _strings(report)
                    if any(word in text.lower() for word in forbidden)})
    check("no string value of the report names a state", not found,
          f"found {found[:3]}")
    check("the report declares that it writes none",
          report["rules"]["state_labels_written"] is False)
    check("the two channels are separate blocks",
          sorted(report["channels"]) == sorted(CHANNELS),
          str(sorted(report["channels"])))
    check("DESCENDING keeps its pre-collapse gap on the record",
          report["channels"]["DESCENDING"]["gap_before_event_days"] > 14,
          f"{report['channels']['DESCENDING']['gap_before_event_days']} d")

    print()
    if failures:
        print(f"FAILED — {len(failures)} check(s): {failures}")
        return 1
    print("all checks passed")
    return 0


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def build_parser():
    parser = argparse.ArgumentParser(
        description=("Pre-collapse trend and late level of the Carolabrücke OCV "
                     "dimensions, per observation channel"))
    parser.add_argument("--table", default=TABLE, help="channel table to read")
    parser.add_argument("--json", default=OUT_JSON)
    parser.add_argument("--md", default=OUT_MD)
    parser.add_argument("--boot", type=int, default=BOOT_DELTA,
                        help=f"moving-block resamples for the delta and median "
                             f"CIs (default {BOOT_DELTA})")
    parser.add_argument("--no-ladder", action="store_true",
                        help="skip the 2/3/4/6-month sensitivity grid")
    parser.add_argument("--selftest", action="store_true",
                        help="run the offline checks and exit")
    parser.add_argument("--verify", action="store_true",
                        help="recompute the JSON from the table and compare")
    parser.add_argument("--quiet", action="store_true")
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    if args.selftest:
        return selftest()
    if args.verify:
        return verify(table=args.table)

    print(f"{LABEL} — pre-collapse trend and late level")
    report = analysis(args.table, boot=args.boot, ladder=not args.no_ladder)
    write_outputs(report, args.json, args.md, quiet=args.quiet)
    if not args.quiet:
        for channel, block in report["channels"].items():
            if not block.get("available"):
                print(f"  {channel[:4]} no pre-collapse overpasses")
                continue
            print(f"  {channel[:4]} n={block['n_pre']} "
                  f"(last pre {block['last'][:10]}, "
                  f"{block['gap_before_event_days']} d before the event)")
            for dim in DIMS:
                cell = block["level"][dim]
                if not cell.get("available"):
                    print(f"    {dim}  late level n/a")
                    continue
                print(f"    {dim}  delta {cell['delta']:+.3f} "
                      f"{cell['delta_ci95']}  MW p {_p(cell['mannwhitney_p'])} "
                      f"holm {_p(cell['holm_p'])}  "
                      f"n {cell['n']['baseline']}/{cell['n']['late_pre']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())









