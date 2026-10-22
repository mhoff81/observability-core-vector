#!/usr/bin/env python3
"""Espoo over two years, without any state — and the question is stability.

Espoo is the event-free control of this project: no collapse, no rebuild, no
damage label (`damage_label` and `structural_state` are empty in every row of
the committed table). What the record does have is **two full years**
(2024-09-05 … 2026-09-07), a regular six-rows-per-month cadence, and two orbits
whose passes are the same days: 69 dates carry both an afternoon ASC and a
morning DESC.

So the figure draws trajectories only: one panel per dimension, two series (ASC
and DESC), points per overpass, no bands, no states, no smoothing. Each panel
carries the trend as a **Theil-Sen line**, drawn solid only when its moving-block
bootstrap interval excludes zero and dashed when it does not — the line is a
statement with an interval, not a smoother.

Three stability tests, pre-registered:

1. **drift** per dimension and orbit — Theil-Sen slope, block-bootstrap CI,
   Mann-Kendall, lag-1, runs, Ljung-Box and a split-half comparison (the ported
   ``espoo_ocv_stats.stationarity_of`` block, reused rather than re-written);
2. **year repeat** — the same calendar month in year 1 against year 2, twelve
   paired months: stability here means *no systematic year difference*, which is
   a different statement from *no trend*;
3. **orbit contrast, paired by date** (n = 69) — the ASC minus DESC difference on
   the same day: the geometry effect with no event in it, and whether *that*
   effect is itself stable over the two years.

What this figure cannot say: Espoo has **no mask layer** — `A` is the pipeline's
`coherence_masked_pixels`, and `D, F, S, P` and `gamma2_mask` do not exist for
this site (see `data/espoo/asset_evidence.json`). The stability statement is
about coherence, phase, displacement and brightness, not about the mask
geometry. A chip-sized version needs a confirmed `base_width_m` first, which is
still the unconfirmed registry default.

Usage
-----
  python3 code/espoo/fig_espoo_ocv_orbits.py
  python3 code/espoo/fig_espoo_ocv_orbits.py --boot 400 --selftest
  python3 code/espoo/fig_espoo_ocv_orbits.py --verify

Outputs
-------
  figures/espoo/fig_espoo_ocv_orbits.jpg   (630 dpi JPEG)
  figures/espoo/fig_espoo_ocv_orbits.md
  data/espoo/fig_espoo_ocv_orbits.json
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

import espoo_ocv_stats as st  # noqa: E402

SCHEMA = "espoo-ocv-orbits-v1"
SITE = "espoo"
LABEL = "Espoo mast (OSM 2233589804, event-free control)"
TABLE = os.path.join(REPO, "data", SITE, "espoo_ocv_channels.csv")
#: 630 dpi JPEG, like the shared orbit family in code/orbits/
FIG_DPI = 630
OUT_FIG = os.path.join(REPO, "figures", SITE, "fig_espoo_ocv_orbits.jpg")
OUT_MD = os.path.join(REPO, "figures", SITE, "fig_espoo_ocv_orbits.md")
OUT_JSON = os.path.join(REPO, "data", SITE, "fig_espoo_ocv_orbits.json")

#: dimension -> (column, axis label, log scale)
DIMS = {
    "gamma2": ("coherence_gamma2", "echo coherence", False),
    "A": ("coherence_masked_pixels", "masked pixels", False),
    "displacement": ("displacement_los_m", "displacement LOS (m)", False),
    "brightness": ("brightness_ratio", "brightness ratio", False),
    "phase_coherence": ("phase_coherence", "phase coherence", False),
    "phase_rms": ("phase_rms_rad", "phase rms (rad)", False),
}
ORBITS = (("ASCENDING", "ASC afternoon", "#2b6cb0"),
          ("DESCENDING", "DESC morning", "#c53030"))
BLOCK = 5
BOOT = 2000
SEED = 13


# ---------------------------------------------------------------------------
# Data and the three stability tests
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
    """One dict per overpass: date, orbit and the six dimensions."""
    table = table or TABLE
    rows = []
    with open(table, newline="") as handle:
        for raw in csv.DictReader(handle):
            when = dt.date.fromisoformat((raw.get("acquisition_date") or "").strip()
                                         or (raw["acquisition_ts"][:10]))
            values = {name: _number(raw, column) for name, (column, _l, _s)
                      in DIMS.items()}
            rows.append({"date": when, "ordinal": when.toordinal(),
                         "year": when.year, "month": when.month,
                         "orbit": raw["orbit_direction"], "values": values})
    rows.sort(key=lambda row: (row["ordinal"], row["orbit"]))
    return rows


def series_for(rows, orbit, dim):
    """``[[day_ordinal, value], ...]`` for one orbit and dimension."""
    out = [[row["ordinal"], row["values"][dim]] for row in rows
           if row["orbit"] == orbit and row["values"][dim] is not None]
    out.sort()
    return out


def drift(rows, boot=BOOT):
    """Test 1: the ported stationarity block, per dimension and orbit."""
    out = {}
    for dim in DIMS:
        out[dim] = {}
        for orbit, name, _colour in ORBITS:
            series = series_for(rows, orbit, dim)
            if len(series) < 12:
                out[dim][orbit] = {"n": len(series), "available": False}
                continue
            block = st.stationarity_of(series, block=BLOCK, n_boot=boot)
            block["available"] = True
            block["orbit_label"] = name
            # A censored series (brightness capped at 1.0, displacement floored at
            # 0.0) forces the Theil-Sen slope to zero and makes its bootstrap
            # interval degenerate. That is "no power", not "stable", so the share
            # of the single most frequent value travels with the row.
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
    """Moving-block bootstrap CI of the median of a short ordered series."""
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
    """Test 2: the same calendar month in the two years.

    A month pair is ``median(year 2) - median(year 1)`` for that month; the sign
    test counts how many of the twelve months move the same way, and the interval
    is a bootstrap over the twelve month differences.
    """
    start = min(row["date"] for row in rows)
    out = {"window_days": 365, "by_dim": {}}
    if (max(row["date"] for row in rows) - start).days < 700:
        return out
    #: campaign year 0 = the first 365 days, year 1 = the next 365. Calendar years
    #: would pair September 2024 with September 2026 and leave one month over.
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
                if len(left) >= 1 and len(right) >= 1:
                    pairs.append({"month": month,
                                  "delta": float(np.median(right) - np.median(left)),
                                  "n": [len(left), len(right)]})
            if len(pairs) < 8:
                dates = [row["date"] for row in rows
                         if row["orbit"] == orbit and row["values"][dim] is not None]
                out["by_dim"][dim][orbit] = {
                    "available": False, "months": len(pairs),
                    "reason": ("only " + str(len(pairs)) + " calendar month(s) "
                               "carry a value in both campaign years"),
                    "date_range": ([min(dates).isoformat(), max(dates).isoformat()]
                                   if dates else None)}
                continue
            deltas = [pair["delta"] for pair in pairs]
            positive = sum(1 for value in deltas if value > 0)
            out["by_dim"][dim][orbit] = {
                "available": True, "orbit_label": name, "months": len(deltas),
                "median_delta": float(np.median(deltas)),
                "ci": bootstrap_median_ci(deltas, n=boot),
                "months_positive": int(positive),
                "months_negative": int(len(deltas) - positive),
                "pairs": pairs,
            }
    return out


def orbit_pairs(rows, boot=BOOT):
    """Test 3: ASC minus DESC on the same day -- the geometry effect, no event."""
    out = {}
    for dim in DIMS:
        by_date = {}
        for row in rows:
            if row["values"][dim] is None:
                continue
            by_date.setdefault(row["date"], {})[row["orbit"]] = row["values"][dim]
        pairs = [(date, values["ASCENDING"], values["DESCENDING"])
                 for date, values in sorted(by_date.items())
                 if "ASCENDING" in values and "DESCENDING" in values]
        if len(pairs) < 12:
            out[dim] = {"available": False, "n": len(pairs)}
            continue
        diffs = [left - right for _date, left, right in pairs]
        positive = sum(1 for value in diffs if value > 0)
        series = [[date.toordinal(), left - right] for date, left, right in pairs]
        stationarity = st.stationarity_of(series, block=BLOCK, n_boot=boot)
        out[dim] = {
            "available": True, "n": len(diffs),
            "median_difference": float(np.median(diffs)),
            "ci": bootstrap_median_ci(diffs, n=boot),
            "share_positive": float(positive / len(diffs)),
            "drift_of_difference_per_year": stationarity["theil_sen_per_year"],
            "drift_ci95": stationarity["slope_ci95_per_year"],
            "drift_excludes_zero": stationarity["slope_excludes_zero"],
            "flags": stationarity["flags"],
        }
    return out


def analyse(table=None, boot=BOOT):
    rows = load(table)
    table = table or TABLE
    by_orbit = {orbit: sum(1 for row in rows if row["orbit"] == orbit)
                for orbit, _name, _colour in ORBITS}
    return {
        "schema": SCHEMA, "site": SITE, "label": LABEL,
        "table": {"path": os.path.relpath(table, REPO),
                  "sha256": st_hashing(table), "n_rows": len(rows)},
        "by_orbit": by_orbit,
        "dates": [min(row["date"] for row in rows).isoformat(),
                  max(row["date"] for row in rows).isoformat()],
        "dimensions": {name: {"column": column, "label": label}
                       for name, (column, label, _log) in DIMS.items()},
        "drift": drift(rows, boot=boot),
        "year_repeat": year_repeat(rows, boot=boot),
        "orbit_pairs": orbit_pairs(rows, boot=boot),
        "drawing_rule": "the Theil-Sen line is solid only where the ported "
                        "stationarity flag says trend (block interval excludes "
                        "zero *and* Mann-Kendall p < 0.05), dashed otherwise",
        "rules": {"block": BLOCK, "boot": int(boot), "seed": SEED,
                  "states": "none: the record carries no damage or structural "
                            "label, so the figure has no state bands",
                  "chip": "none: Espoo has no mask layer, so A is the pipeline's "
                          "coherence_masked_pixels and D/F/S/P do not exist",
                  "phase_n": "phase_coherence covers 80 of 150 rows"},
    }


def st_hashing(path):
    import hashlib
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


# ---------------------------------------------------------------------------
# The figure
# ---------------------------------------------------------------------------
def make_figure(report, rows, fig_path=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.dates as mdates
    import matplotlib.pyplot as plt

    fig_path = fig_path or OUT_FIG
    os.makedirs(os.path.dirname(os.path.abspath(fig_path)), exist_ok=True)
    names = list(DIMS)
    fig, axes = plt.subplots(3, 2, figsize=(11.0, 9.6), dpi=110)
    axes = axes.ravel()
    for index, dim in enumerate(names):
        ax = axes[index]
        label = DIMS[dim][1]
        for orbit, name, colour in ORBITS:
            series = series_for(rows, orbit, dim)
            if not series:
                continue
            dates = [dt.date.fromordinal(int(day)) for day, _v in series]
            values = [value for _d, value in series]
            ax.plot(dates, values, linestyle="none", marker="o", ms=2.6,
                    mfc=colour, mec="white", mew=0.25, alpha=0.85,
                    label=f"{name} (n={len(values)})")
            block = report["drift"][dim].get(orbit) or {}
            slope = block.get("theil_sen_per_year")
            if slope is None:
                continue
            x = np.array([day for day, _v in series], float)
            y = np.array(values, float)
            per_day = slope / 365.25
            intercept = float(np.median(y - per_day * x))
            xs = np.array([x.min(), x.max()])
            flagged = "trend" in (block.get("flags") or [])
            ax.plot([dt.date.fromordinal(int(v)) for v in xs],
                    intercept + per_day * xs, color=colour, lw=1.2,
                    ls="-" if flagged else "--", alpha=0.9, zorder=3)
            ci = block.get("slope_ci95_per_year")
            verdict = "trend" if "trend" in (block.get("flags") or []) else "no trend"
            ax.text(0.02, 0.97,
                    f"{name}: slope {slope:+.3g}/a "
                    f"[{ci[0]:+.3g}, {ci[1]:+.3g}]\n"
                    f"MK p {block['mann_kendall'].get('p')!s:.3} · lag-1 "
                    f"{block.get('lag1_autocorr')!s:.2} · {verdict}",
                    transform=ax.transAxes, va="top", ha="left", fontsize=4.6,
                    color=colour)
        ax.set_title(f"{dim} — {label}", fontsize=7.2)
        ax.tick_params(labelsize=5.0)
        ax.xaxis.set_major_locator(mdates.AutoDateLocator())
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m"))
        ax.grid(alpha=0.18, lw=0.4)
        ax.legend(fontsize=4.6, loc="lower right", frameon=False)
    fig.suptitle(
        f"{report['label']} — two years, two orbits, no states\n"
        f"{report['table']['n_rows']} overpasses, "
        f"{report['dates'][0]} … {report['dates'][1]}; "
        "solid trend = interval excludes zero, dashed = it does not",
        fontsize=7.4)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(fig_path, dpi=FIG_DPI, facecolor="white")
    plt.close(fig)
    return fig_path


# ---------------------------------------------------------------------------
# Report, outputs, checks
# ---------------------------------------------------------------------------
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


def to_markdown(report):
    lines = []
    add = lines.append
    add(f"# {report['label']} — two years, two orbits, no states")
    add("")
    add(f"* table `{report['table']['path']}` (sha256 "
        f"`{report['table']['sha256'][:16]}…`), {report['table']['n_rows']} overpasses")
    add(f"* {report['dates'][0]} … {report['dates'][1]}; orbits: "
        + ", ".join(f"{k} {v}" for k, v in sorted(report['by_orbit'].items())))
    add("* **no states**: the record carries no damage or structural label, so the "
        "figure has no state bands — the axis is time only")
    add("* **no chip**: Espoo has no mask layer; `A` is the pipeline's "
        "`coherence_masked_pixels`, and `D, F, S, P, gamma2_mask` do not exist here")
    add(f"* drawing rule: {report['drawing_rule']}")
    add("")
    add("## 1. Drift inside the two years")
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
            split = (block.get("split_half") or {}).get("delta")
            flags = list(block.get("flags") or [])
            trend = "yes" if "trend" in flags else "no"
            serial = ", ".join(f for f in flags if f != "trend") or "none"
            add(f"| `{dim}` | {name} | {block['n']} | "
                f"{_f(block['theil_sen_per_year'])} | "
                f"{_ci(block['slope_ci95_per_year'])} | "
                f"{_p(block['mann_kendall'].get('p'))} | "
                f"{_f(block.get('lag1_autocorr'), 2)} | {_f(split, 2)} | "
                f"{trend} | {serial} | "
                f"{block.get('top_value_share', 0):.2f} at {block.get('top_value')}"
                + (" **(degenerate)**" if block.get("degenerate") else "") + " |")
    add("")
    add("## 2. Year repeat — the same calendar month in both years")
    add("")
    add("| dimension | orbit | months | median Δ (y2 − y1) | 95% CI | "
        "months + / − or the reason |")
    add("|---|---|---|---|---|---|")
    for dim in DIMS:
        for orbit, name, _colour in ORBITS:
            block = report["year_repeat"]["by_dim"][dim][orbit]
            if not block.get("available"):
                span = block.get("date_range") or ["—", "—"]
                add(f"| `{dim}` | {name} | {block.get('months')} | — | — | "
                    f"{block.get('reason', '')} ({span[0]} … {span[1]}) |")
                continue
            add(f"| `{dim}` | {name} | {block['months']} | "
                f"{_f(block['median_delta'])} | {_ci(block['ci'])} | "
                f"{block['months_positive']} / {block['months_negative']} |")
    add("")
    add("## 3. Orbit contrast, paired by date (no event in it)")
    add("")
    add("| dimension | pairs | median ASC − DESC | 95% CI | share positive | "
        "drift of that difference |")
    add("|---|---|---|---|---|---|")
    for dim in DIMS:
        block = report["orbit_pairs"][dim]
        if not block.get("available"):
            add(f"| `{dim}` | {block.get('n')} | — | — | — | — |")
            continue
        drift = ("drift" if block["drift_excludes_zero"] else "stable")
        add(f"| `{dim}` | {block['n']} | {_f(block['median_difference'])} | "
            f"{_ci(block['ci'])} | {block['share_positive']:.2f} | "
            f"{_f(block['drift_of_difference_per_year'])} "
            f"{_ci(block['drift_ci95'])} → {drift} |")
    add("")
    add("## What this can and cannot say")
    add("")
    add("* **Stability, not a signature.** There is no event here, so nothing in "
        "this file is an effect; it is the baseline the event sites lack.")
    add("* **No mask geometry.** `D, F, S, P` and `gamma2_mask` do not exist for "
        "Espoo (`data/espoo/asset_evidence.json`), so the stability statement is "
        "about coherence, phase, displacement and brightness. A chip-sized "
        "version needs a confirmed `base_width_m` first.")
    add("* **`phase_coherence` and `phase_rms_rad` cover half the rows** (80 of "
        "150), so their `n` is smaller in every table.")
    add("* **`trend = no` is not `stable`**: the serial flags (autocorrelation, "
        "runs, Ljung-Box) are a different statement, and the split-half and the "
        "year-repeat columns are the second and third check. A dimension can pass "
        "one and fail another.")
    add("* **The solid line follows the ported flag rule**, which is stricter than "
        "the interval alone: it needs the block interval to exclude zero *and* "
        "Mann-Kendall p < 0.05 (`espoo_ocv_stats.stationarity_flags`).")
    add("* **A zero slope with a zero-width interval is not stability.** "
        "`brightness_ratio` is capped at 1.0 and `displacement_los_m` floored at "
        "0.0; where that value fills at least half the rows the pairwise slopes "
        "are all zero and the bootstrap can only return zero. Those rows carry "
        "`(degenerate)` in the ties column: no power, not a result.")
    add("")
    return "\n".join(lines)


def write_outputs(report, fig_path=None, md_path=None, json_path=None, quiet=False):
    fig_path = fig_path or OUT_FIG
    md_path = md_path or OUT_MD
    json_path = json_path or OUT_JSON
    for path in (fig_path, md_path, json_path):
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(json_path, "w") as handle:
        json.dump(report, handle, indent=1)
        handle.write(chr(10))
    with open(md_path, "w") as handle:
        handle.write(to_markdown(report))
        handle.write(chr(10))
    if not quiet:
        print(f"  wrote {os.path.relpath(json_path, REPO)}")
        print(f"  wrote {os.path.relpath(md_path, REPO)}")
        print(f"  wrote {os.path.relpath(fig_path, REPO)}")
    return fig_path, md_path, json_path


def verify(report=None, table=None):
    """Recompute the deterministic numbers (n, medians, Theil-Sen slope)."""
    if report is None:
        with open(OUT_JSON) as handle:
            report = json.load(handle)
    table = table or TABLE
    if report["table"]["sha256"] != st_hashing(table):
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
                failures.append(f"{dim}/{orbit}: n {y.size} != {block['n']}")
            if abs(float(np.median(y)) - block["median"]) > 1e-9:
                failures.append(f"{dim}/{orbit}: median")
            if slope is None or abs(slope * 365.25 - (block["theil_sen_per_year"] or 0)) > 1e-9:
                failures.append(f"{dim}/{orbit}: slope")
        pairs = report["orbit_pairs"][dim]
        if pairs.get("available"):
            by_date = {}
            for row in rows:
                if row["values"][dim] is not None:
                    by_date.setdefault(row["date"], {})[row["orbit"]] = row["values"][dim]
            diffs = [values["ASCENDING"] - values["DESCENDING"]
                     for values in by_date.values()
                     if "ASCENDING" in values and "DESCENDING" in values]
            checked += 1
            if abs(float(np.median(diffs)) - pairs["median_difference"]) > 1e-9:
                failures.append(f"{dim}: paired median difference")
    print(f"  {checked} field(s) recomputed from {os.path.relpath(table, REPO)}")
    if failures:
        print(f"FAILED — {len(failures)} mismatch(es): {failures[:5]}")
        return 1
    print("  all recomputed values match the JSON")
    return 0


def _synthetic(kind, n_days=120, seed=5):
    """Rows in ``load()``'s shape: flat, drifting, or a year shift."""
    rng = np.random.default_rng(seed)
    rows = []
    start = dt.date(2024, 9, 5)
    for index in range(n_days):
        day = start + dt.timedelta(days=3 * index)
        drift = 0.004 * index if kind == "drift" else 0.0
        year_shift = (0.5 if (kind == "year" and (day - start).days >= 365)
                      else 0.0)
        for orbit in ("ASCENDING", "DESCENDING"):
            offset = 0.25 if orbit == "ASCENDING" else 0.0
            base = 1.0 + drift + year_shift + offset + rng.normal(0, 0.02)
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

    print("\n1. A flat series is called flat")
    rows = _synthetic("flat", n_days=150)
    block = st.stationarity_of(series_for(rows, "ASCENDING", "gamma2"),
                               block=5, n_boot=300)
    check("no drift flag on a flat series", not block["slope_excludes_zero"],
          f"slope {block['theil_sen_per_year']:+.3g}, CI {block['slope_ci95_per_year']}")

    print("\n2. A planted drift is found")
    rows = _synthetic("drift", n_days=150)
    block = st.stationarity_of(series_for(rows, "ASCENDING", "gamma2"),
                               block=5, n_boot=300)
    check("the interval excludes zero", block["slope_excludes_zero"],
          f"slope {block['theil_sen_per_year']:+.3g}")
    check("at the planted size (~+0.49/a)",
          abs(block["theil_sen_per_year"] - 0.487) < 0.15,
          f"{block['theil_sen_per_year']:+.3f}")

    print("\n3. A year shift shows up in the repeat test, not in the drift")
    rows = _synthetic("year", n_days=260)
    repeat = year_repeat(rows, boot=400)
    cell = repeat["by_dim"]["gamma2"]["ASCENDING"]
    check("the year-repeat test sees it",
          cell.get("available") and cell["median_delta"] > 0.2,
          f"median delta {cell.get('median_delta')}")

    print("\n4. The orbit contrast on the same days")
    rows = _synthetic("flat", n_days=150)
    pairs = orbit_pairs(rows, boot=300)["gamma2"]
    check("the planted offset is recovered",
          abs(pairs["median_difference"] - 0.25) < 0.03,
          f"{pairs['median_difference']:+.3f}")
    check("and it has no drift of its own", not pairs["drift_excludes_zero"],
          f"{pairs['drift_of_difference_per_year']:+.3g}")

    print("\n5. Determinism and the report")
    first = bootstrap_median_ci([1.0, 2.0, 3.0, 4.0, 5.0, 6.0], n=200, seed=7)
    again = bootstrap_median_ci([1.0, 2.0, 3.0, 4.0, 5.0, 6.0], n=200, seed=7)
    check("the bootstrap is deterministic", first == again, str(first))
    report = analyse(boot=100)
    text = to_markdown(report)
    check("the report says there are no states", "no states" in text)
    check("and that there is no chip", "no chip" in text)
    check("the three tables are present",
          all(section in text for section in
              ("1. Drift", "2. Year repeat", "3. Orbit contrast")))

    print()
    if failures:
        print(f"FAILED — {len(failures)} check(s): {failures}")
        return 1
    print("all checks passed")
    return 0


def build_parser():
    parser = argparse.ArgumentParser(
        description=("Espoo over two years without states: drift, year repeat "
                     "and the paired orbit contrast"))
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
    print(f"{LABEL} — two years, two orbits, no states")
    rows = load(args.table)
    report = analyse(args.table, boot=args.boot)
    write_outputs(report, args.fig, args.md, args.json, quiet=args.quiet)
    if not args.quiet:
        for dim in DIMS:
            line = f"  {dim:16s}"
            for orbit, _name, _colour in ORBITS:
                block = report["drift"][dim][orbit]
                if not block.get("available"):
                    line += "  n too small"
                    continue
                slope = block["theil_sen_per_year"]
                text_slope = (f"{slope:+.2e}" if slope and abs(slope) < 1e-3
                              else f"{slope:+.3g}")
                flag = "trend" if "trend" in (block.get("flags") or []) else "flat "
                line += f"  {orbit[:4]} {text_slope}/a {flag}"
            print(line)
    make_figure(report, rows, args.fig)
    if not args.quiet:
        print(f"  wrote {os.path.relpath(args.fig or OUT_FIG, REPO)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
