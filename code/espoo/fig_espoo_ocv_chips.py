#!/usr/bin/env python3
"""Espoo on the **chip** — the OCV vector of the 75 support chips.

This is a second instrument beside ``fig_espoo_ocv_orbits.py``, and the two are
not interchangeable:

* the pipeline figure reads the committed export (``coherence_gamma2``,
  ``coherence_masked_pixels``, ...) over 150 overpasses;
* this one reads ``data/espoo/support/espoo_ocv_channels.csv`` -- 75 chips cut
  with the reference rule (support 12 m base + height layover, guard 30 m, ring
  30-200 m, ``ring_complete: true``) -- and carries ``A, D, F, S, P`` and
  ``gamma2_mask`` for the first time at this site.

They agree only loosely: r(log A, log pipeline mask size) = -0.24 and
r(gamma2_mask, coherence_gamma2) = 0.13. Neither is wrong; they measure
different things.

What the chip series is, structurally
-------------------------------------
* **75 of 150 overpasses.** The mast sits on the far range edge of IW1 in two of
the four relative orbits, so the cut is restricted to the two tracks that place
  it inside the swath: **R160 (ASC, 31 dates)** and **R51 (DESC, 44 dates)**.
  The two tracks have different dates -- only **3 days are shared** -- so the
  orbit contrast is **unpaired** here, unlike the 69 same-day pairs of the
  pipeline series.
* **`A` is orbit-dependent by a factor of seven** (ASC median 76, DESC 528).
* **Two dimensions are functions of `A`**: ``D = A / bbox_area`` exactly, and
  ``gamma2_mask`` carries ``n = A`` in its denominator. A drift in ``D`` or
  ``gamma2_mask`` is therefore partly a drift in the mask size, not an
  independent observation.

The tests are the same three as in the pipeline figure -- drift, year repeat
(two 365-day campaign windows), orbit contrast -- with the pairing dropped and
that stated in the report.

Usage
-----
  python3 code/espoo/fig_espoo_ocv_chips.py
  python3 code/espoo/fig_espoo_ocv_chips.py --selftest
  python3 code/espoo/fig_espoo_ocv_chips.py --verify

Outputs
-------
  figures/espoo/fig_espoo_ocv_chips.jpg   (630 dpi JPEG)
  figures/espoo/fig_espoo_ocv_chips.md
  data/espoo/fig_espoo_ocv_chips.json
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
for path in (HERE, os.path.join(REPO, "code", "lumo")):
    if path not in sys.path:
        sys.path.insert(0, path)

import espoo_ocv_stats as st  # noqa: E402
from lumo_ocv_stats import cliffs_delta  # noqa: E402

SCHEMA = "espoo-ocv-chips-v1"
SITE = "espoo"
LABEL = "Espoo mast (OSM 2233589804) — support chips"
TABLE = os.path.join(REPO, "data", SITE, "support", f"{SITE}_ocv_channels.csv")
#: 630 dpi JPEG, like the shared orbit family in code/orbits/
FIG_DPI = 630
OUT_FIG = os.path.join(REPO, "figures", SITE, "fig_espoo_ocv_chips.jpg")
OUT_MD = os.path.join(REPO, "figures", SITE, "fig_espoo_ocv_chips.md")
OUT_JSON = os.path.join(REPO, "data", SITE, "fig_espoo_ocv_chips.json")

DIMS = {
    "A": ("A", "masked cells (echo area)", False),
    "D": ("D", "density (A / bbox)", False),
    "F": ("F", "fragmentation (components)", False),
    "S": ("S", "shift (centroid to peak)", False),
    "P": ("P", "persistence (pixel frequency)", False),
    "gamma2_mask": ("gamma2_mask", "echo coherence", False),
}
ORBITS = (("ASCENDING", "ASC R160", "#2b6cb0"),
          ("DESCENDING", "DESC R51", "#c53030"))
BLOCK = 5
BOOT = 2000
SEED = 13


# ---------------------------------------------------------------------------
# Data and the three tests
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
    """One dict per chip: date, orbit, and the six chip dimensions."""
    table = table or TABLE
    rows = []
    with open(table, newline="") as handle:
        for raw in csv.DictReader(handle):
            stamp = (raw.get("acquisition_date") or raw.get("acquisition_ts") or "")[:10]
            when = dt.date.fromisoformat(stamp)
            rows.append({"date": when, "ordinal": when.toordinal(),
                         "year": when.year, "month": when.month,
                         "orbit": raw["orbit_direction"],
                         "values": {name: _number(raw, column)
                                    for name, (column, _l, _s) in DIMS.items()}})
    rows.sort(key=lambda row: (row["ordinal"], row["orbit"]))
    return rows


def series_for(rows, orbit, dim):
    out = [[row["ordinal"], row["values"][dim]] for row in rows
           if row["orbit"] == orbit and row["values"][dim] is not None]
    out.sort()
    return out


def drift(rows, boot=BOOT):
    """Test 1: the ported stationarity block plus a ties diagnostic."""
    out = {}
    for dim in DIMS:
        out[dim] = {}
        for orbit, name, _colour in ORBITS:
            series = series_for(rows, orbit, dim)
            if len(series) < 12:
                out[dim][orbit] = {"available": False, "n": len(series)}
                continue
            block = st.stationarity_of(series, block=BLOCK, n_boot=boot)
            block["available"] = True
            block["orbit_label"] = name
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
    """Test 2: the same calendar month in the two 365-day campaign windows."""
    start = min(row["date"] for row in rows)
    out = {"by_dim": {}}
    if (max(row["date"] for row in rows) - start).days < 700:
        return out
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
                if left and right:
                    pairs.append(float(np.median(right) - np.median(left)))
            if len(pairs) < 6:
                out["by_dim"][dim][orbit] = {"available": False,
                                             "months": len(pairs)}
                continue
            positive = sum(1 for value in pairs if value > 0)
            out["by_dim"][dim][orbit] = {
                "available": True, "orbit_label": name, "months": len(pairs),
                "median_delta": float(np.median(pairs)),
                "ci": bootstrap_median_ci(pairs, n=boot),
                "months_positive": int(positive),
                "months_negative": int(len(pairs) - positive),
            }
    return out


def orbit_contrast(rows, boot=BOOT):
    """Test 3, **unpaired**: ASC against DESC as two groups, not by date.

    The track filter that makes the chip possible puts the two orbits on
    different dates (only three days are shared), so there is no same-day pair
    to difference. The contrast is therefore a two-sample comparison, and the
    report says so instead of pretending otherwise.
    """
    out = {}
    for dim in DIMS:
        left = [row["values"][dim] for row in rows
                if row["orbit"] == "ASCENDING" and row["values"][dim] is not None]
        right = [row["values"][dim] for row in rows
                 if row["orbit"] == "DESCENDING" and row["values"][dim] is not None]
        if len(left) < 8 or len(right) < 8:
            out[dim] = {"available": False, "n": [len(left), len(right)]}
            continue
        # ``cliffs_delta(a, b)`` is positive when b > a; pass DESC as the
        # reference so that a positive delta means "ASC above DESC".
        delta = cliffs_delta(right, left)
        rng = np.random.default_rng(SEED)
        draws = np.empty(boot)
        for draw in range(boot):
            a = np.asarray(left)[rng.integers(0, len(left), len(left))]
            b = np.asarray(right)[rng.integers(0, len(right), len(right))]
            draws[draw] = cliffs_delta(b, a)
        test = st.compare_orbit(left, right) or {}
        out[dim] = {
            "available": True, "n": [len(left), len(right)],
            "median": [float(np.median(left)), float(np.median(right))],
            "delta": delta,
            "delta_ci95": [float(np.percentile(draws, 2.5)),
                           float(np.percentile(draws, 97.5))],
            "mannwhitney_p": test.get("mannwhitney_p"),
            "welch_p": test.get("welch_p"),
            "paired": False,
            "shared_dates": None,
        }
    shared = ({row["date"] for row in rows if row["orbit"] == "ASCENDING"}
              & {row["date"] for row in rows if row["orbit"] == "DESCENDING"})
    for dim in out:
        out[dim]["shared_dates"] = len(shared)
    return out


def analyse(table=None, boot=BOOT):
    rows = load(table)
    table = table or TABLE
    import hashlib
    digest = hashlib.sha256(open(table, "rb").read()).hexdigest()
    by_orbit = {orbit: sum(1 for row in rows if row["orbit"] == orbit)
                for orbit, _name, _colour in ORBITS}
    shared = len({row["date"] for row in rows if row["orbit"] == "ASCENDING"}
                 & {row["date"] for row in rows if row["orbit"] == "DESCENDING"})
    return {
        "schema": SCHEMA, "site": SITE, "label": LABEL,
        "table": {"path": os.path.relpath(table, REPO), "sha256": digest,
                  "n_rows": len(rows)},
        "by_orbit": by_orbit, "shared_dates": shared,
        "dates": [min(row["date"] for row in rows).isoformat(),
                  max(row["date"] for row in rows).isoformat()],
        "dimensions": {name: {"column": column, "label": label}
                       for name, (column, label, _s) in DIMS.items()},
        "drift": drift(rows, boot=boot),
        "year_repeat": year_repeat(rows, boot=boot),
        "orbit_contrast": orbit_contrast(rows, boot=boot),
        "drawing_rule": "the Theil-Sen line is solid only where the ported "
                        "stationarity flag says trend, dashed otherwise",
        "rules": {"block": BLOCK, "boot": int(boot), "seed": SEED,
                  "states": "none: the record carries no damage label",
                  "tracks": "ASC = relative orbit 160, DESC = 51; the other "
                            "two tracks put the mast on the far range edge",
                  "pairing": "none: the tracks have 3 shared dates",
                  "tautology": "D = A / bbox_area exactly; gamma2_mask "
                               "contains n = A in its denominator",
                  "instrument": "chip values; not the pipeline columns "
                                "(r(log A, log coherence_masked_pixels) = -0.24)"},
    }


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
    fig, axes = plt.subplots(3, 2, figsize=(11.0, 9.6), dpi=110)
    axes = axes.ravel()
    for index, dim in enumerate(DIMS):
        ax = axes[index]
        for orbit, name, colour in ORBITS:
            series = series_for(rows, orbit, dim)
            if not series:
                continue
            dates = [dt.date.fromordinal(int(day)) for day, _v in series]
            values = [value for _d, value in series]
            ax.plot(dates, values, linestyle="none", marker="o", ms=2.8,
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
            ci = block.get("slope_ci95_per_year") or [float("nan"), float("nan")]
            ax.text(0.02, 0.97,
                    f"{name}: slope {slope:+.3g}/a [{ci[0]:+.3g}, {ci[1]:+.3g}]\n"
                    f"MK p {block['mann_kendall'].get('p')!s:.3} "
                    f"{'trend' if flagged else 'no trend'}",
                    transform=ax.transAxes, va="top", ha="left", fontsize=4.6,
                    color=colour)
        ax.set_title(f"{dim} — {DIMS[dim][1]}", fontsize=7.2)
        ax.tick_params(labelsize=5.0)
        ax.xaxis.set_major_locator(mdates.AutoDateLocator())
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m"))
        ax.grid(alpha=0.18, lw=0.4)
        ax.legend(fontsize=4.6, loc="lower right", frameon=False)
    fig.suptitle(
        f"{report['label']} — {report['table']['n_rows']} chips, two tracks, "
        "no states\n"
        f"{report['dates'][0]} … {report['dates'][1]}; ASC = relative orbit 160, "
        f"DESC = 51; {report['shared_dates']} shared dates\n"
        "solid trend = ported flag rule says trend, dashed = it does not",
        fontsize=7.2)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(fig_path, dpi=FIG_DPI, facecolor="white")
    plt.close(fig)
    return fig_path


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
    add(f"# {report['label']} — the OCV vector of the support chips")
    add("")
    add(f"* table `{report['table']['path']}` (sha256 "
        f"`{report['table']['sha256'][:16]}…`), {report['table']['n_rows']} chips")
    add(f"* {report['dates'][0]} … {report['dates'][1]}; orbits: "
        + ", ".join(f"{k} {v}" for k, v in sorted(report['by_orbit'].items()))
        + f"; **{report['shared_dates']} shared dates**")
    add("* **second instrument**: these are chip values, not the pipeline columns "
        "of `fig_espoo_ocv_orbits.py` (r(log A, log `coherence_masked_pixels`) = "
        "−0.24, r(`gamma2_mask`, `coherence_gamma2`) = 0.13)")
    add("* **two tracks**: ASC = relative orbit 160, DESC = 51; the other two "
        "tracks place the mast on the far range edge, where the 200 m ring runs "
        "past the swath")
    add("* **no pairing**: the tracks have three shared dates, so the orbit "
        "contrast below is a two-sample comparison")
    add("* **tautology**: `D = A / bbox_area` exactly, and `gamma2_mask` carries "
        "`n = A` in its denominator; `A` differs by a factor 7 between the "
        "orbits, so both dimensions are dominated by the mask size")
    add("")
    add("## 1. Drift across the two years")
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
            flags = list(block.get("flags") or [])
            trend = "yes" if "trend" in flags else "no"
            serial = ", ".join(f for f in flags if f != "trend") or "none"
            split = (block.get("split_half") or {}).get("delta")
            add(f"| `{dim}` | {name} | {block['n']} | "
                f"{_f(block['theil_sen_per_year'])} | "
                f"{_ci(block['slope_ci95_per_year'])} | "
                f"{_p(block['mann_kendall'].get('p'))} | "
                f"{_f(block.get('lag1_autocorr'), 2)} | {_f(split, 2)} | "
                f"{trend} | {serial} | "
                f"{block.get('top_value_share', 0):.2f} at {block.get('top_value')}"
                + (" **(degenerate)**" if block.get("degenerate") else "") + " |")
    add("")
    add("## 2. Year repeat — the same calendar month in both campaign years")
    add("")
    add("| dimension | orbit | months | median Δ (y2 − y1) | 95% CI | months + / − |")
    add("|---|---|---|---|---|---|")
    for dim in DIMS:
        for orbit, name, _colour in ORBITS:
            block = report["year_repeat"]["by_dim"][dim][orbit]
            if not block.get("available"):
                add(f"| `{dim}` | {name} | {block.get('months')} | — | — | zu wenige "
                    "Monatspaare |")
                continue
            add(f"| `{dim}` | {name} | {block['months']} | {_f(block['median_delta'])} "
                f"| {_ci(block['ci'])} | {block['months_positive']} / "
                f"{block['months_negative']} |")
    add("")
    add("## 3. Orbit contrast, unpaired (ASC R160 vs DESC R51)")
    add("")
    add("| dimension | n ASC/DESC | median ASC | median DESC | Cliff's δ | 95% CI | "
        "Mann-Whitney p |")
    add("|---|---|---|---|---|---|---|")
    for dim in DIMS:
        block = report["orbit_contrast"][dim]
        if not block.get("available"):
            add(f"| `{dim}` | {block.get('n')} | — | — | — | — | — |")
            continue
        add(f"| `{dim}` | {block['n'][0]}/{block['n'][1]} | "
            f"{_f(block['median'][0], 3, False)} | {_f(block['median'][1], 3, False)} | "
            f"{_f(block['delta'], 2)} | {_ci(block['delta_ci95'], 2)} | "
            f"{_p(block['mannwhitney_p'])} |")
    add("")
    add("## What this can and cannot say")
    add("")
    add("* **It is a different instrument from the pipeline figure.** The two "
        "must not be mixed: they correlate at −0.24 (`A`) and 0.13 (`gamma2`).")
    add("* **The series is halved and selected by track.** 75 of 150 overpasses; "
        "the excluded tracks (ASC R58, DESC R153) are named in the site README.")
    add("* **No paired orbit contrast.** Three shared dates; the contrast is a "
        "two-sample comparison and is reported as such.")
    add("* **`D` and `gamma2_mask` are functions of `A`**, which itself differs "
        "by a factor 7 between the orbits.")
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
    return fig_path, md_path, json_path


def verify(report=None, table=None):
    if report is None:
        with open(OUT_JSON) as handle:
            report = json.load(handle)
    table = table or TABLE
    import hashlib
    if report["table"]["sha256"] != hashlib.sha256(open(table, "rb").read()).hexdigest():
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
                failures.append(f"{dim}/{orbit}: n")
            if abs(float(np.median(y)) - block["median"]) > 1e-9:
                failures.append(f"{dim}/{orbit}: median")
            if slope is None or abs(slope * 365.25 - (block["theil_sen_per_year"] or 0)) > 1e-9:
                failures.append(f"{dim}/{orbit}: slope")
        contrast = report["orbit_contrast"][dim]
        if contrast.get("available"):
            left = [r["values"][dim] for r in rows if r["orbit"] == "ASCENDING"
                    and r["values"][dim] is not None]
            right = [r["values"][dim] for r in rows if r["orbit"] == "DESCENDING"
                     and r["values"][dim] is not None]
            checked += 1
            if abs(cliffs_delta(right, left) - contrast["delta"]) > 1e-9:
                failures.append(f"{dim}: orbit delta")
    print(f"  {checked} field(s) recomputed from {os.path.relpath(table, REPO)}")
    if failures:
        print(f"FAILED — {len(failures)} mismatch(es): {failures[:5]}")
        return 1
    print("  all recomputed values match the JSON")
    return 0


def _synthetic(kind, n_days=80, seed=5):
    rng = np.random.default_rng(seed)
    start = dt.date(2024, 9, 5)
    rows = []
    for index in range(n_days):
        day = start + dt.timedelta(days=12 * index)
        drift = 0.01 * index if kind == "drift" else 0.0
        year = 0.6 if (kind == "year" and (day - start).days >= 365) else 0.0
        for orbit in ("ASCENDING", "DESCENDING"):
            offset = 3.0 if orbit == "ASCENDING" else 0.0
            base = 10.0 + drift + year + offset + rng.normal(0, 0.1)
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

    print("\n1. Flat, drifting, year shift")
    flat = _synthetic("flat")
    block = st.stationarity_of(series_for(flat, "ASCENDING", "A"), block=5, n_boot=300)
    check("a flat series has no trend flag", not block["slope_excludes_zero"],
          f"{block['theil_sen_per_year']:+.3g}")
    drifting = _synthetic("drift")
    block = st.stationarity_of(series_for(drifting, "ASCENDING", "A"),
                               block=5, n_boot=300)
    check("a planted drift is seen", block["slope_excludes_zero"],
          f"{block['theil_sen_per_year']:+.3g}")
    shifted = _synthetic("year")
    repeat = year_repeat(shifted, boot=300)
    cell = repeat["by_dim"]["A"]["ASCENDING"]
    check("the year shift is caught by the repeat test",
          cell.get("available") and cell["median_delta"] > 0.2,
          f"{cell.get('median_delta')}")

    print("\n2. The unpaired orbit contrast")
    flat = _synthetic("flat")
    contrast = orbit_contrast(flat, boot=200)["A"]
    check("the planted ASC offset comes out positive",
          contrast["delta"] > 0.5, f"delta {contrast['delta']:+.2f}")
    check("and it is marked as unpaired", contrast["paired"] is False)

    print("\n3. The real table")
    report = analyse(boot=100)
    check("75 chips are read", report["table"]["n_rows"] == 75,
          str(report["table"]["n_rows"]))
    check("the track split is reported",
          report["by_orbit"].get("ASCENDING") == 31
          and report["by_orbit"].get("DESCENDING") == 44,
          str(report["by_orbit"]))
    check("only three dates are shared", report["shared_dates"] == 3,
          str(report["shared_dates"]))
    text = to_markdown(report)
    check("the report names the second instrument", "second instrument" in text)
    check("and the tautology", "tautology" in text)

    print()
    if failures:
        print(f"FAILED — {len(failures)} check(s): {failures}")
        return 1
    print("all checks passed")
    return 0


def build_parser():
    parser = argparse.ArgumentParser(
        description="Espoo on the support chips: drift, year repeat, unpaired "
                    "orbit contrast")
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
    print(f"{LABEL} — {LABEL and 'the OCV vector of the support chips'}")
    rows = load(args.table)
    report = analyse(args.table, boot=args.boot)
    write_outputs(report, args.fig, args.md, args.json, quiet=args.quiet)
    make_figure(report, rows, args.fig)
    if not args.quiet:
        print(f"  {report['table']['n_rows']} chips, orbits {report['by_orbit']}, "
              f"{report['shared_dates']} shared dates")
        for dim in DIMS:
            line = f"  {dim:11s}"
            for orbit, _name, _colour in ORBITS:
                block = report["drift"][dim][orbit]
                if not block.get("available"):
                    line += "  n too small"
                    continue
                slope = block["theil_sen_per_year"]
                text = (f"{slope:+.2e}" if slope and abs(slope) < 1e-3
                        else f"{slope:+.3g}")
                flag = "trend" if "trend" in (block.get("flags") or []) else "flat "
                line += f"  {orbit[:4]} {text}/a {flag}"
            print(line)
        print(f"  wrote {os.path.relpath(args.fig or OUT_FIG, REPO)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
