#!/usr/bin/env python3
"""SDS for the OCV dimensions on the support-sized chips — Bautzen and LUMO.

The chip changed from a fixed search window to a support-sized, odd,
anchor-centred one (guard 30 m, reference ring 30..200 m). This script asks the
question that change has to answer: **does the channel still separate the
labelled states, and for which state pair?**

The metric is the project's own, not a re-implementation: ``cliffs_delta`` and
``bootstrap_delta_ci`` are imported from ``code/lumo/lumo_ocv_stats.py``, so the
numbers here are the numbers the LUMO reference JSON was built with

    SDS      = 10 * |Cliff's delta|
    delta    = (#{y > x} - #{y < x}) / (n1 * n2)      ties count for neither
    95% CI   = percentile bootstrap, both groups resampled independently,
               10 000 resamples, seed 7 — one interval *per comparison*

Reporting hierarchy (this is the scientific claim, not a formatting choice)
--------------------------------------------------------------------------
``c`` is the **observation channel**, i.e. the orbit. ASC and DESC see different
geometries, so the primary analysis is per channel and per state pair:

1. ``per_channel_per_state``     — **primary**: SDS_j,c(Z1, Z2) for every pair;
2. ``per_channel_max_pairwise``  — compact summary, and it never appears without
   the state pair that produced it, the signed delta, that pair's CI and n;
3. ``per_channel_pooled_damage`` — healthy versus all damage states per channel;
4. ``pooled_all_damage``         — both channels pooled. This is the basis of the
   committed LUMO reference and is kept as a **compatibility pin and context**,
   never as the primary evidence.

Two wording rules are enforced in code, because both are easy to get wrong:

* the interval belongs to **one comparison**. A maximum over several comparisons
  has no interval of its own, so the report says "95% CI of the winning pair"
  and never "CI of the max SDS" (the writer refuses to emit that string);
* a maximum without its winning pair, signed delta, CI and n1/n2 is an error,
  not a blank (``max_pairwise`` raises).

Usage
-----
  python3 code/sds/sds_support.py --site all
  python3 code/sds/sds_support.py --site lumo --boot 2000 --json /tmp/x.json
  python3 code/sds/sds_support.py --selftest

Outputs (written next to the tables they read):
  data/<site>/support/sds_support.json   the four levels, machine-readable
  data/<site>/support/sds_support.md     the report, primary analysis first
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
for path in (os.path.join(REPO, "code", "lumo"),):
    if path not in sys.path:
        sys.path.insert(0, path)

from lumo_ocv_stats import (BOOT_MAIN, bootstrap_delta_ci,  # noqa: E402
                            cliffs_delta)

SCHEMA = "ocv-sds-support-v2"
SEED = 7

#: The five bands of the 0-10 scale, printed beside every score.
BANDS = ((2.0, "very weak"), (4.0, "weak"), (6.0, "moderate"),
         (8.0, "strong"), (float("inf"), "very strong"))

#: Below this many damaged values a score is flagged; the CI is reported anyway
##and the flag travels with every copy of the number.
SMALL_N = 5

#: The wording the writer must never produce: it would claim an interval for a
## selection over comparisons.
FORBIDDEN_PHRASES = ("CI of max", "CI of the max", "confidence interval of the max")

SITES = {
    "lumo": {
        "label": "LUMO lattice tower (Hannover)",
        "table_key": "lumo",
        "channel_column": "orbit_direction",
        "channels": ("ASCENDING", "DESCENDING"),
        "label_column": "damage_label",
        "healthy": "healthy",
        "states": ("DAM 3", "DAM 4", "DAM 6"),
        "record_class": True,
        "dimensions": (("gamma2", "coherence_gamma2"), ("P", "P"), ("D", "D"),
                       ("F", "F"), ("A", "A"), ("S", "S")),
        "committed": "data/lumo/lumo_ocv_channels.csv",
        "support": "data/lumo/support/lumo_ocv_channels.csv",
        "emv_reference": "data/lumo/reference/lumo_echo_mask_vector.json",
    },
    "cts": {
        "label": "Champlain Towers South (target footprint)",
        "table_key": "cts",
        "channel_column": "orbit_direction",
        "channels": ("ASCENDING",),
        "label_column": "state",
        "healthy": "pre-collapse (healthy)",
        "states": ("post-collapse",),
        "record_class": False,
        "dimensions": (("gamma2", "gamma2_mask"), ("P", "P"), ("D", "D"),
                       ("F", "F"), ("A", "A"), ("S", "S")),
        "committed": "data/cts/cts_ocv_channels.csv",
        "support": "data/cts/support/cts_ocv_channels.csv",
        "emv_reference": None,
    },
    "kdlo": {
        "label": "KDLO media tower",
        "table_key": "kdlo",
        "channel_column": "orbit_direction",
        "channels": ("DESCENDING",),
        "label_column": "state",
        "healthy": "pre",
        "states": ("rebuild",),
        # `single` carries n = 1: a state with one overpass has no contrast, and
        # scoring it would turn a single comparison into an SDS.
        "exclude_states": ("single",),
        "record_class": False,
        "dimensions": (("gamma2", "gamma2_mask"), ("P", "P"), ("D", "D"),
                       ("F", "F"), ("A", "A"), ("S", "S")),
        "committed": "data/kdlo/kdlo_ocv_channels.csv",
        "support": "data/kdlo/support/kdlo_ocv_channels.csv",
        "emv_reference": None,
    },
    "carola": {
        "label": "Carolabrücke deck",
        "table_key": "carola",
        "channel_column": "orbit_direction",
        "channels": ("ASCENDING", "DESCENDING"),
        "label_column": "state",
        "healthy": "pre-collapse (healthy)",
        "states": ("post-collapse",),
        "record_class": False,
        "dimensions": (("gamma2", "gamma2_mask"), ("P", "P"), ("D", "D"),
                       ("F", "F"), ("A", "A"), ("S", "S")),
        "committed": "data/carola/carola_ocv_channels.csv",
        "support": "data/carola/support/carola_ocv_channels.csv",
        "emv_reference": None,
    },
    "haltern": {
        "label": "Haltern AV9 wind turbine",
        "table_key": "haltern",
        "channel_column": "orbit_direction",
        "channels": ("ASCENDING", "DESCENDING"),
        "label_column": "state",
        "healthy": "pre-collapse (healthy)",
        "states": ("post-collapse",),
        "record_class": False,
        "dimensions": (("gamma2", "gamma2_mask"), ("P", "P"), ("D", "D"),
                       ("F", "F"), ("A", "A"), ("S", "S")),
        "committed": "data/haltern/haltern_ocv_channels.csv",
        "support": "data/haltern/support/haltern_ocv_channels.csv",
        "emv_reference": None,
    },
    "ywf": {
        "label": "Yeongdeok Unit 21",
        "table_key": "ywf",
        "channel_column": "orbit_direction",
        "channels": ("ASCENDING", "DESCENDING"),
        "label_column": "state",
        "healthy": "pre-collapse (healthy)",
        "states": ("post-collapse",),
        "record_class": False,
        "dimensions": (("gamma2", "gamma2_mask"), ("P", "P"), ("D", "D"),
                       ("F", "F"), ("A", "A"), ("S", "S")),
        "committed": "data/ywf/ywf_ocv_channels.csv",
        "support": "data/ywf/support/ywf_ocv_channels.csv",
        "emv_reference": None,
    },
    "morandi": {
        "label": "Ponte Morandi deck",
        "table_key": "morandi",
        "channel_column": "orbit_direction",
        "channels": ("ASCENDING", "DESCENDING"),
        "label_column": "state",
        "healthy": "pre-collapse (healthy)",
        "states": ("post-collapse",),
        "record_class": False,
        "dimensions": (("gamma2", "gamma2_mask"), ("P", "P"), ("D", "D"),
                       ("F", "F"), ("A", "A"), ("S", "S")),
        "committed": "data/morandi/morandi_ocv_channels.csv",
        "support": "data/morandi/support/morandi_ocv_channels.csv",
        "emv_reference": None,
    },
    "bautzen": {
        "label": "Bautzen OpenLabs bridge",
        "table_key": "bautzen",
        "channel_column": "orbit_direction",
        "channels": ("ASC", "DESC"),
        "label_column": "state",
        "healthy": "RS",
        "states": ("DS1", "DS2"),
        "record_class": False,
        "dimensions": (("gamma2", "gamma2_band_floor_corrected"), ("P", "P"),
                       ("D", "D"), ("F", "F"), ("A", "A"), ("S", "S")),
        "committed": "data/bautzen/bautzen_ocv_channels.csv",
        "support": "data/bautzen/support/bautzen_ocv_channels.csv",
        "emv_reference": None,
    },
}


# ---------------------------------------------------------------------------
# Reading a table, and one comparison
# ---------------------------------------------------------------------------
def read_rows(path, spec):
    """The coherence rows of a channel table, with the empty cell kept empty."""
    with open(path, newline="") as handle:
        rows = list(csv.DictReader(handle))
    if spec.get("record_class"):
        rows = [row for row in rows
                if row.get("record_class") in ("", "coherence")]
    return rows


def values(rows, field, label_column, labels):
    """The finite values of one field under one or more state labels."""
    out = []
    for row in rows:
        if row.get(label_column) not in labels:
            continue
        cell = row.get(field)
        if cell in (None, ""):
            continue
        try:
            number = float(cell)
        except (TypeError, ValueError):
            continue
        if math.isfinite(number):
            out.append(number)
    return out


def band(sds):
    """The band name of a 0-10 score (``None`` stays ``None``)."""
    if sds is None:
        return None
    for edge, label in BANDS:
        if sds < edge:
            return label
    return BANDS[-1][1]


def comparison(healthy, damaged, label, n_boot=BOOT_MAIN, seed=SEED):
    """One SDS: delta, score, that comparison's CI, the counts and the flags."""
    delta = cliffs_delta(healthy, damaged)
    sds = None if delta is None else 10.0 * abs(delta)
    interval = (None if delta is None
                else bootstrap_delta_ci(healthy, damaged, n=n_boot, seed=seed))
    low, high = interval if interval else (None, None)
    return {
        "comparison": label,
        "sds": sds,
        "band": band(sds),
        "delta": delta,
        "delta_sign": (None if delta is None
                       else "damaged > healthy" if delta > 0
                       else "damaged < healthy" if delta < 0 else "tie"),
        "ci": {"low": low, "high": high, "level": 0.95,
               "n_boot": n_boot, "seed": seed,
               "scope": "per comparison (this state pair only)"},
        "n_healthy": len(healthy),
        "n_damaged": len(damaged),
        "small_n": len(damaged) < SMALL_N or len(healthy) < SMALL_N,
        "straddles_zero": (None if low is None else bool(low <= 0.0 <= high)),
        "median_healthy": _median(healthy),
        "median_damaged": _median(damaged),
        "reason": (None if delta is not None else
                   f"{len(healthy)} healthy and {len(damaged)} damaged value(s)"),
    }


def _median(numbers):
    if not numbers:
        return None
    ordered = sorted(numbers)
    middle = len(ordered) // 2
    if len(ordered) % 2:
        return float(ordered[middle])
    return float((ordered[middle - 1] + ordered[middle]) / 2.0)


def _states_in(rows, spec):
    """The state labels the table actually carries, checked against the registry.

    A renamed label must fail loudly: silently scoring an empty set would report
    "no discrimination" for a state that is simply spelled differently.
    """
    found = {row.get(spec["label_column"]) for row in rows}
    found.discard(spec["healthy"])
    found = {value for value in found if value not in (None, "")}
    found -= set(spec.get("exclude_states") or ())
    if found != set(spec["states"]):
        raise ValueError(
            f"{spec['label']}: the table carries {sorted(found)} but the "
            f"registry says {sorted(spec['states'])}")
    return list(spec["states"])


def _channel_rows(rows, spec, channel):
    return [row for row in rows
            if row.get(spec["channel_column"]) == channel]


# ---------------------------------------------------------------------------
# The four levels
# ---------------------------------------------------------------------------
def level_pooled_all_damage(rows, spec, n_boot=BOOT_MAIN, seed=SEED):
    """Both channels pooled, healthy versus every damaged row.

    This is the basis of the committed LUMO EMV reference, so it is the level
    that can be pinned. It is **not** the primary evidence: pooling ASC and DESC
    averages two geometries.
    """
    out = {}
    for dim, field in spec["dimensions"]:
        healthy = values(rows, field, spec["label_column"], (spec["healthy"],))
        damaged = values(rows, field, spec["label_column"], spec["states"])
        out[dim] = comparison(healthy, damaged,
                              f"{spec['healthy']} vs all damaged", n_boot, seed)
    return out


def level_per_channel_pooled_damage(rows, spec, n_boot=BOOT_MAIN, seed=SEED):
    """Healthy versus all damage states, one channel at a time."""
    out = {}
    for channel in spec["channels"]:
        subset = _channel_rows(rows, spec, channel)
        out[channel] = {}
        for dim, field in spec["dimensions"]:
            healthy = values(subset, field, spec["label_column"],
                             (spec["healthy"],))
            damaged = values(subset, field, spec["label_column"],
                             spec["states"])
            out[channel][dim] = comparison(
                healthy, damaged, f"{spec['healthy']} vs all damaged",
                n_boot, seed)
    return out


def level_per_channel_per_state(rows, spec, n_boot=BOOT_MAIN, seed=SEED):
    """**The primary analysis**: one SDS per channel, dimension and state pair.

    Nothing is pooled here: a channel is one observation geometry and a state
    pair is one contrast, so each cell answers its own question.
    """
    states = _states_in(rows, spec)
    out = {}
    for channel in spec["channels"]:
        subset = _channel_rows(rows, spec, channel)
        out[channel] = {}
        for dim, field in spec["dimensions"]:
            healthy = values(subset, field, spec["label_column"],
                             (spec["healthy"],))
            per_state = {}
            for state in states:
                damaged = values(subset, field, spec["label_column"], (state,))
                per_state[state] = comparison(
                    healthy, damaged, f"{spec['healthy']} – {state}",
                    n_boot, seed)
            out[channel][dim] = per_state
    return out



def level_per_channel_max_pairwise(per_state):
    """The compact summary. It carries the winning pair and raises without it.

    A maximum over several comparisons is a *selection*: it has no interval of
    its own, so the interval reported here is copied from the winning comparison
    and labelled as such. An entry missing the pair, the signed delta, the
    interval or the counts is a programming error, not a blank cell.
    """
    out = {}
    for channel, dims in per_state.items():
        out[channel] = {}
        for dim, cells in dims.items():
            scored = [(state, cell) for state, cell in cells.items()
                      if cell["sds"] is not None]
            if not scored:
                out[channel][dim] = {
                    "sds": None, "band": None, "winning_pair": None,
                    "reason": "no state pair could be scored in this channel",
                }
                continue
            scored.sort(key=lambda item: item[1]["sds"], reverse=True)
            state, winner = scored[0]
            entry = {
                "sds": winner["sds"],
                "band": winner["band"],
                "winning_pair": winner["comparison"],
                "state": state,
                "delta": winner["delta"],
                "delta_sign": winner["delta_sign"],
                "ci": dict(winner["ci"]),
                "n_healthy": winner["n_healthy"],
                "n_damaged": winner["n_damaged"],
                "small_n": winner["small_n"],
                "straddles_zero": winner["straddles_zero"],
                "median_healthy": winner["median_healthy"],
                "median_damaged": winner["median_damaged"],
                "runner_up": ({"state": scored[1][0],
                               "sds": scored[1][1]["sds"],
                               "winning_pair": scored[1][1]["comparison"]}
                              if len(scored) > 1 else None),
            }
            entry["ci"]["scope"] = ("of the winning pair (this comparison "
                                    "only) - not of the maximum")
            missing = [key for key in ("winning_pair", "delta", "ci",
                                       "n_healthy", "n_damaged")
                       if entry.get(key) is None]
            if missing:
                raise ValueError(
                    f"max pairwise for {channel}/{dim} would be incomplete: "
                    f"{missing}")
            out[channel][dim] = entry
    return out


def pin_emv_reference(reference_path, pooled, repo):
    """Reproduce the committed per-dimension SDS of the LUMO EMV reference.

    ``sds_pooled_all_damaged`` in that JSON was computed on the committed table
    with every row counted; if this script cannot reproduce it exactly, the
    metric here is not the project's metric and the report is void.
    """
    if not reference_path:
        return {"available": False,
                "reason": "the registry names no reference for this site"}
    path = (reference_path if os.path.isabs(reference_path)
            else os.path.join(repo, reference_path))
    if not os.path.exists(path):
        return {"available": False, "reason": f"{path} is not on this machine"}
    document = json.load(open(path))
    checks = []
    for dim in ("A", "D", "F", "S", "P"):
        wanted = ((document.get("per_dimension") or {}).get(dim)
                  or {}).get("sds_pooled_all_damaged")
        got = (pooled.get(dim) or {}).get("sds")
        same = (wanted is not None and got is not None
                and abs(got - wanted) <= 1e-9)
        checks.append({"dimension": dim, "reference": wanted,
                       "recomputed": got, "ok": bool(same)})
    return {"available": True, "path": os.path.relpath(path, repo),
            "checks": checks, "ok": all(check["ok"] for check in checks)}



# ---------------------------------------------------------------------------
# Assembling a site
# ---------------------------------------------------------------------------
def sha256(path):
    """The file's SHA-256, so a report names the bytes it was built from."""
    import hashlib

    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def analyse_site(site, spec, repo, n_boot=BOOT_MAIN, seed=SEED):
    """All four levels for both tables of one site."""
    result = {"site": site, "label": spec["label"], "tables": {},
              "levels": {}, "primary_analysis": "per_channel_per_state"}
    for kind in ("committed", "support"):
        path = os.path.join(repo, spec[kind])
        rows = read_rows(path, spec)
        per_state = level_per_channel_per_state(rows, spec, n_boot, seed)
        levels = {
            "per_channel_per_state": per_state,
            "per_channel_max_pairwise": level_per_channel_max_pairwise(per_state),
            "per_channel_pooled_damage": level_per_channel_pooled_damage(
                rows, spec, n_boot, seed),
            "pooled_all_damage": level_pooled_all_damage(rows, spec, n_boot,
                                                         seed),
        }
        counts = {}
        for channel in spec["channels"]:
            subset = _channel_rows(rows, spec, channel)
            counts[channel] = {
                state: len([row for row in subset
                            if row.get(spec["label_column"]) == state])
                for state in (spec["healthy"],) + tuple(spec["states"])}
        result["tables"][kind] = {
            "path": os.path.relpath(path, repo),
            "sha256": sha256(path),
            "n_rows": len(rows),
            "counts_by_channel": counts,
        }
        result["levels"][kind] = levels
    result["pin"] = pin_emv_reference(spec.get("emv_reference"),
                                      result["levels"]["committed"]["pooled_all_damage"],
                                      repo)
    return result


def analysis(results, repo, n_boot, seed):
    """The whole report: the metric, its scope, and every site's four levels."""
    return {
        "schema": SCHEMA,
        "generated_by": "code/sds/sds_support.py",
        "metric": {
            "definition": ("SDS = 10 * |Cliff's delta|; delta = (#{y>x} - "
                           "#{y<x}) / (n1*n2), ties count for neither side"),
            "library": "code/lumo/lumo_ocv_stats.py (cliffs_delta)",
            "scale": 10.0,
            "bands": [[edge, name] for edge, name in BANDS
                      if math.isfinite(edge)] + [[None, BANDS[-1][1]]],
            "ci": {"method": ("percentile bootstrap, both groups resampled "
                              "independently"),
                   "n_boot": n_boot, "seed": seed, "level": 0.95,
                   "scope": ("one interval per comparison; a maximum over "
                             "comparisons has no interval of its own")},
        },
        "primary_analysis": ("per_channel_per_state - the observation channel "
                             "is the orbit, and ASC/DESC are never pooled into "
                             "the primary number"),
        "levels": {
            "pooled_all_damage": "compatibility pin and context only",
            "per_channel_pooled_damage": "healthy versus all damage states, per channel",
            "per_channel_per_state": "primary analysis",
            "per_channel_max_pairwise": "compact summary, always with its winning pair",
        },
        "sites": {site: analyse_site(site, spec, repo, n_boot, seed)
                  for site, spec in SITES.items()
                  if site in results},
    }



# ---------------------------------------------------------------------------
# The markdown report
# ---------------------------------------------------------------------------
def _num(value, digits=2, signed=False):
    if value is None:
        return "–"
    text = f"{value:+.{digits}f}" if signed else f"{value:.{digits}f}"
    return text


def _ci(cell, digits=3):
    interval = (cell or {}).get("ci") or {}
    low, high = interval.get("low"), interval.get("high")
    if low is None or high is None:
        return "–"
    return f"[{low:+.{digits}f}, {high:+.{digits}f}]"


def _flags(cell):
    flags = []
    if (cell or {}).get("straddles_zero"):
        flags.append("CI includes 0")
    if (cell or {}).get("small_n"):
        flags.append(f"n<{SMALL_N}")
    return ", ".join(flags) or "–"


def _row(cells):
    return "| " + " | ".join(str(cell) for cell in cells) + " |"


def markdown(report, repo):
    """The report. Primary analysis first; the pool is a pin, not evidence."""
    lines = []
    add = lines.append

    add("# SDS of the OCV dimensions on the support-sized chips")
    add("")
    metric = report["metric"]
    ci_meta = metric["ci"]
    add(f"{len(report['sites'])} sites (" + ", ".join(sorted(report['sites'])) +
        "), recomputed on the **support-sized chip** (odd, anchor-centred, "
        "guard 30 m, reference ring 30–200 m) and beside it on "
        "the committed table. The metric is the project's own: "
        "`SDS = 10 · |Cliff's delta|`, computed with "
        "`code/lumo/lumo_ocv_stats.cliffs_delta`; the interval is a percentile "
        f"bootstrap ({ci_meta['n_boot']:,} resamples, seed {ci_meta['seed']}) "
        "of **one comparison**.")
    add("")
    add("## 0. Reading rules")
    add("")
    add("1. **The primary analysis is per channel and per state pair.** The "
        "observation channel is the orbit; ASC and DESC see different "
        "geometries and are never pooled into the primary number.")
    add("2. **A maximum always travels with its winning pair** — the signed "
        "delta, that pair's interval and its two sample sizes.")
    add("3. **The interval belongs to one comparison.** For a maximum over "
        "several comparisons the report prints the *winning pair's* interval "
        "and says so; it is not an interval for the selection.")
    add("4. **`pooled_all_damage` is a compatibility pin and context**, not the "
        "primary evidence: it is the basis the committed LUMO reference used.")
    add("")

    add("## 1. Primary analysis — per channel, per state pair")
    add("")
    add("Support chip. ★ marks the winning comparison of its "
        "dimension (the one the compact table in §2 carries).")
    add("")
    for site, block in report["sites"].items():
        levels = block["levels"]["support"]["per_channel_per_state"]
        add(f"### {block['label']}")
        add("")
        for channel, dims in levels.items():
            counts = (block["tables"]["support"]["counts_by_channel"]
                      .get(channel) or {})
            add(f"**{channel}** — rows per class: "
                + ", ".join(f"{name} {count}"
                            for name, count in counts.items()))
            add("")
            add(_row(["dimension", "state pair", "SDS", "band", "delta",
                      "95% CI", "n1 / n2", "flags"]))
            add(_row(["---"] * 8))
            for dim, cells in dims.items():
                best = max((cell["sds"] for cell in cells.values()
                            if cell["sds"] is not None), default=None)
                for pair, cell in cells.items():
                    star = " ★" if (best is not None
                                    and cell["sds"] == best) else ""
                    add(_row([f"`{dim}`", cell["comparison"] + star,
                              _num(cell["sds"]),
                              cell["band"] or "–",
                              _num(cell["delta"], 3, signed=True),
                              _ci(cell),
                              f"{cell['n_healthy']} / {cell['n_damaged']}",
                              _flags(cell)]))
            add("")
            add("")

    add("## 2. Compact summary — max pairwise")
    add("")
    add("The paper table. Every row names the state pair that produced the "
        "maximum, its signed delta, **that pair's** interval, and the two "
        "sample sizes — a maximum without them would read as a general "
        "discrimination the number does not support.")
    add("")
    add(_row(["asset", "channel", "dimension", "max SDS", "winning pair",
              "delta", "95% CI (winning pair)", "n1 / n2", "flags"]))
    add(_row(["---"] * 9))
    for site, block in report["sites"].items():
        arms = block["levels"]["support"]["per_channel_max_pairwise"]
        for channel, dims in arms.items():
            for dim, entry in dims.items():
                add(_row([site, channel, f"`{dim}`", _num(entry.get("sds")),
                          entry.get("winning_pair") or "–",
                          _num(entry.get("delta"), 3, signed=True),
                          _ci(entry),
                          (f"{entry['n_healthy']} / {entry['n_damaged']}"
                           if entry.get("n_healthy") is not None else "–"),
                          _flags(entry)]))
    add("")

    add("## 3. Compatibility pin and context — pooled over both channels")
    add("")
    add("**Not the primary evidence.** This is the aggregation the committed "
        "LUMO reference used (`sds_pooled_all_damaged`), so it is reproduced "
        "here as a pin and reported for context. Pooling ASC and DESC mixes two "
        "observation geometries.")
    add("")
    add(_row(["asset", "table", "dimension", "SDS", "delta", "95% CI",
              "n1 / n2", "flags"]))
    add(_row(["---"] * 8))
    for site, block in report["sites"].items():
        for kind in ("committed", "support"):
            pooled = block["levels"][kind]["pooled_all_damage"]
            for dim, cell in pooled.items():
                add(_row([site, kind, f"`{dim}`", _num(cell["sds"]),
                          _num(cell["delta"], 3, signed=True), _ci(cell),
                          f"{cell['n_healthy']} / {cell['n_damaged']}",
                          _flags(cell)]))
    add("")
    for site, block in report["sites"].items():
        pin = block.get("pin") or {}
        if not pin.get("available"):
            add(f"- **{site}**: nothing to pin — "
                f"{pin.get('reason', 'not available')}.")
            continue
        status = "reproduced exactly" if pin.get("ok") else "**MISMATCH**"
        add(f"- **{site}**: the committed EMV reference "
            f"(`{pin['path']}`) is {status} for all five dimensions.")
    add("")

    add("## 4. What the chip size changed")
    add("")
    add("Max pairwise, committed → support, with the winning pair on each side. "
        "A changed pair is the point: the maximum moved to a different "
        "comparison, not just to a different number.")
    add("")
    add(_row(["asset", "channel", "dimension", "committed", "support",
              "Δ SDS", "winning pair committed → support"]))
    add(_row(["---"] * 7))
    for site, block in report["sites"].items():
        old = block["levels"]["committed"]["per_channel_max_pairwise"]
        new = block["levels"]["support"]["per_channel_max_pairwise"]
        for channel, dims in new.items():
            for dim, entry in dims.items():
                before = old.get(channel, {}).get(dim, {})
                delta = (None if entry.get("sds") is None
                         or before.get("sds") is None
                         else entry["sds"] - before["sds"])
                add(_row([site, channel, f"`{dim}`", _num(before.get("sds")),
                          _num(entry.get("sds")), _num(delta, 2, signed=True),
                          f"{before.get('winning_pair') or '–'} → "
                          f"{entry.get('winning_pair') or '–'}"]))
    add("")

    add("## 5. What this does not claim")
    add("")
    add("- SDS measures **separability on the given contrast** — two labelled "
        "groups on one record — never damage.")
    add("- A channel with few damaged overpasses gets a wide interval and the "
        f"`n<{SMALL_N}` flag; the score is printed, the flag is not optional.")
    add("- The two tables are two measurements, not a better and a worse "
        "instrument: the chip changes which cells enter the mask.")
    add("- State, season and orbit are confounded in these records; a "
        "per-channel SDS does not remove that.")
    add("")
    text = "\n".join(lines)
    lowered = text.lower()
    for phrase in FORBIDDEN_PHRASES:
        if phrase.lower() in lowered:
            raise ValueError(
                f"the report would claim an interval for a selection: "
                f"{phrase!r} - use 'CI of the winning pair' instead")
    return text




# ---------------------------------------------------------------------------
# Writing
# ---------------------------------------------------------------------------
def _guard_text(text):
    """Refuse to ship a report that puts an interval on a selection."""
    lowered = text.lower()
    for phrase in FORBIDDEN_PHRASES:
        if phrase.lower() in lowered:
            raise ValueError(
                f"the report would claim an interval for a selection: "
                f"{phrase!r} - use 'CI of the winning pair' instead")
    return text


def write_outputs(report, repo, out_dir=None, quiet=False):
    """Write the JSON and the markdown; return the two paths."""
    target = out_dir or os.path.join(repo, "data", "sds")
    os.makedirs(target, exist_ok=True)
    json_path = os.path.join(target, "sds_support.json")
    md_path = os.path.join(target, "sds_support.md")
    with open(json_path, "w") as handle:
        json.dump(report, handle, indent=1, sort_keys=True)
        handle.write("\n")
    text = _guard_text(markdown(report, repo))
    with open(md_path, "w") as handle:
        handle.write(text)
    if not quiet:
        print(f"  wrote {os.path.relpath(json_path, repo)}")
        print(f"  wrote {os.path.relpath(md_path, repo)}")
    return json_path, md_path



# ---------------------------------------------------------------------------
# Self-test — the metric, the flags and the wording rules
# ---------------------------------------------------------------------------
def selftest():
    """Offline checks: the metric by hand, the flags, the guard, the pin."""
    failures = []

    def check(label, condition, detail=""):
        print(f"  [{'ok  ' if condition else 'FAIL'}] {label}"
              + (f" — {detail}" if detail else ""))
        if not condition:
            failures.append(label)

    print("\n1. Cliff's delta, by hand")
    full = comparison([1.0, 2.0, 3.0], [4.0, 5.0], "test", n_boot=200, seed=7)
    check("complete separation gives delta = +1 and SDS 10",
          full["delta"] == 1.0 and full["sds"] == 10.0,
          f"delta {full['delta']}, SDS {full['sds']}")
    check("and the band is the top one", full["band"] == "very strong")
    reverse = comparison([4.0, 5.0], [1.0, 2.0], "test", n_boot=200, seed=7)
    check("the score takes the size, not the direction",
          reverse["sds"] == 10.0 and reverse["delta"] == -1.0,
          f"delta {reverse['delta']}, SDS {reverse['sds']}")
    check("the direction is reported beside it",
          reverse["delta_sign"] == "damaged < healthy")
    tied = comparison([1.0, 2.0], [1.0, 2.0], "test", n_boot=200, seed=7)
    check("a tie counts for neither side: delta = 0",
          tied["delta"] == 0.0 and tied["sds"] == 0.0,
          f"delta {tied['delta']}")

    print("\n2. the interval is one comparison's, and it is reproducible")
    first = comparison([1.0, 2.0, 3.0, 4.0], [5.0, 6.0], "test", n_boot=500,
                       seed=7)
    second = comparison([1.0, 2.0, 3.0, 4.0], [5.0, 6.0], "test", n_boot=500,
                        seed=7)
    check("same seed, same interval", first["ci"] == second["ci"],
          str(first["ci"]["low"]))
    check("the interval is declared per comparison",
          "per comparison" in first["ci"]["scope"])

    print("\n3. a missing class is None with a reason, never a zero")
    empty = comparison([], [1.0, 2.0], "test", n_boot=100, seed=7)
    check("no healthy values -> no score", empty["sds"] is None)
    check("and it says why", bool(empty["reason"]), empty["reason"])

    print("\n4. the small-n flag travels with the number")
    small = comparison([1.0] * 40, [2.0, 3.0, 4.0, 5.0], "test", n_boot=100,
                       seed=7)
    check("four damaged values set small_n", small["small_n"] is True,
          f"n_damaged {small['n_damaged']}")

    print("\n5. a maximum cannot be shipped without its winning pair")
    per_state = {"ASC": {"F": {"DAM 3": comparison([1.0] * 5, [2.0] * 5,
                                                  "healthy - DAM 3",
                                                  n_boot=100, seed=7)}}}
    summary = level_per_channel_max_pairwise(per_state)
    entry = summary["ASC"]["F"]
    check("the winning pair is named",
          entry["winning_pair"] == "healthy - DAM 3", str(entry["winning_pair"]))
    check("with the signed delta, the interval and both counts",
          entry["delta"] is not None and entry["ci"]["low"] is not None
          and entry["n_healthy"] == 5 and entry["n_damaged"] == 5)
    check("and the interval is scoped to that pair",
          "not of the maximum" in entry["ci"]["scope"])
    broken = {"ASC": {"F": {"DAM 3": dict(per_state["ASC"]["F"]["DAM 3"],
                                          n_damaged=None)}}}
    try:
        level_per_channel_max_pairwise(broken)
        check("an incomplete maximum raises", False, "no exception")
    except ValueError as exc:
        check("an incomplete maximum raises", True, str(exc)[:70])

    print("\n6. the wording guard")
    try:
        _guard_text("the 95% CI of max SDS is wide")
        check("the guard refuses an interval on a selection", False)
    except ValueError:
        check("the guard refuses an interval on a selection", True)
    check("a correct sentence passes",
          _guard_text("95% CI of the winning pair") != "")


    print("\n7. the registry matches the tables")
    for site, spec in SITES.items():
        for kind in ("committed", "support"):
            path = os.path.join(REPO, spec[kind])
            if not os.path.exists(path):
                check(f"{site}/{kind}: table present", False, path)
                continue
            rows = read_rows(path, spec)
            header = set(rows[0]) if rows else set()
            missing = [field for _dim, field in spec["dimensions"]
                       if field not in header]
            check(f"{site}/{kind}: every dimension column exists",
                  not missing, str(missing))
            present = {row.get(spec["channel_column"]) for row in rows}
            check(f"{site}/{kind}: both channels exist",
                  set(spec["channels"]) <= present,
                  str(sorted(value for value in present if value)))
            try:
                _states_in(rows, spec)
                check(f"{site}/{kind}: the states are the registry's", True)
            except ValueError as exc:
                check(f"{site}/{kind}: the states are the registry's", False,
                      str(exc)[:70])

    print("\n8. the LUMO reference pin (committed table)")
    spec = SITES["lumo"]
    rows = read_rows(os.path.join(REPO, spec["committed"]), spec)
    pooled = level_pooled_all_damage(rows, spec, n_boot=100, seed=7)
    pin = pin_emv_reference(spec["emv_reference"], pooled, REPO)
    if not pin.get("available"):
        check("the pin is available", False, pin.get("reason"))
    else:
        worst = max((abs(check["recomputed"] - check["reference"])
                     for check in pin["checks"]
                     if check["recomputed"] is not None
                     and check["reference"] is not None), default=None)
        check("the committed EMV reference is reproduced", pin["ok"],
              f"worst |diff| {worst}")
        for check in pin["checks"]:
            if not check["ok"]:
                print(f"        {check['dimension']}: "
                      f"{check['recomputed']} vs {check['reference']}")

    print()
    if failures:
        print(f"FAILED — {len(failures)} check(s): {failures}")
        return 1
    print(f"all checks passed ({len(SITES)} site(s))")
    return 0


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def build_parser():
    parser = argparse.ArgumentParser(
        description=("SDS of the OCV dimensions on the support-sized chips, "
                     "per observation channel and per state pair"))
    parser.add_argument("--site", default="all",
                        choices=("all",) + tuple(SITES),
                        help="which site(s) to analyse (default: all)")
    parser.add_argument("--out", default=None,
                        help="output directory (default: data/sds)")
    parser.add_argument("--json", default=None,
                        help="write the JSON here instead of --out")
    parser.add_argument("--md", default=None,
                        help="write the markdown here instead of --out")
    parser.add_argument("--boot", type=int, default=BOOT_MAIN,
                        help=f"bootstrap resamples per comparison "
                             f"(default {BOOT_MAIN})")
    parser.add_argument("--seed", type=int, default=SEED,
                        help=f"bootstrap seed (default {SEED})")
    parser.add_argument("--selftest", action="store_true",
                        help="run the offline checks and exit")
    parser.add_argument("--quiet", action="store_true")
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    if args.selftest:
        return selftest()

    wanted = list(SITES) if args.site == "all" else [args.site]
    print(f"SDS on the support-sized chips — {', '.join(wanted)}")
    report = analysis(wanted, REPO, args.boot, args.seed)
    if not args.quiet:
        for site, block in report["sites"].items():
            for kind in ("committed", "support"):
                table = block["tables"][kind]
                print(f"  {site:8s} {kind:9s} {table['n_rows']:4d} rows  "
                      f"{table['sha256'][:12]}")
        for site, block in report["sites"].items():
            pin = block.get("pin") or {}
            if pin.get("available"):
                print(f"  pin {site}: committed EMV "
                      f"{'reproduced' if pin.get('ok') else 'MISMATCH'}")
    out_dir = args.out
    if args.json or args.md:
        target = os.path.dirname(os.path.abspath(args.json or args.md))
        json_path, md_path = write_outputs(report, REPO, target, quiet=args.quiet)
        if args.json:
            os.replace(json_path, args.json)
        if args.md:
            os.replace(md_path, args.md)
    else:
        write_outputs(report, REPO, out_dir, quiet=args.quiet)

    failed = [site for site, block in report["sites"].items()
              if (block.get("pin") or {}).get("available")
              and not (block.get("pin") or {}).get("ok")]
    if failed:
        print(f"REFUSED: the committed reference is not reproduced for "
              f"{failed}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())


