#!/usr/bin/env python3
"""Bautzen: the damage ladder, as a **paired** bridge-minus-control contrast.

Bautzen does not fit the pre/post decomposition of ``ocv_decompose.py``, and the
record itself says why:

* the states are a **dose ladder** ``RS -> DS1 -> DS2``, not a step (RS 28,
  DS1 92, DS2 22);
* the states are **also the calendar**: RS exists only in April/May, DS2 only in
  September/October, so state and season are the same axis and a seasonal or
  temperature control cannot be identified -- there is no temperature column
  either;
* ``gamma2`` is the floor-corrected band coherence, which **can be negative**
  (DS1 median -0.021), so it cannot be logged;
* ``D`` is ``A / bbox_area`` exactly and ``A == n_band`` in every row, so neither
  can serve as evidence *against* the area;
* and the table carries a **control track 500 m east** on the same dates, which
  is the one design feature the other sites do not have.

That last point is the design: every ``(date, orbit)`` pair exists twice, once on
the bridge and once on the control, so the analysis works on the **difference**
bridge minus control. That difference removes the season, the acquisition day and
everything else the two tracks share -- which is exactly what the missing
temperature control would otherwise have to do -- and it is the only contrast
here that can be read as structural.

Design
------
* differences on the transformed scale: ``log`` for ``P, D, F`` (and ``log A`` as
  the covariate), the identity for ``gamma2`` (signed) and ``S``;
* contrasts: ``RS -> DS1``, ``RS -> DS2`` and ``RS -> {DS1, DS2}``, each with and
  without ``log A``;
* moving-block bootstrap inside each state, blocks never crossing the dose step;
* Holm over the targets within a contrast;
* the within-bridge and within-control contrasts are printed beside the
  difference, so the reader sees what the pairing removed.

What this design cannot do: there are **no pseudo events**. The state *is* a
campaign phase, so a placebo inside the record would be a placebo in time only.
Cells here are therefore descriptive/supporting by construction and can never
reach the ``generalisation-grade`` label of the four-site comparison.

Usage
-----
  python3 code/orbits/ocv_dose.py
  python3 code/orbits/ocv_dose.py --selftest
  python3 code/orbits/ocv_dose.py --verify

Outputs
-------
  data/bautzen/support/bautzen_ocv_dose.json
  data/bautzen/support/bautzen_ocv_dose.md
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
for path in (HERE, os.path.join(REPO, "code", "carola"),
             os.path.join(REPO, "code", "lumo")):
    if path not in sys.path:
        sys.path.insert(0, path)

import carola_ocv_stats as stats_mod  # noqa: E402
from carola_ocv_precollapse import holm, sha256_of  # noqa: E402
from ocv_decompose import block_index  # noqa: E402

SCHEMA = "ocv-dose-v1"
SITE = "bautzen"
LABEL = "Bautzen bridge (RS -> DS1 -> DS2, with a control track)"
TABLE = os.path.join(REPO, "data", SITE, "support", f"{SITE}_ocv_channels.csv")
OUT_JSON = os.path.join(REPO, "data", SITE, "support", f"{SITE}_ocv_dose.json")
OUT_MD = os.path.join(REPO, "data", SITE, "support", f"{SITE}_ocv_dose.md")

COLUMNS = {"gamma2": "gamma2_band_floor_corrected",
           "gamma2_raw": "gamma2_band_raw",
           "P": "P", "D": "D", "F": "F", "S": "S", "A": "A"}
TARGETS = ("gamma2", "P", "D", "F", "S")
LOG_TARGETS = ("P", "D", "F")
#: ``RS -> DS1``, ``RS -> DS2``, ``RS -> {DS1, DS2}``.
CONTRASTS = {
    "RS_vs_DS1": ("RS", ("DS1",)),
    "RS_vs_DS2": ("RS", ("DS2",)),
    "RS_vs_damaged": ("RS", ("DS1", "DS2")),
}
BLOCK = 5
BOOT = 2000
SEED = 7


# ---------------------------------------------------------------------------
# Data: pairs of (date, orbit), bridge and control
# ---------------------------------------------------------------------------
def _number(row, key):
    raw = (row.get(key) or "").strip()
    if raw in ("", "NA", "null", "NaN"):
        return None
    try:
        return float(raw)
    except ValueError:
        return None


def _group(series):
    return "bridge" if (series or "").endswith("Br\u00fccke") else "control"


def load(table=None):
    """Every row, keyed by ``(date, orbit)``; bridge and control separated."""
    table = table or TABLE
    rows = []
    with open(table, newline="") as handle:
        for raw in csv.DictReader(handle):
            values = {name: _number(raw, column) for name, column in COLUMNS.items()}
            rows.append({
                "date": dt.date.fromisoformat(raw["acquisition_date"]),
                "key": (raw["acquisition_date"], raw["orbit_direction"]),
                "orbit": raw["orbit_direction"],
                "state": raw["state"],
                "group": _group(raw["series"]),
                "series": raw["series"],
                "values": values,
            })
    rows.sort(key=lambda row: (row["date"], row["orbit"]))
    return rows


def pairs(rows):
    """``{key: {'bridge': row, 'control': row}}`` for the keys that have both."""
    out = {}
    for row in rows:
        bucket = out.setdefault(row["key"], {})
        bucket[row["group"]] = row
    missing = [key for key, bucket in out.items()
               if "bridge" not in bucket or "control" not in bucket]
    complete = {key: bucket for key, bucket in out.items()
                if "bridge" in bucket and "control" in bucket}
    return complete, missing


def transform(target, value):
    """``log`` only where it is defined and meant; ``gamma2`` stays signed."""
    if value is None:
        return None
    if target in LOG_TARGETS:
        return math.log(value) if value > 0 else None
    return value


def log_a(row):
    area = row["values"]["A"]
    return math.log(area) if (area is not None and area > 0) else None


def differences(complete, target):
    """Paired ``g(bridge) - g(control)`` per key, in date order."""
    out = []
    for key in sorted(complete):
        bucket = complete[key]
        bridge, control = bucket["bridge"], bucket["control"]
        left = transform(target, bridge["values"][target])
        right = transform(target, control["values"][target])
        area_left, area_right = log_a(bridge), log_a(control)
        if left is None or right is None or area_left is None or area_right is None:
            continue
        out.append({"key": key, "date": bridge["date"],
                    "orbit": bridge["orbit"], "state": bridge["state"],
                    "d": left - right, "d_area": area_left - area_right,
                    "bridge": left, "control": right})
    return out


# ---------------------------------------------------------------------------
# Fits: the paired contrast, then the same contrast inside each track
# ---------------------------------------------------------------------------
def designs(observations, state_a, states_b, with_area=False):
    used = [o for o in observations
            if o["state"] == state_a or o["state"] in states_b]
    if len(used) < 6:
        return None
    y = np.array([o["d"] for o in used], float)
    post = np.array([1.0 if o["state"] in states_b else 0.0 for o in used], float)
    columns = [np.ones(len(used)), post]
    names = ["intercept", "post"]
    if with_area:
        columns.append(np.array([o["d_area"] for o in used], float))
        names.append("d_logA")
    return {"rows": used, "y": y, "post": post,
            "X": np.column_stack(columns), "names": names}


def level_designs(observations, state_a, states_b, with_area=False):
    """The same contrast on one track's levels (bridge or control alone)."""
    used = [o for o in observations
            if o["state"] == state_a or o["state"] in states_b]
    if len(used) < 6:
        return None
    y = np.array([o["level"] for o in used], float)
    post = np.array([1.0 if o["state"] in states_b else 0.0 for o in used], float)
    columns = [np.ones(len(used)), post]
    names = ["intercept", "post"]
    if with_area:
        columns.append(np.array([o["log_area"] for o in used], float))
        names.append("logA")
    return {"rows": used, "y": y, "post": post,
            "X": np.column_stack(columns), "names": names}


def bootstrap_post(prepared, boot=BOOT, seed=SEED):
    """Moving-block bootstrap inside each state; blocks never cross the step."""
    post = prepared["post"]
    groups = [np.where(post == 0)[0], np.where(post == 1)[0]]
    lengths = [max(2, min(BLOCK, group.size // 3)) for group in groups]
    column = prepared["names"].index("post")
    observed = stats_mod.ols_fit(prepared["y"], prepared["X"], prepared["names"])
    rng = np.random.default_rng(seed)
    draws = np.empty(boot)
    for draw in range(boot):
        index = np.concatenate([group[block_index(group.size, length, rng)]
                                for group, length in zip(groups, lengths)])
        beta, *_ = np.linalg.lstsq(prepared["X"][index], prepared["y"][index],
                                   rcond=None)
        draws[draw] = float(beta[column])
    low, high = np.percentile(draws, [2.5, 97.5])
    share_le = float((draws <= 0).mean())
    share_ge = float((draws >= 0).mean())
    return {"beta": float(observed["beta"][column]),
            "ci": [float(low), float(high)],
            "p_boot": min(1.0, 2.0 * min(share_le, share_ge)),
            "se_boot": float(draws.std(ddof=1)),
            "n": int(observed["n"]),
            "n_a": int((post == 0).sum()), "n_b": int((post == 1).sum()),
            "block": lengths, "boot": int(boot), "seed": int(seed),
            "excludes_null": bool(low > 0 or high < 0)}


# ---------------------------------------------------------------------------
# The report
# ---------------------------------------------------------------------------
def level_observations(complete, target, group):
    out = []
    for key in sorted(complete):
        row = complete[key].get(group)
        if row is None:
            continue
        value = transform(target, row["values"][target])
        area = log_a(row)
        if value is None or area is None:
            continue
        out.append({"key": key, "state": row["state"], "level": value,
                    "log_area": area})
    return out


def analyse(table=None, boot=BOOT):
    rows = load(table)
    complete, missing = pairs(rows)
    states = {}
    for key, bucket in complete.items():
        state = bucket["bridge"]["state"]
        states[state] = states.get(state, 0) + 1
    report = {
        "schema": SCHEMA, "site": SITE, "label": LABEL,
        "table": {"path": os.path.relpath(table or TABLE, REPO),
                  "sha256": sha256_of(table or TABLE), "n_rows": len(rows)},
        "pairs": {"complete": len(complete), "incomplete": len(missing),
                  "by_state": states,
                  "dates": [min(key[0] for key in complete),
                            max(key[0] for key in complete)]},
        "design": "paired bridge-minus-control per (date, orbit), then the state "
                  "contrast on that difference",
        "targets": list(TARGETS), "log_targets": list(LOG_TARGETS),
        "contrasts": {}, "area": {},
        "rules": {"block": BLOCK, "boot": int(boot), "seed": SEED,
                  "placebo": "not applicable: the state is a campaign phase, so a "
                             "pseudo event would be a placebo in time only",
                  "gamma2": "sigma_band_floor_corrected is signed, identity scale",
                  "tautology": "D = A / bbox_area exactly; A == n_band in every row"},
    }
    for name, (state_a, states_b) in CONTRASTS.items():
        block = {"state_a": state_a, "states_b": list(states_b), "targets": {}}
        pvals = []
        for target in TARGETS:
            diffs = differences(complete, target)
            plain = designs(diffs, state_a, states_b, with_area=False)
            adjusted = designs(diffs, state_a, states_b, with_area=True)
            bridge = level_designs(level_observations(complete, target, "bridge"),
                                   state_a, states_b)
            control = level_designs(level_observations(complete, target, "control"),
                                    state_a, states_b)
            cell = {
                "did": bootstrap_post(plain, boot=boot) if plain else None,
                "did_area": bootstrap_post(adjusted, boot=boot) if adjusted else None,
                "bridge": bootstrap_post(bridge, boot=boot) if bridge else None,
                "control": bootstrap_post(control, boot=boot) if control else None,
            }
            block["targets"][target] = cell
            pvals.append((cell["did_area"] or {}).get("p_boot"))
        adjusted_p, family = holm(pvals)
        for target, p in zip(TARGETS, adjusted_p):
            if block["targets"][target]["did_area"] is not None:
                block["targets"][target]["did_area"]["holm_p"] = p
        block["holm_family"] = family
        # the area itself, as the thing that moves
        diffs = differences(complete, "A")
        if diffs:
            used = [o for o in diffs if o["state"] == state_a or o["state"] in states_b]
            if len(used) >= 6:
                post = np.array([1.0 if o["state"] in states_b else 0.0 for o in used])
                X = np.column_stack([np.ones(len(used)), post])
                beta, *_ = np.linalg.lstsq(X, np.array([o["d_area"] for o in used]),
                                           rcond=None)
                pre = [o["d_area"] for o in used if o["state"] == state_a]
                postv = [o["d_area"] for o in used if o["state"] in states_b]
                report["area"][name] = {
                    "did_log_area": float(beta[1]),
                    "ratio": float(math.exp(beta[1])),
                    "n_a": len(pre), "n_b": len(postv),
                    "median_d_a": float(np.median(pre)) if pre else None,
                    "median_d_b": float(np.median(postv)) if postv else None,
                }
        report["contrasts"][name] = block
    return report


# ---------------------------------------------------------------------------
# Markdown
# ---------------------------------------------------------------------------
def _f(value, digits=3, signed=True):
    if value is None:
        return "—"
    return f"{value:+.{digits}f}" if signed else f"{value:.{digits}f}"


def _ci(pair, digits=3):
    if not pair:
        return "—"
    return f"[{pair[0]:+.{digits}f}, {pair[1]:+.{digits}f}]"


def _p(value):
    if value is None:
        return "—"
    return "<0.001" if value < 0.001 else f"{value:.3f}"


def to_markdown(report):
    lines = []
    add = lines.append
    add(f"# {report['label']} — the damage ladder as a paired contrast")
    add("")
    add(f"* table `{report['table']['path']}` (sha256 "
        f"`{report['table']['sha256'][:16]}…`), {report['table']['n_rows']} rows")
    add(f"* {report['pairs']['complete']} complete `(date, orbit)` pairs, "
        f"{report['pairs']['incomplete']} without a partner; "
        f"{report['pairs']['dates'][0]} … {report['pairs']['dates'][1]}")
    add(f"* pairs per state: "
        + ", ".join(f"{k} {v}" for k, v in sorted(report['pairs']['by_state'].items())))
    add("")
    add("## Why the difference, and not the state itself")
    add("")
    add("The states are also the calendar: `RS` exists only in April/May, `DS2` only "
        "in September/October. A state effect and a season effect are therefore the "
        "same axis in this record, and the table carries no temperature. What it "
        "does carry is a **control track 500 m east on the same dates**, so every "
        "`(date, orbit)` exists twice and the analysis works on the difference "
        "bridge minus control. That removes whatever the two tracks share — season, "
        "acquisition day, processing — and leaves the structural part.")
    add("")
    add("`gamma2` is the floor-corrected band coherence and can be **negative**, so it "
        "is analysed on the identity scale; `P, D, F` are logged. `D` is "
        "`A / bbox_area` exactly and `A == n_band` in every row, so neither can be "
        "read as evidence against the area.")
    add("")
    add("## The area itself")
    add("")
    add("| contrast | did log A | ratio | n (RS / damaged) |")
    add("|---|---|---|---|")
    for name, area in report["area"].items():
        add(f"| {name} | {_f(area['did_log_area'])} | {area['ratio']:.2f} | "
            f"{area['n_a']} / {area['n_b']} |")
    add("")
    add("## The contrasts")
    add("")
    for name, block in report["contrasts"].items():
        add(f"### {name} — {block['state_a']} vs "
            + "/".join(block["states_b"]))
        add("")
        add("| target | Did (bridge−control) | + dlog A | within bridge | "
            "within control | Holm |")
        add("|---|---|---|---|---|---|")
        for target in report["targets"]:
            cell = block["targets"][target]
            did = cell.get("did") or {}
            did_area = cell.get("did_area") or {}
            bridge = cell.get("bridge") or {}
            control = cell.get("control") or {}
            add(f"| `{target}` | {_f(did.get('beta'))} {_ci(did.get('ci'))} | "
                f"{_f(did_area.get('beta'))} {_ci(did_area.get('ci'))} | "
                f"{_f(bridge.get('beta'))} {_ci(bridge.get('ci'), 2)} | "
                f"{_f(control.get('beta'))} {_ci(control.get('ci'), 2)} | "
                f"{_p(did_area.get('holm_p'))} |")
        add("")
        add(f"Holm family: {block['holm_family']} targets within this contrast.")
        add("")
    add("## What this can and cannot say")
    add("")
    add("* **No placebo.** The state is a campaign phase, so a pseudo event would be "
        "a placebo in time only. These cells are descriptive/supporting by "
        "construction and never reach the `generalisation-grade` label of the "
        "four-site comparison.")
    add("* **No temperature, no season control.** Both are absorbed by the pairing, "
        "not estimated.")
    add("* **`D` is the area.** `D = A / bbox_area` in every row, so a `D` effect "
        "after `log A` adjustment is a change of the bounding box, not an "
        "independent dimension.")
    add("* **`gamma2` is signed** (floor-corrected) and unlike the other sites not "
        "log-transformable.")
    add("")
    return "\n".join(lines)


def write_outputs(report, json_path=None, md_path=None, quiet=False):
    json_path = json_path or OUT_JSON
    md_path = md_path or OUT_MD
    for path in (json_path, md_path):
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
    return json_path, md_path


# ---------------------------------------------------------------------------
# Self-test: a planted structural step, and a step both tracks share
# ---------------------------------------------------------------------------
def _synthetic(kind, n_rs=14, n_ds2=12, seed=3):
    """Rows in ``load()``'s shape; ``kind`` decides where the step lives."""
    rng = np.random.default_rng(seed)
    rows = []
    start = dt.date(2025, 4, 2)
    index = 0
    for state, count in (("RS", n_rs), ("DS2", n_ds2)):
        for _ in range(count):
            day = start + dt.timedelta(days=index)
            index += 1
            shared = 0.0
            bridge_only = 0.4 if state == "DS2" else 0.0
            if kind == "both":
                shared, bridge_only = bridge_only, 0.0
            base = 10.0 + shared + rng.normal(0, 0.05)
            for group in ("bridge", "control"):
                shift = bridge_only if group == "bridge" else 0.0
                rows.append({
                    "date": day, "key": (day.isoformat(), "ASC"),
                    "orbit": "ASC", "state": state, "group": group,
                    "series": "ASC Br\u00fccke" if group == "bridge" else "ASC Kontrolle",
                    "values": {"P": math.exp(base + shift), "D": math.exp(base + shift),
                               "F": math.exp(base + shift), "S": base + shift,
                               "gamma2": -0.01 + 0.01 * (base - 10.0) + (0.3 * (state == "DS2") if group == "bridge" else 0.0),
                               "A": math.exp(2.0 + 0.1 * (state == "DS2") * (group == "bridge")),
                               "gamma2_raw": math.exp(base + shift)},
                })
    return rows


def selftest():
    failures = []

    def check(label, condition, detail=""):
        print(f"  [{'ok  ' if condition else 'FAIL'}] {label}"
              + (f" — {detail}" if detail else ""))
        if not condition:
            failures.append(label)

    print("\n1. Pairing and the signed gamma2")
    rows = _synthetic("bridge")
    complete, missing = pairs(rows)
    check("every (date, orbit) pairs once", len(complete) == 26 and not missing,
          f"{len(complete)} complete, {len(missing)} missing")
    check("gamma2 keeps a negative value on the identity scale",
          transform("gamma2", -0.02) == -0.02)
    check("a non-positive log target is refused, not logged",
          transform("P", -1.0) is None and transform("P", 0.0) is None)

    print("\n2. A step that exists only at the bridge")
    complete, _ = pairs(_synthetic("bridge"))
    diffs = differences(complete, "P")
    prepared = designs(diffs, "RS", ("DS2",))
    fit = bootstrap_post(prepared, boot=300)
    check("the paired contrast finds it", fit["beta"] > 0.25 and fit["excludes_null"],
          f"beta {fit['beta']:+.3f}, CI {fit['ci']}")
    control = level_designs(level_observations(complete, "P", "control"),
                            "RS", ("DS2",))
    control_fit = bootstrap_post(control, boot=300)
    check("the control track alone does not", not control_fit["excludes_null"],
          f"beta {control_fit['beta']:+.3f}, CI {control_fit['ci']}")

    print("\n3. A step both tracks share")
    complete, _ = pairs(_synthetic("both"))
    fit = bootstrap_post(designs(differences(complete, "P"), "RS", ("DS2",)),
                         boot=300)
    check("the pairing removes it", not fit["excludes_null"],
          f"beta {fit['beta']:+.3f}, CI {fit['ci']}")

    print("\n4. An incomplete pair is reported, never filled")
    broken = _synthetic("bridge")
    broken = [row for row in broken
              if not (row["group"] == "control" and row["key"] == broken[0]["key"])]
    complete, missing = pairs(broken)
    check("the orphan key is named", len(missing) == 1 and len(complete) == 25,
          f"{len(complete)} complete, {len(missing)} missing")

    print("\n5. The report")
    report = analyse(boot=100)
    text = to_markdown(report)
    check("the report states that no placebo is possible",
          "No placebo" in text)
    check("and that D is the area", "D = A / bbox_area" in text)
    check("the contrasts are all present",
          set(report["contrasts"]) == set(CONTRASTS), str(sorted(report["contrasts"])))

    print()
    if failures:
        print(f"FAILED — {len(failures)} check(s): {failures}")
        return 1
    print("all checks passed")
    return 0


def build_parser():
    parser = argparse.ArgumentParser(
        description=("Bautzen: the damage ladder as a paired bridge-minus-control "
                     "contrast"))
    parser.add_argument("--table", default=None)
    parser.add_argument("--json", default=None)
    parser.add_argument("--md", default=None)
    parser.add_argument("--boot", type=int, default=BOOT)
    parser.add_argument("--verify", action="store_true")
    parser.add_argument("--selftest", action="store_true")
    parser.add_argument("--quiet", action="store_true")
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    if args.selftest:
        return selftest()
    if args.verify:
        return verify(table=args.table)
    print(f"{LABEL}")
    report = analyse(args.table, boot=args.boot)
    write_outputs(report, args.json, args.md, quiet=args.quiet)
    if not args.quiet:
        print(f"  {report['pairs']['complete']} pairs, "
              f"{report['pairs']['incomplete']} incomplete")
        for name, block in report["contrasts"].items():
            area = report["area"].get(name) or {}
            print(f"  {name}: did log A {area.get('did_log_area', float('nan')):+.3f}"
                  f" (ratio {area.get('ratio', float('nan')):.2f})")
            for target in report["targets"]:
                cell = block["targets"][target].get("did_area") or {}
                print(f"    {target:7s} Did+A {cell.get('beta', float('nan')):+7.3f} "
                      f"[{cell.get('ci', [float('nan'), float('nan')])[0]:+.3f},"
                      f"{cell.get('ci', [float('nan'), float('nan')])[1]:+.3f}]  "
                      f"holm {_p(cell.get('holm_p'))}")
    return 0



def verify(report=None, table=None):
    """Recompute the deterministic point estimates from the table."""
    if report is None:
        with open(OUT_JSON) as handle:
            report = json.load(handle)
    table = table or TABLE
    if report["table"]["sha256"] != sha256_of(table):
        print("FAILED — the table on disk is not the one this report came from")
        return 1
    rows = load(table)
    complete, missing = pairs(rows)
    failures, checked = [], 0
    if len(complete) != report["pairs"]["complete"]:
        failures.append(f"pairs: {len(complete)} vs {report['pairs']['complete']}")
    if len(missing) != report["pairs"]["incomplete"]:
        failures.append(f"incomplete pairs: {len(missing)} vs "
                        f"{report['pairs']['incomplete']}")
    for name, (state_a, states_b) in CONTRASTS.items():
        block = report["contrasts"][name]
        for target in TARGETS:
            diffs = differences(complete, target)
            for key, prepared in (("did", designs(diffs, state_a, states_b)),
                                  ("did_area", designs(diffs, state_a, states_b, True))):
                cell = block["targets"][target].get(key) or {}
                if prepared is None or not cell:
                    continue
                fit = stats_mod.ols_fit(prepared["y"], prepared["X"],
                                        prepared["names"])
                beta = float(fit["beta"][prepared["names"].index("post")])
                checked += 1
                if abs(beta - cell["beta"]) > 1e-9:
                    failures.append(f"{name}/{target}/{key}: {beta} != {cell['beta']}")
    print(f"  {checked} field(s) recomputed from {os.path.relpath(table, REPO)}")
    if failures:
        print(f"FAILED — {len(failures)} mismatch(es): {failures[:5]}")
        return 1
    print("  all recomputed values match the JSON")
    return 0

if __name__ == "__main__":
    sys.exit(main())
