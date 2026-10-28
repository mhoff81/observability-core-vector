#!/usr/bin/env python3
"""The frozen-reference **prospective** run of one case, read by two readers.

The question
------------
The descriptive plates ask *where in the record the vector moves*. This one asks
what a reader could have said **before** the event, with the reference frozen at a
window a stated rule selected, and with nothing read after the event.

It is one run with **two readers over the same numbers**:

* the **trajectory rule** of this repository -- ``ocv_trajectory.py``'s own frozen
  quantiles, its monthly envelope of the reference months, its persistence
  minimum -- and
* **EAS policy v1.0.0** (``erlking_mvp/eas/policy_v1.json``, digest-pinned in
  ``data/cases/eas_policy_pin.json``), the rule set the live engine applies.

They are reported as a **co-primary pair** and never collapsed into one verdict.
Where they agree the agreement is the finding; where they differ the difference
is a property of the two measurement layers and is reported as one. Neither
reader has precedence, and no sentence here says one of them is right.

Why the window comes from a rule
--------------------------------
Morandi's registry declares **no** reference window: its healthy state is the
whole 198-chip pre-event record, which leaves no held-out class at all. The
window therefore has to come from a rule written down before the read, and that
rule may not look at the event or at any measured quantity:

    the earliest run of twelve consecutive calendar months in which every month
    holds at least four chips

It reads acquisition dates and per-month counts, nothing else. ``resolve_window``
takes no event argument and never calls ``event_day``; ``--selftest`` asserts
that from the function's own source. A healthy-state window that does not end at
least a year before the event is refused outright, so the rule cannot swallow
the months it is supposed to be tested against.

The control keeps its registry window
-------------------------------------
Carola already declares one, and its descriptive plate is published against it.
Reading the control through the same *layer* -- same standardisation, same
envelope construction, same persistence rule -- is what makes the comparison
fair; re-selecting the control's window would move a published number for the
sake of tidiness. The difference in window *source* is printed beside the
comparison rather than hidden by it.

What the run will not do
------------------------
* It re-derives no metric. Every number is one the plate already computes:
  per-chip ``R_B = rms(gamma2, P)``, ``R_G = rms(A, D, F, S)``, the pooled-orbit
  standardisation, the monthly median, the monthly envelope of the frozen
  reference months. The monthly-RMS variant is reported **beside** the per-chip
  one and never replaces it.
* It does not pool cases, and the control is never an input to the case.
* It does not read the post-event tail for anything but the hindsight guards.
* It does not turn a null into a failure: a trajectory positive beside a live
  null is exactly the co-primary finding this run exists to print.

The live reader uses the policy's own words
-------------------------------------------
``NOMINAL``, ``WATCH``, ``INSUFFICIENT`` -- read from the pinned policy, unaltered.
An earlier draft used a private word, ``NOT_EVALUABLE``, for the rule this cut
cannot read; it was dropped because v1.0.0 already names that condition
``INSUFFICIENT`` ("a missing observable the rule needed"). ``DEGRADED_OBSERVABILITY``,
``UNRESOLVED`` and ``CORROBORATED_WATCH`` are unreachable in an observability cut
that carries no observability level and no second stream, and the artefact says
so beside the word instead of letting a reader assume the whole policy spoke.

Outputs
-------
``data/morandi/fig_morandi_ocv_trajectory_prospective.json``
``figures/morandi/fig_morandi_ocv_trajectory_prospective.{md,jpg,sh}``
``data/cases/fig_ocv_prospective_comparison.json`` and its report

The descriptive plate of the same case is **not** touched: the registry keeps its
window, and ``fig_morandi_ocv_trajectory.{json,md,jpg}`` stay byte-identical
(``ocv_trajectory.py --verify`` checks that independently).

Usage:
  python3 code/cases/ocv_prospective.py
  python3 code/cases/ocv_prospective.py --selftest
  python3 code/cases/ocv_prospective.py --verify
"""
from __future__ import annotations

import argparse
import datetime as dt
import inspect
import json
import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, os.pardir, os.pardir))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

import ocv_trajectory as plate  # noqa: E402

SCHEMA = "ocv-prospective-v1"
CASE = "morandi"
CONTROL = "carola"
DATA = os.path.join(REPO, "data")
FIGURES = os.path.join(REPO, "figures")
PIN = os.path.join(DATA, "cases", "eas_policy_pin.json")

#: The rule that selects the reference window: written down before the run, read
#: by ``resolve_window`` alone, and never by a threshold or a state.
WINDOW_RULE = {
    "name": "earliest full cycle",
    "months": 12,
    "min_chips": 4,
    "statement": ("the earliest run of twelve consecutive calendar months in "
                  "which every month holds at least four chips"),
    "reads": ["the acquisition date of every chip",
              "the number of chips in each calendar month"],
    "never_reads": ["the event date", "any measured quantity of any chip"],
}

#: How far the frozen window must stay from the event: a full year of calendar
#: months. The unit is months, not days, because every bin here is a month -- a
#: 365-day rule would refuse a window that is twelve bins away to the day.
GUARD_MONTHS = 12

#: The reader vocabularies. The trajectory words are this run's own, and are
#: carried beside the plate's five-word verdict; the live words are the pinned
#: policy's, verbatim, and nothing else may be printed in that field.
TRAJECTORY_STATES = ("STABLE", "WATCH", "CHANGE", "NOT_ASSESSABLE")
LIVE_STATES = ("NOMINAL", "WATCH", "INSUFFICIENT")

#: How a count is spelled when a sentence opens with it.
WORD_NUMBERS = {1: "One", 2: "Two", 3: "Three", 4: "Four", 5: "Five",
                6: "Six", 7: "Seven", 8: "Eight", 9: "Nine"}


def word_number(n):
    """``5 -> 'five'``, ``12 -> '12'``: a count as a word where one exists.

    Counts above nine keep their digits rather than growing a second spelling
    table, so a page never prints "twelve" where the data says twelve for a
    reason nobody can check.
    """
    word = WORD_NUMBERS.get(n)
    return word.lower() if word else "%d" % n


def sentence_number(n):
    """``word_number`` for a count a sentence opens with: ``5 -> 'Five'``."""
    text = word_number(n)
    return text[:1].upper() + text[1:]

MAGNITUDES = ("R_B", "R_G")
STREAMS = MAGNITUDES + plate.DIMS

#: The two readers, named. Neither is derived from the other.
READERS = ("trajectory", "live")

#: Every guard the run promises, in the order the report prints them. The
#: selftest asserts that this list and the guards actually computed agree in both
#: directions, so a guard cannot be dropped from the page by being forgotten here.
GUARD_ORDER = (
    "window_event_free",
    "buckets_disjoint_and_reconciled",
    "no_pooling",
    "envelope_from_reference_months_only",
    "every_held_out_bin_fielded",
    "state_recomputes_from_stored_numbers",
    "record_is_max_over_streams",
    "control_after_case_no_influence",
    "runs_reported_separately_from_magnitude",
    "magnitudes_at_1e_12",
    "vocabulary_scan",
    "no_failure_language",
    "hindsight",
    "thresholds_pinned",
    "event_is_the_registrys",
    "policy_sha_matches_pin",
)


def month_labels(first, count):
    """``count`` consecutive ``YYYY-MM`` labels from ``first``."""
    year, month = int(first[:4]), int(first[5:7])
    labels = []
    for _ in range(count):
        labels.append("%04d-%02d" % (year, month))
        year, month = (year + 1, 1) if month == 12 else (year, month + 1)
    return labels


def after(label):
    """The calendar month after ``label``."""
    return month_labels(label, 2)[1]


def month_counts(chips):
    """``{YYYY-MM: n chips}`` -- a census of acquisition months, nothing else."""
    counts = {}
    for chip in chips:
        key = plate.month_text(plate.month_key(chip["day"]))
        counts[key] = counts.get(key, 0) + 1
    return counts


def month_index(label):
    """A calendar month as a whole number of months since year 0."""
    return int(label[:4]) * 12 + int(label[5:7]) - 1


def resolve_window(site, chips=None, event=None):
    """The frozen window ``WINDOW_RULE`` selects -- event-free by construction.

    The rule walks the occupied months from the first one and stops at the first
    run of ``months`` consecutive months in which every month holds at least
    ``min_chips`` chips. Every start it rejects is recorded with the months that
    failed, so a reader can re-run the rule by hand.

    ``event`` is **not** an input to the selection. It arrives only to be used
    against the answer: the chosen window must end at least ``GUARD_MONTHS``
    calendar months before the event's own month, or the function refuses. That is
    the one direction in which the event may not travel, and the selftest asserts
    that the rule itself never reads it.
    """
    if chips is None:
        chips, _cut = plate.load_chips(site)
    counts = month_counts(chips)
    occupied = sorted(counts)
    need, floor = WINDOW_RULE["months"], WINDOW_RULE["min_chips"]
    start = occupied[0]
    tried = []
    while True:
        span = month_labels(start, need)
        if span[-1] > occupied[-1]:
            raise SystemExit(
                "%s: the rule finds no run of %d months with >= %d chips"
                % (site, need, floor))
        short = [m for m in span if counts.get(m, 0) < floor]
        tried.append({"start": span[0], "end": span[-1], "short": short,
                      "shorts": [{"month": m, "chips": counts.get(m, 0)}
                                 for m in short]})
        if not short:
            break
        start = after(start)
    window = [span[0], span[-1]]
    resolved = {
        "site": site,
        "rule": dict(WINDOW_RULE),
        "window": window,
        "months": [{"month": m, "chips": counts.get(m, 0)} for m in span],
        "n_months": need,
        "n_chips": sum(counts.get(m, 0) for m in span),
        "first_occupied": occupied[0],
        "last_occupied": occupied[-1],
        "starts_tried": tried,
        "n_starts_tried": len(tried),
        "basis": "acquisition months and per-month chip counts only",
        "event_read_during_selection": False,
    }
    if event:
        event_month = plate.month_text(plate.month_key(dt.date.fromisoformat(event)))
        gap = month_index(event_month) - month_index(window[1])
        last = dt.date(int(window[1][:4]), int(window[1][5:7]), 28)
        resolved["event_month"] = event_month
        resolved["month_gap"] = gap
        resolved["days_before_event"] = (dt.date.fromisoformat(event) - last).days
        resolved["guard_months"] = GUARD_MONTHS
        resolved["guard_passed"] = bool(gap >= GUARD_MONTHS)
        if not resolved["guard_passed"]:
            raise SystemExit(
                "%s: the rule's window ends %d calendar months before the event "
                "month, inside the %d-month guard; refusing to read a held-out "
                "class that reaches into the year the window is tested against"
                % (site, gap, GUARD_MONTHS))
    return resolved


def build_case(site, window=None, rule=None):
    """Build the plate report for ``site`` against a window, then restore.

    ``window=None`` builds the site exactly as its registry declares it -- which
    is how the control is read, so this run cannot move the control.
    """
    if window is None:
        return plate.build(site)
    with plate.window_override(site, tuple(window),
                               rule or WINDOW_RULE["statement"]):
        return plate.build(site)


def partition(report):
    """The frozen buckets of one report, plus the reconciliation that ties them.

    The four audit buckets are the reference state's chosen window, the held-out
    months, the tail after the event and the *prefix* -- reference-state chips the
    frozen window leaves behind it. They are disjoint by construction
    (``chip_class`` is a single-valued function), and they must add up to the cut's
    own kept count, or the partition is wrong and nothing downstream is read.
    """
    classes = report["classes"]
    counts = classes["reconciliation"]
    labels = classes["labels"]
    buckets = {}
    for name in ("reference", "outside", "prefix", "post"):
        cell = classes.get(name) or {}
        buckets[name] = {"label": labels.get(name),
                         "n_chips": cell.get("n_chips", 0),
                         "months": cell.get("months", 0),
                         "window": cell.get("window")}
    buckets["unassigned"] = {"label": "unassigned",
                             "n_chips": counts.get("unassigned", 0),
                             "months": 0, "window": None,
                             "chips": classes.get("unassigned", [])}
    total = sum(buckets[name]["n_chips"] for name in buckets)
    return {"buckets": buckets,
            "reconciliation": dict(counts),
            "total": total,
            "kept": counts.get("kept"),
            "reconciles": bool(total == counts.get("kept")),
            "disjoint": True,
            "rule": ("one chip, one bucket: chip_class returns exactly one of "
                     "reference / outside / prefix / post / unassigned"),
            "empty": [name for name in buckets if not buckets[name]["n_chips"]]}


def phase_audit(report):
    """One row per monthly bin: its class, its census, and its position.

    A row is *held-out* when every chip in its bin is held out -- the plate's own
    definition, not a second one. Each stream of ``plate.DIMS`` is then placed
    against the frozen envelope of the reference months: the block magnitudes
    one-sided (a departure is a rise), the dimensions two-sided (``[q05, q95]``
    is the reference's own range). No threshold is chosen here; every one is read
    out of ``reference_envelope.monthly``.
    """
    monthly = report["reference_envelope"]["monthly"]
    rows = []
    for entry in report["months"]:
        census = {"reference": entry.get("n_reference", 0),
                  "outside": entry.get("n_outside", 0),
                  "prefix": entry.get("n_prefix", 0),
                  "post": entry.get("n_post", 0)}
        held_out = bool(entry["n"] and census["outside"] == entry["n"])
        values, flags = {}, {}
        for key in STREAMS:
            cell = monthly[key]
            value = plate.bin_value(entry, key)
            low, high = cell["q05"], cell["q95"]
            values[key] = value
            if value is None:
                position, fired = None, False
            elif value > high:
                position, fired = "above", True
            elif value < low:
                position, fired = "below", True
            else:
                position, fired = "inside", False
            flags[key] = {"value": value, "q05": low, "q95": high,
                          "position": position, "fired": fired,
                          "two_sided": key in plate.DIMS}
        klass = plate.month_class(entry)
        rows.append({
            "label": entry["label"],
            "month": entry["month"],
            "part": entry.get("part"),
            "class": klass,
            "n": entry["n"],
            "census": census,
            "held_out": held_out,
            "fielded": all(values[key] is not None for key in STREAMS),
            "reason": ("every chip in the bin is held out"
                       if held_out else
                       "before the frozen window" if klass == "prefix" else
                       "on or after the event" if klass == "post" else
                       "inside the frozen reference window"
                       if klass == "ref" else
                       "more than one class in the bin"),
            "values": values,
            "flags": flags,
        })
    held = [row for row in rows if row["held_out"]]
    return {"rows": rows,
            "held_out": held,
            "n_held_out": len(held),
            "unfielded": [row["label"] for row in held if not row["fielded"]],
            "fully_fielded": all(row["fielded"] for row in held),
            "envelope_source": ("the monthly medians of the frozen reference "
                                "months only"),
            "thresholds": {key: {"q05": monthly[key]["q05"],
                                 "q50": monthly[key]["q50"],
                                 "q95": monthly[key]["q95"]}
                           for key in STREAMS}}


def run_of(rows, flags):
    """The longest run of consecutive ``True`` in ``flags`` (rows aligned).

    Consecutiveness is the plate's own ``follows``, so a split bin (a month cut
    in two by the event) counts as one step, not as a gap.
    """
    best, current = [], []
    for row, flag in zip(rows, flags):
        if not flag:
            current = []
            continue
        if current and plate.follows(current[-1], row):
            current = current + [row]
        else:
            current = [row]
        if len(current) > len(best):
            best = list(current)
    return {"length": len(best),
            "months": [row["label"] for row in best],
            "start": best[0]["label"] if best else None,
            "end": best[-1]["label"] if best else None}


def streams(audit):
    """Per-stream runs over the held-out months, from the flags and nothing else."""
    held = audit["held_out"]
    result = {}
    for key in STREAMS:
        flags = [bool(row["flags"][key]["fired"]) for row in held]
        run = run_of(held, flags)
        values = [row["values"][key] for row in held
                  if row["values"][key] is not None]
        result[key] = {
            "n_months": len(held),
            "n_fired": sum(flags),
            "months_fired": [row["label"] for row, flag in zip(held, flags)
                             if flag],
            "longest_run": run,
            "run_length": run["length"],
            "max_abs": max((abs(value) for value in values), default=None),
            "min": min(values) if values else None,
            "max": max(values) if values else None,
            "median": plate.median(values) if values else None,
            "q05": audit["thresholds"][key]["q05"],
            "q95": audit["thresholds"][key]["q95"],
            "two_sided": key in plate.DIMS,
            "direction": ("both" if key in plate.DIMS else "above"),
        }
    return {"streams": result,
            "n_months": len(held),
            "order": list(STREAMS),
            "runs": {key: result[key]["run_length"] for key in STREAMS}}


def state_of_trajectory(record, persist=None):
    """The coarse trajectory state, from run lengths and one boolean.

    ``record`` is numeric only -- ``{"fielded": bool, "runs": {stream: months}}``
    -- so the function cannot know which asset, which month or which label it is
    reading. Rename every series and the answer is the same; that is the check the
    selftest makes.

    The rule, stated: a departure is *sustained* when a stream stays outside the
    frozen envelope for ``persist`` consecutive held-out months. ``CHANGE`` needs a
    sustained block magnitude **and** a sustained dimension, because a departure
    that moves the blocks alone is one plane moving, not the vector;
    ``WATCH`` is either one alone; ``STABLE`` is neither;
    ``NOT_ASSESSABLE`` is a held-out bin that does not carry every stream.
    """
    persist = plate.PERSIST_MIN if persist is None else persist
    if not record.get("fielded"):
        return "NOT_ASSESSABLE"
    runs = record["runs"]
    block = any(runs.get(key, 0) >= persist for key in MAGNITUDES)
    dimension = any(runs.get(key, 0) >= persist for key in plate.DIMS)
    if block and dimension:
        return "CHANGE"
    if block or dimension:
        return "WATCH"
    return "STABLE"


def state_of_live(record):
    """The pinned policy's state, from numbers and booleans only.

    ``record`` carries the longest runs over the policy's own floors, the six gate
    values and which of the policy's rules this cut could read at all. It carries
    no label, no date and no site.

    Precedence, and why: a failed gate is definitive whatever else fired, and so
    is a gate the cut cannot fill -- v1.0.0 calls both ``INSUFFICIENT`` ("too few
    months ... or a missing observable the rule needed"), so both land there and
    ``reason`` says which. Otherwise a firing rule is ``WATCH``, and silence is
    ``NOMINAL`` *over the rules this cut carries* -- never over the rules it cannot
    read, which the artefact states beside the word.
    """
    failed = [gate["gate"] for gate in record["gates"]
              if gate["passed"] is False]
    if failed:
        return "INSUFFICIENT"
    if record["unreadable"]:
        return "INSUFFICIENT"
    if any(rule["fired"] for rule in record["rules"] if rule["applied"]):
        return "WATCH"
    return "NOMINAL"


def live_reason(record, state):
    """The sentence that has to travel with a live state."""
    if state == "INSUFFICIENT":
        failed = [gate["gate"] for gate in record["gates"]
                  if gate["passed"] is False]
        if failed:
            return "quality gate(s) not met: %s" % ", ".join(failed)
        return "the cut cannot fill: %s" % ", ".join(record["unreadable"])
    if state == "WATCH":
        return "rule(s) firing: %s" % ", ".join(
            "%s/%s" % (rule["rule"], rule["quantity"])
            for rule in record["rules"] if rule["fired"] and rule["applied"])
    total = len(record.get("firing_rules") or ())
    return ("no rule this cut can read fires; %s of the policy's %s firing rules "
            "%s not carried by an observability cut"
            % (word_number(len(record["not_applied"])),
               word_number(total) if total else "own",
               "is" if len(record["not_applied"]) == 1 else "are"))


def live_record(report, audit, policy):
    """Apply the pinned policy's observability rules to this cut's numbers.

    Each rule is recorded as ``{rule, value, threshold, fired}`` whatever it did,
    so a reader can see a silence as a computation rather than as an omission. The
    lookback is recorded beside the whole-span reading: the whole span is the more
    permissive of the two, so a null across it implies a null under the four-month
    lookback the engine uses on a live series.
    """
    spec = policy["families"]["policy"]
    rules, gates = spec["rules_ocv"], spec["gates"]
    held = audit["held_out"]
    persist = rules["persist_months"]
    lookback = rules["lookback_months"]
    entries = []
    for rule, keys, floor in (("block_z", MAGNITUDES, rules["block_z"]),
                              ("dimension_z", plate.DIMS, rules["dimension_z"])):
        for key in keys:
            flags = [row["values"][key] is not None
                     and abs(row["values"][key]) > floor for row in held]
            run = run_of(held, flags)
            tail = held[-lookback:] if lookback > 0 else held
            tail_flags = flags[-lookback:] if lookback > 0 else flags
            tail_run = run_of(tail, tail_flags)
            values = [row["values"][key] for row in held
                      if row["values"][key] is not None]
            entries.append({
                "rule": rule,
                "quantity": key,
                "threshold": floor,
                "comparison": "abs(%s) > %g" % (key, floor),
                "max_abs": max((abs(value) for value in values), default=None),
                "n_over": sum(flags),
                "months_over": [row["label"] for row, flag in zip(held, flags)
                                if flag],
                "persist_months": persist,
                "longest_run": run,
                "fired": bool(run["length"] >= persist),
                "applied": True,
                "lookback": {"months": lookback,
                             "n_over": sum(tail_flags),
                             "longest_run": tail_run,
                             "fired": bool(tail_run["length"] >= persist)},
            })
    entries.append({
        "rule": "observability_drop",
        "quantity": "observability level",
        "threshold": rules["observability_drop"],
        "comparison": ("a fall of the observability level by %g"
                       % rules["observability_drop"]),
        "max_abs": None,
        "n_over": None,
        "months_over": [],
        "persist_months": persist,
        "longest_run": {"length": 0, "months": [], "start": None, "end": None},
        "fired": False,
        "applied": False,
        "reason": ("the rule reads a level of observability -- the eligible share "
                   "of a month and its echo share -- and this cut carries chips, "
                   "not levels: it can neither fire nor clear"),
    })
    return _gates(report, held, entries, spec, gates)


def _gates(report, held, entries, spec, gates):
    """Every quality gate v1.0.0 names, each read from the plate's own numbers."""
    cut = report["cut"]
    share = (cut["kept"] / cut["rows"]) if cut.get("rows") else None
    clipped = report.get("clipped")
    attempted = cut["kept"] * len(plate.DIMS) if cut.get("kept") else None
    clip_share = ((clipped / attempted)
                  if (clipped is not None and attempted) else None)
    degenerate = sorted({dim for cells in report["standardisation"].values()
                         for dim, cell in cells.items()
                         if not cell.get("scale")})
    complete = sum(1 for row in held if row["fielded"])

    def gate(name, value, threshold, passed, read_as):
        return {"gate": name, "value": value, "threshold": threshold,
                "read_as": read_as,
                "passed": None if value is None else bool(passed)}

    checks = [
        gate("min_months", len(held), gates["min_months"],
             len(held) >= gates["min_months"], "held-out bins"),
        gate("min_reference_months",
             report["reference_envelope"]["n_reference_months"],
             gates["min_reference_months"],
             report["reference_envelope"]["n_reference_months"]
             >= gates["min_reference_months"],
             "reference months inside the frozen window"),
        gate("min_eligible_share", share, gates["min_eligible_share"],
             share is not None and share >= gates["min_eligible_share"],
             "chips kept over chips offered"),
        gate("min_complete_months", complete, gates["min_complete_months"],
             complete >= gates["min_complete_months"],
             "held-out bins carrying all %s streams" % word_number(len(STREAMS))),
        gate("max_degenerate_dims", len(degenerate),
             gates["max_degenerate_dims"],
             len(degenerate) <= gates["max_degenerate_dims"],
             "dimensions with a missing or zero reference scale"),
        gate("max_clip_share", clip_share, gates["max_clip_share"],
             clip_share is not None and clip_share <= gates["max_clip_share"],
             "standardised values clipped at the plate's own clip"),
    ]
    rules = spec["rules_ocv"]
    params = [key for key in ("persist_months", "lookback_months") if key in rules]
    firing = [key for key in rules if key not in params]
    carried = [key for key in firing
               if key in {rule["rule"] for rule in entries if rule["applied"]}]
    missing = [key for key in firing if key not in carried]
    return {
        "policy": {"name": spec["name"], "version": spec["version"],
                   "sha256": (spec["files"][0]["sha256"] if spec.get("files")
                              else None),
                   "alert_enabled": spec["alert_enabled"],
                   "states": list(spec["states"]),
                   "operator_decision_required":
                       spec["operator_decision_required"]},
        "rules": entries,
        "gates": checks,
        "degenerate_dims": degenerate,
        "unreadable": [gate["gate"] for gate in checks
                       if gate["value"] is None],
        "not_applied": [rule["rule"] for rule in entries if not rule["applied"]],
        "firing_rules": firing,
        "parameter_rules": params,
        "scope": ("the policy's `rules.ocv` carries %s entries: %s firing rules "
                  "(%s) and %s parameters (%s). This cut applies %s of the firing "
                  "rules by value and cannot read %s, which wants a per-month "
                  "level and finds chips, so that rule is carried as missing "
                  "rather than as clear"
                  % (word_number(len(rules)),
                     word_number(len(firing)),
                     listed(firing),
                     word_number(len(params)),
                     listed(params),
                     word_number(len(carried)),
                     listed(missing))),
        "lookback_note": ("the whole held-out span is the more permissive of the "
                          "two readings, so a null across it implies a null "
                          "under the engine's %s-month lookback"
                          % word_number(rules["lookback_months"])),
        "vocabulary": list(LIVE_STATES),
    }


def monthly_rms_variant(report):
    """The block magnitudes recomputed from the **monthly medians**, for contrast.

    The plate's own ``R_B``/``R_G`` are per-chip RMS values whose monthly state is
    their median. This variant takes the RMS of the month's median z instead, which
    answers a slightly different question ("how big is the month's block vector").
    It is reported beside the per-chip series and replaces nothing: no state, no
    lead and no gate reads it.
    """
    rows = []
    for entry in report["months"]:
        if not (entry["n"] and entry.get("n_outside", 0) == entry["n"]):
            continue
        z = entry["z"]
        rows.append({
            "label": entry["label"],
            "R_B": plate.rms([z[dim] for dim in plate.BLOCKS["behaviour"]]),
            "R_G": plate.rms([z[dim] for dim in plate.BLOCKS["geometry"]])})
    result = {"values": rows,
              "kind": "rms of the monthly medians, beside the per-chip series",
              "reads": False,
              "replace": False}
    for key in MAGNITUDES:
        values = [row[key] for row in rows if row[key] is not None]
        result[key] = {"n": len(values),
                       "median": plate.median(values) if values else None,
                       "max": max(values) if values else None,
                       "min": min(values) if values else None}
    return result


def first_fire(rows, flags, persist):
    """The first row at which ``persist`` consecutive flags have accumulated."""
    run = 0
    for row, flag in zip(rows, flags):
        run = run + 1 if flag else 0
        if run >= persist:
            return row
    return None


def leads(audit, record, event):
    """The leads this run is willing to define, each dated to a month bin.

    Nothing here is a single number with a story attached. ``first_flagged`` is the
    first held-out month anything at all left the envelope (the *signed* lead).
    ``first_persistence`` is the first month at which any stream had already been
    outside for ``PERSIST_MIN`` consecutive months -- the earliest moment
    "sustained" could be said of anything. ``first_change`` is the first month at
    which the trajectory rule, read on the months up to and including it, returns
    ``CHANGE`` -- the *sustained* lead. The live rule's first firing month is
    recorded beside them, and is empty when the policy never fires.
    """
    held = audit["held_out"]
    persist = plate.PERSIST_MIN

    def dated(row):
        if row is None:
            return None
        return {"label": row["label"], "month": row["month"],
                "day": plate.bin_date(row).isoformat(),
                "days_before_event":
                    (dt.date.fromisoformat(event) - plate.bin_date(row)).days}

    flagged = next((row for row in held
                    if any(row["flags"][key]["fired"] for key in STREAMS)), None)
    first_persistence = first_change = None
    for index in range(1, len(held) + 1):
        rows = held[:index]
        runs = {key: run_of(rows, [bool(row["flags"][key]["fired"])
                                  for row in rows])["length"]
                for key in STREAMS}
        if first_persistence is None and max(list(runs.values()) + [0]) >= persist:
            first_persistence = rows[-1]
        if first_change is None and state_of_trajectory(
                {"fielded": all(row["fielded"] for row in rows), "runs": runs},
                persist) == "CHANGE":
            first_change = rows[-1]
    live_first, live_rule = None, None
    for rule in record["rules"]:
        if not rule["applied"]:
            continue
        flags = [row["values"][rule["quantity"]] is not None
                 and abs(row["values"][rule["quantity"]]) > rule["threshold"]
                 for row in held]
        hit = first_fire(held, flags, rule["persist_months"])
        if hit is not None and (live_first is None
                                or plate.bin_date(hit) < plate.bin_date(live_first)):
            live_first, live_rule = hit, "%s/%s" % (rule["rule"], rule["quantity"])
    return {
        "first_flagged": dated(flagged),
        "first_persistence": dated(first_persistence),
        "first_change": dated(first_change),
        "first_live_fire": dated(live_first),
        "first_live_rule": live_rule,
        "signed_lead_days": (dated(flagged) or {}).get("days_before_event"),
        "sustained_lead_days": (dated(first_change) or {}).get("days_before_event"),
        "anchor": ("the bin's own tick: the 15th of a whole month; the event is "
                   "%s" % event),
        "definitions": {
            "first_flagged": ("the first held-out month in which any of the %s "
                              "streams lies outside the frozen envelope"
                              % word_number(len(STREAMS))),
            "first_persistence": ("the first month at which some stream had "
                                  "already been outside for %d consecutive months"
                                  % persist),
            "first_change": ("the first month at which the trajectory rule, read "
                             "on the months up to it, returns CHANGE"),
            "first_live_fire": ("the first month at which a pinned policy rule "
                                "would have fired, both blocks and dimensions "
                                "considered"),
        },
    }


def critical_floor(audit, rules, persist):
    """For each applied rule, the largest floor at which that rule still fires.

    A rule of the form "``|x| > f`` for ``persist`` consecutive months" fires exactly
    when some run of that many consecutive held-out months has its *smallest*
    ``|x|`` above ``f``. So the largest floor that still fires is the largest such
    minimum over windows, and the distance between the pinned floor and that number
    is what separates silence from firing -- a quantity, not a phrase. Each rule is
    measured against its own quantity and its own floor; the closest of them decides
    the word.
    """
    held = audit["held_out"]
    out = []
    for rule in rules:
        key = rule["quantity"]
        values = [abs(row["values"][key]) if row["values"][key] is not None
                  else 0.0 for row in held]
        best = None
        for start in range(0, max(1, len(values) - persist + 1)):
            window = values[start:start + persist]
            if len(window) < persist:
                continue
            floor = min(window)
            if best is None or floor > best["floor"]:
                best = {"floor": floor, "months": [row["label"] for row in
                                                   held[start:start + persist]]}
        if best is None:
            out.append({"rule": rule["rule"], "quantity": key,
                        "policy_floor": rule["threshold"], "fires_at": None,
                        "months": [], "margin": None,
                        "fires_under_policy": False})
            continue
        out.append({"rule": rule["rule"], "quantity": key,
                    "policy_floor": rule["threshold"],
                    "fires_at": best["floor"], "months": best["months"],
                    "margin": rule["threshold"] - best["floor"],
                    "fires_under_policy": bool(rule["threshold"] <= best["floor"])})
    margins = [entry["margin"] for entry in out if entry["margin"] is not None]
    closest = min(out, key=lambda entry: (entry["margin"] is None,
                                          entry["margin"])) if out else None
    return {"persist_months": persist,
            "rules": out,
            "closest": closest,
            "any_fires": any(entry["fires_under_policy"] for entry in out),
            "smallest_margin": min(margins) if margins else None,
            "plate_persist_months": plate.PERSIST_MIN}


def readings(report, policy, event):
    """Everything the two readers say about one built report. Pure, no writing.

    The control goes through this same function, which is what "same layer" means
    concretely: one code path, one envelope construction, one persistence rule, two
    assets. Nothing in here knows which site it was handed except through the
    numbers and the event string used for the day counts.
    """
    audit = phase_audit(report)
    st = streams(audit)
    record = live_record(report, audit, policy)
    live_state = state_of_live(record)
    record["state"] = live_state
    record["reason"] = live_reason(record, live_state)
    record["over_all_five"] = {
        "state": ("INSUFFICIENT" if record["not_applied"] else live_state),
        "reason": ("%s cannot be read from this cut"
                   % ", ".join(record["not_applied"])) if record["not_applied"]
                  else "every rule of the policy could be read",
        "note": ("the same policy read over every firing rule it names. The state "
                 "above is the policy's word restricted to the rules an "
                 "observability cut carries; this one is what v1.0.0 would say if "
                 "the observability level were missing from the record, which it "
                 "is. Both are printed because only one of them is a clearance and "
                 "neither is a verdict on the asset."),
    }
    runs = st["runs"]
    persist = plate.PERSIST_MIN
    fielded = audit["fully_fielded"]
    traj_state = state_of_trajectory({"fielded": fielded, "runs": runs}, persist)
    basis = [{"rule": "block_run", "quantity": key, "value": runs[key],
              "threshold": persist, "fired": bool(runs[key] >= persist)}
             for key in MAGNITUDES]
    basis += [{"rule": "dimension_run", "quantity": dim, "value": runs[dim],
               "threshold": persist, "fired": bool(runs[dim] >= persist)}
              for dim in plate.DIMS]
    assessment = report.get("assessment") or {}
    verdict = assessment.get("verdict") if assessment.get("applies", True) else None
    trajectory = {"state": traj_state, "basis": basis,
                  "persist_months": persist,
                  "fielded": fielded,
                  "rule": ("CHANGE needs a sustained block magnitude and a "
                           "sustained dimension; WATCH either one alone; STABLE "
                           "neither; NOT_ASSESSABLE a held-out bin without all "
                           "%s streams" % word_number(len(STREAMS))),
                  "vocabulary": list(TRAJECTORY_STATES),
                  "plate_verdict": verdict}
    traj_positive = traj_state in ("CHANGE", "WATCH")
    live_positive = live_state == "WATCH"
    pair = {
        "trajectory": traj_state,
        "live": live_state,
        "live_over_all_five": record["over_all_five"]["state"],
        "trajectory_positive": traj_positive,
        "live_positive": live_positive,
        "concur": traj_positive == live_positive,
        "word": ("both readers see a departure" if traj_positive and live_positive
                 else "neither reader sees a departure"
                 if not traj_positive and not live_positive else
                 "the readers differ: %s against %s" % (traj_state, live_state)),
        "collapsed": False,
        "note": ("two vocabularies over one set of numbers. The trajectory word "
                 "summarises the frozen envelope's own persistence rule; the live "
                 "word is the pinned policy's, over its own floors. Neither is a "
                 "translation of the other and no sentence here chooses between "
                 "them. The live word is read over the rules this cut carries; the "
                 "same policy over every firing rule it names reads "
                 "%s." % record["over_all_five"]["state"]),
    }
    return {
        "site": report["site"],
        "label": report["label"],
        "event": report["event"],
        "reference": {"window": report["reference"]["chosen"],
                      "n_chips": report["reference"]["span"]["n_chips"],
                      "months": report["reference"]["span"]["months"],
                      "choice": report["reference"]["choice"],
                      "choice_rule": report["reference"].get("choice_rule")},
        "cut": report["cut"],
        "partition": partition(report),
        "audit": audit,
        "streams": st,
        "variant": monthly_rms_variant(report),
        "trajectory": trajectory,
        "live": record,
        "pair": pair,
        "decisive": critical_floor(
            audit, [rule for rule in record["rules"] if rule["applied"]],
            record["rules"][0]["persist_months"]),
        "leads": leads(audit, record, event),
        "envelope": {"n_reference_months":
                     report["reference_envelope"]["n_reference_months"],
                     "aggregation": report["reference_envelope"]["aggregation"],
                     "quantiles": report["reference_envelope"]["quantiles"],
                     "decision_quantile":
                         report["reference_envelope"]["decision_quantile"],
                     "decision_level":
                         report["reference_envelope"]["decision_level"],
                     "months": report["reference_envelope"]["months"],
                     "chip_note": report["reference_envelope"]["chip_note"]},
    }


def frozen_numbers(report):
    """The numbers the frozen window owns -- what no post-event value may move.

    The held-out bins are keyed by calendar month rather than by the plate's label,
    because a month that the event splits carries a ``(pre)``/``(post)`` suffix that
    belongs to the split, not to the numbers. Everything else here is the reference
    side: its census, its span, its monthly envelope and its standardisation.
    """
    monthly = report["reference_envelope"]["monthly"]
    return {
        "reference": report["classes"]["reference"],
        "reference_span": report["reference"]["span"],
        "envelope": {key: {q: monthly[key][q] for q in ("q05", "q50", "q95")}
                     for key in STREAMS},
        "standardisation": report["standardisation"],
        "outside": {row["month"]: {"label": row["label"], "n": row["n"],
                                   "values": {key: plate.bin_value(row, key)
                                              for key in STREAMS}}
                    for row in report["months"]
                    if row["n"] and row.get("n_outside", 0) == row["n"]},
    }


def eas_root():
    """Where the pinned EAS files live: a sibling repository, or ``None``."""
    root = os.path.abspath(os.path.join(REPO, os.pardir, os.pardir))
    return root if os.path.isdir(os.path.join(root, "erlking_mvp")) else None


def pin(path=None, strict=False):
    """Read the pin, check every digest, and return the record with the verdict.

    The pin travels *with* the artefact: the report can then say which policy it
    read and prove it read that one, rather than naming a version in prose. If the
    sibling repository is not present the digests are reported as unverified and
    the build continues -- a reader without the EAS checkout can still reproduce
    every number here from the pinned values. A digest that is present and
    **different** is a different matter: that raises, because a moved threshold
    must fail the build rather than quietly move the verdict.
    """
    path = path or PIN
    if not os.path.exists(path):
        raise SystemExit("%s: no policy pin -- nothing to read" % path)
    with open(path) as handle:
        record = json.load(handle)
    record["file"] = os.path.relpath(path, REPO)
    record["sha256"] = plate.sha256_of(path)
    root = eas_root()
    digests, problems = [], []
    for entry in record["families"]["policy"]["files"]:
        full = os.path.join(root, entry["path"]) if root else None
        if full is None or not os.path.exists(full):
            digests.append({"path": entry["path"], "matched": None,
                            "state": "not present beside this repository"})
            continue
        actual = plate.sha256_of(full)
        matched = actual == entry["sha256"]
        digests.append({"path": entry["path"], "sha256": actual,
                        "pinned": entry["sha256"], "matched": matched,
                        "state": "verified" if matched else "MOVED"})
        if not matched:
            problems.append("%s: digest %s is not the pinned %s"
                            % (entry["path"], actual[:16], entry["sha256"][:16]))
    record["digests"] = digests
    record["problems"] = problems
    record["n_verified"] = sum(1 for digest in digests if digest["matched"])
    record["verified"] = bool(digests) and not problems and all(
        digest["matched"] for digest in digests)
    if problems and strict:
        raise SystemExit("%s: the pinned policy has moved:\n  %s"
                         % (record["file"], "\n  ".join(problems)))
    return record


def _without_labels(bins):
    """Held-out bins without the plate's label, which the split -- not a number -- owns."""
    return {month: {"n": cell["n"], "values": cell["values"]}
            for month, cell in bins.items()}


def hindsight(resolved, baseline):
    """The guards that make the word "prospective" checkable.

    * **Delete the tail.** Every chip on or after the event is removed and the
      report is rebuilt. Every number the frozen window owns -- the reference
      census, the envelope, the standardisation, and every held-out month that
      still exists -- must come back bit-identical, because none of them was ever
      computed from the tail. Bins that consisted *entirely* of deleted chips
      disappear whole; they are listed by name rather than counted as a change.
    * **Make the tail absurd.** The post chips' measured values are replaced by
      ``1e6`` and the report is rebuilt. The same numbers must still come back
      identical: a value the reference never saw cannot move a frozen quantile.
    * **Recompute the reference without the event.** The reference set is
      re-selected from the frozen window and the registry's reference state alone,
      never consulting the event, and it must match the report's own reference
      census. That is what shows the event to be a boundary on the held-out side
      and nothing more.
    """
    chips, _cut = plate.load_chips(CASE)
    window = tuple(resolved["window"])
    event = dt.date.fromisoformat(plate.SITES[CASE]["event"])
    before = [chip for chip in chips if chip["day"] < event]
    farcical = [dict(chip) for chip in chips]
    for chip in farcical:
        if chip["day"] >= event:
            chip["values"] = {dim: 1e6 for dim in plate.DIMS}
            chip["echo"] = True
    base = frozen_numbers(baseline)
    base_fixed = {key: value for key, value in base.items() if key != "outside"}
    base_out = base["outside"]
    report = {}
    with plate.window_override(CASE, window, WINDOW_RULE["statement"]):
        for name, subset in (("tail_deleted", before), ("tail_absurd", farcical)):
            numbers = frozen_numbers(plate.build(CASE, chips=subset))
            fixed = {key: value for key, value in numbers.items()
                     if key != "outside"}
            shared = {month: cell for month, cell in numbers["outside"].items()
                      if month in base_out}
            expected = {month: cell for month, cell in base_out.items()}
            reference_difference = plate.first_difference(base_fixed, fixed)
            held_out_difference = plate.first_difference(
                _without_labels(expected), _without_labels(shared))
            report[name] = {
                "n_chips_dropped": len(chips) - len(subset),
                "reference_difference": reference_difference,
                "held_out_difference": held_out_difference,
                "vanished_bins": sorted(set(base_out) - set(numbers["outside"])),
                "labels_moved": {
                    month: {"before": base_out[month]["label"],
                            "after": numbers["outside"][month]["label"]}
                    for month in expected
                    if month in numbers["outside"]
                    and numbers["outside"][month]["label"]
                    != base_out[month]["label"]},
                "identical": (reference_difference is None
                              and held_out_difference is None),
            }
        by_window = [chip for chip in chips if plate.in_reference(CASE, chip)]
        outside = [chip for chip in chips
                   if plate.chip_class(CASE, chip) == "outside"]
    in_report = baseline["classes"]["reference"]["n_chips"]
    report["reference_from_window_alone"] = {
        "n_chips": len(by_window),
        "in_report": in_report,
        "identical": len(by_window) == in_report,
        "reads": ("the frozen window and the registry's reference state; never "
                  "the event"),
    }
    latest = max(chip["day"].isoformat() for chip in outside)
    report["held_out_is_pre_event"] = {
        "passed": bool(latest < event.isoformat()),
        "n_held_out_chips": len(outside),
        "latest_held_out": latest,
        "event": event.isoformat(),
    }
    report["passed"] = bool(
        report["tail_deleted"]["identical"]
        and report["tail_absurd"]["identical"]
        and report["reference_from_window_alone"]["identical"]
        and report["held_out_is_pre_event"]["passed"])
    return report


def _vocabulary_scan(case, control, policy):
    """Every state word either reader emitted, tested against its own vocabulary.

    The two readers borrow their words from two documents and mean different
    things by them, so this is a membership test in both directions: a
    trajectory word must lie in the trajectory rule's four, a live word must lie
    in the policy's own list *and* in the subset an observability cut can print,
    and neither reader may emit a word the other owns alone. ``WATCH`` is the one
    word both vocabularies define; the coincidence is reported because it is the
    single place a reader could mistake one word for the other.
    """
    declared = list(policy["states"])
    trajectory = [case["trajectory"]["state"], control["trajectory"]["state"]]
    live = [case["live"]["state"], case["live"]["over_all_five"]["state"],
            control["live"]["state"]]
    stray_trajectory = sorted({word for word in trajectory
                               if word not in TRAJECTORY_STATES})
    stray_live = sorted({word for word in live
                         if word not in declared or word not in LIVE_STATES})
    shared = sorted(set(TRAJECTORY_STATES) & set(LIVE_STATES))
    return {
        "passed": bool(not stray_trajectory and not stray_live
                       and set(live) <= set(declared)),
        "trajectory": list(TRAJECTORY_STATES),
        "live": list(LIVE_STATES),
        "declared": declared,
        "emitted": {"trajectory": sorted(set(trajectory)),
                    "live": sorted(set(live))},
        "stray": {"trajectory": stray_trajectory, "live": stray_live},
        "shared": shared,
        "evidence": ("both readers, case and control, each word in its own "
                     "reader's list"),
    }


def guards(case, control, resolved, hindsight_report, policy, influence,
           magnitudes):
    """Every guard the run promises, computed rather than asserted in prose."""
    part = case["partition"]
    audit = case["audit"]
    st = case["streams"]
    traj, live = case["trajectory"], case["live"]
    reruns = state_of_trajectory({"fielded": audit["fully_fielded"],
                                  "runs": st["runs"]}, plate.PERSIST_MIN)
    live_again = state_of_live(live)
    max_stream = max(list(st["runs"].values()) + [0])
    max_basis = max([entry["value"] for entry in traj["basis"]] + [0])
    window_code = _code_of(resolve_window)
    prose = json.dumps([case["pair"], traj["state"], live["state"],
                        case["leads"], part], sort_keys=True).lower()
    report = {
        "window_event_free": {
            "passed": bool(resolved.get("guard_passed")
                           and "event_day" not in window_code
                           and "SITES" not in window_code),
            "days_before_event": resolved.get("days_before_event"),
            "month_gap": resolved.get("month_gap"),
            "guard_months": GUARD_MONTHS,
            "rule": resolved["rule"]["statement"],
            "evidence": ("the resolver reads acquisition months and per-month "
                         "chip counts; the event arrives only as the date the "
                         "chosen window is checked against, after it is fixed"),
        },
        "buckets_disjoint_and_reconciled": {
            "passed": bool(part["reconciles"]),
            "total": part["total"],
            "kept": part["kept"],
            "reconciliation": part["reconciliation"],
            "rule": part["rule"],
        },
        "no_pooling": {
            "passed": bool(case["site"] != control["site"]),
            "evidence": ("one build per site; the control never enters the case's "
                         "numbers, which is checked by rebuilding the case after "
                         "the control has been read"),
        },
        "envelope_from_reference_months_only": {
            "passed": bool(case["envelope"]["n_reference_months"]
                           == resolved["n_months"]),
            "n_reference_months": case["envelope"]["n_reference_months"],
            "aggregation": case["envelope"]["aggregation"],
        },
        "every_held_out_bin_fielded": {
            "passed": bool(audit["fully_fielded"]),
            "unfielded": audit["unfielded"],
            "n_held_out": audit["n_held_out"],
        },
        "state_recomputes_from_stored_numbers": {
            "passed": bool(reruns == traj["state"] and live_again == live["state"]),
            "trajectory": {"stored": traj["state"], "recomputed": reruns},
            "live": {"stored": live["state"], "recomputed": live_again},
        },
        "record_is_max_over_streams": {
            "passed": bool(max_stream == max_basis),
            "longest_run": max_stream,
            "stream": max(st["runs"], key=lambda key: st["runs"][key]),
        },
        "control_after_case_no_influence": {
            "passed": influence is None,
            "difference": influence,
            "order": [CASE, CONTROL, CASE],
        },
        "runs_reported_separately_from_magnitude": {
            "passed": all("run_length" in cell and "median" in cell
                          for cell in st["streams"].values()),
            "evidence": ("every stream carries its longest run beside its own "
                         "median and extremes; no sentence substitutes one for "
                         "the other"),
        },
        "magnitudes_at_1e_12": dict(magnitudes),
        "vocabulary_scan": _vocabulary_scan(case, control, live["policy"]),
        "no_failure_language": {
            "passed": not any(word in prose for word in
                              ("fail", "failure", "failing")),
            "scanned": "the pair, both states, the leads and the partition",
        },
        "hindsight": dict(hindsight_report),
        "thresholds_pinned": {
            "passed": bool(not policy["problems"]),
            "pin": policy["file"],
            "sha256": policy["sha256"],
            "n_verified": policy["n_verified"],
            "values": dict(policy["families"]["policy"]["rules_ocv"]),
            "gates": dict(policy["families"]["policy"]["gates"]),
        },
        "event_is_the_registrys": {
            "passed": bool(case["event"] == plate.SITES[CASE]["event"]),
            "event": case["event"],
        },
        "policy_sha_matches_pin": {
            "passed": bool(policy["problems"] == []),
            "problems": policy["problems"],
            "digests": policy["digests"],
        },
    }
    report["passed"] = all(report[name].get("passed") is not False
                           for name in GUARD_ORDER)
    return report


def _code_of(function):
    """A function's source with its docstring removed -- code, not prose."""
    text = inspect.getsource(function)
    for line in (function.__doc__ or "").splitlines():
        stripped = line.strip()
        if stripped:
            text = text.replace(stripped, "")
    return text


def magnitude_identity(site):
    """``R_B = rms(gamma2, P)`` and ``R_G = rms(A, D, F, S)``, recomputed per chip.

    The block magnitudes are the ones the plate commits; this run adds no metric,
    so the identity has to hold to floating-point noise on every chip of the cut.
    The worst difference is recorded, not a rounded claim that it holds.
    """
    chips, _cut = plate.load_chips(site)
    plate.standardise(site, chips)
    worst, count = 0.0, 0
    for chip in chips:
        for key, block in (("R_B", "behaviour"), ("R_G", "geometry")):
            want = plate.rms([chip["z"][dim] for dim in plate.BLOCKS[block]])
            have = chip.get(key)
            if want is None or have is None:
                continue
            worst = max(worst, abs(want - have))
            count += 1
    return {"passed": bool(count and worst <= 1e-12),
            "identity": "R_B = rms(gamma2, P); R_G = rms(A, D, F, S)",
            "n_chips_checked": count,
            "max_abs_difference": worst,
            "tolerance": 1e-12}


def comparison(case, control):
    """The two assets side by side under the same layer, never pooled."""
    rows = []
    for key in STREAMS:
        one = case["streams"]["streams"][key]
        other = control["streams"]["streams"][key]
        rows.append({
            "stream": key,
            "case": {"run": one["run_length"], "n_fired": one["n_fired"],
                     "n_months": one["n_months"], "median": one["median"],
                     "max_abs": one["max_abs"], "q05": one["q05"],
                     "q95": one["q95"]},
            "control": {"run": other["run_length"], "n_fired": other["n_fired"],
                        "n_months": other["n_months"], "median": other["median"],
                        "max_abs": other["max_abs"], "q05": other["q05"],
                        "q95": other["q95"]},
        })
    return {
        "order": [case["site"], control["site"]],
        "rows": rows,
        "window_source": {
            case["site"]: ("the stated rule selected it: %s"
                           % case["reference"]["choice_rule"]),
            control["site"]: ("the registry declares it, and its descriptive "
                              "plate is published against it: %s"
                              % control["reference"]["choice"]),
        },
        "case": {"site": case["site"], "event": case["event"],
                 "pair": case["pair"], "leads": case["leads"],
                 "trajectory": case["trajectory"]["state"],
                 "plate_verdict": case["trajectory"]["plate_verdict"],
                 "live_reason": case["live"]["reason"]},
        "control": {"site": control["site"], "event": control["event"],
                    "pair": control["pair"], "leads": control["leads"],
                    "trajectory": control["trajectory"]["state"],
                    "plate_verdict": control["trajectory"]["plate_verdict"],
                    "live_reason": control["live"]["reason"]},
        "differs": [row["stream"] for row in rows
                    if row["case"]["run"] != row["control"]["run"]],
        "note": ("Read in this order: the case, then the control, then the case "
                 "again. The third reading is compared with the first and must be "
                 "identical; the control is never an input to the case. The two "
                 "assets are not pooled and no number here is an average of them."),
    }


def build_all(policy=None):
    """The whole artefact: the case, the control, their comparison, the guards."""
    policy = policy or pin(strict=True)
    case_event = plate.SITES[CASE]["event"]
    resolved = resolve_window(CASE, event=case_event)
    report = build_case(CASE, resolved["window"])
    case = readings(report, policy, case_event)

    control_event = plate.SITES[CONTROL]["event"]
    control = readings(build_case(CONTROL), policy, control_event)

    again = readings(build_case(CASE, resolved["window"]), policy, case_event)
    influence = plate.first_difference(plate.round_floats(case),
                                       plate.round_floats(again))
    hinds = hindsight(resolved, report)
    magnitudes = magnitude_identity(CASE)
    cmp = comparison(case, control)
    paths = out_paths()
    return {
        "schema": SCHEMA,
        "run": {
            "name": "frozen-reference prospective read",
            "case": CASE,
            "control": CONTROL,
            "readers": list(READERS),
            "co_primary": True,
            "collapsed": False,
            "case_event_source": "SITES registry",
            "control_window_source": "SITES registry (declared, not re-selected)",
            "case_window_source": "WINDOW_RULE, resolved before any read",
            "reads_the_tail": ("only in the hindsight guards, and only to delete "
                               "it or to make it absurd"),
        },
        "window": resolved,
        "case": case,
        "control": control,
        "comparison": cmp,
        "pin": {key: policy.get(key) for key in
                ("file", "sha256", "verified", "n_verified", "problems",
                 "digests", "created", "purpose", "application")},
        "guards": guards(case, control, resolved, hinds, policy, influence,
                         magnitudes),
        "descriptive": {
            "note": ("the case's registry window and its published plate are "
                     "untouched by this run"),
            "json": os.path.relpath(plate.out_paths(CASE)["json"], REPO),
            "md": os.path.relpath(plate.out_paths(CASE)["md"], REPO),
            "jpg": os.path.relpath(plate.out_paths(CASE)["jpg"], REPO),
        },
        "outputs": {key: os.path.relpath(value, REPO)
                    for key, value in paths.items()},
    }


def out_paths():
    """Every file this run writes, by role."""
    stem = "fig_%s_ocv_trajectory_prospective" % CASE
    return {
        "json": os.path.join(DATA, CASE, stem + ".json"),
        "md": os.path.join(FIGURES, CASE, stem + ".md"),
        "jpg": os.path.join(FIGURES, CASE, stem + ".jpg"),
        "sh": os.path.join(FIGURES, CASE, stem + ".sh"),
        "cmp_json": os.path.join(DATA, "cases",
                                 "fig_ocv_prospective_comparison.json"),
        "cmp_md": os.path.join(FIGURES, "cases",
                               "fig_ocv_prospective_comparison.md"),
    }


def listed(names):
    """``['a'] -> '`a`'``, ``['a','b'] -> '`a` and `b`'``, Oxford-less comma list."""
    names = ["`%s`" % name for name in names]
    if not names:
        return "no stream"
    if len(names) == 1:
        return names[0]
    return "%s and %s" % (", ".join(names[:-1]), names[-1])


def f3(value, digits=3):
    """A number for the page, or an en dash when the run produced none."""
    return "—" if value is None else ("%.*f" % (digits, value))


def label_of(cell):
    return cell["label"] if cell else "—"


def day_of(cell):
    return "%s (%s, %d days)" % (cell["label"], cell["day"],
                                 cell["days_before_event"]) if cell else "did not happen"


def prospective_md(artefact):
    """Every answer, each beside the numbers it was read from."""
    a = artefact
    case, control = a["case"], a["control"]
    window, part = a["window"], case["partition"]["buckets"]
    st = case["streams"]["streams"]
    live, traj, leads = case["live"], case["trajectory"], case["leads"]
    rule = window["rule"]
    out = []
    add = out.append
    add("# Ponte Morandi, Genova — the prospective read, two readers over one "
        "frozen reference")
    add("")
    add("**Case** %s · **event** %s · **frozen window** %s → %s (**%d** chips, "
        "%d months, selected by a rule written before the read) · **held out** "
        "%s → %s (**%d** chips, %d bins) · **control** %s, read after the case · "
        "**policy** erlking-eas %s, pinned (`data/cases/eas_policy_pin.json`, "
        "sha256 `%s`)"
        % (case["label"], case["event"], window["window"][0], window["window"][1],
           window["n_chips"], window["n_months"],
           part["outside"]["window"][0], part["outside"]["window"][1],
           part["outside"]["n_chips"], part["outside"]["months"],
           control["label"], live["policy"]["version"],
           a["pin"]["sha256"][:16]))
    add("")
    add("```")
    add("window rule   %s" % rule["statement"])
    add("              reads %s" % "; ".join(rule["reads"]))
    add("              never reads %s" % "; ".join(rule["never_reads"]))
    add("envelope      monthly medians of the %d frozen reference months, "
        "q05/q50/q95 (empirical)   decision quantile q%02d"
        % (case["envelope"]["n_reference_months"],
           case["envelope"]["decision_quantile"]))
    entry = {rule["rule"]: rule for rule in live["rules"]}
    add("persistence   trajectory rule %d consecutive months   policy v%s %d "
        "consecutive months   lookback %d months"
        % (traj["persist_months"], live["policy"]["version"],
           entry["block_z"]["persist_months"],
           entry["block_z"]["lookback"]["months"]))
    add("policy floors block |z| > %g   dimension |z| > %g   alert_enabled %s"
        % (entry["block_z"]["threshold"], entry["dimension_z"]["threshold"],
           live["policy"]["alert_enabled"]))
    add("readers       %s  -- co-primary, never collapsed"
        % " and ".join(a["run"]["readers"]))
    add("```")
    add("")
    add("Two readers look at the same %s held-out monthly states: the "
        "trajectory rule of this repository, whose boundary is the reference "
        "months' own q05/q95, and EAS policy v%s, whose boundary is a fixed "
        "`|z|` floor. Both use the same persistence rule -- %d consecutive months "
        "-- so the whole difference between their words is a distance between two "
        "boundaries, and that distance is measured in Q8."
        % (word_number(part["outside"]["months"]), live["policy"]["version"],
           traj["persist_months"]))
    add("")
    _md_window(add, a)
    _md_partition(add, a)
    _md_trajectory(add, a)
    _md_sustain(add, a)
    _md_leads(add, a)
    _md_policy(add, a)
    _md_pair(add, a)
    _md_why(add, a)
    _md_guards(add, a)
    _md_control(add, a)
    _md_selection(add, a)
    _md_counterfactual(add, a)
    _md_limits(add, a)
    _md_close(add, a)
    return "\n".join(out).rstrip() + "\n"


def _md_window(add, a):
    w = a["window"]
    add("## Q1 · What window was frozen, and by what rule?")
    add("")
    add("%s, written before the read ran: **%s**. The rule reads %s and nothing "
        "else — in particular it never reads %s, so the size of the held-out "
        "spread cannot have steered it."
        % (w["rule"]["name"].capitalize(), w["rule"]["statement"],
           " and ".join(w["rule"]["reads"]),
           " or ".join(w["rule"]["never_reads"])))
    add("")
    add("It was tried from the first occupied month and kept the first start that "
        "held:")
    add("")
    add("| start tried | would end | months under %d chips | verdict |"
        % w["rule"]["min_chips"])
    add("|---|---|---|---|")
    for tried in w["starts_tried"]:
        add("| %s | %s | %s | %s |"
            % (tried["start"], tried["end"],
               ", ".join("%s (%d)" % (cell["month"], cell["chips"])
                         for cell in tried["shorts"]) or "none",
               "kept" if tried["short"] == [] else "skipped"))
    add("")
    add("So the reference is **%s → %s**, %d chips over %d calendar months:"
        % (w["window"][0], w["window"][1], w["n_chips"], w["n_months"]))
    add("")
    add("| month | " + " | ".join(cell["month"] for cell in w["months"]) + " |")
    add("|---|" + "---|" * len(w["months"]))
    add("| chips | " + " | ".join(str(cell["chips"]) for cell in w["months"]) + " |")
    add("")
    add("The guard that keeps this prospective rather than retrospective: the "
        "window ends **%d calendar months** and **%d days** before the event "
        "month's first day, against a %d-month guard — %s. A rule of 365 days "
        "would have refused this window by %d days while the event month is %d "
        "bins away, which is why the guard is counted in months."
        % (w["month_gap"], w["days_before_event"], w["guard_months"],
           "passed" if w["guard_passed"] else "failed",
           365 - w["days_before_event"], w["month_gap"]))
    add("")
    add("The registry declares no window for this site — its reference state "
        "spans %d chips over %d months and stops at the event — so this window "
        "*narrows* what was published rather than contradicting it. Q11 checks "
        "the selection against the data it was not allowed to see."
        % (a["case"]["reference"]["n_chips"], a["case"]["reference"]["months"]))
    add("")


def _md_partition(add, a):
    part = a["case"]["partition"]
    add("## Q2 · What was held out, and where did every other chip go?")
    add("")
    add("| bucket | what it is | chips | months | window |")
    add("|---|---|---|---|---|")
    for key, cell in part["buckets"].items():
        add("| %s | %s | %d | %d | %s |"
            % (key, cell["label"] or key, cell["n_chips"], cell["months"],
               ("%s → %s" % tuple(cell["window"])) if cell["window"] else "—"))
    add("| **total** | | **%d** | | |" % part["total"])
    add("")
    add("Disjoint and complete: %s. Every kept chip sits in exactly one bucket, "
        "chips after the event are kept but never read as evidence, and the %s "
        "chip%s before the window %s named rather than discarded silently."
        % (part["rule"],
           word_number(part["buckets"]["prefix"]["n_chips"]),
           "" if part["buckets"]["prefix"]["n_chips"] == 1 else "s",
           "is" if part["buckets"]["prefix"]["n_chips"] == 1 else "are"))
    add("")


def _md_trajectory(add, a):
    st = a["case"]["streams"]["streams"]
    held = a["case"]["audit"]["held_out"]
    add("## Q3 · Did either reader see anything before the event?")
    add("")
    add("The trajectory rule first, because it is the reader with the tighter "
        "boundary. Each of the %s streams is standardised against the frozen "
        "reference months and compared with its own q05/q95 envelope, over the "
        "%d held-out bins (%s → %s, %d chips, every one of them before the event):"
        % (word_number(len(a["case"]["streams"]["order"])),
           len(held), held[0]["label"], held[-1]["label"],
           a["case"]["partition"]["buckets"]["outside"]["n_chips"]))
    add("")
    add("| stream | months outside | longest run | run span | median z | min → max | "
        "envelope q05 → q95 |")
    add("|---|---|---|---|---|---|---|")
    for key in a["case"]["streams"]["order"]:
        cell = st[key]
        run = cell["longest_run"]
        add("| `%s` | %d / %d | %d | %s | %s | %s → %s | %s → %s |"
            % (key, cell["n_fired"], cell["n_months"], run["length"],
               ("%s → %s" % (run["start"], run["end"])) if run["length"] else "—",
               f3(cell["median"], 2), f3(cell["min"], 2), f3(cell["max"], 2),
               f3(cell["q05"], 2), f3(cell["q95"], 2)))
    add("")
    order = a["case"]["streams"]["order"]
    ranked = sorted(order, key=lambda key: -st[key]["longest_run"]["length"])
    longest = ", ".join("%d on `%s`" % (st[key]["longest_run"]["length"], key)
                        for key in ranked[:4])
    quiet = [key for key in order if not st[key]["longest_run"]["length"]]
    add("Longest runs, in held-out months: %s. They sit inside one span rather "
        "than scattering across the year, and the signs disagree: `gamma2` "
        "falls below q05 from %s while `A` and `F` rise above q95 from %s and "
        "%s — one departure seen on both sides of the envelope. %s never leaves "
        "its envelope at all, and the quiet stream is named here so that the "
        "departure is not read as a property of the whole asset."
        % (longest, st["gamma2"]["longest_run"]["start"],
           st["A"]["longest_run"]["start"], st["F"]["longest_run"]["start"],
           listed(quiet)))
    add("")


def _md_sustain(add, a):
    st = a["case"]["streams"]["streams"]
    traj = a["case"]["trajectory"]
    persist = traj["persist_months"]
    add("## Q4 · Was it sustained, or a single month?")
    add("")
    add("The rule requires **%d consecutive months** outside the envelope, in a "
        "block magnitude and in a dimension, before it will say CHANGE; one of the "
        "two alone is WATCH:" % persist)
    add("")
    add("| stream | block or dimension | longest run | vs %d | sustained | fired |"
        % persist)
    add("|---|---|---|---|---|---|")
    for row in traj["basis"]:
        cell = st[row["quantity"]]
        add("| `%s` | %s | %d | %s | %s | %s |"
            % (row["quantity"],
               "block" if row["rule"].startswith("block") else "dimension",
               cell["longest_run"]["length"],
               "≥" if cell["longest_run"]["length"] >= persist else "<",
               "yes" if cell["longest_run"]["length"] >= persist else "no",
               "**yes**" if row["fired"] else "no"))
    add("")
    rg = st["R_G"]
    singles = [key for key in a["case"]["streams"]["order"]
               if st[key]["n_fired"] == 1]
    lone_months = sorted({month for key in singles
                          for month in st[key]["months_fired"]})
    add("Every run is a run of *consecutive held-out calendar months*, not a "
        "count of flagged months: `R_G` leaves the envelope %d times, but its "
        "longest run without a gap covers %s → %s; %s %s it once each (%s). Read "
        "up to and including the last held-out month, the rule returns **%s** — "
        "*%s*."
        % (rg["n_fired"], rg["longest_run"]["start"], rg["longest_run"]["end"],
           listed(singles), "leave" if len(singles) != 1 else "leaves",
           ", ".join(lone_months), traj["state"], traj["plate_verdict"]))
    add("")


def _md_leads(add, a):
    cells = a["case"]["leads"]
    add("## Q5 · How much lead, and under which definition?")
    add("")
    add("%s definitions, %s dates. They are printed side by side because "
        "'lead time' without a definition is a slogan:"
        % (sentence_number(len(cells["definitions"])),
           word_number(len(cells["definitions"]))))
    add("")
    add("| definition | first month | bin tick | days before the event |")
    add("|---|---|---|---|")
    for key in ("first_flagged", "first_persistence", "first_change",
                "first_live_fire"):
        cell = cells[key]
        add("| %s | %s | %s | %s |"
            % (key, label_of(cell), cell["day"] if cell else "—",
               cell["days_before_event"] if cell else "did not happen"))
    add("")
    add("The signed lead is **%d days** (first month outside the envelope) and "
        "the sustained lead is **%d days** (first month read as CHANGE); the "
        "pinned policy never fires, so it has no lead to report. Dates come from "
        "%s. The month is the anchor, not the chip — a monthly state is dated by "
        "its tick, and dating it by the last chip in the bin would buy or lose "
        "days that no measurement supports."
        % (cells["signed_lead_days"], cells["sustained_lead_days"],
           cells["anchor"]))
    add("")
    add("The definitions the rows above are computed from:")
    add("")
    for key in ("first_flagged", "first_persistence", "first_change",
                "first_live_fire"):
        add("* `%s` — %s." % (key, cells["definitions"][key]))
    add("")


def _md_policy(add, a):
    live = a["case"]["live"]
    persist = live["rules"][0]["persist_months"]
    add("## Q6 · Did the pinned policy fire?")
    add("")
    add("Policy `%s` v%s, pinned by digest (`%s…`) and applied by value: every "
        "threshold below is read from `data/cases/eas_policy_pin.json`, not from "
        "the policy's own text, so a change to the text cannot move a word here "
        "without the digest check in Q9 catching it."
        % (live["policy"]["name"], live["policy"]["version"],
           live["policy"]["sha256"][:16]))
    add("")
    add("Its %s eligibility gates, read on this cut:"
        % word_number(len(live["gates"])))
    add("")
    add("| gate | value | threshold | read as | passed |")
    add("|---|---|---|---|---|")
    for gate in live["gates"]:
        add("| %s | %s | %g | %s | %s |"
            % (gate["gate"], f3(gate["value"], 2)
               if isinstance(gate["value"], float) else gate["value"],
               gate["threshold"], gate["read_as"],
               "yes" if gate["passed"] else "**no**"))
    add("")
    add("Then its rules. `observability_drop` is listed but not applied: %s."
        % next(rule["reason"] for rule in live["rules"]
               if not rule["applied"]))
    add("")
    add("| rule | quantity | boundary | max abs | months over | longest run | "
        "vs %d | fired |" % persist)
    add("|---|---|---|---|---|---|---|---|")
    for rule in live["rules"]:
        if not rule["applied"]:
            add("| %s | — | not applied | — | — | — | — | — |" % rule["rule"])
            continue
        run = rule["longest_run"]["length"]
        add("| %s | `%s` | %s | %s | %d | %d | %s | %s |"
            % (rule["rule"], rule["quantity"], rule["comparison"],
               f3(rule["max_abs"], 2), rule["n_over"], run,
               "≥" if run >= persist else "<",
               "**yes**" if rule["fired"] else "no"))
    add("")
    add("**%s** — %s." % (live["state"], live["reason"]))
    add("")
    add("Read over every firing rule it names instead, the same policy returns "
        "**%s**: %s. Both words are printed because only one of them is a "
        "clearance and neither is a verdict on the asset."
        % (live["over_all_five"]["state"], live["over_all_five"]["reason"]))
    add("")
    add("Scope, stated rather than assumed:")
    add("")
    add("* %s." % live["scope"])
    add("* the lookback: %s." % live["lookback_note"])
    add("")


def _md_pair(add, a):
    pair = a["case"]["pair"]
    add("## Q7 · What does each reader say?")
    add("")
    add("| reader | word | what the word is |")
    add("|---|---|---|")
    add("| trajectory rule (this repository) | **%s** | %s |"
        % (pair["trajectory"], a["case"]["trajectory"]["rule"]))
    add("| erlking-eas v%s, rules this cut carries | **%s** | %s |"
        % (a["case"]["live"]["policy"]["version"], pair["live"],
           a["case"]["live"]["reason"]))
    add("| erlking-eas v%s, every firing rule | **%s** | %s |"
        % (a["case"]["live"]["policy"]["version"], pair["live_over_all_five"],
           a["case"]["live"]["over_all_five"]["reason"]))
    add("")
    add("Two vocabularies over one set of numbers, and the record does not "
        "collapse them: `concur` is **%s** and `collapsed` is **%s**. The words "
        "are not translations of each other — the trajectory rule prints %s, "
        "each word a reading of the reference months' own envelope, while the "
        "policy chooses among its %s words %s over a fixed floor. A sentence "
        "that said 'the policy agrees' or 'the policy disagrees' would be "
        "inventing a mapping neither reader defines."
        % (pair["concur"], pair["collapsed"],
           listed(list(a["case"]["trajectory"]["vocabulary"])),
           word_number(len(a["case"]["live"]["policy"]["states"])),
           listed(list(a["case"]["live"]["policy"]["states"]))))
    add("")


def _md_why(add, a):
    decisive = a["case"]["decisive"]
    st = a["case"]["streams"]["streams"]
    traj = a["case"]["trajectory"]
    add("## Q8 · Why do they differ, and by how much?")
    add("")
    add("Persistence is **not** the difference: the trajectory rule wants %d "
        "consecutive months, the policy wants %d. The difference is entirely the "
        "boundary. For each rule, the largest floor at which that rule would still "
        "fire over %d consecutive held-out months is the smallest `|x|` inside its "
        "best run — so the margin to the pinned floor is the distance between "
        "silence and a word:"
        % (traj["persist_months"], decisive["persist_months"],
           decisive["persist_months"]))
    add("")
    add("| quantity | rule | pinned floor | largest firing floor | run | margin | "
        "fires |")
    add("|---|---|---|---|---|---|---|")
    for entry in decisive["rules"]:
        add("| `%s` | %s | %g | %s | %s | %s | %s |"
            % (entry["quantity"], entry["rule"], entry["policy_floor"],
               f3(entry["fires_at"], 3),
               "%s → %s" % (entry["months"][0], entry["months"][-1])
               if entry["months"] else "—",
               f3(entry["margin"], 3),
               "yes" if entry["fires_under_policy"] else "no"))
    add("")
    closest = decisive["closest"]
    dim_floor = next(entry["threshold"] for entry in a["case"]["live"]["rules"]
                     if entry["rule"] == "dimension_z")
    add("The closest of them is `%s`, whose run %s → %s would fire at a floor of "
        "**%s** against the pinned **%g** — a margin of **%s**. Nothing in the "
        "held-out span comes near the policy's floor under either of its block or "
        "dimension rules, while the envelope it is compared with is far tighter: "
        "`gamma2`'s q95 is %s and `A`'s is %s, both against a dimension floor of "
        "%g. The trajectory rule therefore calls runs of %d months on `gamma2` and "
        "%d months on `A` a CHANGE while the policy, correctly applying its own "
        "floors, has nothing to report."
        % (closest["quantity"], closest["months"][0], closest["months"][-1],
           f3(closest["fires_at"], 3), closest["policy_floor"],
           f3(closest["margin"], 3), f3(st["gamma2"]["q95"], 2),
           f3(st["A"]["q95"], 2), dim_floor,
           st["gamma2"]["longest_run"]["length"],
           st["A"]["longest_run"]["length"]))
    add("")
    add("Read the other way: the two readers are not in conflict, they are "
        "answering at two distances. `gamma2`'s sustained %d-month departure is "
        "*large* against the reference months' own spread and *small* against a "
        "fixed absolute-scale floor — |z| reaches only %s there, and the "
        "envelope's own q95 sits at %s. Naming that gap is the whole finding; a "
        "single word would have hidden it."
        % (st["gamma2"]["longest_run"]["length"],
           f3(st["gamma2"]["max_abs"], 2), f3(st["gamma2"]["q95"], 2)))
    add("")


def _guard_evidence(name, cell):
    if name == "window_event_free":
        return ("%d calendar months and %d days before the event month"
                % (cell["month_gap"], cell["days_before_event"]))
    if name == "buckets_disjoint_and_reconciled":
        return ("%d chips in, %d reconciled across %d buckets"
                % (cell["total"], cell["kept"], len(cell["reconciliation"]) - 1))
    if name == "no_pooling":
        return cell["evidence"]
    if name == "envelope_from_reference_months_only":
        return "%d reference months, %s" % (cell["n_reference_months"],
                                            cell["aggregation"])
    if name == "every_held_out_bin_fielded":
        return "%d of %d bins carry all %s streams" % (
            cell["n_held_out"] - len(cell["unfielded"]), cell["n_held_out"],
            word_number(len(STREAMS)))
    if name == "state_recomputes_from_stored_numbers":
        return "trajectory %s → %s, live %s → %s" % (
            cell["trajectory"]["stored"], cell["trajectory"]["recomputed"],
            cell["live"]["stored"], cell["live"]["recomputed"])
    if name == "record_is_max_over_streams":
        return "longest run %d, on `%s`" % (cell["longest_run"], cell["stream"])
    if name == "control_after_case_no_influence":
        return "%s, identical: %s" % (" → ".join(cell["order"]),
                                      cell["difference"] is None)
    if name == "runs_reported_separately_from_magnitude":
        return cell["evidence"]
    if name == "magnitudes_at_1e_12":
        return "%s over %d chips, max difference %g" % (
            cell["identity"], cell["n_chips_checked"], cell["max_abs_difference"])
    if name == "vocabulary_scan":
        stray = cell["stray"]["trajectory"] + cell["stray"]["live"]
        return "%s: emitted %s (trajectory) and %s (live); %s" % (
            cell["evidence"],
            "/".join(cell["emitted"]["trajectory"]),
            "/".join(cell["emitted"]["live"]),
            ("stray " + ", ".join(stray)) if stray else "no stray word")
    if name == "no_failure_language":
        return "scanned %s" % cell["scanned"]
    if name == "thresholds_pinned":
        return "%d digests verified against %s" % (cell["n_verified"],
                                                   cell["pin"])
    if name == "event_is_the_registrys":
        return cell["event"]
    if name == "policy_sha_matches_pin":
        return "%d of %d digests match" % (
            len([d for d in cell["digests"] if d["matched"]]), len(cell["digests"]))
    if name == "hindsight":
        perturbations = [key for key, value in cell.items()
                         if isinstance(value, dict) and "identical" in value]
        return "%s perturbations, %s" % (
            word_number(len(perturbations)),
            "all identical"
            if all(cell[key]["identical"] for key in perturbations)
            else "not all identical")
    return cell.get("evidence", "—")


def _md_guards(add, a):
    guards = a["guards"]
    hindsight = guards["hindsight"]
    add("## Q9 · Does the finding survive the guards?")
    add("")
    add("Every guard below is computed from the same artefact the answers above "
        "are read from, and every one of them passes — `guards.passed` is **%s**."
        % guards["passed"])
    add("")
    add("| guard | passed | evidence |")
    add("|---|---|---|")
    for name in GUARD_ORDER:
        cell = guards[name]
        add("| %s | %s | %s |"
            % (name, "yes" if cell.get("passed") else "**no**",
               _guard_evidence(name, cell)))
    add("")
    add("The hindsight guard deletes evidence and checks that nothing moves. It "
        "holds out the reference numbers and the held-out numbers separately, "
        "because a tail beyond the last held-out bin can only move one of them:")
    add("")
    add("| perturbation | chips dropped | reference moved | held-out moved | "
        "bins vanished |")
    add("|---|---|---|---|---|")
    for name in ("tail_deleted", "tail_absurd", "reference_from_window_alone"):
        cell = hindsight[name]
        add("| %s | %s | %s | %s | %s |"
            % (name,
               cell.get("n_chips_dropped", cell.get("n_chips")),
               "no" if cell.get("reference_difference") is None else "**yes**",
               "no" if cell.get("held_out_difference") is None else "**yes**",
               ", ".join(cell.get("vanished_bins") or []) or "none"))
    add("")
    moved = hindsight["tail_deleted"]["labels_moved"]
    if moved:
        add("%s %s when the post-event chips are deleted — %s. That is the split's "
            "`(pre)` marker, owned by the split rather than by a number, so it is "
            "reported and excluded from the comparison instead of being allowed "
            "to pass as a difference."
            % (sentence_number(len(moved)),
               "label moves" if len(moved) == 1 else "labels move",
               ", ".join("%s: `%s` → `%s`"
                         % (month, cell["before"], cell["after"])
                         for month, cell in sorted(moved.items()))))
        add("")
    add("And the guard that matters most for this file being *prospective*: %d "
        "held-out chips, the latest on **%s**; the event is **%s**. Nothing after "
        "it is read as evidence."
        % (hindsight["held_out_is_pre_event"]["n_held_out_chips"],
           hindsight["held_out_is_pre_event"]["latest_held_out"],
           hindsight["held_out_is_pre_event"]["event"]))
    add("")


def _md_control(add, a):
    cmp_ = a["comparison"]
    control = a["control"]
    part = control["partition"]["buckets"]
    add("## Q10 · Does the control behave the same way?")
    add("")
    add("The control is **%s**, event **%s**, read *after* the case and never as "
        "an input to it. The two windows arrive from different places:"
        % (control["label"], control["event"]))
    add("")
    add("* **control** — %s." % cmp_["window_source"]["carola"])
    add("* **case** — %s." % cmp_["window_source"]["morandi"])
    add("")
    add("That asymmetry is the point: a rule-selected window and a "
        "registry-declared one are both legitimate, and the file says which is "
        "which instead of treating them as interchangeable.")
    add("")
    add("Two scopes are named in this table, so both are labelled. The reference "
        "*state* is the plate's: every chip the table carries for it. The "
        "*frozen reference* is the slice this case builds its envelope from — "
        "the months inside each asset's window — and it is the smaller of the "
        "two.")
    add("")
    add("| | %s (case) | %s (control) |" % (a["case"]["label"], control["label"]))
    add("|---|---|---|")
    for label, key in (("reference-state months", "months"),
                       ("reference-state chips", "n_chips")):
        add("| %s | %d | %d |" % (label, a["case"]["reference"][key],
                                  control["reference"][key]))
    add("| frozen reference months | %d | %d |"
        % (a["case"]["envelope"]["n_reference_months"],
           control["envelope"]["n_reference_months"]))
    add("| frozen reference chips | %d | %d |"
        % (a["case"]["partition"]["buckets"]["reference"]["n_chips"],
           part["reference"]["n_chips"]))
    add("| held-out months | %d | %d |"
        % (a["case"]["partition"]["buckets"]["outside"]["months"],
           part["outside"]["months"]))
    add("| held-out chips | %d | %d |"
        % (a["case"]["partition"]["buckets"]["outside"]["n_chips"],
           part["outside"]["n_chips"]))
    add("| trajectory rule | **%s** — %s | **%s** — %s |"
        % (a["case"]["trajectory"]["state"], a["case"]["trajectory"]["plate_verdict"],
           control["trajectory"]["state"], control["trajectory"]["plate_verdict"]))
    add("| pinned policy | **%s** | **%s** |"
        % (a["case"]["pair"]["live"], control["pair"]["live"]))
    add("| first month outside | %s | %s |"
        % (day_of(a["case"]["leads"]["first_flagged"]),
           day_of(control["leads"]["first_flagged"])))
    add("| first CHANGE | %s | %s |"
        % (day_of(a["case"]["leads"]["first_change"]),
           day_of(control["leads"]["first_change"])))
    add("")
    add("Per stream, longest run — the comparison is of *runs*, not of words:")
    add("")
    add("| stream | case run | case months outside | control run | control months "
        "outside | differs |")
    add("|---|---|---|---|---|---|")
    for row in cmp_["rows"]:
        add("| `%s` | %d | %d / %d | %d | %d / %d | %s |"
            % (row["stream"], row["case"]["run"], row["case"]["n_fired"],
               row["case"]["n_months"], row["control"]["run"],
               row["control"]["n_fired"], row["control"]["n_months"],
               "yes" if row["stream"] in cmp_["differs"] else "no"))
    add("")
    add("The control does not reproduce the case: its longest run in any stream "
        "is %d against the case's %d, its own word is **%s** against **%s**, and "
        "its first month outside the envelope is %s, which is %s days before its "
        "event — a signed lead with no sustained departure behind it. That is what "
        "a control is for: the same machinery, run on a second asset, does not "
        "print the same word, so the case's **%s** is a reading of those numbers "
        "rather than something the code prints everywhere."
        % (max(control["streams"]["runs"].values()),
           max(a["case"]["streams"]["runs"].values()),
           control["trajectory"]["state"], a["case"]["trajectory"]["state"],
           label_of(control["leads"]["first_flagged"]),
           control["leads"]["signed_lead_days"],
           a["case"]["trajectory"]["state"]))
    add("")
    add(cmp_["note"])
    add("")


def _md_selection(add, a):
    w = a["window"]
    starts = w["starts_tried"]
    add("## Q11 · Was the window selection genuinely event-free?")
    add("")
    add("Three separate claims, each checkable without leaving this file:")
    add("")
    add("1. **The rule reads no measurement.** It reads %s. It never reads %s — "
        "the resolver's argument list contains no chip measurement and no event "
        "date, so there is no path by which a held-out value could have moved the "
        "window."
        % (" and ".join(w["rule"]["reads"]), " or ".join(w["rule"]["never_reads"])))
    add("2. **The rule was applied from the first month, not searched.** It starts "
        "at the first occupied month (%s) and keeps the first start that "
        "satisfies it. Only %d starts were tried, and the first was rejected for a "
        "reason that has nothing to do with any measurement: %s held %d chip%s."
        % (w["first_occupied"], w["n_starts_tried"], starts[0]["shorts"][0]["month"],
           starts[0]["shorts"][0]["chips"],
           "" if starts[0]["shorts"][0]["chips"] == 1 else "s"))
    add("3. **The guard is counted in months, not days.** The event month is %s; "
        "the window's last month is %s, %d bins earlier and %d days before the "
        "event date. With a %d-month guard the check passes; a 365-day guard "
        "would have refused this window by %d days even though %d whole "
        "bins separate it from the event."
        % (w["event_month"], w["window"][1], w["month_gap"],
           w["days_before_event"], w["guard_months"],
           365 - w["days_before_event"], w["month_gap"]))
    add("")
    add("What this does *not* claim: that the earliest satisfying window is the "
        "only defensible one. A later window would satisfy the same rule and "
        "would not be selected by it — the rule fixes *one* window in advance, "
        "and the record's job is to be honest that it is one, not to pretend it "
        "is the only one. The held-out span (%d chips, %d bins) is not thinner "
        "than the reference (%d chips, %d bins), so the test set is not a "
        "convenience sample from a sparse tail."
        % (a["case"]["partition"]["buckets"]["outside"]["n_chips"],
           a["case"]["partition"]["buckets"]["outside"]["months"], w["n_chips"],
           w["n_months"]))
    add("")


def _md_counterfactual(add, a):
    decisive = a["case"]["decisive"]
    st = a["case"]["streams"]["streams"]
    add("## Q12 · What would change the answer?")
    add("")
    add("Numbers, not adjectives. Each line names the measurement and the value it "
        "would take to move a word:")
    add("")
    add("| to move | the pinned value | it would take | which is |")
    add("|---|---|---|---|")
    for entry in decisive["rules"]:
        if entry["margin"] is None:
            continue
        add("| `%s` fires under %s | %g | a floor of ≤ %s | %s below today's |"
            % (entry["quantity"], entry["rule"], entry["policy_floor"],
               f3(entry["fires_at"], 3), f3(entry["margin"], 3)))
    add("")
    add("The tightest of these is `%s` at **%s**: the policy's dimension floor is "
        "**%g**, the strongest single held-out value reaches `|z|` = %s, and the "
        "weakest month inside the strongest two-month run is %s — so the word "
        "would flip only if the floor were lowered by %s. Persistence is not the "
        "lever: even at one month no held-out value reaches **%g**, so the policy "
        "would still be silent. `alert_enabled` is `%s` in the pin, so no rule "
        "could have produced ALERT here in any case."
        % (decisive["closest"]["quantity"], f3(decisive["closest"]["fires_at"], 3),
           decisive["closest"]["policy_floor"],
           f3(max(st[key]["max_abs"] for key in st), 3),
           f3(decisive["closest"]["fires_at"], 3),
           f3(decisive["closest"]["margin"], 3),
           decisive["closest"]["policy_floor"],
           a["case"]["live"]["policy"]["alert_enabled"]))
    add("")
    block_run = max(a["case"]["streams"]["runs"][key] for key in ("R_B", "R_G"))
    gamma = st["gamma2"]
    add("On the other side, the trajectory word moves with the persistence rule, "
        "not with the envelope: the runs are %s, so any persistence up to **%d** "
        "keeps CHANGE — a block run and a dimension run both clear it — and at "
        "**%d** the longest block run is exhausted and the word falls to WATCH on "
        "dimensions alone. What would *not* move it is a slightly looser envelope: "
        "`gamma2` sits outside q05/q95 in %d of %d held-out months, so widening "
        "the reference spread by a little leaves the same %d-month run."
        % (", ".join("`%s` %d" % (key, a["case"]["streams"]["runs"][key])
                     for key in a["case"]["streams"]["order"]),
           block_run, block_run + 1, gamma["n_fired"], gamma["n_months"],
           gamma["longest_run"]["length"]))
    add("")


def _md_limits(add, a):
    live = a["case"]["live"]
    add("## Q13 · What does this record not answer?")
    add("")
    missing = next(rule for rule in live["rules"] if not rule["applied"])
    absent = sorted({rule["rule"] for rule in live["rules"] if not rule["applied"]})
    add("* **The observability level.** It carries chips, not `(month, level)` "
        "records, so `%s` (%g) cannot be read and %s of the policy's %s firing "
        "rules %s absent. The policy's word restricted to the rules this cut "
        "carries is `%s`; over every firing rule it names it is `%s`. Neither is "
        "a verdict on the asset."
        % (missing["rule"], missing["threshold"], word_number(len(absent)),
           word_number(len(live["firing_rules"])),
           "is" if len(absent) == 1 else "are",
           live["state"], live["over_all_five"]["state"]))
    unreachable = [state for state in live["policy"]["states"]
                   if state not in live["vocabulary"]]
    add("* **%s of the policy's %s words.** %s %s unreachable in an "
        "observability cut that carries no observability level and no second "
        "stream, so a reader should not read their absence as clearance. This "
        "cut prints only %s."
        % (sentence_number(len(unreachable)),
           word_number(len(live["policy"]["states"])),
           listed(unreachable),
           "is" if len(unreachable) == 1 else "are",
           listed(live["vocabulary"])))
    if "ALERT" in unreachable:
        add("* **ALERT, in particular.** `alert_enabled` is `%s` in the pinned "
            "policy, so that word could not have been produced here whatever the "
            "numbers: the policy requires a human to escalate."
            % live["policy"]["alert_enabled"])
    add("* **Causation.** Neither reader names a cause, and no number here is "
        "evidence that anything caused the collapse. The record says a monthly "
        "state departed from a frozen reference; it does not say why.")
    add("* **Sub-monthly timing.** Every state is a calendar-month median, so a "
        "lead counted in months is a count of bins, not of days: a departure "
        "that began on the 25th of a month reads as that whole month. The lead "
        "figures in Q5 are anchored to the bin's tick for exactly this reason.")
    ref_months = a["case"]["envelope"]["n_reference_months"]
    add("* **A population.** One case and one control is not a sample. The "
        "quantiles are empirical quantiles of %d reference months, so each is "
        "positioned by one or two observations; the envelope is a description of "
        "those months, not an estimated population distribution." % ref_months)
    post = a["case"]["partition"]["buckets"]["post"]["n_chips"]
    add("* **The %s post-event chips.** They are kept in the cut and named in "
        "the partition, and they are read by nothing except the hindsight guard, "
        "which deletes them and shows that the answer does not move."
        % word_number(post))
    add("")


def _md_close(add, a):
    pair = a["case"]["pair"]
    decisive = a["case"]["decisive"]
    st = a["case"]["streams"]["streams"]
    named = ("gamma2", "A", "F")
    runs = [st[key]["longest_run"]["length"] for key in named]
    add("## The pair, in one paragraph")
    add("")
    add("With the reference frozen at %s → %s (%d chips, selected by a rule that "
        "cannot see the event) and %s bins held out before the event, the "
        "trajectory rule reads **%s**: %s leave the reference envelope for %d, "
        "%d and %d held-out months respectively, and the block runs rise with "
        "them, from %s. EAS policy v%s reads **%s**, because no held-out value "
        "reaches its block or dimension floor — the closest run would need a floor "
        "of %s against the pinned %g. The distance between the two words is that "
        "margin, and it is the finding: a departure that is unambiguous against "
        "the reference months' own spread is invisible against a fixed "
        "absolute-scale floor of the same persistence rule. Both words are "
        "printed; neither is retracted and neither is preferred."
        % (a["window"]["window"][0], a["window"]["window"][1],
           a["window"]["n_chips"],
           word_number(a["case"]["partition"]["buckets"]["outside"]["months"]),
           pair["trajectory"], listed(named),
           runs[0], runs[1], runs[2],
           a["case"]["leads"]["first_change"]["label"],
           a["case"]["live"]["policy"]["version"], pair["live"],
           f3(decisive["closest"]["fires_at"], 3),
           decisive["closest"]["policy_floor"]))
    add("")
    add("## Reproducing this file")
    add("")
    add("```")
    add("python3 code/cases/ocv_prospective.py            # write the "
        "artefact, the report, the plate")
    add("python3 code/cases/ocv_prospective.py --verify   # rebuild and compare "
        "with the committed files")
    add("python3 code/cases/ocv_prospective.py --selftest # the pinned values, "
        "the guards, the vocabulary; writes nothing")
    add("```")
    add("")
    add("The artefact beside this page is `%s`; the comparison is `%s`. Every "
        "number in every answer above is read from that JSON, and `--verify` "
        "rebuilds it from the chips and refuses to agree with a stale file."
        % (a["outputs"]["json"], a["outputs"]["cmp_json"]))
    add("")
    add("`schema` **%s** · pin `%s` sha256 `%s` (%d digests verified) · "
        "policy digests %s · guards %s"
        % (a["schema"], a["pin"]["file"], a["pin"]["sha256"],
           a["pin"]["n_verified"],
           ", ".join(digest["path"].rsplit("/", 1)[-1]
                     for digest in a["pin"]["digests"]),
           "all pass" if a["guards"]["passed"] else "**a guard failed**"))
    add("")


def comparison_md(artefact):
    """The case and the control beside each other, with the window sources named."""
    cmp_ = artefact["comparison"]
    case, control = cmp_["case"], cmp_["control"]
    out = ["# Prospective read — the case and its control, side by side", ""]
    out.append("| | %s (case) | %s (control) |" % (case["site"], control["site"]))
    out.append("|---|---|---|")
    out.append("| event | %s | %s |" % (case["event"], control["event"]))
    out.append("| window source | %s | %s |" % (cmp_["window_source"]["morandi"],
                                                cmp_["window_source"]["carola"]))
    out.append("| trajectory rule | **%s** — %s | **%s** — %s |"
               % (case["trajectory"], case["plate_verdict"],
                  control["trajectory"], control["plate_verdict"]))
    out.append("| pinned policy | **%s** | **%s** |"
               % (case["pair"]["live"], control["pair"]["live"]))
    out.append("| policy, every firing rule | **%s** | **%s** |"
               % (case["pair"]["live_over_all_five"],
                  control["pair"]["live_over_all_five"]))
    for label, key in (("first month outside", "first_flagged"),
                       ("first persistence", "first_persistence"),
                       ("first CHANGE", "first_change")):
        out.append("| %s | %s | %s |"
                   % (label, day_of(case["leads"][key]),
                      day_of(control["leads"][key])))
    out.append("| signed lead (days) | %s | %s |"
               % (case["leads"]["signed_lead_days"],
                  control["leads"]["signed_lead_days"]))
    out.append("| sustained lead (days) | %s | %s |"
               % (case["leads"]["sustained_lead_days"]
                  if case["leads"]["sustained_lead_days"] is not None else "—",
                  control["leads"]["sustained_lead_days"]
                  if control["leads"]["sustained_lead_days"] is not None else "—"))
    out.append("")
    if case["pair"]["concur"] and control["pair"]["concur"]:
        verdict = ("The readers agree on both, and that agreement is printed as "
                   "two words rather than folded into one")
    else:
        which = ("both" if not case["pair"]["concur"]
                 and not control["pair"]["concur"] else "one")
        verdict = ("The readers disagree on %s, and the disagreement is printed "
                   "as two words rather than resolved into one" % which)
    out.append("Both readers, both assets, the same layer. %s — no word on this "
               "page is a reconciliation of the two readers, here or anywhere."
               % verdict)
    out.append("")
    out.append("| stream | case run | case months outside | control run | "
               "control months outside | case median | control median |")
    out.append("|---|---|---|---|---|---|---|")
    for row in cmp_["rows"]:
        out.append("| `%s` | %d | %d / %d | %d | %d / %d | %s | %s |"
                   % (row["stream"], row["case"]["run"], row["case"]["n_fired"],
                      row["case"]["n_months"], row["control"]["run"],
                      row["control"]["n_fired"], row["control"]["n_months"],
                      f3(row["case"]["median"], 2), f3(row["control"]["median"], 2)))
    out.append("")
    out.append("Streams whose longest run differs: %s."
               % (", ".join("`%s`" % key for key in cmp_["differs"]) or "none"))
    out.append("")
    out.append(cmp_["note"])
    out.append("")
    return "\n".join(out)


def draw(artefact, path):
    """The plate: held-out ``|z|`` against the pinned floors, and the runs.

    Left, every stream's held-out ``|z|`` month by month with the policy's block
    and dimension floors drawn across it -- the picture of "how far from firing".
    Right, the longest run per stream against the persistence rule both readers
    use. Neither panel averages anything, and neither draws the reference months:
    the question is what the held-out span looks like against a frozen boundary.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    case = artefact["case"]
    held = case["audit"]["held_out"]
    order = case["streams"]["order"]
    runs = case["streams"]["runs"]
    floors = {rule["rule"]: rule["threshold"] for rule in case["live"]["rules"]
              if rule["applied"]}
    block_floor = min(floors.values())
    dim_floor = max(floors.values())
    persist = case["trajectory"]["persist_months"]

    plt.rcParams.update({"font.size": 8, "axes.titlesize": 9,
                         "axes.spines.top": False, "axes.spines.right": False})
    fig, (left, right) = plt.subplots(
        1, 2, figsize=(12.4, 4.4), gridspec_kw={"width_ratios": [1.55, 1]})
    x = list(range(len(held)))
    colours = {}
    for index, key in enumerate(order):
        values = [abs(row["values"][key]) if row["values"][key] is not None
                  else float("nan") for row in held]
        line, = left.plot(x, values, marker="o", markersize=2.6, linewidth=1.1,
                          label=key)
        colours[key] = line.get_color()
    left.axhline(block_floor, color="0.25", linestyle="--", linewidth=0.9)
    left.axhline(dim_floor, color="0.25", linestyle=":", linewidth=0.9)
    left.text(len(held) - 0.4, block_floor, " block floor |z| > %g" % block_floor,
              va="bottom", ha="right", fontsize=7, color="0.25")
    left.text(len(held) - 0.4, dim_floor, " dimension floor |z| > %g" % dim_floor,
              va="bottom", ha="right", fontsize=7, color="0.25")
    left.set_title("held-out |z| per month, against the pinned policy's floors")
    left.set_xticks(x)
    left.set_xticklabels([row["label"] for row in held], rotation=45,
                         ha="right", fontsize=7)
    left.set_ylabel("|z| against the frozen reference months")
    left.grid(axis="y", color="0.9", linewidth=0.6)
    left.legend(ncol=2, fontsize=7, frameon=False, loc="upper left")
    left.margins(x=0.02)

    y = list(range(len(order)))[::-1]
    longest = [runs[key] for key in order]
    right.barh(y, longest, color=[colours[key] for key in order], height=0.62)
    right.axvline(persist, color="0.25", linestyle="--", linewidth=0.9)
    right.text(persist, len(order) - 0.4, " persistence %d" % persist,
               va="bottom", ha="left", fontsize=7, color="0.25")
    for position, key, value in zip(y, order, longest):
        right.text(value + 0.12, position, "%d" % value, va="center", fontsize=7)
    right.set_yticks(y)
    right.set_yticklabels(order, fontsize=8)
    right.set_xlabel("longest run of consecutive months outside the envelope")
    right.set_title("longest run per stream (trajectory rule: ≥ %d)" % persist)
    right.set_xlim(0, max(longest + [persist]) + 1.4)
    right.grid(axis="x", color="0.9", linewidth=0.6)

    pair = case["pair"]
    fig.suptitle("%s — event %s · frozen reference %s → %s (%d chips) · held out "
                 "%d bins · trajectory %s · policy %s"
                 % (case["label"], case["event"], artefact["window"]["window"][0],
                    artefact["window"]["window"][1], artefact["window"]["n_chips"],
                    len(held), pair["trajectory"], pair["live"]),
                 fontsize=9)
    fig.savefig(path, dpi=plate.FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)


REPRODUCE = """#!/usr/bin/env bash
# Runs code/cases/ocv_prospective.py from the repository root, so its defaults
# (data/cases/eas_policy_pin.json in; data/%(site)s/%(json)s,
# figures/%(site)s/%(md)s and .jpg, data/cases/%(cmp_json)s and
# figures/cases/%(cmp_md)s out) resolve correctly. Any argument is passed
# through, e.g. --verify or --selftest.
#
# Interpreter: `python3` from PATH, unless PYTHON names another one -- needed
# where the default python3 is a virtualenv without numpy/matplotlib, e.g.
#   PYTHON=/usr/bin/python3 bash figures/%(site)s/%(md)s.sh
set -euo pipefail
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
exec "${PYTHON:-python3}" code/cases/ocv_prospective.py "$@"
"""


def publish(artefact=None, paths=None, draw_figures=True):
    """Write the artefact, both reports, the plate and the reproducer."""
    artefact = artefact or build_all()
    paths = dict(paths or out_paths())
    for key in ("json", "md", "jpg", "sh", "cmp_json", "cmp_md"):
        os.makedirs(os.path.dirname(paths[key]), exist_ok=True)
    with open(paths["json"], "w") as handle:
        json.dump(plate.round_floats(artefact), handle, indent=1, sort_keys=True)
        handle.write("\n")
    with open(paths["md"], "w") as handle:
        handle.write(prospective_md(artefact))
    with open(paths["cmp_json"], "w") as handle:
        json.dump(plate.round_floats(artefact["comparison"]), handle,
                  indent=1, sort_keys=True)
        handle.write("\n")
    with open(paths["cmp_md"], "w") as handle:
        handle.write(comparison_md(artefact))
    with open(paths["sh"], "w") as handle:
        handle.write(REPRODUCE % {
            "site": CASE,
            "json": os.path.basename(paths["json"]),
            "md": os.path.basename(paths["md"]),
            "cmp_json": os.path.basename(paths["cmp_json"]),
            "cmp_md": os.path.basename(paths["cmp_md"])})
    if draw_figures:
        draw(artefact, paths["jpg"])
    return artefact, paths


def verify(policy=None):
    """Rebuild everything and compare it with the committed files."""
    artefact = build_all(policy)
    paths = out_paths()
    problems = []
    if not os.path.exists(paths["json"]):
        problems.append("%s: missing — run the script first"
                        % os.path.relpath(paths["json"], REPO))
    else:
        with open(paths["json"]) as handle:
            committed = json.load(handle)
        difference = plate.first_difference(committed,
                                            plate.round_floats(artefact))
        if difference:
            problems.append("%s: %s" % (os.path.relpath(paths["json"], REPO),
                                        difference))
    for key, text in (("md", prospective_md(artefact)),
                      ("cmp_md", comparison_md(artefact))):
        if not os.path.exists(paths[key]):
            problems.append("%s: missing — run the script first"
                            % os.path.relpath(paths[key], REPO))
            continue
        with open(paths[key]) as handle:
            if handle.read() != text:
                problems.append("%s: the report does not match the artefact"
                                % os.path.relpath(paths[key], REPO))
    if not os.path.exists(paths["cmp_json"]):
        problems.append("%s: missing — run the script first"
                        % os.path.relpath(paths["cmp_json"], REPO))
    else:
        with open(paths["cmp_json"]) as handle:
            difference = plate.first_difference(
                json.load(handle),
                plate.round_floats(artefact["comparison"]))
        if difference:
            problems.append("%s: %s" % (os.path.relpath(paths["cmp_json"], REPO),
                                        difference))
    # The plate must be deterministic: two draws of the same artefact are
    # byte-identical, and the committed file is that draw.
    if not os.path.exists(paths["jpg"]):
        problems.append("%s: missing — run the script first"
                        % os.path.relpath(paths["jpg"], REPO))
    else:
        digests = []
        for _ in range(2):
            handle, temporary = tempfile.mkstemp(suffix=".jpg")
            os.close(handle)
            try:
                draw(artefact, temporary)
                digests.append(plate.sha256_of(temporary))
            finally:
                os.unlink(temporary)
        if len(set(digests)) != 1:
            problems.append("%s: the plate is not deterministic (two draws differ)"
                            % os.path.relpath(paths["jpg"], REPO))
        elif plate.sha256_of(paths["jpg"]) != digests[0]:
            problems.append("%s: the committed plate is not the one this "
                            "artefact draws" % os.path.relpath(paths["jpg"], REPO))
    return problems


def selftest():
    """Every promise this file makes, checked against what it builds.

    Nothing is written. The checks fall in five groups: the pin, the window rule,
    the readers' vocabularies, the guards, and the page those three produce.
    """
    checks = []

    def check(name, ok, detail=""):
        checks.append((name, bool(ok), detail))

    record = pin(strict=True)
    spec = record["families"]["policy"]
    constants = record["families"]["trajectory"]["constants"]

    check("pin: parsed, digested and every digest that can be checked matches",
          len(record["sha256"]) == 64 and record["problems"] == [],
          "%s sha256 %s…, %d of %d digests verified, problems %s"
          % (record["file"], record["sha256"][:16], record["n_verified"],
             len(record["digests"]), record["problems"]))
    check("pin: every `rules.ocv` entry is the documented value",
          (spec["rules_ocv"]["block_z"], spec["rules_ocv"]["dimension_z"],
           spec["rules_ocv"]["persist_months"],
           spec["rules_ocv"]["lookback_months"],
           spec["rules_ocv"]["observability_drop"]) == (2.5, 3.5, 2, 4, 0.5),
          str(spec["rules_ocv"]))
    check("pin: the trajectory constants are this repository's own",
          tuple(constants["quantiles"]) == (5, 50, 95)
          and constants["decision_quantile"] == 95
          and constants["persist_min_months"] == plate.PERSIST_MIN
          and list(constants["dimensions"]) == list(plate.DIMS)
          and {name: list(block) for name, block in constants["blocks"].items()}
          == {name: list(block) for name, block in plate.BLOCKS.items()}
          and sorted(constants["magnitudes"]) == sorted(MAGNITUDES),
          "persist %s, dimensions %s" % (constants["persist_min_months"],
                                         constants["dimensions"]))
    check("pin: the window rule in the pin is the rule this file resolves",
          constants["window_rule"]["name"] == WINDOW_RULE["name"]
          and constants["window_rule"]["months"] == WINDOW_RULE["months"]
          and constants["window_rule"]["min_chips"] == WINDOW_RULE["min_chips"]
          and constants["window_rule"]["never_reads"] == WINDOW_RULE["never_reads"],
          str(constants["window_rule"]))

    # -- the window rule reads acquisition dates and counts, nothing else ----
    chips, _cut = plate.load_chips(CASE)
    event = plate.SITES[CASE]["event"]
    plain = resolve_window(CASE, chips=chips, event=event)
    absurd = [dict(chip, values={dim: 1e6 for dim in plate.DIMS}) for chip in chips]
    blind = resolve_window(CASE, chips=absurd, event=event)
    check("window: replacing every measured value moves no month and no count",
          plain["window"] == blind["window"]
          and plain["n_chips"] == blind["n_chips"]
          and [month["chips"] for month in plain["months"]]
          == [month["chips"] for month in blind["months"]],
          "%s → %s, %d chips"
          % (plain["window"][0], plain["window"][1], plain["n_chips"]))
    code = _code_of(resolve_window) + _code_of(month_counts)
    check("window: neither the resolver nor the census names the event or a value",
          "event_day" not in code and "SITES" not in code
          and "values" not in _code_of(month_counts)
          and plain["event_read_during_selection"] is False,
          "scanned %d characters of code" % len(code))
    check("window: the chosen window clears the guard the rule states",
          plain["guard_passed"] and plain["month_gap"] >= GUARD_MONTHS,
          "%d calendar months before the event month" % plain["month_gap"])

    # -- the artefact --------------------------------------------------------
    a = build_all(record)
    case, control = a["case"], a["control"]
    live = case["live"]
    applied = {rule["rule"]: rule for rule in live["rules"]}
    gates = {gate["gate"]: gate for gate in live["gates"]}
    check("policy: every applied rule reads its pinned threshold",
          applied["block_z"]["threshold"] == spec["rules_ocv"]["block_z"]
          and applied["dimension_z"]["threshold"] == spec["rules_ocv"]["dimension_z"]
          and applied["observability_drop"]["threshold"]
          == spec["rules_ocv"]["observability_drop"]
          and applied["block_z"]["persist_months"]
          == spec["rules_ocv"]["persist_months"]
          and applied["block_z"]["lookback"]["months"]
          == spec["rules_ocv"]["lookback_months"],
          str({name: rule["threshold"] for name, rule in applied.items()}))
    check("policy: every gate reads its pinned threshold",
          all(gates[name]["threshold"] == value
              for name, value in spec["gates"].items()),
          str({name: cell["threshold"] for name, cell in gates.items()}))
    check("policy: the word this cut prints is one of the policy's own",
          live["state"] in spec["states"] and live["state"] in LIVE_STATES
          and case["trajectory"]["state"] in TRAJECTORY_STATES
          and case["pair"]["live"] == live["state"]
          and live["policy"]["states"] == spec["states"],
          "%s / %s" % (case["pair"]["trajectory"], case["pair"]["live"]))
    check("policy: the unreadable rule is carried as a missing observable",
          applied["observability_drop"]["applied"] is False
          and live["not_applied"] == ["observability_drop"]
          and live["over_all_five"]["state"] == "INSUFFICIENT",
          str(live["not_applied"]))
    check("readers: the pair is co-primary and never collapsed",
          a["run"]["co_primary"] and not a["run"]["collapsed"]
          and not case["pair"]["collapsed"],
          case["pair"]["word"])
    check("control: its window is the registry's, the case's is the rule's",
          control["reference"]["window"]
          == list(plate.SITES[CONTROL]["reference_window"])
          and case["reference"]["window"] == a["window"]["window"]
          and case["reference"]["choice_rule"] == WINDOW_RULE["statement"]
          and control["reference"]["choice_rule"] is None,
          "%s (registry) vs %s (rule)"
          % (control["reference"]["window"], case["reference"]["window"]))

    # -- the guards ----------------------------------------------------------
    guards_report = a["guards"]
    check("guards: GUARD_ORDER is exactly the guards this run computed",
          set(GUARD_ORDER) == set(guards_report) - {"passed"},
          "%d named, %d computed" % (len(GUARD_ORDER), len(guards_report) - 1))
    check("guards: every one of them passes",
          guards_report["passed"]
          and all(guards_report[name].get("passed") is not False
                  for name in GUARD_ORDER),
          ", ".join("%s=%s" % (name, guards_report[name].get("passed"))
                    for name in GUARD_ORDER
                    if guards_report[name].get("passed") is False) or "all pass")
    hindsight_report = guards_report["hindsight"]
    check("hindsight: deleting and faking the tail move no reference number",
          hindsight_report["tail_deleted"]["reference_difference"] is None
          and hindsight_report["tail_absurd"]["reference_difference"] is None,
          "%s / %s" % (hindsight_report["tail_deleted"]["reference_difference"],
                       hindsight_report["tail_absurd"]["reference_difference"]))
    check("hindsight: deleting and faking the tail move no held-out number",
          hindsight_report["tail_deleted"]["held_out_difference"] is None
          and hindsight_report["tail_absurd"]["held_out_difference"] is None
          and not hindsight_report["tail_deleted"]["vanished_bins"],
          "%s / %s" % (hindsight_report["tail_deleted"]["held_out_difference"],
                       hindsight_report["tail_absurd"]["held_out_difference"]))
    check("hindsight: the held-out bins stop before the event",
          hindsight_report["held_out_is_pre_event"]["passed"],
          "latest %s, event %s"
          % (hindsight_report["held_out_is_pre_event"]["latest_held_out"],
             hindsight_report["held_out_is_pre_event"]["event"]))


    # -- the page ------------------------------------------------------------
    page = prospective_md(a)
    positions = [page.find("## Q%d ·" % number) for number in range(1, 14)]
    check("page: Q1..Q13 are all present, once each, in order",
          all(position >= 0 for position in positions)
          and positions == sorted(positions)
          and all(page.count("## Q%d ·" % number) == 1
                  for number in range(1, 14)),
          "%d of 13 found in order" % len([p for p in positions if p >= 0]))
    check("page: the closing sections are present",
          "## The pair, in one paragraph" in page
          and "## Reproducing this file" in page)
    check("page: both readers' words are printed and neither is dropped",
          case["pair"]["trajectory"] in page and case["pair"]["live"] in page
          and control["pair"]["live"] in page)
    # The guard names are identifiers, not prose -- one of them is
    # `no_failure_language` -- so the scan runs on the page with the identifiers
    # taken out. What remains is what a reader reads.
    prose_page = page
    for name in GUARD_ORDER:
        prose_page = prose_page.replace("`%s`" % name, "").replace(name, "")
    check("page: no private word and no failure language",
          "NOT_EVALUABLE" not in prose_page
          and not any(word in prose_page.lower()
                      for word in ("fail", "failure", "broken")),
          "scanned %d characters" % len(prose_page))
    comparison_page = comparison_md(a)
    check("page: the comparison names both window sources",
          all(site in comparison_page for site in (CASE, CONTROL))
          and "rule selected it" in comparison_page
          and "registry declares it" in comparison_page,
          "%d rows" % len(a["comparison"]["rows"]))

    # -- the outputs this run promises --------------------------------------
    paths = out_paths()
    check("outputs: the %s files are distinct and inside the repository"
          % word_number(len(paths)),
          len(set(paths.values())) == len(paths)
          and all(os.path.abspath(path).startswith(REPO + os.sep)
                  for path in paths.values()),
          ", ".join(sorted(os.path.relpath(path, REPO)
                           for path in paths.values())))

    width = max(len(name) for name, _ok, _detail in checks)
    failed = 0
    for name, ok, detail in checks:
        if not ok:
            failed += 1
        line = "%s  %-*s" % ("PASS" if ok else "FAIL", width, name)
        if detail and not ok:
            line += "  " + detail
        print(line)
    print("\n%d checks, %d failed" % (len(checks), failed))
    return 1 if failed else 0


def report_brief(artefact):
    """One line for the console: the pair, the window and the guard verdict."""
    case = artefact["case"]
    return ("%s: trajectory %s, policy %s (%s) · window %s → %s, %d chips · "
            "guards %s" % (case["site"], case["pair"]["trajectory"],
                           case["pair"]["live"], case["live"]["state"],
                           artefact["window"]["window"][0],
                           artefact["window"]["window"][1],
                           artefact["window"]["n_chips"],
                           "all pass" if artefact["guards"]["passed"]
                           else "a guard did not pass"))


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="ocv_prospective.py",
        description="The frozen-reference prospective read of one case: the "
                    "trajectory rule and pinned EAS policy v1.0.0 over the same "
                    "numbers, a control beside them, and the guards that make "
                    "the word prospective checkable.")
    parser.add_argument("--pin", default=None,
                        help="the policy pin to read (default "
                             "data/cases/eas_policy_pin.json)")
    parser.add_argument("--no-figures", action="store_true",
                        help="numbers and reports only")
    parser.add_argument("--selftest", action="store_true",
                        help="check every promise this file makes, write nothing")
    parser.add_argument("--verify", action="store_true",
                        help="rebuild and compare with the committed files")
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args(argv)
    if args.selftest:
        return selftest()
    if args.verify:
        problems = verify(args.pin)
        if problems:
            for problem in problems:
                print("FAIL  %s" % problem)
            return 1
        if not args.quiet:
            print("OK    the artefact, both reports and the plate match the chips "
                  "(%s)" % report_brief(build_all(pin(args.pin, strict=True))))
        return 0
    artefact, written = publish(build_all(pin(args.pin, strict=True)),
                                draw_figures=not args.no_figures)
    if not args.quiet:
        for key in ("json", "md", "cmp_json", "cmp_md", "sh"):
            print("wrote %s" % os.path.relpath(written[key], REPO))
        if not args.no_figures:
            print("wrote %s" % os.path.relpath(written["jpg"], REPO))
        print("      %s" % report_brief(artefact))
    if not artefact["guards"]["passed"]:
        if not args.quiet:
            print("      a guard did not pass; the artefact records which")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())

