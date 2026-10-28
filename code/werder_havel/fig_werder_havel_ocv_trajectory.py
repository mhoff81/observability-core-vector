#!/usr/bin/env python3
"""WERDER/HAVEL - the OCV trajectory of the A 10 crossing, across a documented repair.

Figure: ``figures/werder_havel/fig_werder_havel_ocv_trajectory.jpg`` (+ ``.md``)
Result: ``data/werder_havel/fig_werder_havel_ocv_trajectory.json``

The record is the support-sized cut of the Havelbrücke Werder (A 10): **491** target
chips, one per acquisition, 2022-07-02 -> 2025-12-04, each carrying the six registered
OCV dimensions of ``ocv_reconfigure`` -- ``gamma2``, ``P`` (behaviour) and ``A``, ``D``,
``F``, ``S`` (geometry) -- computed once, from the committed mask cache, by
``werder_havel_ocv_channels_csv.py``. ``P`` is a *within-patch* persistence: the target,
reference and held-out patches are never mixed, and the export's own
``coherence_gamma2`` / ``coherence_masked_pixels`` columns are carried as provenance
only, never conflated with the mask rule's ``gamma2`` / ``A``.

The record's one documented state change is the renewal of the joints, the guardrails,
the mastic asphalt and the steel's corrosion protection, **begun 2024-07-01 and finished
2025-09-05** (Autobahn GmbH, via the Step-0 hierarchy document). The two dates are not
written into this file: they are read from the committed AOI reference
(``data/werder_havel/reference/werder_havel_aoi.json``, field ``event``, which also
carries ``from_map: false`` and the note that no version of the OSM way names the
renewal), cross-checked against the external event record when it is present, and
cross-checked against the export's own strata, which are cut on exactly those two dates
(``W2_bestand`` ends 2024-06-29, ``W3_sanierung`` runs 2024-07-03 -> 2025-09-04,
``W4_post`` starts 2025-09-05).

This plate is the **continuous-case** design of ``code/cases/ocv_trajectory.py``
(Carola's), reused rather than re-authored: the same twelve panels (a1-a4, b1-b4, c1-c4),
the same axes logic, the same standardisation and the same decision rule. Two things
differ, and only two:

  * the anchor is an **interval**, not a point: the documented works are drawn as a band
    2024-07-01 -> 2025-09-05 with a dashed edge at each end, instead of a single line;
  * the classes are **four**, because a repair has two sides: ``reference`` (the frozen
    pre-intervention annual cycle), ``outside`` (the held-out pre-intervention year),
    ``during`` (the intervention stratum) and ``post`` (the post stratum). Carola's
    third class is R2's held-out one; here the intervention stratum *is* a class,
    because the question "did the vector move while the works ran?" needs one.

The frozen reference is the record's **first full annual cycle, 2022-07 -> 2023-06**
(inclusive calendar months, 119 chips). That is a stated choice, made before any
statistic was read, and the plate draws what it leaves out: the other 122 pre-intervention
chips (2023-07 -> 2024-06) are a class of their own, drawn in the reference frame, and
standardising nothing. The choice is not tuned afterwards: the reference is the earliest
cycle the record holds, so it cannot be moved towards the result it is read against.

Rules the plate enforces (the template's, verbatim):

  * **R1 -- one cut.** Only the support cut, only its ``target`` role (491 of 991 rows).
  * **R2 -- a chosen window, and its third class drawn.** The window is inclusive
    calendar months; a reference-state chip outside it and strictly before the
    intervention is OUTSIDE, a reference-state chip on or after the start is reported
    as *unassigned* rather than absorbed.
  * **R3 -- per-block standardisation.** ``{gamma2, P}`` against the echo-bearing
    reference chips, ``{A, D, F, S}`` against every reference chip (identical here,
    since all 491 chips carry an echo: ``A >= 3`` everywhere), pooled over the orbits.
  * **R4 -- the transforms.** ``log`` for the five strictly positive dimensions, identity
    for the non-negative shift ``S``; ``ocv_reconfigure.TRANSFORM``'s own rule.
  * **R5 -- the interval.** 2000 circular moving-block draws, block 3, seed 7; monthly
    bins; the decision distribution is the distribution of the **reference months'**
    monthly medians, never of single reference chips.

What the plate answers, in the nine registered questions of ``QUESTIONS`` below: whether
the vector departs while the works run, in which block, along which channel, whether it
had already departed before the start, whether it persists, whether it returns after the
documented end, whether the post state is distinguishable from the held-out
pre-intervention state, whether the pronounced months coincide with weather, and whether
both orbit directions carry it.

What it does **not** claim: no damage, no failure precursor, no collapse language, no
causal attribution to the works, no weather correction, no new significance statistic
and no directional score. An OCV departure is a statement about the observability of the
crossing, not about its safety. The pooled record mixes four pass geometries; the plate
reports the orbit-separated reading beside the pooled one and never substitutes one for
the other.

Usage:
  python3 code/werder_havel/fig_werder_havel_ocv_trajectory.py
  python3 code/werder_havel/fig_werder_havel_ocv_trajectory.py --no-figures
  python3 code/werder_havel/fig_werder_havel_ocv_trajectory.py --selftest
  python3 code/werder_havel/fig_werder_havel_ocv_trajectory.py --verify
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import math
import os
import random
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, os.pardir, os.pardir))
CASES = os.path.join(ROOT, "code", "cases")
sys.path.insert(0, HERE)
sys.path.insert(0, CASES)

import matplotlib  # noqa: E402

matplotlib.use("Agg")

import matplotlib.dates as mdates  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

import ocv_trajectory as tpl  # noqa: E402  the continuous-case template

#: The plate's numbers are the template's numbers. They are *aliased*, not copied:
#: re-declaring ``SEED = 7`` here would let the two plates drift apart silently, and
#: ``--selftest`` could not tell the difference.
median = tpl.median
percentile = tpl.percentile
scale_of = tpl.scale_of
apply_transform = tpl.apply_transform
month_key = tpl.month_key
month_text = tpl.month_text
month_date = tpl.month_date
next_month = tpl.next_month
parse_day = tpl.parse_day
rms = tpl.rms
resampled_medians = tpl.resampled_medians
bootstrap_cell = tpl.bootstrap_cell
first_difference = tpl.first_difference
sha256_of = tpl.sha256_of
round_floats = tpl.round_floats
stdev = tpl.stdev

#: The monthly-bin machinery is generic over ``{"month", "part", ...}`` entries and is
#: reused verbatim: the bin key, the x anchor, the adjacency rule, the quantile block, the
#: envelope position, the persistence runs and the joint-pattern runs all behave exactly
#: as they do on the point-event plate, so a reader can compare the two tables line by
#: line. Only the *binning* below (four classes instead of three, two boundaries) is
#: re-implemented.
bin_key = tpl.bin_key
bin_date = tpl.bin_date
follows = tpl.follows
bin_value = tpl.bin_value
quantile_block = tpl.quantile_block
envelope_position = tpl.envelope_position
persistence_runs = tpl.persistence_runs
directional_runs = tpl.directional_runs

DIMS = tpl.DIMS
BLOCKS = tpl.BLOCKS
BLOCK_OF = tpl.BLOCK_OF
DIM_LABEL = tpl.DIM_LABEL
DIM_MD = tpl.DIM_MD
DIM_LONG = tpl.DIM_LONG
GEOMETRY_PLANE = tpl.GEOMETRY_PLANE
BOOT = tpl.BOOT
BLOCK = tpl.BLOCK
SEED = tpl.SEED
Z_CLIP = tpl.Z_CLIP
QUANTILES = tpl.QUANTILES
DECISION_QUANTILE = tpl.DECISION_QUANTILE
PERSIST_MIN = tpl.PERSIST_MIN
FIGURE_DPI = tpl.FIGURE_DPI
COL_WIDTH = tpl.COL_WIDTH

#: R4, verbatim from the template and from ``ocv_reconfigure``'s rule for the sites it
#: registers: log for every strictly positive dimension, identity for ``S``.
TRANSFORM = {dim: ("identity" if dim == "S" else "log") for dim in DIMS}
CHANNEL_COLUMN = dict(tpl.CHANNEL_COLUMN)

#: The dimensions drawn as the two planes of the plate (the template's grammar).
BEHAVIOUR_PLANE = ("gamma2", "P")

SCHEMA = "werder-havel-ocv-trajectory-v1"
NAME = "fig_werder_havel_ocv_trajectory"
DATA = os.path.join(ROOT, "data", "werder_havel")
FIGDIR = os.path.join(ROOT, "figures", "werder_havel")
TABLE = os.path.join(DATA, "support", "werder_havel_ocv_channels.csv")
AOI = os.path.join(DATA, "reference", "werder_havel_aoi.json")
META = os.path.join(DATA, "support", "werder_havel_ocv_channels_meta.json")
EVENTS_SOURCE = "/home/projects/havel_werder_analysis/havel_werder_events.json"
OUT_JSON = os.path.join(DATA, NAME + ".json")
OUT_MD = os.path.join(FIGDIR, NAME + ".md")
#: The plate is a **JPEG**, like the other companion plates of this repository
#: (``fig_kdlo_ocv_series``, ``fig_lumo_ocv_orbits``, ``fig_bautzen_ocv_orbits``):
#: matplotlib picks the encoder from the extension, and ``figures/**/*.jpg`` is
#: ignored by git exactly as the former ``.png`` was, so the plate stays a
#: regenerated binary and never enters the history.
OUT_JPG = os.path.join(FIGDIR, NAME + ".jpg")

ASSET = "werder"
ASSET_ID = "havelbruecke_werder_a10"
LABEL = "Havelbrücke Werder, A 10"
SITE_NOTE = ("The crossing's three strata are the export's own (``W2_bestand`` / "
             "``W3_sanierung`` / ``W4_post``), cut on the two documented dates of the "
             "2024/2025 renewal.")

REFERENCE_STATE = "W2_bestand"
DURING_STATE = "W3_sanierung"
POST_STATE = "W4_post"
STATE_ORDER = (REFERENCE_STATE, DURING_STATE, POST_STATE)

#: R2's choice, made before any statistic of this plate was computed: the record's
#: **first full annual cycle**, 2022-07 -> 2023-06 inclusive. It is the earliest cycle
#: the record holds, so the choice cannot have been steered by the result it is read
#: against; the year it leaves out is drawn as OUTSIDE and standardises nothing.
REFERENCE_WINDOW = ("2022-07", "2023-06")

CLASSES = ("reference", "outside", "during", "post")
CLASS_LABEL = {
    "reference": "frozen pre-intervention reference",
    "outside": "held-out pre-intervention",
    "during": "intervention period",
    "post": "post-intervention",
}
CLASS_NOTE = {
    "reference": "the frozen reference: the chips the standardisation rests on (R2/R3)",
    "outside": "held out: pre-intervention chips of the reference state that the chosen "
               "window leaves out; drawn in the reference frame, standardising nothing",
    "during": "the export's intervention stratum W3_sanierung: the works interval",
    "post": "the export's post stratum W4_post",
}
#: ``ref`` is the template's label for the reference class; the report uses it for the
#: month-class column, so a reader of both plates reads the same word.
MONTH_CLASS = {"reference": "ref", "outside": "outside", "during": "during", "post": "post"}
COLOUR = {"behaviour": "#2b6cb0", "geometry": "#2f855a", "context": "#6b46c1",
          "asc": "#2b6cb0", "desc": "#dd6b20", "reference": "#4a5568",
          "outside": "#a0aec0", "during": "#b7791f", "post": "#c53030",
          "mixed": "#b7791f", "band": "#e2e8f0", "envelope": "#cbd5e0",
          "interval": "#c53030"}
MARKER = {"reference": "o", "outside": "o", "during": "D", "post": "s", "mixed": "X"}
WEATHER_LABEL = {"temperature_c": "T [C]", "wind_speed_ms": "wind [m/s]",
                 "precipitation_mm": "precipitation [mm/h]",
                 "relative_humidity_pct": "relative humidity [%]"}
WEATHER = ("temperature_c", "precipitation_mm")
ECHO_MIN = tpl.ECHO_MIN
ORBITS = tpl.ORBITS
ORBIT_SHORT = tpl.ORBIT_SHORT
#: The header packer is the template's, imported rather than copied: it measures the
#: lines on the plate's own renderer and spreads the parts over the fewest lines that
#: fit, so the header ends where the plate ends instead of the plate growing to fit it.
packed = tpl.packed
#: The nine registered questions of the repair/intervention assessment, in the order the
#: report answers them. They are a roster, not a narrative: each is answered from the
#: frozen classes and from nothing else, and each carries its own evidence block into
#: the JSON. Q1/Q5/Q6 are the only ones the derived verdict reads.
QUESTIONS = (
    ("onset", "Does the vector depart from the frozen reference while the documented "
              "works run?"),
    ("block", "Which response block departs - behaviour ($\\gamma^2$, P), geometry "
              "(A, D, F, S), or both?"),
    ("channel", "Which channel carries the departure, and in which direction does it "
                "sit in the reference's frame?"),
    ("precursor", "Had the vector already departed before the documented start?"),
    ("persistence", "Is a departure sustained across consecutive months of the "
                    "interval?"),
    ("recovery", "Does the vector move back towards the reference after the "
                 "documented end, or does it stay displaced?"),
    ("distinction", "Is the post-intervention state distinguishable from the "
                    "held-out pre-intervention year?"),
    ("environment", "Do the pronounced months coincide with environmental variation?"),
    ("orbit", "Does the pattern hold in both orbit directions, or is it one pass "
              "geometry?"),
)

#: The only words a verdict may use. The first five are the template's own
#: (``ocv_trajectory.ALLOWED_VERDICTS``), so the two plates can be read side by side;
#: the last two name the one thing this plate can say that a point-event plate cannot:
#: whether the departure outlasted the documented end. No word here names a cause, a
#: safety state or a failure.
ALLOWED_VERDICTS = tuple(tpl.ALLOWED_VERDICTS) + (
    "sustained departure through the intervention interval",
    "sustained departure that outlasts the intervention interval",
)
#: The orthogonal second word: *where* the persistence sits relative to the interval.
ALLOWED_EXTENTS = ("no sustained departure", "confined to the intervention interval",
                   "into the post-intervention period")

#: The registered cross-check of the event metadata: the export's own strata must be cut
#: on the documented dates *inside this resolution*. A chip is a day and a stratum
#: boundary is a date, but the record is sampled, not continuous -- its 490 acquisition
#: gaps are 1 d (77), 2 d (159), 3 d (170), 4 d (80), 5 d (1) and 7 d (3) -- so a boundary
#: can only be located to within one acquisition gap, and the bound is the record's own
#: longest gap rather than an aspirational cadence. The precise work is done by the
#: side-of-the-date test below (each stratum must lie wholly on its own side of both
#: dates); this bound only refuses a *hole* at the transition, i.e. chips missing for
#: longer than the record ever misses them. The observed gaps are printed, so a reader
#: sees how far from the ideal 0 d the export actually sits.
EVENT_TOLERANCE_DAYS = 7

#: The pins: the numbers this plate is allowed to have. They are constants, not choices
#: made after reading the result, and ``--selftest`` asserts every one of them on every
#: run. The structural block is fixed by the record; the numeric block is the record's
#: census, re-derived and compared, never typed from a fit.
PINS = {
    "table": "data/werder_havel/support/werder_havel_ocv_channels.csv",
    "table_rows": 991,
    "cut_rows": 491,
    "reference_window": ("2022-07", "2023-06"),
    "reference_months": 12,
    "classes": {"reference": 119, "outside": 122, "during": 189, "post": 61},
    "unassigned": 0,
    "monthly_bins": 43,
    "event": ("2024-07-01", "2025-09-05"),
    "state_windows": {REFERENCE_STATE: ("2022-07-02", "2024-06-29"),
                      DURING_STATE: ("2024-07-03", "2025-09-04"),
                      POST_STATE: ("2025-09-05", "2025-12-04")},
}


def _fail(message):
    """Abort without writing: the artifact must never outrun its evidence."""
    raise SystemExit(NAME + ": " + message)


def rel(path):
    """A repository-relative path, for every message and every report line."""
    return os.path.relpath(path, ROOT)


def short_digest(path):
    """The first 16 hex digits of a file's sha256 (the family's citation form)."""
    return sha256_of(path)[:16]


def read_json(path):
    with open(path, encoding="utf-8") as handle:
        return json.load(handle)


def load_chips():
    """Every **target** chip of the support cut, with the audit of what it refused.

    The cut is read as it was committed (pipe-delimited, like the rest of the family's
    support tables) and re-verified in two ways before a single statistic is computed:
    the file's sha256 must be the one the generator wrote into its own meta (an
    integrity check, not a warning), and every row must carry a date, a known orbit and
    one of the export's three state windows. A row that fails is *refused* and counted,
    never silently dropped.
    """
    for path in (TABLE, META, AOI):
        if not os.path.exists(path):
            _fail("missing committed input %s" % rel(path))
    meta = read_json(META)
    digest = sha256_of(TABLE)
    if meta.get("table_sha256") != digest:
        _fail("the committed table is not the one its meta describes (%s vs %s)"
              % (digest[:16], str(meta.get("table_sha256"))[:16]))
    with open(TABLE, newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="|"))
    refused, chips, roles = [], [], {}
    for row in rows:
        role = (row.get("role") or "").strip()
        roles[role] = roles.get(role, 0) + 1
        if role != "target":
            continue
        day = parse_day(row.get("acquisition_ts") or row.get("date"))
        if day is None:
            refused.append({"id": row.get("id"), "reason": "no acquisition date"})
            continue
        orbit = (row.get("orbit_direction") or "").strip().upper()
        if orbit not in ORBITS:
            refused.append({"id": row.get("id"),
                            "reason": "orbit_direction %r" % (orbit,)})
            continue
        state = (row.get("state_window") or "").strip()
        if state not in STATE_ORDER:
            refused.append({"id": row.get("id"), "reason": "state_window %r" % (state,)})
            continue
        values = {}
        for dim in DIMS:
            raw = (row.get(CHANNEL_COLUMN.get(dim, dim)) or "").strip()
            values[dim] = float(raw) if raw else None
        weather = {}
        for column in WEATHER:
            raw = (row.get(column) or "").strip()
            weather[column] = float(raw) if raw else None
        chips.append({"id": (row.get("id") or "").strip(), "day": day, "orbit": orbit,
                      "state": state, "track": (row.get("track") or "").strip(),
                      "subswath": (row.get("subswath") or "").strip(),
                      "platform": (row.get("platform") or "").strip(),
                      "month": month_text(month_key(day)),
                      "values": values, "weather": weather,
                      "echo": (values.get("A") or 0.0) >= ECHO_MIN})
    chips.sort(key=lambda chip: (chip["day"], chip["orbit"], chip["track"]))
    audit = {"rows": len(rows), "roles": dict(sorted(roles.items())),
             "kept": len(chips), "refused": refused,
             "table_sha256": digest, "table": rel(TABLE),
             "n_dates": len({chip["day"] for chip in chips})}
    return chips, audit, meta


def event_metadata(chips, meta):
    """The documented interval, from the committed metadata -- never from this file.

    Three statements of the same two dates have to agree before the plate draws: the
    committed AOI reference's ``event`` block (which is where the dates are *read* from,
    and which also records that no version of the OSM way names the renewal), the
    channel table's own meta (written by the generator), and -- when the external
    analysis directory is present -- the operator record both derive from. On top of
    that the export's strata must be cut on those dates: ``W2`` must end before the
    start, ``W3`` must open after it and close before the end, ``W4`` must open at the
    end. That last check is what makes the metadata testable against the data it
    labels, rather than against another document.
    """
    aoi = read_json(AOI)
    block = aoi.get("event") or {}
    start, end = block.get("start"), block.get("end")
    if not start or not end:
        _fail("%s carries no event block" % rel(AOI))
    start, end = dt.date.fromisoformat(start), dt.date.fromisoformat(end)
    if not start < end:
        _fail("the documented interval is not an interval: %s -> %s" % (start, end))
    if block.get("from_map"):
        _fail("the AOI event block claims to come from the map; this plate requires "
              "the external record")
    sources = [{"name": rel(AOI), "role": "committed reference (read)",
                "kind": block.get("kind"), "note": block.get("map_note"),
                "source": block.get("source")}]
    conflicts = []
    meta_event = meta.get("event") or {}
    if (meta_event.get("start"), meta_event.get("end")) != (start.isoformat(),
                                                            end.isoformat()):
        conflicts.append("the channel table's meta names %s -> %s"
                         % (meta_event.get("start"), meta_event.get("end")))
    else:
        sources.append({"name": rel(META), "role": "generator's own provenance",
                        "kind": meta_event.get("source")})
    external = None
    if os.path.exists(EVENTS_SOURCE):
        try:
            record = read_json(EVENTS_SOURCE)
        except (ValueError, OSError):
            record = []
        if not isinstance(record, list) or not record:
            # absent or empty: the cross-check degrades to a note, never to a pass
            external = {"name": EVENTS_SOURCE, "role": "operator record (unreadable)",
                        "entries": 0, "present": True}
            sources.append({"name": EVENTS_SOURCE, "role": "operator record (unreadable)",
                            "entries": 0})
        else:
            dates = {entry.get("label"): entry.get("date") for entry in record
                     if isinstance(entry, dict)}
            found = (dates.get("sanierung_beginn"), dates.get("sanierung_ende"))
            external = {"name": EVENTS_SOURCE, "role": "operator record (cross-check)",
                        "sha256": sha256_of(EVENTS_SOURCE),
                        "begin": found[0], "finish": found[1],
                        "labels": {entry.get("label"): entry.get("source")
                                   for entry in record if isinstance(entry, dict)}}
            if found != (start.isoformat(), end.isoformat()):
                conflicts.append("the operator record names %s -> %s" % (found,))
            else:
                sources.append({"name": EVENTS_SOURCE,
                                "role": "operator record (agrees)",
                                "sha256": external["sha256"]})
    if conflicts:
        _fail("the documented interval disagrees with itself: " + "; ".join(conflicts))

    windows = {}
    for state in STATE_ORDER:
        days = [chip["day"] for chip in chips if chip["state"] == state]
        if not days:
            _fail("no target chip carries the state %s" % state)
        windows[state] = {"first_date": min(days).isoformat(),
                          "last_date": max(days).isoformat(), "n": len(days)}
    problems = []
    if windows[REFERENCE_STATE]["last_date"] >= start.isoformat():
        problems.append("%s runs past the documented start" % REFERENCE_STATE)
    if windows[DURING_STATE]["first_date"] < start.isoformat():
        problems.append("%s opens before the documented start" % DURING_STATE)
    if windows[DURING_STATE]["last_date"] >= end.isoformat():
        problems.append("%s runs past the documented end" % DURING_STATE)
    if windows[POST_STATE]["first_date"] < end.isoformat():
        problems.append("%s opens before the documented end" % POST_STATE)
    pairs = ((REFERENCE_STATE, DURING_STATE), (DURING_STATE, POST_STATE))
    gaps = {}
    for left, right in pairs:
        gap = (dt.date.fromisoformat(windows[right]["first_date"])
               - dt.date.fromisoformat(windows[left]["last_date"])).days
        gaps["%s -> %s" % (left, right)] = gap
        if gap > EVENT_TOLERANCE_DAYS:
            problems.append("%s -> %s gap %d d" % (left, right, gap))
    if problems:
        _fail("the export's strata are not cut on the documented dates: "
              + "; ".join(problems))
    return {"start": start, "end": end, "span_days": (end - start).days,
            "start_text": start.isoformat(), "end_text": end.isoformat(),
            "kind": block.get("kind"), "from_map": bool(block.get("from_map")),
            "map_note": block.get("map_note"), "sources": sources,
            "external": external, "state_windows": windows,
            "tolerance_days": EVENT_TOLERANCE_DAYS, "gaps_days": gaps,
            "check": "The export's strata are cut on the documented dates: " + ", ".join(
                "%s %s -> %s" % (state, windows[state]["first_date"],
                                 windows[state]["last_date"]) for state in STATE_ORDER)
            + " (%s)" % ", ".join("%s %d d" % (pair, gap)
                                  for pair, gap in gaps.items())}


def in_reference_window(month):
    """True when a ``YYYY-MM`` month lies inside R2's frozen window, inclusive."""
    return REFERENCE_WINDOW[0] <= month <= REFERENCE_WINDOW[1]


def assign_classes(chips):
    """R2's partition: four disjoint classes, reconciling to every target chip.

    The three *state* strata are the export's, cut on the documented dates
    (``event_metadata``); the two *pre-intervention* strata are this plate's, cut on the
    frozen window. A chip is placed by its state and, for the reference state only, by
    the month its acquisition falls in -- so a chip can never be in two classes, and the
    window is a statement about months, not about chip order. ``during`` spans the whole
    of ``W3``; no month may hold chips of two states (a documented interval that cuts
    mid-month would be visible immediately, and ``--verify`` asserts it rather than
    hiding it).
    """
    for chip in chips:
        month = chip["month"]
        if chip["state"] == REFERENCE_STATE:
            chip["cls"] = "reference" if in_reference_window(month) else "outside"
        elif chip["state"] == DURING_STATE:
            chip["cls"] = "during"
        elif chip["state"] == POST_STATE:
            chip["cls"] = "post"
        else:  # unreachable: load_chips refuses anything else
            _fail("chip %s carries the unknown state %r" % (chip["id"], chip["state"]))
    members = {name: [chip for chip in chips if chip["cls"] == name] for name in CLASSES}
    problems = []
    if sum(len(group) for group in members.values()) != len(chips):
        problems.append("the classes do not reconcile to the record")
    if len({chip["id"] for chip in chips}) != len(chips):
        problems.append("a chip id repeats")
    for name in CLASSES:
        if not members[name]:
            problems.append("the class %s is empty" % name)
    if problems:
        _fail("the partition is not usable: " + "; ".join(problems))

    months = {}
    for chip in chips:
        months.setdefault(chip["month"], set()).add(chip["cls"])
    spanning = {month: sorted(names) for month, names in months.items()
                if len(names) > 1}
    placed = {}
    for chip in chips:
        placed.setdefault(chip["cls"], []).append(chip["id"])
    disjoint = (sorted(placed) == sorted(CLASSES)
                and sum(len(ids) for ids in placed.values()) == len(chips)
                and len({chip_id for ids in placed.values() for chip_id in ids})
                == len(chips))
    if not disjoint:
        _fail("a chip sits in more than one class, or in none")
    window = {name: sorted({chip["month"] for chip in members[name]})
              for name in CLASSES}
    states = {name: sorted({chip["state"] for chip in members[name]}) for name in CLASSES}
    audit = {
        "rule": "state stratum (the export's, cut on the documented dates) x, for the "
                "reference state only, the frozen window's months",
        "reference_window": list(REFERENCE_WINDOW),
        "classes": {name: {"n": len(members[name]),
                           "months": [window[name][0], window[name][-1],
                                      len(window[name])],
                           "states": states[name],
                           "orbits": {orbit: sum(1 for c in members[name]
                                                 if c["orbit"] == orbit)
                                      for orbit in ORBITS},
                           "echo": sum(1 for c in members[name] if c["echo"]),
                           "n_dates": len({c["day"] for c in members[name]})}
                    for name in CLASSES},
        "total": len(chips),
        "reconciles": True,
        "disjoint": disjoint,
        "months_spanning_classes": spanning,
        "unassigned": 0,
    }
    return members, months, audit


def reference_stats(chips):
    """R3/R4: the centre and scale of every channel, from the **frozen reference only**.

    The reference is a single pooled set. The template splits a site's reference per
    orbit only where the pass geometry *is* the contrast (Espoo); here the contrast is
    the documented intervention, so the orbits are pooled and Q9 checks them afterwards.
    A dimension's usable chips are R3's rule, not a choice: the two behaviour channels
    are standardised on the reference chips that carry an echo (on a chip without one,
    ``gamma2`` and ``P`` have no population to be drawn from), the four geometry
    channels on every reference chip. R4's transform is applied *before* the centre and
    scale are taken, so the z of a chip is ``(t(x) - centre) / scale`` in the frame the
    reference was measured in.
    """
    reference = [chip for chip in chips if chip["cls"] == "reference"]
    if not reference:
        _fail("the frozen reference is empty")
    stats = {}
    for dim in DIMS:
        rule = ("echo-bearing reference chips" if BLOCK_OF[dim] == "behaviour"
                else "all reference chips")
        usable = [chip for chip in reference
                  if chip["values"][dim] is not None
                  and (chip["echo"] or BLOCK_OF[dim] != "behaviour")]
        if not usable:
            _fail("no reference chip carries %s (%s)" % (dim, rule))
        values = [apply_transform(chip["values"][dim], TRANSFORM[dim]) for chip in usable]
        centre, scale, source = scale_of(values)
        stats[dim] = {"transform": TRANSFORM[dim], "rule": rule, "n": len(values),
                      "n_usable_reference": len(usable), "centre": centre,
                      "scale": scale, "scale_source": source,
                      "distinct": len(set(values)),
                      "z_min": (min(values) - centre) / scale,
                      "z_max": (max(values) - centre) / scale}
    return stats


def weather_reference(chips):
    """The robust frame of the record's own co-variates (identity scale, R4's rule)."""
    reference = [chip for chip in chips if chip["cls"] == "reference"]
    out = {}
    for column in WEATHER:
        values = [chip["weather"][column] for chip in reference
                  if chip["weather"].get(column) is not None]
        if not values:
            continue
        centre, scale, source = scale_of(values)
        out[column] = {"transform": "identity", "n": len(values),
                       "distinct": len(set(values)), "centre": centre,
                       "scale": scale, "scale_source": source}
    return out


def standardise(chips, stats, weather):
    """The z of every chip, and the two block magnitudes R_B and R_G.

    No chip is standardised on itself: the centre and the scale come from
    ``reference_stats``, so the held-out pre-intervention year, the intervention and the
    post strata are all read in the frozen reference's frame. Values are clipped at
    ``Z_CLIP`` (the template's rule) and every clip is counted, because a clipped chip
    is no longer a measurement of the tail -- it is a measurement of "further out than
    the frame can say".
    """
    clipped = []
    for chip in chips:
        z = {}
        for dim in DIMS:
            raw = chip["values"][dim]
            if raw is None:
                z[dim] = None
                continue
            value = (apply_transform(raw, stats[dim]["transform"])
                     - stats[dim]["centre"]) / stats[dim]["scale"]
            if value > Z_CLIP or value < -Z_CLIP:
                clipped.append({"id": chip["id"], "month": chip["month"],
                                "dim": dim, "z": value})
                value = max(-Z_CLIP, min(Z_CLIP, value))
            z[dim] = value
        chip["z"] = z
        chip["R_B"] = rms([z[dim] for dim in BLOCKS["behaviour"]])
        chip["R_G"] = rms([z[dim] for dim in BLOCKS["geometry"]])
        weather_z = {}
        for column, cell in weather.items():
            raw = chip["weather"].get(column)
            if raw is None:
                weather_z[column] = None
                continue
            weather_z[column] = max(-Z_CLIP, min(Z_CLIP,
                                                (raw - cell["centre"]) / cell["scale"]))
        chip["weather_z"] = weather_z
    return {"n_clipped": len(clipped), "clipped": clipped[:24],
            "clip_rule": "|z| <= %.1f, counted not hidden" % Z_CLIP,
            "blocks": {name: list(dims) for name, dims in BLOCKS.items()}}


def assign_parts(chips, event):
    """Split the calendar month a documented boundary falls inside, where it must be.

    The bin is a calendar month (R5), but a bin may never mix the two sides of a
    documented boundary: a monthly state that combined, say, the last four days of the
    works with the first weeks after them would be a number about nothing. So the month
    holding a boundary is cut into a ``pre`` and a ``post`` part exactly where the record
    holds chips on both sides -- the template's rule, applied at both ends of an interval
    instead of at one point. A month that holds only one side of the boundary stays a
    plain month. The one configuration that cannot be represented (both boundaries inside
    a single month, i.e. an interval shorter than a month) is refused, not approximated.
    """
    splits, parts = {}, {}
    for boundary, name in ((event["start"], "start"), (event["end"], "end")):
        month = month_text(month_key(boundary))
        days = [chip["day"] for chip in chips if chip["month"] == month]
        if not days:
            continue
        if any(day < boundary for day in days) and any(day >= boundary for day in days):
            if month in splits and splits[month] != boundary:
                _fail("both documented boundaries fall inside %s; a calendar-month bin "
                      "cannot keep them apart" % month)
            splits[month] = boundary
            parts[month] = name
    for chip in chips:
        boundary = splits.get(chip["month"])
        chip["part"] = None if boundary is None else (
            "pre" if chip["day"] < boundary else "post")
        chip["bin"] = month_text(month_key(chip["day"])) + (
            "" if chip["part"] is None else " (%s)" % chip["part"])
    return {"split_months": {month: {"boundary": boundary.isoformat(),
                                     "part": parts[month]}
                             for month, boundary in sorted(splits.items())},
            "rule": "a calendar month is one bin unless a documented boundary falls "
                    "inside it and the record holds chips on both sides"}


def monthly(chips):
    """One row per monthly bin, in date order, with its class census and medians.

    Every bin is drawn from a single class (``assign_parts`` guarantees that no bin
    combines two sides of a boundary, and ``assign_classes`` that no month combines two
    strata), so ``cls`` is well defined even in a month that ``assign_parts`` had to cut:
    there the two parts are two bins of two classes, never one mixed bin. The row
    carries the bin's census, the monthly medians the plate draws, and -- for the
    environment panel -- the monthly medians of the two co-variates.
    """
    groups = {}
    for chip in chips:
        groups.setdefault((chip["month"], chip["part"]), []).append(chip)
    order = {"whole": 0, "pre": 0, "post": 1}
    series = []
    for key, part in sorted(groups, key=lambda item: (item[0], order[item[1] or "whole"])):
        group = groups[(key, part)]
        label = key + ("" if part is None else " (%s)" % part)
        names = {chip["cls"] for chip in group}
        entry = {"month": key, "label": label, "part": part,
                 "anchor": bin_date({"month": key, "part": part}).isoformat(),
                 "cls": names.pop() if len(names) == 1 else "mixed",
                 "n": len(group),
                 "n_dates": len({chip["day"] for chip in group}),
                 "n_ascending": sum(1 for chip in group
                                    if chip["orbit"] == "ASCENDING"),
                 "n_descending": sum(1 for chip in group
                                     if chip["orbit"] == "DESCENDING"),
                 "n_echo": sum(1 for chip in group if chip["echo"]),
                 "z": {}, "weather_z": {}, "sides": {}}
        for name in CLASSES:
            entry["sides"][name] = sum(1 for chip in group if chip["cls"] == name)
        for dim in DIMS:
            values = [chip["z"][dim] for chip in group if chip["z"][dim] is not None]
            entry["z"][dim] = median(values) if values else None
        for name in ("R_B", "R_G"):
            values = [chip[name] for chip in group if chip[name] is not None]
            entry[name] = median(values) if values else None
        for column in group[0]["weather_z"]:
            values = [chip["weather_z"][column] for chip in group
                      if chip["weather_z"].get(column) is not None]
            entry["weather_z"][column] = median(values) if values else None
        series.append(entry)
    return series


def reference_envelope(chips, months):
    """R5: the reference's own distributions, at MONTHLY and at CHIP aggregation.

    The decision distribution is the distribution of the **monthly medians of the frozen
    reference's months** -- never of single reference chips, and never the reference
    maximum. That is the like-with-like rule: a held-out monthly state may only be
    compared with a monthly reference state. The chip-level block is carried for the
    plot's band and is marked descriptive: a single chip is never compared with a
    monthly state anywhere in this analysis. ``q05``/``q50``/``q95`` are empirical type-7
    quantiles of those monthly states (R5).
    """
    def series(key):
        pairs = []
        for entry in months:
            if not entry["n"] or entry["cls"] != "reference":
                continue
            value = bin_value(entry, key)
            if value is not None:
                pairs.append({"month": entry["label"], "value": value})
        return pairs

    monthly_dist = {}
    for key in list(DIMS) + ["R_B", "R_G"]:
        pairs = series(key)
        if len(pairs) < 2:
            continue
        block = quantile_block([pair["value"] for pair in pairs])
        block["values"] = pairs
        block["role"] = ("decision distribution: a magnitude has no sign, so only an "
                         "exceedance above q95 counts"
                         if key in ("R_B", "R_G") else
                         "directional reference distribution: outside [q05, q95] is a "
                         "departure from the reference's own monthly range")
        monthly_dist[key] = block

    reference_chips = [chip for chip in chips if chip["cls"] == "reference"]
    chip_dist = {}
    for key in list(DIMS) + ["R_B", "R_G"]:
        if key in DIMS:
            values = [chip["z"][key] for chip in reference_chips
                      if chip["z"][key] is not None]
        else:
            values = [chip[key] for chip in reference_chips if chip[key] is not None]
        if values:
            chip_dist[key] = quantile_block(values)
    labels = [pair["month"] for pair in monthly_dist.get("R_B", {}).get("values", [])]
    return {"aggregation": "monthly median of each calendar-month bin",
            "quantiles": list(QUANTILES),
            "decision_quantile": DECISION_QUANTILE,
            "decision_level": "monthly",
            "persistence_min_months": PERSIST_MIN,
            "n_reference_months": len(labels),
            "months": labels,
            "monthly": monthly_dist,
            "chip": chip_dist,
            "chip_note": "descriptive only: a single chip is never compared with a "
                         "monthly state in this analysis"}
COMPARE_ORDER = ("gamma2", "P", "R_B", "A", "D", "F", "S", "R_G")
COMPARE_LABEL = {"gamma2": "z(gamma2)", "P": "z(P)", "R_B": "R_B", "A": "z(A)",
                 "D": "z(D)", "F": "z(F)", "S": "z(S)", "R_G": "R_G"}


def class_window(members):
    """``[first, last]`` month of a class, as one interval."""
    months = sorted({chip["month"] for chip in members})
    return [months[0], months[-1]] if months else [None, None]


def class_block(name, members, months, envelope, stats):
    """One held-out class read against the frozen reference, channel by channel.

    For every channel and both block magnitudes the block reports: the median over the
    class's chips in the reference's own scales, the signed departure ``r`` in transformed
    units, R5's seeded circular block interval on ``median(class) - median(reference)``,
    the direction, and where the class's **monthly** states sit in the reference's
    monthly envelope. The bootstrap is the plate's own stream -- ``random.Random(SEED)``
    per class, reference resampled first and then the class -- so a class's interval does
    not depend on which other classes were computed, and reproducing this file's numbers
    needs nothing but the seed.
    """
    reference = [chip for chip in members["reference"]]
    group = members[name]
    bins = [entry for entry in months if entry["cls"] == name and entry["n"]]
    rng = random.Random(SEED)
    result = {
        "applies": bool(group), "class": name, "label": CLASS_LABEL[name],
        "note": CLASS_NOTE[name], "draws": BOOT, "block": BLOCK, "seed": SEED,
        "window": class_window(group), "n_chips": len(group),
        "n_reference": len(reference), "n_intervals": len(bins),
        "intervals": [entry["label"] for entry in bins],
        "n_dates": len({chip["day"] for chip in group}),
        "orbits": {orbit: sum(1 for chip in group if chip["orbit"] == orbit)
                   for orbit in ORBITS},
        "criterion": ("a month is a departure only when its monthly state exceeds the "
                      "monthly reference q95 of the same magnitude (R_B, R_G: a "
                      "magnitude has no sign) or leaves the monthly reference "
                      "[q05, q95] (a signed z); a run of at least %d consecutive such "
                      "months is sustained (R1)" % PERSIST_MIN),
        "dims": {}, "departure_months": {},
    }
    if not group:
        result["reason"] = "the class holds no chip"
        return result
    for key in COMPARE_ORDER:
        if key in DIMS:
            ref_values = [chip["z"][key] for chip in reference
                          if chip["z"][key] is not None]
            values = [chip["z"][key] for chip in group if chip["z"][key] is not None]
        else:
            ref_values = [chip[key] for chip in reference if chip[key] is not None]
            values = [chip[key] for chip in group if chip[key] is not None]
        if not ref_values or not values:
            continue
        entry = bootstrap_cell(ref_values, values, rng)
        one_sided = key in ("R_B", "R_G")
        comparison = envelope_position(bins, key, envelope["monthly"][key],
                                       one_sided=one_sided)
        entry.update({
            "quantity": key, "label": COMPARE_LABEL[key],
            "block": BLOCK_OF[key] if key in DIMS else
                     ("behaviour" if key == "R_B" else "geometry"),
            "n_reference": len(ref_values), "n_class": len(values),
            "value": median(values),
            "monthly_median": comparison["level"],
            "reference_comparison": comparison,
            "one_sided": one_sided,
        })
        if key in DIMS:
            cell = stats[key]
            transformed = [apply_transform(chip["values"][key], cell["transform"])
                           for chip in group if chip["values"][key] is not None]
            entry["r"] = median(transformed) - cell["centre"]
            entry["r_unit"] = "transformed units (%s)" % cell["transform"]
            entry["reference_centre"] = cell["centre"]
        result["dims"][key] = entry
        result["departure_months"][key] = comparison["outside_envelope"]
    result["persistence"] = {
        key: persistence_runs(bins, key, envelope["monthly"][key])
        for key in ("R_B", "R_G") if key in envelope["monthly"]}
    result["directional"] = directional_runs(bins)
    return result


def medians_of(group, keys=COMPARE_ORDER):
    """The median of every registered quantity over one set of chips."""
    out = {}
    for key in keys:
        if key in DIMS:
            numbers = [chip["z"][key] for chip in group if chip["z"][key] is not None]
        else:
            numbers = [chip[key] for chip in group if chip[key] is not None]
        out[key] = median(numbers) if numbers else None
    return out


#: The verdict vocabulary is the template's, word for word (the two plates are read side
#: by side), and so is its derivation: a function of the two block magnitudes'
#: persistence rules and of nothing else.
verdict_word = tpl.verdict_word


def class_comparison(members, months):
    """The eight quantities, once per class -- the report's first table.

    Each cell is the median over the chips of that class, in the reference's own scales:
    the reference column standardises to zero by construction, the held-out column is the
    departure the assessment reads, the intervention column is the same reading over the
    works interval, and the post column is the tail that outlasts it.
    """
    result = {"order": list(COMPARE_ORDER), "labels": dict(COMPARE_LABEL),
              "units": "z (reference scales); R is a per-chip RMS in the same units",
              "classes": {}}
    for name in CLASSES:
        group = members[name]
        bins = [entry for entry in months if entry["cls"] == name and entry["n"]]
        result["classes"][name] = {
            "label": CLASS_LABEL[name], "note": CLASS_NOTE[name],
            "n_chips": len(group), "n_intervals": len(bins),
            "intervals": [entry["label"] for entry in bins],
            "values": medians_of(group)}
    return result


def orbit_block(members, months):
    """Q9: the same quantitities per orbit direction, and the departure's sign agreement.

    The reference is pooled (``reference_stats``), so the pooled median of a class says
    nothing about whether both pass geometries moved. This block splits each class by
    orbit, reports the chip-level medians in the pooled frame, and -- for the three
    held-out classes -- states whether the *departure from the reference* has the same
    sign in both orbits. A pattern that appears in one orbit only is one pass geometry's
    pattern, and the report must be able to say so.
    """
    keys = ("gamma2", "P", "R_B", "R_G")
    base = medians_of(members["reference"], keys)
    out = {"quantities": list(keys),
           "reference": {orbit: {"value": medians_of(
               [chip for chip in members["reference"] if chip["orbit"] == orbit], keys),
               "n": sum(1 for chip in members["reference"] if chip["orbit"] == orbit)}
               for orbit in ORBITS},
           "classes": {}}
    for name in ("outside", "during", "post"):
        per_orbit, separation, agreement = {}, {}, {}
        for orbit in ORBITS:
            group = [chip for chip in members[name] if chip["orbit"] == orbit]
            bins = [entry for entry in months if entry["cls"] == name and entry["n"]
                    and entry["n_ascending" if orbit == "ASCENDING" else "n_descending"]]
            per_orbit[orbit] = {"value": medians_of(group, keys), "n": len(group),
                               "n_bins": len(bins),
                               "bins": [entry["label"] for entry in bins]}
        for key in keys:
            asc = per_orbit["ASCENDING"]["value"].get(key)
            desc = per_orbit["DESCENDING"]["value"].get(key)
            separation[key] = (None if asc is None or desc is None else asc - desc)
            signs = []
            for orbit in ORBITS:
                value = per_orbit[orbit]["value"].get(key)
                if value is None or base.get(key) is None:
                    signs.append(None)
                else:
                    signs.append(value - base[key])
            agreement[key] = (None if any(sign is None for sign in signs)
                              else "same sign in both orbits" if (signs[0] > 0) ==
                              (signs[1] > 0) and (signs[0] != 0 or signs[1] != 0)
                              else "opposite sign in the two orbits")
        out["classes"][name] = {"orbits": per_orbit, "asc_minus_desc": separation,
                                "departure_sign": agreement}
    out["note"] = ("a chip-level split in the pooled frame; the verdict never reads it -- "
                   "a departure that only one orbit carries is reported as one pass "
                   "geometry's departure, not as a site-wide one")
    return out


def class_contrast(members, left, right):
    """R5's interval on ``median(right) - median(left)``, channel by channel.

    Used for the one comparison the two held-out classes need with each other: the
    post-intervention stratum against the held-out pre-intervention year (Q7). The
    bootstrap is a fresh ``random.Random(SEED)`` stream, ``left`` resampled first, so the
    interval is reproducible from the seed alone and never borrows another block's draws.
    """
    rng = random.Random(SEED)
    result = {"left": left, "right": right, "draws": BOOT, "block": BLOCK,
              "seed": SEED, "quantity": "median(%s) - median(%s)" % (right, left),
              "dims": {}}
    for key in COMPARE_ORDER:
        if key in DIMS:
            left_values = [chip["z"][key] for chip in members[left]
                           if chip["z"][key] is not None]
            right_values = [chip["z"][key] for chip in members[right]
                            if chip["z"][key] is not None]
        else:
            left_values = [chip[key] for chip in members[left] if chip[key] is not None]
            right_values = [chip[key] for chip in members[right] if chip[key] is not None]
        if not left_values or not right_values:
            continue
        entry = bootstrap_cell(left_values, right_values, rng)
        entry.update({"n_left": len(left_values), "n_right": len(right_values),
                      "block": BLOCK_OF[key] if key in DIMS else
                      ("behaviour" if key == "R_B" else "geometry")})
        result["dims"][key] = entry
    result["note"] = ("a descriptive interval on a difference between two held-out "
                      "classes; no hypothesis test and no verdict input")
    return result


def environment_context(months, blocks):
    """Q8: descriptive coincidence of the pronounced months with the co-variates.

    The environment is reported *beside* the OCV numbers, never inside them: no weather
    correction, no attribution, no causal claim. A month whose R_B or R_G exceeds the
    reference's monthly q95 is still an OCV statement if it also carries an unusual
    temperature or precipitation; the block reports that they happened at the same time
    and stops there.
    """
    rows = []
    for entry in months:
        if not entry["n"] or entry["cls"] == "reference":
            continue
        row = {"month": entry["label"], "class": entry["cls"],
               "R_B": entry["R_B"], "R_G": entry["R_G"],
               "gamma2": entry["z"]["gamma2"], "P": entry["z"]["P"],
               "weather_z": dict(entry["weather_z"]), "above_q95": {}}
        for name in ("R_B", "R_G"):
            threshold = blocks.get(entry["cls"], {}).get("persistence", {}).get(
                name, {}).get("threshold")
            row["above_q95"][name] = bool(
                threshold is not None and entry[name] is not None
                and entry[name] > threshold)
        rows.append(row)
    pronounced = [row["month"] for row in rows if any(row["above_q95"].values())]
    weather = sorted({column for row in rows for column in row["weather_z"]})
    coincident = [row["month"] for row in rows
                  if any(row["above_q95"].values())
                  and any(value is not None and abs(value) > 1.0
                          for value in row["weather_z"].values())]
    return {"applies": True, "variables": weather, "months": rows,
            "pronounced_months": pronounced, "coincident_months": coincident,
            "coincidence_rule": ("a pronounced month (R_B or R_G above the class's "
                                 "reference q95) that also carries |z| > 1 in a "
                                 "co-variate"),
            "note": "The co-variates are reported beside the OCV quantities, not "
                    "inside them: no weather correction is applied and no causal "
                    "attribution is made"}


def _blocks_word(behaviour, geometry):
    """The two persistence flags, in the verdict's own vocabulary."""
    if behaviour and geometry:
        return "both blocks"
    if behaviour:
        return "the behaviour block ($\\gamma^2$, P)"
    if geometry:
        return "the geometry block (A, D, F, S)"
    return "neither block"


def _sustained(block, key):
    return bool(block.get("persistence", {}).get(key, {}).get("sustained"))


def _run_span(run):
    """How long a persistence run is, as a noun phrase, or that there is no run.

    A run of one month is printed as a single date, not as ``2024-12 → 2024-12``: the
    span is a property of the run, and a one-month run has no span.
    """
    if not run:
        return "no month"
    if run["length"] == 1:
        return "a single month (%s)" % run["start"]
    return "%d consecutive months (%s → %s)" % (run["length"], run["start"], run["end"])


def assessment(members, months, blocks, envelope, orbit, environment, cross, event):
    """The nine registered questions, each answered from the frozen classes alone.

    Nothing here is read from the plot, from the co-variates or from a class that the
    question does not name: Q1-Q6 read the intervention and post blocks, Q4 reads the
    held-out pre-intervention year, Q7 reads the one contrast between the two held-out
    classes, Q8 reads the descriptive environment block, Q9 reads the orbit split. The
    verdict is derived -- ``verdict_word`` on the intervention's two block magnitudes --
    and the extent word says where the persistence sits relative to the interval. Neither
    word names a cause, a safety state or a failure; the plate measures observability.
    """
    during = blocks["during"]
    post = blocks["post"]
    outside = blocks["outside"]
    behaviour = _sustained(during, "R_B")
    geometry = _sustained(during, "R_G")
    verdict = verdict_word(during["persistence"]["R_B"], during["persistence"]["R_G"])
    post_sustained = _sustained(post, "R_B") or _sustained(post, "R_G")
    extent = ("into the post-intervention period" if post_sustained else
              "confined to the intervention interval" if (behaviour or geometry) else
              "no sustained departure")
    if extent not in ALLOWED_EXTENTS or verdict not in ALLOWED_VERDICTS:
        _fail("a verdict word outside the registered vocabulary: %r / %r"
              % (verdict, extent))

    def above(name_block, key):
        run = name_block["persistence"].get(key, {})
        return run.get("n_above_q95", 0), run.get("n_months", 0)

    rb_above, rb_n = above(during, "R_B")
    rg_above, rg_n = above(during, "R_G")
    ob_rb, ob_n = above(outside, "R_B")
    ob_rg, ob_n2 = above(outside, "R_G")

    channels = []
    for key in COMPARE_ORDER:
        cell = during["dims"].get(key)
        if cell is None:
            continue
        channels.append({"quantity": key, "label": cell["label"],
                         "months_outside": len(cell["reference_comparison"]
                                               ["outside_envelope"]),
                         "position": cell["reference_comparison"]["position"],
                         "delta": cell["delta"], "ci": cell["ci"],
                         "excludes_zero": cell["excludes_zero"],
                         "direction": cell["direction"]})
    moved = [cell for cell in channels
             if cell["excludes_zero"] and cell["months_outside"]]
    moved.sort(key=lambda cell: (-cell["months_outside"], -abs(cell["delta"])))

    answers = []

    answers.append({
        "key": "onset", "question": dict(QUESTIONS)["onset"],
        "answer": ("'%s' -- over the %d monthly bins of the documented interval, R_B "
                  "exceeds the reference's monthly q95 in %d of them and R_G in %d. The "
                  "word is a function of those two persistence rules and of nothing else."
                  % (verdict, during["n_intervals"], rb_above, rg_above)),
        "evidence": {"verdict": verdict, "extent": extent,
                     "R_B": {"months_above_q95": rb_above, "n_months": rb_n},
                     "R_G": {"months_above_q95": rg_above, "n_months": rg_n},
                     "n_intervals": during["n_intervals"],
                     "intervals": during["intervals"]},
        "reads": ["blocks.during.persistence.R_B", "blocks.during.persistence.R_G"]})

    answers.append({
        "key": "block", "question": dict(QUESTIONS)["block"],
        "answer": ("The sustained departure sits in %s: R_G is above its monthly q95 "
                   "for %s, R_B for %s, against a threshold of %d consecutive months."
                   % (_blocks_word(behaviour, geometry),
                      _run_span(during["persistence"]["R_G"]["longest_run"]),
                      _run_span(during["persistence"]["R_B"]["longest_run"]),
                      PERSIST_MIN)),
        "evidence": {"behaviour_sustained": behaviour, "geometry_sustained": geometry,
                     "R_B_longest_run": during["persistence"]["R_B"]["longest_run"],
                     "R_G_longest_run": during["persistence"]["R_G"]["longest_run"]},
        "reads": ["blocks.during.persistence"]})

    answers.append({
        "key": "channel", "question": dict(QUESTIONS)["channel"],
        "answer": ("%d of the %d registered quantities leave the reference's monthly "
                   "envelope in at least one interval month while their R5 interval "
                   "excludes zero: %s." % (
                       len(moved), len(channels),
                       ", ".join("%s (%s, %d months)" % (cell["label"], cell["direction"],
                                                         cell["months_outside"])
                                 for cell in moved[:6]) or "none")),
        "evidence": {"channels": channels,
                     "ranked": [cell["quantity"] for cell in moved]},
        "reads": ["blocks.during.dims", "blocks.during.departure_months"]})
    ob_word = verdict_word(outside["persistence"]["R_B"], outside["persistence"]["R_G"])
    answers.append({
        "key": "precursor", "question": dict(QUESTIONS)["precursor"],
        "answer": ("In the held-out pre-intervention year, R_B exceeds the reference's "
                   "monthly q95 in %d of %d months and R_G in %d of %d; that class's own "
                   "verdict is '%s'." % (ob_rb, ob_n, ob_rg, ob_n2, ob_word)),
        "evidence": {"verdict": ob_word, "R_B": outside["persistence"]["R_B"]
                     ["longest_run"], "R_G": outside["persistence"]["R_G"]
                     ["longest_run"], "n_intervals": outside["n_intervals"],
                     "intervals": outside["intervals"]},
        "reads": ["blocks.outside.persistence"]})

    run_b = during["persistence"]["R_B"]["longest_run"]
    run_g = during["persistence"]["R_G"]["longest_run"]

    answers.append({
        "key": "persistence", "question": dict(QUESTIONS)["persistence"],
        "answer": ("R_B: %s; R_G: %s above the reference's monthly q95. Sustained means "
                   "at least %d consecutive such months (R1)."
                   % (_run_span(run_b), _run_span(run_g), PERSIST_MIN)),
        "evidence": {"R_B": during["persistence"]["R_B"],
                     "R_G": during["persistence"]["R_G"],
                     "directional_runs": during["directional"]},
        "reads": ["blocks.during.persistence", "blocks.during.directional"]})

    post_rb = post["dims"]["R_B"]["reference_comparison"]
    post_rg = post["dims"]["R_G"]["reference_comparison"]
    placement = ("sits %s for both R_B and R_G" % post_rb["position"]
                 if post_rb["position"] == post_rg["position"] else
                 "sits %s for R_B and %s for R_G"
                 % (post_rb["position"], post_rg["position"]))
    answers.append({
        "key": "recovery", "question": dict(QUESTIONS)["recovery"],
        "answer": ("After the documented end the post class's monthly median %s; counted "
                   "month by month, %d of %d post bins above the R_B q95 and %d of %d "
                   "above the R_G q95. The extent word, from the post block's own "
                   "persistence rules, is '%s'."
                   % (placement,
                      len(post["persistence"]["R_B"]["months_above_q95"]),
                      post["n_intervals"],
                      len(post["persistence"]["R_G"]["months_above_q95"]),
                      post["n_intervals"], extent)),
        "evidence": {"R_B": {"position": post_rb["position"],
                             "months_above_q95": post["persistence"]["R_B"]
                             ["months_above_q95"]},
                     "R_G": {"position": post_rg["position"],
                             "months_above_q95": post["persistence"]["R_G"]
                             ["months_above_q95"]},
                     "during_R_B": during["dims"]["R_B"]["value"],
                     "during_R_G": during["dims"]["R_G"]["value"],
                     "post_R_B": post["dims"]["R_B"]["value"],
                     "post_R_G": post["dims"]["R_G"]["value"],
                     "extent": extent, "n_intervals": post["n_intervals"]},
        "reads": ["blocks.post.dims", "blocks.post.persistence"]})
    apart = [key for key in COMPARE_ORDER
             if cross["dims"].get(key, {}).get("excludes_zero")]
    answers.append({
        "key": "distinction", "question": dict(QUESTIONS)["distinction"],
        "answer": ("Post minus held-out pre-intervention: the R5 interval excludes zero "
                   "for %d of the %d registered quantities (%s). This is a descriptive "
                   "difference between two held-out classes, not a test."
                   % (len(apart), len(cross["dims"]),
                      ", ".join(COMPARE_LABEL[key] for key in apart) or "none")),
        "evidence": {"quantities": apart, "dims": cross["dims"]},
        "reads": ["cross.post_minus_outside"]})

    answers.append({
        "key": "environment", "question": dict(QUESTIONS)["environment"],
        "answer": ("%d pronounced %s outside the reference (%s); %d of them also "
                   "carry a co-variate departure larger than one robust scale (%s). The "
                   "co-variates are reported beside the OCV quantities: no correction, "
                   "no attribution."
                   % (len(environment["pronounced_months"]),
                      "month" if len(environment["pronounced_months"]) == 1 else "months",
                      ", ".join(environment["pronounced_months"]) or "none",
                      len(environment["coincident_months"]),
                      ", ".join(environment["coincident_months"]) or "none")),
        "evidence": {"pronounced_months": environment["pronounced_months"],
                     "coincident_months": environment["coincident_months"],
                     "months": environment["months"]},
        "reads": ["environment"]})

    during_orbits = orbit["classes"]["during"]["departure_sign"]
    answers.append({
        "key": "orbit", "question": dict(QUESTIONS)["orbit"],
        "answer": ("Splitting the intervention interval by pass direction, the departure "
                   "from the reference has %s. The reference is pooled, so this split is "
                   "the check that the pattern is not one pass geometry's."
                   % "; ".join("%s: %s" % (COMPARE_LABEL[key], during_orbits[key])
                               for key in orbit["quantities"])),
        "evidence": {"during": orbit["classes"]["during"],
                     "outside": orbit["classes"]["outside"],
                     "post": orbit["classes"]["post"],
                     "reference": orbit["reference"]},
        "reads": ["orbit"]})

    return {"applies": True, "verdict": verdict, "extent": extent,
            "verdict_is": ("derived: verdict_word on blocks.during.persistence.R_B/R_G, "
                           "then the extent word from blocks.post.persistence"),
            "vocabulary": {"verdicts": list(ALLOWED_VERDICTS),
                           "extents": list(ALLOWED_EXTENTS)},
            "questions": [{"key": key, "question": text} for key, text in QUESTIONS],
            "answers": answers,
            "interval": [event["start_text"], event["end_text"]]}


def class_colour(cls, colour):
    """The colour a class is drawn in: its block's, or its own class colour."""
    return {"outside": COLOUR["outside"], "during": COLOUR["during"],
            "post": COLOUR["post"]}.get(cls, colour)


def draw(report, path):
    """The 3 x 4 plate. Every ordinate drawn here is a number of ``report``.

    The template's grammar, unchanged: the frozen reference months filled, the held-out
    pre-intervention months hollow grey circles, the intervention months hollow amber
    diamonds, the post months hollow red squares, the frozen reference window a grey
    band, and the monthly reference envelope (q05-q95 of the reference months' monthly
    medians) a light band behind each panel that carries one. Two things differ, and both
    are the interval the record documents: the anchor is a **band** (red, with a dashed
    edge at each end) rather than a single line, and a month never spans two classes --
    ``assign_parts`` splits the month the documented end falls inside -- so a segment is
    drawn only between two bins of the same class.
    """
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "mathtext.fontset": "stix",
        "font.size": 7.4,
        "axes.titlesize": 7.4,
        "axes.labelsize": 6.9,
        "xtick.labelsize": 6.2,
        "ytick.labelsize": 6.2,
        "legend.fontsize": 6.2,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.alpha": 0.22,
        "grid.linewidth": 0.4,
        "figure.dpi": 110,
        "savefig.dpi": FIGURE_DPI,
        "savefig.bbox": "tight",
        "pdf.fonttype": 42,
    })
    months = report["months"]
    xs = [dt.date.fromisoformat(month["anchor"]) for month in months]
    classes = [month["cls"] for month in months]
    envelopes = report["envelope"]["monthly"]
    start = dt.date.fromisoformat(report["event"]["start"])
    end = dt.date.fromisoformat(report["event"]["end"])
    ref_window = report["reference"]["window"]
    fig = plt.figure(figsize=(COL_WIDTH, 0.92 * COL_WIDTH))
    grid = fig.add_gridspec(3, 4, hspace=0.80, wspace=0.62,
                            left=0.068, right=0.975, top=0.875, bottom=0.175)
    where = {}
    for row in range(3):
        for col in range(4):
            where["abc"[row] + str(col + 1)] = fig.add_subplot(grid[row, col])

    def linked(index):
        """True when bin i and bin i + 1 are adjacent in time **and** of one class."""
        if classes[index] != classes[index + 1]:
            return False
        return follows(months[index], months[index + 1])

    def decorate(ax, zero=True, envelope_key=None):
        cell = envelopes.get(envelope_key) if envelope_key else None
        if cell:
            ax.axhspan(cell["q05"], cell["q95"], color=COLOUR["envelope"],
                       alpha=0.55, zorder=0, lw=0)
        if ref_window:
            ax.axvspan(dt.date.fromisoformat(ref_window[0] + "-01"),
                       dt.date.fromisoformat(ref_window[1] + "-28"),
                       color=COLOUR["band"], alpha=0.75, zorder=0, lw=0)
        ax.axvspan(start, end, color=COLOUR["interval"], alpha=0.10, zorder=0, lw=0)
        for boundary in (start, end):
            ax.axvline(boundary, color=COLOUR["interval"], lw=0.8,
                       ls=(0, (4, 2)), zorder=1)
        if zero:
            ax.axhline(0.0, color="#4a5568", lw=0.6, ls=":", zorder=1)
        ax.set_xlim(dt.date.fromisoformat(ref_window[0] + "-01") - dt.timedelta(days=25),
                    max(xs) + dt.timedelta(days=25))
        ax.xaxis.set_major_locator(mdates.AutoDateLocator(minticks=3, maxticks=6))
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m"))
        for label in ax.get_xticklabels():
            label.set_rotation(30)
            label.set_horizontalalignment("right")

    def title(ax, ident, text):
        ax.set_title("%s  ·  %s" % (ident, text), loc="left", pad=3.0, fontsize=7.2)

    def series(ax, values, ident, text, colour, ylabel, envelope_key=None):
        """One monthly series: solid inside the reference, dashed outside it."""
        decorate(ax, envelope_key=envelope_key)
        for i in range(len(xs) - 1):
            if not linked(i):
                continue
            style = "-" if classes[i] == "reference" else (0, (3, 2))
            ax.plot(xs[i:i + 2], values[i:i + 2],
                    color=class_colour(classes[i], colour), lw=1.0, ls=style,
                    zorder=2)
        for i, (x, value) in enumerate(zip(xs, values)):
            if value is None:
                continue
            ax.plot([x], [value], marker=MARKER[classes[i]], ms=2.6,
                    mfc=colour if classes[i] == "reference" else "white",
                    mec=class_colour(classes[i], colour), mew=0.7, zorder=3)
        ax.set_ylabel(ylabel)
        title(ax, ident, text)
        return ax

    def plane(ax, ident, px, py, cx, cy, text):
        """The two dimensions of a block, month by month -- no date axis here."""
        ax.axhline(0.0, color="#4a5568", lw=0.5, ls=":")
        ax.axvline(0.0, color="#4a5568", lw=0.5, ls=":")
        for i in range(len(px) - 1):
            if linked(i) and px[i] is not None and py[i] is not None \
                    and px[i + 1] is not None and py[i + 1] is not None:
                ax.plot(px[i:i + 2], py[i:i + 2], color="#a0aec0", lw=0.6, zorder=1)
        for i in range(len(px)):
            if px[i] is None or py[i] is None:
                continue
            ax.plot([px[i]], [py[i]], marker=MARKER[classes[i]], ms=2.8,
                    mfc=COLOUR["reference"] if classes[i] == "reference" else "white",
                    mec=class_colour(classes[i], COLOUR["reference"]), mew=0.7,
                    zorder=3)
        ax.set_xlabel(cx)
        ax.set_ylabel(cy)
        title(ax, ident, text)
        return ax
    def line_only(ax, values, colour):
        for i in range(len(xs) - 1):
            if not linked(i):
                continue
            if values[i] is None or values[i + 1] is None:
                continue
            ax.plot(xs[i:i + 2], values[i:i + 2], color=colour, lw=1.0, zorder=2)

    dim_series = {dim: [month["z"][dim] for month in months] for dim in DIMS}
    block_series = {name: [month[name] for month in months] for name in ("R_B", "R_G")}

    # --- row a: the behaviour block ---------------------------------------
    series(where["a1"], dim_series["gamma2"], "a1",
           "z(%s)" % DIM_LABEL["gamma2"], COLOUR["behaviour"],
           "z(%s)" % DIM_LABEL["gamma2"], "gamma2")
    series(where["a2"], dim_series["P"], "a2", "z(P)", COLOUR["behaviour"], "z(P)", "P")
    plane(where["a3"], "a3", dim_series["gamma2"], dim_series["P"],
          "z(%s)" % DIM_LABEL["gamma2"], "z(P)", "behaviour plane")
    series(where["a4"], block_series["R_B"], "a4", "R_B", COLOUR["behaviour"],
           "R_B", "R_B")

    # --- row b: the geometry block ----------------------------------------
    series(where["b1"], dim_series["A"], "b1", "z(A)", COLOUR["geometry"], "z(A)", "A")
    series(where["b2"], dim_series["F"], "b2", "z(F)", COLOUR["geometry"], "z(F)", "F")
    plane(where["b3"], "b3", dim_series["A"], dim_series["F"], "z(A)", "z(F)",
          "geometry plane")
    series(where["b4"], block_series["R_G"], "b4", "R_G", COLOUR["geometry"],
           "R_G", "R_G")

    # --- row c: the record's own context ----------------------------------
    ax = where["c1"]
    columns = list(months[0]["weather_z"]) if months else []
    decorate(ax, zero=True)
    for index, column in enumerate(columns[:2]):
        line_only(ax, [month["weather_z"].get(column) for month in months],
                  [COLOUR["context"], COLOUR["asc"]][index])
    ax.set_ylabel("z (own reference)")
    for index, column in enumerate(columns[:2]):
        ax.text(0.03, 0.96 - 0.11 * index, WEATHER_LABEL.get(column, column),
                transform=ax.transAxes, fontsize=5.4, va="top",
                color=[COLOUR["context"], COLOUR["asc"]][index],
                bbox=dict(facecolor="white", alpha=0.7, lw=0, pad=1.0))
    title(ax, "c1", "environmental context")

    ax = where["c2"]
    decorate(ax, zero=False)
    ax.bar(xs, [month["n_ascending"] for month in months], width=20,
           color=COLOUR["asc"], lw=0)
    ax.bar(xs, [month["n_descending"] for month in months], width=20,
           bottom=[month["n_ascending"] for month in months], color=COLOUR["desc"], lw=0)
    ax.set_ylabel("chips per month")
    ax.set_ylim(0, max([month["n"] for month in months] or [1]) * 1.32)
    ax.legend(handles=[Patch(facecolor=COLOUR["asc"], label="ASC"),
                       Patch(facecolor=COLOUR["desc"], label="DESC")],
              loc="upper left", frameon=False, handlelength=1.0,
              borderaxespad=0.4, labelspacing=0.35, fontsize=5.2)
    title(ax, "c2", "cadence")

    ax = where["c3"]
    decorate(ax, zero=False)
    echo = [month["n_echo"] for month in months]
    plain = [month["n"] - month["n_echo"] for month in months]
    ax.bar(xs, echo, width=20, color=COLOUR["behaviour"], lw=0)
    ax.bar(xs, plain, width=20, bottom=echo, facecolor="white",
           edgecolor=COLOUR["behaviour"], lw=0.6)
    ax.set_ylabel("chips per month")
    ax.set_ylim(0, max([month["n"] for month in months] or [1]) * 1.12)
    title(ax, "c3", "observation availability")

    ax = where["c4"]
    grid_values = [[(month["z"][dim] if month["z"][dim] is not None else float("nan"))
                    for month in months] for dim in DIMS]
    cmap = plt.get_cmap("RdBu_r").copy()
    cmap.set_bad("#f7fafc")
    image = ax.imshow(grid_values, aspect="auto", cmap=cmap, vmin=-3.0, vmax=3.0,
                      interpolation="nearest")
    ax.set_yticks(range(len(DIMS)))
    ax.set_yticklabels([DIM_LABEL[dim] for dim in DIMS])
    step = max(1, int(math.ceil(len(months) / 5.0)))
    ticks = list(range(0, len(months), step))
    ax.set_xticks(ticks)
    ax.set_xticklabels([months[i]["label"] for i in ticks], rotation=30, ha="right")
    ax.grid(False)
    for boundary in (start, end):
        column = sum(1 for i, month in enumerate(months)
                     if dt.date.fromisoformat(month["anchor"]) < boundary)
        ax.axvline(column - 0.5, color=COLOUR["interval"], lw=0.8, ls=(0, (4, 2)))
    bar = fig.colorbar(image, ax=ax, fraction=0.05, pad=0.02)
    bar.set_label("z", fontsize=6.2)
    bar.ax.tick_params(labelsize=5.8)
    #: Two lines, not one: this panel is the narrowest of the four columns and its
    #: single-line title is wider than the panel, which -- with ``bbox_inches="tight"``
    #: -- widened the whole file to fit the title. The break is where the panel ends.
    ax.set_title("c4  ·  trajectory matrix" + chr(10) + "z per dimension × month",
                 loc="left", pad=3.0, fontsize=7.2, linespacing=1.15)
    # --- the plate's own header and legend --------------------------------
    ref = report["reference"]
    classes_audit = report["classification"]["classes"]
    head = "%s — OCV trajectory on the support-sized chip" % report["label"]
    parts = [("n = %d chips (%s → %s) · works %s → %s"
              % (report["cut"]["kept"], report["cut"]["window"][0],
                 report["cut"]["window"][1], report["event"]["start"],
                 report["event"]["end"])),
             "frozen pre-intervention reference: %d chips, %s → %s"
             % (ref["n_chips"], ref["window"][0], ref["window"][1]),
             "held-out pre-intervention: %d chips, %s → %s"
             % (classes_audit["outside"]["n"], ref["outside_window"][0],
                ref["outside_window"][1]),
             "intervention %s: %d chips, %s → %s"
             % (DURING_STATE, classes_audit["during"]["n"],
                classes_audit["during"]["months"][0],
                classes_audit["during"]["months"][1]),
             "post %s: %d chips, %s → %s"
             % (POST_STATE, classes_audit["post"]["n"],
                classes_audit["post"]["months"][0],
                classes_audit["post"]["months"][1])]
    #: The header is packed to the panel block's own width, measured from the block's
    #: leftmost ink (column one's y label) to the last column's right spine, and it starts
    #: at that same left edge -- so it reads as part of the block rather than as a caption
    #: laid over it. Written as one line the header is wider than the plate, and
    #: ``bbox_inches="tight"`` would then grow the *file* to fit the header instead of the
    #: header to fit the plate.
    renderer = fig.canvas.get_renderer()
    panel_left = min(ax.get_tightbbox(renderer).x0 for ax in where.values())
    panel_left /= fig.bbox.width
    panel_right = max(ax.get_position().x1 for ax in where.values())
    fig.suptitle(packed(fig, [head] + parts, 7.8,
                        (panel_right - panel_left) * COL_WIDTH),
                 fontsize=7.8, x=panel_left, ha="left", y=0.995, va="top",
                 linespacing=1.35)
    #: The legend is the plate's only key, so it has to fit inside the column: seven
    #: long labels on one row would be wider than the plate and ``bbox_inches="tight"``
    #: would then widen the *file* to fit the key. Two rows of four, short labels.
    handles = [
        Line2D([], [], color=COLOUR["reference"], marker="o", ms=3.2, lw=0.9,
               label="frozen reference month"),
        Patch(facecolor=COLOUR["band"], alpha=0.75, label="reference window"),
        Patch(facecolor=COLOUR["envelope"], alpha=0.55,
              label="reference envelope (q05–q95)"),
        Line2D([], [], color=COLOUR["outside"], marker="o", ms=3.2, lw=0.9,
               ls=(0, (3, 2)), mfc="white", label="held-out pre-intervention month"),
        Line2D([], [], color=COLOUR["during"], marker="D", ms=3.2, lw=0.9,
               ls=(0, (3, 2)), mfc="white", label="intervention month"),
        Line2D([], [], color=COLOUR["post"], marker="s", ms=3.2, lw=0.9,
               ls=(0, (3, 2)), mfc="white", label="post-intervention month"),
        Patch(facecolor=COLOUR["interval"], alpha=0.18, label="works interval"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False,
               bbox_to_anchor=(0.5, 0.045), columnspacing=0.9, handletextpad=0.45,
               handlelength=1.1, labelspacing=0.4)
    fig.savefig(path, dpi=FIGURE_DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return path


def number(value, digits=3):
    """``value`` with ``digits`` significant digits, or an en dash."""
    if value is None:
        return "–"
    if value == int(value) and abs(value) < 1e12:
        return "%d" % int(value)
    return "%.*g" % (digits, value)


def signed(value, digits=2):
    """A signed number with ``digits`` decimals, or an en dash."""
    return "–" if value is None else "%+.*f" % (digits, value)


def build():
    """The whole report of this asset -- every number the plate draws.

    One pass, one order: the cut is verified, the documented interval is read and checked
    against the export's own strata, the four classes are cut, the reference's frame is
    measured, every chip is standardised in that frame, the bins are built, and only then
    is anything read against the reference. Nothing below re-opens a class, and no
    statistic is computed on a chip that standardised anything.
    """
    chips, audit, meta = load_chips()
    if not chips:
        _fail("the support cut holds no target chip")
    event = event_metadata(chips, meta)
    members, class_months, classification = assign_classes(chips)
    stats = reference_stats(chips)
    weather = weather_reference(chips)
    clipping = standardise(chips, stats, weather)
    strips = assign_parts(chips, event)
    classification["split_months"] = {month: cell["boundary"]
                                      for month, cell in strips["split_months"].items()}
    for month, names in classification["months_spanning_classes"].items():
        if month not in classification["split_months"]:
            _fail("%s holds chips of two classes but no documented boundary falls "
                  "inside it (%s)" % (month, ", ".join(names)))
    moments = monthly(chips)
    if any(month["n"] == 0 or month["cls"] == "mixed" for month in moments):
        _fail("a monthly bin is empty or spans two classes")
    envelope = reference_envelope(chips, moments)
    blocks = {name: class_block(name, members, moments, envelope, stats)
              for name in ("outside", "during", "post")}
    comparison = class_comparison(members, moments)
    orbit = orbit_block(members, moments)
    environment = environment_context(moments, blocks)
    cross = class_contrast(members, "outside", "post")
    verdict = assessment(members, moments, blocks, envelope, orbit, environment, cross,
                         event)

    reference = {
        "window": list(REFERENCE_WINDOW), "state": REFERENCE_STATE,
        "label": CLASS_LABEL["reference"], "note": CLASS_NOTE["reference"],
        "n_chips": len(members["reference"]),
        "n_months": classification["classes"]["reference"]["months"][2],
        "n_intervals": sum(1 for month in moments if month["cls"] == "reference"),
        "outside_window": class_window(members["outside"]),
        "chosen": None,
        "choice": "stated choice, frozen before any result of this plate was read",
        "policy": "pooled over the orbits",
        "standardises": "every held-out chip, through the centre and scale of these "
                        "chips alone",
    }
    states = {"reference": REFERENCE_STATE, "during": DURING_STATE, "post": POST_STATE}
    report = {
        "schema": SCHEMA, "name": NAME, "asset": ASSET, "asset_id": ASSET_ID,
        "label": LABEL, "note": SITE_NOTE,
        "table": audit["table"], "sha256": audit["table_sha256"],
        "inputs": [{"path": rel(path), "sha256": sha256_of(path), "role": role}
                   for path, role in ((TABLE, "the cut this plate reads"),
                                      (META, "the generator's provenance"),
                                      (AOI, "the documented interval"))]
        + ([{"path": EVENTS_SOURCE, "sha256": sha256_of(EVENTS_SOURCE),
             "role": "the external event record (cross-check)"}]
           if os.path.exists(EVENTS_SOURCE) else []),
        "cut": {"rows": audit["rows"], "kept": audit["kept"],
                "excluded": audit["rows"] - audit["kept"],
                "roles": audit["roles"], "refused": audit["refused"],
                "n_dates": audit["n_dates"],
                "window": [chips[0]["day"].isoformat(),
                           chips[-1]["day"].isoformat()]},
        "event": {
            "start": event["start_text"], "end": event["end_text"],
            "span_days": event["span_days"], "kind": event["kind"],
            "from_map": event["from_map"], "map_note": event["map_note"],
            "sources": event["sources"], "external": event["external"],
            "state_windows": event["state_windows"], "states": states,
            "tolerance_days": event["tolerance_days"], "check": event["check"],
            "gaps_days": event["gaps_days"],
            "reading": "read from %s and cross-checked against the generator's meta and "
                       "the export's strata; never typed into this file" % rel(AOI),
        },
        "reference": reference,
        "classification": classification,
        "bins": strips,
        "months": moments,
        "standardisation": {
            "formula": "z = (t(x) - centre) / scale, centre and scale the median and "
                       "1.4826 * MAD of the frozen reference (R3/R4)",
            "transform": dict(TRANSFORM),
            "channels": stats, "weather": weather, "clipping": clipping,
            "reference": reference,
            "rule": "no held-out chip standardises anything, and no weather column "
                    "enters an OCV channel",
        },
        "envelope": envelope,
        "blocks": blocks,
        "comparison": comparison,
        "orbit": orbit,
        "environment": environment,
        "cross": {"post_minus_outside": cross},
        "assessment": verdict,
        "bootstrap": {"draws": BOOT, "block": BLOCK, "seed": SEED,
                      "kind": "circular moving-block, one resample per class"},
        "figure": {
            "file": rel(OUT_JPG), "dpi": FIGURE_DPI,
            "size_mm": [round(COL_WIDTH * 25.4, 1),
                        round(0.92 * COL_WIDTH * 25.4, 1)],
            "panels": [
                {"id": "a1", "block": "behaviour", "shows": "monthly z(gamma2)"},
                {"id": "a2", "block": "behaviour", "shows": "monthly z(P)"},
                {"id": "a3", "block": "behaviour",
                 "shows": "behaviour plane, z(gamma2) x z(P)"},
                {"id": "a4", "block": "behaviour", "shows": "monthly R_B"},
                {"id": "b1", "block": "geometry", "shows": "monthly z(A)"},
                {"id": "b2", "block": "geometry", "shows": "monthly z(F)"},
                {"id": "b3", "block": "geometry",
                 "shows": "geometry plane, z(A) x z(F)"},
                {"id": "b4", "block": "geometry", "shows": "monthly R_G"},
                {"id": "c1", "block": "context",
                 "shows": "monthly co-variates, in their own reference frame"},
                {"id": "c2", "block": "context",
                 "shows": "cadence, chips per bin by orbit"},
                {"id": "c3", "block": "context", "shows": "echo-bearing share per bin"},
                {"id": "c4", "block": "context",
                 "shows": "trajectory matrix, z per dimension and bin"},
            ],
        },
        "provenance": {
            "script": rel(os.path.abspath(__file__)),
            "script_sha256": sha256_of(os.path.abspath(__file__)),
            "template": os.path.relpath(tpl.__file__, ROOT),
            "template_sha256": sha256_of(os.path.abspath(tpl.__file__)),
            "python": "%d.%d.%d" % sys.version_info[:3],
            "matplotlib": matplotlib.__version__,
            "seed": SEED, "deterministic": True,
            "note": "no wall-clock value enters this report: two runs are identical",
        },
    }
    return report


def _md_head(report, add):
    """Title, the site line and the numbers that define the whole reading."""
    event = report["event"]
    ref = report["reference"]
    add("# %s — OCV trajectory on the support-sized chip" % report["label"])
    add("")
    add("**Asset** `%s` · **documented interval** %s → %s (%d d) · **reference state** "
        "`%s` (%s) · **table** `%s` (sha256 `%s`) · **plate** `%s`"
        % (report["asset_id"], event["start"], event["end"], event["span_days"],
           ref["state"], ref["policy"], report["table"], report["sha256"][:16],
           report["figure"]["file"]))
    add("")
    add("```")
    add("standardisation  z = (t(x) - centre) / scale; centre = median_ref, "
        "scale = 1.4826 * MAD_ref, |z| <= %g" % Z_CLIP)
    add("                 t = log for gamma2, P, A, D, F; identity for S (R4)")
    add("blocks           behaviour = %s   geometry = %s   (ocv_reconfigure.BLOCKS)"
        % (list(BLOCKS["behaviour"]), list(BLOCKS["geometry"])))
    add("reference per block  behaviour: echo-bearing reference chips; geometry: every "
        "reference chip (R3)")
    add("bootstrap        %d circular moving-block draws, block %d, seed %d (R5)"
        % (BOOT, BLOCK, SEED))
    add("bins             calendar months, split at a documented boundary where the "
        "record holds chips on both sides")
    add("```")
    add("")
    add("%s" % report["note"].replace("``", "`"))
    add("")


def _md_window(report, add):
    """The cut, the four classes and the reconciliation."""
    cut = report["cut"]
    classes = report["classification"]["classes"]
    add("## The record and its cut")
    add("")
    add("The support cut holds **%d** of the table's **%d** rows (%d refused, %d not a "
        "target chip), over **%d** distinct acquisition dates from %s to %s. The plate "
        "draws **%d** monthly bins."
        % (cut["kept"], cut["rows"], len(cut["refused"]),
           cut["rows"] - cut["kept"] - len(cut["refused"]), cut["n_dates"],
           cut["window"][0], cut["window"][1], len(report["months"])))
    add("")
    add("| class | label in this report | state | chips | ASC | DESC | echo-bearing | "
        "months |")
    add("| --- | --- | --- | --- | --- | --- | --- | --- |")
    for name in CLASSES:
        cell = classes[name]
        add("| `%s` | %s | `%s` | **%d** | %d | %d | %d | %s → %s (%d) |"
            % (name, CLASS_LABEL[name], "/".join(cell["states"]), cell["n"],
               cell["orbits"]["ASCENDING"], cell["orbits"]["DESCENDING"], cell["echo"],
               cell["months"][0], cell["months"][1], cell["months"][2]))
    add("")
    add("The classes are disjoint by construction -- a chip is placed by the export's "
        "state window and, for the reference state only, by the month its acquisition "
        "falls in -- so **%d** kept = **%d** reference + **%d** held-out "
        "pre-intervention + **%d** intervention + **%d** post, with **%d** unassigned."
        % (report["classification"]["total"], classes["reference"]["n"],
           classes["outside"]["n"], classes["during"]["n"], classes["post"]["n"],
           report["classification"]["unassigned"]))
    add("")
    add("The frozen reference is **%s → %s** (inclusive calendar months, %d bins, %d "
        "chips): the record's first full annual cycle, fixed before any result of this "
        "plate was read. The remaining %d chips of the same state (%s → %s) are drawn "
        "hollow and grey as the held-out pre-intervention year; they standardise "
        "nothing."
        % (report["reference"]["window"][0], report["reference"]["window"][1],
           report["reference"]["n_intervals"], report["reference"]["n_chips"],
           classes["outside"]["n"], report["reference"]["outside_window"][0],
           report["reference"]["outside_window"][1]))
    add("")
    if report["bins"]["split_months"]:
        add("A month that holds chips on both sides of a documented boundary is split, "
            "so no bin mixes two sides of the interval: %s. A month is otherwise one "
            "bin, and no line is drawn across a boundary."
            % "; ".join("`%s` split at %s (%s)" % (month, cell["boundary"], cell["part"])
                        for month, cell in report["bins"]["split_months"].items()))
    else:
        add("No month holds chips on both sides of a documented boundary, so every bin "
            "is a plain calendar month and no line is drawn across a boundary.")
    add("")


def _md_reference(report, add):
    """R3/R4: the frame every held-out number is measured in."""
    add("## The reference (R3 / R4)")
    add("")
    add("**%s** — %d chips, %s → %s. Each channel is standardised on the reference "
        "alone." % (report["reference"]["policy"], report["reference"]["n_chips"],
                    report["reference"]["window"][0],
                    report["reference"]["window"][1]))
    add("")
    add("| dimension | block | standardised against | n | t | centre | scale | from | "
        "z range (before the clip) |")
    add("| --- | --- | --- | --- | --- | --- | --- | --- | --- |")
    channels = report["standardisation"]["channels"]
    for dim in DIMS:
        cell = channels[dim]
        add("| `%s` | %s | %s | %d | %s | %s | %s | %s | %s … %s |"
            % (dim, BLOCK_OF[dim], cell["rule"], cell["n"], cell["transform"],
               number(cell["centre"], 4), number(cell["scale"], 4),
               cell["scale_source"], signed(cell["z_min"]), signed(cell["z_max"])))
    add("")
    add("A scale of zero is refused by `scale_of`, so no channel is standardised by a "
        "spread that is not its own. The geometry block is standardised on every "
        "reference chip; the behaviour block only on the reference chips that carry an "
        "echo, because on a chip without one a squared coherence and a scatterer count "
        "have no population to be drawn from (R3).")
    add("")
    clipped = report["standardisation"]["clipping"]
    if clipped["n_clipped"]:
        add("**%d** value(s) lie beyond |z| = %g and are drawn at the clip, never "
            "dropped; they are listed in `standardisation.clipping`, and the range in "
            "the table above is the raw one, before the clip."
            % (clipped["n_clipped"], Z_CLIP))
    else:
        add("No value lies beyond |z| = %g, so no ordinate of the plate is a clipped "
            "measurement." % Z_CLIP)
    add("")
    add("### Environmental context (panel c1)")
    add("")
    add("Standardised on the same reference chips, on their own (identity) scale, and "
        "reported beside the OCV numbers -- never inside them.")
    add("")
    add("| column | n | centre | scale | from |")
    add("| --- | --- | --- | --- | --- |")
    for column, cell in report["standardisation"]["weather"].items():
        add("| `%s` | %d | %s | %s | %s |"
            % (column, cell["n"], number(cell["centre"], 4),
               number(cell["scale"], 4), cell["scale_source"]))
    add("")


def _md_envelope(report, add):
    """R5: the decision distribution is the reference *months*, not the chips."""
    envelope = report["envelope"]
    add("## The reference distribution and the monthly envelope (R5)")
    add("")
    add("The decision distribution is built at **monthly** aggregation: each of the "
        "**%d** reference months contributes one state (the median of its chips), and "
        "the envelope is the empirical q05/q50/q95 of those states. A held-out month is "
        "compared only with this monthly distribution -- never with a single reference "
        "chip and never with the reference maximum, which is printed below for "
        "transparency and is **not** a decision boundary. With %d monthly states the "
        "quantiles are empirical order statistics, not confidence bounds."
        % (envelope["n_reference_months"], envelope["n_reference_months"]))
    add("")
    add("| quantity | monthly reference states | q05 | q50 | q95 | min | max | read as |")
    add("| --- | --- | --- | --- | --- | --- | --- | --- |")
    for key in list(DIMS) + ["R_B", "R_G"]:
        cell = envelope["monthly"].get(key)
        if not cell:
            continue
        add("| `%s` | %d | %s | %s | %s | %s | %s | %s |"
            % (key, cell["n"], signed(cell["q05"]), signed(cell["q50"]),
               signed(cell["q95"]), signed(cell["min"]), signed(cell["max"]),
               "above q95 = unusual magnitude" if key in ("R_B", "R_G")
               else "outside [q05, q95] = signed departure"))
    add("")
    add("The chip-level distributions are stored beside them in the JSON "
        "(`envelope.chip`) and are **descriptive only**: nothing in the held-out reading "
        "compares a monthly state with a single chip, and the chip band is the light "
        "band behind the plate's panels only.")
    add("")


def _md_months(report, add):
    """The trajectory table: the same numbers the plate draws."""
    add("## The monthly trajectory")
    add("")
    add("| bin | class | n | ASC | DESC | echo | " +
        " | ".join("z(%s)" % DIM_MD[dim] for dim in DIMS) + " | R_B | R_G |")
    add("| --- | --- | --- | --- | --- | --- | " + " | ".join("---" for _ in DIMS) +
        " | --- | --- |")
    for month in report["months"]:
        add("| %s | %s | %d | %d | %d | %d | %s | %s | %s |"
            % (month["label"], month["cls"], month["n"], month["n_ascending"],
               month["n_descending"], month["n_echo"],
               " | ".join(signed(month["z"][dim]) for dim in DIMS),
               signed(month["R_B"]), signed(month["R_G"])))
    add("")
    add("A block magnitude is computed **per chip** (the RMS of that chip's standardised "
        "dimensions) and the table and the plate show the **median of those chip "
        "values** -- so `R_B` is the typical joint displacement of a single overpass, "
        "and it is *not* the RMS of the two columns beside it: a bin whose chips scatter "
        "widely in opposite directions has a near-zero median in each dimension and a "
        "large magnitude per chip.")
    add("")


def _md_classes(report, add):
    """The four classes side by side, then each held-out class against the reference."""
    comparison = report["comparison"]
    add("## The frozen reference, the held-out year, the interval and the post period")
    add("")
    add("Each cell is the median over the chips of that class, in the reference's own "
        "scales. The reference column is zero by construction (the standardisation rests "
        "on its median); it is printed so the other columns are read against the frame "
        "they are measured in.")
    add("")
    add("| quantity | reference | held-out pre-intervention | intervention | post |")
    add("| --- | --- | --- | --- | --- |")
    for key in comparison["order"]:
        cells = [comparison["classes"][name]["values"].get(key)
                 for name in CLASSES]
        add("| `%s` | %s | %s | %s | %s |"
            % (key, *(signed(value) for value in cells)))
    add("")
    counts = ", ".join("%s **%d** chips / **%d** bins"
                       % (name, comparison["classes"][name]["n_chips"],
                          comparison["classes"][name]["n_intervals"])
                       for name in CLASSES)
    add("n = " + counts + ".")
    add("")
    add("### Each held-out class against the reference")
    add("")
    add("The interval is R5's seeded circular block bootstrap on "
        "`median(class) - median(reference)`; `outside` counts the class's monthly "
        "states beyond the reference envelope (above q95 for the two magnitudes, outside "
        "[q05, q95] for a signed z).")
    add("")
    #: The reference column is zero for every signed dimension by construction, so it is
    #: not evidence of a mass at zero; only the held-out classes are.
    zeros = [key for key in comparison["order"]
             if any(comparison["classes"][name]["values"].get(key) == 0
                    for name in CLASSES if name != "reference")]
    add("A bound printed as `+0.000` is not a small number rounded to nothing: it is a "
        "draw in which the two medians coincide exactly, which is what a quantity with "
        "a mass at zero does%s. Such an interval reads `excludes 0 = no`, and that is "
        "the honest reading of it."
        % (" (`%s` sit%s at zero in at least one class here)"
           % ("`, `".join(zeros), "s" if len(zeros) == 1 else "")
           if zeros else ""))
    add("")
    for name in ("outside", "during", "post"):
        block = report["blocks"][name]
        add("**%s** — %d chips, %d bins, %s → %s." % (
            CLASS_LABEL[name], block["n_chips"], block["n_intervals"],
            block["window"][0], block["window"][1]))
        add("")
        add("| quantity | median | monthly state | delta | interval | excludes 0 | "
            "months outside | direction |")
        add("| --- | --- | --- | --- | --- | --- | --- | --- |")
        for key in COMPARE_ORDER:
            cell = block["dims"].get(key)
            if cell is None:
                continue
            low, high = cell["ci"]
            add("| `%s` | %s | %s | %s | [%s, %s] | %s | %d | %s |"
                % (key, signed(cell["value"], 3), signed(cell["monthly_median"], 3),
                   signed(cell["delta"], 3), signed(low, 3), signed(high, 3),
                   "yes" if cell["excludes_zero"] else "no",
                   len(cell["reference_comparison"]["outside_envelope"]),
                   cell["direction"]))
        add("")
    add("The bootstrap resamples the reference first and then the class, from a fresh "
        "`random.Random(%d)` per class, so a class's interval does not depend on which "
        "other classes were computed." % SEED)
    add("")


def _md_assessment(report, add):
    """The nine registered questions, the verdict and the extent word."""
    verdict = report["assessment"]
    add("## Repair/intervention assessment")
    add("")
    add("The nine questions below are the registered ones (`QUESTIONS` in the "
        "generator), answered in that order and only from the frozen classes and the "
        "registered rules; the word in parentheses under each heading is the question's "
        "key in the JSON.")
    add("")
    add("The three statements are kept apart throughout:")
    add("")
    add("* **Directional departure** -- a signed monthly median leaves the monthly "
        "reference envelope; a statement about direction, not about size.")
    add("* **Departure magnitude** -- `R_B` or `R_G` above the monthly reference q95; a "
        "statement about the joint size of the displacement, with no direction.")
    add("* **Sustained departure** -- a run of at least **%d** consecutive months that "
        "meet the criterion, decided by the fixed rule and not by the shape of the "
        "curve." % PERSIST_MIN)
    add("")
    for index, answer in enumerate(verdict["answers"], start=1):
        add("**Q%d · %s** (`%s`)"
            % (index, answer["question"].replace("$\\gamma^2$", "γ²"),
               answer["key"]))
        add("")
        add(answer["answer"])
        add("")
    add("**Verdict — %s.** The extent word is **%s**. The verdict is derived from the "
        "intervention interval's two block magnitudes' persistence rules alone "
        "(`verdict_word` on `blocks.during.persistence.R_B` / `R_G`) and the extent word "
        "from the post block's persistence, so neither reads the co-variates, the orbit "
        "split or any interval. It names no damage, no failure precursor, no repair "
        "effect and no cause; the plate measures the observability of the crossing."
        % (verdict["verdict"], verdict["extent"]))
    add("")


def _md_context(report, add):
    """Q8's table, Q9's split, Q7's contrast and the plate's own panel list."""
    orbit = report["orbit"]
    add("## The orbit split (Q9)")
    add("")
    add("The reference is pooled, so the pooled median of a class says nothing about "
        "whether both pass geometries moved. This is a chip-level split in the pooled "
        "frame, and the verdict never reads it: a departure that only one orbit carries "
        "is reported as one pass geometry's departure, not as a site-wide one.")
    add("")
    add("| class | orbit | n | chips' medians: " +
        " | ".join("`%s`" % key for key in orbit["quantities"]) + " |")
    add("| --- | --- | --- | " + " | ".join("---" for _ in orbit["quantities"]) + " |")
    for name in ("outside", "during", "post"):
        for direction in ORBITS:
            cell = orbit["classes"][name]["orbits"][direction]
            add("| %s | %s | %d | %s |"
                % (name, ORBIT_SHORT.get(direction, direction), cell["n"],
                   " | ".join(signed(cell["value"].get(key), 3)
                              for key in orbit["quantities"])))
    add("")
    add("Sign of the departure from the reference, per orbit: " + "; ".join(
        "%s -- %s" % (name, ", ".join("`%s` %s" % (key, orbit["classes"][name]
                                                    ["departure_sign"][key])
                                      for key in orbit["quantities"]))
        for name in ("outside", "during", "post")) + ".")
    add("")
    add("## The post state against the held-out year (Q7)")
    add("")
    cross = report["cross"]["post_minus_outside"]
    add("`%s`, with the same seeded block bootstrap; a descriptive interval on a "
        "difference between two held-out classes, not a test, and not a verdict input."
        % cross["quantity"])
    add("")
    add("| quantity | median(post) - median(held-out) | interval | excludes 0 |")
    add("| --- | --- | --- | --- |")
    for key in COMPARE_ORDER:
        cell = cross["dims"].get(key)
        if cell is None:
            continue
        low, high = cell["ci"]
        add("| `%s` | %s | [%s, %s] | %s |"
            % (key, signed(cell["delta"], 3), signed(low, 3), signed(high, 3),
               "yes" if cell["excludes_zero"] else "no"))
    add("")
    environment = report["environment"]
    add("## The environmental context (Q8)")
    add("")
    add("The co-variates are reported beside the OCV quantities, not inside them: no "
        "weather correction is applied and no causal attribution is made. A pronounced "
        "month is one whose `R_B` or `R_G` exceeds its class's reference q95.")
    add("")
    add("| bin | class | R_B | R_G | z(γ²) | z(P) | " +
        " | ".join("`%s`" % column for column in environment["variables"]) +
        " | above q95 |")
    add("| --- | --- | --- | --- | --- | --- | " +
        " | ".join("---" for _ in environment["variables"]) + " | --- |")
    for row in environment["months"]:
        flagged = ", ".join(name for name, hit in row["above_q95"].items() if hit)
        add("| %s | %s | %s | %s | %s | %s | %s | %s |"
            % (row["month"], row["class"], signed(row["R_B"], 3),
               signed(row["R_G"], 3), signed(row["gamma2"]), signed(row["P"]),
               " | ".join(signed(row["weather_z"].get(column))
                          for column in environment["variables"]),
               flagged or "none"))
    add("")
    add("Pronounced months: %s. Of these, %s also carry a co-variate departure larger "
        "than one robust scale. This is a temporal coincidence, not an explanation."
        % (", ".join(environment["pronounced_months"]) or "none",
           ", ".join(environment["coincident_months"]) or "none"))
    add("")
    add("## The panels")
    add("")
    add("| panel | block | what it draws |")
    add("| --- | --- | --- |")
    for panel in report["figure"]["panels"]:
        add("| **%s** | %s | %s |" % (panel["id"], panel["block"], panel["shows"]))
    add("")
    add("The plate is set at the column's own width -- `figsize = (%.3f in, %.3f in)`, "
        "%d dpi, so %.1f mm x %.1f mm -- and saved with `bbox_inches=\"tight\"`, which "
        "crops the figure to its ink rather than to that rectangle, so the delivered file "
        "is a few millimetres wider than the setting (the colour bar's ticks) and shorter "
        "than it (the empty margin the grid leaves below the key). The header is packed to "
        "the panel block's own width, measured from column one's y label to column four's "
        "right spine, and starts at that same left edge; panel c4's title is broken in "
        "two. Either, left as one line, would be wider than the block, and the tight box "
        "would then widen the *file* to hold the text instead of the text to fit the block."
        % (COL_WIDTH, 0.92 * COL_WIDTH, FIGURE_DPI,
           report["figure"]["size_mm"][0], report["figure"]["size_mm"][1]))
    add("")


def _md_event(report, add):
    """The documented interval and the three statements that had to agree."""
    event = report["event"]
    add("## The documented interval")
    add("")
    add("The interval **%s → %s** (%d d) is read from `%s` and never typed into the "
        "plate. Its stated kind is *%s*, and it does **not** come from the map: %s"
        % (event["start"], event["end"], event["span_days"], rel(AOI), event["kind"],
           (event["map_note"] or "no map note is recorded.").replace("``", "`")))
    add("")
    add("| source | role | outcome |")
    add("| --- | --- | --- |")
    outcome = {
        "committed reference (read)": "read -- the interval is taken from here",
        "generator's own provenance": "agrees with the block above",
        "operator record (agrees)": "agrees on both dates -- cross-check passed",
        "operator record (unreadable)": "unavailable -- recorded as a note, not a pass",
    }
    for source in event["sources"]:
        add("| `%s` | %s | %s |"
            % (source["name"], source["role"],
               outcome.get(source["role"], source["role"])))
    add("")
    external = event["external"]
    if external and external.get("role") == "operator record (unreadable)":
        add("The external operator record exists but holds **no entry**, so the "
            "cross-check it was meant to provide is **unavailable**. The plate records "
            "that as a note and does not claim the cross-check as passed: the interval "
            "rests on the committed AOI block, on the generator's own meta and on the "
            "export's strata, not on that file.")
        add("")
    elif external:
        add("The external operator record is not taken on trust: it is read as data "
            "(`%s`, sha256 `%s`), its two entries are matched by label to the two "
            "documented dates, and any disagreement is a failure of the run, not a "
            "note. It names `sanierung_beginn` %s and `sanierung_ende` %s -- the same "
            "two dates, to the day."
            % (external["name"], external["sha256"][:16],
               external["begin"], external["finish"]))
        add("")
    add("The metadata is also checked against the data it labels, not only against "
        "another document: the export's strata must be cut on these dates, within a "
        "tolerance of **%d d** -- the record's own longest acquisition gap, because a "
        "boundary can only be located to within one gap. The observed gaps are %s."
        % (event["tolerance_days"],
           ", ".join("%s %d d" % (pair, gap)
                     for pair, gap in event["gaps_days"].items())))
    add("")
    add("| state | first chip | last chip | chips |")
    add("| --- | --- | --- | --- |")
    for state in STATE_ORDER:
        cell = event["state_windows"][state]
        add("| `%s` | %s | %s | %d |"
            % (state, cell["first_date"], cell["last_date"], cell["n"]))
    add("")
    add("%s." % event["check"])
    add("")


def _md_limits(report, add):
    """What the plate does not claim, the pins and how to reproduce it."""
    add("## What the plate does not claim")
    add("")
    add("* **No damage, failure or safety verdict.** Every ordinate is a deviation from "
        "the crossing's *own* reference, in that reference's robust scales; two plates "
        "are not comparable in magnitude, only in shape and in where the movement sits.")
    add("* **No attribution to the works.** The interval is the record's; a departure "
        "that coincides with it in time is not evidence that the works caused it, and no "
        "weather correction is applied.")
    add("* **No significance on the curve.** The intervals are differences of medians, "
        "not bands on the curve: the monthly medians are drawn with their `n` and "
        "without one.")
    add("* **No directional score.** `R_B` and `R_G` keep their frozen definitions, and "
        "the verdict never reads the orbit split, the co-variates or any interval.")
    add("* **The reference maximum is not a boundary.** The envelope is the q05/q50/q95 "
        "of the monthly reference states; the maximum is printed only so a reader can "
        "see how far q95 sits from it.")
    add("* **No interpolation.** A boundary, a hole and a class change break the line: "
        "nothing is drawn across them, and no bin mixes two classes.")
    add("* **A month is a bin, not a date.** Its median can mix both pass directions; "
        "panel c2 says how many chips of each the bin holds.")
    add("")
    add("## Pins and reproduction")
    add("")
    add("| pin | value |")
    add("| --- | --- |")
    for key, value in sorted(PINS.items()):
        add("| `%s` | `%s` |" % (key, value))
    add("")
    add("`--verify` re-derives every number of this report from `%s` and fails if the "
        "JSON, the report and the plate disagree with the table; `--selftest` checks the "
        "rules and the pins above without writing anything. The bootstrap is seeded "
        "(%d, block %d, %d draws), no wall-clock value enters any output, and the plate "
        "is drawn twice and compared byte for byte, so the JSON, the report and the "
        "figure reproduce deterministically."
        % (report["table"], SEED, BLOCK, BOOT))
    add("")
    add("```bash")
    add("bash figures/werder_havel/fig_werder_havel_ocv_trajectory.sh            # "
        "JSON, report, plate")
    add("bash figures/werder_havel/fig_werder_havel_ocv_trajectory.sh --no-figures # "
        "numbers only")
    add("bash figures/werder_havel/fig_werder_havel_ocv_trajectory.sh --verify      # "
        "re-derive and compare")
    add("bash figures/werder_havel/fig_werder_havel_ocv_trajectory.sh --selftest    # "
        "rules and pins, nothing written")
    add("```")
    add("")
    add("The generator of the cut is "
        "`code/werder_havel/werder_havel_ocv_channels_csv.py`; the template this plate "
        "follows is `%s` (sha256 `%s`), imported rather than copied, so the seed, the "
        "transforms, the block definitions and the monthly-bin machinery cannot drift "
        "between the two plates."
        % (report["provenance"]["template"],
           report["provenance"]["template_sha256"][:16]))
    add("")


def report_md(report):
    """The report -- the same numbers as the plate, never a new one."""
    out = []
    add = out.append
    _md_head(report, add)
    _md_event(report, add)
    _md_window(report, add)
    _md_reference(report, add)
    _md_envelope(report, add)
    _md_months(report, add)
    _md_classes(report, add)
    _md_assessment(report, add)
    _md_context(report, add)
    _md_limits(report, add)
    return "\n".join(out).rstrip("\n") + "\n"


def _rules(report):
    """The registered rules, re-derived from one built report.

    Returns ``(name, ok, detail)`` triples. Nothing here repeats a number typed into
    ``build()``: every check recomputes what it asserts from the report's own leaves, so a
    silent change in the analysis fails the check instead of moving it.
    """
    rules = []

    def rule(name, ok, detail=""):
        rules.append((name, bool(ok), detail))

    classes = report["classification"]["classes"]
    months = report["months"]
    env = report["envelope"]
    stats = report["standardisation"]["channels"]
    reference_echo = classes["reference"]["echo"]
    reference_n = classes["reference"]["n"]

    rule("R1 the cut is the support cut's target role only",
         report["cut"]["kept"] + report["cut"]["excluded"] == report["cut"]["rows"],
         "%d kept + %d excluded = %d rows"
         % (report["cut"]["kept"], report["cut"]["excluded"], report["cut"]["rows"]))
    rule("R2 the four classes are disjoint",
         report["classification"]["disjoint"]
         and set(report["classification"]["months_spanning_classes"])
         <= set(report["classification"]["split_months"]),
         "every chip in exactly one class; a month holds two classes only where a "
         "documented boundary splits it (%s)"
         % (", ".join("%s at %s" % (month, boundary) for month, boundary
                      in sorted(report["classification"]["split_months"].items()))
            or "no split month"))
    rule("R2 the classes reconcile to every kept chip, with nothing unassigned",
         sum(classes[name]["n"] for name in CLASSES) == report["classification"]["total"]
         == report["cut"]["kept"] and report["classification"]["unassigned"] == 0,
         " + ".join("%s %d" % (name, classes[name]["n"]) for name in CLASSES)
         + " = %d" % report["cut"]["kept"])
    rule("R2 no bin mixes two classes or the two sides of a boundary",
         all(month["cls"] in CLASSES and month["n"] > 0 for month in months)
         and all(sum(month["sides"].values()) == month["n"] for month in months)
         and all(sum(1 for count in month["sides"].values() if count) == 1
                 for month in months),
         "%d bins, split months %s"
         % (len(months), sorted(report["bins"]["split_months"]) or "none"))
    rule("R3 behaviour is standardised on the echo-bearing reference chips only",
         stats["gamma2"]["n"] == reference_echo == stats["P"]["n"]
         and stats["gamma2"]["rule"] == "echo-bearing reference chips",
         "%d of %d reference chips carry an echo" % (reference_echo, reference_n))
    rule("R3 geometry is standardised on every reference chip",
         all(stats[dim]["n"] == reference_n
             and stats[dim]["rule"] == "all reference chips"
             for dim in BLOCKS["geometry"]),
         ", ".join("%s %d" % (dim, stats[dim]["n"]) for dim in BLOCKS["geometry"]))
    rule("R4 the transform is log for five channels and identity for S",
         report["standardisation"]["transform"] == TRANSFORM,
         ", ".join("%s %s" % (dim, TRANSFORM[dim]) for dim in DIMS))
    rule("R4 no channel is standardised by a scale of zero",
         all(stats[dim]["scale"] > 0 for dim in DIMS)
         and all(cell["scale"] > 0
                 for cell in report["standardisation"]["weather"].values()),
         "smallest scale %s"
         % number(min([stats[dim]["scale"] for dim in DIMS]
                      + [cell["scale"] for cell in
                         report["standardisation"]["weather"].values()]), 3))

    reference_bins = [month for month in months
                      if month["cls"] == "reference" and month["n"]]
    rule("R5 the envelope is built from the reference months, not the chips",
         env["decision_level"] == "monthly"
         and env["n_reference_months"] == len(reference_bins)
         and env["months"] == [month["label"] for month in reference_bins]
         and all({pair["month"] for pair in cell["values"]} == set(env["months"])
                 for cell in env["monthly"].values()),
         "%d monthly reference states, %d reference bins"
         % (env["n_reference_months"], len(reference_bins)))
    rule("R5 the chip-level distribution is marked descriptive only",
         "descriptive" in env["chip_note"] and "never compared" in env["chip_note"],
         env["chip_note"])
    rule("R5 the bootstrap is the registered one, in every block",
         (report["bootstrap"]["draws"], report["bootstrap"]["block"],
          report["bootstrap"]["seed"]) == (BOOT, BLOCK, SEED)
         and all(report["blocks"][name]["seed"] == SEED
                 and report["blocks"][name]["draws"] == BOOT
                 and report["blocks"][name]["block"] == BLOCK
                 for name in ("outside", "during", "post"))
         and (report["cross"]["post_minus_outside"]["seed"],
              report["cross"]["post_minus_outside"]["block"]) == (SEED, BLOCK),
         "%d draws, block %d, seed %d" % (BOOT, BLOCK, SEED))
    rule("Q9 both pass geometries are present in every held-out class",
         all(all(cell["n"] > 0 for cell in
                 report["orbit"]["classes"][name]["orbits"].values())
             for name in ("outside", "during", "post")),
         ", ".join("%s %s" % (name, "/".join(
             "%s %d" % (ORBIT_SHORT[orbit],
                        report["orbit"]["classes"][name]["orbits"][orbit]["n"])
             for orbit in ORBITS))
             for name in ("outside", "during", "post")))
    rule("Q9 the departure's sign is stated for every registered quantity",
         all(set(report["orbit"]["classes"][name]["departure_sign"])
             == set(report["orbit"]["quantities"])
             for name in ("outside", "during", "post")),
         ", ".join(sorted({report["orbit"]["classes"]["during"]["departure_sign"][key]
                           for key in report["orbit"]["quantities"]})))
    event = report["event"]
    rule("the documented interval is read from the committed metadata, not typed",
         (event["start"], event["end"]) == PINS["event"]
         and event["from_map"] is False
         and event["span_days"] == (dt.date.fromisoformat(event["end"])
                                    - dt.date.fromisoformat(event["start"])).days,
         "%s -> %s (%d d), from_map %s"
         % (event["start"], event["end"], event["span_days"], event["from_map"]))
    rule("the export's strata are cut on those two dates",
         all((event["state_windows"][state]["first_date"],
              event["state_windows"][state]["last_date"]) == PINS["state_windows"][state]
             for state in STATE_ORDER),
         ", ".join("%s %s -> %s" % (state, event["state_windows"][state]["first_date"],
                                    event["state_windows"][state]["last_date"])
                   for state in STATE_ORDER))
    external = event["external"]
    if external is None:
        external_ok, external_note = True, "absent -- no cross-check is claimed"
    elif external.get("role") == "operator record (unreadable)":
        external_ok = any(source["role"] == "operator record (unreadable)"
                          for source in event["sources"])
        external_note = "empty, and reported as a note rather than a pass"
    else:
        external_ok = ((external.get("begin"), external.get("finish"))
                       == (event["start"], event["end"]))
        external_note = "present, read as data, and agreeing on both dates"
    rule("an external event record can only agree or be reported as a note",
         external_ok, external_note)
    verdict = report["assessment"]
    during = report["blocks"]["during"]
    rule("the verdict is one of the registered words, and is derived",
         verdict["verdict"] in ALLOWED_VERDICTS
         and verdict["extent"] in ALLOWED_EXTENTS
         and verdict["verdict"] == verdict_word(during["persistence"]["R_B"],
                                                during["persistence"]["R_G"]),
         "%r / %r" % (verdict["verdict"], verdict["extent"]))
    rule("the pins are the record's own census",
         report["cut"]["rows"] == PINS["table_rows"]
         and report["cut"]["kept"] == PINS["cut_rows"]
         and {name: classes[name]["n"] for name in CLASSES} == PINS["classes"]
         and report["reference"]["n_months"] == PINS["reference_months"]
         and len(months) == PINS["monthly_bins"]
         and list(report["reference"]["window"]) == list(PINS["reference_window"]),
         "%d/%d rows, %s, %d bins"
         % (report["cut"]["kept"], report["cut"]["rows"],
            "/".join(str(classes[name]["n"]) for name in CLASSES), len(months)))
    return rules


def selftest():
    """The rules and the pins, on a fresh build; nothing is written."""
    started = time.time()
    report = build()
    again = build()
    rules = _rules(report)
    rules.append(("the JSON reproduces deterministically",
                  json.dumps(round_floats(report), sort_keys=True)
                  == json.dumps(round_floats(again), sort_keys=True),
                  "two builds compared field by field"))
    rules.append(("the report reproduces deterministically",
                  report_md(report) == report_md(again),
                  "%d characters" % len(report_md(report))))
    width = max(len(name) for name, _ok, _detail in rules)
    failed = 0
    for name, ok, detail in rules:
        if not ok:
            failed += 1
        line = "%s  %-*s" % ("PASS" if ok else "FAIL", width, name)
        if detail:
            line += "  " + detail
        print(line)
    print("\n%d checks, %d failed (%.1f s, nothing written)"
          % (len(rules), failed, time.time() - started))
    return 1 if failed else 0


def publish(paths=None, draw_figures=True):
    """Write the JSON, the report and -- unless asked not to -- the plate."""
    report = build()
    paths = dict(paths or {"json": OUT_JSON, "md": OUT_MD, "jpg": OUT_JPG})
    for key in ("json", "md"):
        os.makedirs(os.path.dirname(os.path.abspath(paths[key])), exist_ok=True)
    with open(paths["json"], "w", encoding="utf-8") as handle:
        json.dump(round_floats(report), handle, indent=1, sort_keys=True)
        handle.write("\n")
    with open(paths["md"], "w", encoding="utf-8") as handle:
        handle.write(report_md(report))
    if draw_figures:
        os.makedirs(os.path.dirname(os.path.abspath(paths["jpg"])), exist_ok=True)
        draw(report, paths["jpg"])
    return report, paths


def verify():
    """Re-derive everything and compare it with the committed JSON, report and plate.

    Three comparisons, in order of what a reader can check: the numbers (the JSON is
    rebuilt and diffed field by field), the prose (the report must be the string this
    table produces), and the plate (drawn twice, compared byte for byte, and then
    compared with the committed file). On top of that every registered rule is re-run on
    the fresh report, so a JSON that matches a *wrong* analysis still fails.
    """
    report = build()
    problems = []
    for name, ok, detail in _rules(report):
        if not ok:
            problems.append("%s (%s)" % (name, detail))
    if not os.path.exists(OUT_JSON):
        problems.append("%s: missing -- run the script first"
                        % os.path.relpath(OUT_JSON, ROOT))
    else:
        with open(OUT_JSON, encoding="utf-8") as handle:
            committed = json.load(handle)
        difference = first_difference(committed, round_floats(report))
        if difference:
            problems.append("%s: %s" % (os.path.relpath(OUT_JSON, ROOT), difference))
    if not os.path.exists(OUT_MD):
        problems.append("%s: missing -- run the script first"
                        % os.path.relpath(OUT_MD, ROOT))
    else:
        with open(OUT_MD, encoding="utf-8") as handle:
            committed_text = handle.read()
        if committed_text != report_md(report):
            problems.append("%s: the report does not match the table"
                            % os.path.relpath(OUT_MD, ROOT))
    if not os.path.exists(OUT_JPG):
        problems.append("%s: missing -- run the script first"
                        % os.path.relpath(OUT_JPG, ROOT))
    else:
        digests = []
        for _ in range(2):
            handle, temporary = tempfile.mkstemp(suffix=".jpg")
            os.close(handle)
            try:
                draw(report, temporary)
                digests.append(sha256_of(temporary))
            finally:
                os.unlink(temporary)
        if len(set(digests)) != 1:
            problems.append("%s: the plate is not deterministic (two draws differ)"
                            % os.path.relpath(OUT_JPG, ROOT))
        elif sha256_of(OUT_JPG) != digests[0]:
            problems.append("%s: the committed plate is not the one this report draws"
                            % os.path.relpath(OUT_JPG, ROOT))
    return problems


def report_brief(report):
    """One line: the cut, the four classes and the two derived words."""
    classes = report["classification"]["classes"]
    verdict = report["assessment"]
    return ("%d of %d rows, %s; reference %s %d (%s -> %s); works %s -> %s; "
            "%s / %s"
            % (report["cut"]["kept"], report["cut"]["rows"],
               " / ".join("%s %d" % (name, classes[name]["n"]) for name in CLASSES),
               report["reference"]["state"], report["reference"]["n_chips"],
               report["reference"]["window"][0], report["reference"]["window"][1],
               report["event"]["start"], report["event"]["end"],
               verdict["verdict"], verdict["extent"]))


def main(argv=None):
    started = time.time()
    parser = argparse.ArgumentParser(
        prog=NAME + ".py",
        description="The OCV trajectory of the Havelbruecke Werder (A 10) across the "
                    "documented 2024/2025 renewal: four classes, a frozen annual "
                    "reference, the works interval as a band.")
    parser.add_argument("--json", default=OUT_JSON, help="path of the result JSON")
    parser.add_argument("--md", default=OUT_MD, help="path of the report")
    parser.add_argument("--jpg", default=OUT_JPG, help="path of the plate (JPEG)")
    parser.add_argument("--no-figures", action="store_true",
                        help="write the JSON and the report only (no matplotlib output)")
    parser.add_argument("--selftest", action="store_true",
                        help="check the rules and the pins, write nothing")
    parser.add_argument("--verify", action="store_true",
                        help="re-derive and compare with the committed files")
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args(argv)

    if args.selftest:
        return selftest()

    if args.verify:
        problems = verify()
        if problems:
            for problem in problems:
                print("FAIL  %s" % problem)
            print("\n%d problem(s): the committed files are not what this table "
                  "produces" % len(problems))
            return 1
        if not args.quiet:
            print("OK    %s" % report_brief(build()))
            print("OK    the JSON, the report and the plate reproduce from the table")
        return 0

    report, written = publish({"json": args.json, "md": args.md, "jpg": args.jpg},
                              draw_figures=not args.no_figures)
    if not args.quiet:
        print("%s: %s" % (NAME, report_brief(report)))
        print("  " + os.path.relpath(written["json"], ROOT))
        print("  " + os.path.relpath(written["md"], ROOT))
        if not args.no_figures:
            print("  " + os.path.relpath(written["jpg"], ROOT))
    print("%s (%.1f s)" % ("DONE" if not args.no_figures else "DONE (numbers only)",
                           time.time() - started))
    return 0


if __name__ == "__main__":
    sys.exit(main())

