#!/usr/bin/env python3
"""OCV reconfiguration against a controlled state change — three registered pairs.

Two questions, one metric. **Does the magnitude of OCV reconfiguration grow with
an independently defined state dose, and is the direction of that reconfiguration
the same in both assets?** The quantity is the distance from each experiment's
*own* reference state; raw OCV values are never compared between two assets.

The dose pair is LUMO and Bautzen; the bridge pair is Carolabrücke and Ponte
Morandi, read across a collapse instead of a dose (amendment 4); the tower pair is
Haltern and KDLO, read across a collapse and a rebuild (amendment 5). The design is
fixed in ``code/sds/ocv_reconfigure_spec.md`` and was written before this module
ran. Short version:

* **dose** — LUMO's dose is measured, not asserted: the first-higher natural
  frequency of the tower drops 0 % (healthy) / 7.41 % (DAM 3 and DAM 4) / 12.96 %
  (DAM 6), so the ordered test has three levels and DAM 3 and DAM 4 are one level
  (``data/lumo/reference/lumo_state_frequencies.json``). Bautzen's ladder
  ``RS < DS1 < DS2`` is ordinal only. The bridges and the towers have no dose at
  all: two states, one binary event, and the ordered test reduces to its
  two-group case;
* **response** — ``r = median(g(X_state)) - median(g(X_ref))`` on the transformed
  scale ``g`` (``log`` for positive dimensions, identity for ``S`` and Bautzen's
  floor-corrected ``gamma2``), standardised by the **MAD of the reference state**
  with a fixed fallback chain;
* **Bautzen** is analysed on the **paired bridge-minus-control difference**, because
  its states are also its calendar;
* **blocks** — ``mask = [A, F]`` and ``other = [gamma2, P, S]``. ``D = A/bbox_area``
  holds in 178 of 178 LUMO rows, so the literal ``[A, D, F]`` block is a
  sensitivity analysis, not the primary one;
* **channels** — a record carries the orbit directions it was actually acquired in.
  KDLO was acquired DESCENDING only, so its run has one channel, one H1 cell pair
  per block and a two-test Holm family. A channel a record does not carry is
  *declared* (``channels_absent``), never filled and never counted as a null;
* **H1** — Jonckheere-Terpstra on the *per-unit* magnitudes (the reference's median
  magnitude is 0 by construction, so a test on state medians would be degenerate),
  Holm over the asset's own family — two cells per channel it carries — and a
  moving-block bootstrap for the intervals;
* **H2** — cosine between the standardised response vectors, block by block, with a
  random-direction null; γ² is reported with and without, because it is a different
  quantity in the two experiments. Descriptive in every run, never a hypothesis.
  It needs both assets on the same channel, so a pair whose records share only one
  orbit direction carries its cosines on that one;
* **excluded rows** — an acquisition whose state the registry excludes is removed
  before any statistic and *counted* in the audit, never folded into a neighbouring
  state: KDLO's ``single`` is the one acquisition at the collapse itself, and one
  date can carry no median of its own.

A cell whose ``n`` cannot fill the bootstrap block is flagged wherever it is quoted
(``✱`` on the plates, a named list in the report): its interval then carries the
reference's uncertainty only.

Usage
-----
  python3 code/sds/ocv_reconfigure.py
  python3 code/sds/ocv_reconfigure.py --sites carola,morandi
  python3 code/sds/ocv_reconfigure.py --sites haltern,kdlo
  python3 code/sds/ocv_reconfigure.py --sites kdlo
  python3 code/sds/ocv_reconfigure.py --selftest
  python3 code/sds/ocv_reconfigure.py --verify

Outputs
-------
  data/sds/<stem>.json, <stem>.md, <stem>_table.csv
  figures/sds/<figure stem>.jpg + .pdf + .eps   (170 mm wide, 630 dpi raster)
  figures/sds/<figure stem>.md                  (the panel guide + the plate audit)
where ``<stem>`` is ``ocv_reconfigure`` for the dose pair, ``ocv_carola_morandi``
for the bridge pair, ``ocv_haltern_kdlo`` for the tower pair, and ``ocv_<site>``
for a single asset.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import statistics
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
for _path in (os.path.join(REPO, "code", "lumo"), os.path.join(REPO, "code")):
    if _path not in sys.path:
        sys.path.insert(0, _path)

from lumo_ocv_stats import _jt_statistic  # noqa: E402  (verbatim port, see module)
import figure_style as style             # noqa: E402  (the one typography and width)

SCHEMA = "ocv-reconfigure-v1"
SEED = 7
BOOT = 2000
PERM = 20000
BLOCK_BOOT = 3
FIGURE_DPI = 630

#: orbit_direction -> the observation channel used everywhere in the reports.
CHANNEL = {"ASCENDING": "ASC", "DESCENDING": "DESC", "ASC": "ASC", "DESC": "DESC"}

DIMS = ("A", "D", "F", "P", "S", "gamma2")

#: dimension -> column, per site (the two exports name gamma2 differently).
COLUMNS = {
    "lumo": {"A": "A", "D": "D", "F": "F", "P": "P", "S": "S",
             "gamma2": "coherence_gamma2"},
    "bautzen": {"A": "A", "D": "D", "F": "F", "P": "P", "S": "S",
                "gamma2": "gamma2_band_floor_corrected"},
    "carola": {"A": "A", "D": "D", "F": "F", "P": "P", "S": "S",
               "gamma2": "gamma2_mask"},
    "morandi": {"A": "A", "D": "D", "F": "F", "P": "P", "S": "S",
                "gamma2": "gamma2_mask"},
    # Both tower chips carry `gamma2_mask` and also Haltern's asset-point
    # `gamma2_export`. The masked coherence is the one the window rule defines (the
    # mask is the support-sized chip's own), so `gamma2_export` is left unused: it is
    # a single point of the scene, not the site's masked coherence.
    "haltern": {"A": "A", "D": "D", "F": "F", "P": "P", "S": "S",
                "gamma2": "gamma2_mask"},
    "kdlo": {"A": "A", "D": "D", "F": "F", "P": "P", "S": "S",
             "gamma2": "gamma2_mask"},
}

#: The scale each dimension is analysed on. `S` is a non-negative distance and
#: Bautzen's gamma2 is floor-corrected and signed, so both stay on identity.
TRANSFORM = {
    "lumo": {"A": "log", "D": "log", "F": "log", "P": "log", "gamma2": "log",
             "S": "identity"},
    "bautzen": {"A": "log", "D": "log", "F": "log", "P": "log", "gamma2": "identity",
                "S": "identity"},
    "carola": {"A": "log", "D": "log", "F": "log", "P": "log",
               "gamma2": "log", "S": "identity"},
    # Morandi's export names the same coherence channel `gamma2_mask`, and it is
    # positive in every row, so it stays on the log scale like Carola's.
    "morandi": {"A": "log", "D": "log", "F": "log", "P": "log",
                "gamma2": "log", "S": "identity"},
    # The tower chips: every dimension is positive in every row of both, so all six
    # except `S` are analysed on the log scale.
    "haltern": {"A": "log", "D": "log", "F": "log", "P": "log",
                "gamma2": "log", "S": "identity"},
    "kdlo": {"A": "log", "D": "log", "F": "log", "P": "log",
             "gamma2": "log", "S": "identity"},
}

#: The two response blocks (see the dict below); see the 2026-10-22 amendment.
#: How a dimension is written when it is a *label* (the figures), never as a key.
DIM_DISPLAY = {"gamma2": r"$\gamma^2$", "A": "A", "D": "D", "F": "F",
               "P": "P", "S": "S"}


#: The two response blocks, in the wording the paper uses.  The keys keep the
#: pre-registration's spelling -- see ``BLOCK_DISPLAY`` for what a reader sees.
#:
#: * ``behaviour`` -- cross-acquisition echo behavior: `gamma2` is the consistency
#:   of the echo across acquisitions, `P` its persistence/continuity across
#:   acquisitions. `P` belongs here and not in the geometry block, even though its
#:   computation depends on valid echo detections.
#: * ``geometry`` -- spatial echo-mask geometry: `A` extent (how much is
#:   observable), `D` occupancy, `F` fragmentation (how spatially organised the
#:   mask is), `S` the spatial shift. `D` is kept out of the figure's axes because
#:   it partly inherits `A` algebraically, and `S` is reported as a supplementary
#:   coordinate because a shift describes position, not morphology -- both stay in
#:   the block magnitude.
BLOCKS = {
    "behaviour": ("gamma2", "P"),
    "geometry": ("A", "D", "F", "S"),
}

#: How a block is *printed*.  The keys above are the pre-registration's tokens and
#: stay as they are: the other generators index ``BLOCKS`` and this module keys its
#: report, its H1/H2 ids and its plate audit on them.  The manuscript is set in
#: American English, so every name that reaches a reader says ``behavior`` -- no
#: key, no field name in the JSON is routed through this table, only labels.
BLOCK_DISPLAY = {"behaviour": "behavior", "geometry": "geometry"}

#: The two dimensions drawn as the plane of the geometry block; `D` and `S` are
#: annotated beside it (see ``BLOCKS``).
GEOMETRY_PLANE = ("A", "F")

#: Where a single-table asset keeps its state, orbit and date.
UNIT_SHAPE = {
    "lumo": {"state": "damage_label", "orbit": "orbit_direction",
             "date": "acquisition_date"},
    "carola": {"state": "state", "orbit": "orbit_direction",
               "date": "acquisition_ts"},
    "morandi": {"state": "state", "orbit": "orbit_direction",
                "date": "acquisition_ts"},
    "haltern": {"state": "state", "orbit": "orbit_direction",
                "date": "acquisition_ts"},
    "kdlo": {"state": "state", "orbit": "orbit_direction",
             "date": "acquisition_ts"},
}

REFERENCE = {"lumo": "healthy", "bautzen": "RS",
             "carola": "pre-collapse (healthy)",
             "morandi": "pre-collapse (healthy)",
             "haltern": "pre-collapse (healthy)",
             "kdlo": "pre"}

STATES = {"lumo": ["healthy", "DAM 3", "DAM 4", "DAM 6"],
          "bautzen": ["RS", "DS1", "DS2"],
          "carola": ["pre-collapse (healthy)", "post-collapse"],
          "morandi": ["pre-collapse (healthy)", "post-collapse"],
          "haltern": ["pre-collapse (healthy)", "post-collapse"],
          "kdlo": ["pre", "rebuild"]}

#: Rows the registry removes before any statistic, with the reason it gives. An
#: excluded state is *counted* in the audit and never folded into a neighbour. KDLO's
#: `single` is the one acquisition at the collapse itself (2022-12-23): one date
#: cannot carry a median of its own, and it belongs to neither the standing nor the
#: rebuilt tower. Registering it would create an n = 1 cell; dropping it silently
#: would hide a decision, so it is excluded here and reported.
EXCLUDED = {
    "kdlo": {"single": "the collapse itself, 2022-12-23, and one date cannot carry a "
                       "median"},
}

#: How the two-state contrast is called in a verdict, and how a single asset's title
#: reads. KDLO's second state is a rebuild, and its report says so.
CONTRAST = {"carola": "pre-to-post", "morandi": "pre-to-post",
            "haltern": "pre-to-post", "kdlo": "pre-to-rebuild"}
SITE_HEADLINE = {"carola": "across a collapse", "morandi": "across a collapse",
                 "haltern": "across a collapse",
                 "kdlo": "across a collapse and a rebuild"}

#: Dose levels for the ordered test. LUMO's dose is measured (percent drop of the
#: first-higher mode); DAM 3 and DAM 4 share it to one FFT bin. Bautzen has no
#: measured dose, so its levels are the ordinal ladder.
LEVELS = {
    "lumo": [("L0", ("healthy",), 0.0),
             ("L1", ("DAM 3", "DAM 4"), 7.41),
             ("L2", ("DAM 6",), 12.96)],
    "bautzen": [("RS", ("RS",), 0.0), ("DS1", ("DS1",), 1.0), ("DS2", ("DS2",), 2.0)],
    # Carola has two states and nothing in between: the ordered test collapses to a
    # one-sided rank test between them, and dR is the post-minus-pre step. Morandi is
    # the same shape: 2018-08-14 is the collapse, and no ladder exists.
    "carola": [("pre", ("pre-collapse (healthy)",), 0.0),
               ("post", ("post-collapse",), 1.0)],
    "morandi": [("pre", ("pre-collapse (healthy)",), 0.0),
                ("post", ("post-collapse",), 1.0)],
    # The towers carry the same two-level shape. Haltern's `post` is the fallen
    # turbine; KDLO's second level is the *rebuilt* mast, which is why it is called
    # `rebuild` and not `post` — the level names travel from here into the reports.
    "haltern": [("pre", ("pre-collapse (healthy)",), 0.0),
                ("post", ("post-collapse",), 1.0)],
    "kdlo": [("pre", ("pre",), 0.0), ("rebuild", ("rebuild",), 1.0)],
}

#: The state used for the cross-experiment direction (and the middle one).
LARGEST = {"lumo": "DAM 6", "bautzen": "DS2",
           "carola": "post-collapse", "morandi": "post-collapse",
           "haltern": "post-collapse", "kdlo": "rebuild"}
MIDDLE = {"lumo": ("DAM 3", "DAM 4"), "bautzen": ("DS1",),
          "carola": ("pre-collapse (healthy)",),
          "morandi": ("pre-collapse (healthy)",),
          "haltern": ("pre-collapse (healthy)",),
          "kdlo": ("pre",)}

#: Below this many units a cell is marked thin and must be quoted with its n.
SMALL_N = 6

#: Below this many reference units the reference MAD is itself too unstable to
#: standardise against, and the pooled MAD takes over (amendment 2026-10-22).
SCALE_MIN_REF = 10

TABLE = {"lumo": os.path.join(REPO, "data", "lumo", "lumo_ocv_channels.csv"),
         "bautzen": os.path.join(REPO, "data", "bautzen",
                                 "bautzen_ocv_channels.csv"),
         # Carola is read on the support-sized chip, where only Girder A survives
         # the cut (436 chips of the committed 1 717).
         "carola": os.path.join(REPO, "data", "carola", "support",
                                 "carola_ocv_channels.csv"),
         # Morandi is read on its support-sized chip too (207 of the committed 232),
         # which is what makes the two bridges a comparison of one window rule.
         "morandi": os.path.join(REPO, "data", "morandi", "support",
                                 "morandi_ocv_channels.csv"),
         # The towers are read on the same rule's chip: Haltern 107 of 120 committed
         # overpasses, KDLO 60 of 60 (every one of its acquisitions lies inside the
         # support).
         "haltern": os.path.join(REPO, "data", "haltern", "support",
                                 "haltern_ocv_channels.csv"),
         "kdlo": os.path.join(REPO, "data", "kdlo", "support",
                              "kdlo_ocv_channels.csv")}

#: The pair the module was built for; ``--sites`` may select any registered asset.
DEFAULT_SITES = ("lumo", "bautzen")

#: The two collapsed bridges, read with the same metric and against their own
#: pre-collapse reference. A pair, not a new hypothesis class: amendment 4.
BRIDGE_PAIR = ("carola", "morandi")

#: The two towers, read against their own pre-collapse reference: one came down and
#: stayed down, the other came down and was rebuilt. Amendment 5.
TOWER_PAIR = ("haltern", "kdlo")

#: The pairs that carry a descriptive H2; every other selection is per-asset only.
PAIRS = (DEFAULT_SITES, BRIDGE_PAIR, TOWER_PAIR)

#: The pairs whose assets each have exactly two states, so H2 reduces to the one
#: largest-state cosine per channel and block.
TWO_STATE_PAIRS = (BRIDGE_PAIR, TOWER_PAIR)

FREQUENCIES = os.path.join(REPO, "data", "lumo", "reference",
                           "lumo_state_frequencies.json")


def pair_of(sites):
    """The registered pair a run stands for, or ``None`` for a single asset."""
    sites = tuple(sites)
    for pair in PAIRS:
        if set(sites) == set(pair):
            return pair
    return None


# ---------------------------------------------------------------------------
# Numeric helpers — medians and the robust scales, in one place
# ---------------------------------------------------------------------------
def _number(value):
    """A CSV cell as a float, or None when it is empty or not a number."""
    text = (value or "").strip() if isinstance(value, str) else value
    if text in (None, "", "nan", "NaN", "None"):
        return None
    try:
        return float(text)
    except (TypeError, ValueError):
        return None


def median(values):
    clean = [v for v in values if v is not None]
    return statistics.median(clean) if clean else None


def mad(values, centre=None):
    """Median absolute deviation about the median (unscaled)."""
    clean = [v for v in values if v is not None]
    if not clean:
        return None
    centre = median(clean) if centre is None else centre
    return statistics.median([abs(v - centre) for v in clean])


def _quantile(values, q):
    clean = sorted(v for v in values if v is not None)
    if not clean:
        return None
    if len(clean) == 1:
        return clean[0]
    pos = (len(clean) - 1) * q
    low = int(math.floor(pos))
    high = min(low + 1, len(clean) - 1)
    return clean[low] + (clean[high] - clean[low]) * (pos - low)


def robust_scale(values):
    """MAD -> IQR/1.349 -> SD -> None, with the step that produced it.

    The chain is fixed in the pre-registration; which step was used travels with
    every cell instead of being inferred later.
    """
    clean = [v for v in values if v is not None]
    if len(clean) < 2:
        return None, "too_few"
    value = mad(clean)
    if value:
        return value, "mad"
    spread = _quantile(clean, 0.75) - _quantile(clean, 0.25)
    if spread:
        return spread / 1.349, "iqr"
    sd = statistics.pstdev(clean)
    if sd:
        return sd, "sd"
    return None, "degenerate"


def transform(site, dim, value):
    """``log`` where the pre-registration says log, identity otherwise.

    A non-positive value makes the log undefined; that unit is dropped for that
    dimension only, and the drop is counted.
    """
    if value is None:
        return None
    if TRANSFORM[site][dim] == "log":
        return math.log(value) if value > 0 else None
    return value


def percentile_ci(draws, low=2.5, high=97.5):
    clean = [v for v in draws if v is not None]
    if len(clean) < 20:
        return None
    return [_quantile(clean, low / 100.0), _quantile(clean, high / 100.0)]


# ---------------------------------------------------------------------------
# The two records as one list of units
# ---------------------------------------------------------------------------
def _is_bridge(series):
    return "Brücke" in series and "Kontrolle" not in series


def load_units(site, path=None):
    """One unit per observation: an overpass (LUMO) or a bridge/control pair.

    A unit carries its state, channel, date and the dimensions on the
    transformed scale; ``None`` marks a dimension the log made undefined, and
    the audit counts those drops. A row whose state the registry excludes is
    removed here, before any statistic, and counted under ``excluded`` — the
    registry decides the vocabulary, so a state that is neither registered nor
    excluded stops the run instead of being folded into a neighbour.
    """
    path = path or TABLE[site]
    columns = COLUMNS[site]
    audit = {"dropped": {}, "rows": 0, "pairs": 0, "unpaired": 0, "excluded": {}}
    drop = EXCLUDED.get(site) or {}
    units = []
    if site != "bautzen":
        # one row per observation (a coherence row for LUMO, a chip for Carola)
        with open(path, newline="") as handle:
            rows = list(csv.DictReader(handle))
        if site == "lumo":
            rows = [r for r in rows if r.get("record_class") == "coherence"]
        shape = UNIT_SHAPE[site]
        audit["rows"] = len(rows)
        for row in rows:
            state = row[shape["state"]]
            if state in drop:
                audit["excluded"][state] = audit["excluded"].get(state, 0) + 1
                continue
            if state not in STATES[site]:
                raise SystemExit(
                    f"{site}: state {state!r} is neither registered in STATES nor "
                    f"excluded in EXCLUDED — the registry decides the vocabulary, so "
                    f"the run stops here")
            values = {}
            for dim in DIMS:
                raw = _number(row.get(columns[dim]))
                values[dim] = transform(site, dim, raw)
                if raw is not None and values[dim] is None:
                    audit["dropped"][dim] = audit["dropped"].get(dim, 0) + 1
            units.append({"state": state,
                          "channel": CHANNEL.get(row[shape["orbit"]], "?"),
                          "date": row[shape["date"]], "values": values})
    else:
        with open(path, newline="") as handle:
            rows = list(csv.DictReader(handle))
        audit["rows"] = len(rows)
        buckets = {}
        for row in rows:
            key = (row["acquisition_date"], row["orbit_direction"])
            side = "bridge" if _is_bridge(row["series"]) else "control"
            buckets.setdefault(key, {})[side] = row
        for key in sorted(buckets):
            pair = buckets[key]
            if "bridge" not in pair or "control" not in pair:
                audit["unpaired"] += 1
                continue
            bridge, control = pair["bridge"], pair["control"]
            values = {}
            for dim in DIMS:
                left = transform(site, dim, _number(bridge.get(columns[dim])))
                right = transform(site, dim, _number(control.get(columns[dim])))
                values[dim] = None if (left is None or right is None) else left - right
                if values[dim] is None:
                    audit["dropped"][dim] = audit["dropped"].get(dim, 0) + 1
            units.append({"state": bridge["state"],
                          "channel": CHANNEL.get(bridge["orbit_direction"], "?"),
                          "date": key[0], "values": values})
        audit["pairs"] = len(units)
    units.sort(key=lambda u: (u["channel"], u["date"]))
    return units, audit


def units_of(units, channel, states):
    """The units of one channel and one or more states, in date order."""
    wanted = set(states)
    return [u for u in units if u["channel"] == channel and u["state"] in wanted]


def channels_present(units):
    """The orbit directions a record was actually acquired in, in report order.

    A record that was never acquired ASCENDING has no ASC units, so every cell,
    test and panel that would quote one is *declared* absent rather than filled
    with an empty cell — an empty cell would read as a null result.
    """
    return tuple(channel for channel in ("ASC", "DESC")
                 if any(u["channel"] == channel for u in units))


def channels_of(report, site):
    """The channels one asset's record carries, as stored in the report."""
    store = report["experiments"][site]
    return tuple(store.get("channels_present") or store["channels"])


def shared_channels(report, pair):
    """The channels both assets of a pair carry: what a cosine can be built on."""
    left, right = pair
    return tuple(channel for channel in channels_of(report, left)
                 if channel in channels_of(report, right))


def comparison_channels(report, pair):
    """The channels a pair's comparison section actually carries."""
    h2 = report["hypotheses"]["H2"]
    return tuple(channel for channel in shared_channels(report, pair)
                 if channel in h2)




# ---------------------------------------------------------------------------
# Reference, responses, magnitudes
# ---------------------------------------------------------------------------
def reference_stats(units, site, channel, dim):
    """The reference centre and the scale actually used, with both candidates.

    The 2026-10-22 amendment: when the reference state has fewer than
    ``SCALE_MIN_REF`` units, its MAD is estimated from too few points to
    standardise against, so the pooled MAD of the same asset takes over and
    ``scale_used`` records which one it was.
    """
    ref = units_of(units, channel, [REFERENCE[site]])
    values = [u["values"][dim] for u in ref if u["values"][dim] is not None]
    pooled = [u["values"][dim] for u in units
              if u["channel"] == channel and u["values"][dim] is not None]
    reference_scale, reference_kind = robust_scale(values)
    pooled_scale, pooled_kind = robust_scale(pooled)
    if len(values) >= SCALE_MIN_REF or not pooled_scale:
        scale, kind, used = reference_scale, reference_kind, "mad"
    else:
        scale, kind, used = pooled_scale, pooled_kind, "pooled_mad"
    return {"n": len(values), "median": median(values), "mad": mad(values),
            "scale": scale, "scale_kind": kind, "scale_used": used,
            "reference_scale": reference_scale, "reference_kind": reference_kind,
            "pooled_scale": pooled_scale, "pooled_kind": pooled_kind}


def _cell(values, ref):
    """One state's median, response and standardised response against ``ref``."""
    centre = median(values)
    r = None if (centre is None or ref["median"] is None) else centre - ref["median"]
    z = None if (r is None or not ref["scale"]) else r / ref["scale"]
    z_pooled = None if (r is None or not ref["pooled_scale"]) else r / ref["pooled_scale"]
    return {"n": len(values), "median": centre, "r": r, "z": z,
            "z_pooled": z_pooled, "thin": len(values) < SMALL_N}


def responses(units, site, channels=None):
    """``r`` and ``z`` for every channel, dimension, state and dose level."""
    out = {}
    for channel in (channels or channels_present(units)):
        block = {"dims": {}, "states": {}, "levels": {}}
        for dim in DIMS:
            ref = reference_stats(units, site, channel, dim)
            block["dims"][dim] = ref
            for state in STATES[site]:
                rows = units_of(units, channel, [state])
                values = [u["values"][dim] for u in rows if u["values"][dim] is not None]
                block["states"].setdefault(state, {})[dim] = _cell(values, ref)
            for name, states, dose in LEVELS[site]:
                rows = units_of(units, channel, states)
                values = [u["values"][dim] for u in rows if u["values"][dim] is not None]
                block["levels"].setdefault(name, {"dose": dose, "states": list(states),
                                                  "dims": {}})
                block["levels"][name]["dims"][dim] = _cell(values, ref)
        out[channel] = block
    return out


def magnitude(z_by_dim, dims):
    """Root mean square of the available standardised responses."""
    values = [z_by_dim[d] for d in dims if z_by_dim.get(d) is not None]
    if not values:
        return None, 0
    return math.sqrt(sum(v * v for v in values) / len(values)), len(values)


def z_of(cells, dims):
    """``{dim: z}`` from a ``{dim: cell}`` mapping."""
    return {d: (cells.get(d) or {}).get("z") for d in dims}


def unit_magnitudes(units, site, channel, dims, refs):
    """One magnitude per unit, against the fixed reference centre and scale."""
    out = []
    for unit in units_of(units, channel, STATES[site]):
        z = {}
        for dim in dims:
            ref = refs[dim]
            value = unit["values"][dim]
            z[dim] = None if (value is None or not ref["scale"]) else \
                (value - ref["median"]) / ref["scale"]
        value, k = magnitude(z, dims)
        out.append({"state": unit["state"], "R": value, "k": k})
    return out


def jonckheere(groups, n_perm=PERM, seed=SEED):
    """Jonckheere-Terpstra on ordered groups; one-sided p for an increase.

    ``groups`` is a list of value sequences in dose order. The statistic is the
    repo's verbatim port (``lumo_ocv_stats._jt_statistic``); the null permutes
    the labels and keeps the group sizes.
    """
    arrays = [np.sort(np.asarray([v for v in g if v is not None], float)) for g in groups]
    arrays = [a for a in arrays if a.size]
    if len(arrays) < 2:
        return None
    observed = _jt_statistic(arrays)
    sizes = [a.size for a in arrays]
    labels = np.concatenate([[i] * n for i, n in enumerate(sizes)])
    values = np.concatenate(arrays)
    rng = np.random.default_rng(seed)
    at_least = 0
    for _ in range(n_perm):
        shuffled = rng.permutation(values)
        draw = [np.sort(shuffled[labels == i]) for i in range(len(sizes))]
        if _jt_statistic(draw) >= observed:
            at_least += 1
    return {"statistic": float(observed), "n_groups": len(sizes),
            "n_by_group": sizes, "p_increasing": (at_least + 1) / (n_perm + 1),
            "n_perm": n_perm}



# ---------------------------------------------------------------------------
# Bootstrap — moving blocks inside a state, reference re-estimated every draw
# ---------------------------------------------------------------------------
def _resample(rows, block, rng):
    n = len(rows)
    if not n:
        return []
    block = max(1, min(block, n))
    starts = rng.integers(0, n - block + 1, int(math.ceil(n / block)))
    out = []
    for start in starts:
        out.extend(rows[int(start):int(start) + block])
    return out[:n]


def _z_vector(sample, dims, stats):
    """One state's ``z`` vector from a resample and the draw's reference stats."""
    vector = {}
    for dim in dims:
        values = [u["values"][dim] for u in sample if u["values"][dim] is not None]
        centre, scale = stats[dim]
        value = median(values)
        vector[dim] = None if (value is None or centre is None or not scale) \
            else (value - centre) / scale
    return vector


def _draw_stats(sample, dims):
    stats = {}
    for dim in dims:
        values = [u["values"][dim] for u in sample if u["values"][dim] is not None]
        stats[dim] = (median(values), robust_scale(values)[0])
    return stats


def bootstrap_magnitudes(units, site, channel, dims, boot=BOOT, seed=SEED,
                         block=BLOCK_BOOT):
    """Percentile intervals for ``R`` per level and for the top step ``dR``."""
    ref_rows = units_of(units, channel, [REFERENCE[site]])
    reference_states = {REFERENCE[site]}
    level_rows = [(name, units_of(units, channel, states), set(states) == reference_states)
                  for name, states, _dose in LEVELS[site]]
    names = [row[0] for row in level_rows]
    rng = np.random.default_rng(seed)
    draws = {name: [] for name in names}
    deltas = []
    for _ in range(boot):
        sample = _resample(ref_rows, block, rng)
        stats = _draw_stats(sample, dims)
        current = {}
        for name, rows, is_reference in level_rows:
            # The reference level is drawn from the very sample its centre and
            # scale come from, so its z is 0 by construction and its interval is
            # [0, 0] -- its sampling noise belongs to dR, where both levels share
            # one draw and it cancels (amendment 2026-10-22, item 3).
            sample_rows = sample if is_reference else _resample(rows, block, rng)
            vector = _z_vector(sample_rows, dims, stats)
            value, _k = magnitude(vector, dims)
            current[name] = value
            if value is not None:
                draws[name].append(value)
        if len(names) >= 2:
            low, high = current.get(names[-2]), current.get(names[-1])
            if low is not None and high is not None:
                deltas.append(high - low)
    return {"ci": {name: percentile_ci(values) for name, values in draws.items()},
            "by_construction": [name for name, _r, ref in level_rows if ref],
            "delta_ci": percentile_ci(deltas),
            "delta_label": f"{names[-1]}-{names[-2]}" if len(names) >= 2 else None,
            "n_draws": boot}


def bootstrap_z(units, site, channel, dims, states, boot=BOOT, seed=SEED,
                block=BLOCK_BOOT):
    """Bootstrap draws of one state's ``z`` vector (reference re-estimated)."""
    ref_rows = units_of(units, channel, [REFERENCE[site]])
    rows = units_of(units, channel, states)
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(boot):
        stats = _draw_stats(_resample(ref_rows, block, rng), dims)
        out.append(_z_vector(_resample(rows, block, rng), dims, stats))
    return out


# ---------------------------------------------------------------------------
# Direction — cosine between the two experiments' standardised responses
# ---------------------------------------------------------------------------
def cosine(left, right, dims):
    """Cosine over the dimensions defined in both vectors, and how many."""
    pairs = [(left[d], right[d]) for d in dims
             if left.get(d) is not None and right.get(d) is not None]
    if len(pairs) < 2:
        return None, 0
    dot = sum(a * b for a, b in pairs)
    norm = math.sqrt(sum(a * a for a, _b in pairs)) * math.sqrt(
        sum(b * b for _a, b in pairs))
    if not norm:
        return None, len(pairs)
    return dot / norm, len(pairs)


def sign_pattern(vector, dims):
    """``'A+ F-'`` — the direction in words, dimension by dimension."""
    parts = []
    for dim in dims:
        value = vector.get(dim)
        if value is None:
            parts.append(f"{dim}·")
        elif value == 0:
            parts.append(f"{dim}=0")
        else:
            parts.append(f"{dim}{'+' if value > 0 else '-'}")
    return " ".join(parts)


def sign_agreement(left, right, dims):
    """How many dimensions point the same way, and how many were comparable."""
    agree = total = 0
    for dim in dims:
        a, b = left.get(dim), right.get(dim)
        if a is None or b is None or a == 0 or b == 0:
            continue
        total += 1
        if (a > 0) == (b > 0):
            agree += 1
    return agree, total


def _random_abs_cosines(k, n=PERM, seed=SEED):
    rng = np.random.default_rng(seed)
    for _ in range(n):
        a = rng.normal(size=k)
        b = rng.normal(size=k)
        yield abs(float(a @ b) / (np.linalg.norm(a) * np.linalg.norm(b)))


def null_percentile(observed, k, n=PERM, seed=SEED):
    """Share of random |cos| at least as large as the observed value."""
    count = sum(1 for value in _random_abs_cosines(k, n, seed) if value >= observed)
    return (count + 1) / (n + 1)



# ---------------------------------------------------------------------------
# Holm and the report
# ---------------------------------------------------------------------------
def holm(pvalues):
    """Holm-adjusted p per label (monotone, capped at 1)."""
    items = sorted(((k, v) for k, v in pvalues.items() if v is not None),
                   key=lambda kv: kv[1])
    total = len(items)
    out, running = {}, 0.0
    for index, (key, value) in enumerate(items):
        running = max(running, (total - index) * value)
        out[key] = min(1.0, running)
    return out


def _sha256(path):
    import hashlib
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _block_cell(z, dims):
    value, k = magnitude(z, dims)
    return {"value": value, "k": k, "dims": list(dims),
            "direction": sign_pattern(z, dims), "z": z}


def _vector(report, site, channel, kind, key, dims):
    cell = report["experiments"][site]["channels"][channel]
    source = cell["states"][key]["dims"] if kind == "state" else \
        cell["levels"][key]["dims"]
    return z_of(source, dims)


def compute(sites=None, boot=BOOT, perm=PERM, seed=SEED):
    """Every number the table, the markdown and the figure quote.

    ``sites`` selects the assets: the dose pair by default (LUMO and Bautzen), the
    bridge pair with ``("carola", "morandi")``, or a single asset such as
    ``("carola",)``. A registered pair carries H2; a single asset has no partner
    and says so.
    """
    sites = tuple(sites or DEFAULT_SITES)
    for site in sites:
        if site not in TABLE:
            raise SystemExit(f"unknown site: {site}")
    with open(FREQUENCIES) as handle:
        frequencies = json.load(handle)
    report = {
        "schema": SCHEMA, "sites": list(sites),
        "generated_by": "code/sds/ocv_reconfigure.py",
        "spec": "code/sds/ocv_reconfigure_spec.md",
        "rules": {"boot": boot, "perm": perm, "seed": seed,
                  "block_boot": BLOCK_BOOT, "small_n": SMALL_N,
                  "reference": dict(REFERENCE),
                  "blocks": {k: list(v) for k, v in BLOCKS.items()},
                  "transform": {k: dict(v) for k, v in TRANSFORM.items()},
                  "levels": {k: [[n, list(st), d] for n, st, d in v]
                             for k, v in LEVELS.items()}},
        "dose": {"applies": "lumo" in sites,
                 "source": os.path.relpath(FREQUENCIES, REPO),
                 "sha256": _sha256(FREQUENCIES),
                 "definition": frequencies.get("dose_definition"),
                 "levels": frequencies.get("levels"),
                 "per_state_pct": frequencies.get("dose")},
        "experiments": {},
        "hypotheses": {"H1": {}, "H2": {}},
    }
    units_by_site = {}
    for site in sites:
        units, audit = load_units(site)
        units_by_site[site] = units
        block = responses(units, site)
        channels = channels_present(units)
        absent = {channel: ("no %s acquisition"
                            % ("ASCENDING" if channel == "ASC" else "DESCENDING"))
                  for channel in ("ASC", "DESC") if channel not in channels}
        store = {"table": os.path.relpath(TABLE[site], REPO),
                 "sha256": _sha256(TABLE[site]), "audit": audit,
                 "channels_present": list(channels), "channels_absent": absent,
                 "channels": {}}
        for channel in channels:
            cell = {"n_by_state": {s: len(units_of(units, channel, [s]))
                                   for s in STATES[site]},
                    "dims": block[channel]["dims"], "states": {}, "levels": {},
                    "blocks": {}}
            for state in STATES[site]:
                cells = block[channel]["states"][state]
                z = z_of(cells, DIMS)
                cell["states"][state] = {
                    "n_units": len(units_of(units, channel, [state])), "dims": cells,
                    "blocks": {name: _block_cell(z, dims)
                               for name, dims in BLOCKS.items()}}
            for name, states, dose in LEVELS[site]:
                cells = block[channel]["levels"][name]["dims"]
                z = z_of(cells, DIMS)
                cell["levels"][name] = {
                    "states": list(states), "dose": dose,
                    "n_units": len(units_of(units, channel, states)), "dims": cells,
                    "blocks": {b: _block_cell(z, d) for b, d in BLOCKS.items()}}
            for name, dims in BLOCKS.items():
                refs = {d: block[channel]["dims"][d] for d in dims}
                per_unit = unit_magnitudes(units, site, channel, dims, refs)
                groups = [[x["R"] for x in per_unit if x["state"] in states]
                          for _n, states, _d in LEVELS[site]]
                cell["blocks"][name] = {
                    "dims": list(dims), "n_by_level": [len(g) for g in groups],
                    "jonckheere": jonckheere(groups, n_perm=perm, seed=seed),
                    "bootstrap": bootstrap_magnitudes(units, site, channel, dims,
                                                      boot=boot, seed=seed)}
            store["channels"][channel] = cell
        report["experiments"][site] = store
    _hypotheses(report, units_by_site, boot=boot, seed=seed)
    report["table"] = table_rows(report)
    report["bootstrap_caveats"] = bootstrap_caveats(report)
    report["reading"] = reading(report)
    report["predictions"] = predictions(report)
    return report


def _partner_hint(sites):
    """Which registered pair contains these assets, for a skipped-H2 note."""
    wanted = set(sites)
    for pair in PAIRS:
        if wanted <= set(pair):
            others = [site for site in pair if site not in wanted]
            if len(others) == 1:
                return ("`%s` is its partner in the registered pair `--sites %s`, "
                        "and a cosine needs both" % (EXPERIMENT[others[0]],
                                                     ",".join(pair)))
    return ("no registered pair contains these assets, so there is nothing to "
            "compare them with")


def _hypotheses_two_state(h2, report, units_by_site, pair, boot=BOOT, seed=SEED):
    """H2 for a pair of assets that each have exactly two states.

    Only ``largest`` exists — both assets have a single response state — and only
    the ``post|post`` (or, for the towers, ``post|rebuild``) cell of the 2 x 2
    carries a response, so no ``middle`` and no full ``matrix`` are computed. The
    cosine stays descriptive, as in the dose pair: it asks whether the two assets'
    responses point the same way, and it never compares a magnitude.

    A cosine needs a vector on each side, so it runs on the channels *both* records
    carry. KDLO was acquired DESCENDING only, which leaves one channel, two cells
    and a two-test Holm family; the ASCENDING channel is declared under
    ``not_computed`` instead of being filled with an empty cell.
    """
    left_site, right_site = pair
    h2["rules"] = PAIR_NOTE[pair]["h2_rules"]
    h2["pair"] = list(pair)
    h2["not_computed"] = {
        "middle": "both assets have exactly two states, so there is no "
                  "intermediate state to compare",
        "matrix": "2 x 2, of which only the "
                  f"{LEVELS[left_site][-1][0]}|{LEVELS[right_site][-1][0]} cell "
                  "carries a response"}
    shared = shared_channels(report, pair)
    h2["shared_channels"] = list(shared)
    for channel in ("ASC", "DESC"):
        if channel not in shared:
            h2["not_computed"][channel] = (
                "no %s acquisition in both records, so no cosine is carried on "
                "that channel" % ("ASCENDING" if channel == "ASC" else "DESCENDING"))
    for channel in shared:
        h2[channel] = {}
        for block_name in ("behaviour", "geometry"):
            dims = BLOCKS[block_name]
            left = _vector(report, left_site, channel, "state",
                           LARGEST[left_site], dims)
            right = _vector(report, right_site, channel, "state",
                            LARGEST[right_site], dims)
            value, k = cosine(left, right, dims)
            draws_left = bootstrap_z(units_by_site[left_site], left_site, channel,
                                     dims, (LARGEST[left_site],),
                                     boot=boot, seed=seed)
            draws_right = bootstrap_z(units_by_site[right_site], right_site, channel,
                                      dims, (LARGEST[right_site],),
                                      boot=boot, seed=seed)
            draws = [item for item, _k in
                     (cosine(a, b, dims) for a, b in zip(draws_left, draws_right))
                     if item is not None]
            agree, total = sign_agreement(left, right, dims)
            cell = {left_site: {"kind": "state", "key": LARGEST[left_site]},
                    right_site: {"kind": "state", "key": LARGEST[right_site]},
                    "dims": list(dims), "k": k, "cosine": value,
                    "ci": percentile_ci(draws), "sign_agreement": [agree, total],
                    left_site + "_direction": sign_pattern(left, dims),
                    right_site + "_direction": sign_pattern(right, dims),
                    left_site + "_z": left, right_site + "_z": right}
            if value is not None:
                cell["p_null"] = null_percentile(abs(value), k, seed=seed)
            if block_name == "behaviour":
                # a cosine needs a vector; a single dimension is a sign comparison
                for dim in BLOCKS["behaviour"]:
                    one_left = _vector(report, left_site, channel, "state",
                                       LARGEST[left_site], (dim,))
                    one_right = _vector(report, right_site, channel, "state",
                                        LARGEST[right_site], (dim,))
                    cell["only_" + dim] = {
                        left_site + "_z": one_left.get(dim),
                        right_site + "_z": one_right.get(dim),
                        "same_sign": (None if one_left.get(dim) is None or
                                      one_right.get(dim) is None
                                      else (one_left[dim] > 0) ==
                                      (one_right[dim] > 0))}
            h2[channel][block_name] = {"largest": cell}



def _hypotheses(report, units_by_site, boot=BOOT, seed=SEED):
    """H1 (Holm over the eight primary cells) and H2 (cosines, with intervals)."""
    raw = {}
    for site, store in report["experiments"].items():
        for channel, cell in store["channels"].items():
            for name in ("behaviour", "geometry"):
                test = cell["blocks"][name].get("jonckheere") or {}
                raw[f"{site}/{channel}/{name}"] = test.get("p_increasing")
    # Holm per asset (amendment 2026-10-22, item 1): four tests each, because the
    # two assets no longer carry a joint claim.
    sites = tuple(report.get("sites") or DEFAULT_SITES)
    adjusted = {}
    for site in sites:
        family = {k: v for k, v in raw.items() if k.split("/")[0] == site}
        for key, value in holm(family).items():
            adjusted[key] = value
    report["hypotheses"]["H1_family"] = {
        site: sorted(k for k in raw if k.split("/")[0] == site)
        for site in sites}
    for site, store in report["experiments"].items():
        for channel, cell in store["channels"].items():
            for name in ("behaviour", "geometry"):
                key = f"{site}/{channel}/{name}"
                entry = cell["blocks"][name]
                entry["p_holm"] = adjusted.get(key)
                test = entry.get("jonckheere") or {}
                report["hypotheses"]["H1"][key] = {
                    "site": site, "channel": channel, "block": name,
                    "dims": entry["dims"],
                    "p_increasing": test.get("p_increasing"),
                    "p_holm": adjusted.get(key), "statistic": test.get("statistic"),
                    "R_by_level": {level: cell["levels"][level]["blocks"][name]["value"]
                                   for level in cell["levels"]},
                    "R_ci_by_level": entry["bootstrap"]["ci"],
                    "delta_ci": entry["bootstrap"]["delta_ci"],
                    "delta_label": entry["bootstrap"]["delta_label"],
                    "n_by_level": entry["n_by_level"]}

    h2 = report["hypotheses"]["H2"]
    pair = pair_of(sites)
    if pair is None:
        h2["rules"] = ("not computed: a single-asset run has no comparison "
                       "partner, so no cosine is carried — "
                       + _partner_hint(sites))
        h2["skipped"] = True
        return
    if pair in TWO_STATE_PAIRS:
        _hypotheses_two_state(h2, report, units_by_site, pair, boot=boot, seed=seed)
        return
    h2["rules"] = ("cosine of the standardised response vectors, each against its own "
                   "reference; magnitudes are never compared across experiments")
    comparisons = {"largest": {"lumo": ("state", LARGEST["lumo"], (LARGEST["lumo"],)),
                               "bautzen": ("state", LARGEST["bautzen"],
                                           (LARGEST["bautzen"],))},
                   "middle": {"lumo": ("level", "L1", tuple(MIDDLE["lumo"])),
                              "bautzen": ("level", "DS1", tuple(MIDDLE["bautzen"]))}}
    for channel in ("ASC", "DESC"):
        h2[channel] = {}
        for block_name in ("behaviour", "geometry"):
            dims = BLOCKS[block_name]
            cell = {}
            for label, sides in comparisons.items():
                left = _vector(report, "lumo", channel, sides["lumo"][0],
                               sides["lumo"][1], dims)
                right = _vector(report, "bautzen", channel, sides["bautzen"][0],
                                sides["bautzen"][1], dims)
                value, k = cosine(left, right, dims)
                draws_left = bootstrap_z(units_by_site["lumo"], "lumo", channel, dims,
                                         sides["lumo"][2], boot=boot, seed=seed)
                draws_right = bootstrap_z(units_by_site["bautzen"], "bautzen", channel,
                                          dims, sides["bautzen"][2], boot=boot, seed=seed)
                draws = [item for item, _k in
                         (cosine(a, b, dims) for a, b in zip(draws_left, draws_right))
                         if item is not None]
                agree, total = sign_agreement(left, right, dims)
                cell[label] = {
                    "lumo": {"kind": sides["lumo"][0], "key": sides["lumo"][1]},
                    "bautzen": {"kind": sides["bautzen"][0],
                                "key": sides["bautzen"][1]},
                    "dims": list(dims), "k": k, "cosine": value,
                    "ci": percentile_ci(draws), "sign_agreement": [agree, total],
                    "lumo_direction": sign_pattern(left, dims),
                    "bautzen_direction": sign_pattern(right, dims),
                    "lumo_z": left, "bautzen_z": right}
                if value is not None:
                    cell[label]["p_null"] = null_percentile(abs(value), k, seed=seed)
            if block_name == "behaviour":
                # the two dimensions of the block, side by side: a cosine needs a
                # vector, so a single dimension is a sign comparison
                for dim in BLOCKS["behaviour"]:
                    left = _vector(report, "lumo", channel, "state",
                                   LARGEST["lumo"], (dim,))
                    right = _vector(report, "bautzen", channel, "state",
                                    LARGEST["bautzen"], (dim,))
                    cell[f"only_{dim}"] = {
                        "lumo_z": left.get(dim), "bautzen_z": right.get(dim),
                        "same_sign": (None if left.get(dim) is None or
                                      right.get(dim) is None
                                      else (left[dim] > 0) == (right[dim] > 0))}
            matrix = {}
            for lumo_state in STATES["lumo"]:
                for bautzen_state in STATES["bautzen"]:
                    left = _vector(report, "lumo", channel, "state", lumo_state, dims)
                    right = _vector(report, "bautzen", channel, "state",
                                    bautzen_state, dims)
                    value, k = cosine(left, right, dims)
                    agree, total = sign_agreement(left, right, dims)
                    matrix[f"{lumo_state}|{bautzen_state}"] = {
                        "cosine": value, "k": k, "sign_agreement": [agree, total]}
            cell["matrix"] = matrix
            h2[channel][block_name] = cell


# ---------------------------------------------------------------------------
# The table, the reading, the predictions
# ---------------------------------------------------------------------------
EXPERIMENT = {"lumo": "LUMO", "bautzen": "Bautzen",
              "carola": "Carolabrücke", "morandi": "Ponte Morandi",
              "haltern": "Haltern", "kdlo": "KDLO"}

#: What one row of an asset's table *is*, for the sentences that name it.
TABLE_NOUN = {"lumo": "coherence", "bautzen": "bridge", "carola": "bridge",
              "morandi": "bridge", "haltern": "tower", "kdlo": "tower"}

#: The prose a registered pair owns, so that no sentence has to be typed twice and
#: nothing about one pair is asserted for another. The towers are read across a
#: collapse *and a rebuild*, and only the second asset has a rebuilt state.
PAIR_NOTE = {
    BRIDGE_PAIR: {
        "headline": "across a collapse",
        "noun": "bridges",
        "unit": "chips",
        "amendment": 4,
        "h2_rules": "descriptive: cosine of the two bridges' standardised "
                    "post-collapse response vectors, each against its own "
                    "pre-collapse reference; no magnitude is compared and the two "
                    "assets share no orbit, date or atmosphere",
        "window": "Both chips were masked by the same rule (`30 %` of the peak and "
                  "`5 ×` the ring median, 200 m ring, 30 m guard) in "
                  "`erlking_mvp/tools/ocv_table.py`; the comparison inherits that "
                  "rule rather than re-deriving it.",
        "event": "",
    },
    TOWER_PAIR: {
        "headline": "across a collapse and a rebuild",
        "noun": "towers",
        "unit": "overpasses",
        "amendment": 5,
        "h2_rules": "descriptive: cosine of the two towers' standardised response "
                    "vectors at the collapsed state — Haltern's post-collapse, "
                    "KDLO's rebuild — each against its own pre-collapse reference; "
                    "no magnitude is compared and the two assets share no orbit, "
                    "date or atmosphere",
        "window": "Both chips were masked by the same rule (`30 %` of the peak and "
                  "`5 ×` the ring median, 200 m ring, 30 m guard) in "
                  "`erlking_mvp/tools/ocv_table.py`, the same rule the two bridges "
                  "were cut with; the comparison inherits that rule rather than "
                  "re-deriving it.",
        "event": " — Haltern's turbine came down on 2021-09-29, and KDLO's mast "
                 "collapsed between 2022-12-11 and 2022-12-23 and was rebuilt "
                 "afterwards (its second state is a *rebuild*, not a ruin)",
    },
}


def bootstrap_caveats(report):
    """Levels whose ``n`` does not exceed the bootstrap block length.

    A moving-block resample of a cell with ``n <= BLOCK_BOOT`` draws that cell as a
    single block on every iteration, so that level's interval carries the
    *reference's* uncertainty only — the cell itself contributes a fixed number.
    The list is recorded in the report rather than tolerated silently; the refusal
    in the spec ("no thin cell quoted without its n") is what it serves.
    """
    cells = []
    for site in (report.get("sites") or DEFAULT_SITES):
        for channel in channels_of(report, site):
            blocks = report["experiments"][site]["channels"][channel]["blocks"]
            for name in ("behaviour", "geometry"):
                sizes = blocks[name].get("n_by_level") or []
                for index, (level, _states, _dose) in enumerate(LEVELS[site]):
                    count = sizes[index] if index < len(sizes) else None
                    if count is not None and count <= BLOCK_BOOT:
                        cells.append({"site": site, "channel": channel,
                                      "block": name, "level": level,
                                      "n": count, "block_covers_cell": True})
    return {"block": BLOCK_BOOT, "small_n": SMALL_N, "cells": cells}


def table_rows(report):
    """The compact table the pre-registration asks for, one row per cell."""
    rows = []
    for site in (report.get("sites") or DEFAULT_SITES):
        store = report["experiments"][site]
        for channel in channels_of(report, site):
            cell = store["channels"][channel]
            for state in STATES[site]:
                entry = cell["states"][state]
                echo, geometry = (entry["blocks"]["behaviour"],
                                  entry["blocks"]["geometry"])
                rows.append({
                    "Experiment": EXPERIMENT[site], "State": state,
                    "Channel": channel, "n": entry["n_units"],
                    "R_behavior": echo["value"], "R_geometry": geometry["value"],
                    "Direction": (f"{echo['direction']} | "
                                  f"{geometry['direction']}"),
                    "Scale": cell["dims"]["A"]["scale_used"],
                    "Support": "thin" if entry["n_units"] < SMALL_N else "ok"})
    return rows


def reading_by_site(report):
    """One verdict per asset and per block, with the bound each interval licenses."""
    out = {}
    h1 = report["hypotheses"]["H1"]
    for site in (report.get("sites") or DEFAULT_SITES):
        channels = channels_of(report, site)
        block = {name: {ch: h1[f"{site}/{ch}/{name}"] for ch in channels}
                 for name in ("behaviour", "geometry")}
        summary = {}
        for name, cells in block.items():
            significant = [ch for ch, c in cells.items()
                           if (c["p_holm"] or 1.0) < 0.05]
            monotone = {}
            bounds = {}
            for channel in channels:
                values = [v for v in cells[channel]["R_by_level"].values()
                          if v is not None]
                monotone[channel] = (all(b >= a for a, b in zip(values, values[1:]))
                                     if len(values) > 1 else None)
                pair = cells[channel]["delta_ci"]
                bounds[channel] = {"delta_ci": pair,
                                   "excludes_increase_above": (None if not pair
                                                               else pair[1]),
                                   "excludes_decrease_below": (None if not pair
                                                               else pair[0])}
            ordered = len(LEVELS[site]) > 2
            if significant:
                what = ("a dose ordering" if ordered else
                        f"the {CONTRAST[site]} contrast")
                verdict = (f"{what} in the {BLOCK_DISPLAY[name]} block survives the "
                           f"asset's own Holm family on " + ", ".join(significant))
            elif any(monotone.values()):
                rising = [ch for ch, flag in monotone.items() if flag]
                where = ("both channels" if len(rising) == len(channels) == 2
                         else " ".join(rising) + " only")
                verdict = (f"no {('dose ordering' if ordered else 'state contrast')} "
                           f"in the {BLOCK_DISPLAY[name]} block survives the asset's "
                           f"own Holm family; "
                           f"the point estimate of the step is positive on " + where)
            else:
                verdict = (f"no {('dose ordering' if ordered else 'state contrast')} "
                           f"in the {BLOCK_DISPLAY[name]} block survives the asset's "
                           f"own Holm family, "
                           f"and the point estimate of the step is not monotone")
                verdict += (" on either channel" if len(channels) > 1
                            else " on the one channel this record carries")
            summary[name] = {
                "verdict": verdict, "significant": significant,
                "monotone_by_channel": monotone, "delta_ci": bounds,
                "widest_upper_bound": max((c["delta_ci"] or [None, None])[1] or 0.0
                                          for c in cells.values())}
        verdicts = {name: item["verdict"] for name, item in summary.items()}
        out[site] = {
            "verdicts": verdicts, "blocks": summary,
            "n_ref": {ch: report["experiments"][site]["channels"][ch]["dims"]["A"]["n"]
                      for ch in channels},
            "scale_used": {ch: report["experiments"][site]["channels"][ch]["dims"]["A"]
                                       ["scale_used"] for ch in channels}}
    return out


def reading(report):
    """Per-asset verdicts plus — for the pair — the descriptive comparison."""
    h2 = report["hypotheses"]["H2"]
    pair = pair_of(report.get("sites") or DEFAULT_SITES)
    left_site, right_site = pair or (None, None)
    comparison = {}
    if not h2.get("skipped"):
        for channel in comparison_channels(report, pair):
            comparison[channel] = {}
            for block_name in ("behaviour", "geometry"):
                item = h2[channel][block_name]["largest"]
                comparison[channel][block_name] = {
                    "cosine": item["cosine"], "ci": item["ci"],
                    "sign_agreement": item["sign_agreement"],
                    "p_null": item.get("p_null"),
                    left_site + "_direction": item[left_site + "_direction"],
                    right_site + "_direction": item[right_site + "_direction"]}
    role = ("not computed: a single-asset run has no comparison partner"
            if h2.get("skipped") else
            "descriptive: the same metric applied to two assets, not a "
            "transferability test")
    return {"per_asset": reading_by_site(report), "comparison": comparison,
            "comparison_role": role}


def predictions(report):
    """P1-P3, re-stated per block in the 2026-10-22 amendment (block split).

    The three predictions belong to LUMO (they are statements about the measured
    frequency dose), so a run without the tower has none.
    """
    out = []
    if "lumo" not in (report.get("sites") or DEFAULT_SITES):
        return out

    # P1 -- the geometry block of the tower rises with the measured dose
    details, holds = [], []
    for channel in ("ASC", "DESC"):
        levels = report["hypotheses"]["H1"][f"lumo/{channel}/geometry"]["R_by_level"]
        low, mid, high = levels.get("L0"), levels.get("L1"), levels.get("L2")
        holds.append(all(v is not None for v in (low, mid, high)) and high > mid > 0)
        details.append(f"{channel} {_num(low)} < {_num(mid)} < {_num(high)}")
    out.append({"key": "P1",
                "verdict": "held" if all(holds) else
                           ("partly" if any(holds) else "did not hold"),
                "statement": "the tower's **geometry** magnitude rises with the dose, "
                             "L0 < L1 < L2",
                "detail": "; ".join(details)})

    # P2 -- DAM 3 and DAM 4 share one dose, so the same-dose gap stays inside the
    # cross-dose step. Evaluated per block, because the amendment made the two
    # blocks answer differently.
    details, verdicts = [], []
    for name in ("geometry", "behaviour"):
        holds = []
        for channel in ("ASC", "DESC"):
            states = report["experiments"]["lumo"]["channels"][channel]["states"]
            three = states["DAM 3"]["blocks"][name]["value"]
            four = states["DAM 4"]["blocks"][name]["value"]
            levels = report["hypotheses"]["H1"][f"lumo/{channel}/{name}"]["R_by_level"]
            step = (None if levels.get("L1") is None or levels.get("L2") is None
                    else abs(levels["L2"] - levels["L1"]))
            gap = None if (three is None or four is None) else abs(three - four)
            holds.append(gap is not None and step is not None and gap <= step)
            details.append(f"{BLOCK_DISPLAY[name]}/{channel}: "
                           f"|R(DAM 3) - R(DAM 4)| = {_num(gap)} "
                           f"against a cross-dose step of {_num(step)}")
        verdicts.append("held" if all(holds) else
                        ("partly" if any(holds) else "did not hold"))
    out.append({"key": "P2",
                "verdict": "; ".join(f"{BLOCK_DISPLAY[name]} {verdict}"
                                     for name, verdict in
                                     zip(("geometry", "behaviour"), verdicts))
                           if len(set(verdicts)) > 1 else verdicts[0],
                "statement": "DAM 3 and DAM 4 share one dose, so their magnitudes "
                             "agree (the same-dose gap is not larger than the "
                             "cross-dose step) — asked per block",
                "detail": "; ".join(details)})

    # P3 -- the behavior block need not follow the dose
    behavior = [report["hypotheses"]["H1"][f"{site}/{channel}/behaviour"]["p_holm"]
                for site in (report.get("sites") or DEFAULT_SITES)
                for channel in ("ASC", "DESC")]
    following = [p for p in behavior if p is not None and p < 0.05]
    out.append({"key": "P3",
                "verdict": "consistent" if not following else
                           "the behavior block follows the dose too",
                "statement": "the γ²/P behavior block need not follow the dose",
                "detail": f"{len(following)} of {len(behavior)} behavior-block tests "
                          f"reach Holm < 0.05"})
    return out


def _num(value, digits=3):
    return "–" if value is None else f"{value:.{digits}g}"


def _ci(pair, digits=3):
    if not pair:
        return "–"
    return f"[{pair[0]:.{digits}g}, {pair[1]:.{digits}g}]"


def _p(value):
    if value is None:
        return "–"
    return "< 0.001" if value < 0.001 else f"{value:.3f}"



def _caveat_sentence(report):
    """One sentence naming the cells whose interval cannot resolve their own cell.

    An empty caveat list produces no sentence at all: a plate with no such cell
    must not carry a sentence saying so.
    """
    caveats = report["bootstrap_caveats"]
    labels = []
    for cell in caveats["cells"]:
        label = (f"{EXPERIMENT[cell['site']]} {cell['channel']} "
                 f"`{cell['level']}` (n = {cell['n']})")
        if label not in labels:
            labels.append(label)
    if not labels:
        return ""
    return ("**Bootstrap caveat:** the moving-block resample runs on blocks of "
            f"{caveats['block']} units, so a cell with `n` at or below that is drawn "
            "as a single block and its interval reflects the reference's uncertainty "
            "only — here: " + "; ".join(labels) + ".")


#: One clause per asset for the opening paragraph; the pair text is these clauses
#: joined, so it is generated rather than typed twice.
REFERENCE_NOTE = {
    "lumo": "LUMO against `healthy` on the coherence rows",
    "bautzen": "Bautzen against `RS` on the paired bridge-minus-control difference",
    "carola": "Carola against `pre-collapse (healthy)` on the support-sized chip, "
              "with the post-collapse state as the single contrast",
    "morandi": "Morandi against `pre-collapse (healthy)` on the support-sized chip, "
               "with the post-collapse state as the single contrast",
    "haltern": "Haltern against `pre-collapse (healthy)` on the support-sized chip, "
               "with the post-collapse state as the single contrast",
    "kdlo": "KDLO against `pre` on the support-sized chip, with the rebuilt state as "
            "the single contrast and one orbit direction only",
}


#: Small counts, spelled out: "of four tests", not "of 4 tests".
_WORDS = {1: "one", 2: "two", 3: "three", 4: "four", 5: "five", 6: "six"}


def markdown(report):
    """Per asset one verdict; across a registered pair one descriptive comparison."""
    lines = []
    add = lines.append
    dose = report["dose"]
    reading = report["reading"]
    sites = tuple(report.get("sites") or DEFAULT_SITES)
    pair = pair_of(sites)
    left_site, right_site = pair or (None, None)
    left_name = EXPERIMENT[left_site] if left_site else ""
    right_name = EXPERIMENT[right_site] if right_site else ""
    if pair == DEFAULT_SITES:
        add("# OCV reconfiguration across controlled states — two assets, "
            "two verdicts")
    elif pair:
        add(f"# OCV reconfiguration {PAIR_NOTE[pair]['headline']} — {left_name} vs. "
            f"{right_name}, two verdicts")
    elif len(STATES[sites[0]]) == 2:
        add(f"# OCV reconfiguration {SITE_HEADLINE.get(sites[0], 'across a collapse')} "
            f"— {EXPERIMENT[sites[0]]}, two states")
    else:
        add(f"# OCV reconfiguration across the controlled states — "
            f"{EXPERIMENT[sites[0]]}")
    add("")
    add("*generated by `code/sds/ocv_reconfigure.py`; the design is fixed in "
        "`code/sds/ocv_reconfigure_spec.md`, including the 2026-10-22 amendment that "
        "separates the two assets"
        + (f" and amendment {PAIR_NOTE[pair]['amendment']} that registers this pair*"
           if pair and pair != DEFAULT_SITES else "*"))
    add("")
    families = {site: 2 * len(channels_of(report, site)) for site in sites}
    if len(set(families.values())) == 1:
        family = f"of {_WORDS.get(families[sites[0]], families[sites[0]])} tests"
    else:
        family = ("— two cells per channel, so "
                  + " and ".join(f"{_WORDS.get(count, count)} for {EXPERIMENT[site]}"
                                 for site, count in families.items()) + " —")
    add("Each asset is measured against **its own reference state**, in units of that "
        "reference's own variability"
        + " — " + ", ".join(REFERENCE_NOTE[site] for site in sites)
        + ". Assets share a metric, not a scale, so their magnitudes are never "
        "compared. Each carries its **own** Holm family " + family + " and its own "
        "verdict; no joint label is issued"
        + (", and the comparison section is a description, not a hypothesis."
           if pair else "."))
    add("")
    if "lumo" in sites:
        add("## What the dose is")
        add("")
        add(f"LUMO's dose is measured (`{dose['source']}`, sha256 "
            f"`{str(dose['sha256'])[:12]}`): {dose['definition']}")
        add("")
        add("| level | states | dose (frequency drop) |")
        add("| --- | --- | ---: |")
        for level in dose["levels"] or []:
            add(f"| {level['level']} | {' / '.join(level['states'])} | "
                f"{level['dose_pct']:g} % |")
        add("")
        add("DAM 3 and DAM 4 sit one FFT bin apart and are therefore one dose level: the "
            "physics does not license a four-step ladder. Bautzen has **no** measured "
            "frequency per state — its ladder `RS < DS1 < DS2` is ordinal, and no dose "
            "regression is run on it.")
        add("")
    else:
        add("## There is no dose")
        add("")
        # A two-state pair whose levels are literally `pre`/`post` keeps the wording
        # the bridge amendment used. KDLO's second level is a *rebuild*, so the
        # sentence names the reference instead of a level that does not exist.
        plain_states = {LEVELS[site][-1][0] for site in sites} == {"post"}
        add(("Both assets have" if pair else "This asset has")
            + " exactly two states, separated by the collapse date"
            + (PAIR_NOTE[pair]["event"] if pair else "")
            + ", so there is no " "ladder to order and no frequency drop to regress "
            "against. The ordered "
            "test reduces to its two-group case — one one-sided rank test between "
            + ("`pre` and `post`" if plain_states else
               "the reference and the other state")
            + " per block"
            + (" and per asset — and `ΔR` is the post-minus-pre step." if pair and plain_states
               else " and per asset — and `ΔR` is the step." if pair
               else ", and `ΔR` is the post-minus-pre step." if plain_states
               else ", and `ΔR` is the step."))
        add("")
        counted = [channel for channel in ("ASC", "DESC")
                   if any(channel in channels_of(report, site) for site in sites)]
        add("| experiment | state | "
            + " | ".join(f"n {channel}" for channel in counted) + " |")
        add("| --- | --- |" + " ---: |" * len(counted))
        for site in sites:
            for index, (_level, state_names, _dose_value) in enumerate(LEVELS[site]):
                counts = []
                for channel in counted:
                    if channel not in channels_of(report, site):
                        counts.append("–")
                        continue
                    sizes = (report["experiments"][site]["channels"][channel]["blocks"]
                             ["geometry"].get("n_by_level") or [])
                    counts.append(sizes[index] if index < len(sizes) else "–")
                mark = " (reference)" if state_names == (REFERENCE[site],) else ""
                add(f"| {EXPERIMENT[site]} | {state_names[0]}{mark} | "
                    + " | ".join(str(count) for count in counts) + " |")
        add("")
        absent = [(site, channel) for site in sites
                  for channel in (report["experiments"][site].get("channels_absent")
                                  or {})]
        if absent:
            add("A `–` is a channel the record does not contain, not a cell that came "
                "out empty: "
                + "; ".join(f"{EXPERIMENT[site]} was never acquired "
                            f"{'ASCENDING' if channel == 'ASC' else 'DESCENDING'}"
                            for site, channel in absent)
                + ". Such a channel carries no cell, no test and no panel, and it "
                "contributes nothing to the asset's Holm family.")
            add("")
        add("LUMO's measured frequency ladder is not used here — "
            + ("not by either asset" if pair else "not by this asset")
            + "; it is kept under `dose` in the JSON, flagged `applies: false`, and "
            "quoted in `data/sds/ocv_reconfigure.md` where it is the dose.")
        add("")
    for site in (report.get("sites") or DEFAULT_SITES):
        block = reading["per_asset"][site]
        add(f"## {EXPERIMENT[site]}")
        add("")
        add(f"* reference `{REFERENCE[site]}`, n_ref = "
            + " / ".join(f"{block['n_ref'][channel]} ({channel})"
                         for channel in channels_of(report, site))
            + "; scale used = "
            + " / ".join(f"`{block['scale_used'][channel]}` ({channel})"
                         for channel in channels_of(report, site)) + ".")
        add("")
        add("| State | Channel | n | R_behavior | R_geometry | Direction | Scale | Support |")
        add("| --- | --- | ---: | ---: | ---: | --- | --- | --- |")
        for row in report["table"]:
            if row["Experiment"] != EXPERIMENT[site]:
                continue
            add(f"| {row['State']} | {row['Channel']} | {row['n']}"
                + ("✱" if isinstance(row["n"], int) and row["n"] <= BLOCK_BOOT else "")
                + f" | "
                f"{_num(row['R_behavior'])} | {_num(row['R_geometry'])} | "
                f"{row['Direction']} | `{row['Scale']}` | {row['Support']} |")
        add("")
        for channel, note in (report["experiments"][site].get("channels_absent")
                              or {}).items():
            add(f"*{channel} is absent from this record ({note}), so it has no rows in "
                f"the table above, no H1 cells below and no panels on the plate: "
                f"`–`, `not_computed` and an empty panel, never a null.*")
            add("")
        add("| channel | block | R by level | ΔR, top step | p (one-sided) | Holm "
            "in the asset's family |")
        add("| --- | --- | --- | --- | ---: | ---: |")
        for channel in channels_of(report, site):
            for name in ("behaviour", "geometry"):
                cell = report["hypotheses"]["H1"][f"{site}/{channel}/{name}"]
                levels = " < ".join(f"{key} {_num(value)}"
                                    for key, value in cell["R_by_level"].items())
                interval = cell["delta_ci"] or [None, None]
                add(f"| {channel} | `{BLOCK_DISPLAY[name]}` "
                    f"[{', '.join(cell['dims'])}] | {levels} | "
                    f"{_num(interval[0])} … {_num(interval[1])} | "
                    f"{_p(cell['p_increasing'])} | {_p(cell['p_holm'])} |")
        add("")
        for name in ("behaviour", "geometry"):
            item = block["blocks"][name]
            add(f"**Verdict — {EXPERIMENT[site]}, {BLOCK_DISPLAY[name]}**: "
                f"{item['verdict']}. Its own "
                f"interval still covers a top-step change of up to "
                f"{_num(item['widest_upper_bound'], 2)} R-units, so what the asset "
                f"**excludes** for this block is bounded by that number, and nothing "
                f"smaller.")
            add("")
        add("")
        if site == "lumo":
            add("The tower's two visits per state make one further statement "
                "possible, and it is the strongest single result here: "
                "**the two states at the same measured dose differ more than the two "
                "dose levels do** (P2 below), so the response is not a function of "
                "the frequency drop.")
            add("")
        elif site == "bautzen":
            add("The bridge's states are its calendar; the paired "
                "bridge-minus-control difference removes everything the two tracks "
                "share (season, acquisition day, processing) and nothing more.")
            add("")
        elif site == "carola":
            add("Carola has two states and no control track: `pre` and `post` are "
                "separated by the collapse date itself, so this is a before/after "
                "contrast in time, not a controlled comparison. The table is the "
                "support-sized chip of `Girder A` alone.")
            add("")
        elif site == "morandi":
            audit = report["experiments"]["morandi"]["audit"]
            units = audit.get("pairs") or audit.get("rows")
            sizes = {channel: (report["experiments"]["morandi"]["channels"][channel]
                               ["blocks"]["geometry"]["n_by_level"])
                     for channel in ("ASC", "DESC")}
            add(f"Morandi has two states and no control track either: `pre` and "
                f"`post` are separated by the collapse date, so this is a "
                f"before/after contrast in time. The table is the support-sized chip "
                f"({units} chips), one platform track, both orbit directions. The "
                f"post-collapse state holds {sizes['ASC'][1]} ASC and "
                f"{sizes['DESC'][1]} DESC chips against {sizes['ASC'][0]} and "
                f"{sizes['DESC'][0]} pre-collapse ones: **p = 1.000 here is not "
                f"evidence of no change — it is what this many chips can resolve.**")
            add("")
        elif site == "haltern":
            audit = report["experiments"]["haltern"]["audit"]
            sizes = {channel: (report["experiments"]["haltern"]["channels"][channel]
                               ["blocks"]["geometry"]["n_by_level"])
                     for channel in channels_of(report, site)}
            after = " and ".join(f"{sizes[channel][1]} {channel}"
                                 for channel in channels_of(report, site))
            before = " and ".join(f"{sizes[channel][0]} {channel}"
                                  for channel in channels_of(report, site))
            thin = [channel for channel in channels_of(report, site)
                    if sizes[channel][1] <= BLOCK_BOOT]
            below = [channel for channel in channels_of(report, site)
                     if sizes[channel][1] < SMALL_N]
            add(f"Haltern has two states and no control track: the turbine came down "
                f"on 2021-09-29, and `pre-collapse (healthy)` and `post-collapse` are "
                f"separated by that date alone, so this is a before/after contrast in "
                f"time, not a controlled comparison. The table is the support-sized "
                f"chip of the turbine mast ({audit.get('rows')} of the 120 committed "
                f"overpasses). The post-collapse state holds {after} units against "
                f"{before} pre-collapse ones"
                + (f"; {' and '.join(thin)} is at or below the bootstrap block of "
                   f"{BLOCK_BOOT}" if thin else "")
                + (f"; {', '.join(below)} sits below the six-unit threshold the table "
                   f"calls `thin`" if below else "")
                + (". No contrast survives the asset's own family, so the intervals "
                   "above are what these data give."
                   if not [1 for channel in channels_of(report, site)
                           for name in ("behaviour", "geometry")
                           if (report["hypotheses"]["H1"][f"{site}/{channel}/{name}"]
                               ["p_holm"] or 1.0) < 0.05] else "."))
            add("")
        elif site == "kdlo":
            audit = report["experiments"]["kdlo"]["audit"]
            channels = channels_of(report, site)
            sizes = {channel: (report["experiments"]["kdlo"]["channels"][channel]
                               ["blocks"]["geometry"]["n_by_level"])
                     for channel in channels}
            excluded = audit.get("excluded") or {}
            cell = report["experiments"]["kdlo"]["channels"][channels[0]]
            largest = cell["states"][LARGEST["kdlo"]]
            block = largest["blocks"]["geometry"]["z"]
            shares = {dim: (block.get(dim) or 0.0) ** 2 for dim in BLOCKS["geometry"]}
            leading = max(shares, key=shares.get)
            lead_share = shares[leading] / sum(shares.values()) if sum(shares.values()) else 0.0
            others = max(abs(block[dim]) for dim in BLOCKS["geometry"] if dim != leading)
            add(f"KDLO has two states and no control track either: the mast collapsed "
                f"between 2022-12-11 and 2022-12-23 and was rebuilt, so `pre` and "
                f"`rebuild` are separated by that event alone — a before/after "
                f"contrast in time in which the second state is a *new* structure, "
                f"not a damaged one. The table is the support-sized chip, which here "
                f"is the whole record ({audit.get('pairs') or audit.get('rows')} "
                f"acquisitions)"
                + (", and the record carries no ASCENDING acquisition, so this asset "
                   "has one channel, two H1 cells per block and a two-test Holm "
                   "family." if len(channels) == 1 else ".")
                + " The reference holds "
                + " and ".join(f"{sizes[channel][0]} {channel}" for channel in channels)
                + " units and the rebuilt state "
                + " and ".join(f"{sizes[channel][1]} {channel}" for channel in channels)
                + " units."
                + "".join(f" The {'one' if count == 1 else count} acquisition"
                          f"{'' if count == 1 else 's'} in state `{state}` "
                          f"{'is' if count == 1 else 'are'} excluded before any "
                          f"statistic — {EXCLUDED[site][state]} — and counted in the "
                          f"audit, so neither state inherits it."
                          for state, count in excluded.items())
                + (f" The geometry magnitude of the rebuilt state (`R` = "
                   f"{_num(largest['blocks']['geometry']['value'], 3)}) is carried by "
                   f"`{leading}` (z = {block[leading]:+.1f}, "
                   f"{lead_share * 100:.1f} % of the block's squared magnitude): the "
                   f"reference's `{leading}` median is "
                   f"{_num(cell['dims'][leading]['median'], 3)} (on the "
                   f"`{TRANSFORM['kdlo'][leading]}` scale) with a MAD of "
                   f"{_num(cell['dims'][leading]['scale'], 2)}, the rebuilt state's is "
                   f"{_num(largest['dims'][leading]['median'], 3)}, and the other "
                   f"three dimensions stay within {_num(others, 2)} z. `R` is an RMS "
                   f"over four dimensions, so one coordinate can carry it — which is "
                   f"why `{leading}` is annotated beside the plane and not drawn on "
                   f"it." if lead_share > 0.5 else ""))
            add("")
        else:
            levels = " and ".join(f"`{name}`" for name, _states, _dose in LEVELS[site])
            add(f"{EXPERIMENT[site]} has two states and no control track: {levels} "
                f"are separated by the event date alone, so this is a before/after "
                f"contrast in time, not a controlled comparison.")
            add("")
        applicable = [item for item in report["predictions"]
                      if (site == "lumo" and item["key"] in ("P1", "P2"))
                      or (site == "bautzen" and item["key"] == "P3")]
        if applicable:
            add("**Pre-registered predictions that belong to this asset:**")
            add("")
            for item in applicable:
                add(f"* **{item['key']}** — {item['statement']}: "
                    f"*{item['verdict']}*. {item['detail']}")
            add("")
    if not pair:
        add("## No comparison section")
        add("")
        add("A single-asset run carries no cosine: " + _partner_hint(sites) + ".")
        add("")
        partners = [item for item in PAIRS if set(sites) <= set(item)]
        if partners:
            add(f"The comparison is reported in `data/sds/{output_stem(partners[0])}.md`, "
                f"written by `--sites {','.join(partners[0])}`.")
            add("")
    if pair:
        add("## Comparison — the same metric, two assets (descriptive)")
        add("")
        add("Not a hypothesis and not a test: the cosine of the two standardised "
            "response vectors at the largest state says whether the *same metric*, applied "
            "twice to two assets, answers with the same direction. Magnitudes are not "
            "compared — each asset's `R` lives on its own scale.")
        add("")
        add(f"| channel | block | cosine | 95 % CI | sign agreement | null p | "
            f"{left_name} | {right_name} |")
        add("| --- | --- | ---: | --- | --- | ---: | --- | --- |")
        for channel in comparison_channels(report, pair):
            for block_name in ("behaviour", "geometry"):
                cell = report["hypotheses"]["H2"][channel][block_name]
                item = cell["largest"]
                add(f"| {channel} | `{BLOCK_DISPLAY[block_name]}` "
                    f"[{', '.join(item['dims'])}] | {_num(item['cosine'])} | "
                    f"{_ci(item['ci'])} | "
                    f"{item['sign_agreement'][0]}/{item['sign_agreement'][1]} | "
                    f"{_p(item.get('p_null'))} | {item[left_site + '_direction']} | "
                    f"{item[right_site + '_direction']} |")
                for dim in ("gamma2", "P"):
                    only = cell.get(f"only_{dim}")
                    if not only:
                        continue
                    add(f"| {channel} | `{dim}` alone | — | – | – | – | "
                        f"{_num(only[left_site + '_z'], 2)} | "
                        f"{_num(only[right_site + '_z'], 2)} |")
        add("")
        add("`null p` is the share of random directions in 3 (or 2) dimensions whose "
            "absolute cosine is at least as large — a guard against a coarse direction "
            "space, not a significance claim. "
            + ("γ² is a different quantity in the two experiments, "
               if pair == DEFAULT_SITES else
               "γ² is the same exported quantity in the two chips but measured on "
               "different structures at different dates, ")
            + "which is why the two dimensions of the behavior block are "
            "also printed one by one (`z` per asset). The geometry block carries four "
            "dimensions while its plane shows two: `R_geometry` is the RMS over "
            "`A, D, F, S`, and `D` (which partly inherits `A`) and `S` (position, not "
            "morphology) are annotated beside the plane rather than drawn as axes.")
        add("")
    if pair in TWO_STATE_PAIRS:
        add("### Two verdicts, read together")
        add("")
        h1 = report["hypotheses"]["H1"]
        resolved_sites = []
        for site in pair:
            cells = [(channel, block_name, h1[f"{site}/{channel}/{block_name}"])
                     for channel in channels_of(report, site)
                     for block_name in ("behaviour", "geometry")]
            level = LEVELS[site][-1][0]
            parts = "; ".join(
                f"{channel} {BLOCK_DISPLAY[block_name]} `R` = "
                f"{_num(cell['R_by_level'].get(level))}, "
                f"Holm {_p(cell['p_holm'])}" for channel, block_name, cell in cells)
            resolved = [f"{channel} {BLOCK_DISPLAY[block_name]}"
                        for channel, block_name, cell in cells
                        if cell["p_holm"] is not None and cell["p_holm"] < 0.05]
            if resolved:
                resolved_sites.append(site)
            add(f"* **{EXPERIMENT[site]}** — {parts}. "
                + ("Resolves " + ", ".join(resolved) + " at Holm < 0.05."
                   if resolved else "Resolves no cell at Holm < 0.05."))
        add("")
        cells = [(channel, block_name,
                  report["hypotheses"]["H2"][channel][block_name]["largest"])
                 for channel in comparison_channels(report, pair)
                 for block_name in ("behaviour", "geometry")]
        # A shared direction needs a positive cosine *and* the signs to agree
        # dimension by dimension: a cosine that one dimension carries alone is not a
        # direction the two assets share.
        shared = [item for item in cells
                  if (item[2]["cosine"] or 0) > 0
                  and item[2]["sign_agreement"][0] == item[2]["sign_agreement"][1]]
        if shared:
            phrases = []
            for channel, block_name, item in shared:
                low, high = item["ci"] or (None, None)
                covers = (low is None or high is None or low <= 0.0 <= high)
                phrases.append(
                    f"{channel} {BLOCK_DISPLAY[block_name]}, cosine {_num(item['cosine'])} "
                    f"{_ci(item['ci'])}, {item['sign_agreement'][0]}/"
                    f"{item['sign_agreement'][1]} signs equal "
                    f"(`{item[left_site + '_direction']}` in {left_name}, "
                    f"`{item[right_site + '_direction']}` in {right_name}), null p "
                    f"{_p(item.get('p_null'))}"
                    + (" — an interval that still covers zero, so the agreement is in "
                       "the signs and not in the size" if covers else
                       " — an interval that excludes zero"))
            add(f"Of the {len(cells)} cells, "
                + ("exactly one shares a direction: " if len(shared) == 1
                   else f"{len(shared)} share a direction: ")
                + "; ".join(phrases) + ".")
            add("")
        elif cells:
            channel, block_name, item = max(
                cells, key=lambda cell: (cell[2]["sign_agreement"][0]
                                         / max(1, cell[2]["sign_agreement"][1])))
            agree, total = item["sign_agreement"]
            add(f"None of the {len(cells)} cells shares a direction: in every one the "
                f"cosine is not positive or the signs disagree in at least one "
                f"dimension. The closest is {channel} {BLOCK_DISPLAY[block_name]}, cosine "
                f"{_num(item['cosine'])} {_ci(item['ci'])} with {agree}/{total} signs "
                f"equal (null p {_p(item.get('p_null'))}) — a cosine that does not "
                f"carry the block's dimensions with it.")
            add("")
        sizes = {channel: report["experiments"][right_site]["channels"][channel]
                 ["blocks"]["geometry"]["n_by_level"]
                 for channel in channels_of(report, right_site)}
        post = sum(sum(sizes[channel][1:]) for channel in sizes)
        pre = sum(sum(sizes[channel][:1]) for channel in sizes)
        unit = PAIR_NOTE[pair]["unit"]
        lead = ("both reorganised measurably, each against its own reference and on "
                "its own scale" if len(resolved_sites) == 2 else
                "one of them reorganised measurably and the other did not"
                if resolved_sites else
                "neither reorganised measurably at this n")
        add(f"So the comparison does not say that the two {PAIR_NOTE[pair]['noun']} "
            f"reorganised alike. It says that {lead}"
            + (", and that where the two directions do agree, the agreement rests on "
               "signs rather than on magnitudes, on a null that is too coarse to call "
               f"significant, and on {post} {unit} after the event against {pre} "
               "before it in the second asset."
               if shared else
               ", and that no cell has both a positive cosine and agreeing signs — "
               f"the second asset holds {post} {unit} after the event against {pre} "
               "before it.")
            + " That is the answer these data support, and it is not the same as *no "
              "change*.")
        add("")
    add("### Units dropped, per dimension")
    add("")
    add("| experiment | units | " + " | ".join(DIMS) + " | unpaired |")
    add("| --- | ---: |" + " ---: |" * (len(DIMS) + 1))
    for site in (report.get("sites") or DEFAULT_SITES):
        audit = report["experiments"][site]["audit"]
        dropped = audit.get("dropped") or {}
        cells = " | ".join(str(dropped.get(dim, 0)) for dim in DIMS)
        add(f"| {EXPERIMENT[site]} | {audit.get('pairs') or audit.get('rows')} | "
            f"{cells} | {audit.get('unpaired', 0)} |")
    add("")
    dropped_total = sum(sum((report["experiments"][site]["audit"].get("dropped") or
                             {}).values()) for site in sites)
    noun = TABLE_NOUN[pair[0]] if pair else TABLE_NOUN[sites[0]]
    add("A unit is dropped for one dimension only when the transformed value is "
        "undefined. "
        + ("LUMO loses one `gamma2` row (a zero); Bautzen loses a *pair* "
           "whenever either side of the bridge-control difference is missing or zero."
           if pair == DEFAULT_SITES else
           (f"Neither {noun} table loses a unit for any dimension, so every `R` above "
            "is computed on the full cell." if len(sites) > 1 else
            f"The {noun} table loses no unit for any dimension, so every `R` above is "
            "computed on the full cell.") if not dropped_total else
           f"The {noun} table{'s' if len(sites) > 1 else ''} lose"
           f"{'s' if len(sites) == 1 else ''} {dropped_total} unit-dimension(s) in "
           "total.")
        + ("" if not report.get("bootstrap_caveats", {}).get("cells") else
           " " + _caveat_sentence(report)))
    add("")
    excluded = [(site, state, count) for site in sites
                for state, count in (report["experiments"][site]["audit"].get("excluded")
                                     or {}).items()]
    if excluded:
        add("**Excluded before any statistic:** "
            + "; ".join(f"{EXPERIMENT[site]} `{state}` — {count} acquisition"
                        f"{'' if count == 1 else 's'}: {EXCLUDED[site][state]}; counted "
                        "here and in the JSON audit and never folded into a "
                        "neighbouring state" for site, state, count in excluded)
            + ".")
        add("")
    add("## What this does not claim")
    add("")
    add("* No magnitude comparison between the assets: `R` is scaled by each "
        "reference's own variation — and, where the reference is smaller than ten "
        "units, by that asset's pooled variation — so the numbers live on different "
        "units by construction.")
    add("* No joint verdict. Each asset's Holm family, and therefore each asset's "
        "verdict and bound, is its own.")
    if pair == DEFAULT_SITES:
        severity = ("LUMO's ordering comes from the measured natural frequencies, "
                    "Bautzen's from the documented ladder.")
    elif len(sites) == 1:
        one = " and ".join(f"`{name}`" for name, _states, _dose in LEVELS[sites[0]])
        severity = (f"this asset has no ordered state ladder — {one} are one binary "
                    "event, and `ΔR` is its size, not a dose response.")
    else:
        one = ("`pre` and `post`"
               if {name for site in sites for name, _s, _d in LEVELS[site]}
               == {"pre", "post"} else "the two states")
        severity = (f"neither {noun} has an ordered state ladder — {one} are one "
                    "binary event, and `ΔR` is its size, not a dose response.")
    add("* No severity derived from an OCV response: " + severity)
    add("* No failure case enters the construction of this metric; projecting them "
        "afterwards is a separate step with these definitions frozen.")
    add("")
    if "carola" in sites:
        add("* **One girder, not three.** The support-sized chip keeps only "
            "`Girder A` (436 chips of the committed 1 717); `Girder B` and "
            "`Girder C` are not in this table, so nothing here speaks for them.")
        add("* **Two states, no control.** Carola's `pre` and `post` are separated "
            "by the collapse date and the table carries no control track, so the "
            "post signature is confounded with the calendar — unlike Bautzen, "
            "whose bridge-minus-control difference removes what the two tracks "
            "share.")
        add("* **No ordering claim.** With two states there is no dose ladder to "
            "order: the ordered test reduces to a one-sided rank test between "
            "`pre` and `post`, and only its two-group version is reported.")
    if "morandi" in sites:
        add("* **One platform track.** Morandi's support chip is the counterweight "
            "of a single platform track, so nothing here speaks for the tower as a "
            "whole, and no second track exists to control for what the two dates "
            "share.")
        sizes = {channel: report["experiments"]["morandi"]["channels"][channel]
                 ["blocks"]["geometry"]["n_by_level"]
                 for channel in channels_of(report, "morandi")}
        thin = [f"{channel} (n = {sizes[channel][1]})"
                for channel in channels_of(report, "morandi")
                if sizes[channel][1] <= BLOCK_BOOT]
        add("* **A thin post-collapse cell.** It holds "
            + " and ".join(f"{sizes[channel][1]} {channel}"
                           for channel in channels_of(report, "morandi")) + " units"
            + (f"; that cell is at or below the bootstrap block of {BLOCK_BOOT} on "
               f"{', '.join(thin)}, so its interval cannot resolve it — a null there "
               "is a statement about what this many chips can see, not about the "
               "structure." if thin else "."))
        add("* **No ordering claim.** `pre` and `post` are one binary event, so only "
            "the two-group version of the ordered test is reported.")
    if "haltern" in sites:
        sizes = {channel: report["experiments"]["haltern"]["channels"][channel]
                 ["blocks"]["geometry"]["n_by_level"]
                 for channel in channels_of(report, "haltern")}
        thin = [f"{channel} (n = {sizes[channel][1]})"
                for channel in channels_of(report, "haltern")
                if sizes[channel][1] <= BLOCK_BOOT]
        below = [channel for channel in channels_of(report, "haltern")
                 if sizes[channel][1] < SMALL_N]
        add("* **One section of one mast, no control.** The support-sized chip is the "
            "`turbine-mast-section` of a single AV9 turbine and the table carries no "
            "control track, so the post-collapse signature is confounded with the "
            "calendar — unlike Bautzen, whose bridge-minus-control difference removes "
            "what the two tracks share.")
        add("* **A thin post-collapse cell.** It holds "
            + " and ".join(f"{sizes[channel][1]} {channel}"
                           for channel in channels_of(report, "haltern"))
            + " units against "
            + " and ".join(f"{sizes[channel][0]} {channel}"
                           for channel in channels_of(report, "haltern"))
            + " pre-collapse ones"
            + (f", and {', '.join(thin)} is at or below the bootstrap block of "
               f"{BLOCK_BOOT}, so its interval carries the reference's uncertainty "
               "only" if thin else "")
            + (f"; {', '.join(below)} sits below the six-unit threshold the `Support` "
               "column calls `thin`" if below else "")
            + ". With this many post-collapse overpasses a null is a statement about "
            "what the data can resolve, not about the structure.")
        add("* **`gamma2_mask`, not `gamma2_export`.** The behavior block uses the "
            "masked coherence the window rule defines; the asset-point export the same "
            "table carries is a single point of the scene, not the site's masked "
            "coherence, and is not used.")
        add("* **No ordering claim.** `pre` and `post` are one binary event, so only "
            "the two-group version of the ordered test is reported.")
    if "kdlo" in sites:
        channels = channels_of(report, "kdlo")
        excluded = report["experiments"]["kdlo"]["audit"].get("excluded") or {}
        add("* **A rebuild, not a ruin.** KDLO's second state is the mast as rebuilt, "
            "not the mast as collapsed, so a response there is the difference between "
            "two structures — the standing one and its replacement — and not a damage "
            "state. Nothing here reads severity off the size of `ΔR`.")
        add("* **One channel, half the tests.** "
            + ("The record carries no ASCENDING acquisition, so every cell is a DESC "
               "cell: two H1 tests for the asset instead of four, and two cosine cells "
               "instead of four. The missing channel is declared in the JSON "
               "(`channels_absent`) and drawn as an empty panel; it is never a null."
               if len(channels) == 1 else
               "Both orbit directions are present, so the asset carries the full "
               "four-cell family."))
        add("* **One acquisition excluded, and counted.** "
            + ("; ".join(f"`{state}` — {count} acquisition"
                         f"{'' if count == 1 else 's'} ({EXCLUDED['kdlo'][state]})"
                         for state, count in excluded.items())
               if excluded else "none") + ".")
        add("* **Support chip = whole record.** All 60 acquisitions of the record lie "
            "inside the support, so this table is the record rather than a window on "
            "it — the opposite of Haltern's 107 of 120, and the reason the two counts "
            "are printed rather than compared.")
    if pair in TWO_STATE_PAIRS:
        add(f"* **One metric, two events.** The two {PAIR_NOTE[pair]['noun']} "
            + ("collapse for different reasons on different dates with different "
               "structures" if pair == BRIDGE_PAIR else
               "fail for different reasons on different dates with different "
               "structures, and only one of them was rebuilt")
            + "; the shared element is the measurement and the window rule, not the "
            "event. A shared cosine direction is not a transfer claim, and no "
            "magnitude crosses between the two.")
        add("* **Same window rule, verified once.** " + PAIR_NOTE[pair]["window"])
    return "\n".join(lines)


def write_outputs(report, json_path=None, md_path=None, csv_path=None, quiet=False):
    """JSON, the report, and the compact table — all from the same object."""
    stem = output_stem(report.get("sites") or DEFAULT_SITES)
    json_path = json_path or os.path.join(REPO, "data", "sds", f"{stem}.json")
    md_path = md_path or os.path.join(REPO, "data", "sds", f"{stem}.md")
    csv_path = csv_path or os.path.join(REPO, "data", "sds", f"{stem}_table.csv")
    os.makedirs(os.path.dirname(os.path.abspath(json_path)), exist_ok=True)
    with open(json_path, "w") as handle:
        json.dump(report, handle, indent=1, sort_keys=True, default=str)
        handle.write("\n")
    with open(md_path, "w") as handle:
        handle.write(markdown(report))
    fields = ["Experiment", "State", "Channel", "n", "R_behavior", "R_geometry",
              "Direction", "Scale", "Support"]
    with open(csv_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in report["table"]:
            writer.writerow({key: (f"{value:.6g}" if isinstance(value, float) else value)
                             for key, value in row.items()})
    if not quiet:
        for path in (json_path, md_path, csv_path):
            print(f"  {os.path.relpath(path, REPO)}")
    return json_path, md_path, csv_path


# ---------------------------------------------------------------------------
# The plate — eight z-planes on one 170 mm canvas (amendment 6)
# ---------------------------------------------------------------------------
#: The channels, in the order an asset's two rows are drawn: ASC above DESC.
CHANNEL_ORDER = ("ASC", "DESC")

#: The two columns, in the order ``BLOCKS`` states them.
PLATE_COLUMNS = ("behaviour", "geometry")

#: The canvas.  ``style.FINAL_WIDTH_MM`` is 170 mm -- the text width every
#: figure of this repository is set in.  The pre-registered 2 x 4 grid (rows =
#: experiment, columns = channel x block) could not be kept at that width: four
#: columns inside 170 mm leave 30 mm per square z-plane, and the longest state
#: label is 24.8 mm at the ladder's 6.2 pt annotation floor, so the width rule
#: would be met and the typography rule broken.  Amendment 6 therefore keeps
#: *rows = experiment* and *columns = (channel x block)* as the reading order but
#: nests them: two columns of blocks, two rows per experiment (ASC, then DESC).
#: Eight panels in two columns is four rows, and a row can afford 44 mm of plane.
PLATE_LEFT_MM = 14.0          # y tick labels + the rotated y label
PLATE_RIGHT_MM = 3.0          # a tick at the axis end is labelled half outside it
PLATE_COL_GAP_MM = 11.0       # the right column's y ticks + its label
PLATE_PLANE_MM = 44.0         # the square z-plane: four 6.2 pt labels need it
PLATE_ROW_GAP_MM = 12.0       # x tick labels + x label + the next row's title
PLATE_HEADER_MM = 15.0        # the one-line header, and the ✱ line under it
#: Where the two header lines start, measured from the top of the canvas.  The
#: ✱ line is set on the label rung, so the band has to hold the header (3.6 mm),
#: the line itself (2.9 mm) and 1.5 mm of air above the first row of titles --
#: which is why the header sits at 2 mm and not at 3.
PLATE_TEXT_MM = 2.0
PLATE_LEGEND_MM = 6.2
PLATE_FOOT_MM = 10.0          # x tick labels + x label under the last row
PLATE_DPI = FIGURE_DPI        # 630, as the spec fixes it
#: The state label inside a panel: the ladder's annotation floor, not a size of
#: its own.  Nothing on a plate is written smaller than this.
STATE_PT = style.FS["annotation"]
#: The inset strip (R per level with its interval).  Its height is a tag line, one
#: row per level and the interval's ticks -- a fixed 15 mm box for three levels
#: was a fifth of the plane spent on white space, and it is the plane the labels
#: need.  ``INSET_LABEL_MM`` is the width the level names take on the left,
#: ``INSET_TICK_MM`` the height the interval ticks take at the bottom and
#: ``INSET_TAG_MM`` the height the ``R`` tag takes above the rows; all three are
#: subtracted before anything is drawn, so the footprint that is reserved for the
#: strip is the box the ink actually needs.
INSET_WIDTH_MM = 28.0
INSET_ROW_MM = 3.4             # one row per level: a 6.8 pt label plus its gap
INSET_LABEL_MM = 10.0          # `rebuild` is 8.6 mm at the tick rung
INSET_TICK_MM = 3.6
INSET_TAG_MM = 3.2
#: The drawn content must stay this far inside the canvas: the raster is written
#: from the canvas, so a glyph at the edge would be clipped -- a tight bbox would
#: have hidden that by cropping, which is the defect it must not be.
MARGIN_MIN_MM = 0.3
#: Clearances the guards hold, in millimetres.
LABEL_CLEAR_MM = 0.0          # boxes may touch, never overlap
MARKER_CLEAR_MM = 0.8         # a label keeps this far off a marker
PANEL_CLEAR_MM = 0.4          # a label stays this far inside its own slot
INSET_CLEAR_MM = 1.0          # the R strip clears labels and markers by this
HEADER_CLEAR_MM = 1.5         # the header clears the first row of titles


def plate_layout(rows):
    """The canvas and the grid fractions of a plate with *rows* z-plane rows.

    Pure arithmetic on the constants above, so ``--selftest`` can check the
    canvas without matplotlib: two slots, the gutters and the gap add up to
    exactly ``FINAL_WIDTH_MM``, and the height is the header, the planes, the
    gaps and the foot, with nothing left over for a tight bbox to steal.
    """
    width_mm = style.FINAL_WIDTH_MM
    columns = len(PLATE_COLUMNS)
    slot_mm = (width_mm - PLATE_LEFT_MM - PLATE_RIGHT_MM
               - (columns - 1) * PLATE_COL_GAP_MM) / columns
    height_mm = (PLATE_HEADER_MM + rows * PLATE_PLANE_MM
                 + (rows - 1) * PLATE_ROW_GAP_MM + PLATE_FOOT_MM)
    return {
        "width_mm": width_mm, "height_mm": height_mm,
        "rows": rows, "columns": columns, "slot_mm": slot_mm,
        "plane_mm": PLATE_PLANE_MM,
        #: A z-plane is square (``aspect="equal"``); the slot is wider than tall,
        #: so the x limits are stretched by exactly this factor and the plane
        #: keeps its units.  The lateral margin is where the state labels sit.
        "x_stretch": slot_mm / PLATE_PLANE_MM,
        "figsize": (width_mm / 25.4, height_mm / 25.4),
        "grid": {
            "left": PLATE_LEFT_MM / width_mm,
            "right": 1.0 - PLATE_RIGHT_MM / width_mm,
            "bottom": PLATE_FOOT_MM / height_mm,
            "top": 1.0 - PLATE_HEADER_MM / height_mm,
            "wspace": PLATE_COL_GAP_MM / slot_mm,
            "hspace": PLATE_ROW_GAP_MM / PLATE_PLANE_MM,
        },
    }


def _magnitude_strip(panel_ax, layout, entry, site):
    """R per level with its bootstrap interval, as a strip inside the panel.

    The strip is its own axes (``inset_axes``), so its tick labels are on the
    ladder like every other tick label.  Its height follows the levels -- one row
    each -- because the plane is what the state labels have to live in, and a box
    sized for the worst case would cost every two-level asset a third of its
    pane.  It carries a 6.2 pt ``R`` tag of its own: an axis label would cost the
    strip two thirds of its height.
    """
    levels = [name for name, _states, _dose in LEVELS[site]]
    ci = entry["bootstrap"]["ci"]
    width = INSET_WIDTH_MM / layout["slot_mm"]
    height = ((INSET_TAG_MM + INSET_TICK_MM + len(levels) * INSET_ROW_MM)
              / layout["plane_mm"])
    return {"levels": levels, "ci": {name: ci.get(name) for name in levels},
            "size": [width, height]}


def _draw_magnitude_strip(inset, strip):
    """Draw one prepared strip; the numbers are the ones the JSON holds."""
    for index, name in enumerate(strip["levels"]):
        pair = strip["ci"].get(name)
        if not pair:
            continue
        inset.plot([pair[0], pair[1]], [index, index], color="0.35", lw=0.8)
        inset.plot([(pair[0] + pair[1]) / 2.0], [index], marker="|", ms=4.0,
                   color="#c53030")
    inset.set_yticks(range(len(strip["levels"])))
    inset.set_yticklabels(strip["levels"])
    inset.set_ylim(-0.55, len(strip["levels"]) - 0.45)   # one row, one level
    # R is 0 at the reference level by construction, so its marker must not sit
    # on the axis itself where it would read as a tick
    inset.margins(x=0.08)
    inset.tick_params(labelsize=style.FS["tick"], length=1.5, pad=1.0)
    # the tag is the strip's title, not a mark inside it: at (0.02, 0.98) it sat on
    # the top row's own interval, and next to the level names it read as one of
    # them -- the top right says what the ticks measure without joining a row
    inset.set_title("R", loc="right", pad=1.0, fontsize=STATE_PT, color="0.35")
    for spine in ("top", "right"):
        inset.spines[spine].set_visible(False)
    for spine in ("top", "right"):
        inset.spines[spine].set_visible(False)


#: How a state is written on a plate. The two-state assets use the level name,
#: because `pre-collapse (healthy)` next to `post-collapse` doubles the label
#: length without adding anything the section text does not already say.
STATE_LABEL = {"pre-collapse (healthy)": "pre", "post-collapse": "post"}


def _empty_panel(ax, label, note):
    """A panel for a channel the record does not contain: empty, and says so.

    It is drawn, not skipped: the plate keeps its two rows per asset so the
    assets stay comparable, and the panel states that the cell does not exist.
    An empty cell with an axis and no note would read as a null result.
    """
    ax.set_xlim(-1.0, 1.0)
    ax.set_ylim(-1.0, 1.0)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.set_title(label, loc="left", color="0.55")
    ax.text(0.5, 0.56, "—", transform=ax.transAxes, ha="center", va="center",
            fontsize=13.0, color="0.7")
    ax.text(0.5, 0.42, note, transform=ax.transAxes, ha="center", va="center",
            fontsize=STATE_PT, color="0.45", wrap=True)
    ax.text(0.5, 0.30, "no cell, no test, no interval", transform=ax.transAxes,
            ha="center", va="center", fontsize=STATE_PT, color="0.45")


def _panel(ax, report, site, channel, block_name, layout, fig):
    """One cell of a plate: arrows from the reference, its labels, its R strip.

    Returns what was drawn -- the points, the label texts, the strip -- so the
    audit can pin it and :func:`assert_figure_matches_report` can compare it
    with the numbers the JSON holds.  Nothing here rounds a value: the label
    shows ``_num(value, 2)`` and the point *is* the JSON's ``z``.
    """
    store = report["experiments"][site]
    dims = GEOMETRY_PLANE if block_name == "geometry" else BLOCKS[block_name]
    title = (f"{EXPERIMENT[site]} · {channel} · {BLOCK_DISPLAY[block_name]} "
             f"[{', '.join(DIM_DISPLAY.get(dim, dim) for dim in dims)}]")
    panel = {"site": site, "channel": channel, "block": block_name,
             "title": title, "dims": list(dims), "empty": False,
             "points": {}, "labels": {}, "note": None, "strip": None,
             "inset": None, "limits": None}
    if channel not in store["channels"]:
        note = (store.get("channels_absent") or {}).get(channel, "no acquisition")
        _empty_panel(ax, title, note)
        panel.update(empty=True, absent=note, limits=[[-1.0, 1.0], [-1.0, 1.0]])
        return panel
    cell = store["channels"][channel]
    span = 0.3
    marks = []
    for state in STATES[site]:
        mark = cell["states"][state]["blocks"][block_name]
        x, y = mark["z"].get(dims[0]), mark["z"].get(dims[1])
        if x is None or y is None:
            continue
        marks.append((state, mark, x, y))
        span = max(span, abs(x), abs(y))
    # The plane keeps its units (aspect="equal"); the slot is wider than tall, so
    # the x limits are stretched by exactly slot/plane and the extra width is the
    # margin the state labels are placed in -- outside the data, not over it.
    limit = span * 1.35
    x_lim = limit * layout["x_stretch"]
    ax.set_xlim(-x_lim, x_lim)
    ax.set_ylim(-limit, limit)
    ax.set_aspect("equal")
    ax.axhline(0.0, color="0.75", lw=0.5, ls=":")
    ax.axvline(0.0, color="0.75", lw=0.5, ls=":")
    ax.tick_params(labelsize=style.FS["tick"], length=2.0)
    ax.set_xlabel(f"z({DIM_DISPLAY.get(dims[0], dims[0])})",
                  fontsize=style.FS["label"], labelpad=1.5)
    ax.set_ylabel(f"z({DIM_DISPLAY.get(dims[1], dims[1])})",
                  fontsize=style.FS["label"], labelpad=1.5)
    ax.set_title(title, loc="left", pad=3.0, fontsize=style.FS["title"])
    renderer = fig.canvas.get_renderer()
    to_mm = 25.4 / fig.dpi
    reserved = []
    for state, mark, x, y in marks:
        colour = "#2b6cb0" if state == REFERENCE[site] else "#c53030"
        ax.annotate("", xy=(x, y), xytext=(0.0, 0.0),
                    arrowprops={"arrowstyle": "->", "lw": 1.0, "color": colour,
                                "alpha": 0.9, "shrinkA": 0, "shrinkB": 0})
        ax.plot([x], [y], marker="o", ms=3.0, color=colour, zorder=3)
        # a marker is a box no label may sit on
        reserved.append(_grow(ax.lines[-1].get_window_extent(renderer),
                              MARKER_CLEAR_MM / to_mm))
        count = cell["states"][state].get("n_units")
        thin = "✱" if count is not None and count <= BLOCK_BOOT else ""
        text = (f"{STATE_LABEL.get(state, state)} "
                f"({_num(mark['value'], 2)}; n={count}{thin})")
        panel["points"][state] = {
            "x": x, "y": y,
            "z": {dim: mark["z"].get(dim)
                  for dim in ("A", "D", "F", "S", "P", "gamma2")},
            "value": mark["value"], "n": count, "thin": bool(thin),
            "reference": state == REFERENCE[site], "label": text}
    # D and S stay in the block magnitude but are not on the plane: D partly
    # inherits A algebraically, and S describes position, not morphology.  The
    # note is placed before the labels, so they are laid out around it.
    if block_name == "geometry":
        mark = cell["states"][LARGEST[site]]["blocks"]["geometry"]["z"]
        note = " ".join(f"{dim}{mark[dim]:+.2f}" for dim in ("D", "S")
                        if mark.get(dim) is not None)
        if note:
            artist = ax.text(0.995, 0.995, note, transform=ax.transAxes,
                             ha="right", va="top", fontsize=STATE_PT,
                             color="0.35",
                             bbox={"facecolor": "white", "edgecolor": "none",
                                   "alpha": 0.75, "pad": 0.6})
            reserved.append(_grow(artist.get_window_extent(renderer),
                                  LABEL_CLEAR_MM / to_mm))
            panel["note"] = note
    # The strip goes in before the labels, and is re-placed if they do not fit
    # around it: it is the coarser box -- 28 mm by one row per level -- while a
    # label can step a millimetre at a time.  The windows are tried in order of
    # their clearance and the first one that still leaves the cell's own labels
    # room is kept, so neither element is drawn over the other.
    strip = _magnitude_strip(ax, layout, cell["blocks"][block_name], site)
    inset = footprint = labels = None
    if any(strip["ci"].values()):
        for rect, footprint in _inset_windows(ax, fig, strip, reserved):
            trial = reserved + [_grow(footprint, INSET_CLEAR_MM / to_mm)]
            labels = _lay_out_labels(ax, fig, marks, panel, trial, limit, x_lim)
            if labels is not None:
                inset = rect
                break
    else:
        labels = _lay_out_labels(ax, fig, marks, panel, reserved, limit, x_lim)
    if labels is None:
        raise RuntimeError(
            "the %s panel holds %d state labels and the R strip nowhere: the "
            "plane is too full for both -- widen it rather than shrink the type"
            % (title, len(marks)))
    if inset is not None:
        _backdrop(ax, footprint, fig)
        _draw_magnitude_strip(ax.inset_axes(inset), strip)
        panel["strip"] = dict(strip["ci"])
        panel["inset"] = [round(float(value), 4) for value in inset]
    panel["limits"] = [[-x_lim, x_lim], [-limit, limit]]
    return panel


def _backdrop(ax, footprint, fig):
    """A hairline box behind the whole strip: plot, level names and ticks.

    The inset's own frame is the plotting box only, so without this the level
    names and the interval ticks float in the plane as loose text.  The backdrop
    is drawn on the parent axes, under the inset and under the parent's spines,
    clipped to the slot, and grown by exactly the ``INSET_CLEAR_MM`` the
    reservation uses -- the box no state label may enter is the box the reader
    sees.
    """
    from matplotlib.patches import Rectangle

    frame = ax.get_window_extent(fig.canvas.get_renderer())
    box = _grow(footprint, INSET_CLEAR_MM / (25.4 / fig.dpi))
    ax.add_patch(Rectangle(
        ((box.x0 - frame.x0) / frame.width, (box.y0 - frame.y0) / frame.height),
        box.width / frame.width, box.height / frame.height,
        transform=ax.transAxes, facecolor="white", edgecolor="0.8", lw=0.4,
        zorder=2.1, clip_on=True))


def _lay_out_labels(ax, fig, marks, panel, reserved, limit, x_lim):
    """Place every state label, or none: the artists, or ``None``.

    Being able to *try* a layout is what lets the strip take a window the labels
    can live with: a failed attempt removes its own artists, so the panel is
    handed back untouched and the next window is tried.  A half-labelled panel
    is not a layout to fall back on.
    """
    site = panel["site"]
    placed, drawn = [], []
    for state, _mark, x, y in marks:
        colour = "#2b6cb0" if state == REFERENCE[site] else "#c53030"
        text = panel["points"][state]["label"]
        # the label leaves its marker on the side the other states are not on,
        # and steps off that line until it clears every box already placed
        others = [other_x for other_state, _m, other_x, _y in marks
                  if other_state != state]
        ha = "right" if others and x >= sum(others) / len(others) else "left"
        found = _place_state_label(ax, fig, x, y, text, colour, ha,
                                   placed, reserved, limit, x_lim)
        if found is None:
            for artist in drawn:
                artist.remove()
            return None
        artist, box = found
        placed.append(box)
        drawn.append(artist)
        panel["labels"][state] = text
    return drawn


def _grow(box, pad):
    """*box* grown by *pad* display units on every side."""
    return type(box).from_extents(box.x0 - pad, box.y0 - pad,
                                  box.x1 + pad, box.y1 + pad)


def _place_state_label(ax, fig, x, y, text, colour, ha_first, placed, reserved,
                       limit, x_lim):
    """Put *text* beside ``(x, y)`` where its measured box clears the panel.

    The box is the renderer's, not an estimate from the character count: a
    6.2 pt ``rebuild (0.56; n=43)`` is 24.8 mm wide and a 44 mm plane's lateral
    margin is 13.5 mm, so whether a label fits on one side of its marker is a
    fact about the text, and the search runs on that fact.  Candidate positions
    step away from the marker in millimetres -- up and down, on either side --
    and the first box that is inside the slot and clear of every other box wins.
    ``None`` means this text fits nowhere; the caller decides whether that stops
    the build or only this window of the panel.
    """
    renderer = fig.canvas.get_renderer()
    frame = ax.get_window_extent(renderer)
    to_mm = 25.4 / fig.dpi
    mm_x = 2.0 * x_lim / (frame.width * to_mm)      # data units per millimetre
    mm_y = 2.0 * limit / (frame.height * to_mm)
    clear_mm = PANEL_CLEAR_MM / to_mm
    steps = (0.0, 1.7, -1.7, 3.4, -3.4, 5.1, -5.1, 6.8, -6.8, 8.5, -8.5,
             10.2, -10.2)
    for ha in (ha_first, "left" if ha_first == "right" else "right"):
        dx = -1.0 if ha == "right" else 1.0
        for step in steps:
            artist = ax.text(x + dx * mm_x, y + step * mm_y, text,
                             fontsize=STATE_PT, color=colour, va="bottom",
                             ha=ha,
                             bbox={"facecolor": "white", "edgecolor": "none",
                                   "alpha": 0.75, "pad": 0.6})
            box = _grow(artist.get_window_extent(renderer),
                        LABEL_CLEAR_MM / to_mm)
            inside = (box.x0 >= frame.x0 + clear_mm
                      and box.x1 <= frame.x1 - clear_mm
                      and box.y0 >= frame.y0 + clear_mm
                      and box.y1 <= frame.y1 - clear_mm)
            if inside and not any(box.overlaps(other)
                                  for other in placed + reserved):
                return artist, box
            artist.remove()
    return None


def _inset_windows(ax, fig, strip, reserved):
    """The windows the R strip could take, freest first.

    The strip is placed *before* the state labels, because it is the coarser box:
    26 mm by one row per level has to find a window the data leaves empty, while
    a label can step a millimetre at a time.  The four corners and the four edge
    midpoints are scored by the distance from the strip's footprint to the
    nearest marker or note; the caller takes them in that order and keeps the
    first window whose footprint still leaves the labels room for themselves.

    Each window is ``(rect, footprint)``: ``rect`` in axes fractions for
    ``inset_axes``, ``footprint`` in display units for the caller's reservation.
    A window with less than ``INSET_CLEAR_MM`` of clear space is not offered: the
    strip is an annotation and the marker is the datum, so the strip moves, not
    the datum.
    """
    from matplotlib.transforms import Bbox

    renderer = fig.canvas.get_renderer()
    frame = ax.get_window_extent(renderer)
    to_mm = 25.4 / fig.dpi
    width, height = strip["size"]                  # fractions of the slot
    margin = 0.02
    lo_x, hi_x = margin, 1.0 - margin - width
    lo_y, hi_y = margin, 1.0 - margin - height
    mid_x = 0.5 * (lo_x + hi_x)
    mid_y = 0.5 * (lo_y + hi_y)
    windows = {"bl": (lo_x, lo_y), "br": (hi_x, lo_y),
               "tl": (lo_x, hi_y), "tr": (hi_x, hi_y),
               "bottom": (mid_x, lo_y), "top": (mid_x, hi_y),
               "left": (lo_x, mid_y), "right": (hi_x, mid_y)}

    def footprint(where):
        x0 = frame.x0 + where[0] * frame.width
        y0 = frame.y0 + where[1] * frame.height
        return Bbox.from_extents(x0, y0, x0 + width * frame.width,
                                 y0 + height * frame.height)

    def clearance(where):
        box = footprint(where)
        gaps = []
        for other in reserved:
            dx = max(other.x0 - box.x1, box.x0 - other.x1, 0.0)
            dy = max(other.y0 - box.y1, box.y0 - other.y1, 0.0)
            gaps.append(math.hypot(dx, dy))
        return min(gaps) if gaps else float("inf")

    order = sorted(windows, key=lambda name: -clearance(windows[name]))
    # the plotting box inside the footprint: the level names take their width on
    # the left, the interval ticks their height at the bottom and the R tag its
    # height above the top row
    pad_x = INSET_LABEL_MM / (frame.width * to_mm)
    pad_y = INSET_TICK_MM / (frame.height * to_mm)
    pad_top = INSET_TAG_MM / (frame.height * to_mm)
    found = []
    for name in order:
        where = windows[name]
        if clearance(where) * to_mm < INSET_CLEAR_MM:
            break
        found.append(([where[0] + pad_x, where[1] + pad_y,
                       width - pad_x, height - pad_y - pad_top],
                      footprint(where)))
    return found


def _thin_legend(report, sites):
    """The ✱ legend as one line, naming the cells that carry the mark.

    It is a second header line, on the label rung of the ladder, not part of the
    title: the title says what the plate shows, the legend explains a mark.
    """
    thin = []
    for site in sites:
        for channel in channels_of(report, site):
            cell = report["experiments"][site]["channels"][channel]
            for state in STATES[site]:
                count = cell["states"][state].get("n_units")
                if count is not None and count <= BLOCK_BOOT:
                    label = STATE_LABEL.get(state, state)
                    thin.append(f"{EXPERIMENT[site]} {channel} {label} n={count}")
    if not thin:
        return ""
    return ("✱ n ≤ %d (%s): resampled as one block, so its interval carries the "
            "reference's uncertainty only" % (BLOCK_BOOT, ", ".join(thin)))


def _thin_rule():
    """The same rule without the cell names, for a plate whose list is too long.

    One line at 7.6 pt holds about 150 characters, and a named list grows with
    every thin cell: on the pair whose full length it is 216 mm the names are what
    is dropped, not the sentence -- the guide still lists the cells, and the mark
    is on the cell in the panel.
    """
    return ("✱ n ≤ %d: the cell is one resampling block, so its interval carries "
            "the reference's uncertainty only" % BLOCK_BOOT)


def _legend_line(fig, named):
    """The ✱ line that fits the canvas: the cells by name, else only the rule."""
    if not named:
        return ""
    renderer = fig.canvas.get_renderer()
    to_mm = 25.4 / fig.dpi
    room = style.FINAL_WIDTH_MM - 1.5        # the line starts 1 mm in
    for text in (named, _thin_rule()):
        probe = fig.text(0.006, 0.5, text, fontsize=style.FS["label"])
        width = probe.get_window_extent(renderer).width * to_mm
        probe.remove()
        if width <= room:
            return text
    raise RuntimeError(
        "the ✱ line cannot be set on one %.0f mm row even without the cell names"
        % style.FINAL_WIDTH_MM)


def _plate(report, sites, figdir, stem, headline):
    """Draw one plate at 170 mm, guard it, write jpg + pdf + eps; return the audit.

    The canvas is authored, not measured from the content: 170 mm is the print
    width and the panels are the ones ``plate_layout`` sizes.  A tight bbox would
    re-cut this after the fact and shrink every glyph with it, which is the
    defect the guards below exist to make impossible to hit silently.
    """
    import matplotlib
    matplotlib.use("Agg")
    plt = style.apply(dpi=PLATE_DPI)

    layout = plate_layout(len(sites) * len(CHANNEL_ORDER))
    fig, slots = plt.subplots(layout["rows"], layout["columns"],
                              figsize=layout["figsize"], squeeze=False)
    fig.subplots_adjust(**layout["grid"])
    drawn = {"stem": stem, "sites": list(sites), "layout": layout,
             "panel_order": [], "panels": {}}
    titles = {}
    for row, site in enumerate(sites):
        for offset, channel in enumerate(CHANNEL_ORDER):
            line = row * len(CHANNEL_ORDER) + offset
            for column, block_name in enumerate(PLATE_COLUMNS):
                ax = slots[line][column]
                key = f"{site}/{channel}/{block_name}"
                drawn["panels"][key] = _panel(ax, report, site, channel,
                                              block_name, layout, fig)
                titles[key] = ax._left_title          # `ax.title` is the centre one
                drawn["panel_order"].append(key)
    header = fig.suptitle(headline, x=0.006, ha="left", va="top",
                          y=1.0 - PLATE_TEXT_MM / layout["height_mm"],
                          fontsize=style.FS["suptitle"])
    thin = _legend_line(fig, _thin_legend(report, sites))
    legend = None
    if thin:
        legend = fig.text(0.006, 1.0 - PLATE_LEGEND_MM / layout["height_mm"],
                          thin, ha="left", va="top",
                          fontsize=style.FS["label"], color="0.35")
    # --- guards -------------------------------------------------------------
    # These are the checks a tight bbox hides: with the canvas as the export box
    # a glyph is either inside it or cut off, and a panel that grew by half a
    # point collides instead of pushing the frame wider.
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    to_mm = 25.4 / fig.dpi
    header_box = header.get_window_extent(renderer)
    rows_here = [("header", header_box)]
    if legend is not None:
        rows_here.append(("✱ legend", legend.get_window_extent(renderer)))
    for name, span in rows_here:
        if (span.x0 < 0 or span.x1 > fig.bbox.width
                or span.y0 < 0 or span.y1 > fig.bbox.height):
            raise RuntimeError("the %s leaves the canvas" % name)
        if (fig.bbox.width - span.x1) * to_mm < MARGIN_MIN_MM:
            raise RuntimeError(
                "the %s stops %.1f mm from the right edge of the canvas: "
                "shorten it" % (name, (fig.bbox.width - span.x1) * to_mm))
    if legend is not None:
        gap = (header_box.y0 - legend.get_window_extent(renderer).y1) * to_mm
        if gap < 0.4:
            raise RuntimeError("the header and the ✱ line touch (%.1f mm)" % gap)
    first_rows = len(CHANNEL_ORDER) * len(PLATE_COLUMNS)
    top = max(titles[key].get_window_extent(renderer).y1
              for key in drawn["panel_order"][:first_rows])
    low = rows_here[-1][1].y0
    if (low - top) * to_mm < HEADER_CLEAR_MM:
        raise RuntimeError("the header clears the first row of titles by only "
                           "%.1f mm" % ((low - top) * to_mm))
    tight = fig.get_tightbbox(renderer)                       # inches
    edges = {"left": tight.x0, "right": layout["figsize"][0] - tight.x1,
             "bottom": tight.y0, "top": layout["figsize"][1] - tight.y1}
    worst = min(edges.values()) * 25.4
    if worst < MARGIN_MIN_MM:
        raise RuntimeError(
            "the drawn content is %.2f mm from the canvas edge (min %.1f mm): "
            "re-cut the grid -- do not crop the export" % (worst, MARGIN_MIN_MM))
    assert_figure_matches_report(report, drawn)
    translucent = style.alpha_artists(fig)     # counted before they are blended
    out = os.path.join(figdir, f"{stem}.jpg")
    style.save(fig, out, dpi=PLATE_DPI, eps="blend")
    plt.close(fig)
    return plate_audit(report, drawn, out, tight, edges, translucent, headline,
                       thin)


def plate_audit(report, drawn, fig_path, tight, edges, translucent, headline,
                thin):
    """The machine-readable record of a plate, written into the guide.

    It answers what a reader of the file cannot see by looking: how wide each of
    the three containers really is, what happened to the alpha in the EPS, which
    numbers each panel pins, and where the labels and the strip ended up.
    Nothing is recomputed -- these are the values that were drawn.
    """
    stem = os.path.splitext(fig_path)[0]
    layout = drawn["layout"]
    audit = {
        "schema": "sds-plate-audit-v1",
        "figure": os.path.relpath(fig_path, REPO),
        "headline": headline,
        "thin_legend": thin or None,
        "target_width_mm": style.FINAL_WIDTH_MM,
        "target_height_mm": round(layout["height_mm"], 2),
        "figsize_inches": [round(layout["figsize"][0], 6),
                           round(layout["figsize"][1], 6)],
        "native_dpi": PLATE_DPI,
        #: Keyed by the ``figure_style`` role, so a check can read its band.
        "final_font_sizes_pt": {
            "suptitle": style.FS["suptitle"],
            "title": style.FS["title"],
            "label": style.FS["label"],
            "tick": style.FS["tick"],
            "annotation": STATE_PT,
        },
        "native_width_mm": {}, "native_height_mm": {},
        "latex_display_width_mm": style.FINAL_WIDTH_MM,
        "latex_scaling": "the file is 170 mm wide -- included with "
                         "`width=\\linewidth`; no trim, clip or resizebox",
        "layout": {
            "rows": layout["rows"], "columns": layout["columns"],
            "plane_mm": layout["plane_mm"],
            "slot_mm": round(layout["slot_mm"], 2),
            "x_stretch": round(layout["x_stretch"], 4),
            "grid": {key: round(value, 5)
                     for key, value in layout["grid"].items()},
            "reading_order": "rows = experiment x channel (ASC above DESC), "
                             "columns = block ("
                             + ", ".join(BLOCK_DISPLAY[name]
                                         for name in PLATE_COLUMNS) + ")",
        },
        "bounds": {
            "export": "canvas-authorised (figure_style.save); bbox_inches is "
                      "never 'tight'",
            "canvas_mm": [round(layout["width_mm"], 2),
                          round(layout["height_mm"], 2)],
            "content_mm": [round(tight.width * 25.4, 1),
                           round(tight.height * 25.4, 1)],
            "edge_mm": {key: round(value * 25.4, 2)
                        for key, value in edges.items()},
            "minimum_edge_mm": MARGIN_MIN_MM,
        },
        "eps": {"alpha_artists": len(translucent), "strategy": "blend",
                "background": "white", "rasterised_artists": 0,
                "vector_text": "Type 42 (pdf/ps fonttype 42)"},
        "panel_order": list(drawn["panel_order"]),
        "panels": {},
    }
    for path in (fig_path, stem + ".pdf", stem + ".eps"):
        info = style.measure(path)
        audit["native_width_mm"][info["format"]] = round(info["mm"], 2)
        audit["native_height_mm"][info["format"]] = round(info["height_mm"], 2)
    for key in drawn["panel_order"]:
        panel = drawn["panels"][key]
        audit["panels"][key] = {
            "title": panel["title"], "empty": panel["empty"],
            "dims": panel["dims"], "note": panel["note"],
            "limits": [[round(panel["limits"][0][0], 4),
                        round(panel["limits"][0][1], 4)],
                       [round(panel["limits"][1][0], 4),
                        round(panel["limits"][1][1], 4)]],
            "states": {state: {"z": {dim: (None if value is None
                                           else round(value, 6))
                                     for dim, value in point["z"].items()},
                               "x": round(point["x"], 6),
                               "y": round(point["y"], 6),
                               "R": (None if point["value"] is None
                                     else round(point["value"], 6)),
                               "n": point["n"], "thin": point["thin"],
                               "reference": point["reference"],
                               "label": point["label"]}
                       for state, point in panel["points"].items()},
            "labels": panel["labels"],
            "strip": ({name: [round(pair[0], 6), round(pair[1], 6)]
                       for name, pair in panel["strip"].items() if pair}
                      if panel["strip"] else None),
            "strip_inset": panel["inset"],
        }
    return audit


def assert_figure_matches_report(report, drawn):
    """The plate is the report: same points, same labels, same intervals.

    Runs on every build, so a layout change can never quietly move a value: the
    drawn marker has to be the JSON's ``z``, the drawn label the JSON's ``R``
    and ``n``, and the drawn strip the JSON's interval.  A missing panel is a
    defect too -- an empty cell is drawn, and says that it is empty.
    """
    problems = []
    for key in drawn["panel_order"]:
        site, channel, block_name = key.split("/")
        store = report["experiments"][site]
        panel = drawn["panels"][key]
        if channel not in store["channels"]:
            if not panel["empty"]:
                problems.append(f"{key}: drawn with data, the record has none")
            continue
        if panel["empty"]:
            problems.append(f"{key}: drawn empty, the record has the cell")
            continue
        cell = store["channels"][channel]
        for state in STATES[site]:
            mark = cell["states"][state]["blocks"][block_name]
            point = panel["points"].get(state)
            if point is None:
                problems.append(f"{key}/{state}: not drawn")
                continue
            for dim in panel["dims"]:
                if point["z"][dim] != mark["z"].get(dim):
                    problems.append(
                        f"{key}/{state}: z({dim}) {point['z'][dim]!r} != "
                        f"{mark['z'].get(dim)!r}")
            if point["n"] != cell["states"][state].get("n_units"):
                problems.append(f"{key}/{state}: n {point['n']} != "
                                f"{cell['states'][state].get('n_units')}")
            if (point["value"] is None) != (mark["value"] is None):
                problems.append(f"{key}/{state}: R {point['value']!r} != "
                                f"{mark['value']!r}")
            elif (point["value"] is not None
                  and abs(point["value"] - mark["value"]) > 1e-12):
                problems.append(f"{key}/{state}: R {point['value']!r} != "
                                f"{mark['value']!r}")
            if _num(mark["value"], 2) not in point["label"]:
                problems.append(f"{key}/{state}: the label does not quote R "
                                f"({point['label']!r})")
            if f"n={point['n']}" not in point["label"]:
                problems.append(f"{key}/{state}: the label does not quote n "
                                f"({point['label']!r})")
        strip = cell["blocks"][block_name]["bootstrap"]["ci"]
        for level, _states, _dose in LEVELS[site]:
            pair = strip.get(level)
            got = (panel["strip"] or {}).get(level)
            if (pair is None) != (got is None):
                problems.append(f"{key}/{level}: the strip {got} != {pair}")
            elif pair is not None and any(abs(a - b) > 1e-12
                                          for a, b in zip(got, pair)):
                problems.append(f"{key}/{level}: the strip {got} != {pair}")
    if problems:
        raise RuntimeError("the plate does not match the report:\n  "
                           + "\n  ".join(problems))




#: One line per asset's plate header.  On the plate the header states *what is
#: shown*; the long description (chip, window rule, which states) lives in the
#: guide, because a 170 mm line at 9 pt holds about 85 characters and a header
#: that has to be shrunk to fit is a header nobody can read.
PLATE_HEADLINE = {
    "lumo": "OCV reconfiguration — LUMO tower, against its own reference state",
    "bautzen": "OCV reconfiguration — Bautzen bridge, against its own ordinal "
               "state ladder",
    "carola": "OCV reconfiguration — Carolabrücke, pre- against post-collapse",
    "morandi": "OCV reconfiguration — Ponte Morandi, pre- against post-collapse",
    "haltern": "OCV reconfiguration — Haltern turbine, standing against fallen",
    "kdlo": "OCV reconfiguration — KDLO mast, standing against rebuilt",
}

#: The juxtaposition header per registered pair; its stem comes from the site tuple.
#: The default pair's header names the states as *controlled* and stops there: a
#: shared scale is what the plate does **not** have, so the title must not read
#: like a common dose axis the two assets both sit on.
PAIR_HEADLINE = {
    DEFAULT_SITES: "OCV reconfiguration across controlled LUMO and Bautzen states",
    BRIDGE_PAIR: "OCV reconfiguration across a collapse — Carolabrücke and Ponte "
                 "Morandi (no shared scale)",
    TOWER_PAIR: "OCV reconfiguration across a collapse and a rebuild — Haltern and "
                "KDLO (no shared scale)",
}


def make_figure(report, figdir=None):
    """One juxtaposition plate for a registered pair, one plate per asset otherwise.

    A pair writes **one** plate -- the juxtaposition, both assets side by side.
    The per-asset companions of earlier revisions are no longer written: the
    figure set of this repository is a fixed list, and the pair's plate already
    carries each asset's panels.  An asset named alone on ``--sites`` still gets
    its own plate.

    Returns one audit per plate (see :func:`plate_audit`), each carrying the
    path of the raster it was measured from, so the guide can print the widths
    and the caller can name the files without guessing.
    """
    figdir = figdir or os.path.join(REPO, "figures", "sds")
    os.makedirs(figdir, exist_ok=True)
    sites = tuple(report.get("sites") or DEFAULT_SITES)
    pair = pair_of(sites)
    if pair:
        plates = [(figure_stem(pair), tuple(pair), PAIR_HEADLINE[pair])]
    else:
        plates = [(f"fig_sds_ocv_reconfigure_{site}", (site,),
                   PLATE_HEADLINE[site]) for site in sites]
    return [_plate(report, subset, figdir, stem, headline)
            for stem, subset, headline in plates]


def figure_markdown(report, plates=None):
    """The panel guide: a pair's juxtaposition, or one plate per asset.

    *plates* are the audits :func:`make_figure` returned.  When they are passed,
    the guide also carries the plate record -- the measured width of all three
    containers, the rungs that were used and, per panel, the numbers that were
    drawn -- so the guide and the file cannot drift apart silently.
    """
    sites = tuple(report.get("sites") or DEFAULT_SITES)
    stem = output_stem(sites)
    pair = pair_of(sites)
    lines = ["# `fig_sds_ocv_reconfigure*.jpg` — panel guide", "",
             f"*generated by `code/sds/ocv_reconfigure.py`; the numbers behind every "
             f"arrow are in `data/sds/{stem}.json`*", "",
             ("One plate — the juxtaposition of both assets:" if pair
              else "One plate — the asset standing alone:"), "",
             "| file | role |", "| --- | --- |"]
    if pair:
        lines.append(f"| `{figure_stem(pair)}.jpg` | the juxtaposition — rows = "
                     f"experiment × channel (ASC above DESC), columns = block "
                     f"({' and '.join(BLOCK_DISPLAY[name] for name in PLATE_COLUMNS)}); "
                     f"every panel has its own plane |")
    else:
        for site in sites:
            lines.append(f"| `fig_sds_ocv_reconfigure_{site}.jpg` | "
                         f"{EXPERIMENT[site]} alone ({_panels_note(report, site)}) |")
    lines += ["",
             "In every panel the origin is that asset's **reference state** and the "
             "arrows end at its standardised responses "
             "(`z = (median − median_ref) / scale`), so the geometry is "
             "within-asset and no scale is shared between assets.", ""]
    for site in (report.get("sites") or DEFAULT_SITES):
        block = report["reading"]["per_asset"][site]
        if site == "lumo":
            levels = ", ".join(f"`{name}` ({_num(dose, 2)} % frequency drop)"
                               for name, _states, dose in LEVELS[site])
        elif len(STATES[site]) > 2:
            levels = ", ".join(f"`{name}` (rank {_num(dose, 2)})"
                               for name, _states, dose in LEVELS[site])
        else:
            levels = " and ".join(
                f"`{name}` ({'reference' if states == (REFERENCE[site],) else 'the single contrast'})"
                for name, states, _dose in LEVELS[site])
        lines.append(f"* **{EXPERIMENT[site]}** — reference `{REFERENCE[site]}` "
                     + "(n_ref "
                     + " / ".join(f"{block['n_ref'][channel]} {channel}"
                                  for channel in channels_of(report, site))
                     + "), scale "
                     + " / ".join(f"`{block['scale_used'][channel]}`"
                                  for channel in channels_of(report, site))
                     + f", levels {levels}; {BLOCK_DISPLAY['behaviour']} "
                     f"`[{'`, `'.join(BLOCKS['behaviour'])}]`, {BLOCK_DISPLAY['geometry']} "
                     f"`[{'`, `'.join(BLOCKS['geometry'])}]` (plane "
                     f"`[{', '.join(GEOMETRY_PLANE)}]`).")
    lines += ["",
              "The inset in each panel repeats the block magnitude `R` per level with "
              "its moving-block 95 % interval — the quantity H1 tests. The reference "
              "level sits at the origin by construction and carries no interval.", "",
              "Every state label carries its cell's `n`; a cell at or below the "
              "bootstrap block is marked ✱ and its interval then reflects the "
              "reference's uncertainty only.", ""]
    if report.get("bootstrap_caveats", {}).get("cells"):
        lines += [_caveat_sentence(report), ""]
    if pair in TWO_STATE_PAIRS:
        lines += [f"Both rows come from the same window rule in "
                  f"`erlking_mvp/tools/ocv_table.py` (30 % of the peak and 5 × the "
                  f"ring median, 200 m ring, 30 m guard), which is what makes the "
                  f"juxtaposition a comparison of one rule rather than of two.", ""]
    empty = [f"{EXPERIMENT[site]} {channel}"
             for site in sites
             for channel in (report["experiments"][site].get("channels_absent")
                             or {})]
    if empty:
        lines += ["A panel for a channel the record does not contain is drawn empty "
                  "and says so — " + ", ".join(empty) + " — because the plate keeps "
                  "its two rows per asset so that the assets stay comparable. It "
                  "carries no arrow, no inset and no ✱: there is no cell to mark.", ""]
    lines += ["Nothing is clipped: the axis limits follow the largest coordinate in "
              "the panel.", ""]
    lines += _plate_section(plates)
    return "\n".join(lines)


def _plate_section(plates):
    """The plate record: one 170 mm canvas per file, measured after it was written."""
    if not plates:
        return []
    lines = ["## The plate", "",
             "Every plate is authored at 170 mm — the text width of "
             "`code/figure_style.py` — and written from that canvas as `.jpg` "
             f"({PLATE_DPI} dpi), `.pdf` and `.eps`. It is never cropped to its "
             "content: `figure_style.save()` refuses a layout that does not fit, "
             "measures all three containers afterwards and aborts on a deviation "
             "of more than 0.5 mm.", "",
             "| rung | size | used for |", "| --- | --- | --- |",
             f"| header | {style.FS['suptitle']:.1f} pt | the one-line headline, "
             "left-aligned |",
             f"| panel title | {style.FS['title']:.1f} pt | experiment · channel · "
             "block · the plane's pair |",
             f"| axis label | {style.FS['label']:.1f} pt | `z(·)` on both axes |",
             f"| tick label | {style.FS['tick']:.1f} pt | the plane's ticks and the "
             "inset's levels |",
             f"| annotation | {STATE_PT:.1f} pt | the state labels, the `D`/`S` note "
             "and the inset's `R` |", "",
             "The plane is square (`aspect=\"equal\"`): the axes are stretched "
             "horizontally by `slot / plane`, so the units stay equal in both "
             "directions and the added width is the margin the state labels are "
             "placed in — outside the data, never over it.", ""]
    for audit in plates:
        layout = audit["layout"]
        lines += [f"### `{os.path.basename(audit['figure'])}`", "",
                  f"* {layout['rows']} rows × {layout['columns']} columns, plane "
                  f"{layout['plane_mm']:.0f} mm square in a {layout['slot_mm']:.2f} mm "
                  f"slot, canvas {audit['bounds']['canvas_mm'][0]:.1f} × "
                  f"{audit['bounds']['canvas_mm'][1]:.1f} mm.",
                  f"* measured width: "
                  + ", ".join(f"{name} {value:.2f} mm"
                              for name, value in audit["native_width_mm"].items())
                  + f" (target {audit['target_width_mm']:.1f} mm).",
                  f"* the drawn content stops "
                  f"{audit['bounds']['edge_mm']['left']:.2f} / "
                  f"{audit['bounds']['edge_mm']['right']:.2f} / "
                  f"{audit['bounds']['edge_mm']['bottom']:.2f} / "
                  f"{audit['bounds']['edge_mm']['top']:.2f} mm from the left, "
                  f"right, bottom and top edge.",
                  f"* the EPS carries the {audit['eps']['alpha_artists']} "
                  f"translucent artist(s) blended onto white; text stays vector "
                  f"(Type 42).", ""]
    lines += ["The audit below is the machine-readable record of the files that were "
              "written — the widths from the containers themselves, and per panel "
              "the `z` that was drawn, its `R`, its `n` and its interval. "
              "`--verify` re-reads it and re-derives every one of those numbers "
              "from `data/sds/` and the records.", ""]
    for audit in plates:
        lines += ["```json", json.dumps(audit, indent=1, ensure_ascii=False), "```",
                  ""]
    return lines


def _panels_note(report, site):
    """``4 panels``, or how many of them are empty and why."""
    absent = report["experiments"][site].get("channels_absent") or {}
    if not absent:
        return "4 panels"
    return ("4 panels, %d empty — %s" % (2 * len(absent), " and ".join(
        f"no {channel}" for channel in absent)))


def _synthetic(site, effect=0.0, n=40, seed=3):
    """Units with a planted dose effect in `A` and `F`, on the analysis scale."""
    rng = np.random.default_rng(seed)
    states = STATES[site]
    units = []
    for index, state in enumerate(states):
        dose = index / max(1, len(states) - 1)
        for _ in range(n):
            area = math.exp(2.0 + effect * dose + rng.normal(scale=0.25))
            frag = max(1.0, 5.0 + 6.0 * effect * dose + rng.normal(scale=0.6))
            units.append({
                "state": state, "channel": "ASC",
                "date": f"2024-01-{index:02d}",
                "values": {"A": math.log(area), "D": math.log(area / 400.0),
                           "F": math.log(frag), "P": 0.02 + rng.normal(scale=0.002),
                           "S": abs(rng.normal(scale=3.0)),
                           "gamma2": 0.3 + 0.05 * rng.normal()}})
    return units


def selftest():
    """The behaviors the metric must have, on synthetic units."""
    failures = []
    counted = {"checks": 0}

    def check(label, condition, detail=""):
        counted["checks"] += 1
        if not condition:
            failures.append(f"{label}: {detail}")

    # 1./2. the ordered test on the geometry block: it must find a planted
    # progression and must not manufacture one. Three seeds each, because a single
    # draw of a permutation p is itself noisy.
    def ordered_p(effect, seed):
        units = _synthetic("lumo", effect=effect, seed=seed)
        refs = {d: reference_stats(units, "lumo", "ASC", d)
                for d in BLOCKS["geometry"]}
        per_unit = unit_magnitudes(units, "lumo", "ASC", BLOCKS["geometry"], refs)
        groups = [[x["R"] for x in per_unit if x["state"] in states]
                  for _n, states, _d in LEVELS["lumo"]]
        return jonckheere(groups, n_perm=2000, seed=seed)["p_increasing"]

    planted = [ordered_p(2.0, seed) for seed in (3, 5, 7)]
    check("a planted progression is detected in at least 2 of 3 seeds",
          sum(1 for p in planted if p < 0.01) >= 2, str(planted))
    null = [ordered_p(0.0, seed) for seed in (11, 12, 13)]
    check("a null draws p > 0.05 in at least 2 of 3 seeds",
          sum(1 for p in null if p > 0.05) >= 2, str(null))

    # 3. the scale chain names its step
    scale, kind = robust_scale([1.0])
    check("too few values is named", scale is None and kind == "too_few",
          f"{scale} {kind}")
    scale, kind = robust_scale([2.0, 2.0, 2.0, 2.0])
    check("a constant series is named", scale is None and kind == "degenerate",
          f"{scale} {kind}")

    # 4. the cosine behaves — identical, mirrored, orthogonal
    same, _k = cosine({"A": 1.0, "F": 2.0}, {"A": 1.0, "F": 2.0}, ("A", "F"))
    mirrored, _k = cosine({"A": 1.0, "F": 2.0}, {"A": -1.0, "F": -2.0}, ("A", "F"))
    orthogonal, _k = cosine({"A": 1.0, "F": 2.0}, {"A": 2.0, "F": -1.0}, ("A", "F"))
    check("cosine of identical vectors is 1", abs(same - 1.0) < 1e-12, str(same))
    check("cosine of mirrored vectors is -1", abs(mirrored + 1.0) < 1e-12, str(mirrored))
    check("cosine of orthogonal vectors is 0", abs(orthogonal) < 1e-12, str(orthogonal))

    # 5. the sign pattern marks a dropped dimension instead of inventing one
    pattern = sign_pattern({"A": 1.0, "F": None}, ("A", "F"))
    check("a dropped dimension is marked", pattern == "A+ F·", pattern)

    # 6. zeros are dropped for the log scale only
    check("a zero survives on the identity scale", transform("lumo", "S", 0.0) == 0.0)
    check("a zero is dropped on the log scale", transform("lumo", "A", 0.0) is None)

    # 7. the reference level carries no interval, by construction
    units = _synthetic("lumo", effect=1.0, n=12, seed=5)
    boot = bootstrap_magnitudes(units, "lumo", "ASC", BLOCKS["behaviour"], boot=120, seed=7)
    reference_name = LEVELS["lumo"][0][0]
    check("the reference level has a [0, 0] interval",
          boot["ci"].get(reference_name) == [0.0, 0.0] and
          reference_name in boot["by_construction"], str(boot["ci"]))

    # 8. the scale rule: below ten reference units the pooled MAD takes over
    small = _synthetic("bautzen", effect=0.0, n=2)
    large = _synthetic("lumo", effect=0.0, n=40)
    check("a small reference switches to the pooled MAD",
          reference_stats(small, "bautzen", "ASC", "A")["scale_used"] == "pooled_mad",
          str(reference_stats(small, "bautzen", "ASC", "A")["scale_used"]))
    check("a large reference keeps its own MAD",
          reference_stats(large, "lumo", "ASC", "A")["scale_used"] == "mad",
          str(reference_stats(large, "lumo", "ASC", "A")["scale_used"]))

    # 9. Holm is monotone and capped
    adjusted = holm({"a": 0.001, "b": 0.04, "c": 0.03})
    check("Holm adjusts and caps",
          adjusted["a"] == 0.003 and adjusted["c"] == 0.06 and adjusted["b"] == 0.06,
          str(adjusted))

    # 10. the pair machinery: a registered pair is recognised, a lone asset is not
    check("a registered pair is recognised in either order",
          pair_of(("bautzen", "lumo")) == DEFAULT_SITES and
          pair_of(("morandi", "carola")) == BRIDGE_PAIR,
          str(pair_of(("morandi", "carola"))))
    check("a single asset has no pair and keeps its own stem",
          pair_of(("carola",)) is None and output_stem(("carola",)) == "ocv_carola"
          and figure_stem(BRIDGE_PAIR) == "fig_sds_ocv_reconfigure_carola_morandi"
          and output_stem(DEFAULT_SITES) == "ocv_reconfigure",
          output_stem(BRIDGE_PAIR))

    # 11. a cell at or below the bootstrap block is flagged, a larger one is not
    stub = {"sites": ["bautzen"],
            "experiments": {"bautzen": {"channels": {
                channel: {"blocks": {name: {"n_by_level": [BLOCK_BOOT, 9, 9]}
                                     for name in BLOCKS}}
                for channel in ("ASC", "DESC")}}}}
    caveats = bootstrap_caveats(stub)
    check("a cell at or below the block is flagged, a larger one is not",
          len(caveats["cells"]) == 4 and
          all(cell["n"] == BLOCK_BOOT and cell["level"] == "RS"
              for cell in caveats["cells"]), str(caveats["cells"]))

    # 12. the registry is complete and self-consistent across all six assets
    registry = {"COLUMNS": COLUMNS, "TRANSFORM": TRANSFORM, "REFERENCE": REFERENCE,
                "STATES": STATES, "LEVELS": LEVELS, "LARGEST": LARGEST,
                "MIDDLE": MIDDLE, "EXPERIMENT": EXPERIMENT, "TABLE_NOUN": TABLE_NOUN}
    gaps = {name: sorted(set(TABLE) ^ set(mapping))
            for name, mapping in registry.items() if set(TABLE) ^ set(mapping)}
    # Bautzen is the one asset with no single-row table -- it is read on the paired
    # bridge-minus-control difference -- so `UNIT_SHAPE` names every asset but that
    # one, and the gap is asserted rather than tolerated.
    if (set(TABLE) - {"bautzen"}) ^ set(UNIT_SHAPE):
        gaps["UNIT_SHAPE"] = sorted((set(TABLE) - {"bautzen"}) ^ set(UNIT_SHAPE))
    check("every registered asset appears in every registry table", not gaps, str(gaps))
    broken = [site for site in TABLE
              if REFERENCE[site] not in STATES[site]
              or LARGEST[site] not in STATES[site]
              or any(state not in STATES[site] for state in MIDDLE[site])
              or LEVELS[site][0][1] != (REFERENCE[site],)
              or LEVELS[site][-1][1] != (LARGEST[site],)
              or set(DIMS) - set(COLUMNS[site]) or set(DIMS) - set(TRANSFORM[site])]
    check("the reference is the first level and the largest state the last",
          not broken, str(broken))
    check("an excluded state is never also a registered state",
          set(EXCLUDED) <= set(TABLE)
          and all(state not in STATES[site]
                  for site, drop in EXCLUDED.items() for state in drop),
          str(EXCLUDED))

    # 13. the tower pair is registered and keeps its own stem
    check("the tower pair is recognised in either order and named once",
          pair_of(("kdlo", "haltern")) == TOWER_PAIR
          and output_stem(TOWER_PAIR) == "ocv_haltern_kdlo"
          and figure_stem(TOWER_PAIR) == "fig_sds_ocv_reconfigure_haltern_kdlo",
          output_stem(TOWER_PAIR))

    # 14. the two records load as the registry says they do
    kdlo_units, kdlo_audit = load_units("kdlo")
    haltern_units, haltern_audit = load_units("haltern")
    check("KDLO's excluded acquisition is counted and not loaded",
          len(kdlo_units) == 59 and kdlo_audit["excluded"] == {"single": 1}
          and channels_present(kdlo_units) == ("DESC",),
          f"{len(kdlo_units)} units, {kdlo_audit['excluded']}, "
          f"{channels_present(kdlo_units)}")
    check("Haltern loads both channels and drops no unit",
          len(haltern_units) == 107
          and channels_present(haltern_units) == ("ASC", "DESC")
          and not any(haltern_audit["dropped"].values()),
          f"{len(haltern_units)} units, {haltern_audit['dropped']}")

    # 15. channels come from the units, not from the registry's expectation
    one_track = [dict(unit, channel="DESC") for unit in _synthetic("kdlo", n=2)]
    check("a record's channels are read off its units",
          channels_present(one_track) == ("DESC",)
          and channels_present(one_track + [dict(one_track[0], channel="ASC")])
          == ("ASC", "DESC"), str(channels_present(one_track)))

    # 16. an absent channel carries no caveat cell, and an empty list no sentence
    single = {"sites": ["kdlo"],
              "experiments": {"kdlo": {"channels_present": ["DESC"],
                                       "channels": {"DESC": {"blocks": {
                                           name: {"n_by_level": [2, 9]}
                                           for name in BLOCKS}}}}}}
    caveats = bootstrap_caveats(single)
    check("an absent channel contributes no caveat cells",
          {cell["channel"] for cell in caveats["cells"]} == {"DESC"}
          and len(caveats["cells"]) == 2, str(caveats["cells"]))
    check("an empty caveat list yields no sentence",
          _caveat_sentence({"bootstrap_caveats": {"block": BLOCK_BOOT,
                                                  "small_n": SMALL_N,
                                                  "cells": []}}) == "")

    # 17. a pair's cosine runs only where both records carry the channel
    two_records = {"experiments": {"haltern": {"channels_present": ["ASC", "DESC"]},
                                   "kdlo": {"channels_present": ["DESC"]}},
                   "hypotheses": {"H2": {"DESC": {}}}}
    check("a pair's cosine runs on the channels both records carry",
          shared_channels(two_records, TOWER_PAIR) == ("DESC",)
          and comparison_channels(two_records, TOWER_PAIR) == ("DESC",),
          str(shared_channels(two_records, TOWER_PAIR)))

    # 18. the plate's arithmetic: the slots, the gutters and the gap are the canvas
    for rows in (1, 2, 4):
        plate = plate_layout(rows)
        filled = (PLATE_LEFT_MM + plate["columns"] * plate["slot_mm"]
                  + (plate["columns"] - 1) * PLATE_COL_GAP_MM + PLATE_RIGHT_MM)
        check(f"a {rows}-row plate's slots fill the canvas exactly",
              abs(filled - style.FINAL_WIDTH_MM) < 1e-9, f"{filled!r}")
        check(f"a {rows}-row plate's plane keeps its units",
              abs(plate["slot_mm"] / plate["x_stretch"] - plate["plane_mm"]) < 1e-9
              and plate["x_stretch"] > 1.0,
              f"{plate['slot_mm']!r} / {plate['x_stretch']!r}")
        check(f"a {rows}-row plate's height is the sum of its bands",
              abs(plate["height_mm"] - (PLATE_HEADER_MM + PLATE_FOOT_MM
                                        + rows * plate["plane_mm"]
                                        + (rows - 1) * PLATE_ROW_GAP_MM)) < 1e-9,
              f"{plate['height_mm']!r}")
        check(f"a {rows}-row plate's grid stays inside the canvas",
              plate["grid"]["left"] < plate["grid"]["right"]
              and plate["grid"]["bottom"] < plate["grid"]["top"]
              and 0.0 < plate["grid"]["wspace"] < 0.5,
              str(plate["grid"]))
    check("the typography floor is the annotation rung",
          STATE_PT == style.FS["annotation"] >= 6.0
          and all(style.FS[rung] > style.FS["annotation"]
                  for rung in ("tick", "label", "title", "suptitle")),
          f"{STATE_PT} {style.FS}")

    for line in failures:
        print(f"  SELFTEST FAIL {line}")
    print(f"  selftest: {counted['checks']} checks, {len(failures)} failure(s)")
    return 0 if not failures else 1



def verify(report):
    """Re-derive every table cell from the two records and compare with the JSON."""
    problems = []
    want_flags = set()
    for site in (report.get("sites") or DEFAULT_SITES):
        units, audit = load_units(site)
        block = responses(units, site)
        store = report["experiments"][site]
        channels = channels_present(units)
        if list(channels) != list(channels_of(report, site)):
            problems.append(f"{site}: channels_present {store.get('channels_present')} "
                            f"vs {list(channels)}")
        if sorted(store.get("channels_absent") or {}) != [
                channel for channel in ("ASC", "DESC") if channel not in channels]:
            problems.append(f"{site}: channels_absent "
                            f"{sorted(store.get('channels_absent') or {})} does not "
                            f"match the record {list(channels)}")
        if audit.get("excluded") != ((store.get("audit") or {}).get("excluded") or {}):
            problems.append(f"{site}: excluded {audit.get('excluded')} vs "
                            f"{(store.get('audit') or {}).get('excluded')}")
        for channel in channels:
            cell = store["channels"][channel]
            for state in STATES[site]:
                cells = block[channel]["states"][state]
                count = len(units_of(units, channel, [state]))
                if count != cell["states"][state]["n_units"]:
                    problems.append(f"{site}/{channel}/{state}: n {count} vs "
                                    f"{cell['states'][state]['n_units']}")
                for name, dims in BLOCKS.items():
                    value, _k = magnitude(z_of(cells, dims), dims)
                    want = cell["states"][state]["blocks"][name]["value"]
                    if (value is None) != (want is None):
                        problems.append(f"{site}/{channel}/{state}/{name}: "
                                        f"{value} vs {want}")
                    elif value is not None and abs(value - want) > 1e-9:
                        problems.append(f"{site}/{channel}/{state}/{name}: "
                                        f"{value!r} vs {want!r}")
            for name, dims in BLOCKS.items():
                refs = {d: block[channel]["dims"][d] for d in dims}
                per_unit = unit_magnitudes(units, site, channel, dims, refs)
                groups = [[x["R"] for x in per_unit if x["state"] in states]
                          for _n, states, _d in LEVELS[site]]
                if cell["blocks"][name]["n_by_level"] != [len(g) for g in groups]:
                    problems.append(f"{site}/{channel}/{name}: level sizes differ")
                total = len(units_of(units, channel,
                                     [state for _n, states, _d in LEVELS[site]
                                      for state in states]))
                if sum(len(g) for g in groups) != total:
                    problems.append(f"{site}/{channel}/{name}: the levels hold "
                                    f"{sum(len(g) for g in groups)} units, the record "
                                    f"has {total} in those states")
            # the caveat list must name exactly the cells whose n cannot fill a block
            for name in BLOCKS:
                for level, states, _dose in LEVELS[site]:
                    if len(units_of(units, channel, list(states))) <= BLOCK_BOOT:
                        want_flags.add((site, channel, name, level))
        # the Holm value must follow from the asset's own family -- two cells per
        # channel the record carries, so two for a single-track asset and four for a
        # record with both orbit directions
        family = {k: v["p_increasing"] for k, v in report["hypotheses"]["H1"].items()
                  if k.split("/")[0] == site}
        expected = 2 * len(channels)
        if len(family) != expected:
            problems.append(f"{site}: Holm family has {len(family)} tests, not "
                            f"{expected} ({len(channels)} channel(s))")
        recomputed = holm(family)
        for key, value in recomputed.items():
            stored = report["hypotheses"]["H1"][key]["p_holm"]
            if stored is None or abs(stored - value) > 1e-12:
                problems.append(f"{key}: p_holm {stored} vs {value}")
    flagged = {(cell["site"], cell["channel"], cell["block"], cell["level"])
               for cell in report.get("bootstrap_caveats", {}).get("cells") or []}
    if flagged != want_flags:
        problems.append(f"bootstrap caveats {sorted(flagged)} vs expected "
                        f"{sorted(want_flags)}")
    return problems


#: The containers every plate is written in, in the order the guide lists them.
PLATE_FORMATS = ("jpg", "pdf", "eps")


def verify_plates(report, plates, figure_md):
    """Re-check the plates against the files on disk and against the report.

    Three separate claims are tested, because each can fail on its own:

    *every container is 170 mm wide* -- measured from the raster's own DPI and
    from the PDF/EPS bounding box, not from the number the build printed;
    *the guide quotes the plate that was built* -- the JSON blocks in the
    markdown are compared field by field with the fresh audits, so a stale or
    hand-edited guide is a failure rather than a silent inconsistency;
    *the quoted numbers are the report's* -- every ``z``, ``R``, ``n`` and
    interval in the guide is looked up in the JSON report and must agree.
    """
    problems = []
    for audit in plates:
        stem = os.path.splitext(os.path.join(REPO, audit["figure"]))[0]
        for fmt in PLATE_FORMATS:
            path = f"{stem}.{fmt}"
            rel = os.path.relpath(path, REPO)
            if not os.path.exists(path):
                problems.append(f"{rel}: missing")
                continue
            try:
                info = style.measure(path)
            except Exception as exc:  # noqa: BLE001  (any parse failure is a defect)
                problems.append(f"{rel}: cannot measure ({exc})")
                continue
            width = info.get("mm")
            if width is None:
                problems.append(f"{rel}: no width in the container")
            elif abs(width - style.FINAL_WIDTH_MM) > style.WIDTH_TOL_MM:
                problems.append(f"{rel}: {width:.2f} mm, not "
                                f"{style.FINAL_WIDTH_MM:.1f} mm")
            elif abs(width - audit["native_width_mm"].get(fmt, float("nan"))) > 1e-6:
                problems.append(f"{rel}: {width:.2f} mm, the audit says "
                                f"{audit['native_width_mm'].get(fmt):.2f} mm")
    problems += _verify_guide(report, plates, figure_md)
    return problems


def _verify_guide(report, plates, figure_md):
    """The guide's audit blocks: present, current and faithful to the report."""
    rel = os.path.relpath(figure_md, REPO)
    if not os.path.exists(figure_md):
        return [f"{rel}: missing"]
    with open(figure_md) as handle:
        text = handle.read()
    quoted = []
    for chunk in re.findall(r"```json\n(.*?)\n```", text, re.S):
        try:
            quoted.append(json.loads(chunk))
        except json.JSONDecodeError as exc:
            return [f"{rel}: an audit block does not parse ({exc})"]
    if len(quoted) != len(plates):
        return [f"{rel}: {len(quoted)} audit block(s) for {len(plates)} plate(s)"]
    problems = []
    for audit, stored in zip(plates, quoted):
        name = os.path.basename(audit["figure"])
        if stored.get("schema") != audit["schema"]:
            problems.append(f"{rel}/{name}: schema {stored.get('schema')!r}")
            continue
        for field in ("figure", "target_width_mm", "native_width_mm", "layout",
                      "final_font_sizes_pt", "bounds", "eps"):
            if stored.get(field) != audit[field]:
                problems.append(f"{rel}/{name}: {field} has drifted from the plate "
                                f"that was built "
                                f"({_first_difference(stored.get(field), audit[field])})")
        if stored.get("panel_order") != audit["panel_order"]:
            problems.append(f"{rel}/{name}: the panels are not the plate's")
            continue
        problems += _verify_quoted_panels(report, name, stored)
    return problems


def _verify_quoted_panels(report, name, audit):
    """Every number the guide quotes must still be the report's number."""
    problems = []
    for key in audit["panel_order"]:
        site, channel, block_name = key.split("/")
        panel = (audit.get("panels") or {}).get(key) or {}
        store = report["experiments"][site]
        if channel not in store["channels"]:
            if not panel.get("empty"):
                problems.append(f"{name}/{key}: quoted with data, the record has none")
            continue
        if panel.get("empty"):
            problems.append(f"{name}/{key}: quoted empty, the record has the cell")
            continue
        cell = store["channels"][channel]
        for state in STATES[site]:
            mark = cell["states"][state]["blocks"][block_name]
            point = (panel.get("states") or {}).get(state)
            if point is None:
                problems.append(f"{name}/{key}/{state}: not quoted")
                continue
            for dim in panel["dims"]:
                if not _close(point["z"].get(dim), mark["z"].get(dim)):
                    problems.append(f"{name}/{key}/{state}: z({dim}) "
                                    f"{point['z'].get(dim)} vs {mark['z'].get(dim)}")
            if point["n"] != cell["states"][state].get("n_units"):
                problems.append(f"{name}/{key}/{state}: n {point['n']} vs "
                                f"{cell['states'][state].get('n_units')}")
            if not _close(point["R"], mark["value"]):
                problems.append(f"{name}/{key}/{state}: R {point['R']} vs "
                                f"{mark['value']}")
            if _num(mark["value"], 2) not in point["label"]:
                problems.append(f"{name}/{key}/{state}: the label does not quote R "
                                f"({point['label']!r})")
            if f"n={point['n']}" not in point["label"]:
                problems.append(f"{name}/{key}/{state}: the label does not quote n "
                                f"({point['label']!r})")
        strip = cell["blocks"][block_name]["bootstrap"]["ci"]
        for level, _states, _dose in LEVELS[site]:
            want = strip.get(level)
            got = (panel.get("strip") or {}).get(level)
            if (want is None) != (got is None):
                problems.append(f"{name}/{key}/{level}: the strip {got} vs {want}")
            elif want is not None and not all(_close(a, b)
                                              for a, b in zip(got, want)):
                problems.append(f"{name}/{key}/{level}: the strip {got} vs {want}")
    return problems


def _close(a, b, tol=1e-6):
    """Equal after the rounding the audit applies, or both absent."""
    if a is None or b is None:
        return a is None and b is None
    return abs(a - b) <= tol


def _first_difference(got, want, path=""):
    """A short, addressable description of where two audit values part company."""
    if isinstance(got, dict) and isinstance(want, dict):
        for key in list(want) + [k for k in got if k not in want]:
            if got.get(key) != want.get(key):
                return _first_difference(got.get(key), want.get(key),
                                         f"{path}.{key}" if path else str(key))
        return f"{path or '<root>'}: equal"
    return f"{path or '<root>'}: {got!r} != {want!r}"


def build_parser():
    parser = argparse.ArgumentParser(
        description=("OCV reconfiguration across controlled states — the "
                     "LUMO/Bautzen pair by default, the two collapsed bridges with "
                     "--sites carola,morandi, the two towers with --sites "
                     "haltern,kdlo, or a single asset"))
    parser.add_argument("--sites", default=",".join(DEFAULT_SITES),
                        help="comma-separated assets: lumo,bautzen (default), "
                             "carola,morandi (Carolabrücke and Ponte Morandi), "
                             "haltern,kdlo (Haltern and KDLO), or a single one such "
                             "as kdlo")
    parser.add_argument("--boot", type=int, default=BOOT)
    parser.add_argument("--perm", type=int, default=PERM)
    parser.add_argument("--seed", type=int, default=SEED)
    parser.add_argument("--json", default=None)
    parser.add_argument("--md", default=None)
    parser.add_argument("--csv", default=None)
    parser.add_argument("--figdir", default=None)
    parser.add_argument("--selftest", action="store_true")
    parser.add_argument("--verify", action="store_true",
                        help="re-derive every table cell and compare with the JSON")
    parser.add_argument("--quiet", action="store_true")
    return parser


def output_stem(sites):
    """The file stem of a run: a registered pair keeps the plain name."""
    sites = tuple(sites)
    if set(sites) == set(DEFAULT_SITES):
        return "ocv_reconfigure"
    return "ocv_" + "_".join(sites)


def figure_stem(sites):
    """The plate stem: `fig_sds_ocv_reconfigure`, plus the assets beyond the pair."""
    sites = tuple(sites)
    base = "fig_sds_ocv_reconfigure"
    if set(sites) == set(DEFAULT_SITES):
        return base
    return base + "_" + "_".join(sites)


def main(argv=None):
    args = build_parser().parse_args(argv)
    if args.selftest:
        return selftest()
    sites = tuple(part.strip() for part in args.sites.split(",") if part.strip())
    stem = output_stem(sites)
    committed_path = args.json or os.path.join(REPO, "data", "sds",
                                               f"{stem}.json")
    report = compute(sites=sites, boot=args.boot, perm=args.perm, seed=args.seed)
    if args.verify:
        with open(committed_path) as handle:
            committed = json.load(handle)
        problems = verify(committed)
        for problem in problems:
            print(f"  VERIFY {problem}")
        print(f"  verify: {len(problems)} problem(s)")
        if problems:
            return 1
    json_path, md_path, csv_path = write_outputs(report, args.json, args.md,
                                                 args.csv, quiet=args.quiet)
    plates = make_figure(report, args.figdir)
    figure_md_name = ("fig_sds_ocv_reconfigure.md" if stem == "ocv_reconfigure"
                      else f"fig_sds_ocv_reconfigure_{'_'.join(sites)}.md")
    figure_md = os.path.join(os.path.dirname(
        os.path.join(REPO, plates[0]["figure"])), figure_md_name)
    with open(figure_md, "w") as handle:
        handle.write(figure_markdown(report, plates))
    if args.verify:
        problems = verify_plates(report, plates, figure_md)
        for problem in problems:
            print(f"  VERIFY {problem}")
        print(f"  verify: {len(problems)} plate problem(s)")
        if problems:
            return 1
    if not args.quiet:
        for plate in plates:
            print(f"  {plate['figure']}")
        print(f"  {os.path.relpath(figure_md, REPO)}")
        for site, block in report["reading"]["per_asset"].items():
            for name, verdict in block["verdicts"].items():
                print(f"  {EXPERIMENT[site]} / {BLOCK_DISPLAY[name]}: {verdict}")
        print(f"  comparison: {report['reading']['comparison_role']}")
        for item in report["predictions"]:
            print(f"  {item['key']}: {item['verdict']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
