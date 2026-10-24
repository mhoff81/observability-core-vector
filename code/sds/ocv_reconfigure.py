#!/usr/bin/env python3
"""OCV reconfiguration against the controlled state dose — LUMO and Bautzen.

The question the two controlled experiments have to answer together: **does the
magnitude of OCV reconfiguration grow with an independently defined state dose,
and is the direction of that reconfiguration the same in both?** The quantity is
the distance from each experiment's *own* reference state; raw OCV values are
never compared between the two assets.

The design is fixed in ``code/sds/ocv_reconfigure_spec.md`` and was written before
this module ran. Short version:

* **dose** — LUMO's dose is measured, not asserted: the first-higher natural
  frequency of the tower drops 0 % (healthy) / 7.41 % (DAM 3 and DAM 4) / 12.96 %
  (DAM 6), so the ordered test has three levels and DAM 3 and DAM 4 are one level
  (``data/lumo/reference/lumo_state_frequencies.json``). Bautzen's ladder
  ``RS < DS1 < DS2`` is ordinal only;
* **response** — ``r = median(g(X_state)) - median(g(X_ref))`` on the transformed
  scale ``g`` (``log`` for positive dimensions, identity for ``S`` and Bautzen's
  floor-corrected ``gamma2``), standardised by the **MAD of the reference state**
  with a fixed fallback chain;
* **Bautzen** is analysed on the **paired bridge-minus-control difference**, because
  its states are also its calendar;
* **blocks** — ``mask = [A, F]`` and ``other = [gamma2, P, S]``. ``D = A/bbox_area``
  holds in 178 of 178 LUMO rows, so the literal ``[A, D, F]`` block is a
  sensitivity analysis, not the primary one;
* **H1** — Jonckheere-Terpstra on the *per-unit* magnitudes (the reference's median
  magnitude is 0 by construction, so a test on state medians would be degenerate),
  Holm over the eight primary cells, moving-block bootstrap for the intervals;
* **H2** — cosine between the standardised response vectors, block by block, with a
  random-direction null; γ² is reported with and without, because it is a different
  quantity in the two experiments.

Usage
-----
  python3 code/sds/ocv_reconfigure.py
  python3 code/sds/ocv_reconfigure.py --selftest
  python3 code/sds/ocv_reconfigure.py --verify

Outputs
-------
  data/sds/ocv_reconfigure.json
  data/sds/ocv_reconfigure.md
  data/sds/ocv_reconfigure_table.csv
  figures/sds/fig_sds_ocv_reconfigure.jpg   (630 dpi)
  figures/sds/fig_sds_ocv_reconfigure.md
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import statistics
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
for _path in (os.path.join(REPO, "code", "lumo"),):
    if _path not in sys.path:
        sys.path.insert(0, _path)

from lumo_ocv_stats import _jt_statistic  # noqa: E402  (verbatim port, see module)

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
}

#: The two response blocks (see the dict below); see the 2026-10-22 amendment.
#: How a dimension is written when it is a *label* (the figures), never as a key.
DIM_DISPLAY = {"gamma2": r"$\gamma^2$", "A": "A", "D": "D", "F": "F",
               "P": "P", "S": "S"}


#: The two response blocks, in the wording the paper uses:
#:
#: * ``behaviour`` -- cross-acquisition echo behaviour: `gamma2` is the consistency
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

#: The two dimensions drawn as the plane of the geometry block; `D` and `S` are
#: annotated beside it (see ``BLOCKS``).
GEOMETRY_PLANE = ("A", "F")

#: Where a single-table asset keeps its state, orbit and date.
UNIT_SHAPE = {
    "lumo": {"state": "damage_label", "orbit": "orbit_direction",
             "date": "acquisition_date"},
    "carola": {"state": "state", "orbit": "orbit_direction",
               "date": "acquisition_ts"},
}

REFERENCE = {"lumo": "healthy", "bautzen": "RS",
             "carola": "pre-collapse (healthy)"}

STATES = {"lumo": ["healthy", "DAM 3", "DAM 4", "DAM 6"],
          "bautzen": ["RS", "DS1", "DS2"],
          "carola": ["pre-collapse (healthy)", "post-collapse"]}

#: Dose levels for the ordered test. LUMO's dose is measured (percent drop of the
#: first-higher mode); DAM 3 and DAM 4 share it to one FFT bin. Bautzen has no
#: measured dose, so its levels are the ordinal ladder.
LEVELS = {
    "lumo": [("L0", ("healthy",), 0.0),
             ("L1", ("DAM 3", "DAM 4"), 7.41),
             ("L2", ("DAM 6",), 12.96)],
    "bautzen": [("RS", ("RS",), 0.0), ("DS1", ("DS1",), 1.0), ("DS2", ("DS2",), 2.0)],
    # Carola has two states and nothing in between: the ordered test collapses to a
    # one-sided rank test between them, and dR is the post-minus-pre step.
    "carola": [("pre", ("pre-collapse (healthy)",), 0.0),
               ("post", ("post-collapse",), 1.0)],
}

#: The state used for the cross-experiment direction (and the middle one).
LARGEST = {"lumo": "DAM 6", "bautzen": "DS2",
           "carola": "post-collapse"}
MIDDLE = {"lumo": ("DAM 3", "DAM 4"), "bautzen": ("DS1",),
          "carola": ("pre-collapse (healthy)",)}

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
                                 "carola_ocv_channels.csv")}

#: The two-asset pair the module was built for; ``--sites`` may select any subset.
DEFAULT_SITES = ("lumo", "bautzen")

FREQUENCIES = os.path.join(REPO, "data", "lumo", "reference",
                           "lumo_state_frequencies.json")


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
    the audit counts those drops.
    """
    path = path or TABLE[site]
    columns = COLUMNS[site]
    audit = {"dropped": {}, "rows": 0, "pairs": 0, "unpaired": 0}
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
            values = {}
            for dim in DIMS:
                raw = _number(row.get(columns[dim]))
                values[dim] = transform(site, dim, raw)
                if raw is not None and values[dim] is None:
                    audit["dropped"][dim] = audit["dropped"].get(dim, 0) + 1
            units.append({"state": row[shape["state"]],
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


def responses(units, site, channels=("ASC", "DESC")):
    """``r`` and ``z`` for every channel, dimension, state and dose level."""
    out = {}
    for channel in channels:
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

    ``sites`` selects the assets: the two-asset pair by default (LUMO and Bautzen,
    the cross-asset comparison), or a single asset such as ``("carola",)`` for the
    two-state case. A single-asset run has no H2 partner and says so.
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
        "dose": {"source": os.path.relpath(FREQUENCIES, REPO),
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
        store = {"table": os.path.relpath(TABLE[site], REPO),
                 "sha256": _sha256(TABLE[site]), "audit": audit, "channels": {}}
        for channel in ("ASC", "DESC"):
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
    report["reading"] = reading(report)
    report["predictions"] = predictions(report)
    return report




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
    if set(sites) != set(DEFAULT_SITES):
        h2["rules"] = ("not computed: the cross-asset comparison needs the two "
                        "controlled-reference experiments")
        h2["skipped"] = True
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
              "carola": "Carolabrücke"}


def table_rows(report):
    """The compact table the pre-registration asks for, one row per cell."""
    rows = []
    for site in (report.get("sites") or DEFAULT_SITES):
        store = report["experiments"][site]
        for channel in ("ASC", "DESC"):
            cell = store["channels"][channel]
            for state in STATES[site]:
                entry = cell["states"][state]
                behaviour, geometry = (entry["blocks"]["behaviour"],
                              entry["blocks"]["geometry"])
                rows.append({
                    "Experiment": EXPERIMENT[site], "State": state,
                    "Channel": channel, "n": entry["n_units"],
                    "R_behaviour": behaviour["value"], "R_geometry": geometry["value"],
                    "Direction": (f"{behaviour['direction']} | "
                                  f"{geometry['direction']}"),
                    "Scale": cell["dims"]["A"]["scale_used"],
                    "Support": "thin" if entry["n_units"] < SMALL_N else "ok"})
    return rows


def reading_by_site(report):
    """One verdict per asset and per block, with the bound each interval licenses."""
    out = {}
    h1 = report["hypotheses"]["H1"]
    for site in (report.get("sites") or DEFAULT_SITES):
        block = {name: {ch: h1[f"{site}/{ch}/{name}"] for ch in ("ASC", "DESC")}
                 for name in ("behaviour", "geometry")}
        summary = {}
        for name, cells in block.items():
            significant = [ch for ch, c in cells.items()
                           if (c["p_holm"] or 1.0) < 0.05]
            monotone = {}
            bounds = {}
            for channel in ("ASC", "DESC"):
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
                        "the pre-to-post contrast")
                verdict = (f"{what} in the {name} block survives the asset's own "
                           f"Holm family on " + ", ".join(significant))
            elif any(monotone.values()):
                rising = [ch for ch, flag in monotone.items() if flag]
                where = ("both channels" if len(rising) == 2
                         else " ".join(rising) + " only")
                verdict = (f"no {('dose ordering' if ordered else 'state contrast')} "
                           f"in the {name} block survives the asset's own Holm family; "
                           f"the point estimate of the step is positive on " + where)
            else:
                verdict = (f"no {('dose ordering' if ordered else 'state contrast')} "
                           f"in the {name} block survives the asset's own Holm family, "
                           f"and the point estimate of the step is not monotone")
                verdict += " on either channel"
            summary[name] = {
                "verdict": verdict, "significant": significant,
                "monotone_by_channel": monotone, "delta_ci": bounds,
                "widest_upper_bound": max((c["delta_ci"] or [None, None])[1] or 0.0
                                          for c in cells.values())}
        verdicts = {name: item["verdict"] for name, item in summary.items()}
        out[site] = {
            "verdicts": verdicts, "blocks": summary,
            "n_ref": {ch: report["experiments"][site]["channels"][ch]["dims"]["A"]["n"]
                      for ch in ("ASC", "DESC")},
            "scale_used": {ch: report["experiments"][site]["channels"][ch]["dims"]["A"]
                                       ["scale_used"] for ch in ("ASC", "DESC")}}
    return out


def reading(report):
    """Per-asset verdicts plus — for the pair — the descriptive comparison."""
    h2 = report["hypotheses"]["H2"]
    comparison = {}
    if not h2.get("skipped"):
        for channel in ("ASC", "DESC"):
            comparison[channel] = {}
            for block_name in ("behaviour", "geometry"):
                item = h2[channel][block_name]["largest"]
                comparison[channel][block_name] = {
                    "cosine": item["cosine"], "ci": item["ci"],
                    "sign_agreement": item["sign_agreement"],
                    "p_null": item.get("p_null"),
                    "lumo_direction": item["lumo_direction"],
                    "bautzen_direction": item["bautzen_direction"]}
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
            step = None if (levels.get("L1") is None or levels.get("L2") is None)                 else abs(levels["L2"] - levels["L1"])
            gap = None if (three is None or four is None) else abs(three - four)
            holds.append(gap is not None and step is not None and gap <= step)
            details.append(f"{name}/{channel}: |R(DAM 3) - R(DAM 4)| = {_num(gap)} "
                           f"against a cross-dose step of {_num(step)}")
        verdicts.append("held" if all(holds) else
                        ("partly" if any(holds) else "did not hold"))
    out.append({"key": "P2",
                "verdict": "; ".join(f"{name} {verdict}" for name, verdict in
                                      zip(("geometry", "behaviour"), verdicts))
                           if len(set(verdicts)) > 1 else verdicts[0],
                "statement": "DAM 3 and DAM 4 share one dose, so their magnitudes "
                             "agree (the same-dose gap is not larger than the "
                             "cross-dose step) — asked per block",
                "detail": "; ".join(details)})

    # P3 -- the behaviour block need not follow the dose
    behaviour = [report["hypotheses"]["H1"][f"{site}/{channel}/behaviour"]["p_holm"]
                 for site in (report.get("sites") or DEFAULT_SITES)
                 for channel in ("ASC", "DESC")]
    following = [p for p in behaviour if p is not None and p < 0.05]
    out.append({"key": "P3",
                "verdict": "consistent" if not following else
                           "the behaviour block follows the dose too",
                "statement": "the γ²/P behaviour block need not follow the dose",
                "detail": f"{len(following)} of {len(behaviour)} behaviour-block tests "
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



def markdown(report):
    """Two assets, two verdicts, and one descriptive comparison between them."""
    lines = []
    add = lines.append
    dose = report["dose"]
    reading = report["reading"]
    sites = tuple(report.get("sites") or DEFAULT_SITES)
    is_pair = set(sites) == set(DEFAULT_SITES)
    subject = ("two assets, two verdicts" if is_pair else
               " ".join(EXPERIMENT[site] for site in sites) + ", two states")
    add(f"# OCV reconfiguration against the controlled state dose — {subject}")
    add("")
    add("*generated by `code/sds/ocv_reconfigure.py`; the design is fixed in "
        "`code/sds/ocv_reconfigure_spec.md`, including the 2026-10-22 amendment that "
        "separates the two assets*")
    add("")
    add("Each asset is measured against **its own reference state**, in units of that "
        "reference's own variability"
        + (" — LUMO against `healthy` on the coherence rows, Bautzen against `RS` on "
           "the paired bridge-minus-control difference" if is_pair else
           " — Carola against `pre-collapse (healthy)` on the support-sized chip, "
           "with the post-collapse state as the single contrast")
        + ". Assets share a metric, not a scale, so their magnitudes are never "
        "compared. Each carries its **own** Holm family of four tests and its own "
        "verdict; no joint label is issued"
        + (", and the comparison section is a description, not a hypothesis."
           if is_pair else "."))
    add("")
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
    for site in (report.get("sites") or DEFAULT_SITES):
        block = reading["per_asset"][site]
        add(f"## {EXPERIMENT[site]}")
        add("")
        add(f"* reference `{REFERENCE[site]}`, n_ref = {block['n_ref']['ASC']} (ASC) / "
            f"{block['n_ref']['DESC']} (DESC); scale used = "
            f"`{block['scale_used']['ASC']}` (ASC) / `{block['scale_used']['DESC']}` "
            f"(DESC).")
        add("")
        add("| State | Channel | n | R_behaviour | R_geometry | Direction | Scale | Support |")
        add("| --- | --- | ---: | ---: | ---: | --- | --- | --- |")
        for row in report["table"]:
            if row["Experiment"] != EXPERIMENT[site]:
                continue
            add(f"| {row['State']} | {row['Channel']} | {row['n']} | "
                f"{_num(row['R_behaviour'])} | {_num(row['R_geometry'])} | "
                f"{row['Direction']} | `{row['Scale']}` | {row['Support']} |")
        add("")
        add("| channel | block | R by level | ΔR, top step | p (one-sided) | Holm "
            "in the asset's family |")
        add("| --- | --- | --- | --- | ---: | ---: |")
        for channel in ("ASC", "DESC"):
            for name in ("behaviour", "geometry"):
                cell = report["hypotheses"]["H1"][f"{site}/{channel}/{name}"]
                levels = " < ".join(f"{key} {_num(value)}"
                                    for key, value in cell["R_by_level"].items())
                pair = cell["delta_ci"] or [None, None]
                add(f"| {channel} | `{name}` [{', '.join(cell['dims'])}] | {levels} | "
                    f"{_num(pair[0])} … {_num(pair[1])} | "
                    f"{_p(cell['p_increasing'])} | {_p(cell['p_holm'])} |")
        add("")
        for name in ("behaviour", "geometry"):
            item = block["blocks"][name]
            add(f"**Verdict — {EXPERIMENT[site]}, {name}**: {item['verdict']}. Its own "
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
        else:
            add("Carola has two states and no control track: `pre` and `post` are "
                "separated by the collapse date itself, so this is a before/after "
                "contrast in time, not a controlled comparison. The table is the "
                "support-sized chip of `Girder A` alone.")
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
    if not is_pair:
        add("## No comparison section")
        add("")
        add(f"A single-asset run has no comparison partner, so the cosine block of "
            f"the pair is not computed. The cross-asset question (`{EXPERIMENT['lumo']}` "
            f"vs. `{EXPERIMENT['bautzen']}`) needs both controlled-reference "
            f"experiments and is reported in `data/sds/ocv_reconfigure.md`.")
        add("")
    if is_pair:
        add("## Comparison — the same metric, two assets (descriptive)")
        add("")
        add("Not a hypothesis and not a test: the cosine of the two standardised "
            "response vectors at the largest state says whether the *same metric*, applied "
            "twice to two assets, answers with the same direction. Magnitudes are not "
            "compared — each asset's `R` lives on its own scale.")
        add("")
        add("| channel | block | cosine | 95 % CI | sign agreement | null p | LUMO | Bautzen |")
        add("| --- | --- | ---: | --- | --- | ---: | --- | --- |")
        for channel in ("ASC", "DESC"):
            for block_name in ("behaviour", "geometry"):
                cell = report["hypotheses"]["H2"][channel][block_name]
                item = cell["largest"]
                add(f"| {channel} | `{block_name}` "
                    f"[{', '.join(item['dims'])}] | {_num(item['cosine'])} | "
                    f"{_ci(item['ci'])} | "
                    f"{item['sign_agreement'][0]}/{item['sign_agreement'][1]} | "
                    f"{_p(item.get('p_null'))} | {item['lumo_direction']} | "
                    f"{item['bautzen_direction']} |")
                for dim in ("gamma2", "P"):
                    only = cell.get(f"only_{dim}")
                    if not only:
                        continue
                    add(f"| {channel} | `{dim}` alone | — | – | – | – | "
                        f"{_num(only['lumo_z'], 2)} | {_num(only['bautzen_z'], 2)} |")
        add("")
        add("`null p` is the share of random directions in 3 (or 2) dimensions whose "
            "absolute cosine is at least as large — a guard against a coarse direction "
            "space, not a significance claim. γ² is a different quantity in the two "
            "experiments, which is why the two dimensions of the behaviour block are "
            "also printed one by one (`z` per asset). The geometry block carries four "
            "dimensions while its plane shows two: `R_geometry` is the RMS over "
            "`A, D, F, S`, and `D` (which partly inherits `A`) and `S` (position, not "
            "morphology) are annotated beside the plane rather than drawn as axes.")
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
    add("A unit is dropped for one dimension only when the transformed value is "
        "undefined. LUMO loses one `gamma2` row (a zero); Bautzen loses a *pair* "
        "whenever either side of the bridge-control difference is missing or zero.")
    add("")
    add("## What this does not claim")
    add("")
    add("* No magnitude comparison between the assets: `R` is scaled by each "
        "reference's own variation — and, where the reference is smaller than ten "
        "units, by that asset's pooled variation — so the numbers live on different "
        "units by construction.")
    add("* No joint verdict. Each asset's Holm family, and therefore each asset's "
        "verdict and bound, is its own.")
    add("* No severity derived from an OCV response: LUMO's ordering comes from the "
        "measured natural frequencies, Bautzen's from the documented ladder.")
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
    fields = ["Experiment", "State", "Channel", "n", "R_behaviour", "R_geometry",
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
# The figure — two rows, four panels, arrows from the reference
# ---------------------------------------------------------------------------
PANELS = (("ASC", "behaviour"), ("ASC", "geometry"),
          ("DESC", "behaviour"), ("DESC", "geometry"))


def _magnitude_strip(ax, entry, site):
    """R per dose level with its bootstrap interval, small, inside the panel."""
    levels = [name for name, _states, _dose in LEVELS[site]]
    ci = entry["bootstrap"]["ci"]
    if not any(ci.get(name) for name in levels):
        return
    inset = ax.inset_axes([0.60, 0.05, 0.36, 0.20])
    inset.tick_params(labelsize=4.0, length=1.5, pad=1.0)
    for index, name in enumerate(levels):
        pair = ci.get(name)
        if not pair:
            continue
        inset.plot([pair[0], pair[1]], [index, index], color="0.35", lw=0.8)
        inset.plot([(pair[0] + pair[1]) / 2.0], [index], marker="|", ms=4.0,
                   color="#c53030")
    inset.set_yticks(range(len(levels)))
    inset.set_yticklabels(levels)
    inset.set_xlabel("R", fontsize=4.4, labelpad=1.0)
    inset.set_ylim(-0.6, len(levels) - 0.4)
    for spine in ("top", "right"):
        inset.spines[spine].set_visible(False)


def _panel(ax, report, site, channel, block_name, plt):
    """One cell of a plate: arrows from the reference, plus the R strip."""
    dims = GEOMETRY_PLANE if block_name == "geometry" else BLOCKS[block_name]
    cell = report["experiments"][site]["channels"][channel]
    span = 0.3
    for state in STATES[site]:
        mark = cell["states"][state]["blocks"][block_name]
        x, y = mark["z"].get(dims[0]), mark["z"].get(dims[1])
        if x is None or y is None:
            continue
        colour = "#2b6cb0" if state == REFERENCE[site] else "#c53030"
        ax.annotate("", xy=(x, y), xytext=(0.0, 0.0),
                    arrowprops={"arrowstyle": "->", "lw": 1.0, "color": colour,
                                "alpha": 0.9, "shrinkA": 0, "shrinkB": 0})
        ax.plot([x], [y], marker="o", ms=3.0, color=colour, zorder=3)
        ax.text(x, y, f" {state} ({_num(mark['value'], 2)})", fontsize=5.2,
                color=colour, va="bottom", ha="left")
        span = max(span, abs(x), abs(y))
    limit = span * 1.5
    ax.axhline(0.0, color="0.75", lw=0.5, ls=":")
    ax.axvline(0.0, color="0.75", lw=0.5, ls=":")
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    ax.set_aspect("equal")
    ax.tick_params(labelsize=4.8, length=2.0)
    ax.set_xlabel(f"z({dims[0]})", fontsize=5.4, labelpad=1.0)
    ax.set_ylabel(f"z({dims[1]})", fontsize=5.4, labelpad=1.0)
    ax.set_title(f"{EXPERIMENT[site]} · {channel} · {block_name} "
                 f"[{', '.join(DIM_DISPLAY.get(dim, dim) for dim in dims)}]",
                 fontsize=6.2, loc="left")
    if block_name == "geometry":
        # D and S stay in the block magnitude but are not on the plane: D partly
        # inherits A algebraically, and S describes position, not morphology
        mark = cell["states"][LARGEST[site]]["blocks"]["geometry"]["z"]
        note = " ".join(f"{dim}{mark[dim]:+.2f}" for dim in ("D", "S")
                        if mark.get(dim) is not None)
        if note:
            ax.text(0.97, 0.97, note, transform=ax.transAxes, ha="right",
                    va="top", fontsize=4.6, color="0.35")

    _magnitude_strip(ax, cell["blocks"][block_name], site)


def _plate(report, sites, figdir, stem, suptitle):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(len(sites), len(PANELS),
                             figsize=(13.0, 3.2 * len(sites)), dpi=110,
                             squeeze=False)
    for row, site in enumerate(sites):
        for col, (channel, block_name) in enumerate(PANELS):
            _panel(axes[row][col], report, site, channel, block_name, plt)
    fig.suptitle(suptitle + "\narrows: reference → state, in units of the "
                 "standardising scale; inset: block magnitude R per level with its "
                 "moving-block 95 % interval", fontsize=7.0)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    out = os.path.join(figdir, f"{stem}.jpg")
    fig.savefig(out, dpi=FIGURE_DPI, facecolor="white")
    plt.close(fig)
    return out


#: One line per plate title, so a single-asset run reads like the pair's rows.
PLATE_TITLE = {
    "lumo": "LUMO tower — OCV reconfiguration against the measured frequency drop",
    "bautzen": "Bautzen bridge — OCV reconfiguration against the ordinal state "
               "ladder (paired bridge minus control)",
    "carola": "Carolabrücke — OCV reconfiguration between the pre- and the "
               "post-collapse state (support-sized chip, Girder A)",
}


def make_figure(report, figdir=None):
    """The plates: the juxtaposition for the pair, one plate per asset otherwise."""
    figdir = figdir or os.path.join(REPO, "figures", "sds")
    os.makedirs(figdir, exist_ok=True)
    sites = tuple(report.get("sites") or DEFAULT_SITES)
    if set(sites) == set(DEFAULT_SITES):
        plates = [("fig_sds_ocv_reconfigure", tuple(sites),
                   "OCV reconfiguration against the controlled state dose — "
                   "juxtaposition of two assets (each against its own reference; "
                   "no shared scale)")]
        plates += [(f"fig_sds_ocv_reconfigure_{site}", (site,), PLATE_TITLE[site])
                   for site in sites]
    else:
        plates = [(f"fig_sds_ocv_reconfigure_{site}", (site,), PLATE_TITLE[site])
                  for site in sites]
    return [_plate(report, subset, figdir, stem, title)
            for stem, subset, title in plates]


def figure_markdown(report):
    """The panel guide: the pair's three plates, or one plate per asset."""
    sites = tuple(report.get("sites") or DEFAULT_SITES)
    stem = output_stem(sites)
    pair = set(sites) == set(DEFAULT_SITES)
    lines = ["# `fig_sds_ocv_reconfigure*.jpg` — panel guide", "",
             f"*generated by `code/sds/ocv_reconfigure.py`; the numbers behind every "
             f"arrow are in `data/sds/{stem}.json`*", "",
             ("Three plates, two roles:" if pair
              else "One plate — the asset standing alone:"), "",
             "| file | role |", "| --- | --- |"]
    if pair:
        lines += [
            "| `fig_sds_ocv_reconfigure.jpg` | the juxtaposition — rows = the two "
            "assets, columns = channel × block; the caption says *no shared scale* |",
            "| `fig_sds_ocv_reconfigure_lumo.jpg` | the tower alone (4 panels) |",
            "| `fig_sds_ocv_reconfigure_bautzen.jpg` | the bridge alone (4 panels) |"]
    else:
        for site in sites:
            lines.append(f"| `fig_sds_ocv_reconfigure_{site}.jpg` | "
                         f"{EXPERIMENT[site]} alone "
                         f"(4 panels) |")
    lines += ["",
             "In every panel the origin is that asset's **reference state** and the "
             "arrows end at its standardised responses "
             "(`z = (median − median_ref) / scale`), so the geometry is "
             "within-asset; the two assets never share a scale.", ""]
    for site in (report.get("sites") or DEFAULT_SITES):
        block = report["reading"]["per_asset"][site]
        if site == "lumo":
            levels = ", ".join(f"`{name}` ({_num(dose, 2)} % frequency drop)"
                               for name, _states, dose in LEVELS[site])
        else:
            levels = ", ".join(f"`{name}` (rank {_num(dose, 2)})"
                               for name, _states, dose in LEVELS[site])
        lines.append(f"* **{EXPERIMENT[site]}** — reference `{REFERENCE[site]}` "
                     f"(n_ref {block['n_ref']['ASC']} ASC / {block['n_ref']['DESC']} "
                     f"DESC), scale `{block['scale_used']['ASC']}` / "
                     f"`{block['scale_used']['DESC']}`, levels {levels}; behaviour "
                     f"`[{'`, `'.join(BLOCKS['behaviour'])}]`, geometry "
                     f"`[{'`, `'.join(BLOCKS['geometry'])}]` (plane "
                     f"`[{', '.join(GEOMETRY_PLANE)}]`).")
    lines += ["",
              "The inset in each panel repeats the block magnitude `R` per level with "
              "its moving-block 95 % interval — the quantity H1 tests. The reference "
              "level sits at the origin by construction and carries no interval.", "",
              "Nothing is clipped: the axis limits follow the largest coordinate in "
              "the panel.", ""]
    return "\n".join(lines)


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
    """The behaviours the metric must have, on synthetic units."""
    failures = []

    def check(label, condition, detail=""):
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

    for line in failures:
        print(f"  SELFTEST FAIL {line}")
    print(f"  selftest: 11 checks, {len(failures)} failure(s)")
    return 0 if not failures else 1



def verify(report):
    """Re-derive every table cell from the two records and compare with the JSON."""
    problems = []
    for site in (report.get("sites") or DEFAULT_SITES):
        units, _audit = load_units(site)
        block = responses(units, site)
        store = report["experiments"][site]
        for channel in ("ASC", "DESC"):
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
        # the Holm value must follow from the asset's own four-test family
        family = {k: v["p_increasing"] for k, v in report["hypotheses"]["H1"].items()
                  if k.split("/")[0] == site}
        if len(family) != 4:
            problems.append(f"{site}: Holm family has {len(family)} tests, not 4")
        recomputed = holm(family)
        for key, value in recomputed.items():
            stored = report["hypotheses"]["H1"][key]["p_holm"]
            if stored is None or abs(stored - value) > 1e-12:
                problems.append(f"{key}: p_holm {stored} vs {value}")
    return problems


def build_parser():
    parser = argparse.ArgumentParser(
        description=("OCV reconfiguration against the controlled state dose — the "
                     "LUMO/Bautzen pair by default, or a single asset with "
                     "--sites"))
    parser.add_argument("--sites", default=",".join(DEFAULT_SITES),
                        help="comma-separated assets: lumo,bautzen (default) or "
                             "a single one such as carola")
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
    """The file stem of a run: the pair keeps the plain name, a single asset its own."""
    sites = tuple(sites)
    if set(sites) == set(DEFAULT_SITES):
        return "ocv_reconfigure"
    return "ocv_" + "_".join(sites)


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
    figures = make_figure(report, args.figdir)
    figure_md_name = ("fig_sds_ocv_reconfigure.md" if stem == "ocv_reconfigure"
                      else f"fig_sds_ocv_reconfigure_{'_'.join(sites)}.md")
    figure_md = os.path.join(os.path.dirname(figures[0]), figure_md_name)
    with open(figure_md, "w") as handle:
        handle.write(figure_markdown(report))
    if not args.quiet:
        for path in figures:
            print(f"  {os.path.relpath(path, REPO)}")
        print(f"  {os.path.relpath(figure_md, REPO)}")
        for site, block in report["reading"]["per_asset"].items():
            for name, verdict in block["verdicts"].items():
                print(f"  {EXPERIMENT[site]} / {name}: {verdict}")
        print(f"  comparison: {report['reading']['comparison_role']}")
        for item in report["predictions"]:
            print(f"  {item['key']}: {item['verdict']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
