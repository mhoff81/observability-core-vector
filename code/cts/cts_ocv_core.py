#!/usr/bin/env python3
"""CTS OCV-Paper — data layer, constants and the committed reference join.

Everything the other modules need, in one place: paths, the analysis channels,
the mask vector, the site constants, the collapse event, the loaders of the
committed tables **and the difference-in-differences reference files**.

The site
--------
Champlain Towers South (Surfside, Miami-Dade, FL, 25.87306 N 80.12083 W) is a
12-storey reinforced-concrete condominium that partially collapsed in the early
hours of 2021-06-24 (01:22 EDT = 05:22 UTC); the standing remnant was brought
down on 2021-07-05. The export stores one 80x80 complex chip per acquisition per
**track** — here four Sentinel-1 observation footprints of the same orbit
(ASC rel-48 IW3):

    cts_asc_48_iw3   the tower itself          role ``target``
    ctn_asc_48_iw3   Champlain Towers North    role ``control_ctn``  (165 m N)
    cte_asc_48_iw3   Champlain Towers East      role ``control_cte``  (~90 m N)
    beach_asc_48_iw3 the beach east of the site  role ``reference_beach``

see ``data/cts/cts_tracks.txt``. The collapse divides the record into the two
states

    pre-collapse (healthy)   date <  2021-06-24
    post-collapse            date >= 2021-06-24

Unlike the single-structure sites, the three towers share identical strata,
dates and counts, so the package can do what a single site cannot: a
**difference-in-differences** (DiD) that contrasts the *target* against a
*control* tower, both of which see the same atmosphere and the same orbit. That
layer is first class here (``cts_ocv_did.py``) and is pinned against the three
committed ``cts_reference_did_*.json`` files.

The post state is **bare sand**: its intensity is moisture-driven, so a change
in the cleared-lot window is *not* evidence about the collapse.

The two mask layers
-------------------
  * ``echo mask`` — the whole chip (the deck as the window sees it), the layer
    of the 6D vector ``x = [gamma2, P, D, A, F, S]``;
  * ``per-track deck mask`` — the same six dimensions, but per track
    (``GIRDERS``), computed on that track's own chip.

Both come from ``cts_ocv_masks`` and therefore from the committed window cache;
the per-date measurements are carried from ``cts_measurements_full.txt``.

Reference files
---------------
  data/cts/reference/cts_reference_did_ctn.json             target vs. CTN
  data/cts/reference/cts_reference_did_cte.json             target vs. CTE
  data/cts/reference/cts_reference_did_placebo_ctn_vs_cte.json  placebo CTN vs. CTE
and the provenance inputs ``cts_measurements_full.txt``, ``cts_tracks.txt``,
``cts_segments.txt``, ``cts_events.json``, ``cts_asset.json``,
``cts_coherence_states.json``, ``cts_acquisition_census.{json,md}`` and
``cts_phase0_audit.json``.
"""
from __future__ import annotations

import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import cts_ocv_masks as masks  # noqa: E402
import cts_ocv_stats as st     # noqa: E402

DATA = os.path.join(HERE, os.pardir, os.pardir, "data", "cts")
FIGDIR = os.path.join(HERE, os.pardir, os.pardir, "figures", "cts")
REFERENCE = os.path.join(DATA, "reference")

CACHE_PATH = masks.CACHE_PATH
MANIFEST_PATH = masks.MANIFEST_PATH
MEAS_PATH = masks.MEAS_PATH
TRACKS_PATH = masks.TRACKS_PATH
# A "segment" is a track here (the girder-equivalent), exactly as in the Morandi
# package; ``cts_segments.txt`` is carried as provenance only.
SEGMENTS_PATH = TRACKS_PATH
REF_SEGMENTS = os.path.join(DATA, "cts_segments.txt")

CSV_PATH = os.path.join(DATA, "cts_ocv_channels.csv")
CSV_META_PATH = os.path.join(DATA, "cts_ocv_channels_meta.json")

# The committed difference-in-differences references (the CTS-specific layer).
REF_DID = {
    "did_ctn": os.path.join(REFERENCE, "cts_reference_did_ctn.json"),
    "did_cte": os.path.join(REFERENCE, "cts_reference_did_cte.json"),
    "did_placebo_ctn_vs_cte":
        os.path.join(REFERENCE, "cts_reference_did_placebo_ctn_vs_cte.json"),
}
# The committed provenance inputs of the site (Step 0 / Phase 0 chain).
REF_CENSUS = os.path.join(DATA, "cts_acquisition_census.json")
REF_EVENTS = os.path.join(DATA, "cts_events.json")
REF_ASSET = os.path.join(DATA, "cts_asset.json")
REF_COHERENCE_STATES = os.path.join(DATA, "cts_coherence_states.json")
REF_PHASE0 = os.path.join(DATA, "cts_phase0_audit.json")

SITE = "Champlain Towers South, Surfside (FL)"
SITE_SHORT = "CTS"
COLLAPSE_DATE = "2021-06-24"
DEMOLITION_DATE = "2021-07-05"
EVENT = COLLAPSE_DATE

STATES = list(st.STATE_ORDER)
PRE, POST = STATES
STATE_SHORT = dict(st.STATE_SHORT)

OCV = ["gamma2", "P", "D"]
FEATURES = ["gamma2", "P", "D", "A", "F", "S"]
DIMS = list(FEATURES)

# ---------------------------------------------------------------------------
# Track (the girder-equivalent) and role tables
# ---------------------------------------------------------------------------
# The four Sentinel-1 footprints of the record; ``GIRDERS`` keeps the LUMO name
# (and ``by_girder`` / ``by_track`` the helper names) so the figure scripts read
# identically to the other sites.
GIRDERS = ["cts_asc_48_iw3", "ctn_asc_48_iw3", "cte_asc_48_iw3", "beach_asc_48_iw3"]
TRACK_SHORT = {"cts_asc_48_iw3": "cts", "ctn_asc_48_iw3": "ctn",
               "cte_asc_48_iw3": "cte", "beach_asc_48_iw3": "beach"}
TRACK_ROLE = {"cts_asc_48_iw3": "target", "ctn_asc_48_iw3": "control_ctn",
              "cte_asc_48_iw3": "control_cte", "beach_asc_48_iw3": "reference_beach"}
ROLES = ["target", "control_ctn", "control_cte", "reference_beach"]

# The DiD contrast: (target role, reference role) -> the committed reference key.
DID_CONTRASTS = [
    ("target", "control_ctn", "did_ctn"),
    ("target", "control_cte", "did_cte"),
    ("control_ctn", "control_cte", "did_placebo_ctn_vs_cte"),
]

# The nine channels the site's DiD reference reports, in the committed order.
DID_CHANNELS = ["phase_coherence", "phase_snr_db", "phase_scatterer_count",
                "displacement_los_m", "coherence_gamma2", "coherence",
                "coherence_masked_pixels", "intensity", "brightness_ratio"]
DID_ORBITS = ["ASCENDING", "DESCENDING"]
# The site's degenerate channel: it is the *fixed* 640-px quantile mask, constant
# for every row, so its OLS interaction is pure floating-point noise. Pinned
# structurally (did/n/cells/deltas only), never on did_t / did_p.
DEGENERATE_CHANNELS = ["coherence_masked_pixels"]

# The site's exported fixed quantile deck-mask size (coherence_masked_pixels).
QUANTILE_MASK_PIXELS = 640

CHANNEL_AXIS = {
    "gamma2": "gamma^2  (masked coherence)",
    "P": "P  (mask persistence)",
    "D": "D  (mask density A / bbox)",
    "A": "A  (echo pixels n_masked)",
    "F": "F  (8-connected components)",
    "S": "S  (||centroid - peak||, px)",
    "temperature_c": "temperature (C)",
    "precipitation_mm": "precip (mm)",
    "shift_x": "shift x (px)",
    "shift_y": "shift y (px)",
    "intensity": "intensity (whole window)",
    "coherence_gamma2": "gamma^2 (pipeline export)",
    "phase_coherence": "phase coherence",
    "phase_scatterer_count": "phase scatterers",
    "coherence_masked_pixels": "mask pixels (fixed 640 px)",
    "brightness_ratio": "brightness ratio",
}
CHANNEL_TITLE = {
    "gamma2": "gamma2", "P": "P", "D": "D", "A": "A", "F": "F", "S": "S",
    "temperature_c": "temperature (C)", "precipitation_mm": "precip (mm)",
    "shift_x": "shift x (px)", "shift_y": "shift y (px)",
    "intensity": "intensity", "coherence_gamma2": "gamma2 (export)",
    "phase_coherence": "phase coherence",
    "phase_scatterer_count": "phase scatterers",
    "coherence_masked_pixels": "mask pixels",
    "brightness_ratio": "brightness ratio",
}

WEATHER = ["temperature_c", "precipitation_mm"]
CONTROLS = ["temperature_c", "precipitation_mm", "shift_x", "shift_y"]

MODELS = {
    "gamma2": ["gamma2"],
    "PD": ["P", "D"],
    "gamma2PD": ["gamma2", "P", "D"],
    "x_6d": ["gamma2", "P", "D", "A", "F", "S"],
    "weather": ["temperature_c", "precipitation_mm"],
}
HEADLINE_PAIR = "x_6d"
LAMBDA_GRID = [0.0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9]
HEADLINE_LAMBDA = 0.5
RNG_SEED = 7
N_BOOT = 10000
N_PERM = 1000

MASK_COLUMNS = ["A", "D", "F", "S", "P", "gamma2_mask", "coherence_masked_pixels",
                "mask_bbox_area", "mask_row_span", "mask_col_span",
                "mask_centroid_row", "mask_centroid_col", "mask_peak_row",
                "mask_peak_col", "mask_n_components", "mask_largest_component",
                "mask_fragmentation"]

# The measurement schema (``cts_measurements_full.txt``, 28 columns, 712 rows:
# 178 dates x 4 tracks). ``id`` is the pipeline key; ``date``/``acquisition_ts``
# the day; ``track``/``role`` the footprint; and the numeric channels below.
SOURCE_TEXT = ["asset_id", "asset_name", "asset_type", "request_id", "date",
               "acquisition_ts", "pass_label", "orbit_direction", "subswath",
               "platform", "state_window", "role", "track"]
SOURCE_NUMERIC = ["relative_orbit", "intensity", "coherence_gamma2", "coherence",
                  "phase_coherence", "coherence_masked_pixels",
                  "phase_scatterer_count", "brightness_ratio", "phase_snr_db",
                  "displacement_los_m", "temperature_c", "precipitation_mm",
                  "shift_y", "shift_x"]
SOURCE_COLUMNS = ["id"] + SOURCE_TEXT + SOURCE_NUMERIC
DERIVED_COLUMNS = ["segment", "track_index", "state", "pre_collapse", "day",
                   "month", "year", "season", "echo_mode"]
CSV_COLUMNS = SOURCE_COLUMNS + DERIVED_COLUMNS + MASK_COLUMNS


def _num(v):
    """CSV/JSON cell -> float or None."""
    if v is None or v == "" or v == "null" or v == "None":
        return None
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def _bool(v):
    return str(v) in ("1", "True", "true")


def day_of(row):
    """Acquisition day (YYYY-MM-DD) of a CSV row."""
    d = row.get("day")
    return d if d else (row.get("acquisition_ts") or "")[:10]


def state_of(day):
    """The collapse alphabet: ``pre-collapse (healthy)`` before the event."""
    return PRE if str(day) < COLLAPSE_DATE else POST


def echo_mode_of(a):
    """Echo-mode bucket of a masked-pixel count (the LUMO/KDLO convention)."""
    a = _num(a)
    if a is None:
        return None
    if a <= 10:
        return "compact"
    if a < 25:
        return "intermediate"
    return "distributed"


def track_index_of(track):
    """0/1/2/3 index of the track (the girder-request index analogue)."""
    try:
        return GIRDERS.index(track)
    except ValueError:
        return None


def role_of_track(track):
    """The role a track plays (``target`` / ``control_ctn`` / ...)."""
    return TRACK_ROLE.get(track)


# ---------------------------------------------------------------------------
# Loaders of the committed tables
# ---------------------------------------------------------------------------
def _read_pipe(path):
    """A pipe-delimited table (``#`` comment lines skipped) as a list of dicts."""
    out = []
    with open(path) as fh:
        header = None
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("|")
            if header is None:
                header = parts
                continue
            if len(parts) == len(header):
                out.append(dict(zip(header, parts)))
    return out


def load_csv(path=CSV_PATH):
    """The generated channel table, numeric cells as floats."""
    with open(path) as fh:
        rows = list(csv.DictReader(fh))
    for r in rows:
        for c in SOURCE_NUMERIC:
            r[c] = _num(r.get(c))
        for c in MASK_COLUMNS:
            r[c] = _num(r.get(c))
        r["gamma2"] = r["gamma2_mask"]
        r["gamma2_export"] = r.get("coherence_gamma2")
        r["track_index"] = _num(r.get("track_index"))
        r["pre_collapse"] = _bool(r.get("pre_collapse"))
    return rows


def load_csv_meta(path=CSV_META_PATH):
    with open(path) as fh:
        return json.load(fh)


def load_cache(path=CACHE_PATH):
    """The committed window cache (the mask layer's input)."""
    return masks.load_cache(path)


def load_manifest(path=MANIFEST_PATH):
    return masks.load_manifest(path)


def load_measurements(path=MEAS_PATH):
    """The committed measurement extract, keyed ``(track, date)``.

    The table carries, per acquisition and per track, the whole-window channels
    (``intensity``, ``coherence_gamma2``, ``phase_coherence``,
    ``coherence_masked_pixels``, ``phase_scatterer_count``, ``brightness_ratio``,
    ``phase_snr_db``, ``displacement_los_m``), the registration shifts
    (``shift_x``/``shift_y``), the Open-Meteo weather (``temperature_c``,
    ``precipitation_mm``) and the site's state window (``state_window``).
    ``coherence_gamma2`` is re-exposed as ``gamma2_export`` (the pipeline's
    *own* spatial coherence, kept for provenance only — the mask layer
    recomputes its own ``gamma2``). ``coherence_masked_pixels`` is the site's
    *fixed 640-px quantile* mask and is **not** the echo mask this package
    recomputes.
    """
    out = {}
    for d in _read_pipe(path):
        if not d.get("track"):
            continue
        d["gamma2_export"] = d.get("coherence_gamma2")
        d["status"] = d.get("state_window", "")
        d["group"] = d.get("role", "")
        out[(d["track"], d["date"])] = d
    return out


def load_segments(path=SEGMENTS_PATH):
    """``track -> {label, orbit, ...}`` (the segment/track table)."""
    return masks.load_tracks(path)


def load_tracks(path=TRACKS_PATH):
    """``track -> {label, role, orbit, ...}``."""
    return masks.load_tracks(path)


def load_reference(path):
    with open(path) as fh:
        return json.load(fh)


# ---------------------------------------------------------------------------
# Small aggregate helpers the figure script shares with the generator
# ---------------------------------------------------------------------------
def row_brief(r):
    """The columns of one chip that enter every figure and every pin."""
    return {
        "id": r.get("id"), "day": day_of(r), "state": r.get("state"),
        "segment": r.get("segment"), "track": r.get("segment"),
        "role": r.get("role"), "platform": r.get("platform"),
        "state_window": r.get("state_window"),
        "gamma2": _num(r.get("gamma2_mask")), "P": _num(r.get("P")),
        "D": _num(r.get("D")), "A": _num(r.get("A")), "F": _num(r.get("F")),
        "S": _num(r.get("S")),
        "temperature_c": _num(r.get("temperature_c")),
        "precipitation_mm": _num(r.get("precipitation_mm")),
        "shift_x": _num(r.get("shift_x")), "shift_y": _num(r.get("shift_y")),
        "coherence_gamma2": _num(r.get("coherence_gamma2")),
        "intensity": _num(r.get("intensity")),
        "brightness_ratio": _num(r.get("brightness_ratio")),
        "phase_coherence": _num(r.get("phase_coherence")),
        "phase_scatterer_count": _num(r.get("phase_scatterer_count")),
        "coherence_masked_pixels": _num(r.get("coherence_masked_pixels")),
    }


def summ(rows, key):
    """``st._stats`` over one channel — the CSV's own summary shape."""
    return st._stats([r.get(key) for r in rows])


def by_girder(rows):
    """Rows per track label, in ``GIRDERS`` order (unknown labels appended)."""
    keys = list(GIRDERS)
    for r in rows:
        s = r.get("segment")
        if s and s not in keys:
            keys.append(s)
    return {g: [r for r in rows if r.get("segment") == g] for g in keys}


# Alias the track wording to the girder helper name.
by_track = by_girder


def feature_column(k):
    """The CSV column that carries feature ``k`` (``gamma2`` is stored as
    ``gamma2_mask``; the pipeline's own ``gamma2_export`` lives beside it)."""
    return "gamma2_mask" if k == "gamma2" else k


def summary_by_state(rows):
    """``{state_short: {feature: st._stats}}`` — the metadata's summary block."""
    return {STATE_SHORT[s]: {k: summ([r for r in rows if r.get("state") == s],
                                     feature_column(k))
                             for k in FEATURES}
            for s in STATES}
