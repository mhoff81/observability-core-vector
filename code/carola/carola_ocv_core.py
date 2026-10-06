#!/usr/bin/env python3
"""Carola OCV-Paper — data layer, constants and the committed reference join.

Everything the other modules need, in one place: paths, the analysis channels,
the mask vector, the site constants, the collapse event and the loaders of the
committed tables.

The site
--------
Carola Bruecke (Dresden, 51.0545 N 13.7466 E) is a three-girder composite
bridge. The export stores one 80x80 complex chip per acquisition per **girder
request** (segments 0/1/2 -> Girder A/B/C, see ``data/carola/carola_segments.txt``),
plus 7x7 chips for a wrongly typed tower asset that this package excludes. The
collapse of 2024-09-11 divides the record into the two states of the site

    pre-collapse (healthy)   date <  2024-09-11
    post-collapse            date >= 2024-09-11

(``analyze_carola_coherence.STATE_LABELS`` / ``COLLAPSE_DATE``). Both states live
on the *same* bridge; the record straddles the event, so season, weather,
traffic and pipeline generation all change with the state — which is why every
contrast below carries a co-variate control beside it.

The two mask layers
-------------------
  * ``echo mask`` — the whole chip (the bridge deck as the window sees it), the
    layer of the 6D vector ``x = [gamma2, P, D, A, F, S]``;
  * ``per-girder deck mask`` — the same six dimensions, but per girder
    (``GIRDERS``), computed on that girder's own chip.

Both come from ``carola_ocv_masks`` and therefore from the committed window
cache; the whole-bridge values are pinned against
``carola_coherence_states.json`` and ``carola_echo_mask_gamma2.json``.

Reference files (``data/carola/reference/``)
--------------------------------------------
  carola_coherence_states.json    per-window gamma2 / n_masked of the first
                                  committed analysis (1,454 stored windows)
  carola_echo_mask_gamma2.json    the newer recomputation: state table, segment
                                  table, nested OLS and SDS of the mask vectors
  carola_bridge_osm.json          the committed bridge geometry (pinned
                                  structurally: it is a map extract, not
                                  rederivable from the record)
"""
from __future__ import annotations

import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import carola_ocv_masks as masks  # noqa: E402
import carola_ocv_stats as st     # noqa: E402

DATA = os.path.join(HERE, os.pardir, os.pardir, "data", "carola")
FIGDIR = os.path.join(HERE, os.pardir, os.pardir, "figures", "carola")
REFERENCE = os.path.join(DATA, "reference")

CACHE_PATH = masks.CACHE_PATH
MANIFEST_PATH = masks.MANIFEST_PATH
MEAS_PATH = masks.MEAS_PATH
INDEX_PATH = masks.INDEX_PATH
SEGMENTS_PATH = masks.SEGMENTS_PATH

CSV_PATH = os.path.join(DATA, "carola_ocv_channels.csv")
CSV_META_PATH = os.path.join(DATA, "carola_ocv_channels_meta.json")
REF_STATES = os.path.join(REFERENCE, "carola_coherence_states.json")
REF_MASK = os.path.join(REFERENCE, "carola_echo_mask_gamma2.json")
REF_OSM = os.path.join(REFERENCE, "carola_bridge_osm.json")

SITE = "Carola Bruecke, Dresden"
SITE_SHORT = "Carola"
COLLAPSE_DATE = "2024-09-11"
EVENT = COLLAPSE_DATE

STATES = list(st.STATE_ORDER)
PRE, POST = STATES
STATE_SHORT = dict(st.STATE_SHORT)

OCV = ["gamma2", "P", "D"]
FEATURES = ["gamma2", "P", "D", "A", "F", "S"]
DIMS = list(FEATURES)
CHANNEL_AXIS = {
    "gamma2": "gamma^2  (masked coherence)",
    "P": "P  (mask persistence)",
    "D": "D  (mask density A / bbox)",
    "A": "A  (echo pixels n_masked)",
    "F": "F  (8-connected components)",
    "S": "S  (||centroid - peak||, px)",
    "wind_speed_ms": "wind (m/s)",
    "temperature_c": "temperature (C)",
    "brightness_ratio": "brightness ratio",
    "measured_frequency_hz": "freq (Hz)",
    "coherence": "coherence",
}
CHANNEL_TITLE = {
    "gamma2": "gamma2", "P": "P", "D": "D", "A": "A", "F": "F", "S": "S",
    "wind_speed_ms": "wind (m/s)", "temperature_c": "temperature (C)",
    "brightness_ratio": "brightness ratio", "measured_frequency_hz": "freq (Hz)",
    "coherence": "coherence", "precipitation_mm": "precip (mm)",
    "wind_gust_ms": "gust (m/s)",
}

WEATHER = ["wind_speed_ms", "temperature_c", "precipitation_mm", "wind_gust_ms"]
CONTROLS = ["wind_speed_ms", "temperature_c", "brightness_ratio",
            "measured_frequency_hz", "coherence"]

GIRDERS = ["Girder A", "Girder B", "Girder C"]
BRIDGE_ASSETS = ["Carola R", "Carolabr", "Carola RPT A-full"]

MODELS = {
    "gamma2": ["gamma2"],
    "PD": ["P", "D"],
    "gamma2PD": ["gamma2", "P", "D"],
    "x_6d": ["gamma2", "P", "D", "A", "F", "S"],
    "weather": ["wind_speed_ms", "temperature_c", "precipitation_mm", "wind_gust_ms"],
}
HEADLINE_PAIR = "x_6d"
LAMBDA_GRID = [0.0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9]
HEADLINE_LAMBDA = 0.5
RNG_SEED = 7
N_BOOT = 10000
N_PERM = 1000

MASK_COLUMNS = ["A", "D", "F", "S", "P", "gamma2_mask", "gamma2_export",
                "coherence_masked_pixels", "mask_bbox_area", "mask_row_span",
                "mask_col_span", "mask_centroid_row", "mask_centroid_col",
                "mask_peak_row", "mask_peak_col", "mask_n_components",
                "mask_largest_component", "mask_fragmentation"]

SOURCE_NUMERIC = ["wind_speed_ms", "temperature_c", "precipitation_mm",
                  "wind_gust_ms", "relative_humidity_pct", "incidence_angle_deg",
                  "intensity", "brightness_ratio", "measured_frequency_hz",
                  "baseline_frequency_hz", "frequency_drop_pct",
                  "displacement_los_m", "coherence", "coherence_masked_pixels",
                  "coherence_gamma2", "phase_coherence", "phase_snr_db",
                  "coherent_sum_amplitude"]
SOURCE_TEXT = ["asset_id", "asset_name", "request_id", "segment_index",
               "acquisition_ts", "pass_label", "orbit_direction",
               "condition_label", "campaign_label", "traffic_load_label",
               "status"]
SOURCE_COLUMNS = ["id"] + SOURCE_TEXT + SOURCE_NUMERIC
DERIVED_COLUMNS = ["segment", "state", "pre_collapse", "day", "month", "year",
                   "season", "echo_mode", "girder_index"]
CSV_COLUMNS = SOURCE_COLUMNS + DERIVED_COLUMNS + MASK_COLUMNS


def _num(v):
    """CSV/JSON cell -> float or None."""
    if v is None or v == "" or v == "null" or v == "None":
        return None
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


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


# ---------------------------------------------------------------------------
# Loaders of the committed tables
# ---------------------------------------------------------------------------
def load_csv(path=CSV_PATH):
    """The generated channel table, numeric cells as floats."""
    with open(path) as fh:
        rows = list(csv.DictReader(fh))
    for r in rows:
        for c in SOURCE_NUMERIC:
            r[c] = _num(r.get(c))
        for c in ("A", "D", "F", "S", "P", "gamma2_mask", "gamma2_export",
                  "coherence_masked_pixels", "mask_bbox_area", "mask_row_span",
                  "mask_col_span", "mask_centroid_row", "mask_centroid_col",
                  "mask_n_components", "mask_largest_component",
                  "mask_fragmentation"):
            r[c] = _num(r.get(c))
        r["gamma2"] = r["gamma2_mask"]
        r["girder_index"] = _num(r.get("girder_index"))
        r["pre_collapse"] = str(r.get("pre_collapse")) in ("1", "True", "true")
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
    """The committed v2 measurement extract as ``id`` -> column dict."""
    out = {}
    with open(path) as fh:
        lines = fh.read().splitlines()
    if not lines:
        return out
    keys = lines[0].split("|")
    for line in lines[1:]:
        p = line.split("|")
        if len(p) != len(keys):
            continue
        d = dict(zip(keys, p))
        if d.get("id"):
            out[d["id"]] = d
    return out


def load_segments(path=SEGMENTS_PATH):
    """``(asset_id, segment_index) -> {label, fundamental_hz}``."""
    out = {}
    if not os.path.isfile(path):
        return out
    with open(path) as fh:
        lines = fh.read().splitlines()
    keys = lines[0].split("|")
    for line in lines[1:]:
        p = line.split("|")
        if len(p) < 5:
            continue
        d = dict(zip(keys, p))
        try:
            sidx = int(p[2])
        except ValueError:
            continue
        out[(p[0], sidx)] = {"label": p[3],
                             "fundamental_hz": _num(d.get("fundamental_hz")),
                             "asset_name": p[1]}
    return out


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
        "segment": r.get("segment"), "pass": r.get("pass_label"),
        "orbit": r.get("orbit_direction"), "request_id": r.get("request_id"),
        "gamma2": _num(r.get("gamma2_mask")), "P": _num(r.get("P")),
        "D": _num(r.get("D")), "A": _num(r.get("A")), "F": _num(r.get("F")),
        "S": _num(r.get("S")),
        "wind_speed_ms": _num(r.get("wind_speed_ms")),
        "temperature_c": _num(r.get("temperature_c")),
        "brightness_ratio": _num(r.get("brightness_ratio")),
        "measured_frequency_hz": _num(r.get("measured_frequency_hz")),
        "coherence": _num(r.get("coherence")),
    }


def summ(rows, key):
    """``st._stats`` over one channel — the CSV's own summary shape."""
    return st._stats([r.get(key) for r in rows])


def by_girder(rows):
    """Rows per girder label, in ``GIRDERS`` order (unknown labels appended)."""
    keys = list(GIRDERS)
    for r in rows:
        s = r.get("segment")
        if s and s not in keys:
            keys.append(s)
    return {g: [r for r in rows if r.get("segment") == g] for g in keys}


def feature_column(k):
    """The CSV column that carries feature ``k`` (``gamma2`` is stored as
    ``gamma2_mask``; the pipeline's own ``coherence_gamma2`` lives beside it)."""
    return "gamma2_mask" if k == "gamma2" else k


def summary_by_state(rows):
    """``{state_short: {feature: st._stats}}`` — the metadata's summary block."""
    return {STATE_SHORT[s]: {k: summ([r for r in rows if r.get("state") == s],
                                     feature_column(k))
                             for k in FEATURES}
            for s in STATES}
