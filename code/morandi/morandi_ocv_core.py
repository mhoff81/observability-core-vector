#!/usr/bin/env python3
"""Morandi OCV-Paper — data layer, constants and the committed reference join.

Everything the other modules need, in one place: paths, the analysis channels,
the mask vector, the site constants, the collapse event and the loaders of the
committed tables.

The site
--------
Ponte Morandi / Viadotto Polcevera (Genova, 44.4243 N 8.8906 E) is a three-span
balanced-cantilever viaduct. The export stores one 80x80 complex chip per
acquisition per **track** (``A_asc`` ASC rel-15 IW1 VV, ``A_des`` DESC rel-66
IW1 VV, see ``data/morandi/morandi_tracks.txt``) — the two Sentinel-1
observation geometries of the deck, the girder-equivalents of the Carola
package. The collapse of 2018-08-14 divides the record into the two states

    pre-collapse (healthy)   date <  2018-08-14
    post-collapse            date >= 2018-08-14

The record straddles the event, so season, weather and pipeline generation all
change with the state — which is why every contrast below carries a co-variate
control beside it. Only 11 chips survive after the event (see
``code/morandi/README.md``): the contrast is reported as an **honest null**, not
as damage.

The two mask layers
-------------------
  * ``echo mask`` — the whole chip (the deck as the window sees it), the layer
    of the 6D vector ``x = [gamma2, P, D, A, F, S]``;
  * ``per-track deck mask`` — the same six dimensions, but per track
    (``GIRDERS``), computed on that track's own chip.

Both come from ``morandi_ocv_masks`` and therefore from the committed window
cache; the per-date measurements are pinned against
``data/morandi/reference/morandi_registration.json``.

Reference files (``data/morandi/reference/``)
--------------------------------------------
  morandi_registration.json       the committed co-registration/QC table: per
                                  track n_dates/n_pre/n_post, first/last/master
                                  day, the fixed 640-px deck mask size and the
                                  per-date whole-window metrics
  morandi_deck_channels.json      the site's fixed quantile deck-mask table
                                  (640 px per track), pinned structurally
  morandi_geometry.json           the viaduct geometry (assumption grid, see
                                  morandi_geometry.md; not rederivable from the
                                  record)
"""
from __future__ import annotations

import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import morandi_ocv_masks as masks  # noqa: E402
import morandi_ocv_stats as st     # noqa: E402

DATA = os.path.join(HERE, os.pardir, os.pardir, "data", "morandi")
FIGDIR = os.path.join(HERE, os.pardir, os.pardir, "figures", "morandi")
REFERENCE = os.path.join(DATA, "reference")

CACHE_PATH = masks.CACHE_PATH
MANIFEST_PATH = masks.MANIFEST_PATH
MEAS_PATH = masks.MEAS_PATH
TRACKS_PATH = masks.TRACKS_PATH
SEGMENTS_PATH = TRACKS_PATH            # alias: a "segment" is a track here

CSV_PATH = os.path.join(DATA, "morandi_ocv_channels.csv")
CSV_META_PATH = os.path.join(DATA, "morandi_ocv_channels_meta.json")
DECK_CSV_PATH = os.path.join(DATA, "morandi_deck_channels.csv")
REF_STATES = os.path.join(REFERENCE, "morandi_registration.json")
REF_MASK = os.path.join(REFERENCE, "morandi_deck_channels.json")
REF_OSM = os.path.join(REFERENCE, "morandi_geometry.json")

SITE = "Ponte Morandi (Polcevera), Genova"
SITE_SHORT = "Morandi"
COLLAPSE_DATE = "2018-08-14"
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
    "precipitation_mm": "precip (mm)",
    "wind_gust_ms": "gust (m/s)",
    "nearest_lag_h": "weather lag (h)",
}
CHANNEL_TITLE = {
    "gamma2": "gamma2", "P": "P", "D": "D", "A": "A", "F": "F", "S": "S",
    "wind_speed_ms": "wind (m/s)", "temperature_c": "temperature (C)",
    "precipitation_mm": "precip (mm)", "wind_gust_ms": "gust (m/s)",
    "nearest_lag_h": "weather lag (h)",
}

WEATHER = ["wind_speed_ms", "temperature_c", "precipitation_mm", "wind_gust_ms"]
CONTROLS = ["wind_speed_ms", "temperature_c", "precipitation_mm", "wind_gust_ms",
            "nearest_lag_h", "temperature_prev6h_c"]

# The tracks play the role of Carola's girders; the name ``GIRDERS`` (and the
# helper ``by_girder``) is kept so the figure scripts read identically.
GIRDERS = ["A_asc", "A_des"]

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

MASK_COLUMNS = ["A", "D", "F", "S", "P", "gamma2_mask", "coherence_masked_pixels",
                "mask_bbox_area", "mask_row_span", "mask_col_span",
                "mask_centroid_row", "mask_centroid_col", "mask_peak_row",
                "mask_peak_col", "mask_n_components", "mask_largest_component",
                "mask_fragmentation"]

SOURCE_NUMERIC = ["nearest_lag_h", "temperature_c", "temperature_prev6h_c",
                  "dT_6h_c", "precipitation_mm", "rain_sum_24h_mm",
                  "rain_sum_72h_mm", "rain_sum_7d_mm", "wind_speed_ms",
                  "wind_gust_ms", "wind_mean_6h_ms", "relative_humidity_pct",
                  "shortwave_wm2", "snow_depth_m",
                  "intensity", "gamma2_export", "phase_R", "phase_R_null",
                  "coh_masked", "coh_vs_master", "shift_abs_px",
                  "coherence_masked_pixels",
                  "deck_intensity", "deck_gamma2_export", "deck_phase_R",
                  "dphi_deck_minus_ref"]
SOURCE_TEXT = ["track", "orbit_direction", "platform", "group",
               "acquisition_ts", "status"]
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
    """0/1 index of the track (the girder-request index analogue)."""
    try:
        return GIRDERS.index(track)
    except ValueError:
        return None


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
        r["gamma2_export"] = r.get("gamma2_export")
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
    """The committed weather/QC measurement extract, keyed ``(track, date)``.

    The weather table carries the whole-window per-date metrics (``intensity``,
    ``gamma2``, ``phase_R``, ``coh_vs_master``, ``coh_masked``, ``shift_abs_px``)
    beside the Sentinel-1 relative-orbit weather columns. ``gamma2`` is
    re-exposed as ``gamma2_export`` (the pipeline's *own* spatial coherence,
    kept for provenance only — the mask layer recomputes its own ``gamma2``),
    ``t_utc`` as ``acquisition_ts`` and ``orbit`` as ``orbit_direction``.
    """
    out = {}
    try:
        reg = {t["tag"]: t for t in load_reference(REF_STATES).get("tracks", [])}
    except (OSError, KeyError, ValueError):
        reg = {}
    for d in _read_pipe(path):
        if not d.get("track"):
            continue
        d["gamma2_export"] = d.get("gamma2")
        d["acquisition_ts"] = d.get("t_utc")
        d["orbit_direction"] = d.get("orbit")
        d["status"] = ""
        # the site's fixed quantile deck mask size (640 px per track) — export
        # provenance only; the package recomputes the echo mask.
        t = reg.get(d["track"])
        if t is not None:
            d["coherence_masked_pixels"] = t.get("n_masked_pixels")
        out[(d["track"], d["date"])] = d
    return out


def load_deck_channels(path=DECK_CSV_PATH):
    """The site's fixed quantile deck-mask table, keyed ``(track, date)``.

    Carries ``deck_intensity``, ``deck_gamma2`` (-> ``deck_gamma2_export``),
    ``deck_phase_R`` and ``dphi_deck_minus_ref``. This is the *exported* deck
    channel, computed by the site with a fixed 640-px (10 %) quantile mask; the
    package's own per-track deck mask is recomputed from the echo rule and does
    **not** use this table (see ``code/morandi/README.md``).
    """
    out = {}
    for d in _read_pipe(path):
        if not d.get("track"):
            continue
        d["deck_gamma2_export"] = d.get("deck_gamma2")
        out[(d["track"], d["date"])] = d
    return out


def load_segments(path=SEGMENTS_PATH):
    """``track -> {label, orbit, ...}`` (the track/segment table)."""
    return masks.load_tracks(path)


def load_tracks(path=TRACKS_PATH):
    """``track -> {label, orbit, ...}``."""
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
        "orbit": r.get("orbit_direction"), "platform": r.get("platform"),
        "gamma2": _num(r.get("gamma2_mask")), "P": _num(r.get("P")),
        "D": _num(r.get("D")), "A": _num(r.get("A")), "F": _num(r.get("F")),
        "S": _num(r.get("S")),
        "wind_speed_ms": _num(r.get("wind_speed_ms")),
        "temperature_c": _num(r.get("temperature_c")),
        "precipitation_mm": _num(r.get("precipitation_mm")),
        "wind_gust_ms": _num(r.get("wind_gust_ms")),
        "coherence": _num(r.get("coh_masked")),
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

