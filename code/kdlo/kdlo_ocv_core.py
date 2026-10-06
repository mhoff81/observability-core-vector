#!/usr/bin/env python3
"""KDLO OCV-Paper — data layer and constants.

Reads only committed files:

  * ``data/kdlo/kdlo_ocv_channels.csv``      (generated, see
    ``code/kdlo/kdlo_ocv_channels_csv.py``)
  * ``data/kdlo/kdlo_ocv_channels_meta.json``(provenance of the generator)
  * the four committed reference JSONs in ``data/kdlo/reference/``

No live API, no Postgres, no burst cache, no third-party imports beyond numpy.
Field names, orderings and aggregate definitions are the ones of the KDLO
project scripts (``kdlo_tower_analysis/fig_kdlo_channels.py``): ``summ`` is
verbatim, the echo-mode split and the collapse pair are verbatim.

What is different from the LUMO/Bautzen packages, and why it is written down
instead of papered over:

  * **state = epoch.** LUMO/Bautzen states are damage labels; a KDLO state is a
    *time window* (``pre`` = standing mast 2022-06-02..2022-12-11, ``rebuild`` =
    reconstruction 2023-05-04..2024-09-19). One acquisition (2022-12-23) lies
    between the two windows and belongs to neither: it is labelled ``single``
    and excluded from every state statistic.
  * **the mask layer was empty, and is not any more.** The pipeline stored no
    echo-mask geometry — only the size ``A = coherence_masked_pixels``. The mask
    of every overpass is now recomputed from the committed full-dwell strip
    cache (``data/kdlo/strips/``, see ``kdlo_ocv_masks``) with the pipeline's own
    rule, so ``P``, ``D``, ``F`` and ``S`` exist and the OCV is
    ``[gamma2, P, D]`` over the 6D echo-mask vector
    ``[gamma2, P, D, A, F, S]`` — the shape of the LUMO package, with ``A``
    pinned bit-exactly to the pipeline's own count and ``gamma2`` to its own
    coherence.
  * **one geometry.** Only S1A descending relative orbit 12 at ~12:43 UTC exists
    over this site before 2025-04; there is no ascending/descending split.
  * **"after the rebuild" does not exist.** The record ends 2024-09-19, i.e.
    *during* the construction (1,542 ft reached 2024-08-19). The second state is
    therefore "during the rebuild", not "repaired".
"""
from __future__ import annotations

import csv
import datetime as dt
import hashlib
import json
import os

import numpy as np

import kdlo_ocv_stats as st

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, os.pardir, os.pardir, "data", "kdlo")
REF = os.path.join(DATA, "reference")
FIGDIR = os.path.join(HERE, os.pardir, os.pardir, "figures", "kdlo")

CSV_PATH = os.path.join(DATA, "kdlo_ocv_channels.csv")
CSV_META_PATH = os.path.join(DATA, "kdlo_ocv_channels_meta.json")
DB_EXTRA_PATH = os.path.join(DATA, "kdlo_db_extra.json")
RAW_PATHS = {1: os.path.join(DATA, "fig_kdlo_1_raw.json"),
             2: os.path.join(DATA, "fig_kdlo_2_raw.json")}

# ---------------------------------------------------------------------------
# Pipeline constants — must stay identical to the KDLO backend / project scripts
# ---------------------------------------------------------------------------
PRED_MIN_COHERENCE = 0.15          # pipeline coherence floor
COMPACT_MASK_PX = 10               # <= 10 masked px  -> compact point echo
MASK_CLUTTER_PX = 25               # >= 25 masked px  -> distributed clutter
# LUMO reference medians (documentation only, never recomputed here)
LUMO_ASC_G2_MEDIAN = 0.0617
LUMO_DESC_G2_MEDIAN = 0.0261
LUMO_HEALTHY_G2_MEDIAN = 0.0379

# The two figure windows of the project record (verbatim from PARTS)
PARTS = {
    1: {"key": "part1", "window": ("2022-06-01", "2022-12-23"),
        "nominal_pass": "afternoon", "collapse_start": "2022-12-11",
        "reference_from": None},
    2: {"key": "part2", "window": ("2023-05-01", "2024-09-30"),
        "nominal_pass": "afternoon", "collapse_start": None,
        "reference_from": 1},
}
PART_ORDER = [1, 2]

# ---------------------------------------------------------------------------
# States = epochs. The windows are the committed sampling dates of the two
# parts minus the one acquisition that belongs to neither.
# ---------------------------------------------------------------------------
SINGLETON_DAY = "2022-12-23"       # the only acquisition between the windows
STATE_WINDOWS = {"pre": ("2022-06-02", "2022-12-11"),
                 "rebuild": ("2023-05-04", "2024-09-19")}
STATE_LABEL = {"pre": "pre-collapse (standing mast, 2022-06-02..2022-12-11)",
               "rebuild": "during the rebuild (2023-05-04..2024-09-19)"}
STATE_LIST = list(st.STATE_ORDER)

# ---------------------------------------------------------------------------
# Channels
# ---------------------------------------------------------------------------
# The OCV proper (what is asked to carry the separation) ...
OCV = ["gamma2", "P", "D"]
# ... and the controls (everything else the record carries).
CONTROLS = ["phase_coherence", "phase_snr_db", "phase_rms_rad",
            "phase_detected_frequency_hz", "frequency_drop_pct",
            "displacement_los_m", "sub_aperture_modulation",
            "coherent_sum_amplitude", "intensity", "incidence_angle_deg",
            "wind_speed_ms", "wind_gust_ms"]
WEATHER = ["temperature_c", "wind_speed_ms", "snow_depth_m",
           "precipitation_mm"]
# Channels the record carries but that are empty for this asset (0 of 60
# measurements) — kept in the CSV so that "empty" is a reproducible fact.
EMPTY_CHANNELS = ["ndvi", "measured_frequency_hz"]
FIELDS = ("coherence_gamma2", "coherence_masked_pixels", "phase_coherence",
          "phase_snr_db", "phase_rms_rad")
CHANNEL_AXIS = {"gamma2": "gamma^2  [1]", "A": "A  [masked px]",
                "P": "P  [1]", "D": "D  [1]", "F": "F  [components]",
                "S": "S  [px]",
                "phase_coherence": "phase coherence  [1]",
                "phase_snr_db": "phase SNR  [dB]",
                "phase_rms_rad": "phase RMS  [rad]",
                "intensity": "intensity  [1]",
                "incidence_angle_deg": "incidence angle  [deg]",
                "wind_speed_ms": "wind speed  [m/s]",
                "wind_gust_ms": "wind gust  [m/s]",
                "temperature_c": "temperature  [deg C]",
                "phase_detected_frequency_hz": "phase freq  [Hz]"}
CHANNEL_TITLE = {"gamma2": "interferometric coherence gamma^2",
                 "A": "echo-mask size (clutter discriminator)",
                 "P": "echo-mask persistence",
                 "D": "echo-mask density",
                 "F": "echo-mask fragmentation",
                 "S": "echo-mask centroid<->peak shift",
                 "phase_coherence": "tower phase coherence",
                 "phase_snr_db": "tower phase SNR",
                 "phase_rms_rad": "tower phase RMS activity",
                 "intensity": "backscatter intensity",
                 "incidence_angle_deg": "incidence angle",
                 "wind_speed_ms": "wind speed at acquisition",
                 "wind_gust_ms": "wind gust at acquisition",
                 "temperature_c": "air temperature at acquisition",
                 "phase_detected_frequency_hz": "phase-detected frequency"}

# Columns only Postgres carries (read once, vendored into
# ``data/kdlo/kdlo_db_extra.json`` by ``kdlo_ocv_channels_csv.py --fetch-db-extra``)
DB_EXTRA_COLUMNS = (
    "phase_snr_db", "phase_observable", "phase_unobservable_reason",
    "snow_depth_m", "snowfall_cm", "precipitation_mm", "wind_gust_ms",
    "intensity", "incidence_angle_deg",
)

# Source columns: the exact key order of the cached API rows
# (``fig_kdlo_{part}_raw.json`` -> rows[*]).
SOURCE_COLUMNS = [
    "id", "acquisition_ts", "pass_label", "orbit_direction",
    "displacement_los_m", "brightness_ratio", "measured_frequency_hz",
    "baseline_frequency_hz", "frequency_drop_pct", "wind_speed_ms",
    "temperature_c", "traffic_load_label", "status",
    "phase_detected_frequency_hz", "phase_coherence", "phase_status",
    "coherence_gamma2", "coherence_masked_pixels", "ndvi", "ndvi_scene_date",
    "ndvi_confound_flagged", "sub_aperture_modulation",
    "coherent_sum_amplitude", "registration_status", "phase_rms_rad",
    "condition_label", "campaign_label",
]
# Derived columns the generator prepends (all recomputable from the source row)
DERIVED_COLUMNS = ["part", "state", "month", "year", "season", "echo_mode"]
# Mask-geometry columns the generator appends. They are recomputed from the
# committed full-dwell strip cache (``data/kdlo/strips/``) by
# ``code/kdlo/kdlo_ocv_masks.py`` with the pipeline's own mask rule, so the echo
# mask stops being a number the record merely asserts: ``A`` must reproduce
# ``coherence_masked_pixels`` and ``gamma2_mask`` must reproduce
# ``coherence_gamma2`` (both checked by ``kdlo_ocv_channels_csv.py --verify``).
MASK_COLUMNS = ["A", "D", "F", "S", "P",
                "mask_rows", "mask_cols", "bbox_area", "row_span", "col_span",
                "centroid_r", "centroid_c", "peak_r", "peak_c",
                "peak_intensity", "gamma2_mask"]
# The generated CSV = derived + source + the vendored DB columns + the masks
CSV_COLUMNS = (DERIVED_COLUMNS + SOURCE_COLUMNS + list(DB_EXTRA_COLUMNS)
               + MASK_COLUMNS)


REFERENCE_FILES = {
    "kdlo_1_channels": "fig_kdlo_1_channels.json",
    "kdlo_2_channels": "fig_kdlo_2_channels.json",
    "s1_availability": "kdlo_s1_availability.json",
    "fem_frequencies": "kdlo_fem_frequencies.json",
}
# The two reference files that are *not* recomputed from the record: they are
# CDSE-catalogue and FEM-solver output. They are pinned structurally instead —
# see ``availability_checks`` and ``fem_checks`` below.
PIN_ONLY_FILES = ("s1_availability", "fem_frequencies")

# ---------------------------------------------------------------------------
# Models / statistics constants
# ---------------------------------------------------------------------------
# The 6D echo-mask vector, in the feature order of the LUMO package so that the
# two packages stay comparable: gamma2 first, then P and D (the OCV), then the
# remaining mask dimensions A, F, S.
FEATURES = ["gamma2", "P", "D", "A", "F", "S"]
DIMS = ["A", "D", "F", "S", "P"]
DIM_LABEL = {
    "A": "Area (n_masked, number of echo pixels)",
    "D": "Density (n_masked / bbox_area)",
    "F": "Fragmentation (number of connected components, 8-neighbourhood)",
    "S": "Shift (distance mask centroid <-> peak pixel, px)",
    "P": "Persistence (mean pixel frequency across all dates inside the own mask)",
}
MODELS = {"gamma2": ["gamma2"], "PDS": ["P", "D", "S"],
          "gamma2PD": ["gamma2", "P", "D"], "x_6d": FEATURES}
MODEL_LABEL = {
    "gamma2": "gamma2 (1D)",
    "PDS": "[P, D, S] (mask-morphology control)",
    "gamma2PD": "[gamma2, P, D] (headline, 3 dims)",
    "x_6d": "[gamma2, P, D, A, F, S] (full 6D vector)",
}
HEADLINE_PAIR = "gamma2PD"
BASELINE = "gamma2"
SECONDARY_PAIR = "PDS"
LAMBDA_GRID = [0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0]
HEADLINE_LAMBDA = 0.5
RNG_SEED = 7
N_BOOT = 10000
N_BOOT_STRATA = 2000
N_PERM_2CLASS = 1000
N_PERM_PAIR = 1000
# n of the two epoch states in the committed record (16 + 43 = 59 of 60 rows;
# the 60th is the singleton 2022-12-23, which belongs to neither state)
N_BY_STATE_EXPECTED = {"pre": 16, "rebuild": 43}


# ---------------------------------------------------------------------------
# Primitives
# ---------------------------------------------------------------------------
def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _num(v):
    """Float or None — verbatim ``fig_kdlo_channels._num`` (bools are not numbers)."""
    if v is None or isinstance(v, bool):
        return None
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def _int(v):
    x = _num(v)
    return None if x is None else int(x)


def _text(v):
    if v is None:
        return None
    s = str(v).strip()
    return s or None


def day_of(row):
    """``YYYY-MM-DD`` of an acquisition — the key of the whole package."""
    ts = row.get("acquisition_ts")
    return None if ts is None else str(ts)[:10]


def state_of(day):
    """Epoch of one acquisition day: ``pre`` / ``rebuild`` / ``single`` / None.

    ``single`` is the one acquisition between the two windows (2022-12-23): it is
    neither pre-collapse nor rebuild and is therefore never part of a state
    statistic. ``None`` means the day lies outside the record.
    """
    for lab in STATE_LIST:
        lo, hi = STATE_WINDOWS[lab]
        if lo <= day <= hi:
            return lab
    return "single" if day == SINGLETON_DAY else None


def echo_mode_of(n_px):
    """Pipeline echo-mode discriminator: compact / intermediate / distributed."""
    if n_px is None:
        return None
    if n_px <= COMPACT_MASK_PX:
        return "compact"
    if n_px >= MASK_CLUTTER_PX:
        return "distributed"
    return "intermediate"


# ---------------------------------------------------------------------------
# Loaders
# ---------------------------------------------------------------------------
def load_raw_parts(paths=None):
    """The two committed raw part files: ``({1: [...], 2: [...]}, provenance)``."""
    paths = paths or RAW_PATHS
    by_part, prov = {}, {}
    for p in PART_ORDER:
        with open(paths[p]) as fh:
            doc = json.load(fh)
        by_part[p] = doc.get("rows", [])
        prov[p] = {
            "file": os.path.relpath(paths[p], os.path.dirname(DATA)),
            "sha256": sha256(paths[p]),
            "asset_id": doc.get("asset_id"),
            "window": doc.get("window"),
            "request_ids": doc.get("request_ids"),
            "request_ids_used": doc.get("request_ids_used"),
            "n_rows": len(doc.get("rows", [])),
        }
    return by_part, prov


def load_db_extra(path=DB_EXTRA_PATH):
    """The vendored one-time Postgres extraction (read-only)."""
    with open(path) as fh:
        return json.load(fh)


def db_extra_rows(path=DB_EXTRA_PATH):
    doc = load_db_extra(path)
    return doc.get("rows", {}), {
        "file": os.path.relpath(path, os.path.dirname(DATA)),
        "sha256": sha256(path),
        "source": doc.get("source"),
        "columns": doc.get("columns"),
        "n_rows": len(doc.get("rows", {})),
        "fetched_at": doc.get("fetched_at"),
    }



def merge_rows(by_part=None, extra=None):
    """The 60 acquisitions, merged with the vendored DB columns, in part order.

    Row keys: the 27 source columns + the 9 DB columns + the derived
    ``part``/``state``/``month``/``year``/``season``/``echo_mode``.

    The DB columns are looked up by the measurement ``id`` (the Postgres primary
    key), not by date or request: the map is keyed by ``id`` so that a duplicated
    or shifted acquisition can never silently receive another row's SNR/weather.
    """
    if by_part is None:
        by_part, _ = load_raw_parts()
    if extra is None:
        extra, _ = db_extra_rows()
    out = []
    for p in PART_ORDER:
        for src in sorted(by_part[p], key=lambda r: str(r.get("acquisition_ts"))):
            row = dict(src)
            day = day_of(row)
            info = extra.get(str(row.get("id"))) or {}
            for col in DB_EXTRA_COLUMNS:
                row[col] = info.get(col)
            row["part"] = f"part{p}"
            row["state"] = state_of(day)
            row["month"] = int(day[5:7])
            row["year"] = int(day[:4])
            row["season"] = st.season_of(row)
            row["echo_mode"] = echo_mode_of(_num(row.get("coherence_masked_pixels")))
            out.append(row)
    return out


def _csv_cell(v):
    """One CSV cell. Floats use ``repr`` so that reading back is bit-exact."""
    if v is None:
        return ""
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, float):
        return repr(v)
    return str(v)


def write_csv(rows, path=CSV_PATH):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=CSV_COLUMNS, lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({c: _csv_cell(r.get(c)) for c in CSV_COLUMNS})
    return path


TEXT_SOURCE_COLUMNS = ("id", "acquisition_ts", "pass_label", "orbit_direction",
                       "status", "traffic_load_label", "phase_status",
                       "ndvi_scene_date", "registration_status",
                       "condition_label", "campaign_label")


def load_csv(path=CSV_PATH):
    """The 60 merged acquisitions in committed order (all states)."""
    rows = []
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh):
            row = dict(r)
            for col in SOURCE_COLUMNS:
                v = (r.get(col) or "").strip()
                row[col] = (v or None) if col in TEXT_SOURCE_COLUMNS else _num(v)
            for col in DB_EXTRA_COLUMNS:
                v = (r.get(col) or "").strip()
                row[col] = (v or None) if col in ("phase_observable",
                                                  "phase_unobservable_reason") \
                    else _num(v)
            row["part"] = r["part"]
            row["state"] = r.get("state") or None
            row["month"] = _int(r.get("month"))
            row["year"] = _int(r.get("year"))
            row["season"] = r.get("season") or None
            row["echo_mode"] = r.get("echo_mode") or None
            row["gamma2"] = _num(r.get("coherence_gamma2"))
            row["A"] = _num(r.get("coherence_masked_pixels"))
            # Echo-mask layer: recomputed from the committed strip cache by the
            # generator (A is the pipeline's own count, gamma2_mask its own
            # coherence — both verified there).
            for col in MASK_COLUMNS:
                row[col] = _num((r.get(col) or "").strip())
            for col in ("mask_rows", "mask_cols", "bbox_area", "row_span",
                        "col_span", "peak_r", "peak_c"):
                row[col] = _int(r.get(col))
            rows.append(row)
    return rows


def load_reference():
    out = {}
    for name, fname in REFERENCE_FILES.items():
        with open(os.path.join(REF, fname)) as fh:
            out[name] = json.load(fh)
    return out


def load_csv_meta(path=CSV_META_PATH):
    with open(path) as fh:
        return json.load(fh)



# ---------------------------------------------------------------------------
# Aggregates — verbatim from ``fig_kdlo_channels.py`` (summ/stats/_echo_modes/
# _row_brief), so that the committed channel JSONs can be reproduced exactly
# ---------------------------------------------------------------------------
def summ(rows, field):
    """Verbatim ``fig_kdlo_channels.summ``: n/median/p25/p75/min/max."""
    vals_ = [v for v in (_num(r.get(field)) for r in rows) if v is not None]
    return {"n": len(vals_),
            "median": float(np.median(vals_)) if vals_ else None,
            "p25": float(np.percentile(vals_, 25)) if vals_ else None,
            "p75": float(np.percentile(vals_, 75)) if vals_ else None,
            "min": float(np.min(vals_)) if vals_ else None,
            "max": float(np.max(vals_)) if vals_ else None}


def row_brief(r):
    """Verbatim ``fig_kdlo_channels._row_brief``."""
    return {
        "day": day_of(r),
        "pass_label": r.get("pass_label"),
        "gamma2": _num(r.get("coherence_gamma2")),
        "masked_pixels": _num(r.get("coherence_masked_pixels")),
        "temperature_c": _num(r.get("temperature_c")),
        "wind_speed_ms": _num(r.get("wind_speed_ms")),
        "phase_coherence": _num(r.get("phase_coherence")),
    }


def echo_modes_of(rows):
    """Verbatim ``fig_kdlo_channels._echo_modes`` (compact / distributed buckets)."""
    px = [(day_of(r), _num(r.get("coherence_masked_pixels")),
           _num(r.get("coherence_gamma2"))) for r in rows]

    def bucket(lo, hi):
        sel = [(d, g) for d, p, g in px if p is not None and lo <= p <= hi]
        gs = [g for _d_, g in sel if g is not None]
        return {"n": len(sel), "days": sorted(d for d, _ in sel),
                "gamma2_median": float(np.median(gs)) if gs else None,
                "mask_range_px": [lo, hi]}

    return {"compact": bucket(0, COMPACT_MASK_PX),
            "distributed": bucket(MASK_CLUTTER_PX, 10 ** 6)}


def phase_observable_of(rows):
    """Verbatim ``fig_kdlo_channels`` phase-observability counters."""
    return {
        "observable": sum(1 for r in rows
                          if str(r.get("phase_observable")) == "t"),
        "not_observable": sum(1 for r in rows
                              if str(r.get("phase_observable")) == "f"),
        "reasons": sorted({str(r.get("phase_unobservable_reason"))
                           for r in rows
                           if r.get("phase_unobservable_reason")}),
    }


def part_stats(rows, part_no, prov=None, reference=None):
    """Verbatim ``fig_kdlo_channels.stats`` for one part (rows already filtered).

    ``prov`` adds the request-id provenance of the part's raw file; ``reference``
    is the already-computed stats dict of the part this window refers to (part 2
    carries part 1 as its intact reference).
    """
    part = PARTS[part_no]
    out = {
        "part": part["key"],
        "window": {"start": part["window"][0], "end": part["window"][1]},
        "n_measurements": len(rows),
        "date_first": min((day_of(r) for r in rows), default=None),
        "date_last": max((day_of(r) for r in rows), default=None),
        "passes": sorted({r.get("pass_label") for r in rows
                          if r.get("pass_label")}),
        "alternate_pass_days": sorted(
            {day_of(r) for r in rows
             if r.get("pass_label")
             and r.get("pass_label") != part.get("nominal_pass", "afternoon")}),
        "orbits": sorted({r.get("orbit_direction") for r in rows
                          if r.get("orbit_direction")}),
        "lumo_reference": {"ascending_gamma2": LUMO_ASC_G2_MEDIAN,
                           "descending_gamma2": LUMO_DESC_G2_MEDIAN,
                           "healthy_gamma2": LUMO_HEALTHY_G2_MEDIAN},
        "pred_min_coherence": PRED_MIN_COHERENCE,
    }
    for field in FIELDS:
        out[field] = summ(rows, field)
    out["weather"] = {f: summ(rows, f) for f in
                      ("temperature_c", "wind_speed_ms", "snow_depth_m",
                       "precipitation_mm")}
    out["phase_observable"] = phase_observable_of(rows)
    out["echo_modes"] = echo_modes_of(rows)

    cut = part.get("collapse_start")
    if cut:
        pre = [r for r in rows if day_of(r) <= cut]
        post = [r for r in rows if day_of(r) > cut]
        if pre and post:
            a = max(pre, key=day_of)
            b = min(post, key=day_of)
            out["collapse_pair"] = {
                cut: {"last_intact": row_brief(a), "first_collapsed": row_brief(b)},
                "series_context": {
                    "gamma2_median": out["coherence_gamma2"]["median"],
                    "gamma2_min": out["coherence_gamma2"]["min"],
                    "gamma2_max": out["coherence_gamma2"]["max"],
                    "n": out["n_measurements"],
                },
            }
    if prov is not None:
        out["request_ids"] = prov[part_no]["request_ids"]
    if reference is not None:
        out["intact_reference"] = {
            "part": str(part["reference_from"]),
            "n": reference["n_measurements"],
            "window": reference["window"],
            "gamma2": reference["coherence_gamma2"],
        }
    return out


def rows_of_part(rows, part_no):
    return [r for r in rows if r.get("part") == PARTS[part_no]["key"]]


def series(rows, field):
    """``{orbit: [(day, value), ...]}`` sorted by day — verbatim ``series``."""
    out = {}
    for r in rows:
        v = _num(r.get(field))
        if v is None:
            continue
        orbit = r.get("orbit_direction") or "DESCENDING"
        out.setdefault(orbit, []).append((day_of(r), v))
    for k in out:
        out[k].sort()
    return out


def state_rows(rows):
    """Only the rows that belong to one of the two epoch states."""
    return [r for r in rows if r.get("state") in st.STATE_ORDER]


def n_by_state(rows):
    counts = {lab: 0 for lab in STATE_LIST}
    for r in rows:
        if r.get("state") in counts:
            counts[r["state"]] += 1
    return counts


def feature_dataset(rows, features=None):
    """X, y of the complete cases, in committed order."""
    features = features or FEATURES
    complete = [r for r in rows if all(r.get(f) is not None for f in features)]
    X = np.array([[float(r[f]) for f in features] for r in complete])
    y = np.array([r["state"] for r in complete], dtype=object)
    return X, y, complete


# ---------------------------------------------------------------------------
# Reference-file checks
#
# Two of the four committed reference files are *not* rederivable from the
# record: ``kdlo_s1_availability.json`` is a CDSE-catalogue query result and
# ``kdlo_fem_frequencies.json`` is an eigenmode-solver output. This package can
# therefore not recompute their leaves. What it *can* do — and does — is
#
#   * rederive every field that follows from the 60 acquisitions (the
#     acquisition-day grid, the burst/product accounting, the pipeline
#     request/day/anomaly counts, the collapse bracket ...), and
#   * check the arithmetic and cross-references inside the files, which is what
#     makes silent tampering impossible even for the leaves we must pin.
#
# Every entry below is a ``(path, got, want)`` triple; ``got`` is the value
# computed here, ``want`` the committed one, so a mismatch is always visible as
# "the record says X, the reference file says Y".
# ---------------------------------------------------------------------------
def days_of(rows):
    return sorted({d for d in (day_of(r) for r in rows) if d})


def add_days(day, n):
    return (dt.date.fromisoformat(day) + dt.timedelta(days=n)).isoformat()


def grid_missing_days(days, cycle=12):
    """Days a ``cycle``-day repeat grid expects between measurements but that
    the record does not contain (S1A repeat cycle over this site)."""
    have, out = set(days), []
    for a, b in zip(days, days[1:]):
        gap = (dt.date.fromisoformat(b) - dt.date.fromisoformat(a)).days
        for k in range(1, max(gap // cycle, 1)):
            d = add_days(a, k * cycle)
            if d not in have and d not in out:
                out.append(d)
    return out


def availability_checks(avail, rows, prov):
    """Recomputed + internal-consistency checks of the S1 availability file."""
    out = []

    def chk(path, got, want):
        out.append((path, got, want))

    site = avail["site"]
    chk("av.site.latitude_deg_in_range", -90.0 <= site["latitude"] <= 90.0, True)
    chk("av.site.longitude_deg_in_range", -180.0 <= site["longitude"] <= 180.0, True)
    chk("av.site.osm_node_positive_int",
        isinstance(site["osm_node"], int) and site["osm_node"] > 0, True)
    chk("av.site.osm_offset_m_non_negative", site["osm_offset_m"] >= 0, True)
    chk("av.site.locality_non_empty", bool(site["locality"].strip()), True)
    chk("av.repeat_cycle_days", avail["repeat_cycle_days"], 12)

    # --- the collapse bracket, rederived from the part-1 days -------------
    p1_days = days_of(rows_of_part(rows, 1))
    cut = PARTS[1]["collapse_start"]
    chk("av.collapse.last_before", max(d for d in p1_days if d <= cut),
        avail["collapse"]["last_before"])
    chk("av.collapse.first_after", min(d for d in p1_days if d > cut),
        avail["collapse"]["first_after"])
    chk("av.collapse.month", avail["collapse"]["last_before"][:7],
        avail["collapse"]["month"])
    chk("av.collapse.bracket_days",
        (dt.date.fromisoformat(avail["collapse"]["first_after"])
         - dt.date.fromisoformat(avail["collapse"]["last_before"])).days,
        avail["collapse"]["bracket_days"])
    chk("av.collapse.bracket_days_eq_repeat_cycle",
        avail["collapse"]["bracket_days"], avail["repeat_cycle_days"])
    chk("av.collapse.singleton_is_first_after",
        avail["collapse"]["first_after"], SINGLETON_DAY)

    # --- the same bracket, as the catalogue names it ----------------------
    br = avail["bracket"]
    lb, fa = br["last_before"][0], br["first_after"][0]
    chk("av.bracket.last_before[0].day", lb["day"], avail["collapse"]["last_before"])
    chk("av.bracket.first_after[0].day", fa["day"], avail["collapse"]["first_after"])
    chk("av.bracket.bracket_days", br["bracket_days"], avail["collapse"]["bracket_days"])
    chk("av.bracket.last_before[0].orbit_direction", lb["orbit_direction"],
        sorted({r.get("orbit_direction") for r in rows_of_part(rows, 1)})[0])
    chk("av.bracket.first_after[0].orbit_direction", fa["orbit_direction"],
        sorted({r.get("orbit_direction") for r in rows_of_part(rows, 1)})[0])
    chk("av.bracket.relative_orbit", lb["relative_orbit"], fa["relative_orbit"])
    chk("av.bracket.bursts_per_product", [lb["bursts"], fa["bursts"]], [2, 2])
    chk("av.bracket.platform_in_product_name",
        [lb["name"].split("_")[0], fa["name"].split("_")[0]],
        [lb["platform"], fa["platform"]])
    chk("av.bracket.names_are_slice_products",
        all(x["name"].startswith("S1") and x["name"].endswith(".SAFE")
            for x in (lb, fa)), True)
    chk("av.bracket.days_in_name",
        [x["day"].replace("-", "") for x in (lb, fa)],
        [x["name"].split("_")[5][:8] for x in (lb, fa)])
    chk("av.bracket.online", [lb["online"], fa["online"]], [True, True])

    # --- the nullcheck window must sit inside the bracket -----------------
    nc = avail["nullcheck"]
    chk("av.nullcheck.window_inside_bracket",
        avail["collapse"]["first_after"] > nc["window"]["end"] > nc["window"]["start"]
        > avail["collapse"]["last_before"], True)
    for key in ("products_all", "products_iw_slc", "bursts"):
        chk(f"av.nullcheck.counts.{key}", len(nc[key]), nc["counts"][key])
    chk("av.nullcheck.available",
        any(v > 0 for v in nc["counts"].values()), nc["available"])

    # --- per-part catalogue accounting, rederived from the record ---------
    for p in PART_ORDER:
        fig = avail[f"figure_{p}"]
        prows = rows_of_part(rows, p)
        days = days_of(prows)
        chk(f"av.figure_{p}.name", f"fig_kdlo_{p}_channels", fig["name"])
        chk(f"av.figure_{p}.window", [fig["window"]["start"], fig["window"]["end"]],
            list(prov[p]["window"]))
        chk(f"av.figure_{p}.n_products", len(prows), fig["n_products"])
        chk(f"av.figure_{p}.acquisition_days", days, fig["acquisition_days"])
        chk(f"av.figure_{p}.n_acquisition_days", len(days), fig["n_acquisition_days"])
        chk(f"av.figure_{p}.one_product_per_day", len(days), len(prows))
        chk(f"av.figure_{p}.products_per_day", {d: 1 for d in days},
            fig["products_per_day"])
        chk(f"av.figure_{p}.bursts_per_day", {d: 2 for d in days},
            fig["bursts_per_day"])
        chk(f"av.figure_{p}.n_bursts", sum(fig["bursts_per_day"].values()),
            fig["n_bursts"])
        chk(f"av.figure_{p}.bursts_per_product", fig["n_bursts"], 2 * fig["n_products"])
        chk(f"av.figure_{p}.swaths", fig["swaths"], ["IW1"])
        chk(f"av.figure_{p}.polarisations", fig["polarisations"], ["VH", "VV"])
        chk(f"av.figure_{p}.platforms_sum", sum(fig["platforms"].values()),
            fig["n_products"])
        chk(f"av.figure_{p}.platforms_keys_subset",
            sorted(set(fig["platforms"]) - {"S1A", "S1B", "S1C", "S1D"}), [])
        chk(f"av.figure_{p}.orbits_sum", sum(fig["orbits"].values()), fig["n_products"])
        chk(f"av.figure_{p}.one_geometry", len(fig["orbits"]), 1)
        chk(f"av.figure_{p}.orbit_direction", list(fig["orbits"])[0].split("|")[0],
            sorted({r.get("orbit_direction") for r in prows})[0])
        chk(f"av.figure_{p}.relative_orbit", int(list(fig["orbits"])[0].split("|")[1]),
            lb["relative_orbit"])
        chk(f"av.figure_{p}.missing_repeat_days", grid_missing_days(days),
            fig["missing_repeat_days"])
        chk(f"av.figure_{p}.days_inside_window",
            all(fig["window"]["start"] <= d <= fig["window"]["end"] for d in days), True)

    # --- the long catalogue window: every count must telescope ------------
    ep = avail["epochs"]
    chk("av.epochs.window.start_before_end",
        ep["window"]["start"] < ep["window"]["end"], True)
    chk("av.epochs.covers_both_figure_windows",
        ep["window"]["start"] <= prov[1]["window"][0]
        and ep["window"]["end"] >= prov[2]["window"][1], True)
    chk("av.epochs.covers_all_measured_days",
        all(ep["window"]["start"] <= d <= ep["window"]["end"] for d in days_of(rows)), True)
    chk("av.epochs.n_products_covers_record", ep["n_products"] >= len(rows), True)
    chk("av.epochs.platforms_sum", sum(ep["platforms"].values()), ep["n_products"])
    chk("av.epochs.orbits_sum", sum(ep["orbits"].values()), ep["n_products"])
    chk("av.epochs.archive_online", ep["archive"]["online"], ep["n_products"])
    chk("av.epochs.platforms_keys_subset",
        sorted(set(ep["platforms"]) - {"S1A", "S1B", "S1C", "S1D"}), [])
    chk("av.epochs.orbits_include_record_geometry",
        f"{lb['orbit_direction']}|{lb['relative_orbit']}" in ep["orbits"], True)
    chk("av.epochs.per_year_n_sum", sum(v["n"] for v in ep["per_year"].values()),
        ep["n_products"])
    for year, blk in ep["per_year"].items():
        chk(f"av.epochs.per_year.{year}.platforms_sum",
            sum(blk["platforms"].values()), blk["n"])
        chk(f"av.epochs.per_year.{year}.orbits_sum",
            sum(blk["orbits"].values()), blk["n"])
        chk(f"av.epochs.per_year.{year}.utc_hours_sum",
            sum(blk["utc_hours"].values()), blk["n"])
    for col in ("platforms", "orbits"):
        chk(f"av.epochs.{col}_sum_over_years",
            {k: sum(b[col].get(k, 0) for b in ep["per_year"].values()) for k in ep[col]},
            ep[col])

    # --- archive availability --------------------------------------------
    on = avail["online"]
    chk("av.online.window.start_before_end",
        on["window"]["start"] < on["window"]["end"], True)
    chk("av.online.n_offline", on["n_offline"], len(on["offline"]))
    chk("av.online.needs_lta_restore", on["n_offline"] > 0, on["needs_lta_restore"])
    chk("av.online.n_offline_le_products", on["n_offline"] <= on["n_products"], True)

    # --- the pipeline's own accounting of the record ----------------------
    pl = avail["pipeline"]
    chk("av.pipeline.status", pl["status"], "ok")
    chk("av.pipeline.asset_id", pl["asset_id"], prov[1]["asset_id"])
    chk("av.pipeline.asset_id_eq_part2", pl["asset_id"], prov[2]["asset_id"])
    chk("av.pipeline.asset_name_non_empty", bool(pl["asset_name"].strip()), True)
    chk("av.pipeline.measured_in_figure_windows", pl["measured_in_figure_windows"],
        len(rows))
    chk("av.pipeline.measured_total",
        pl["measured_in_figure_windows"] + pl["measured_outside_figure_windows"],
        pl["measured_total"])
    for p in PART_ORDER:
        key = PARTS[p]["key"]
        pp = pl["parts"][key]
        prows = rows_of_part(rows, p)
        days = days_of(prows)
        chk(f"av.pipeline.{key}.window", pp["window"], list(prov[p]["window"]))
        chk(f"av.pipeline.{key}.n_requests", pp["n_requests"],
            len(prov[p]["request_ids"]))
        chk(f"av.pipeline.{key}.n_requests_ge_used",
            pp["n_requests"] >= len(prov[p]["request_ids_used"]), True)
        chk(f"av.pipeline.{key}.expected_days", pp["expected_days"], len(days))
        chk(f"av.pipeline.{key}.measured_days", pp["measured_days"], len(days))
        chk(f"av.pipeline.{key}.measured_days_eq_n_products",
            pp["measured_days"], len(prows))
        chk(f"av.pipeline.{key}.missing_days", pp["missing_days"], [])
        chk(f"av.pipeline.{key}.metadata_anomaly_days", pp["metadata_anomaly_days"],
            sorted({day_of(r) for r in prows if r.get("pass_label")
                    and r["pass_label"] != PARTS[p]["nominal_pass"]}))
        chk(f"av.pipeline.{key}.days_inside_window",
            all(pp["window"][0] <= d <= pp["window"][1] for d in days), True)
    chk("av.pipeline.parts.expected_days_total",
        sum(pl["parts"][PARTS[p]["key"]]["expected_days"] for p in PART_ORDER),
        pl["measured_in_figure_windows"])
    chk("av.pipeline.parts.windows_are_the_record_windows",
        [pl["parts"][PARTS[p]["key"]]["window"] for p in PART_ORDER],
        [list(prov[p]["window"]) for p in PART_ORDER])

    return out




def fem_checks(fem):
    """Arithmetic + cross-reference checks of the FEM eigenmode file.

    The solver output itself is not reproducible here (no solver in this
    package); what is checked is that the file's own numbers hang together, that
    the mast model is the *same* mast in all three paths, and that the honest
    caveats the KDLO analysis recorded are still present verbatim.
    """
    out = []

    def chk(path, got, want):
        out.append((path, got, want))

    a = fem["assumptions"]
    chk("fem.asset_id_non_empty", bool(fem["asset_id"].strip()), True)
    chk("fem.site_non_empty", bool(fem["site"].strip()), True)
    chk("fem.solver", fem["solver"],
        "tower_backend::ec_solver::solver_engine::solve_eigenmodes")
    chk("fem.solver_solves_eigenmodes", fem["solver"].endswith("solve_eigenmodes"),
        True)
    chk("fem.payload_source_names_the_project",
        "kdlo_tower_analysis" in fem["payload_source"], True)

    # --- tier table: a closed radius chain summing to the 519.7 m mast ----
    tiers = a["tier_table"]
    chk("fem.assumptions.tier_table.n_tiers", len(tiers), 6)
    chk("fem.assumptions.tier_table.height_sum_m",
        round(sum(t["height_m"] for t in tiers), 4), 519.7)
    chk("fem.assumptions.tier_table.radius_chain",
        [t["top_radius_m"] for t in tiers[:-1]],
        [t["base_radius_m"] for t in tiers[1:]])
    chk("fem.assumptions.tier_table.radii_shrink",
        all(t["top_radius_m"] < t["base_radius_m"] for t in tiers), True)
    chk("fem.assumptions.tier_table.leg_count",
        sorted({t["leg_count"] for t in tiers}), [3])
    chk("fem.assumptions.tier_table.profiles_non_empty",
        all(t["leg_profile"].strip() and t["bracing_profile"].strip()
            for t in tiers), True)
    chk("fem.assumptions.tier_table.heights_positive",
        all(t["height_m"] > 0 for t in tiers), True)

    # --- assumption grid -------------------------------------------------
    chk("fem.assumptions.guy_levels_m", a["guy_levels_m"],
        [60.0, 120.0, 240.0, 360.0, 470.0])
    chk("fem.assumptions.guy_levels_ascending",
        len(set(a["guy_levels_m"])) == len(a["guy_levels_m"]) == 5
        and a["guy_levels_m"] == sorted(a["guy_levels_m"]), True)
    chk("fem.assumptions.guy_levels_below_mast_top",
        max(a["guy_levels_m"]) < sum(t["height_m"] for t in tiers), True)
    chk("fem.assumptions.anchor_ratios", a["anchor_ratios"], [0.5, 0.65, 0.8])
    chk("fem.assumptions.cable_areas_m2", a["cable_areas_m2"],
        [0.0002, 0.00035, 0.0005])
    chk("fem.assumptions.cable_youngs_modulus_pa", a["cable_youngs_modulus_pa"],
        1.6e11)
    chk("fem.assumptions.solidity_ratio", a["solidity_ratio"], 0.25)
    chk("fem.assumptions.wind_pressure_n_m2", a["wind_pressure_n_m2"], 800.0)
    chk("fem.assumptions.element_discretization_applied",
        a["element_discretization_applied"], False)
    chk("fem.assumptions.element_discretization_m", a["element_discretization_m"],
        5.0)
    chk("fem.assumptions.guy_data_source_is_an_assumption",
        "assumption" in a["guy_data_source"], True)

    # --- path A: one cantilever per tier, stored at 4 significant digits --
    pa = fem["path_a_per_tier"]
    chk("fem.path_a.n_tiers", len(pa), len(tiers))
    chk("fem.path_a.segment_index", [e["segment_index"] for e in pa],
        list(range(len(pa))))
    chk("fem.path_a.nodes_per_tier", sorted({e["nodes"] for e in pa}), [6])
    chk("fem.path_a.elements_per_tier", sorted({e["elements"] for e in pa}), [9])
    chk("fem.path_a.labels", [e["label"] for e in pa],
        [f"Mast 0-520 m-tier-{i}" for i in range(len(pa))])
    chk("fem.path_a.stored_is_rounded_fundamental",
        [e["stored_fundamental_hz"] for e in pa],
        [round(e["fundamental_hz"], 4) for e in pa])
    chk("fem.path_a.fundamentals_positive",
        all(e["fundamental_hz"] > 0 for e in pa), True)
    chk("fem.path_a.mast_label_rounds_model_height",
        ("520" in pa[0]["label"])
        and round(sum(t["height_m"] for t in tiers), 1) == 519.7, True)

    # --- path B: the whole mast as one model (production modelling path) --
    pb = fem["path_b_whole_mast"]
    chk("fem.path_b.description_non_empty", bool(pb["description"].strip()), True)
    chk("fem.path_b.description_names_all_tiers",
        "all 6 tiers" in pb["description"], True)
    chk("fem.path_b.elements_is_the_whole_tier_table",
        pb["elements"], 9 * len(pa))
    chk("fem.path_b.f_hz.n_modes", len(pb["f_hz"]), 6)
    chk("fem.path_b.f_hz.strictly_ascending",
        pb["f_hz"] == sorted(pb["f_hz"]) and len(set(pb["f_hz"])) == 6, True)
    chk("fem.path_b.f_hz.positive", all(v > 0 for v in pb["f_hz"]), True)
    chk("fem.path_b.f1_within_self_supporting_span",
        min(e["f1_self_supporting_hz"] for e in fem["path_c_guyed"]) <= pb["f_hz"][0]
        <= max(e["f1_self_supporting_hz"] for e in fem["path_c_guyed"]), True)
    chk("fem.path_b.lumped_mass_t_positive", pb["lumped_mass_t"] > 0, True)
    chk("fem.path_b.nodes_positive", pb["nodes"] > 0, True)
    chk("fem.path_b.nodes_lt_elements", pb["nodes"] < pb["elements"], True)

    # --- path C: the 3x3 guy/anchor/cable grid ---------------------------
    pc = fem["path_c_guyed"]
    combos = [(r, c) for r in a["anchor_ratios"] for c in a["cable_areas_m2"]]
    chk("fem.path_c.n_combinations", len(pc), len(combos))
    chk("fem.path_c.anchor_ratios", [e["anchor_ratio"] for e in pc],
        [r for r, _ in combos])
    chk("fem.path_c.cable_areas_m2", [e["cable_area_m2"] for e in pc],
        [c for _, c in combos])
    chk("fem.path_c.grid_is_complete",
        sorted({(e["anchor_ratio"], e["cable_area_m2"]) for e in pc}),
        sorted(set(combos)))
    chk("fem.path_c.total_guys", sorted({e["total_guys"] for e in pc}),
        [3 * len(a["guy_levels_m"])])
    chk("fem.path_c.active_guys", sorted({e["active_guys"] for e in pc}), [10])
    chk("fem.path_c.active_guys_le_total",
        all(e["active_guys"] <= e["total_guys"] for e in pc), True)
    chk("fem.path_c.mast_elements", sorted({e["mast_elements"] for e in pc}),
        [9 * len(pa)])
    chk("fem.path_c.n_modes_le_nodes",
        max(len(e["f_hz"]) for e in pc) <= min(e["nodes"] for e in pc), True)
    chk("fem.path_c.n_modes", sorted({len(e["f_hz"]) for e in pc}), [6])
    chk("fem.path_c.f_hz.strictly_ascending",
        all(e["f_hz"] == sorted(e["f_hz"]) and len(set(e["f_hz"])) == 6 for e in pc),
        True)
    chk("fem.path_c.f1_raw_is_first_mode",
        [e["f1_raw_hz"] for e in pc], [e["f_hz"][0] for e in pc])
    chk("fem.path_c.guys_stiffen_the_mast",
        all(e["f1_self_supporting_hz"] < e["f1_raw_hz"] < e["f1_all_cables_active_hz"]
            for e in pc), True)
    chk("fem.path_c.tension_only_iterations",
        sorted({e["tension_only_iterations"] for e in pc}), [2])
    chk("fem.path_c.sway_mode_index_all_null",
        [e["sway_mode_index"] for e in pc], [None] * len(pc))
    chk("fem.path_c.f1_sway_hz_all_null",
        [e["f1_sway_hz"] for e in pc], [None] * len(pc))
    chk("fem.path_c.anchor_inclination_deg_in_range",
        all(0.0 < e["anchor_inclination_from_horizontal_deg"] < 90.0 for e in pc),
        True)
    chk("fem.path_c.anchor_inclination_by_ratio",
        {r: sorted({e["anchor_inclination_from_horizontal_deg"] for e in pc
                    if e["anchor_ratio"] == r}) for r in a["anchor_ratios"]},
        {0.5: [63.43494882292201], 0.65: [56.97613244420336],
         0.8: [51.34019174590991]})
    incl = [sorted({e["anchor_inclination_from_horizontal_deg"] for e in pc
                    if e["anchor_ratio"] == r})[0] for r in a["anchor_ratios"]]
    chk("fem.path_c.inclination_falls_with_ratio", incl,
        sorted(incl, reverse=True))

    # --- the honest caveats must survive ---------------------------------
    quirks = fem["known_model_quirks"]
    chk("fem.known_model_quirks.n", len(quirks), 4)
    chk("fem.known_model_quirks.non_empty_strings",
        all(isinstance(s, str) and s.strip() for s in quirks), True)
    chk("fem.known_model_quirks.note_missing_guys",
        any("guy cables" in s for s in quirks), True)
    chk("fem.known_model_quirks.note_element_discretization",
        any("element_discretization" in s for s in quirks), True)
    chk("fem.known_model_quirks.note_per_tier_cantilever",
        any("fixed at its own base" in s for s in quirks), True)
    chk("fem.known_model_quirks.note_antenna_mass",
        any("200 kg" in s for s in quirks), True)
    chk("fem.known_model_quirks.note_tier_scope",
        any("regenerate_segments" in s for s in quirks), True)

    return out

