#!/usr/bin/env python3
"""OCV-Paper (Bautzen bridge site) — data layer and constants.

Reads only ``data/bautzen/bautzen_ocv_channels.csv`` (the per-date table of the
four rect series), the two committed reference JSONs in
``data/bautzen/reference/`` and the hourly environment series
``data/bautzen/openmeteo_bautzen_2025.json``. No burst cache, no
network: the four rect text files committed next to the CSV are the only raw
input and are enough to rebuild the whole site offline (much smaller than the
LUMO burst cache, which is why they are committed here).

The site, mirroring the LUMO package (only the mask layer is swapped, see
``bautzen_ocv_masks.py``):

  series   ASC/DESC bridge and ASC/DESC control area 500 m east of the bridge
           (the control area rules out that the deck band is a processing
           artefact)
  states   RS (reference) / DS1 / DS2, assigned by *overpass time*, not by date:
           both state boundaries fall in the middle of a day (ASC flies at
           18:51 local, DESC at 07:08 local), so a date-only rule would assign
           two dates wrongly -- DESC 12 May (07:08, before the intervention) is
           RS, ASC 29 Sep (18:51, after it) is DS2.
  channels gamma2 (deck band, floor-corrected), A, D, F, S, P of the deck-edge
           mask, exactly the OCV definition of the LUMO package

Series labels are the committed German spellings of the reference JSON (they are
the keys of that file and therefore of the pin block); ``SERIES_EN`` holds the
English display names used in the report.
"""
from __future__ import annotations

import csv
import datetime as dt
import hashlib
import json
import os

import numpy as np

import bautzen_ocv_masks as masks
import bautzen_ocv_stats as st

# Layout: this package lives in ``code/bautzen/``; its inputs and outputs are in
# ``data/bautzen/`` and its figures in ``figures/bautzen/``. This mirrors the
# LUMO package in ``code/lumo/`` (data/lumo/, figures/lumo/) — see the root
# README.
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, os.pardir, os.pardir, "data", "bautzen")
REF = os.path.join(DATA, "reference")
FIGDIR = os.path.join(HERE, os.pardir, os.pardir, "figures", "bautzen")

CSV_PATH = os.path.join(DATA, "bautzen_ocv_channels.csv")
CSV_META_PATH = os.path.join(DATA, "bautzen_ocv_channels_meta.json")
METEO_PATH = os.path.join(DATA, "openmeteo_bautzen_2025.json")

SITE = "bautzen_openbridge"

# ── The four rect series (committed keys of the reference JSON) ─────────────
SERIES = {
    "ASC Brücke": "bautzen_rects_asc_iw2.txt",
    "DESC Brücke": "bautzen_rects_desc_iw3.txt",
    "ASC Kontrolle 500 m östlich": "bautzen_rects_asc_ctrl_500m_east.txt",
    "DESC Kontrolle 500 m östlich": "bautzen_rects_ctrl_500m_east.txt",
}
SERIES_ORDER = list(SERIES)
SERIES_EN = {
    "ASC Brücke": "ASC bridge",
    "DESC Brücke": "DESC bridge",
    "ASC Kontrolle 500 m östlich": "ASC control 500 m east",
    "DESC Kontrolle 500 m östlich": "DESC control 500 m east",
}
BRIDGE_BY_GEOM = {"ASC": "ASC Brücke", "DESC": "DESC Brücke"}
CONTROL_BY_GEOM = {"ASC": "ASC Kontrolle 500 m östlich",
                   "DESC": "DESC Kontrolle 500 m östlich"}

# ── State windows (from datasheet.txt + acquisition times of the manifests) ─
CUT_DS1 = dt.datetime(2025, 5, 12, 10, 0)    # DS1 implemented, 10:00-12:15 local
CUT_DS2 = dt.datetime(2025, 9, 29, 10, 0)    # DS2, 10:00-12:15 local
OVERPASS_LOCAL = {"ASC": (18, 51), "DESC": (7, 8)}   # from the manifest times
STATES = ["RS", "DS1", "DS2"]
STATE_ORDER = st.STATE_ORDER

# Documented state windows (local time, CEST) — the assignment rule of state_of().
STATE_WINDOWS = {
    "RS": "until 2025-05-12 10:00 local (before DS1 was installed)",
    "DS1": "2025-05-12 10:00 .. 2025-09-29 10:00 local",
    "DS2": "from 2025-09-29 10:00 local",
}


# ── Channels ────────────────────────────────────────────────────────────────
DIMS = ["A", "D", "F", "S", "P"]
DIM_LABEL = {
    "A": "Area (n_masked, number of deck-edge pixels)",
    "D": "Density (n_masked / bbox_area)",
    "F": "Fragmentation (number of connected components, 8-neighbourhood)",
    "S": "Shift (distance mask centroid <-> deck band peak, px)",
    "P": "Persistence (mean pixel frequency across the dates of the series "
         "inside the own mask)",
}
FEATURES = ["gamma2", "P", "D", "A", "F", "S"]
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

LAMBDA_GRID = [0.0, 0.25, 0.5, 0.75, 1.0]
HEADLINE_LAMBDA = 0.5
N_PERM_3CLASS = 1000
N_PERM_PAIR = 500
N_BOOT = st.N_BOOT
RNG_SEED = st.RNG_SEED
BOOT_SEED = st.BOOT_SEED

# ── Environment covariates (mirror of the origin's section 5c) ──────────────
ENV_VARS = (("temperature_2m", "T", "°C"), ("relative_humidity_2m", "RH", "%"))

REFERENCE_FILES = {
    "deck_edge_mask": "bautzen_deck_edge_mask.json",
    "deck_edge_state": "bautzen_deck_edge_state.json",
}

# The two reference JSONs are German-language outputs of the bridge project:
# their keys are the German series names, and ``config.*`` as well as
# ``metrics[*].gamma2.note`` / ``.floor_correction`` carry German documentation
# strings. This package is written in English, so every comparison of such a
# string goes through the following explicit and **exhaustive** map: it lists
# every German string that occurs in ``bautzen_deck_edge_state.json`` (the mask
# JSON holds no strings but the series names), and a string that is missing
# from the map would show up as a pin failure instead of being silently
# accepted. Numbers and keys are always compared byte for byte; only these doc
# strings are translated.
LABEL_EN = {
    # config.mask_source, config.gamma2_detail.*
    "analyze_bautzen_deck_edge_mask.py (unverändert übernommen)":
        "analyze_bautzen_deck_edge_mask.py (taken over unchanged)",
    "|sum(z)|^2 / (N * sum(|z|^2)) ueber die Maskenpixel":
        "|sum(z)|^2 / (N * sum(|z|^2)) over the mask pixels",
    "(N*gamma2 - 1)/(N - 1), inkohaerentes Nullmodell (zirkular-Gauss, "
    "unabhaengige Looks); keine Klemmung, negative Werte bleiben erhalten":
        "(N*gamma2 - 1)/(N - 1), incoherent null model (circular Gaussian, "
        "independent looks); no clamping, negative values are kept",
    "komplementaerer Kanal, nicht unabhaengig":
        "complementary channel, not independent",
    "Deckbandmaske = strukturell attribuiert (Deckkanten-Aussage)":
        "deck-band mask = structurally attributed (deck-edge statement)",
    "Haus-/Peakmaske = signalbasiert (|z|^2 >= 0.3*Peak and >= 5.0*Median, "
    "min 2 px), NICHT strukturell attribuiert -> keine Deckkanten-Aussage":
        "house/peak mask = signal-based (|z|^2 >= 0.3*peak and >= 5.0*median, "
        "min 2 px), NOT structurally attributed -> no deck-edge statement",
    # config.environment_source.* and gamma2_environment[*].interpolation
    "linear zwischen den Nachbarstunden des Überflugzeitpunkts "
    "(keine Stundenrundung)":
        "linear between the neighbouring hours of the overpass time "
        "(no hour rounding)",
    "linear zwischen den Nachbarstunden der Überflugzeit (18:51 lokal)":
        "linear between the neighbouring hours of the overpass time "
        "(18:51 local)",
    "linear zwischen den Nachbarstunden der Überflugzeit (07:08 lokal)":
        "linear between the neighbouring hours of the overpass time "
        "(07:08 local)",
    "Gegenprobe zur Saisondeutung: Zustand und Jahreszeit sind konfundiert, "
    "T/RH sind der einfachste Saison-Ersatz":
        "counter-check for the seasonal reading: state and season are "
        "confounded, T/RH are the simplest season proxy",
    # metrics[*].gamma2 -- the two strings state_metrics() states itself
    "gamma2_band_* = strukturell attribuierte Deckbandmaske "
    "(Deckkanten-Aussage); gamma2_house_* = signalbasierte Haus-/Peakmaske, "
    "nicht strukturell attribuiert (keine Deckkanten-Aussage)":
        masks.GAMMA2_NOTE,
    "(N*raw - 1)/(N - 1) unter dem inkohaerenten Nullmodell; negative Werte "
    "bleiben erhalten und bedeuten 'unter dem endlichen-N-Rauschboden'":
        masks.GAMMA2_FLOOR_CORRECTION,
}


def translate(value):
    """Translate a committed doc string (recursively) through ``LABEL_EN``.

    Non-strings, and strings that are not in the map (dates, keys, file names,
    numbers written as text), are returned unchanged — an unmapped German doc
    string therefore does not vanish, it mismatches and becomes a pin failure.
    """
    if isinstance(value, dict):
        return {k: translate(v) for k, v in value.items()}
    if isinstance(value, list):
        return [translate(v) for v in value]
    if isinstance(value, str):
        return LABEL_EN.get(value, value)
    return value

# Columns the generator appends to the enumerated table (duplicate check while
# reading: a column that appears twice would silently overwrite a value).
NEW_COLS = ["A", "D", "F", "S", "P", "P_frac_ge_50pct", "P_frac_ge_30pct",
            "mask_rows", "mask_cols", "bbox_area", "row_span", "col_span",
            "gamma2_band_raw", "gamma2_band_floor_corrected", "n_band",
            "gamma2_house_raw", "gamma2_house_floor_corrected", "n_house",
            "mask_pixels"]

# The column that carries the site's coherence channel of the OCV.
GAMMA2_COLUMN = "gamma2_band_floor_corrected"


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _num(x):
    x = (x or "").strip()
    return float(x) if x else None


def _int(x):
    v = _num(x)
    return None if v is None else int(v)


def _bool(x):
    x = (x or "").strip().lower()
    return None if not x else x == "true"


def pixels_to_mask(text, shape):
    """``"r:c;r:c;..."`` (row-major) -> full-chip bool mask of ``shape``.

    Inverse of the ``mask_pixels`` column of the generated table: it lets the
    figure pipeline rebuild the per-state persistence maps (and therefore
    ``P``, the unions and the Jaccard overlaps) from the CSV alone, without the
    rect files. Pixel order is irrelevant, the mask is a set.
    """
    m = np.zeros(tuple(int(s) for s in shape), dtype=bool)
    for tok in (text or "").split(";"):
        tok = tok.strip()
        if not tok:
            continue
        r, c = tok.split(":")
        m[int(r), int(c)] = True
    return m


INT_COLS = ["peak_row", "peak_row_offset", "n_masked", "n_masked_house",
            "n_components", "largest_component", "band_peak_r", "band_peak_c",
            "bbox_area", "row_span", "col_span", "mask_rows", "mask_cols",
            "n_band", "n_house"]
FLOAT_COLS = DIMS + ["band_contrast", "centroid_row", "centroid_col",
                     "P_frac_ge_50pct", "P_frac_ge_30pct",
                     "gamma2_band_raw", GAMMA2_COLUMN,
                     "gamma2_house_raw", "gamma2_house_floor_corrected"]


def load_csv(path=CSV_PATH):
    """The four rect series, date by date, in committed file order."""
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        fields = reader.fieldnames or []
        dupes = sorted({c for c in fields if fields.count(c) > 1})
        if dupes:
            raise ValueError(f"{path}: duplicate columns {dupes}")
        rows = []
        for r in reader:
            if not (r.get("acquisition_date") or "").strip():
                continue
            row = {
                "series": r["series"],
                "rect_file": r.get("rect_file") or SERIES[r["series"]],
                "date": r["acquisition_date"],
                "month": r["month"],
                "state": r["state"],
                # The ported LUMO estimators address the state through
                # `damage_label`; this alias keeps those ports unchanged.
                "damage_label": r["state"],
                "orbit": r["orbit_direction"],
                "has_mask": _bool(r.get("has_mask")),
            }
            for col in INT_COLS:
                row[col] = _int(r.get(col))
            for col in FLOAT_COLS:
                row[col] = _num(r.get(col))
            # The OCV channel of this site is the floor-corrected deck-band
            # gamma^2; the raw estimator and the house/peak variant stay in the
            # row as audit columns (see the report).
            row["gamma2"] = row[GAMMA2_COLUMN]
            # The mask itself, as committed by the generator (see pixels_to_mask):
            # it is what the report panel A draws and what the per-state
            # persistence maps are rebuilt from.
            row["mask_pixels"] = r.get("mask_pixels") or ""
            row["mask"] = pixels_to_mask(row["mask_pixels"],
                                         (row["mask_rows"], row["mask_cols"]))
            rows.append(row)
    return rows


def load_csv_meta(path=CSV_META_PATH):
    with open(path) as fh:
        return json.load(fh)


def load_reference():
    out = {}
    for name, fname in REFERENCE_FILES.items():
        with open(os.path.join(REF, fname)) as fh:
            out[name] = json.load(fh)
    return out


def load_meteo(path=METEO_PATH):
    """Hourly environment series -> (local timestamps, {tag: float array})."""
    if not os.path.exists(path):
        return None
    with open(path) as fh:
        d = json.load(fh)
    h = d.get("hourly") or {}
    if "time" not in h:
        return None
    ts = [dt.datetime.strptime(t, "%Y-%m-%dT%H:%M") for t in h["time"]]
    cols = {tag: np.asarray([np.nan if v is None else float(v)
                             for v in h.get(name, [None] * len(ts))], float)
            for name, tag, _u in ENV_VARS}
    return ts, cols


def acq_dt(date, geom):
    """Acquisition instant (local) of one date from date + overpass time."""
    y, m, d = (int(p) for p in date.split("-"))
    return dt.datetime(y, m, d, *OVERPASS_LOCAL[geom])


def state_of(date, geom):
    """RS / DS1 / DS2 -- time-correct, not date-based."""
    t = acq_dt(date, geom)
    if t < CUT_DS1:
        return "RS"
    if t < CUT_DS2:
        return "DS1"
    return "DS2"


def geom_of(label):
    return "ASC" if label.startswith("ASC") else "DESC"


def is_bridge(label):
    return "Brücke" in label


def month_of(date):
    return date[:7]


def by_series(rows):
    return {lab: [r for r in rows if r["series"] == lab] for lab in SERIES}


def n_by_state(rows):
    return {s: int(sum(1 for r in rows if r["state"] == s)) for s in STATES}


def feature_dataset(rows, series=None):
    """X, y of the complete cases (all six channels present) of one series."""
    sub = [r for r in rows if series is None or r["series"] == series]
    complete = [r for r in sub if all(r.get(f) is not None for f in FEATURES)]
    X = np.array([[float(r[f]) for f in FEATURES] for r in complete])
    y = np.array([r["state"] for r in complete], dtype=object)
    return X, y, complete

