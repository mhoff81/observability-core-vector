#!/usr/bin/env python3
"""OCV-Paper — data layer.

Reads only ``data/lumo/lumo_ocv_channels.csv`` (the copy of ``lumo_channels.csv``
extended with A/D/F/S/P) plus the five committed reference JSONs in
``data/lumo/reference/``. No burst cache, no third-party imports — but the same
ordering, the same field names and the same conventions as the LUMO project
analysis scripts (see ``lumo_ocv_stats.py``).
"""
from __future__ import annotations

import csv
import hashlib
import json
import os

import numpy as np

import lumo_ocv_stats as st

# Layout: this package lives in ``code/lumo/``; its inputs and outputs are in
# ``data/lumo/`` and its figures in ``figures/lumo/``. The per-site subfolders
# keep the LUMO package and the Bautzen package (``code/bautzen/``) apart — see
# the root README.
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, os.pardir, os.pardir, "data", "lumo")
REF = os.path.join(DATA, "reference")
FIGDIR = os.path.join(HERE, os.pardir, os.pardir, "figures", "lumo")

CSV_PATH = os.path.join(DATA, "lumo_ocv_channels.csv")
CSV_META_PATH = os.path.join(DATA, "lumo_ocv_channels_meta.json")

# The EMV vector and the feature order are taken from the origin
# (analyze_lumo_echo_mask_vector.DIM_LABEL, analyze_lumo_emv_multivariate_cv.FEATURES,
# analyze_lumo_emv_pairwise_triplets_cv.MODELS) so that numbers and labels stay
# comparable with the committed results.
DIMS = ["A", "D", "F", "S", "P"]
DIM_LABEL = {
    "A": "Area (n_masked, number of echo pixels)",
    "D": "Density (n_masked / bbox_area)",
    "F": "Fragmentation (number of connected components, 8-neighbourhood)",
    "S": "Shift (distance mask centroid <-> peak pixel, px)",
    "P": "Persistence (mean pixel frequency across all dates inside the own mask)",
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

# The five reference JSONs are German-language copies of the LUMO project
# outputs (their `per_dimension[*].label`, `emv_definition[*]`, `models[*]` and
# `per_model[*].label` strings are German). This repo is written in English, so
# the pin block compares those label strings through the following explicit and
# exhaustive map: every committed German label is listed here, and a label that
# is missing from the map would show up as a pin failure. Numbers and keys are
# always compared byte-for-byte; only these doc strings are translated.
LABEL_EN = {
    "Area (n_masked, Anzahl Echo-Pixel)":
        "Area (n_masked, number of echo pixels)",
    "Fragmentation (Anzahl zusammenhängender Komponenten, 8-Nachbarschaft)":
        "Fragmentation (number of connected components, 8-neighbourhood)",
    "Shift (Distanz Maskenschwerpunkt <-> Peak-Pixel, px)":
        "Shift (distance mask centroid <-> peak pixel, px)",
    "Persistence (mittlere Pixel-Häufigkeit über alle Termine innerhalb der "
    "eigenen Maske)":
        "Persistence (mean pixel frequency across all dates inside the own mask)",
    "[P, D, S] (Maskenmorphologie-Kontrolle)":
        "[P, D, S] (mask-morphology control)",
    "[gamma2, P, D] (Hauptwahl x3)":
        "[gamma2, P, D] (headline, 3 dims)",
    "[gamma2, P, D, A, F, S] (voller 6D-Vektor)":
        "[gamma2, P, D, A, F, S] (full 6D vector)",
}

LAMBDA_GRID = [0.0, 0.25, 0.5, 0.75, 1.0]
HEADLINE_LAMBDA = 0.5
N_PERM_4CLASS = 1000
N_PERM_PAIR = 500
N_BOOT = 2000
RNG_SEED = 11

STATE_ORDER = st.STATE_ORDER
SHORT = st.STATE_SHORT
# The committed `lumo_tower_coherence_states.json` uses the short form of
# `analyze_lumo_tower_coherence.STATE_SHORT` ("DAM0(healthy)") — the pin block
# against that file needs exactly that spelling.
SHORT_ALC = {"healthy": "DAM0(healthy)", "DAM 3": "DAM3",
             "DAM 4": "DAM4", "DAM 6": "DAM6"}

REFERENCE_FILES = {
    "tower_coherence_states": "lumo_tower_coherence_states.json",
    "gamma2_pairs": "lumo_gamma2_pairs.json",
    "echo_mask_vector": "lumo_echo_mask_vector.json",
    "emv_multivariate_cv": "lumo_emv_multivariate_cv.json",
    "emv_pairwise_triplets_cv": "lumo_emv_pairwise_triplets_cv.json",
}

# Columns appended by the generator (duplicate check while reading)
NEW_COLS = ["A", "D", "F", "S", "P",
            "mask_rows", "mask_cols", "bbox_area",
            "centroid_r", "centroid_c", "peak_r", "peak_c"]


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


def load_csv(path=CSV_PATH):
    """The 178 coherence overpasses in committed order."""
    with open(path, newline="") as fh:
        rows = []
        for r in csv.DictReader(fh):
            if not (r.get("coherence_gamma2") or "").strip():
                continue
            rows.append({
                "burst_id": r["burst_id"],
                "date": r["acquisition_date"],
                "month": r["month"],
                "orbit": r["orbit_direction"],
                "damage_label": r["damage_label"],
                "wind_speed_ms": _num(r.get("wind_speed_ms")),
                "gamma2": _num(r.get("coherence_gamma2")),
                "n_masked": _int(r.get("coherence_masked_pixels")),
                "peak_intensity": _num(r.get("peak_intensity")),
                "A": _num(r.get("A")),
                "D": _num(r.get("D")),
                "F": _num(r.get("F")),
                "S": _num(r.get("S")),
                "P": _num(r.get("P")),
                "mask_rows": _int(r.get("mask_rows")),
                "mask_cols": _int(r.get("mask_cols")),
                "bbox_area": _int(r.get("bbox_area")),
                "centroid_r": _num(r.get("centroid_r")),
                "centroid_c": _num(r.get("centroid_c")),
                "peak_r": _int(r.get("peak_r")),
                "peak_c": _int(r.get("peak_c")),
            })
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


def n_by_state(rows):
    return {SHORT[l]: int(sum(1 for r in rows if r["damage_label"] == l))
            for l in STATE_ORDER}


def feature_dataset(rows):
    """X, y of the complete cases (all 6 features present) — 177 rows.

    Verbatim ``mcv.build_dataset``: only complete-case rows, in committed
    order, feature order ``FEATURES``.
    """
    complete = [r for r in rows if all(r.get(f) is not None for f in FEATURES)]
    X = np.array([[float(r[f]) for f in FEATURES] for r in complete])
    y = np.array([r["damage_label"] for r in complete], dtype=object)
    return X, y, complete
