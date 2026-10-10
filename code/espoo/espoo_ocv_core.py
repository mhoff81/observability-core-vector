#!/usr/bin/env python3
"""Espoo OCV-Paper — data layer.

Reads only ``data/espoo/espoo_ocv_channels.csv`` (the copy of the committed site
record table ``espoo_channels.csv`` extended with the derived bookkeeping columns
``A, date, month, doy, phase_available`` by ``espoo_ocv_channels_csv.py``) plus
the three committed reference JSONs in ``data/espoo/reference/``. No burst cache,
no network, no third-party imports beyond numpy — but the same ordering, the same
field names and the same conventions as the LUMO package (see
``espoo_ocv_stats.py``).

The site is **event-free**: ``espoo_channels.csv`` carries no
``damage_label`` / ``structural_state`` / ``condition_label`` on any of its 150
rows, and the upstream ``espoo_mast_observability.py`` says so in its own verdict
(``"This site has no ground-truth state label"``). The two groups this package
splits into are therefore the **orbit geometries** of the record — ASCENDING and
DESCENDING — and ``state_of()`` is the whole of the state definition.
"""
from __future__ import annotations

import csv
import datetime as dt
import hashlib
import json
import os

import numpy as np

import espoo_ocv_stats as st

# Layout: this package lives in ``code/espoo/``; its inputs and outputs are in
# ``data/espoo/`` and its figures in ``figures/espoo/``.
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, os.pardir, os.pardir, "data", "espoo")
REF = os.path.join(DATA, "reference")
FIGDIR = os.path.join(HERE, os.pardir, os.pardir, "figures", "espoo")

# The committed site record table this package copies (the raw input of the
# generator) and the extended table the figure pipeline actually reads.
RAW_CSV_PATH = os.path.join(DATA, "espoo_channels.csv")
CSV_PATH = os.path.join(DATA, "espoo_ocv_channels.csv")
CSV_META_PATH = os.path.join(DATA, "espoo_ocv_channels_meta.json")

REFERENCE_FILES = {
    "observability": "espoo_mast_observability.json",
    "stationarity": "espoo_phase_stationarity.json",
    "channels": "fig_espoo_channels.json",
}

# ---------------------------------------------------------------------------
# The observability core vector of this site
# ---------------------------------------------------------------------------
# The record delivers exactly one mask-geometry quantity, ``A`` =
# ``coherence_masked_pixels`` — the size in px of the coherent echo mask the
# pipeline found around the mast. The export carries no mask raster
# (``peak_intensity``, ``sub_aperture_brightness``, ``dwell_s``, ``scatterers``
# and every ``peak_*`` column are empty on all 150 rows), so there is **no
# echo-mask layer** of this package and no D/F/S/P to compute. ``A`` is copied,
# not recomputed, and the generator proves the copy is the identity row by row.
DIMS = ["A"]
DIM_LABEL = {
    "A": "Area (A = coherence_masked_pixels, number of echo pixels)",
}
# The vector the figure pipeline models: the two measured quantities of the
# observability core vector of this site, plus the two co-variates panel D asks
# about.
FEATURES = ["gamma2", "A"]
COVARIATES = ["wind_speed_ms", "temperature_c"]
MODELS = {"gamma2": ["gamma2"], "A": ["A"], "gamma2A": ["gamma2", "A"],
          "x_cov": ["gamma2", "A", "wind_speed_ms", "temperature_c"]}
MODEL_LABEL = {
    "gamma2": "gamma2 (1D)",
    "A": "A (1D, the mask size alone)",
    "gamma2A": "[gamma2, A] (headline, 2 dims)",
    "x_cov": "[gamma2, A, wind, temperature] (headline + co-variates)",
}
HEADLINE_PAIR = "gamma2A"
BASELINE = "gamma2"
SECONDARY_PAIR = "A"

# Every phase-ladder quantity of the record is optional and intermittent: the
# phase scan covers **80 of the 150** acquisitions (53.3 %), ``phase_snr_db``
# none of them. ``PHASE_DIMS`` is therefore reported as *availability*, never as
# an observability channel — see panel E4 and §A4 of the report.
PHASE_DIMS = ["phase_coherence", "phase_rms_rad", "phase_snr_db"]

# The two observability metrics the whole package is built on. ``COHERENCE_COLS``
# maps the OCV name to the committed record column so that the pin block can show
# the identity.
COHERENCE_COLS = {"gamma2": "coherence_gamma2", "A": "coherence_masked_pixels"}

SHORT = st.STATE_SHORT
STATE_ORDER = st.STATE_ORDER

LAMBDA_GRID = [0.0, 0.25, 0.5, 0.75, 1.0]
HEADLINE_LAMBDA = 0.5
N_PERM_2CLASS = 1000
N_PERM_PAIR = 1000
N_BOOT = 10000
N_BOOT_STRATA = 2000
RNG_SEED = 7

# The out-of-time split of panel D: the record is cut in half by acquisition
# order and each half is scored with the model fitted on the other one.
OOT_FRACTION = 0.5


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


def _ts(x):
    x = (x or "").strip()
    return dt.datetime.fromisoformat(x) if x else None


def _day(ts):
    """The acquisition day as an ordinal (the day, not the second)."""
    return ts.date().toordinal() if ts else None


def state_of(row):
    """The state of a row: its **orbit geometry**, nothing else.

    The record has no ground-truth label of any kind, so this one-liner is the
    entire state definition of the Espoo package.
    """
    return row["orbit"]


def load_csv(path=CSV_PATH):
    """The 150 acquisitions in committed order."""
    with open(path, newline="") as fh:
        rows = []
        for r in csv.DictReader(fh):
            if not (r.get("coherence_gamma2") or "").strip():
                continue
            ts = _ts(r.get("acquisition_ts"))
            rows.append({
                "id": r["id"],
                "ts": ts,
                "date": r.get("acquisition_date") or (ts.date().isoformat() if ts else None),
                "month": r.get("month") or (ts.strftime("%Y-%m") if ts else None),
                "day": _day(ts),
                "doy": _int(r.get("doy")),
                "pass": r["pass_label"],
                "orbit": r["orbit_direction"],
                "wind_speed_ms": _num(r.get("wind_speed_ms")),
                "temperature_c": _num(r.get("temperature_c")),
                "traffic_load_label": r.get("traffic_load_label"),
                "status": r.get("status"),
                "brightness_ratio": _num(r.get("brightness_ratio")),
                "displacement_los_m": _num(r.get("displacement_los_m")),
                "frequency_drop_pct": _num(r.get("frequency_drop_pct")),
                "gamma2": _num(r.get("coherence_gamma2")),
                "A": _num(r.get("A")),
                "n_masked": _int(r.get("A")),
                "phase_coherence": _num(r.get("phase_coherence")),
                "phase_rms_rad": _num(r.get("phase_rms_rad")),
                "phase_snr_db": _num(r.get("phase_snr_db")),
                "phase_available": _int(r.get("phase_available")) or 0,
                "coherent_sum_amplitude": _num(r.get("coherent_sum_amplitude")),
                "registration_status": r.get("registration_status"),
            })
    for r in rows:
        r["state"] = state_of(r)
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
    return {SHORT[l]: int(sum(1 for r in rows if r["state"] == l))
            for l in STATE_ORDER}


def feature_dataset(rows, feats=None):
    """X, y of the complete cases of ``feats`` (default ``FEATURES``), in order.

    Verbatim ``emv.build_dataset``: only complete-case rows, in committed order,
    feature order as given.
    """
    feats = list(feats if feats is not None else FEATURES)
    complete = [r for r in rows if all(r.get(f) is not None for f in feats)]
    X = np.array([[float(r[f]) for f in feats] for r in complete])
    y = np.array([r["state"] for r in complete], dtype=object)
    return X, y, complete


def row_brief(r):
    """The handful of fields the month companion plots or labels."""
    return {"id": r["id"], "date": r["date"], "month": r["month"],
            "day": r["day"], "doy": r["doy"], "orbit": r["orbit"],
            "state": r["state"], "pass": r["pass"], "gamma2": r["gamma2"],
            "A": r["A"], "phase_coherence": r["phase_coherence"],
            "phase_rms_rad": r["phase_rms_rad"],
            "phase_available": r["phase_available"],
            "wind_speed_ms": r["wind_speed_ms"],
            "temperature_c": r["temperature_c"]}


# ---------------------------------------------------------------------------
# The site
# ---------------------------------------------------------------------------
# The asset the record belongs to, verbatim from the committed upstream report
# (``espoo_mast_observability.json``: ``provenance`` and ``site_context``). The
# package splits the record into the two **orbit geometries** rather than into
# states, so there is no event date here: ``EVENT`` is ``None`` on purpose, and
# every figure says so.
SITE = "Espoo Kurttila mast"
SITE_SHORT = "Espoo"
ASSET_ID = "c4ac1eb8-63e6-481d-997e-4bfe90405d69"
ASSET_NAME = "Espoo Kurttila Mast (OSM 2233589804)"
OSM_NODE = 2233589804
EVENT = None
STATES = STATE_ORDER
# The site's own definition of the two groups, kept explicit so the figure
# script never has to guess what the "state" of a row means.
GROUP_LABEL = {
    "ASCENDING": "ASCENDING — ascending orbit, afternoon pass",
    "DESCENDING": "DESCENDING — descending orbit, morning pass",
}


def brief():
    """Read-only brief of the committed layer (``python3 code/espoo/espoo_ocv_core.py``)."""
    rows = load_csv()
    ref = load_reference()
    meta = load_csv_meta()
    print(f"espoo_ocv_channels.csv  : {len(rows)} acquisitions, "
          f"{rows[0]['date']} .. {rows[-1]['date']}")
    print(f"  derived columns       : {', '.join(meta['new_columns'])}")
    print(f"  by orbit              : {n_by_state(rows)}")
    for k, v in sorted(meta["derived_checks"].items()):
        print(f"  check {k:<28}: {v}")
    obs = ref["observability"]
    print(f"orbit test (committed)  : welch_t {obs['orbit_test']['welch_t']:.4f}, "
          f"p {obs['orbit_test']['welch_p']:.3e}")
    print(f"verdict                 : {obs['verdict']['classes']}")
    sta = ref["stationarity"]
    n_stat, n_grp = st.stationarity_summary(sta)
    print(f"stationarity            : {n_stat} of {n_grp} committed series stationary")
    print(f"primary series          : {sta['primary']}")
    print(f"phase ladder present    : {sum(1 for r in rows if r['phase_available'])}"
          f"/{len(rows)} acquisitions")
    return 0


if __name__ == "__main__":
    raise SystemExit(brief())
