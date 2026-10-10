#!/usr/bin/env python3
"""OCV-Paper — figures A–E, JSON result and English report (Espoo).

Question: the Espoo Kurttila mast record is **event-free** — 150 acquisitions
between 2024-09-05 and 2026-09-07, with no damage of any kind in the window and
no ``damage_label`` / ``structural_state`` / ``condition_label`` on any row. The
committed upstream verdict says so itself: *"This site has no ground-truth state
label"*. So there is no severity axis to regress against, and the only contrast
the record supports is the one between its two **orbit geometries**:

    ASCENDING   ascending orbit, flown in the afternoon pass (n = 70)
    DESCENDING  descending orbit, flown in the morning pass (n = 80)

The site therefore plays the role of the LUMO/Morandi packages' **event-free
control**: the LUMO package asks whether an OCV separates four *damage* states,
Morandi asks the same across a collapse, and this package asks what the same
machinery reports when there is demonstrably nothing to find. Everything the
observability core vector ``OCV = [gamma2, A]`` of this site measures — γ² and
``A`` = ``coherence_masked_pixels`` — is a property of the *geometry*, and the
question is how large that geometry effect is against the record's own
no-event controls.

Approach: the committed table the site's own generator wrote
(``data/espoo/espoo_ocv_channels.csv``, a verified copy of the committed
``espoo_channels.csv`` plus the derived bookkeeping columns), the same
definitions, the same seeds and the same procedures as ``espoo_ocv_stats.py`` —
a verbatim port of the upstream ``espoo_mast_observability.py`` and
``espoo_phase_stationarity.py``. Every number is recomputed here and then pinned
field by field against the three committed reference JSONs in
``data/espoo/reference/``: ``espoo_mast_observability.json``,
``espoo_phase_stationarity.json`` and ``fig_espoo_channels.json``. Any deviation
aborts the script with exit code != 0. On top of that, ``--write-reference``
defines this package's own derivation lock ``espoo_ocv_findings.json`` once; every
later run re-derives it and pins it field by field too. The claim of this repo is
therefore not "we assert", but **"we reproduce and read off"**.

Figures:
  A  raw distributions of the two measured channels per orbit (n = 150)
  B  in-sample effect sizes ASC vs DESC: Cliff's delta + bootstrap CI, SDS
  C  2-class LDA (ASC vs DESC), leave-one-out + label permutation, chance = 0.5
  D  the headline pair across the four feature sets + the out-of-time split
  E  the no-event controls: E1 stationarity, E2 month/season, E3 wind and
     temperature, E4 phase-ladder availability

Usage:
  python3 code/espoo/fig_espoo_ocv_paper.py                   # recompute everything
  python3 code/espoo/fig_espoo_ocv_paper.py --quick           # no permutation nulls
  python3 code/espoo/fig_espoo_ocv_paper.py --write-reference # (once) define the lock
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import math
import os
import platform
import sys
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

try:                                    # scipy is optional in espoo_ocv_stats
    import scipy  # noqa: F401
except ImportError:                     # pragma: no cover
    pass

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import espoo_ocv_core as core   # noqa: E402
import espoo_ocv_stats as st    # noqa: E402

FIGDIR = core.FIGDIR
OUT_JSON = os.path.join(core.DATA, "fig_espoo_ocv_paper.json")
OUT_MD = os.path.join(FIGDIR, "fig_espoo_ocv_paper.md")
REF_FINDINGS = os.path.join(core.REF, "espoo_ocv_findings.json")

# The two states of this package are the orbit geometries, so the contrast runs
# ASC (reference) -> DESC (group), exactly as the committed
# ``espoo_mast_observability.py`` orders its Welch test.
REF_STATE = core.STATE_ORDER[0]     # ASCENDING
GRP_STATE = core.STATE_ORDER[1]     # DESCENDING

# The repository-wide figure style (identical to the LUMO/Morandi/CTS/YWF paper
# scripts and to this site's months companion): DejaVu Sans, 183 mm column width,
# 600 dpi PNG, no top/right spines.
plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "mathtext.fontset": "stix",
    "font.size": 8.5,
    "axes.titlesize": 9.0,
    "axes.labelsize": 8.5,
    "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5,
    "legend.fontsize": 7.5,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.alpha": 0.25,
    "grid.linewidth": 0.4,
    "figure.dpi": 110,
    "savefig.dpi": 600,
    "savefig.bbox": "tight",
    "pdf.fonttype": 42,
})

# 183 mm two-column width (the whole repository uses this one figure size).
MM = 1.0 / 25.4
COL_WIDTH = 183 * MM

STATE_COLORS = {"ASCENDING": "#4C72B0", "DESCENDING": "#DD8452"}
STATE_LABEL = {"ASCENDING": "ASC (afternoon pass)",
               "DESCENDING": "DESC (morning pass)"}
CHANNEL_COLORS = {"gamma2": "#4C72B0", "A": "#8172B3",
                  "phase_coherence": "#55A868", "phase_rms_rad": "#937860",
                  "wind_speed_ms": "#64B5CD", "temperature_c": "#C44E52"}
CHANNEL_TITLE = {"gamma2": "coherence  gamma^2", "A": "mask area  A  [px]",
                 "phase_coherence": "phase coherence",
                 "phase_rms_rad": "phase RMS  [rad]",
                 "wind_speed_ms": "wind  [m/s]",
                 "temperature_c": "temperature  [degC]"}
CHANNEL_AXIS = {"gamma2": "gamma^2  [1]", "A": "A = echo pixels  [px]",
                "phase_coherence": "phase coherence  [1]",
                "phase_rms_rad": "phase RMS  [rad]",
                "wind_speed_ms": "wind speed  [m/s]",
                "temperature_c": "air temperature  [degC]"}

# The four feature sets panels C and D compare (the site's own model registry).
FEATSET_ORDER = ["gamma2", "A", "gamma2A", "x_cov"]
SHORT_MODEL = {"gamma2": "gamma2", "A": "A", "gamma2A": "gamma2+A",
               "x_cov": "gamma2+A+wx"}
MODEL_COLORS = {"gamma2": "#8C8C8C", "A": "#8172B3", "gamma2A": "#4C72B0",
                "x_cov": "#55A868"}
MODEL_MARKERS = {"gamma2": "o", "A": "s", "gamma2A": "^", "x_cov": "D"}

# 2-class chance level — LUMO/Morandi use 1/4 for their four states, this site 1/2.
CHANCE = 0.5

# The channels panel B puts an in-sample effect size on. The two phase-ladder
# columns are reported too, but only over the acquisitions that carry them.
EFFECT_CHANNELS = ["gamma2", "A", "phase_coherence", "phase_rms_rad"]
CONTROL_CHANNELS = ["wind_speed_ms", "temperature_c"]

# Pin comparison rules: exact equality for everything deterministic (delta, SDS,
# medians, counters, CI endpoints, sha256); 1e-6 for numbers that come out of a
# scipy test or a random stream. Deviations above the tolerance => exit != 0.
FLOAT_TOL = 1e-6
TOL_KEYS = {"welch_t", "welch_p", "mannwhitney_p", "mannwhitney_u", "p", "rho",
            "tau", "z", "S", "Q", "r2", "amplitude", "peak_day_of_year",
            "residual_std", "theil_sen_per_year", "lag1_autocorr",
            "lag1_detrended", "p_accuracy", "p_balanced_accuracy",
            "p_decreasing", "p_increasing", "p_two_sided",
            "accuracy_null_mean", "balanced_accuracy_null_mean",
            "accuracy_null_p95", "balanced_accuracy_null_p95",
            "perm_p", "perm_null_mean_sds", "perm_null_p95_sds",
            "diff_mean", "ci95", "slope_ci95_per_year"}


# ---------------------------------------------------------------------------
# Small shared helpers
# ---------------------------------------------------------------------------
def save_fig(fig, name, out_dir=None):
    """Write one 600 dpi PNG (the binaries stay out of git)."""
    out_dir = out_dir or FIGDIR
    png = os.path.join(out_dir, name + ".png")
    fig.savefig(png)
    plt.close(fig)
    print(f"  {os.path.relpath(png)}")
    return png


def panel(ax, letter, text):
    ax.set_title(f"({letter}) {text}", loc="left", fontsize=8.0)


def sha256_file(path):
    return core.sha256(path)


def rel(path):
    """Path relative to the repository root, for the JSON/markdown artifacts."""
    return os.path.relpath(os.path.abspath(path),
                           os.path.abspath(os.path.join(HERE, os.pardir, os.pardir)))


def json_safe(o):
    """Non-finite floats become ``None`` so the artifacts stay standard JSON."""
    if isinstance(o, float):
        return o if math.isfinite(o) else None
    if isinstance(o, np.floating):
        return float(o) if math.isfinite(float(o)) else None
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, dict):
        return {k: json_safe(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [json_safe(v) for v in o]
    return o


def json_default(o):
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.floating):
        return float(o)
    if isinstance(o, np.bool_):
        return bool(o)
    return str(o)


def _jitter(n, width, seed):
    rng = np.random.default_rng(seed)
    return rng.uniform(-width, width, n)


def _vals(rows, key):
    """The defined values of one channel over a row list."""
    return [r[key] for r in rows if r.get(key) is not None]


def _rows(rows, state):
    return [r for r in rows if r["state"] == state]


def _median(vals):
    """The generator's median (``gp._median``), reproduced bit for bit.

    ``espoo_ocv_channels_csv._median`` is a verbatim port of the upstream helper
    and sorts its input before picking the central value. For even n it averages
    the two central values, i.e. it agrees with ``np.median`` on every block of
    ``fig_espoo_channels.json`` (checked, 0 mismatch over the 9 pooled/per-orbit
    blocks: 150/70/80 rows and 80/39/41 phase rows). It is kept as the replica
    rather than replaced by ``np.median`` so that the pin compares the committed
    file against the *same rule that wrote it*, not against a lookalike.
    """
    vals = sorted(v for v in vals if v is not None)
    if not vals:
        return None
    n = len(vals)
    return vals[n // 2] if n % 2 else 0.5 * (vals[n // 2 - 1] + vals[n // 2])


# The five channel blocks of the committed ``fig_espoo_channels.json``: the
# committed key stem and the alias this package uses for the same column.
CHANNEL_REF_KEYS = [("coherence_gamma2", "gamma2"),
                    ("coherence_masked_pixels", "A"),
                    ("phase_coherence", "phase_coherence"),
                    ("phase_rms_rad", "phase_rms_rad"),
                    ("phase_snr_db", "phase_snr_db")]


def channels_block(rows):
    """The committed channel statistics, recomputed from the CSV.

    Verbatim ``espoo_ocv_channels_csv.channel_stats``: pooled and per-orbit
    ``{n, median}`` per channel, plus the monthly composition, the date range
    and the LUMO reference scale the verdict is read against.
    """
    def block(key, sel=None):
        vals = _vals([r for r in rows if sel is None or r["orbit"] == sel], key)
        return {"n": len(vals), "median": _median(vals)}

    out = {"n_acquisitions": len(rows)}
    for committed, alias in CHANNEL_REF_KEYS:
        out[alias] = block(alias)
        for orbit, suffix in (("ASCENDING", "_ascending"),
                              ("DESCENDING", "_descending")):
            out[alias + suffix] = block(alias, orbit)
    out["monthly"] = {m: {o: len([r for r in rows if r["month"] == m
                                  and r["orbit"] == o])
                          for o in core.STATE_ORDER}
                      for m in sorted({r["month"] for r in rows})}
    out["date_first"] = min(r["date"] for r in rows)
    out["date_last"] = max(r["date"] for r in rows)
    out["lumo_reference"] = {"ascending_gamma2": st.LUMO_ASC_G2_MEDIAN,
                             "descending_gamma2": st.LUMO_DESC_G2_MEDIAN}
    return out


# ---------------------------------------------------------------------------
# Analysis blocks (all of them built on espoo_ocv_stats)
# ---------------------------------------------------------------------------
def delta_pair(rows, key, ref=REF_STATE, grp=GRP_STATE,
               n_boot=core.N_BOOT, seed=core.RNG_SEED):
    """``st.delta_block`` of one channel between the two orbit geometries."""
    a = _vals(_rows(rows, ref), key)
    b = _vals(_rows(rows, grp), key)
    return st.delta_block(a, b, f"{core.SHORT[ref]}_vs_{core.SHORT[grp]}",
                          n_boot=n_boot, seed=seed)


def lda_block(rows, features, shrinkage=core.HEADLINE_LAMBDA, n_perm=0,
              seed=core.RNG_SEED, lambda_grid=None):
    """2-class LDA (ASC vs DESC): LOO accuracy, optional lambda grid and
    label-permutation null over one feature set."""
    X, y, complete = core.feature_dataset(rows, features)
    classes = list(core.STATE_ORDER)
    out = {"n": int(X.shape[0]), "n_features": len(features),
           "features": list(features), "classes": classes,
           "shrinkage": shrinkage, "chance_balanced_accuracy": CHANCE}
    for lam in (lambda_grid if lambda_grid is not None else [shrinkage]):
        _, acc, bal, conf = st.loo_cv(X, y, classes, lam)
        out.setdefault("sensitivity_grid", {})[str(lam)] = {
            "accuracy": acc, "balanced_accuracy": bal}
    _, acc, bal, conf = st.loo_cv(X, y, classes, shrinkage)
    out["accuracy"] = acc
    out["balanced_accuracy"] = bal
    out["confusion"] = conf
    if n_perm:
        out["permutation"] = st.permutation_test(X, y, classes, shrinkage,
                                                 n_perm=n_perm, seed=seed)
    return out


def oof_scores_for(rows, features, shrinkage=core.HEADLINE_LAMBDA,
                   class_a=REF_STATE, class_b=GRP_STATE):
    """``(out-of-fold Fisher scores, labels)`` for one feature set."""
    X, y, _complete = core.feature_dataset(rows, features)
    return st.oof_scores(X, y, class_a, class_b, shrinkage), y


def oof_block(rows, features, shrinkage=core.HEADLINE_LAMBDA, n_perm=0,
              n_boot=core.N_BOOT, seed=core.RNG_SEED, q=None, y=None,
              baseline_scores=None, baseline=None,
              class_a=REF_STATE, class_b=GRP_STATE):
    """Out-of-fold Fisher scores + SDS for one feature set, plus the label
    permutation null of that SDS and the paired bootstrap against a baseline."""
    if q is None:
        q, y = oof_scores_for(rows, features, shrinkage, class_a, class_b)
    d, sds, ref_q, grp_q = st.delta_sds(q, y, class_a, class_b)
    out = {"n": int(len(q)), "n_features": len(features),
           "features": list(features), "comparison": f"{class_a}_vs_{class_b}",
           "delta": d, "sds": sds, "shrinkage": shrinkage,
           "score_median": ([float(np.median(ref_q)), float(np.median(grp_q))]
                            if ref_q and grp_q else None)}
    rng = np.random.default_rng(seed)
    if n_perm:
        p, null_mean, null_p95 = st.permutation_p(
            core.feature_dataset(rows, features)[0], y, class_a, class_b,
            shrinkage, sds, n_perm, rng)
        out["permutation"] = {"n_perm": n_perm, "perm_p": p,
                              "perm_null_mean_sds": null_mean,
                              "perm_null_p95_sds": null_p95, "seed": seed}
    if baseline_scores is not None:
        out["paired_vs_baseline"] = dict(
            baseline=baseline,
            **st.paired_bootstrap_diff(q, baseline_scores, y, class_a, class_b,
                                       n_boot, rng))
    return out


def oot_block(rows, features, shrinkage=core.HEADLINE_LAMBDA,
              class_a=REF_STATE, class_b=GRP_STATE, frac=core.OOT_FRACTION):
    """Out-of-**time** split of panel D.

    The complete cases are cut in half by committed acquisition order and each
    half is scored with the model fitted on the other one — an honest
    generalisation check that shares no row between fit and score.
    """
    X, y, complete = core.feature_dataset(rows, features)
    n = X.shape[0]
    k = int(round(n * frac))
    out = {"n": n, "n_first": k, "n_second": n - k,
           "split": [(complete[0]["date"] if k else None),
                     (complete[k]["date"] if k < n else None)]}
    q = np.empty(n)
    q[:] = np.nan
    for fit, score in ((slice(0, k), slice(k, n)), (slice(k, n), slice(0, k))):
        mean, std, w = st.fit_linear_score(X[fit], y[fit], class_a, class_b,
                                           shrinkage)
        for i in range(score.start, score.stop):
            q[i] = st.score_row(X[i], mean, std, w)
    d, sds, ref_q, grp_q = st.delta_sds(q, y, class_a, class_b)
    out.update({"delta": d, "sds": sds, "shrinkage": shrinkage,
                "comparison": f"{class_a}_vs_{class_b}",
                "score_median": ([float(np.median(ref_q)), float(np.median(grp_q))]
                                 if ref_q and grp_q else None)})
    return out


# ---------------------------------------------------------------------------
# The four no-event controls of panel E
# ---------------------------------------------------------------------------
def strata_block(rows, key_fn, n_boot=core.N_BOOT_STRATA):
    """Per-stratum composition + the ASC->DESC effect size of gamma2/A.

    ``espoo_ocv_stats.strata`` is the LUMO version and references a ``healthy``
    state this site does not have, so the same rule is applied to the two orbit
    geometries here.
    """
    out = {}
    for skey in sorted({key_fn(r) for r in rows if key_fn(r) is not None}):
        sub = [r for r in rows if key_fn(r) == skey]
        counts = {core.SHORT[l]: len(_rows(sub, l)) for l in core.STATE_ORDER}
        entry = {"n": len(sub), "n_by_state": counts,
                 "median_gamma2_by_state": {
                     core.SHORT[l]: (float(np.median(_vals(_rows(sub, l), "gamma2")))
                                     if _vals(_rows(sub, l), "gamma2") else None)
                     for l in core.STATE_ORDER}}
        for key in ("gamma2", "A"):
            entry[key] = (delta_pair(sub, key, n_boot=n_boot)
                          if min(counts.values()) >= 3 else None)
        out[str(skey)] = entry
    return out


def _spearman(rows, key, cov):
    """``(n, rho, p)`` of one channel against one co-variate."""
    pairs = [(r[cov], r[key]) for r in rows
             if r.get(cov) is not None and r.get(key) is not None]
    if len(pairs) < 3 or st.sps is None:
        return {"n": len(pairs), "rho": None, "p": None}
    x = np.asarray([a for a, _b in pairs], float)
    y = np.asarray([b for _a, b in pairs], float)
    rho, p = st.sps.spearmanr(x, y)
    return {"n": int(x.size), "rho": float(rho), "p": float(p)}


def covariate_block(rows):
    """E3 — the wind/temperature controls.

    The weather columns are constant over long stretches (they are the upstream
    campaign defaults), so they are reported as *strata* as well as as
    correlations: the geometry contrast is recomputed **inside** each distinct
    wind and temperature value. Only a difference that survives that can be
    blamed on the geometry rather than on the weather.
    """
    out = {"spearman": {}, "by_wind": {}, "by_temperature": {}, "matched_cells": {}}
    for key in ("gamma2", "A"):
        out["spearman"][key] = {cov: _spearman(rows, key, cov)
                                for cov in CONTROL_CHANNELS}
    for cov, slot in (("wind_speed_ms", "by_wind"),
                      ("temperature_c", "by_temperature")):
        table = {}
        for value in sorted({r[cov] for r in rows if r.get(cov) is not None}):
            sub = [r for r in rows if r.get(cov) == value]
            counts = {core.SHORT[l]: len(_rows(sub, l)) for l in core.STATE_ORDER}
            entry = {"value": value, "n": len(sub), "n_by_state": counts}
            for key in ("gamma2", "A"):
                entry[key] = (delta_pair(sub, key, n_boot=core.N_BOOT_STRATA)
                              if min(counts.values()) >= 3 else None)
            table[str(value)] = entry
        out[slot] = table
    cells = {}
    for r in rows:
        if r.get("wind_speed_ms") is None or r.get("temperature_c") is None:
            continue
        cells.setdefault((r["wind_speed_ms"], r["temperature_c"]), []).append(r)
    for (w, t), sub in sorted(cells.items()):
        counts = {core.SHORT[l]: len(_rows(sub, l)) for l in core.STATE_ORDER}
        out["matched_cells"][f"{w:g} m/s / {t:g} degC"] = {
            "n": len(sub), "n_by_state": counts,
            "gamma2": (delta_pair(sub, "gamma2", n_boot=core.N_BOOT_STRATA)
                       if min(counts.values()) >= 3 else None)}
    return out


def phase_availability_block(rows):
    """E4 — the phase ladder: how much of the record carries phase at all."""
    avail = [r for r in rows if r["phase_available"]]
    out = {
        "n_total": len(rows),
        "n_available": len(avail),
        "share": len(avail) / len(rows) if rows else None,
        "n_by_state": {core.SHORT[l]: len(_rows(avail, l))
                       for l in core.STATE_ORDER},
        "n_by_state_total": core.n_by_state(rows),
        "n_snr_db_available": sum(1 for r in rows
                                  if r.get("phase_snr_db") is not None),
        "by_month": {},
        "channels": {k: {core.SHORT[l]: st._stats(_vals(_rows(rows, l), k))
                         for l in core.STATE_ORDER} for k in core.PHASE_DIMS},
    }
    for month in sorted({r["month"] for r in rows}):
        sub = [r for r in rows if r["month"] == month]
        got = [r for r in sub if r["phase_available"]]
        out["by_month"][month] = {
            "n": len(sub), "n_phase": len(got),
            "share": len(got) / len(sub) if sub else None,
            "n_by_state": {core.SHORT[l]: len(_rows(got, l))
                           for l in core.STATE_ORDER}}
    return out


def stationarity_block(sta):
    """E1 — the stationarity layer of ``espoo_phase_stationarity.json``.

    Every one of the 14 committed series is recomputed from its own committed
    ``series_by_group`` entries; the recomputed blocks are pinned field by field
    against the committed ``groups`` blocks, and the file-level claims of
    ``interpretation`` are pinned against the recomputed summary.
    """
    groups, meta = st.phase_series_block(sta)
    n_stat, n_grp = st.stationarity_summary(sta)
    by_orbit = {l: [b["median"] for k, b in groups.items() if k.endswith(l)]
                for l in core.STATE_ORDER}
    # the file-level date range, derived from the committed series themselves
    # rather than copied from the file that states it
    days = [d for series in (sta["series_by_group"] or {}).values()
            for d, _v in series]
    date_range = ([dt.date.fromordinal(int(min(days))).isoformat(),
                   dt.date.fromordinal(int(max(days))).isoformat()]
                  if days else None)
    return {
        "method": sta["method"],
        "metric": sta["metric"],
        "params": sta["params"],
        "n_groups": meta["n_groups"],
        "n_rows": meta["n_rows"],
        "n_by_orbit": meta["n_by_orbit"],
        "date_range": date_range,
        "n_stationary": n_stat,
        "n_series": n_grp,
        "flagged": st.flagged_groups(sta),
        "primary": sta["primary"],
        "primary_block": groups[sta["primary"]],
        "primary_series": sta["series_by_group"][sta["primary"]],
        "median_of_group_medians": {l: float(np.median(v))
                                    for l, v in by_orbit.items() if v},
        "per_group": {k: {"n": b["n"], "median": b["median"],
                          "theil_sen_per_year": b["theil_sen_per_year"],
                          "mann_kendall_p": b["mann_kendall"]["p"],
                          "split_half_delta": (b["split_half"] or {}).get("delta"),
                          "seasonal_r2": (b["seasonal_harmonic"] or {}).get("r2"),
                          "flags": b["flags"], "verdict": b["verdict"]}
                      for k, b in sorted(groups.items())},
        "groups": groups,
    }


# ---------------------------------------------------------------------------
# Compute — one pass over the committed layer
# ---------------------------------------------------------------------------
def compute(quick=False):
    """Every number the JSON, the report and the four pin targets need."""
    t0 = time.time()
    rows = core.load_csv()
    meta = core.load_csv_meta()
    sta = core.load_reference()["stationarity"]
    n_perm = 0 if quick else core.N_PERM_2CLASS

    X2 = core.feature_dataset(rows, core.FEATURES)[0]
    Xp = core.feature_dataset(rows, core.PHASE_DIMS[:2])[0]

    data = {
        "csv": rel(core.CSV_PATH), "csv_sha256": sha256_file(core.CSV_PATH),
        "source_csv": rel(core.RAW_CSV_PATH),
        "source_csv_sha256": sha256_file(core.RAW_CSV_PATH),
        "csv_meta": rel(core.CSV_META_PATH),
        "csv_meta_sha256": sha256_file(core.CSV_META_PATH),
        "generator": meta["generator"],
        "n_acquisitions": len(rows),
        "n_by_state": core.n_by_state(rows),
        "date_first": min(r["date"] for r in rows),
        "date_last": max(r["date"] for r in rows),
        "n_months": len({r["month"] for r in rows}),
        "new_columns": list(meta["new_columns"]),
        "filled_columns": list(meta["filled_columns"]),
        "derived_checks": dict(meta["derived_checks"]),
        "n_complete_2d": int(X2.shape[0]),
        "n_complete_phase2": int(Xp.shape[0]),
        "n_phase_available": sum(1 for r in rows if r["phase_available"]),
        "phase_share": (sum(1 for r in rows if r["phase_available"]) / len(rows)),
        "n_snr_db_available": sum(1 for r in rows
                                  if r.get("phase_snr_db") is not None),
        "months": sorted({r["month"] for r in rows}),
        "orbit": st.orbit_block(rows),
        "channels": channels_block(rows),
        "lumo_reference": {"healthy_g2_median": st.LUMO_HEALTHY_G2_MEDIAN,
                           "healthy_nmasked_median": st.LUMO_HEALTHY_NMASKED_MEDIAN,
                           "asc_g2_median": st.LUMO_ASC_G2_MEDIAN,
                           "desc_g2_median": st.LUMO_DESC_G2_MEDIAN,
                           "clutter_floor": st.LUMO_CLUTTER_FLOOR},
        "constants": {"coherence_peak_frac": st.COHERENCE_PEAK_FRAC,
                      "coherence_median_mult": st.COHERENCE_MEDIAN_MULT,
                      "coherence_min_masked": st.COHERENCE_MIN_MASKED,
                      "pred_min_coherence": st.PRED_MIN_COHERENCE},
    }

    # --- B: in-sample effect sizes ------------------------------------------
    effects = {k: delta_pair(rows, k) for k in EFFECT_CHANNELS}
    effects.update({k: delta_pair(rows, k, n_boot=core.N_BOOT_STRATA)
                    for k in CONTROL_CHANNELS})
    # the headline pair's out-of-fold SDS is the one panel D draws against
    q_base, _y_base = oof_scores_for(rows, core.MODELS[core.BASELINE])
    oof = {name: oof_block(rows, feats, n_perm=n_perm,
                           baseline_scores=q_base if name != core.BASELINE else None,
                           baseline=core.BASELINE)
           for name, feats in core.MODELS.items()}
    oot = {name: oot_block(rows, feats) for name, feats in core.MODELS.items()}

    # --- C: 2-class LDA over the four feature sets ---------------------------
    models = {}
    for name, feats in core.MODELS.items():
        models[name] = lda_block(rows, feats, n_perm=n_perm,
                                 lambda_grid=core.LAMBDA_GRID)
    headline = models[core.HEADLINE_PAIR]
    secondary = models[core.SECONDARY_PAIR]

    # --- E: the four no-event controls --------------------------------------
    controls = {
        "stationarity": stationarity_block(sta),
        "by_season": strata_block(rows, st.season_of),
        "by_orbit": strata_block(rows, lambda r: r.get("orbit")),
        "by_month_n": {m: {"n": len([r for r in rows if r["month"] == m]),
                           "n_by_state": {
                               core.SHORT[l]: len([r for r in rows
                                                   if r["month"] == m
                                                   and r["state"] == l])
                               for l in core.STATE_ORDER}}
                       for m in sorted({r["month"] for r in rows})},
        "covariates": covariate_block(rows),
        "phase_availability": phase_availability_block(rows),
    }

    return {
        "rows": rows, "meta": meta, "data": data, "effects": effects,
        "models": models, "oof": oof, "oot": oot, "controls": controls,
        "headline_permutation": (headline.get("permutation") if not quick else None),
        "secondary_permutation": (secondary.get("permutation") if not quick else None),
        "quick": quick, "n_perm": n_perm, "n_rows": len(rows),
        "runtime_s": time.time() - t0,
    }


# ---------------------------------------------------------------------------
# Rendering — figures A–E, 183 mm column width, 600 dpi PNG (espoo style)
# ---------------------------------------------------------------------------
def fig_A(res, out_dir=None):
    """Raw distributions of the measured channels per orbit geometry."""
    rows = res["rows"]
    fig, axes = plt.subplots(1, 3, figsize=(COL_WIDTH, 0.34 * COL_WIDTH))

    for i, key in enumerate(("gamma2", "A")):
        ax = axes[i]
        for j, s in enumerate(core.STATE_ORDER):
            vals = _vals(_rows(rows, s), key)
            bp = ax.boxplot([vals], positions=[j], widths=0.5, showfliers=False,
                            patch_artist=True,
                            medianprops=dict(color="white", lw=1.1))
            bp["boxes"][0].set(facecolor=STATE_COLORS[s], alpha=0.65,
                               edgecolor=STATE_COLORS[s])
            ax.plot(j + _jitter(len(vals), 0.15, 11 + i), vals, ".", ms=1.7,
                    alpha=0.35, color=STATE_COLORS[s], markeredgewidth=0)
            ax.text(j, 0.99, f"n={len(vals)}\nmed {np.median(vals):.4g}",
                    transform=ax.get_xaxis_transform(), ha="center", va="top",
                    fontsize=6.2, color=STATE_COLORS[s])
        if key == "A":
            ax.set_yscale("log")
        ax.set_xticks(range(len(core.STATE_ORDER)))
        ax.set_xticklabels([core.SHORT[s] for s in core.STATE_ORDER])
        ax.set_ylabel(CHANNEL_AXIS[key])
        panel(ax, "abc"[i], CHANNEL_TITLE[key])

    # (c) the phase ladder, where it exists at all
    ax = axes[2]
    avail = [r for r in rows if r["phase_available"]]
    for j, s in enumerate(core.STATE_ORDER):
        vals = _vals(_rows(avail, s), "phase_coherence")
        ax.plot(j + _jitter(len(vals), 0.15, 42), vals, ".", ms=2.0, alpha=0.45,
                color=STATE_COLORS[s], markeredgewidth=0)
        if vals:
            ax.plot([j - 0.2, j + 0.2], [np.median(vals)] * 2, "-", lw=1.4,
                    color=STATE_COLORS[s])
            ax.text(j, 0.99, f"n={len(vals)}\nmed {np.median(vals):.3f}",
                    transform=ax.get_xaxis_transform(), ha="center", va="top",
                    fontsize=6.2, color=STATE_COLORS[s])
    ax.set_xticks(range(len(core.STATE_ORDER)))
    ax.set_xticklabels([core.SHORT[s] for s in core.STATE_ORDER])
    ax.set_ylabel(CHANNEL_AXIS["phase_coherence"])
    panel(ax, "c", f"phase ladder, the {len(avail)} of {len(rows)} acquisitions that carry it")

    counts = res["data"]["n_by_state"]
    fig.suptitle("(Espoo) two event-free years, one contrast: the orbit geometry — "
                 f"{res['n_rows']} acquisitions, "
                 f"ASC = {counts.get('ASC', 0)}, DESC = {counts.get('DESC', 0)}",
                 fontsize=9.0, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    return save_fig(fig, "fig_espoo_ocv_paper_A", out_dir)


def fig_B(res, out_dir=None):
    """In-sample effect sizes ASC -> DESC: Cliff's delta + CI, SDS."""
    fig, axes = plt.subplots(1, 2, figsize=(COL_WIDTH, 0.40 * COL_WIDTH),
                             gridspec_kw={"width_ratios": [1.05, 1.0]})

    # (a) one row per channel of the observability core vector (+ controls)
    ax = axes[0]
    order = EFFECT_CHANNELS + CONTROL_CHANNELS
    ys = np.arange(len(order))
    for i, k in enumerate(order):
        b = res["effects"][k]
        if not b or b.get("delta") is None:
            ax.text(0.0, i, "n/a", fontsize=6.2, va="center", ha="center",
                    color="#AA2222")
            continue
        ax.plot([b["delta"]], [i], "o", ms=5.0, color=CHANNEL_COLORS[k])
        if b.get("delta_ci95"):
            lo, hi = b["delta_ci95"]
            ax.plot([lo, hi], [i, i], "-", lw=1.4, color=CHANNEL_COLORS[k])
        excl = b.get("ci_excludes_zero")
        ax.text(1.02, i, "CI!=0" if excl else "CI contains 0", fontsize=6.2,
                va="center", color="#444444" if excl else "#AA2222",
                transform=ax.get_yaxis_transform())
    ax.axvline(0.0, color="#888888", lw=0.9)
    ax.set_yticks(ys)
    ax.set_yticklabels([CHANNEL_TITLE[k] + ("  (control)" if k in CONTROL_CHANNELS else "")
                        for k in order], fontsize=7.0)
    ax.invert_yaxis()
    ax.set_xlabel("Cliff's delta  (DESC - ASC)")
    panel(ax, "a", "pooled, all acquisitions (bootstrap 95 % CI)")

    # (b) the same effect size inside each season — the no-event season control
    ax = axes[1]
    seasons = sorted(res["controls"]["by_season"])
    width = 0.38
    for ki, k in enumerate(("gamma2", "A")):
        for si, sk in enumerate(seasons):
            b = res["controls"]["by_season"][sk].get(k)
            if not b or b.get("delta") is None:
                continue
            x = si + (ki - 0.5) * width
            ax.bar(x, b["delta"], width=width, color=CHANNEL_COLORS[k],
                   alpha=0.80 if ki else 0.55,
                   label=CHANNEL_TITLE[k] if si == 0 else None)
    ax.axhline(0.0, color="#888888", lw=0.9)
    ax.set_xticks(range(len(seasons)))
    ax.set_xticklabels([s.split("_")[0] for s in seasons], rotation=10)
    ax.set_ylabel("Cliff's delta  (DESC - ASC)")
    ax.legend(frameon=False)
    panel(ax, "b", "per season (strata of the same contrast)")

    fig.suptitle("(Espoo) the geometry effect is real and it is the largest thing "
                 "in the record", fontsize=9.0, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    return save_fig(fig, "fig_espoo_ocv_paper_B", out_dir)


def fig_C(res, out_dir=None):
    """2-class LDA (ASC vs DESC): LOO accuracy + label-permutation null."""
    fig, axes = plt.subplots(1, 2, figsize=(COL_WIDTH, 0.40 * COL_WIDTH),
                             gridspec_kw={"width_ratios": [1.25, 1.0]})

    # (a) balanced accuracy per feature set against the 0.5 chance level
    ax = axes[0]
    xs = np.arange(len(FEATSET_ORDER))
    accs = [res["models"][m]["balanced_accuracy"] for m in FEATSET_ORDER]
    ax.bar(xs, accs, width=0.62,
           color=[MODEL_COLORS[m] for m in FEATSET_ORDER], alpha=0.80)
    for i, m in enumerate(FEATSET_ORDER):
        blk = res["models"][m]
        p = (blk.get("permutation") or {}).get("p_balanced_accuracy")
        txt = f"p={p:.3g}" if p is not None else "no null (--quick)"
        ax.text(i, min(accs[i] + 0.03, 0.99), txt, ha="center", va="bottom",
                fontsize=6.2, color="#444444")
        ax.text(i, 0.02, f"n={blk['n']}", ha="center", va="bottom", fontsize=6.0,
                color="#FFFFFF")
    ax.axhline(CHANCE, color="#888888", ls="--", lw=0.9)
    ax.text(len(FEATSET_ORDER) - 0.45, CHANCE + 0.015, "chance 0.5", fontsize=6.2,
            color="#666666", ha="right")
    ax.set_xticks(xs)
    ax.set_xticklabels([SHORT_MODEL[m] for m in FEATSET_ORDER], rotation=12,
                       ha="right")
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("LOO balanced accuracy")
    panel(ax, "a", "LDA ASC vs DESC, leave-one-out")

    # (b) headline confusion matrix
    ax = axes[1]
    classes = list(core.STATE_ORDER)
    conf_d = res["models"][core.HEADLINE_PAIR]["confusion"]
    conf = np.asarray([[conf_d[a][b] for b in classes] for a in classes],
                      dtype=float)
    ax.imshow(conf, cmap="Blues", vmin=0, vmax=max(conf.max(), 1))
    for i in range(conf.shape[0]):
        for j in range(conf.shape[1]):
            ax.text(j, i, f"{int(conf[i, j])}", ha="center", va="center",
                    fontsize=7.0,
                    color="#FFFFFF" if conf[i, j] > 0.6 * conf.max() else "#222222")
    ax.set_xticks(range(len(core.STATE_ORDER)))
    ax.set_xticklabels([core.SHORT[s] for s in core.STATE_ORDER])
    ax.set_yticks(range(len(core.STATE_ORDER)))
    ax.set_yticklabels([core.SHORT[s] for s in core.STATE_ORDER])
    ax.set_xlabel("predicted")
    ax.set_ylabel("true (orbit)")
    ax.grid(False)
    panel(ax, "b", f"headline {SHORT_MODEL[core.HEADLINE_PAIR]}, LOO counts")

    perm = res["headline_permutation"]
    tail = ("" if perm is None else
            f"; label permutation n={perm['n_perm']}, "
            f"p(balanced)={perm['p_balanced_accuracy']:.3g}")
    fig.suptitle("(Espoo) can the orbit geometry be read back from the two measured "
                 "channels?" + tail, fontsize=9.0, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    return save_fig(fig, "fig_espoo_ocv_paper_C", out_dir)


def fig_D(res, out_dir=None):
    """The headline pair across the four feature sets, incl. the out-of-time split."""
    fig, axes = plt.subplots(1, 2, figsize=(COL_WIDTH, 0.40 * COL_WIDTH))

    # (a) out-of-fold SDS per feature set, with the permutation null of the SDS
    ax = axes[0]
    xs = np.arange(len(FEATSET_ORDER))
    for i, m in enumerate(FEATSET_ORDER):
        blk = res["oof"][m]
        ax.plot([i], [blk["sds"]], MODEL_MARKERS[m], ms=5.5,
                color=MODEL_COLORS[m], label=SHORT_MODEL[m])
        blk_p = blk.get("paired_vs_baseline")
        if blk_p and blk_p.get("ci95"):
            lo, hi = blk_p["ci95"]
            ax.plot([i, i], [lo + blk["sds"], hi + blk["sds"]], "-", lw=1.2,
                    color=MODEL_COLORS[m], alpha=0.8)
        null = blk.get("permutation")
        if null:
            ax.plot([i - 0.22, i + 0.22],
                    [null["perm_null_p95_sds"]] * 2, "--", lw=1.0,
                    color="#AA2222")
    ax.axhline(0.0, color="#888888", lw=0.9)
    ax.set_xticks(xs)
    ax.set_xticklabels([SHORT_MODEL[m] for m in FEATSET_ORDER])
    ax.set_ylabel("SDS = 10 |delta|  of the out-of-fold score")
    ax.plot([], [], "--", lw=1.0, color="#AA2222", label="permutation null p95")
    ax.legend(frameon=False, ncol=1)
    panel(ax, "a", "out-of-fold SDS (bars = paired CI vs gamma2)")

    # (b) in-sample vs out-of-time
    ax = axes[1]
    width = 0.38
    for ki, (src, lab) in enumerate((("oof", "out-of-fold (LOO)"),
                                     ("oot", "out-of-time (first/second half)"))):
        vals = [res[src][m]["sds"] for m in FEATSET_ORDER]
        ax.bar(xs + (ki - 0.5) * width, vals, width=width,
               color=["#4C72B0", "#DD8452"][ki], alpha=0.85, label=lab)
    ax.axhline(0.0, color="#888888", lw=0.9)
    ax.set_xticks(xs)
    ax.set_xticklabels([SHORT_MODEL[m] for m in FEATSET_ORDER], rotation=12,
                       ha="right")
    ax.set_ylabel("SDS = 10 |delta|")
    ax.legend(frameon=False)
    panel(ax, "b", f"split at {res['oot'][core.HEADLINE_PAIR]['split'][1]}")

    fig.suptitle(f"(Espoo) {core.MODEL_LABEL[core.HEADLINE_PAIR]} is not better than "
                 f"{core.MODEL_LABEL[core.BASELINE]} — the second channel adds nothing",
                 fontsize=9.0, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    return save_fig(fig, "fig_espoo_ocv_paper_D", out_dir)


def fig_E(res, out_dir=None):
    """The four no-event controls: stationarity, season, weather, phase."""
    fig, axes = plt.subplots(2, 2, figsize=(COL_WIDTH, 0.72 * COL_WIDTH))

    # (a) E1 — the stationarity layer: per-series median, flagged series marked
    ax = axes[0][0]
    groups = res["controls"]["stationarity"]["groups"]
    keys = sorted(groups)
    for i, k in enumerate(keys):
        b = groups[k]
        s = k.rsplit(" · ", 1)[-1]
        ax.plot([b["median"]], [i], "o", ms=4.2, color=STATE_COLORS[s],
                markerfacecolor="none" if b["flags"] else STATE_COLORS[s],
                markeredgewidth=1.2)
    for i, k in enumerate(keys):
        if groups[k]["flags"]:
            ax.text(0.02, i, f"  {k}  ({', '.join(groups[k]['flags'])})", va="center",
                    fontsize=6.2, color="#AA2222",
                    transform=ax.get_yaxis_transform())
    ax.set_yticks(range(len(keys)))
    ax.set_yticklabels([k.split(" · ")[0] + " · " + k.split(" · ")[1].replace("_", " ")
                        for k in keys], fontsize=6.0)
    ax.invert_yaxis()
    ax.set_xlabel("median of the committed series")
    n_stat = res["controls"]["stationarity"]["n_stationary"]
    n_ser = res["controls"]["stationarity"]["n_series"]
    panel(ax, "a", f"E1 stationarity: {n_stat} of {n_ser} series stationary")

    # (b) E2 — month/season composition, the acquisition-order control
    ax = axes[0][1]
    months = sorted(res["controls"]["by_month_n"])
    width = 0.8 / len(core.STATE_ORDER)
    for si, s in enumerate(core.STATE_ORDER):
        counts = [res["controls"]["by_month_n"][m]["n_by_state"][core.SHORT[s]]
                  for m in months]
        ax.bar(np.arange(len(months)) + (si - 0.5) * width + width / 2, counts,
               width=width, color=STATE_COLORS[s], alpha=0.85,
               label=core.SHORT[s])
    ax.set_xticks(range(len(months)))
    ax.set_xticklabels([m[2:] for m in months], rotation=90, fontsize=5.6)
    ax.set_ylabel("acquisitions")
    ax.legend(frameon=False, ncol=2)
    panel(ax, "b", f"E2 composition over {len(months)} months (no event to cut on)")

    # (c) E3 — the weather control: the same effect size inside each wind value
    ax = axes[1][0]
    by_wind = res["controls"]["covariates"]["by_wind"]
    wvals = sorted(by_wind, key=float)
    ys, labs = [], []
    for v in wvals:
        entry = by_wind[v]
        b = entry.get("gamma2")
        if not b or b.get("delta") is None:
            continue
        ys.append(len(ys))
        labs.append(f"{float(v):g} m/s\n({entry['n']} acq)")
        ax.plot([b["delta"]], [ys[-1]], "o", ms=4.5, color=CHANNEL_COLORS["gamma2"])
        if b.get("delta_ci95"):
            lo, hi = b["delta_ci95"]
            ax.plot([lo, hi], [ys[-1], ys[-1]], "-", lw=1.3,
                    color=CHANNEL_COLORS["gamma2"])
    pooled = res["effects"]["gamma2"]["delta"]
    ax.axvline(pooled, color="#444444", ls="--", lw=0.9)
    if ys:
        # The label rides on the top of the axes: this site has a single wind
        # stratum, so a data-coordinate offset would land outside the axis.
        ax.text(pooled, 0.97, " pooled", fontsize=6.2, color="#444444",
                va="top", transform=ax.get_xaxis_transform())
    ax.axvline(0.0, color="#888888", lw=0.9)
    ax.set_yticks(ys)
    ax.set_yticklabels(labs, fontsize=6.0)
    # Explicit, inverted limits: one stratum must still give a visible band.
    ax.set_ylim(len(ys) - 0.5 if ys else 0.5, -0.5)
    ax.set_xlabel("Cliff's delta of gamma^2  (DESC - ASC)")
    corr = res["controls"]["covariates"]["spearman"]["gamma2"]["wind_speed_ms"]
    # Two lines: on one line this title reaches into panel (d)'s y label.
    panel(ax, "c", f"E3 wind strata\n(pooled rho={corr['rho']:+.3f}, p={corr['p']:.3f})")

    # (d) E4 — phase-ladder availability
    ax = axes[1][1]
    ph = res["controls"]["phase_availability"]
    months = sorted(ph["by_month"])
    shares = [ph["by_month"][m]["share"] for m in months]
    ns = [ph["by_month"][m]["n"] for m in months]
    ax.bar(range(len(months)), shares, width=0.72, color="#8C8C8C", alpha=0.85)
    for i, (m, s) in enumerate(zip(months, shares)):
        ax.text(i, min(s + 0.03, 1.0), f"{ns[i]}", ha="center", va="bottom",
                fontsize=5.6, color="#444444")
    ax.set_xticks(range(len(months)))
    ax.set_xticklabels([m[2:] for m in months], rotation=90, fontsize=5.6)
    ax.set_ylim(0, 1.15)
    ax.set_ylabel("share of acquisitions with phase")
    panel(ax, "d", f"E4 phase ladder: {ph['n_available']} of {ph['n_total']} acquisitions\n"
                   f"({100 * ph['share']:.1f} %), {ph['n_snr_db_available']} with SNR")

    fig.suptitle("(Espoo) the no-event controls — stationarity, season, weather and "
                 "coverage: nothing here is a damage signal", fontsize=9.0, x=0.01,
                 ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.945))
    return save_fig(fig, "fig_espoo_ocv_paper_E", out_dir)


def fig_paper(res, out_dir=None):
    """All five panels, A–E."""
    return [fn(res, out_dir) for fn in (fig_A, fig_B, fig_C, fig_D, fig_E)]


# ---------------------------------------------------------------------------
# The pin block: every committed number is re-derived and checked; a deviation
# above the tolerance makes the script exit non-zero.
# ---------------------------------------------------------------------------
def compare(prefix, got, want, problems, key=None):
    """Field-by-field comparison; ``want`` (committed) dictates the structure.

    Numbers whose key is in ``TOL_KEYS`` are compared with ``FLOAT_TOL`` (they
    come out of scipy or a random stream); everything else must be exactly equal.
    Non-finite values are normalised to ``None`` on both sides first, so a
    committed ``null`` and a recomputed ``nan`` count as agreement.
    """
    tol = FLOAT_TOL if key in TOL_KEYS else 0.0
    got, want = json_safe(got), json_safe(want)
    if isinstance(want, dict):
        if not isinstance(got, dict):
            problems.append(f"{prefix}: got <{type(got).__name__}> != want <dict>")
            return
        for k in want:
            if k not in got:
                problems.append(f"{prefix}.{k}: missing in the recomputed value")
                continue
            compare(f"{prefix}.{k}", got[k], want[k], problems, k)
        for k in got:
            if k not in want:
                problems.append(f"{prefix}.{k}: extra key in the recomputed value")
        return
    if isinstance(want, list):
        if not isinstance(got, (list, tuple)) or len(got) != len(want):
            problems.append(f"{prefix}: got {got!r} != want {want!r}")
            return
        for i, (g, w) in enumerate(zip(got, want)):
            compare(f"{prefix}[{i}]", g, w, problems, key)
        return
    if isinstance(want, float) or isinstance(got, float):
        if got is None or want is None:
            ok = got is None and want is None
        else:
            ok = abs(float(got) - float(want)) <= tol
    else:
        ok = got == want
    if not ok:
        problems.append(f"{prefix}: got {got!r} != want {want!r}")


def _check(problems, name, got, want, tol=0.0):
    if isinstance(want, (int, float)) and not isinstance(want, bool):
        ok = got is not None and abs(float(got) - float(want)) <= tol
    else:
        ok = got == want
    if not ok:
        problems.append(f"{name}: got {got!r} != want {want!r}")
    return ok


# The five fields of the committed stationarity blocks that this package cannot
# re-derive: the series carry a day and a value, not the acquisition metadata the
# upstream scan recorded, so dwell/sub-aperture/Nyquist counts and the dwell
# control are compared for **internal consistency** rather than recomputed. The
# report says so; nothing here pretends otherwise.
STA_NOT_RECOMPUTED = ("dwell_s", "n_sub", "nyquist_hz", "n_observable",
                      "dwell_control")

# What ``dwell_control`` must contain: the short form when the dwell never
# moved, the full confounder block when it did. Checked per series, so a
# truncated or over-stuffed control cannot pass unnoticed.
DWELL_CONTROL_KEYS = {
    False: ("dwell_s", "varies"),
    True: ("dwell_s", "dwell_vs_time", "metric_vs_dwell", "modal_dwell",
           "modal_subset_mk_p", "modal_subset_n", "modal_subset_trend_per_year",
           "n_sub", "residual_trend_mk_p", "residual_trend_per_year", "varies"),
}


def pin_stationarity(res, sta):
    """Section 4 — the committed stationarity layer, series by series."""
    problems = []
    blk = res["controls"]["stationarity"]
    groups, committed = blk["groups"], sta["groups"]

    _check(problems, "sta.n_groups", blk["n_groups"], 14)
    _check(problems, "sta.n_groups", blk["n_groups"], len(committed))
    _check(problems, "sta.n_rows", blk["n_rows"], sta["n_rows"])
    compare("sta.n_by_orbit", blk["n_by_orbit"], sta["n_by_orbit"], problems)
    _check(problems, "sta.date_range", blk["date_range"], sta["date_range"])

    for key in sorted(committed):
        if key not in groups:
            problems.append(f"sta.groups.{key}: missing in the recomputed blocks")
            continue
        want = {k: v for k, v in committed[key].items()
                if k not in STA_NOT_RECOMPUTED}
        compare(f"sta.groups.{key}", groups[key], want, problems)
        # --- the committed-only fields (see STA_NOT_RECOMPUTED) --------------
        # These four are not re-derivable from the committed series, so they are
        # pinned for internal consistency against the block they live in.
        g = committed[key]
        _check(problems, f"sta.{key}.dwell_s.type", isinstance(g["dwell_s"], list),
               True)
        _check(problems, f"sta.{key}.n_sub.type", isinstance(g["n_sub"], list),
               True)
        _check(problems, f"sta.{key}.nyquist_hz.type",
               isinstance(g["nyquist_hz"], list), True)
        _check(problems, f"sta.{key}.n_observable.type",
               isinstance(g["n_observable"], int), True)
        ctl = g["dwell_control"]
        # the control block is exhaustive iff the dwell actually moved
        _check(problems, f"sta.{key}.dwell_control.varies", ctl["varies"],
               len(set(g["dwell_s"])) > 1)
        _check(problems, f"sta.{key}.dwell_control.keys",
               sorted(ctl), sorted(DWELL_CONTROL_KEYS[ctl.get("varies") is True]))
        # the control repeats the dwell list, rounded to four decimals
        _check(problems, f"sta.{key}.dwell_control.dwell_s", ctl["dwell_s"],
               [round(v, 4) for v in g["dwell_s"]])
        if ctl["varies"]:
            # ... and the sub-slot counts of the series it summarises
            _check(problems, f"sta.{key}.dwell_control.n_sub", ctl["n_sub"],
                   g["n_sub"])
            _check(problems, f"sta.{key}.dwell_control.modal_subset_within_n",
                   ctl["modal_subset_n"] <= g["n"], True)

    # --- the file's own summary, and the prose that states it -------------
    _check(problems, "sta.n_stationary", blk["n_stationary"], 13)
    _check(problems, "sta.n_series", blk["n_series"], 14)
    _check(problems, "sta.flagged", blk["flagged"],
           ["8-18 Hz · full_burst · ASCENDING"])
    _check(problems, "sta.interpretation.len", len(sta["interpretation"]), 7)
    _check(problems, "sta.groups.per_series", len(committed), 14)
    # every number the committed prose quotes, against the recomputed block it
    # describes (the two exceptions read the committed file, and are marked)
    prim = blk["primary_block"]
    _check(problems, "prose.0.n_stationary", blk["n_stationary"], 13)
    _check(problems, "prose.0.n_flagged", blk["n_series"] - blk["n_stationary"], 1)
    _check(problems, "prose.0.flagged", blk["flagged"],
           ["8-18 Hz · full_burst · ASCENDING"])
    _check(problems, "prose.0.flag_names",
           committed["8-18 Hz · full_burst · ASCENDING"]["flags"],
           ["trend", "split_shift"])
    _check(problems, "prose.1.n", prim["n"], 69)
    _check(problems, "prose.1.median", prim["median"], 0.2474, 5e-5)
    _check(problems, "prose.1.trend", prim["theil_sen_per_year"], 0.0041, 5e-5)
    _check(problems, "prose.1.ci_lo", prim["slope_ci95_per_year"][0], -0.0287, 5e-5)
    _check(problems, "prose.1.ci_hi", prim["slope_ci95_per_year"][1], 0.0334, 5e-5)
    _check(problems, "prose.1.mk_p", prim["mann_kendall"]["p"], 0.7442, 5e-5)
    _check(problems, "prose.1.split_delta", prim["split_half"]["delta"], 0.082, 1e-3)
    _check(problems, "prose.1.harmonic_amplitude",
           prim["seasonal_harmonic"]["amplitude"], 0.0154, 5e-5)
    _check(problems, "prose.1.harmonic_r2", prim["seasonal_harmonic"]["r2"], 0.037,
           5e-4)
    _check(problems, "prose.1.verdict", prim["verdict"], "stationary")
    flag = groups["8-18 Hz · full_burst · ASCENDING"]
    _check(problems, "prose.2.trend", flag["theil_sen_per_year"], -0.0030, 5e-5)
    _check(problems, "prose.2.mk_p", flag["mann_kendall"]["p"], 0.0081, 5e-5)
    _check(problems, "prose.2.split_delta", flag["split_half"]["delta"], -0.425,
           5e-4)
    fc = committed["8-18 Hz · full_burst · ASCENDING"]["dwell_control"]
    _check(problems, "prose.2.rho_metric_dwell", fc["metric_vs_dwell"]["rho"],
           0.083, 5e-4)
    _check(problems, "prose.2.p_metric_dwell", fc["metric_vs_dwell"]["p"],
           0.4975, 5e-4)
    _check(problems, "prose.2.rho_dwell_time", fc["dwell_vs_time"]["rho"],
           -0.100, 5e-4)
    _check(problems, "prose.2.residual_trend", fc["residual_trend_per_year"],
           -0.0029, 5e-5)
    _check(problems, "prose.2.residual_mk_p", fc["residual_trend_mk_p"],
           0.0184, 5e-5)
    _check(problems, "prose.2.modal_subset_n", fc["modal_subset_n"], 25)
    _check(problems, "prose.2.modal_subset_trend",
           fc["modal_subset_trend_per_year"], -0.0061, 5e-5)
    mm = blk["median_of_group_medians"]
    _check(problems, "prose.3.asc_median", mm["ASCENDING"], 0.1728, 5e-5)
    _check(problems, "prose.3.desc_median", mm["DESCENDING"], 0.1443, 5e-5)
    _check(problems, "prose.3.ratio", mm["ASCENDING"] / mm["DESCENDING"], 1.20, 5e-3)
    n_zero = sum(1 for g in committed.values() if g["n_observable"] == 0)
    _check(problems, "prose.4.n_zero_observable (committed file)", n_zero, 11)
    slot = sta["slot_persistence"]
    _check(problems, "prose.5.slot_rank_corr.ASC",
           slot["ASCENDING"]["adjacent_rank_correlation_median"], 0.018, 5e-4)
    _check(problems, "prose.5.slot_rank_corr.DESC",
           slot["DESCENDING"]["adjacent_rank_correlation_median"], 0.030, 5e-4)
    _check(problems, "prose.5.columns_never_above_gate.ASC (committed file)",
           slot["ASCENDING"]["columns_never_above_gate"], slot["modal_columns"])
    # the two primary series of the table (panel E1 / the month companion)
    _check(problems, "sta.primary.median", prim["median"], 0.2474, 5e-5)
    _check(problems, "sta.primary.desc.median",
           groups["3-8 Hz · strip400_pipeline · DESCENDING"]["median"], 0.2165,
           5e-5)
    return problems


def pin_block(res, findings_ref=None):
    """``res`` (recomputed) against the three committed site references.

    Returns the list of deviations — empty means every committed number of this
    site reproduced field by field.
    """
    problems = []

    def check(name, got, want, tol=0.0):
        return _check(problems, name, got, want, tol)

    data, meta = res["data"], res["meta"]
    ref = core.load_reference()
    obs, sta, chan = ref["observability"], ref["stationarity"], ref["channels"]

    # --- 1) the generated channel table ---------------------------------
    check("data.csv", data["csv"], "data/espoo/espoo_ocv_channels.csv")
    check("data.source_csv", data["source_csv"], "data/espoo/espoo_channels.csv")
    check("data.generator", data["generator"],
          "code/espoo/espoo_ocv_channels_csv.py")
    check("data.n_acquisitions", data["n_acquisitions"], 150)
    check("data.n_acquisitions", data["n_acquisitions"],
          meta["derived_checks"]["rows_read"])
    check("data.n_by_state.ASC", data["n_by_state"]["ASC"], 70)
    check("data.n_by_state.DESC", data["n_by_state"]["DESC"], 80)
    check("data.date_first", data["date_first"], "2024-09-05")
    check("data.date_last", data["date_last"], "2026-09-07")
    check("data.n_months", data["n_months"], 25)
    check("data.new_columns", data["new_columns"], ["A", "doy", "phase_available"])
    check("data.filled_columns", data["filled_columns"],
          ["month", "acquisition_date"])
    # the digests the generator recorded, plus the digest of its metadata file:
    # a tampered table or metadata file cannot masquerade as the committed one
    check("data.csv_sha256", data["csv_sha256"], meta["output_sha256"])
    check("data.source_csv_sha256", data["source_csv_sha256"],
          meta["source_sha256"])
    for k in ("csv_sha256", "source_csv_sha256", "csv_meta_sha256"):
        check(f"data.{k}.length", len(data[k]), 64)
    compare("data.derived_checks", data["derived_checks"], {
        "A_equals_committed_masked_pixels": True,
        "channel_stats_verified": 235,
        "coherence_rows": 150,
        "columns_appended": 3,
        "columns_filled_in_place": 2,
        "copied_fields_verified_against_committed_json": 4050,
        "loaded_by_core": 150,
        "n_months": 25,
        "rows_left_empty": 0,
        "rows_read": 150}, problems)
    # a three-way cross-check of the record's composition — the generator's own
    # metadata, the committed channel reference, and the table itself
    compare("data.meta_monthly_vs_channels_monthly", meta["monthly"],
            chan["monthly"], problems)
    compare("data.monthly_vs_rows", data["channels"]["monthly"],
            meta["monthly"], problems)
    check("data.months", data["months"], sorted(meta["monthly"]))

    # --- 2) the committed observability layer ---------------------------
    compare("obs.per_pass", data["orbit"]["per_pass"], obs["per_pass"], problems)
    compare("obs.per_orbit", data["orbit"]["per_orbit"], obs["per_orbit"], problems)
    compare("obs.orbit_test", data["orbit"]["orbit_test"], obs["orbit_test"],
            problems)
    compare("obs.verdict.classes", data["orbit"]["classes"],
            obs["verdict"]["classes"], problems)
    compare("obs.constants", data["constants"], obs["constants"], problems)
    compare("obs.lumo_reference", data["lumo_reference"], obs["lumo_reference"],
            problems)
    # the request provenance: 25 monthly cadence slots asked for, 371 candidate
    # scenes in the catalogue, 150 acquisitions kept by the pipeline
    prov, cat = obs["provenance"], obs["catalog"]
    check("obs.provenance.total_acquisitions (cadence slots)",
          prov["total_acquisitions"], len(cat["by_month"]))
    check("obs.provenance.asset_id", prov["asset_id"], core.ASSET_ID)
    check("obs.provenance.asset_name", prov["asset_name"], core.ASSET_NAME)
    check("obs.provenance.time_start", prov["time_start"], "2024-09-01")
    check("obs.provenance.time_end", prov["time_end"], "2026-09-14")
    check("obs.provenance.cadence", prov["cadence"], "monthly")
    check("obs.provenance.passes", prov["passes"], ["morning", "afternoon"])
    check("obs.provenance.status", prov["status"], "done")
    check("obs.catalog.total (candidates before selection)", cat["total"], 371)
    check("obs.catalog.by_track_sum", sum(cat["by_track"].values()), cat["total"])
    check("obs.catalog.by_month_sum", sum(cat["by_month"].values()), cat["total"])
    check("obs.catalog.by_month.n", len(cat["by_month"]),
          prov["total_acquisitions"])
    check("obs.catalog.by_track.n", len(cat["by_track"]), 13)
    check("obs.acquisitions.n", len(obs["acquisitions"]), 150)
    check("obs.verdict.notes.n", len(obs["verdict"]["notes"]), 5)
    check("obs.verdict.headline.orbit_stratified",
          "orbit-stratified" in obs["verdict"]["headline"], True)
    check("obs.verdict.notes.no_state_label",
          any("no ground-truth state label" in n.lower()
              for n in obs["verdict"]["notes"]), True)
    # the headline claim of the site, quoted verbatim in the committed notes
    po = data["orbit"]["per_orbit"]
    check("orbit.ASC.gamma2_median", po["ASCENDING"]["gamma2_median"], 0.1746, 5e-5)
    check("orbit.DESC.gamma2_median", po["DESCENDING"]["gamma2_median"], 0.0258, 5e-5)
    check("orbit.ASC.n_masked_median", po["ASCENDING"]["n_masked_median"], 6.0)
    check("orbit.DESC.n_masked_median", po["DESCENDING"]["n_masked_median"], 11.0)
    check("orbit_test.n_a", data["orbit"]["orbit_test"]["n_a"], 70)
    check("orbit_test.n_b", data["orbit"]["orbit_test"]["n_b"], 80)
    check("orbit_test.welch_t", data["orbit"]["orbit_test"]["welch_t"], 7.31, 5e-3)
    check("orbit_test.welch_p", data["orbit"]["orbit_test"]["welch_p"], 2.7e-10,
          2e-11)

    # --- 3) the committed channel reference -----------------------------
    for committed, alias in CHANNEL_REF_KEYS:
        for suffix in ("", "_ascending", "_descending"):
            compare(f"chan.{committed}{suffix}",
                    data["channels"][alias + suffix], chan[committed + suffix],
                    problems)
    check("chan.n_acquisitions", data["channels"]["n_acquisitions"],
          chan["n_acquisitions"])
    check("chan.date_first", data["channels"]["date_first"], chan["date_first"])
    check("chan.date_last", data["channels"]["date_last"], chan["date_last"])
    compare("chan.lumo_reference", data["channels"]["lumo_reference"],
            chan["lumo_reference"], problems)
    # the LUMO reference scale is not a number of this file — it is the constant
    # the classification layer of espoo_ocv_stats reads, so pin it against that
    check("chan.lumo_reference.ascending_gamma2",
          chan["lumo_reference"]["ascending_gamma2"], st.LUMO_ASC_G2_MEDIAN)
    check("chan.lumo_reference.descending_gamma2",
          chan["lumo_reference"]["descending_gamma2"], st.LUMO_DESC_G2_MEDIAN)
    check("obs.lumo_reference.asc_g2_median",
          obs["lumo_reference"]["asc_g2_median"], st.LUMO_ASC_G2_MEDIAN)
    check("obs.lumo_reference.desc_g2_median",
          obs["lumo_reference"]["desc_g2_median"], st.LUMO_DESC_G2_MEDIAN)

    # --- 4) the committed stationarity layer ----------------------------
    problems.extend(pin_stationarity(res, sta))
    check("sta.primary", res["controls"]["stationarity"]["primary"],
          "3-8 Hz · strip400_pipeline · ASCENDING")

    if findings_ref is not None:
        compare("findings", findings(res), findings_ref, problems)
    return problems


# ---------------------------------------------------------------------------
# The package's own derivation lock, the artifact and the report
# ---------------------------------------------------------------------------
def findings(res):
    """This package's own reference (``data/espoo/reference/espoo_ocv_findings.json``).

    Everything in here follows from the committed tables alone — no runtime, no
    environment, no path that depends on the caller — so the file is a
    *definition* of this package's analysis layer: written once by
    ``--write-reference``, re-derived and pinned field by field on every later
    full run. The three committed site references above stay the upstream pins.
    """
    data = res["data"]
    return {
        "site": core.SITE,
        "site_short": core.SITE_SHORT,
        "generated_by": "code/espoo/fig_espoo_ocv_paper.py",
        "state_order": list(core.STATE_ORDER),
        "states": dict(core.SHORT),
        "comparison": f"{core.SHORT[REF_STATE]}_vs_{core.SHORT[GRP_STATE]}",
        "seed": core.RNG_SEED,
        "n_boot": core.N_BOOT,
        "n_boot_strata": core.N_BOOT_STRATA,
        "n_perm_2class": core.N_PERM_2CLASS,
        "shrinkage": core.HEADLINE_LAMBDA,
        "lambda_grid": list(core.LAMBDA_GRID),
        "features": list(core.FEATURES),
        "phase_dims": list(core.PHASE_DIMS),
        "models": {k: list(v) for k, v in core.MODELS.items()},
        "headline_pair": core.HEADLINE_PAIR,
        "baseline": core.BASELINE,
        "secondary_pair": core.SECONDARY_PAIR,
        "oot_fraction": core.OOT_FRACTION,
        "chance_balanced_accuracy": CHANCE,
        "table": {
            "csv": data["csv"], "csv_sha256": data["csv_sha256"],
            "source_csv": data["source_csv"],
            "source_csv_sha256": data["source_csv_sha256"],
            "csv_meta": data["csv_meta"],
            "csv_meta_sha256": data["csv_meta_sha256"],
            "generator": data["generator"],
            "n_acquisitions": data["n_acquisitions"],
            "n_by_state": data["n_by_state"],
            "date_first": data["date_first"], "date_last": data["date_last"],
            "n_months": data["n_months"],
            "new_columns": data["new_columns"],
            "filled_columns": data["filled_columns"],
            "derived_checks": data["derived_checks"],
            "n_phase_available": data["n_phase_available"],
            "phase_share": data["phase_share"],
            "n_snr_db_available": data["n_snr_db_available"],
            "n_complete_2d": data["n_complete_2d"],
        },
        "orbit": data["orbit"],
        "channels": data["channels"],
        "lumo_reference": data["lumo_reference"],
        "constants": data["constants"],
        "effects": json_safe(res["effects"]),
        "models_fit": json_safe(res["models"]),
        "oof": json_safe(res["oof"]),
        "oot": json_safe(res["oot"]),
        "controls": json_safe(res["controls"]),
    }


def build_json(res, problems, findings_ref, wrote_reference):
    """The result JSON — every recomputed number, without the per-row rows."""
    data = res["data"]
    headline = res["models"][core.HEADLINE_PAIR]
    return {
        "site": core.SITE,
        "site_short": core.SITE_SHORT,
        "generated_by": "code/espoo/fig_espoo_ocv_paper.py",
        "state_order": list(core.STATE_ORDER),
        "states": dict(core.SHORT),
        "comparison": f"{core.SHORT[REF_STATE]}_vs_{core.SHORT[GRP_STATE]}",
        "event": None,
        "event_note": ("the record is event-free: no damage_label, no structural_state, "
                       "no condition_label on any of the 150 rows and no event in the "
                       "window; the two groups are the orbit geometries"),
        "rng_seed": core.RNG_SEED,
        "n_boot": core.N_BOOT,
        "n_boot_strata": core.N_BOOT_STRATA,
        "n_perm_2class": core.N_PERM_2CLASS,
        "shrinkage": core.HEADLINE_LAMBDA,
        "chance_balanced_accuracy": CHANCE,
        "quick": res["quick"],
        "n_rows": res["n_rows"],
        "state_counts": data["n_by_state"],
        "table": {k: data[k] for k in
                  ("csv", "csv_sha256", "source_csv", "source_csv_sha256",
                   "csv_meta", "csv_meta_sha256", "generator",
                   "n_acquisitions", "date_first", "date_last", "n_months",
                   "new_columns", "filled_columns", "derived_checks",
                   "n_phase_available", "phase_share", "n_snr_db_available",
                   "n_complete_2d", "n_complete_phase2")},
        "orbit": data["orbit"],
        "channels": data["channels"],
        "lumo_reference": data["lumo_reference"],
        "constants": data["constants"],
        "effects": json_safe(res["effects"]),
        "models": json_safe(res["models"]),
        "oof": json_safe(res["oof"]),
        "oot": json_safe(res["oot"]),
        "controls": json_safe(res["controls"]),
        "headline": {
            "model": core.HEADLINE_PAIR,
            "balanced_accuracy": headline["balanced_accuracy"],
            "accuracy": headline["accuracy"],
            "permutation": headline.get("permutation"),
        },
        "pins": {"ok": not problems, "n_problems": len(problems),
                 "problems": problems,
                 "findings_pinned": findings_ref is not None,
                 "findings_written": bool(wrote_reference)},
        "figures": [f"figures/espoo/fig_espoo_ocv_paper_{L}.png"
                    for L in "ABCDE"],
    }


def _fmt(x, nd=4):
    if x is None:
        return "--"
    if isinstance(x, float):
        return f"{x:.{nd}g}"
    return str(x)


def md_table(header, rows):
    out = ["| " + " | ".join(header) + " |",
           "|" + "|".join(["---"] * len(header)) + "|"]
    out += ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]
    return "\n".join(out)


def report(res, problems, findings_ref):
    """The markdown report written next to the PNGs."""
    data = res["data"]
    stt = res["controls"]["stationarity"]
    obs = data["orbit"]
    eff = res["effects"]
    po = obs["per_orbit"]
    md = []
    A = md.append

    A(f"# {core.SITE} — OCV-Paper, figures A–E")
    A("")
    A(f"*{core.SITE_SHORT} · event-free control site · "
      f"{data['n_acquisitions']} acquisitions · {data['date_first']} … "
      f"{data['date_last']} · generated by `code/espoo/fig_espoo_ocv_paper.py`*")
    A("")
    A("## The question")
    A("")
    A("This record has **no ground-truth state label** — no `damage_label`, no "
      "`structural_state`, no `condition_label` on any of its "
      f"{data['n_acquisitions']} rows — and the committed upstream verdict says so "
      "in its own words. There is therefore no severity axis to regress against. "
      "The only contrast the record supports is the one between its two **orbit "
      "geometries**:")
    A("")
    A(f"* **{core.GROUP_LABEL['ASCENDING']}** (n = {data['n_by_state'].get('ASC', 0)})")
    A(f"* **{core.GROUP_LABEL['DESCENDING']}** (n = {data['n_by_state'].get('DESC', 0)})")
    A("")
    A("The observability core vector of the site is "
      f"`OCV = [{', '.join(core.FEATURES)}]` — interferometric coherence gamma^2 "
      "and `A` = `coherence_masked_pixels`, the size of the coherent echo mask. "
      "The question of this figure set is not *whether an event can be seen* "
      "(there is none declared and none in the window) but **how large the "
      "geometry effect is against the record's own no-event controls** — which is "
      "exactly the baseline any damage claim at another site has to beat.")
    A("")
    A("## What each panel shows")
    A("")
    A("### A — raw distributions")
    A("")
    rows = []
    for k, col in (("gamma2", "gamma2_median"), ("A", "n_masked_median")):
        rows.append([CHANNEL_TITLE[k]]
                    + [_fmt(po[s][col]) for s in core.STATE_ORDER])
    A(md_table(["channel", f"ASC (n={po['ASCENDING']['n']})",
                f"DESC (n={po['DESCENDING']['n']})"], rows))
    A("")
    A(f"The two distributions barely overlap: gamma^2 sits at "
      f"{_fmt(po['ASCENDING']['gamma2_median'])} on ASC against "
      f"{_fmt(po['DESCENDING']['gamma2_median'])} on DESC, while the mask is "
      f"*smaller* on ASC ({_fmt(po['ASCENDING']['n_masked_median'])} px) than on "
      f"DESC ({_fmt(po['DESCENDING']['n_masked_median'])} px). A compact, highly "
      "coherent echo against a wide, incoherent one — a geometry signature, not a "
      "damage signature. Panel (c) adds the phase ladder, which exists for only "
      f"{data['n_phase_available']} of {data['n_acquisitions']} acquisitions "
      f"({100 * data['phase_share']:.1f} %) and for {data['n_snr_db_available']} "
      "of them carries an SNR at all.")
    A("")
    A("### B — in-sample effect sizes")
    A("")
    rows = []
    for k in EFFECT_CHANNELS + CONTROL_CHANNELS:
        b = eff[k]
        rows.append([f"`{k}`" + (" *(control)*" if k in CONTROL_CHANNELS else ""),
                     _fmt((b or {}).get("delta")),
                     (f"[{_fmt(b['delta_ci95'][0])}, {_fmt(b['delta_ci95'][1])}]"
                      if (b or {}).get("delta_ci95") else "--"),
                     _fmt((b or {}).get("sds")),
                     _fmt((b or {}).get("n", [None, None])[0])])
    A(md_table(["channel", "Cliff's delta (DESC − ASC)", "bootstrap 95 % CI",
                "SDS", "n ASC"], rows))
    A("")
    A(f"The orbit test on gamma^2: Welch t = {_fmt(obs['orbit_test']['welch_t'])}, "
      f"p = {obs['orbit_test']['welch_p']:.3g}, Mann-Whitney p = "
      f"{obs['orbit_test']['mannwhitney_p']:.3g}. The geometry contrast is the "
      "largest effect in the record — larger than any season stratum (panel B b) "
      "and larger than either weather co-variate (panel E c).")
    A("")
    A("### C — is the geometry learnable from the two measured channels?")
    A("")
    rows = []
    for m in FEATSET_ORDER:
        blk = res["models"][m]
        perm = blk.get("permutation") or {}
        rows.append([f"`{m}` — {SHORT_MODEL[m]}", _fmt(blk["accuracy"]),
                     _fmt(blk["balanced_accuracy"]),
                     _fmt(perm.get("p_accuracy")),
                     _fmt(perm.get("p_balanced_accuracy")), _fmt(blk["n"])])
    A(md_table(["model", "accuracy", "balanced accuracy", "p (accuracy)",
                "p (balanced)", "n"], rows))
    A("")
    A(f"Chance is **{CHANCE}** (two classes, balanced accuracy). Leave-one-out "
      "over the whole record — the labels are never shuffled, so the number is "
      "honest; the permutation null only quantifies how surprising it is.")
    A("")

    A("### D — the headline pair across the four feature sets, and out of time")
    A("")
    rows = []
    for m in FEATSET_ORDER:
        o = res["oof"][m]
        t = res["oot"][m]
        pb = o.get("paired_vs_baseline") or {}
        perm = o.get("permutation") or {}
        rows.append([f"`{m}` — {SHORT_MODEL[m]}",
                     _fmt(o["delta"]), _fmt(o["sds"]),
                     _fmt(perm.get("perm_p")),
                     (f"[{_fmt(pb['ci95'][0])}, {_fmt(pb['ci95'][1])}]"
                      if pb.get("ci95") else "--"),
                     _fmt(t["sds"])])
    A(md_table(["model", "delta (out-of-fold)", "SDS", "p (permutation of SDS)",
                "paired CI of the SDS difference vs gamma^2", "SDS (out-of-time)"],
               rows))
    A("")
    oot = res["oot"][core.HEADLINE_PAIR]
    A(f"The out-of-time split cuts the record in half by acquisition order "
      f"({oot['n_first']} / {oot['n_second']}, split at {oot['split'][1]}) and "
      "scores each half with the model fitted on the other one. The second "
      "channel and the two weather co-variates add nothing that survives the "
      "split: the honest reading is that the *orbit* is separable and the "
      "*geometry* is all the vector measures.")
    A("")
    A("### E — the four no-event controls")
    A("")
    A("**E1 stationarity.** "
      f"{stt['n_stationary']} of {stt['n_series']} committed phase-coherence "
      "series are stationary by the tests the site itself applies (Theil-Sen "
      "slope with a moving-block bootstrap CI, Mann-Kendall, first/second-half "
      "split, annual harmonic, lag-1 autocorrelation, runs test and Ljung-Box on "
      "the detrended series). The one series carrying a flag is "
      f"`{stt['flagged'][0] if stt['flagged'] else '--'}`"
      + (f" ({', '.join(stt['groups'][stt['flagged'][0]]['flags'])})"
         if stt["flagged"] else "") + ". A stationary series is the *precondition* "
      "for a persistence channel, not a persistence channel.")
    A("")
    rows = []
    for k in sorted(stt["groups"]):
        b = stt["groups"][k]
        rows.append([f"`{k}`", _fmt(b["n"]), _fmt(b["median"]),
                     _fmt(b["theil_sen_per_year"]), _fmt(b["mann_kendall"]["p"]),
                     _fmt((b["split_half"] or {}).get("delta")),
                     ", ".join(b["flags"]) or "—"])
    A(md_table(["series", "n", "median", "trend /yr", "MK p", "split delta",
                "flags"], rows))
    A("")
    A("Five of the committed fields of these blocks — `dwell_s`, `n_sub`, "
      "`nyquist_hz`, `n_observable` and `dwell_control` — are **not** re-derivable "
      "from the committed series (a series is a list of `[day, value]` pairs), so "
      "the pin block compares them for internal consistency and says so; every "
      "other field of every one of the "
      f"{stt['n_series']} series is recomputed and pinned field by field.")
    A("")
    A("**E2 month and season.** There is no event to cut on, so the record is "
      "split by season and by month instead. Every season reproduces the same "
      "sign and roughly the same size, which is what a geometry effect looks "
      "like.")
    A("")
    rows = []
    for sk in sorted(res["controls"]["by_season"]):
        e = res["controls"]["by_season"][sk]
        g = e.get("gamma2") or {}
        rows.append([f"`{sk}`", _fmt(e["n"]), _fmt(g.get("delta")),
                     _fmt(e["A"].get("delta") if e.get("A") else None)])
    A(md_table(["season", "n", "delta gamma^2", "delta A"], rows))
    A("")

    A("**E3 wind and temperature.** The weather columns of this export are "
      "constant over long stretches, so the contrast is not *corrected* for them — "
      "it is recomputed **inside** each distinct value, which is the strongest "
      "form of the control available here.")
    A("")
    rows = []
    for key, cov_key in (("gamma2", "wind_speed_ms"), ("gamma2", "temperature_c"),
                         ("A", "wind_speed_ms"), ("A", "temperature_c")):
        sp = res["controls"]["covariates"]["spearman"][key][cov_key]
        rows.append([f"`{key}` vs `{cov_key}`", _fmt(sp["n"]), _fmt(sp["rho"]),
                     _fmt(sp["p"])])
    A(md_table(["pair", "n", "Spearman rho", "p"], rows))
    A("")
    cov = res["controls"]["covariates"]
    rows = []
    for label, slot in (("wind", "by_wind"), ("temperature", "by_temperature")):
        for v, e in sorted(cov[slot].items(), key=lambda kv: float(kv[0])):
            g = e.get("gamma2") or {}
            if g.get("delta") is None:
                continue
            rows.append([f"{label} = {_fmt(float(v), 6)}", _fmt(e["n"]),
                         f"{e['n_by_state']['ASC']} / {e['n_by_state']['DESC']}",
                         _fmt(g["delta"]),
                         (f"[{_fmt(g['delta_ci95'][0])}, {_fmt(g['delta_ci95'][1])}]"
                          if g.get("delta_ci95") else "--")])
    A(md_table(["stratum", "n", "ASC / DESC", "delta gamma^2", "bootstrap 95 % CI"],
               rows))
    A("")
    A(f"Pooled delta gamma^2 = {_fmt(eff['gamma2']['delta'])}. Every stratum with "
      "enough acquisitions on both sides keeps the same sign; none inverts it, so "
      "weather does not explain the geometry effect.")
    A("")
    A("**E4 phase-ladder availability.**")
    ph = res["controls"]["phase_availability"]
    A("")
    A(f"`phase_coherence` and `phase_rms_rad` exist for {ph['n_available']} of "
      f"{ph['n_total']} acquisitions ({100 * ph['share']:.1f} %), "
      f"`phase_snr_db` for {ph['n_snr_db_available']}. The ladder is therefore "
      "reported as *availability*, never as a third observability channel: a "
      "channel that is missing on half the record cannot carry a claim about the "
      "record.")
    A("")
    rows = []
    for m in sorted(ph["by_month"]):
        e = ph["by_month"][m]
        rows.append([m, _fmt(e["n"]), _fmt(e["n_phase"]),
                     f"{100 * e['share']:.0f} %" if e["share"] is not None else "--",
                     f"{e['n_by_state']['ASC']} / {e['n_by_state']['DESC']}"])
    A(md_table(["month", "acquisitions", "with phase", "share", "ASC / DESC"], rows))
    A("")

    A("## What this control says")
    A("")
    A("At a site with **no event at all**, the observability core vector still "
      "splits its record into two groups that are almost perfectly separable — "
      "and the split is the **orbit geometry**, not a state of the structure. "
      "That is the point of including Espoo in this repository:")
    A("")
    A(f"1. an effect size of {_fmt(obs['orbit_test']['welch_t'])} (Welch t) and a "
      f"leave-one-out balanced accuracy of "
      f"{_fmt(res['models'][core.HEADLINE_PAIR]['balanced_accuracy'])} look like "
      "*damage* under any classifier-first reading;")
    A("2. the record's own no-event controls (stationarity, season, weather, "
      "coverage) show every one of those numbers is geometry;")
    A("3. so at every other site of this repository the orbit-geometry contrast "
      "has to be reported **beside** the damage contrast before the latter can be "
      "read — which is exactly what the LUMO, Morandi, YWF and Carola packages "
      "do.")
    A("")
    A("## Pins")
    A("")
    A(f"Every committed number of this site is re-derived and compared: "
      + ("**all checks passed**." if not problems
         else f"**{len(problems)} check(s) FAILED**."))
    A("")
    A("| target | what is pinned |")
    A("|---|---|")
    A("| `espoo_ocv_channels.csv` | sha256 against the generator's own "
      "`output_sha256`, plus the source digest; `A` copied, not recomputed; "
      "150 rows, 70 ASC / 80 DESC, 25 months, `2024-09-05` … `2026-09-07` |")
    A("| `espoo_mast_observability.json` | `per_pass`, `per_orbit`, `orbit_test`, "
      "`verdict.classes`, `constants`, `lumo_reference`, provenance and catalogue "
      "— field by field |")
    A("| `fig_espoo_channels.json` | the 15 pooled/per-orbit `{n, median}` "
      "blocks, the 25-month composition and the LUMO reference scale |")
    A("| `espoo_phase_stationarity.json` | all 14 series: median, IQR, Theil-Sen "
      "slope + bootstrap CI, Mann-Kendall, lag-1, runs test, Ljung-Box, "
      "half-split, annual harmonic, monthly means, flags and verdict — plus every "
      "number the file's own prose quotes |")
    A("| `espoo_ocv_findings.json` | this package's own derivation lock, "
      + ("re-derived and pinned field by field."
         if findings_ref is not None
         else "not present yet — write it once with `--write-reference`."))
    A("")
    if problems:
        A("```")
        for p in problems[:40]:
            A(p)
        if len(problems) > 40:
            A(f"... and {len(problems) - 40} more")
        A("```")
        A("")
    A("## Reproduce")
    A("")
    A("```bash")
    A("python3 code/espoo/espoo_ocv_channels_csv.py       # rebuild + verify the table")
    A("python3 code/espoo/espoo_ocv_core.py               # read-only brief of the layer")
    A("python3 code/espoo/fig_espoo_ocv_paper.py          # figures A–E + this report")
    A("python3 code/espoo/fig_espoo_ocv_paper.py --quick  # without the permutation nulls")
    A("python3 code/espoo/fig_espoo_ocv_months.py         # the month companion")
    A("```")
    A("")
    A(f"*environment: python {platform.python_version()}, numpy {np.__version__}, "
      f"matplotlib {matplotlib.__version__}"
      + (", `--quick`" if res["quick"] else "") + "*")
    A("")
    A("*No wall-clock value appears in this file, so two runs of the same code "
      "produce identical bytes.*")
    return "\n".join(md) + "\n"



def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--figdir", default=FIGDIR)
    ap.add_argument("--json", default=OUT_JSON)
    ap.add_argument("--md", default=OUT_MD)
    ap.add_argument("--reference", default=REF_FINDINGS)
    ap.add_argument("--quick", action="store_true",
                    help="skip the label-permutation nulls (fast; no findings pin)")
    ap.add_argument("--write-reference", action="store_true",
                    help="(re)define data/espoo/reference/espoo_ocv_findings.json")
    ap.add_argument("--no-figures", action="store_true",
                    help="write only the JSON, the report and the pins")
    args = ap.parse_args(argv)

    t0 = time.time()
    print(f"computing — {core.SITE}, {core.SHORT[REF_STATE]} -> "
          f"{core.SHORT[GRP_STATE]} (no event) …")
    res = compute(quick=args.quick)
    print(f"  {res['n_rows']} acquisitions ({res['data']['n_by_state']}), "
          f"{res['runtime_s']:.1f} s")

    findings_ref = None
    if args.write_reference:
        print("  --write-reference: (re)defining the package's own lock")
    elif not args.quick and os.path.exists(args.reference):
        with open(args.reference) as fh:
            findings_ref = json.load(fh)

    problems = pin_block(res, findings_ref)
    if problems:
        print(f"PIN FAILED ({len(problems)} problem(s)):")
        for p in problems:
            print(f"  - {p}")
        return 2
    print("  pins OK (findings lock: "
          + ("pinned)" if findings_ref is not None else "not pinned)"))

    if args.write_reference:
        os.makedirs(os.path.dirname(os.path.abspath(args.reference)), exist_ok=True)
        with open(args.reference, "w") as fh:
            json.dump(json_safe(findings(res)), fh, indent=1, sort_keys=True,
                      default=json_default)
        print(f"  {os.path.relpath(args.reference)}")

    if not args.no_figures:
        os.makedirs(args.figdir, exist_ok=True)
        fig_paper(res, args.figdir)

    os.makedirs(os.path.dirname(os.path.abspath(args.json)), exist_ok=True)
    with open(args.json, "w") as fh:
        json.dump(build_json(res, problems, findings_ref, args.write_reference), fh,
                  indent=1, sort_keys=True, default=json_default)
    print(f"  {os.path.relpath(args.json)}")

    os.makedirs(os.path.dirname(os.path.abspath(args.md)), exist_ok=True)
    with open(args.md, "w") as fh:
        fh.write(report(res, problems, findings_ref))
    print(f"  {os.path.relpath(args.md)}")
    print(f"DONE ({time.time() - t0:.1f} s)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
