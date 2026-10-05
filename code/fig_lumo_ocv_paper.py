#!/usr/bin/env python3
"""OCV-Paper — figures A–E, JSON result and English report.

Question: does the observability core vector
``OCV = [gamma2, P, D]`` (three measured quantities of the LUMO mast: interferometric
coherence gamma^2, persistence P and density D of the echo mask) discriminate the four
LUMO structural states (healthy, DAM3, DAM4, DAM6) — or are the states
indistinguishable with respect to these three quantities?

Approach: the same committed 178-overpass sample, the same definitions, the same
seeds and the same procedure as the LUMO project analysis scripts
(``lumo_ocv_stats.py`` is a verbatim port). All numbers are recomputed here —
and then pinned against the five committed reference JSONs: any deviation aborts
the script with exit code != 0. The claim of this repo is therefore not
"we assert", but "we reproduce and read off".

Figures:
  A  distributions of the six OCV components per state (raw data, n=178)
  B  in-sample effect sizes: Cliff's delta + bootstrap CI, SDS (=10|delta|)
  C  4-class LDA, leave-one-out, lambda grid, label permutation, confusion matrix
  D  out-of-fold scores of the six state pairs (Fisher LDA per fold)
  E  controls: severity trend, season/orbit strata, mask brightness, wind

Usage:
  python3 code/fig_lumo_ocv_paper.py            # recompute everything (~8-12 min)
  python3 code/fig_lumo_ocv_paper.py --quick    # without the two permutation tests
"""
from __future__ import annotations

import argparse
import json
import os
import platform
import sys
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

try:                                    # scipy is optional in lumo_ocv_stats
    import scipy
except ImportError:                     # pragma: no cover
    scipy = None

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import lumo_ocv_core as core   # noqa: E402
import lumo_ocv_stats as st    # noqa: E402

FIGDIR = core.FIGDIR
OUT_JSON = os.path.join(core.DATA, "fig_lumo_ocv_paper.json")
OUT_MD = os.path.join(FIGDIR, "fig_lumo_ocv_paper.md")

# ---------------------------------------------------------------------------
# Style (identical to the espoo figures: DejaVu Sans + STIX, 183 mm column width)
# ---------------------------------------------------------------------------
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
MM = 1.0 / 25.4
COL_WIDTH = 183 * MM
STATE_COLORS = {"healthy": "#4C72B0", "DAM 3": "#DD8452",
                "DAM 4": "#C44E52", "DAM 6": "#55A868"}
STATE_MARKERS = {"healthy": "o", "DAM 3": "s", "DAM 4": "^", "DAM 6": "D"}
DIM_COLORS = {"gamma2": "#4C72B0", "A": "#8172B3", "D": "#DD8452",
              "F": "#937860", "S": "#DA8BC3", "P": "#55A868"}
MODEL_COLORS = {"gamma2": "#8C8C8C", "PDS": "#DD8452",
                "gamma2PD": "#4C72B0", "x_6d": "#55A868"}
MODEL_MARKERS = {"gamma2": "o", "PDS": "s", "gamma2PD": "^", "x_6d": "D"}
CHANCE = 1.0 / 4.0
PAIR_KEYS = ["healthy_vs_DAM3", "healthy_vs_DAM4", "healthy_vs_DAM6",
             "DAM3_vs_DAM4", "DAM3_vs_DAM6", "DAM4_vs_DAM6"]
PAIRS = [("healthy", "DAM 3"), ("healthy", "DAM 4"), ("healthy", "DAM 6"),
         ("DAM 3", "DAM 4"), ("DAM 3", "DAM 6"), ("DAM 4", "DAM 6")]
SHORT_MODEL = {"gamma2": "gamma2", "PDS": "P+D+S", "gamma2PD": "gamma2+P+D",
               "x_6d": "6 dims"}


def save_fig(fig, name, out_dir=None):
    """Write one 600 dpi PNG (the espoo set keeps the binary out of git)."""
    out_dir = out_dir or FIGDIR
    png = os.path.join(out_dir, name + ".png")
    fig.savefig(png)
    plt.close(fig)
    print(f"  {os.path.relpath(png, HERE)}")
    return png


# ---------------------------------------------------------------------------
# Compute
# ---------------------------------------------------------------------------
def coh_blocks(rows):
    """The blocks of ``analyze_lumo_tower_coherence`` rebuilt verbatim:
    ``per_state``, ``tests``, ``per_orbit``, ``desc_tests``, ``wind_control``."""
    by_state = st.by_state(rows)
    per_state = {}
    for lab, recs in by_state.items():
        row = {"n": len(recs)}
        if recs:
            g2 = [r["gamma2"] for r in recs]
            npx = [r["n_masked"] for r in recs]
            wind = [r["wind_speed_ms"] for r in recs if r["wind_speed_ms"] is not None]
            row["gamma2"] = st._stats(g2)
            row["n_masked"] = st._stats(npx)
            row["wind_speed_ms"] = st._stats(wind)
            dates = sorted({r["date"] for r in recs if r.get("date")})
            row["date_range"] = [dates[0], dates[-1]] if dates else None
        per_state[core.SHORT_ALC[lab]] = row

    tests = {}
    for other in st.DAMAGED:
        x = [r["gamma2"] for r in by_state["healthy"]]
        y = [r["gamma2"] for r in by_state[other]]
        if len(x) >= 3 and len(y) >= 3:
            tests[f"healthy_vs_{core.SHORT_ALC[other]}"] = st.welch_mw_test(x, y)

    orbit_stats = {}
    for orb in ("ASCENDING", "DESCENDING"):
        orbit_stats[orb] = {}
        for lab in st.STATE_ORDER:
            sub = [r["gamma2"] for r in by_state[lab] if r["orbit"] == orb]
            orbit_stats[orb][core.SHORT_ALC[lab]] = st._stats(sub)

    desc_tests = {}
    for other in st.DAMAGED:
        x = [r["gamma2"] for r in by_state["healthy"] if r["orbit"] == "DESCENDING"]
        y = [r["gamma2"] for r in by_state[other] if r["orbit"] == "DESCENDING"]
        if len(x) >= 3 and len(y) >= 3:
            desc_tests[f"desc_healthy_vs_{core.SHORT_ALC[other]}"] = st.welch_mw_test(x, y)

    wind_control = {}
    if st.sps is not None:
        hw = [(r["wind_speed_ms"], r["gamma2"]) for r in by_state["healthy"]
              if r["wind_speed_ms"] is not None]
        if len(hw) >= 10:
            wx = np.asarray([p[0] for p in hw])
            gy = np.asarray([p[1] for p in hw])
            rho, p = st.sps.spearmanr(wx, gy)
            wind_control["healthy_gamma2_vs_wind"] = {
                "n": len(hw), "spearman_rho": float(rho), "p": float(p)}
        for other in ("DAM 3", "DAM 6"):
            x = [r["gamma2"] for r in by_state["healthy"]
                 if r["wind_speed_ms"] is not None and r["wind_speed_ms"] <= 4.0]
            y = [r["gamma2"] for r in by_state[other]
                 if r["wind_speed_ms"] is not None and r["wind_speed_ms"] <= 4.0]
            if len(x) >= 3 and len(y) >= 3:
                wind_control[f"windmatched_le4_healthy_vs_{core.SHORT_ALC[other]}"] = \
                    st.welch_mw_test(x, y)
    return {"per_state": per_state, "tests": tests, "per_orbit": orbit_stats,
            "desc_tests": desc_tests, "wind_control": wind_control}


def compute(quick=False):
    """Recompute all numbers of figures A–E (and raw values for the pin block)."""
    rows = core.load_csv()
    ref = core.load_reference()
    coh = ref["tower_coherence_states"]
    res = {"rows": rows, "ref": ref}

    X, y, complete = core.feature_dataset(rows)
    res["data"] = {
        "csv": "data/lumo_ocv_channels.csv",
        "csv_meta": "data/lumo_ocv_channels_meta.json",
        "csv_sha256": core.sha256(core.CSV_PATH),
        "n_overpasses": len(rows),
        "n_by_state": core.n_by_state(rows),
        "n_complete_6d": int(X.shape[0]),
        "n_excluded_incomplete": int(len(rows) - X.shape[0]),
        "n_by_state_complete_6d": {
            core.SHORT[l]: int(sum(1 for r in complete if r["damage_label"] == l))
            for l in st.STATE_ORDER},
        "n_P_present": int(sum(1 for r in rows if r["P"] is not None)),
        "mask_shape_majority": [400, 11],
        "n_masked_median_by_state": {
            core.SHORT[l]: float(np.median([r["n_masked"] for r in st.by_state(rows)[l]]))
            for l in st.STATE_ORDER},
        "n_masked_by_state_values": {
            core.SHORT[l]: sorted(int(r["n_masked"]) for r in st.by_state(rows)[l])
            for l in st.STATE_ORDER},
    }

    # --- Figure A: raw distributions per state -------------------------
    per_state = {}
    for key in ["gamma2"] + core.DIMS:
        vals_by_state = {}
        for lab in st.STATE_ORDER:
            v = st.vals(rows, lab, key)
            vals_by_state[core.SHORT[lab]] = st._stats(v)
        per_state[key] = vals_by_state
    res["fig_a"] = {
        "descriptives": per_state,
        "values": {key: {core.SHORT[l]: st.vals(rows, l, key) for l in st.STATE_ORDER}
                   for key in ["gamma2"] + core.DIMS},
        "dim_label": {"gamma2": "gamma2 (coherence)", **core.DIM_LABEL},
    }

    # --- Figure B: in-sample effect sizes ------------------------------
    fig_b = {"gamma2": st.pairs_block(rows, key="gamma2")}
    per_dim, per_dim_full = {}, {}
    for dim in core.DIMS:
        sub = [r for r in rows if r.get(dim) is not None]
        pb = st.pairs_block(sub, key=dim)
        per_dim_full[dim] = pb
        per_dim[dim] = {
            "label": core.DIM_LABEL[dim], "n": len(sub),
            "sds_pooled_all_damaged": pb["sds_pooled_all_damaged"],
            "healthy_vs": {s: {k: b.get(k) for k in ("delta", "sds", "median")}
                           for s, b in pb["healthy_vs"].items()},
            "all_pairs": {k: {kk: b.get(kk) for kk in ("delta", "sds", "median")}
                          for k, b in pb["all_pairs"].items()},
        }
    ranking = sorted(core.DIMS, key=lambda d: -(per_dim[d]["sds_pooled_all_damaged"] or 0.0))
    fig_b["per_dimension"] = per_dim
    fig_b["per_dimension_full"] = per_dim_full
    fig_b["ranking_by_pooled_sds"] = ranking
    fig_b["emv_definition"] = {d: core.DIM_LABEL[d] for d in core.DIMS}
    fig_b["sds_table"] = {
        key: {"pooled": (fig_b["gamma2"] if key == "gamma2" else per_dim_full[key])["sds_pooled_all_damaged"],
              "max_single_group": (fig_b["gamma2"] if key == "gamma2" else per_dim_full[key])["max_single_group"]["sds"],
              "healthy_vs": {s: b.get("sds") for s, b in
                             (fig_b["gamma2"] if key == "gamma2" else per_dim_full[key])["healthy_vs"].items()}}
        for key in ["gamma2"] + core.DIMS}
    res["fig_b"] = fig_b

    # --- Figure E1: trend test + strata (gamma2) -----------------------
    mono = st.monotonicity(rows, fig_b["gamma2"], key="gamma2", n_perm=20000)
    res["fig_e"] = {
        "monotonicity": mono,
        "by_season": st.strata(rows, st.season_of, key="gamma2"),
        "by_orbit": st.strata(rows, lambda r: r.get("orbit"), key="gamma2"),
        "coh": coh_blocks(rows),
        "peak_intensity_median_by_state": {
            core.SHORT[l]: float(np.median([r["peak_intensity"] for r in st.by_state(rows)[l]]))
            for l in st.STATE_ORDER},
        "wind_speed_median_by_state": {
            core.SHORT[l]: float(np.median([r["wind_speed_ms"] for r in st.by_state(rows)[l]]))
            for l in st.STATE_ORDER},
    }
    # --- Figure C: 4-class LDA (LOO, lambda grid, label permutation) ----
    col = {f: i for i, f in enumerate(core.FEATURES)}
    model_specs = {"x_6d": core.FEATURES}
    model_specs.update({f: [f] for f in core.FEATURES})
    grid = {}
    for name, feats in model_specs.items():
        idx = [col[f] for f in feats]
        grid[name] = {}
        for lam in core.LAMBDA_GRID:
            _, acc, bal, _conf = st.loo_cv(X[:, idx], y, core.STATE_ORDER, lam)
            grid[name][str(lam)] = {"accuracy": acc, "balanced_accuracy": bal}
    headline = {}
    for name in ("x_6d", "gamma2"):
        if quick:
            blk = dict(ref["emv_multivariate_cv"]["headline_permutation_test"][name])
            blk["recomputed"] = False
        else:
            idx = [col[f] for f in model_specs[name]]
            blk = st.permutation_test(X[:, idx], y, core.STATE_ORDER, core.HEADLINE_LAMBDA)
            blk["recomputed"] = True
        headline[name] = blk
    _, acc6, bal6, conf6 = st.loo_cv(X, y, core.STATE_ORDER, core.HEADLINE_LAMBDA)
    res["fig_c"] = {
        "sensitivity_grid": grid, "headline_permutation_test": headline,
        "headline_confusion_matrix_6d": conf6,
        "headline_accuracy_6d": acc6, "headline_balanced_accuracy_6d": bal6,
        "lambda_grid": core.LAMBDA_GRID, "headline_lambda": core.HEADLINE_LAMBDA,
        "chance_balanced_accuracy": CHANCE, "permutation_recomputed": not quick,
    }

    # --- Figure D: out-of-fold binary state pairs ------------------------
    pairs_res, q_by_pair = {}, {}
    for class_a, class_b in PAIRS:
        pair_key = f"{st.STATE_SHORT[class_a]}_vs_{st.STATE_SHORT[class_b]}"
        sub = np.isin(y, [class_a, class_b])
        Xp_full, yp = X[sub], y[sub]
        n_a, n_b = int((yp == class_a).sum()), int((yp == class_b).sum())
        per_model, q_by_model = {}, {}
        for name, feats in core.MODELS.items():
            idx = [col[f] for f in feats]
            Xp = Xp_full[:, idx]
            q = st.oof_scores(Xp, yp, class_a, class_b, core.HEADLINE_LAMBDA)
            q_by_model[name] = q
            delta, sds, ref_q, grp_q = st.delta_sds(q, yp, class_a, class_b)
            ci = (st.bootstrap_delta_ci(ref_q, grp_q, n=core.N_BOOT, seed=core.RNG_SEED)
                  if delta is not None else None)
            if quick:
                b_ref = ref["emv_pairwise_triplets_cv"]["pairs"][pair_key]["per_model"][name]
                p_val, null_mean, null_p95, p_recomputed = (
                    b_ref["perm_p"], b_ref["perm_null_mean_sds"],
                    b_ref["perm_null_p95_sds"], False)
            else:
                p_val, null_mean, null_p95 = st.permutation_p(
                    Xp, yp, class_a, class_b, core.HEADLINE_LAMBDA, sds,
                    core.N_PERM_PAIR, np.random.default_rng(core.RNG_SEED + 1))
                p_recomputed = True
            per_model[name] = {
                "label": core.MODEL_LABEL[name], "n": [n_a, n_b],
                "delta": delta, "sds": sds,
                "delta_ci95": list(ci) if ci else None,
                "ci_excludes_zero": bool(ci[0] * ci[1] > 0) if ci else None,
                "perm_p": p_val, "perm_null_mean_sds": null_mean,
                "perm_null_p95_sds": null_p95, "perm_recomputed": p_recomputed,
            }
        paired = {}
        for name in (core.HEADLINE_PAIR, core.SECONDARY_PAIR, "x_6d"):
            paired[name] = st.paired_bootstrap_diff(
                q_by_model[name], q_by_model[core.BASELINE], yp, class_a, class_b,
                core.N_BOOT, np.random.default_rng(core.RNG_SEED + 2))
        pairs_res[pair_key] = {"class_a": st.STATE_SHORT[class_a],
                               "class_b": st.STATE_SHORT[class_b], "n": [n_a, n_b],
                               "per_model": per_model, "paired_diff_vs_gamma2": paired}
        q_by_pair[pair_key] = {"class_a": class_a, "class_b": class_b,
                               "y": list(yp), "q": q_by_model}
    res["fig_d"] = {"pairs": pairs_res, "q_by_pair": q_by_pair}
    return res, X, y, complete


# ---------------------------------------------------------------------------
# Pin block: compares every recomputed number with the committed reference
# ---------------------------------------------------------------------------
def pin_all(res, quick=False):
    """All pins against the five committed reference JSONs."""
    ref, rows = res["ref"], res["rows"]
    coh, gp = ref["tower_coherence_states"], ref["gamma2_pairs"]
    emv, mcv, tri = (ref["echo_mask_vector"], ref["emv_multivariate_cv"],
                     ref["emv_pairwise_triplets_cv"])
    meta = core.load_csv_meta()
    b = coh["bursts"]
    P = Pins()

    # 1) the committed 178-overpass sample, row by row ------------------------
    P.cmp("coh.bursts.burst_id", [r["burst_id"] for r in rows], [x["burst_id"] for x in b])
    P.cmp("coh.bursts.date", [r["date"] for r in rows], [x["date"] for x in b])
    P.cmp("coh.bursts.month", [r["month"] for r in rows], [x["month"] for x in b])
    P.cmp("coh.bursts.orbit", [r["orbit"] for r in rows], [x["orbit"] for x in b])
    P.cmp("coh.bursts.damage_label", [r["damage_label"] for r in rows],
          [x["damage_label"] for x in b])
    P.cmp("coh.bursts.gamma2", [r["gamma2"] for r in rows], [x["gamma2"] for x in b])
    P.cmp("coh.bursts.n_masked", [r["n_masked"] for r in rows], [x["n_masked"] for x in b])
    P.cmp("coh.bursts.peak_intensity", [r["peak_intensity"] for r in rows],
          [x["peak_intensity"] for x in b])
    P.cmp("csv.A_equals_committed_n_masked", [r["A"] for r in rows],
          [float(x["n_masked"]) for x in b])
    P.cmp("coh.params", meta["mask_definition"], coh["params"])
    P.cmp("coh.counts", {"n_bursts_cached": meta["n_cache_dirs_matched"],
                         "n_bursts_total": len(rows),
                         "n_skipped": meta["n_cache_skipped"]},
          {"n_bursts_cached": coh["n_bursts_cached"],
           "n_bursts_total": coh["n_bursts_total"], "n_skipped": coh["n_skipped"]})
    P.cmp("coh.per_state", res["fig_e"]["coh"]["per_state"], coh["per_state"])
    P.cmp("coh.tests", res["fig_e"]["coh"]["tests"], coh["tests"])
    P.cmp("coh.per_orbit", res["fig_e"]["coh"]["per_orbit"], coh["per_orbit"])
    P.cmp("coh.desc_tests", res["fig_e"]["coh"]["desc_tests"], coh["desc_tests"])
    P.cmp("coh.wind_control", res["fig_e"]["coh"]["wind_control"], coh["wind_control"])

    # 2) gamma2 pairs, trend test and strata (lumo_gamma2_pairs.json) ---------
    P.cmp("gp.primary", res["fig_b"]["gamma2"], gp["primary"])
    P.cmp("gp.monotonicity", res["fig_e"]["monotonicity"], gp["monotonicity"])
    P.cmp("gp.controls.by_season", res["fig_e"]["by_season"], gp["controls"]["by_season"])
    P.cmp("gp.controls.by_orbit", res["fig_e"]["by_orbit"], gp["controls"]["by_orbit"])
    P.cmp("gp.controls.n_masked_median_by_state", res["data"]["n_masked_median_by_state"],
          gp["controls"]["n_masked_median_by_state"])
    P.cmp("gp.controls.peak_intensity_median_by_state",
          res["fig_e"]["peak_intensity_median_by_state"],
          gp["controls"]["peak_intensity_median_by_state"])
    P.cmp("gp.params", {"peak_frac": mask_definition("peak_frac"),
                        "median_mult": mask_definition("median_mult"),
                        "n_perm": 20000, "seed": 7, "primary_weighting": "overpass"},
          gp["params"])

    # 3) EMV aggregates (lumo_echo_mask_vector.json) --------------------------
    P.cmp("emv.per_dimension", res["fig_b"]["per_dimension"], emv["per_dimension"])
    P.cmp("emv.ranking_by_pooled_sds", res["fig_b"]["ranking_by_pooled_sds"],
          emv["ranking_by_pooled_sds"])
    P.cmp("emv.emv_definition", res["fig_b"]["emv_definition"], emv["emv_definition"])
    P.cmp("emv.meta.n_by_state", res["data"]["n_by_state"], emv["meta"]["n_by_state"])
    P.cmp("emv.meta.persistence", meta["persistence"], emv["meta"]["persistence"])
    P.cmp("emv.meta.cross_check_n_masked_vs_committed",
          meta["cross_check_n_masked_vs_committed"],
          emv["meta"]["cross_check_n_masked_vs_committed"])

    # 4) 4-class CV (lumo_emv_multivariate_cv.json) ---------------------------
    P.cmp("mcv.sensitivity_grid", res["fig_c"]["sensitivity_grid"], mcv["sensitivity_grid"])
    P.cmp("mcv.headline_accuracy_6d", res["fig_c"]["headline_accuracy_6d"],
          mcv["headline_accuracy_6d"])
    P.cmp("mcv.headline_balanced_accuracy_6d", res["fig_c"]["headline_balanced_accuracy_6d"],
          mcv["headline_balanced_accuracy_6d"])
    P.cmp("mcv.headline_confusion_matrix_6d", res["fig_c"]["headline_confusion_matrix_6d"],
          mcv["headline_confusion_matrix_6d"])
    P.cmp("mcv.lambda_grid", res["fig_c"]["lambda_grid"], mcv["lambda_grid"])
    P.cmp("mcv.headline_lambda", res["fig_c"]["headline_lambda"], mcv["headline_lambda"])
    P.cmp("mcv.features", list(core.FEATURES), mcv["features"])
    P.cmp("mcv.classes", list(core.STATE_ORDER), mcv["classes"])
    P.cmp("mcv.meta.n_total_committed", res["data"]["n_overpasses"],
          mcv["meta"]["n_total_committed"])
    P.cmp("mcv.meta.n_complete", res["data"]["n_complete_6d"], mcv["meta"]["n_complete"])
    P.cmp("mcv.meta.n_excluded_incomplete", res["data"]["n_excluded_incomplete"],
          mcv["meta"]["n_excluded_incomplete"])
    P.cmp("mcv.meta.n_by_state", res["data"]["n_by_state_complete_6d"],
          mcv["meta"]["n_by_state"])
    if not quick:
        P.cmp("mcv.headline_permutation_test", res["fig_c"]["headline_permutation_test"],
              mcv["headline_permutation_test"])

    # 5) out-of-fold pair separation (lumo_emv_pairwise_triplets_cv.json) -----
    P.cmp("tri.models", core.MODEL_LABEL, tri["models"])
    P.cmp("tri.shrinkage", core.HEADLINE_LAMBDA, tri["shrinkage"])
    P.cmp("tri.n_perm", core.N_PERM_PAIR, tri["n_perm"])
    P.cmp("tri.n_boot", core.N_BOOT, tri["n_boot"])
    P.cmp("tri.meta.n_total_committed", res["data"]["n_overpasses"],
          tri["meta"]["n_total_committed"])
    P.cmp("tri.meta.n_complete", res["data"]["n_complete_6d"], tri["meta"]["n_complete"])
    P.cmp("tri.meta.n_by_state", res["data"]["n_by_state_complete_6d"],
          tri["meta"]["n_by_state"])
    for pk in PAIR_KEYS:
        got, want = res["fig_d"]["pairs"][pk], tri["pairs"][pk]
        for k in ("class_a", "class_b", "n"):
            P.cmp(f"tri.pairs.{pk}.{k}", got[k], want[k])
        for model in core.MODELS:
            g = dict(got["per_model"][model])
            g.pop("perm_recomputed", None)
            w = dict(want["per_model"][model])
            if quick:
                for k in ("perm_p", "perm_null_mean_sds", "perm_null_p95_sds"):
                    w.pop(k, None)
            P.cmp(f"tri.pairs.{pk}.per_model.{model}", g, w)
        P.cmp(f"tri.pairs.{pk}.paired_diff_vs_gamma2", got["paired_diff_vs_gamma2"],
              want["paired_diff_vs_gamma2"])
    return P


def mask_definition(key):
    return core.load_csv_meta()["mask_definition"][key]
FLOAT_TOL = 1e-6
TOL_KEYS = {"welch_t", "welch_p", "mannwhitney_p", "perm_p", "p", "rho",
            "p_accuracy", "p_balanced_accuracy", "p_decreasing", "p_increasing",
            "p_two_sided", "accuracy_null_mean", "balanced_accuracy_null_mean",
            "accuracy_null_p95", "balanced_accuracy_null_p95",
            "perm_null_mean_sds", "perm_null_p95_sds", "diff_mean", "ci95"}


class Pins:
    """Collects pin comparisons (recomputed vs. committed) and reports deviations.

    Comparison rule: exact equality for everything deterministic (delta, SDS,
    medians, counters, CI endpoints, gamma2); 1e-6 for numbers that arise from
    scipy tests or permutation p-values (library version / normalisation).
    Deviations above the tolerance => exit code != 0.
    """

    def __init__(self):
        self.checks = []
        self.fails = []
        self.extra_keys = []
        self.label_translations = 0

    def add(self, path, got, want, tol):
        # Committed German doc strings (see core.LABEL_EN) are compared in their
        # English translation; the count is reported in the summary so that no
        # translation can silently hide a genuine deviation.
        if isinstance(got, str) and got in core.LABEL_EN:
            got = core.LABEL_EN[got]
            self.label_translations += 1
        if isinstance(want, str) and want in core.LABEL_EN:
            want = core.LABEL_EN[want]
            self.label_translations += 1
        exact = bool(got == want)
        dev = 0.0
        if not exact:
            numeric = (isinstance(got, (int, float)) and isinstance(want, (int, float))
                       and not isinstance(got, bool) and not isinstance(want, bool))
            dev = abs(float(got) - float(want)) if numeric else None
        ok = exact or (dev is not None and dev <= tol)
        self.checks.append({"path": path, "kind": "exact" if exact else "tolerance",
                            "abs_deviation": dev, "tolerance": tol, "ok": ok})
        if not ok:
            self.fails.append({"path": path, "recomputed": repr(got),
                               "committed": repr(want), "abs_deviation": dev})

    def cmp(self, prefix, got, want, key=None):
        """Recursive comparison: `want` (committed) dictates structure and expectation."""
        tol = FLOAT_TOL if key in TOL_KEYS else 0.0
        if isinstance(want, dict):
            if not isinstance(got, dict):
                self.add(prefix, f"<{type(got).__name__}>", "<dict>", 0.0)
                return
            for k in want:
                if k not in got:
                    self.add(f"{prefix}.{k}", "<missing>", "<present>", 0.0)
                    continue
                self.cmp(f"{prefix}.{k}", got[k], want[k], k)
            for k in got:
                if k not in want:
                    self.extra_keys.append(f"{prefix}.{k}")
            return
        if isinstance(want, list):
            glen = len(got) if hasattr(got, "__len__") else -1
            if not isinstance(got, (list, tuple)) or glen != len(want):
                self.add(prefix, f"<list len {glen}>", f"<list len {len(want)}>", 0.0)
                return
            for i, (g, w) in enumerate(zip(got, want)):
                self.cmp(f"{prefix}[{i}]", g, w, key)
            return
        self.add(prefix, got, want, tol)

    def summary(self):
        devs = [c["abs_deviation"] for c in self.checks
                if c["abs_deviation"] not in (None, 0.0)]
        return {
            "n_checks": len(self.checks),
            "n_exact": sum(1 for c in self.checks if c["kind"] == "exact"),
            "n_within_tolerance": sum(1 for c in self.checks if c["kind"] == "tolerance"),
            "n_failed": len(self.fails),
            "max_abs_deviation": max(devs) if devs else 0.0,
            "label_translations": self.label_translations,
            "float_tolerance": FLOAT_TOL,
            "tolerance_keys": sorted(TOL_KEYS),
            "failures": self.fails,
            "recomputed_keys_not_in_reference": sorted(set(self.extra_keys)),
        }




# ---------------------------------------------------------------------------
# Rendering — figures A–E, 183 mm column width, 600 dpi PNG (espoo style)
# ---------------------------------------------------------------------------
DIM_AXIS = {"gamma2": "gamma^2", "A": "A  [px]", "D": "D  [1]",
            "F": "F  [components]", "S": "S  [px]", "P": "P  [1]"}
DIM_TITLE = {"gamma2": "interferometric coherence gamma^2", "A": "mask area",
             "D": "mask density", "F": "mask fragmentation",
             "S": "centroid<->peak shift", "P": "mask persistence"}


def panel(ax, letter, text):
    ax.set_title(f"({letter}) {text}", fontsize=8.0, loc="left")


def jitter(n, width, seed):
    rng = np.random.default_rng(seed)
    return rng.uniform(-width, width, n)


def state_strip(ax, values_by_state, offset=0.0, span=0.16, seed0=100, rot=0,
                fs=6.6):
    """Raw per-overpass dots + median bar for each of the four states."""
    for i, lab in enumerate(st.STATE_ORDER):
        v = np.asarray(values_by_state[core.SHORT[lab]], float)
        if v.size == 0:
            continue
        x = offset + i + jitter(v.size, span, seed0 + i)
        ax.plot(x, v, STATE_MARKERS[lab], ms=2.4, mfc="none",
                mec=STATE_COLORS[lab], mew=0.5, alpha=0.6, ls="none")
        ax.hlines(float(np.median(v)), offset + i - 0.30, offset + i + 0.30,
                  color=STATE_COLORS[lab], lw=1.7)
    ax.set_xlim(offset - 0.55, offset + 3.55)
    ax.set_xticks([offset + i for i in range(4)])
    ax.set_xticklabels([core.SHORT[l] for l in st.STATE_ORDER], fontsize=fs,
                       rotation=rot, ha="right" if rot else "center",
                       rotation_mode="anchor" if rot else None)


def fig_a(res, out_dir=None):
    """Figure A — raw distributions of the six OCV components per state."""
    vals_ = res["fig_a"]["values"]
    keys = ["gamma2"] + core.DIMS
    states = [core.SHORT[l] for l in st.STATE_ORDER]
    fig, axes = plt.subplots(2, 3, figsize=(COL_WIDTH, 0.80 * COL_WIDTH))
    for ax, key in zip(axes.ravel(), keys):
        state_strip(ax, vals_[key], seed0=10 + keys.index(key))
        if key == "A":
            ax.set_yscale("log")
        elif key == "P":
            ax.set_ylim(0.0, 1.02)
        ax.set_ylabel(DIM_AXIS[key], fontsize=7.0)
        panel(ax, "abcdef"[keys.index(key)], DIM_TITLE[key])
    n_by_state = res["data"]["n_by_state"]
    fig.suptitle(
        "Figure A — OCV components per LUMO structural state  "
        f"({res['data']['n_overpasses']} coherence overpasses; "
        + ", ".join(f"{s} n={n_by_state[s]}" for s in states)
        + ")\n"
        "dots = single overpasses, bar = median;  input: data/lumo_ocv_channels.csv",
        fontsize=7.6)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    return save_fig(fig, "fig_lumo_ocv_paper_A", out_dir)


def fig_b(res, out_dir=None):
    """Figure B — in-sample effect sizes (Cliff's delta, bootstrap CI, SDS)."""
    fb = res["fig_b"]
    tab = fb["sds_table"]
    rank = sorted(["gamma2"] + core.DIMS,
                  key=lambda k: -(tab[k]["pooled"] or 0.0))
    fig, axes = plt.subplots(1, 2, figsize=(COL_WIDTH, 0.40 * COL_WIDTH),
                             width_ratios=[1.0, 1.05])
    ax = axes[0]
    y = np.arange(len(rank))
    ax.barh(y, [tab[k]["pooled"] or 0.0 for k in rank], height=0.62,
            color=[DIM_COLORS[k] for k in rank])
    for i, k in enumerate(rank):
        ax.text((tab[k]["pooled"] or 0.0) + 0.15, i, f"{tab[k]['pooled']:.2f}",
                va="center", fontsize=6.4)
    ax.set_yticks(y)
    ax.set_yticklabels([DIM_AXIS[k] for k in rank], fontsize=7.0)
    ax.invert_yaxis()
    ax.set_xlim(0.0, max([tab[k]["pooled"] or 0.0 for k in rank] + [1.0]) * 1.22)
    ax.set_xlabel("SDS, healthy vs. all damaged pooled", fontsize=7.0)
    panel(ax, "a", "which component separates at all")

    ax = axes[1]
    keys = ["gamma2"] + core.DIMS
    for i, k in enumerate(keys):
        for j, lab in enumerate(st.DAMAGED):
            sds = ((tab[k]["healthy_vs"] or {}).get(core.SHORT[lab]))
            if sds is None:
                continue
            ax.plot(i + (j - 1) * 0.22, sds, STATE_MARKERS[lab], ms=4.2,
                    color=STATE_COLORS[lab],
                    label=core.SHORT[lab] if i == 0 else None)
    ax.set_xticks(range(len(keys)))
    ax.set_xticklabels([DIM_AXIS[k] for k in keys], fontsize=6.6, rotation=18,
                       ha="right")
    ax.set_ylabel("SDS, healthy vs. one state", fontsize=7.0)
    ax.legend(fontsize=6.2, ncol=3, loc="upper right")
    panel(ax, "b", "per damage state (SDS = 10 |Cliff's delta|)")
    sgn = fb["gamma2"]["sign_decomposition"]
    fig.suptitle(
        "Figure B — in-sample effect sizes of the OCV components "
        f"(n={res['data']['n_overpasses']} overpasses; SDS 0 = indistinguishable, "
        "SDS >= 2 = large)\n"
        "pooled = healthy vs. all damaged overpasses at once; gamma^2 cancels "
        "there, because its three contributions have opposite signs "
        "(Cliff's delta "
        + ", ".join(f"{k} {sgn[k]:+.2f}" for k in ("DAM3", "DAM4", "DAM6"))
        + ")", fontsize=7.2)
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    return save_fig(fig, "fig_lumo_ocv_paper_B", out_dir)


def fig_c(res, out_dir=None):
    """Figure C — 4-class LDA: lambda sweep, confusion matrix, permutation null."""
    fc = res["fig_c"]
    grid = fc["sensitivity_grid"]
    lambdas = [str(x) for x in fc["lambda_grid"]]
    fig, axes = plt.subplots(1, 3, figsize=(COL_WIDTH, 0.38 * COL_WIDTH),
                             width_ratios=[1.15, 1.0, 1.0])
    ax = axes[0]
    others = [k for k in grid if k not in ("x_6d", "gamma2")]
    for name in others:
        ax.plot(range(len(lambdas)), [grid[name][l]["balanced_accuracy"]
                                      for l in lambdas], "-", lw=0.9,
                color="0.75", zorder=1)
    for name in ("gamma2", "x_6d"):
        ax.plot(range(len(lambdas)), [grid[name][l]["balanced_accuracy"]
                                      for l in lambdas], "-o", ms=3.2, lw=1.4,
                color=MODEL_COLORS[name],
                label={"gamma2": "gamma2 (1D)",
                       "x_6d": "6D (all channels)"}[name], zorder=3)
    ax.axhline(fc["chance_balanced_accuracy"], color="0.35", ls=":", lw=1.0)
    ax.text(len(lambdas) - 0.5, fc["chance_balanced_accuracy"] + 0.015,
            "chance", ha="right", fontsize=6.2, color="0.35")
    ax.set_xticks(range(len(lambdas)))
    ax.set_xticklabels([f"{float(l):g}" for l in lambdas], fontsize=6.6)
    ax.set_xlabel("shrinkage lambda", fontsize=7.0)
    ax.set_ylabel("balanced accuracy (LOO)", fontsize=7.0)
    ax.set_ylim(0.0, 0.42)
    ax.legend(fontsize=6.0, loc="upper right")
    panel(ax, "a", "6D vs gamma^2 vs single components")

    ax = axes[1]
    conf = fc["headline_confusion_matrix_6d"]
    classes = list(conf)
    M = np.array([[conf[a][b] for b in classes] for a in classes], float)
    im = ax.imshow(M, cmap="Blues", vmin=0.0)
    for i in range(len(classes)):
        for j in range(len(classes)):
            ax.text(j, i, f"{int(M[i, j])}", ha="center", va="center",
                    fontsize=6.2,
                    color="white" if M[i, j] > M.max() * 0.6 else "black")
    ax.set_xticks(range(len(classes)))
    ax.set_xticklabels([core.SHORT[c] for c in classes], fontsize=6.4)
    ax.set_yticks(range(len(classes)))
    ax.set_yticklabels([core.SHORT[c] for c in classes], fontsize=6.4)
    ax.set_xlabel("predicted", fontsize=7.0)
    ax.set_ylabel("true", fontsize=7.0)
    ax.grid(False)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03).ax.tick_params(labelsize=6)
    panel(ax, "b", f"LOO confusion, 6D, lambda={fc['headline_lambda']:g}")

    ax = axes[2]
    # The committed permutation test stores only its summary statistics (mean and
    # p95 of the null), so the panel is drawn as a number line rather than as a
    # histogram of a distribution this repository does not have.
    blk = fc["headline_permutation_test"]["x_6d"]
    obs = blk["balanced_accuracy_observed"]
    nm = blk["balanced_accuracy_null_mean"]
    n95 = blk["balanced_accuracy_null_p95"]
    ch = fc["chance_balanced_accuracy"]
    ax.set_xlim(min(obs, ch, nm) - 0.07, max(n95, ch) + 0.07)
    ax.set_ylim(0.0, 1.0)
    ax.set_yticks([])
    ax.axhline(0.5, color="0.88", lw=6, solid_capstyle="butt", zorder=0)
    for val, lab, color, ls, lw, up, ha in (
            (obs, "observed", MODEL_COLORS["x_6d"], "-", 1.8, True, "center"),
            (nm, "null mean", "0.45", "--", 1.2, False, "right"),
            (n95, "null p95", "0.45", ":", 1.2, False, "left"),
            (ch, "chance", "#DD8452", "-", 1.2, True, "center")):
        ax.vlines(val, 0.5 - 0.09, 0.5 + 0.09, color=color, lw=lw, ls=ls,
                  zorder=3)
        ax.text(val, 0.5 + (0.12 if up else -0.12), lab, fontsize=5.8,
                color=color, ha=ha, va="bottom" if up else "top")
    ax.text(0.5, 0.04,
            f"observed {obs:.3f}  vs.  chance {ch:.3f}\n"
            f"null mean {nm:.3f}, null p95 {n95:.3f}, "
            f"p = {blk['p_balanced_accuracy']:.3f} "
            f"({blk.get('n_perm', core.N_PERM_4CLASS)} permutations)",
            transform=ax.transAxes, ha="center", va="bottom", fontsize=5.8)
    ax.set_xlabel("balanced accuracy", fontsize=7.0)
    tag = "recomputed" if blk.get("recomputed") else "copied (--quick)"
    panel(ax, "c", f"label permutation, p={blk['p_balanced_accuracy']:.3f} [{tag}]")
    fig.suptitle(
        "Figure C — 4-class LDA (LOO) over the six OCV components: "
        "can the state be recovered at all?", fontsize=7.6)
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    return save_fig(fig, "fig_lumo_ocv_paper_C", out_dir)


def fig_d(res, out_dir=None):
    """Figure D — out-of-fold scores of the six state pairs (Fisher LDA per fold)."""
    fd = res["fig_d"]["pairs"]
    qb = res["fig_d"]["q_by_pair"]
    models = list(core.MODELS)
    fig, axes = plt.subplots(2, 3, figsize=(COL_WIDTH, 0.85 * COL_WIDTH))
    for ax, pk in zip(axes.ravel(), PAIR_KEYS):
        d = fd[pk]
        q = qb[pk]
        ya, yb = q["class_a"], q["class_b"]
        for mi, model in enumerate(models):
            sc = np.asarray(q["q"][model], float)
            for gi, cls in enumerate((ya, yb)):
                m = np.asarray([v == cls for v in q["y"]])
                ax.plot(mi + (gi - 0.5) * 0.38 + jitter(int(m.sum()), 0.07, 7 * mi + gi),
                        sc[m], STATE_MARKERS[cls], ms=2.3, mfc="none",
                        mec=STATE_COLORS[cls], mew=0.5, alpha=0.7, ls="none")
            sds = d["per_model"][model]["sds"]
            ax.text(mi, ax.get_ylim()[0], f"{sds:.1f}", ha="center", va="bottom",
                    fontsize=6.0, color=MODEL_COLORS[model])
        ax.set_xticks(range(len(models)))
        ax.set_xticklabels([SHORT_MODEL.get(m, m) for m in models],
                           fontsize=6.2, rotation=20, ha="right")
        ax.axhline(0.0, color="0.6", lw=0.8, ls=":")
        p_ = d["per_model"][core.HEADLINE_PAIR]["perm_p"]
        panel(ax, pk.replace("_vs_", " vs "),
              f"n={d['n'][0]}/{d['n'][1]}\n"
              f"p({SHORT_MODEL[core.HEADLINE_PAIR]})={p_:.3f}")
    fig.suptitle(
        "Figure D — out-of-fold Fisher-LDA scores per state pair "
        "(labels: SDS of the headline model [gamma2, P, D])\n"
        "coloured markers = class_a (o/s) vs class_b (^/D) as in the panel title",
        fontsize=7.6)
    fig.tight_layout(rect=(0, 0, 1, 0.91))
    return save_fig(fig, "fig_lumo_ocv_paper_D", out_dir)


def fig_e(res, out_dir=None):
    """Figure E — controls: severity trend, strata, mask/wind covariates."""
    fe = res["fig_e"]
    fig, axes = plt.subplots(1, 3, figsize=(COL_WIDTH, 0.38 * COL_WIDTH))
    ax = axes[0]
    g2 = res["fig_a"]["values"]["gamma2"]
    state_strip(ax, g2, seed0=31, rot=25, fs=6.2)
    mono = fe["monotonicity"]
    sp = mono.get("spearman_severity_vs_value") or {}
    jt = mono["jonckheere"]
    ax.set_ylabel("gamma^2", fontsize=7.0)
    ax.text(0.32, 0.98,
            f"Spearman rho(severity) = {sp.get('rho', float('nan')):.3f}  "
            f"(p = {sp.get('p', float('nan')):.2f})\n"
            f"Jonckheere p_up = {jt['p_increasing']:.3f}  "
            f"p_down = {jt['p_decreasing']:.3f}\n"
            f"claim |d3| < |d4| < |d6|: "
            f"{'holds' if mono['hypothesis']['holds'] else 'does not hold'}",
            transform=ax.transAxes, fontsize=5.8, va="top")
    panel(ax, "a", "severity trend (gamma^2)")

    ax = axes[1]
    base = res["fig_b"]["sds_table"]["gamma2"]["healthy_vs"]
    entries = [("all", None)]
    entries += [(s, fe["by_season"][s]) for s in fe["by_season"]]
    entries += [(s, fe["by_orbit"][s]) for s in fe["by_orbit"]]
    xs = np.arange(len(entries))
    for lab in st.DAMAGED:
        key = core.SHORT[lab]
        ys = []
        for _name, blk in entries:
            sds = base.get(key) if blk is None else (blk.get(key) or {}).get("sds")
            ys.append(np.nan if sds is None else sds)
        ax.plot(xs, np.asarray(ys, float), STATE_MARKERS[lab] + "-", ms=3.6,
                lw=1.1, color=STATE_COLORS[lab], label=key)
    ax.axvline(0.5, color="0.8", lw=0.8)
    ax.axvline(0.5 + len(fe["by_season"]), color="0.8", lw=0.8)
    ax.set_xticks(xs)
    ax.set_xticklabels([n.replace("_", "\n") for n, _b in entries], fontsize=5.2,
                       rotation=35, ha="right")
    ax.set_ylabel("SDS, healthy vs. state", fontsize=7.0)
    panel(ax, "b", "strata controls (gamma^2)")

    ax = axes[2]
    states = [core.SHORT[l] for l in st.STATE_ORDER]
    x = np.arange(len(states))
    ax.bar(x - 0.19, [fe["peak_intensity_median_by_state"][s] / 1e6 for s in states],
           width=0.36, color="#8172B3", label="median peak intensity")
    ax.set_ylabel("median peak intensity  [1e6]", fontsize=7.0)
    ax.set_xticks(x)
    ax.set_xticklabels(states, fontsize=6.4, rotation=20, ha="right")
    ax2 = ax.twinx()
    ax2.plot(x + 0.19, [fe["wind_speed_median_by_state"][s] for s in states],
             "D-", ms=3.4, lw=1.2, color="#55A868")
    ax2.set_ylabel("median wind speed  [m/s]", fontsize=7.0, color="#55A868")
    ax2.tick_params(axis="y", labelsize=6.4, colors="#55A868")
    ax2.grid(False)
    wc = fe["coh"]["wind_control"]
    ax.margins(y=0.34)
    ax.text(0.02, 0.99,
            "healthy vs. (wind <= 4 m/s):\n"
            + "\n".join(
                f"  {k.split('_vs_')[1]}: p={v['welch_p']:.3f}  "
                f"n={v['n'][0]}/{v['n'][1]}"
                for k, v in wc.items() if k.startswith("windmatched_")),
            transform=ax.transAxes, fontsize=5.4, va="top")
    panel(ax, "c", "mask brightness and wind covariates")
    fig.suptitle(
        "Figure E — controls: is anything beyond the state driving the separation?\n"
        "markers: orange squares = DAM3, red triangles = DAM4, green diamonds = DAM6",
        fontsize=7.6)
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    return save_fig(fig, "fig_lumo_ocv_paper_E", out_dir)


# ---------------------------------------------------------------------------
# Markdown report (no wall-clock info: two runs must give identical bytes)
# ---------------------------------------------------------------------------
def md_table(header, rows):
    head = "| " + " | ".join(header) + " |"
    rule = "|" + "|".join(["---"] * len(header)) + "|"
    body = ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]
    return "\n".join([head, rule] + body)


def fmt(x, nd=3):
    if x is None:
        return "n/a"
    if isinstance(x, bool):
        return "yes" if x else "no"
    if isinstance(x, (int, np.integer)):
        return str(int(x))
    return f"{float(x):.{nd}f}"


def report(res, pins, quick):
    """Deterministic Markdown report of figures A–E plus the pin verdict."""
    d = res["data"]
    fb, fc, fe = res["fig_b"], res["fig_c"], res["fig_e"]
    states = [core.SHORT[l] for l in st.STATE_ORDER]
    L = []
    A = L.append

    def dsgn(x):
        """Signed Cliff's delta for the prose: '+0.139', '-0.053', '0.000'."""
        return f"{x:+.3f}" if x else f"{x:.3f}"

    A("# OCV observability package — do [gamma^2, P, D] separate the LUMO states?")
    A("")
    A("Figures A–E of `research/observability-core-vector/`. This file is generated "
      "by `code/fig_lumo_ocv_paper.py` and deliberately contains no wall-clock "
      "timestamp, so two runs produce identical bytes (the figure PNGs are "
      "deterministic for the same reason).")
    A("")
    A("## 1. Question and design")
    A("")
    A("The **observability core vector** is "
      "OCV = [gamma^2, P, D] = [interferometric coherence, mask persistence, mask "
      "density] — one radiometric and two mask-derived coordinates that are "
      "expected to respond monotonically to progressing structural damage. All "
      "other extracted channels are treated as *controls*: A (mask area), "
      "F (mask fragmentation) and S (centroid<->peak shift) are mask-morphology "
      "quantities, i.e. geometric proxies rather than observability of the "
      "structure itself.")
    A("")
    A("The question of this package is therefore not *\"is there a difference?\"* "
      "(with n=178 overpasses some difference is always detectable) but:")
    A("")
    A("1. **Does the core vector carry the separation?** — Figure B: how large is "
      "the effect per component, and how do P and D compare with gamma^2?")
    A("2. **Is a state recovery possible at all?** — Figure C: 4-class LDA with "
      "leave-one-out cross-validation and a label-permutation null.")
    A("3. **Is the separation state-specific?** — Figure D: six binary state pairs, "
      "out-of-fold scores, four feature sets.")
    A("4. **Is it the state that drives it?** — Figure E: severity trend, season / "
      "orbit strata, mask brightness and wind covariates.")
    A("")
    A("## 2. Data")
    A("")
    A(md_table(["quantity", "value"], [
        ["overpass table", f"`{d['csv']}` (+ `{d['csv_meta']}`)"],
        ["sha256 of the table", f"`{d['csv_sha256']}`"],
        ["overpasses (coherence bursts)", d["n_overpasses"]],
        ["states", ", ".join(f"{s} n={d['n_by_state'][s]}" for s in states)],
        ["complete 6-channel rows", d["n_complete_6d"]],
        ["excluded (incomplete)", d["n_excluded_incomplete"]],
        ["rows with a P value", d["n_P_present"]],
        ["mask geometry", f"{d['mask_shape_majority'][0]} x "
                          f"{d['mask_shape_majority'][1]} px "
                          "(majority shape of the 177 masks)"],
        ["median masked pixels A per state",
         ", ".join(f"{s} {d['n_masked_median_by_state'][s]:.0f}" for s in states)],
    ]))
    A("")
    A("Mask definition (all three thresholds in one rule, identical to the LUMO "
      "pipeline): pixels above `0.30 x peak`, above `5 x median`, at least 2 "
      "pixels; P is the fraction of the 177 cached majority masks in which the "
      "pixel belongs to the mask (one outlier overpass excluded), A its pixel "
      "count, D = A / bounding-box area, F = 8-connected components, "
      "S = distance centroid <-> coherence peak.")
    A("")
    A("## 3. Result at a glance")
    A("")
    tab, rank = fb["sds_table"], fb["ranking_by_pooled_sds"]
    g2 = fb["gamma2"]
    g2_sign = g2["sign_decomposition"]
    A(f"- **{rank[0]} leads** the ranking of the five mask dimensions by pooled SDS "
      "(healthy vs. all damaged overpasses): "
      + ", ".join(f"{k} {fmt(tab[k]['pooled'], 2)}" for k in rank)
      + f". gamma^2 reaches only {fmt(tab['gamma2']['pooled'], 2)} for the same "
      "pooled comparison — the pooling cancels it, because its three per-state "
      "contrasts have different signs (Cliff's delta "
      + ", ".join(f"{k} {dsgn(g2_sign[k])}" for k in ("DAM3", "DAM4", "DAM6"))
      + ") and hence offset each other. Its best single-state effect "
      f"(SDS {fmt(g2['max_single_group']['sds'], 2)} against "
      f"{g2['max_single_group']['group']}) is of the same order as the mask "
      "channels, but restricted to the mildest state.")
    pair_p = [(pk, res["fig_d"]["pairs"][pk]["per_model"][core.HEADLINE_PAIR]["perm_p"])
              for pk in PAIR_KEYS]
    sig = [f"{pk.replace('_vs_', ' vs ')} (p={fmt(p, 3)})"
           for pk, p in pair_p if p is not None and p < 0.05]
    A(f"- **4-class LDA (6 channels, LOO, lambda="
      f"{fc['headline_lambda']:g})**: balanced accuracy "
      f"{fmt(fc['headline_balanced_accuracy_6d'])} against a chance level of "
      f"{fmt(fc['chance_balanced_accuracy'])} — no state recovery; the "
      f"label-permutation test gives p = "
      f"{fmt(fc['headline_permutation_test']['x_6d']['p_balanced_accuracy'])}.")
    mono = fe["monotonicity"]["hypothesis"]
    A(f"- **Severity must be monotone in abs(delta)**: "
      "abs(delta(healthy,DAM3)) < abs(delta(healthy,DAM4)) < abs(delta(healthy,DAM6)) "
      f"is **{'met' if mono.get('holds') else 'NOT met'}** for gamma^2 "
      f"(observed order {mono.get('order_observed', 'n/a')}).")
    A("- **The pairwise test that does hold** (out-of-fold scores, headline "
      "vector [gamma^2, P, D], permutation p < 0.05): "
      + (", ".join(sig) if sig else "none") + ".")
    A("")
    A("## 4. Figures")
    A("")
    A("### Figure A — raw distributions per state "
      "(`figures/fig_lumo_ocv_paper_A.png`)")
    A("")
    A("Dots are single overpasses, the bar is the median; the input is the "
      "overpass table of section 2, nothing else.")
    A("")
    A(md_table(["component, median"] + states, [
        [DIM_TITLE[key]] + [fmt(res["fig_a"]["descriptives"][key][s]["median"])
                            for s in states]
        for key in ["gamma2"] + core.DIMS]))
    A("")
    A("### Figure B — in-sample effect sizes "
      "(`figures/fig_lumo_ocv_paper_B.png`)")
    A("")
    A("SDS = 10 |Cliff's delta| (0 = identical distributions, 2 = large "
      "separation, 4 = complete separation); *pooled* compares healthy with all "
      "damaged overpasses together. Panel (a) ranks all six components by that "
      "pooled SDS, panel (b) shows the three per-state values.")
    A("")
    fam = {k: fb["per_dimension_full"][k]["sign_decomposition"] for k in core.DIMS}
    collapse = [k for k in core.DIMS if tab[k]["pooled"] < 0.5]
    stable = [k for k in core.DIMS if k not in collapse]
    A("Pooling is the sensitive part of the comparison and the reason why gamma^2 "
      "is *not* the leading component in panel (a): its three per-state contrasts "
      "have different signs (Cliff's delta "
      + ", ".join(f"{k} {dsgn(g2_sign[k])}" for k in ("DAM3", "DAM4", "DAM6"))
      + ") at comparable magnitude, so their group-size-weighted mean collapses to "
      f"{dsgn(g2['pooled_all_damaged']['delta'])} (SDS "
      f"{fmt(tab['gamma2']['pooled'], 2)}) — smaller than any of its three parts. "
      "Among the mask channels only "
      + ", ".join(f"`{k}` (deltas "
                  + " / ".join(dsgn(fam[k][s])
                               for s in ("DAM3", "DAM4", "DAM6")) + ")"
                  for k in collapse)
      + " collapses the same way; "
      + ", ".join(f"`{k}`" for k in stable)
      + " keep a pooled SDS of "
      f"{fmt(min(tab[k]['pooled'] for k in stable), 2)}-"
      f"{fmt(max(tab[k]['pooled'] for k in stable), 2)}: their DAM3 and DAM4 "
      "contrasts are negative in "
      f"{len([k for k in stable if fam[k]['DAM3'] < 0 and fam[k]['DAM4'] < 0])} "
      f"of {len(stable)} cases, so the dominant direction survives the pooling.")
    A("")
    A("The one-against-one view of panel (b) is therefore the honest reading of "
      f"gamma^2: it reaches SDS {fmt(g2['max_single_group']['sds'], 2)}, but only "
      f"against {g2['max_single_group']['group']}, the mildest state, and even "
      "there the sign flips between the orbits (ASC "
      f"{dsgn(res['fig_e']['by_orbit']['ASCENDING'][g2['max_single_group']['group']]['delta'])}, "
      "DESC "
      f"{dsgn(res['fig_e']['by_orbit']['DESCENDING'][g2['max_single_group']['group']]['delta'])}).")
    A("")
    A(md_table(["component", "SDS pooled", "SDS max single state",
                "SDS vs. DAM3", "SDS vs. DAM4", "SDS vs. DAM6"], [
        [DIM_TITLE[key], fmt(tab[key]["pooled"], 2),
         fmt(tab[key]["max_single_group"], 2)]
        + [fmt(tab[key]["healthy_vs"].get(k), 2) for k in ("DAM3", "DAM4", "DAM6")]
        for key in ["gamma2"] + core.DIMS]))
    A("")
    A("### Figure C — 4-class LDA, leave-one-out "
      "(`figures/fig_lumo_ocv_paper_C.png`)")
    A("")
    A("Balanced accuracy over the shrinkage grid (uniform priors, standardised "
      "features, identical fold construction as the committed pipeline):")
    A("")
    grid = fc["sensitivity_grid"]
    A(md_table(["feature set", f"lambda -> balanced accuracy"], [
        [core.MODEL_LABEL.get(k, k),
         ", ".join(f"{float(l):g}:{fmt(grid[k][l]['balanced_accuracy'], 3)}"
                   for l in map(str, fc["lambda_grid"]))]
        for k in ["x_6d", "gamma2", "P", "D", "A", "F", "S"] if k in grid]))
    A("")
    conf = fc["headline_confusion_matrix_6d"]
    cls = list(conf)
    A("Headline 6-channel model at lambda = "
      f"{fc['headline_lambda']:g}: accuracy {fmt(fc['headline_accuracy_6d'])}, "
      f"balanced accuracy {fmt(fc['headline_balanced_accuracy_6d'])}; "
      "confusion matrix (rows true, columns predicted):")
    A("")
    A(md_table(["true \\ predicted"] + [core.SHORT[c] for c in cls],
               [[core.SHORT[a]] + [conf[a][b] for b in cls] for a in cls]))
    A("")
    perm = fc["headline_permutation_test"]
    A(md_table(["model", "balanced acc. observed", "null mean", "null p95",
                "permutation p", "recomputed here"], [
        [core.MODEL_LABEL.get(k, k), fmt(v["balanced_accuracy_observed"]),
         fmt(v["balanced_accuracy_null_mean"]), fmt(v["balanced_accuracy_null_p95"]),
         fmt(v["p_balanced_accuracy"]), v.get("recomputed", False)]
        for k, v in perm.items()]))
    A("")
    A("### Figure D — out-of-fold scores of the six state pairs "
      "(`figures/fig_lumo_ocv_paper_D.png`)")
    A("")
    A("Scores come from a Fisher LDA (shrinkage "
      f"{core.HEADLINE_LAMBDA:g}) fitted without the test point; SDS is computed "
      "on these out-of-fold scores, and the permutation p is that of the "
      "headline model [gamma^2, P, D].")
    A("")
    prows = []
    for pk in PAIR_KEYS:
        pr = res["fig_d"]["pairs"][pk]
        prows.append([pk.replace("_vs_", " vs "),
                      f"{pr['n'][0]}/{pr['n'][1]}"]
                     + [fmt(pr["per_model"][m]["sds"], 2) for m in core.MODELS]
                     + [fmt(pr["per_model"][core.HEADLINE_PAIR]["perm_p"]),
                        fmt(pr["per_model"][core.HEADLINE_PAIR]["delta_ci95"][0], 2)
                        + " .. "
                        + fmt(pr["per_model"][core.HEADLINE_PAIR]["delta_ci95"][1], 2)])
    A(md_table(["pair", "n a/b"]
               + [SHORT_MODEL[m] for m in core.MODELS]
               + ["perm p (headline)", "delta CI95 (headline)"], prows))
    A("")
    A("Paired bootstrap of SDS(headline) - SDS(gamma^2 only) on the same "
      "resampled indices (positive = the core vector adds over coherence):")
    A("")
    A(md_table(["pair", "diff mean", "CI95", "CI excludes 0"], [
        [pk.replace("_vs_", " vs "),
         fmt(res["fig_d"]["pairs"][pk]["paired_diff_vs_gamma2"][core.HEADLINE_PAIR]["diff_mean"], 2),
         fmt(res["fig_d"]["pairs"][pk]["paired_diff_vs_gamma2"][core.HEADLINE_PAIR]["ci95"][0], 2)
         + " .. "
         + fmt(res["fig_d"]["pairs"][pk]["paired_diff_vs_gamma2"][core.HEADLINE_PAIR]["ci95"][1], 2),
         res["fig_d"]["pairs"][pk]["paired_diff_vs_gamma2"][core.HEADLINE_PAIR]["ci_excludes_zero"]]
        for pk in PAIR_KEYS]))
    A("")
    A("### Figure E — controls: severity, strata, mask brightness, wind "
      "(`figures/fig_lumo_ocv_paper_E.png`)")
    A("")
    m = fe["monotonicity"]
    hy, sp = m["hypothesis"], m.get("spearman_severity_vs_value") or {}
    jt = m["jonckheere"]
    A("Severity trend of gamma^2 (severity ranks "
      + ", ".join(f"{k}={v}" for k, v in jt["severity_ranks"].items()) + "):")
    A("")
    A(md_table(["test", "value"], [
        ["abs(delta) healthy vs. DAM3 / DAM4 / DAM6",
         " / ".join(fmt(hy["abs_delta"][k], 3) for k in ("DAM3", "DAM4", "DAM6"))],
        ["claim abs(d3) < abs(d4) < abs(d6)", hy.get("holds", "n/a")],
        ["observed order", hy.get("order_observed", "n/a")],
        ["Spearman rho(severity, gamma^2)",
         f"{fmt(sp.get('rho'))} (p = {fmt(sp.get('p'), 2)})"],
        ["Jonckheere-Terpstra p (increasing)", fmt(jt["p_increasing"])],
        ["Jonckheere-Terpstra p (decreasing)", fmt(jt["p_decreasing"])],
        ["permutation draws", jt["n_perm"]],
    ]))
    A("")
    A("Strata controls (SDS of gamma^2, healthy vs. state; `n/a` = fewer than 3 "
      "overpasses on one side of the comparison, i.e. the stratum cannot be "
      "evaluated):")
    A("")
    srows = []
    for name, blk in ([("all overpasses", None)]
                      + [(s, fe["by_season"][s]) for s in fe["by_season"]]
                      + [(s, fe["by_orbit"][s]) for s in fe["by_orbit"]]):
        row = [name]
        for k in ("DAM3", "DAM4", "DAM6"):
            if blk is None:
                row.append(fmt(tab["gamma2"]["healthy_vs"].get(k), 2))
            else:
                b = blk.get(k) or {}
                row.append(fmt(b.get("sds"), 2) if b.get("delta") is not None
                           else f"n/a ({b.get('note', 'n/a')})")
        srows.append(row)
    A(md_table(["stratum", "SDS vs. DAM3", "SDS vs. DAM4", "SDS vs. DAM6"], srows))
    A("")
    A("Covariate controls (overpass level; *windmatched* restricts both sides to "
      "wind speed <= 4 m/s):")
    A("")
    A(md_table(["control test", "n", "Welch p", "Mann-Whitney p", "Spearman"], [
        [k.replace("_", " "), v.get("n"), fmt(v.get("welch_p")),
         fmt(v.get("mannwhitney_p")),
         (f"rho={fmt(v.get('spearman_rho'))} p={fmt(v.get('p'), 2)}"
          if "spearman_rho" in v else "n/a")]
        for k, v in fe["coh"]["wind_control"].items()]))
    A("")
    A(md_table(["state", "median peak intensity", "median masked pixels A",
                "median wind speed [m/s]"], [
        [s, fmt(fe["peak_intensity_median_by_state"][s], 1),
         fmt(d["n_masked_median_by_state"][s], 1),
         fmt(fe["wind_speed_median_by_state"][s], 2)] for s in states]))
    A("")
    A("## 5. Verification against the committed LUMO results")
    A("")
    A("The package recomputes every number from scratch and compares it with the "
      "five committed reference JSONs in `data/reference/` (copies of the LUMO "
      "project outputs):")
    A("")
    A(md_table(["reference file", "pinned blocks"], [
        ["`" + core.REFERENCE_FILES[k] + "`", blk]
        for k, blk in [
            ("tower_coherence_states",
             "burst table, params, counts, per_state, tests, per_orbit, "
             "desc_tests, wind_control"),
            ("gamma2_pairs",
             "primary pairs block, monotonicity, season/orbit controls, mask "
             "controls, params"),
            ("echo_mask_vector",
             "per_dimension (gamma^2 + A/D/F/S/P), ranking, definitions, meta"),
            ("emv_multivariate_cv",
             "lambda grid, headline accuracies, confusion matrix, feature order, "
             "meta, permutation test"),
            ("emv_pairwise_triplets_cv",
             "models, shrinkage, n_perm, per-pair SDS / CI / permutation "
             "blocks")]]))
    A("")
    A(md_table(["pin result", "value"], [
        ["checks", pins["n_checks"]],
        ["exact matches", pins["n_exact"]],
        ["matches within tolerance", pins["n_within_tolerance"]],
        ["failures", pins["n_failed"]],
        ["max |deviation|", pins["max_abs_deviation"]],
        ["label strings translated (German -> English)", pins["label_translations"]],
        ["float tolerance", pins["float_tolerance"]],
        ["blocks only recomputed here",
         (f"{len(pins['recomputed_keys_not_in_reference'])} keys the pin block "
          "recomputes that the references do not store (e.g. "
          + ", ".join(pins["recomputed_keys_not_in_reference"][:3]) + ", ...)"
          if pins["recomputed_keys_not_in_reference"] else "none")],
    ]))
    A("")
    if pins["n_failed"]:
        A("Failing pins (first 10): "
          + ", ".join("`" + str(f["key"]) + "`" for f in pins["failures"][:10]))
        A("")
    A(f"Quick mode: **{'yes' if quick else 'no'}** — in quick mode the two "
      "label-permutation tests (long runtime) are not recomputed but copied from "
      "the committed references and marked as such in the tables above. All other "
      "numbers are always recomputed.")
    A("")
    A("## 6. Reproduce")
    A("")
    A("```bash")
    A("cd research/observability-core-vector")
    A("python3 code/lumo_ocv_channels_csv.py    # rebuild data/lumo_ocv_channels.csv")
    A("python3 code/fig_lumo_ocv_paper.py      # figures A-E, JSON and this report")
    A("bash figures/fig_lumo_ocv_paper.sh      # the same, via the shell wrapper")
    A("python3 code/fig_lumo_ocv_paper.py --quick   # skip the permutation tests")
    A("```")
    A("")
    A("Determinism: all random draws use fixed seeds "
      f"(RNG seed {core.RNG_SEED}, {core.N_BOOT} bootstrap draws, "
      f"{core.N_PERM_4CLASS} 4-class and {core.N_PERM_PAIR} pairwise permutation "
      "draws) and none of the outputs contains a timestamp, so repeated runs "
      "produce identical JSON, Markdown and PNG files.")
    A("")
    A("## 7. Caveats")
    A("")
    A("- The overpass table is a **cached** sample: the mask geometry comes from "
      + f"{d['n_overpasses']} coherence bursts of one LUMO cache directory, and "
      "one outlier overpass was excluded from the P computation. The package can "
      "recompute all statistics from that table but not re-run the "
      "interferometric processing.")
    A("- The states are **label intervals of one infrastructure section**, not "
      "independently surveyed plots, and n differs strongly per state ("
      + ", ".join(f"{s} n={d['n_by_state'][s]}" for s in states)
      + "): balanced accuracy and uniform priors are used because of this, and "
      "the per-state SDS values rest on small groups.")
    A("- Figure D scores are **out-of-fold but not spatially or temporally "
      "separated**: bursts of the same date appear in training and test, so the "
      "numbers answer \"is a state recoverable at all\", not \"does it generalise "
      "to unseen dates\".")
    A("- Wind and season strata thin out the sample per cell; strata with fewer "
      "than 3 overpasses on one side are reported as not evaluable instead of "
      "being silently dropped.")
    A("- P is defined through the masks of the same overpass set (mean pixel "
      "frequency inside the own mask) and is thus not independent of the mask "
      "threshold rule; P, A, D, F and S all share that one mask, so A/D/F/S "
      "control the *morphology* of the same mask rather than being independent "
      "measurements.")
    A("")
    return "\n".join(L) + "\n"


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
def json_default(o):
    """JSON fallback for the numpy scalars/arrays coming out of the statistics."""
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.floating):
        return float(o)
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    raise TypeError(f"not JSON serialisable: {type(o)}")


def build_json(res, pins_summary, quick):
    """Self-contained result JSON: everything the figures plot, plus the pins."""
    return {
        "meta": {
            "generator": "code/fig_lumo_ocv_paper.py",
            "quick": bool(quick),
            "csv": res["data"]["csv"],
            "csv_sha256": res["data"]["csv_sha256"],
            "n_overpasses": res["data"]["n_overpasses"],
            "figures": ["figures/fig_lumo_ocv_paper_%s.png" % k for k in "ABCDE"],
            "report": "figures/fig_lumo_ocv_paper.md",
            "versions": {
                "python": platform.python_version(),
                "numpy": np.__version__,
                "matplotlib": matplotlib.__version__,
                "scipy": getattr(scipy, "__version__", None) if scipy else None,
            },
            "seeds": {"rng_seed": core.RNG_SEED, "n_boot": core.N_BOOT,
                      "n_perm_4class": core.N_PERM_4CLASS,
                      "n_perm_pair": core.N_PERM_PAIR},
            "note": "no timestamps on purpose: repeated runs must be identical",
        },
        "data": res["data"],
        "fig_a": res["fig_a"],
        "fig_b": res["fig_b"],
        "fig_c": res["fig_c"],
        "fig_d": res["fig_d"],
        "fig_e": res["fig_e"],
        "pins": pins_summary,
    }


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Recompute figures A-E of the OCV observability package.")
    ap.add_argument("--quick", action="store_true",
                    help="copy the two long label-permutation tests from the "
                         "committed references instead of recomputing them")
    ap.add_argument("--figdir", default=FIGDIR,
                    help="output directory of the PNG figures")
    ap.add_argument("--json", default=OUT_JSON, help="path of the result JSON")
    ap.add_argument("--md", default=OUT_MD, help="path of the Markdown report")
    args = ap.parse_args(argv)

    t0 = time.perf_counter()
    print("OCV figure package — recomputing figures A-E "
          f"({'quick' if args.quick else 'full'} mode) ...")
    res, _X, _y, _complete = compute(quick=args.quick)
    print(f"  compute: {time.perf_counter() - t0:.1f}s")

    P = pin_all(res, quick=args.quick)
    s = P.summary()
    print(f"  pins: {s['n_checks']} checks, {s['n_exact']} exact, "
          f"{s['n_failed']} failures, max |dev| {s['max_abs_deviation']}, "
          f"labels translated {s['label_translations']}")

    os.makedirs(args.figdir, exist_ok=True)
    print("  writing figures:")
    for fn in (fig_a, fig_b, fig_c, fig_d, fig_e):
        fn(res, args.figdir)

    with open(args.json, "w", encoding="utf-8") as fh:
        json.dump(build_json(res, s, args.quick), fh, indent=1, sort_keys=True,
                  default=json_default)
        fh.write("\n")
    print(f"  {os.path.relpath(args.json, HERE)}")

    with open(args.md, "w", encoding="utf-8") as fh:
        fh.write(report(res, s, args.quick))
    print(f"  {os.path.relpath(args.md, HERE)}")

    print(f"done in {time.perf_counter() - t0:.1f}s, "
          f"{s['n_failed']} pin failures")
    return 0 if s["n_failed"] == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())

