#!/usr/bin/env python3
"""Carola OCV-Paper — figures A–E, JSON result and English report.

Question: does the observability core vector
``OCV = [gamma2, P, D]`` of the Carola bridge (Dresden) — the interferometric
coherence gamma^2 of its bridge-deck echo, plus the persistence ``P`` and the
density ``D`` of that echo mask — separate the two states of the site, the
bridge **before** the 2024-09-11 collapse (``pre-collapse (healthy)``) and the
bridge **after** it (``post-collapse``)? And does the *per-girder deck mask* —
the same six dimensions on Girder A / B / C only, the layer this package adds —
carry the contrast, or does the whole-bridge window alone?

The full mask vector is the 6D ``x = [gamma2, P, D, A, F, S]`` (``A`` = masked
pixels, ``F`` = 8-connected components, ``S`` = ||mask centroid - peak||). Both
mask layers are recomputed here from the committed window cache
``data/carola/carola_windows_mask_cache.txt`` with the mask rule of the analysis
project (``0.30 x peak`` and ``5 x np.median``, ``min_n_masked = 2``); see
``code/carola/carola_ocv_masks.py``.

Approach: the committed 1,717-chip channel table (``carola_ocv_channels.csv``,
built by ``carola_ocv_channels_csv.py`` from the cache, the committed
measurement extract and the segments file), the same state/collapse definitions
as the Carola project, and ``carola_ocv_stats.py`` as a verbatim port of the
LUMO statistics plus a port of the site's own correlation/regression battery.
Every number is recomputed here and then pinned against the committed reference
JSONs; a deviation above the tolerance aborts the script with exit code != 0.

What this package claims, and what it does not
----------------------------------------------
  * **two states on one bridge.** ``pre-collapse (healthy)`` and
    ``post-collapse`` are the two sides of the collapse date of the *same*
    structure, so season, weather, traffic and pipeline generation all move with
    the state. Every contrast is therefore reported beside a co-variate control
    (weather, orbit, girder, mask brightness) — never as damage.
  * **two mask layers.** The echo mask is the whole chip; the deck mask is the
    girder's own chip. ``A`` of neither layer is the pipeline's own
    ``coherence_masked_pixels``: the committed export holds the per-girder
    re-extraction while that column is the pre-fix asset-point chip, which is
    why the site's analysis recomputed the mask from the payloads (and this
    package re-derives the mask from the committed cache).
  * **the mask layer rests on a committed cache.** The 182 MB chip payloads do
    not belong in the repository; the cache carries the echo mask of every chip
    plus the sufficient statistics, and the build re-verifies the mask rule on
    it (``--check``). The pins are the rest of the argument.

Pinned reference files (``data/carola/reference/``):
  carola_coherence_states.json    per-window gamma2 / n_masked, reproduced
                                  exactly for the 925 of its 1,454 windows that
                                  the current export still holds
  carola_echo_mask_gamma2.json    state table, segment table, SDS and nested
                                  OLS — reproduced number by number
  carola_bridge_osm.json          the committed bridge geometry, pinned
                                  structurally (a map extract, not rederivable)

Figures:
  A  raw distributions of the six channels per state (echo mask, all 1,717 chips)
  B  in-sample effect sizes: Cliff's delta + bootstrap CI and SDS for the six
     channels, per state and per girder
  C  2-class LDA (pre vs. post), leave-one-out, lambda grid, confusion matrix,
     label permutation — plus the weather-only competitor
  D  out-of-fold scores of the five models and the pairwise SDS contrast
  E  the evidence against reading the contrast as damage: girder / orbit /
     season strata, the mask-size co-variate and the reference verdicts

Usage:
  python3 code/carola/fig_carola_ocv_paper.py            # recompute everything
  python3 code/carola/fig_carola_ocv_paper.py --quick    # no permutation null
"""
from __future__ import annotations

import argparse
import json
import os
import pickle
import platform
import sys
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import carola_ocv_core as core   # noqa: E402
import carola_ocv_masks as masks  # noqa: E402
import carola_ocv_stats as st    # noqa: E402

FIGDIR = core.FIGDIR
OUT_JSON = os.path.join(core.DATA, "fig_carola_ocv_paper.json")
OUT_MD = os.path.join(FIGDIR, "fig_carola_ocv_paper.md")

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

# The collapse alphabet of the site, short for the figures: the Carola project
# calls the two states "healthy" and "damaged"; this package prefers the
# chronological names and keeps the mapping in one place.
RKEY = {core.PRE: "healthy", core.POST: "damaged"}
SHORT = {core.PRE: "pre\n(healthy)", core.POST: "post\n(collapse)"}
STATE_COLORS = {core.PRE: "#4C72B0", core.POST: "#C44E52"}
GIRDER_COLORS = {"Girder A": "#4C72B0", "Girder B": "#55A868", "Girder C": "#C44E52"}
CHANNEL_COLORS = {"gamma2": "#4C72B0", "P": "#55A868", "D": "#8172B3",
                  "A": "#8172B3", "F": "#64B5CD", "S": "#DA8BC3",
                  "wind_speed_ms": "#64B5CD", "temperature_c": "#C44E52",
                  "brightness_ratio": "#937860", "measured_frequency_hz": "#8C8C8C",
                  "coherence": "#DA8BC3"}
MODEL_COLORS = {"gamma2": "#8C8C8C", "PD": "#DD8452", "gamma2PD": "#4C72B0",
                "x_6d": "#55A868", "weather": "#8172B3"}
MODEL_MARKERS = {"gamma2": "o", "PD": "s", "gamma2PD": "^", "x_6d": "D",
                 "weather": "v"}
SHORT_MODEL = {"gamma2": "gamma2", "PD": "P+D", "gamma2PD": "gamma2+P+D",
               "x_6d": "6 dims", "weather": "weather"}
CHANCE = 0.5
# The nine repair targets of the 429 damage (request-id prefixes) — the
# reference's ``REPAIR_TARGETS``, used only to reproduce its per-segment flag.
REPAIR_TARGETS = ["90a3e881", "ccd664a3", "997ecabe", "4ec6636a", "cc31d61b",
                  "1e1fce7e", "031a8e5b", "b8c52dfc", "c45e1e4c"]
# The committed cross-reference of the analysis project (``COMMITTED`` in
# ``analyze_carola_echo_mask_gamma2.py``, taken there from
# ``carola_coherence_states.md``, the old 1,454-window basis). The pin block
# below also re-derives it from the committed states JSON.
COMMITTED_REFERENCE = {
    "source": "carola_coherence_states.md (1454 stored windows)",
    "healthy": {"n": 880, "n_masked_median": 13.0, "gamma2_median": 0.1228},
    "damaged": {"n": 567, "n_masked_median": 6.0, "gamma2_median": 0.1728},
}
CHANNELS = list(core.FEATURES) + ["wind_speed_ms", "temperature_c",
                                  "brightness_ratio", "measured_frequency_hz",
                                  "coherence"]


def save_fig(fig, name, out_dir=None):
    """Write one 600 dpi PNG (the binaries stay out of git)."""
    out_dir = out_dir or FIGDIR
    png = os.path.join(out_dir, name + ".png")
    fig.savefig(png)
    plt.close(fig)
    print(f"  {os.path.relpath(png)}")
    return png


# ---------------------------------------------------------------------------
# Row views
# ---------------------------------------------------------------------------
def mask_rows(csv_rows):
    """The rows as this package analyses them (full-precision gamma2)."""
    out = []
    for r in csv_rows:
        d = core.row_brief(r)
        d["acquisition_ts"] = r.get("acquisition_ts")
        d["orbit_direction"] = r.get("orbit_direction")
        d["request_id"] = r.get("request_id") or ""
        d["month"] = r.get("month")
        d["day"] = core.day_of(r)
        d["echo_mode"] = r.get("echo_mode")
        d["A_export"] = core._num(r.get("coherence_masked_pixels"))
        out.append(d)
    return out


def ref_rows(csv_rows):
    """The same rows in the *analysis project's* convention.

    ``analyze_carola_echo_mask_gamma2.mask_metrics`` rounds gamma2 to four
    decimals, and every aggregate block of its committed JSON is computed from
    the rounded value — so the pinned numbers must be too.
    """
    out = []
    for r in csv_rows:
        day = core.day_of(r)
        d = {
            "gamma2": round(core._num(r["gamma2_mask"]), 4),
            "coherence_masked_pixels": core._num(r["A"]),
            "asset": r.get("asset_name"),
            "segment": r.get("segment"),
            "date": day,
            "orbit_direction": r.get("orbit_direction"),
            "state": "damaged" if day >= core.COLLAPSE_DATE else "healthy",
            "acquisition_ts": r.get("acquisition_ts"),
            "request_id": r.get("request_id") or "",
            # the analysis project groups by the 8-char prefix (``request_short``)
            "request_short": (r.get("request_id") or "?")[:8],
        }
        d["repair_target"] = any(str(d["request_id"]).startswith(p)
                                 for p in REPAIR_TARGETS)
        out.append(d)
    return out


# ---------------------------------------------------------------------------
# Compute
# ---------------------------------------------------------------------------
def _med(v):
    v = [x for x in v if x is not None]
    return float(np.median(v)) if v else None


def _feature_matrix(rows, keys):
    return np.array([[r.get(k) if r.get(k) is not None else 0.0 for k in keys]
                     for r in rows], dtype=float)


def state_table(rows):
    """The reference's ``state_table``: healthy / damaged / pooled / by_orbit."""
    out = {}
    for state in ("healthy", "damaged"):
        sel = [r for r in rows if r["state"] == state]
        masks = [r["coherence_masked_pixels"] for r in sel]
        sb = st.pair_block(sel, "coherence_masked_pixels", "gamma2")
        out[state] = {
            "n": len(sel),
            "n_masked_median": _med(masks),
            "n_masked_range": [min(masks), max(masks)] if masks else None,
            "gamma2_median": _med([r["gamma2"] for r in sel]),
            "rho_gamma2_vs_n_masked": sb,
            "date_range": ([min(r["date"] for r in sel), max(r["date"] for r in sel)]
                           if sel else None),
        }
    out["pooled"] = st.pair_block(rows, "coherence_masked_pixels", "gamma2")
    out["by_orbit"] = {o: st.pair_block([r for r in rows
                                         if r["orbit_direction"] == o],
                                        "coherence_masked_pixels", "gamma2")
                       for o in ("ASCENDING", "DESCENDING")}
    out["committed_reference"] = COMMITTED_REFERENCE
    return out


def group_table(rows, key):
    """The reference's ``group_table``: one block per group value."""
    groups = {}
    for r in rows:
        groups.setdefault(r.get(key) or "?", []).append(r)
    out = {}
    for k, recs in sorted(groups.items(), key=lambda kv: -len(kv[1])):
        masks = [r["coherence_masked_pixels"] for r in recs]
        sb = st.pair_block(recs, "coherence_masked_pixels", "gamma2")
        out[str(k)] = {
            "n": len(recs), "rho": sb["rho"], "perm_p": sb["perm_p"],
            "gamma2_median": _med([r["gamma2"] for r in recs]),
            "n_masked_median": _med(masks),
            "n_masked_range": [min(masks), max(masks)] if masks else None,
            "n_healthy": sum(1 for r in recs if r["state"] == "healthy"),
            "n_damaged": sum(1 for r in recs if r["state"] == "damaged"),
            "repair_target": any(r["repair_target"] for r in recs),
        }
    return out


def sds_block(rows):
    """The reference's ``sds_block``: Cliff's delta + SDS + scenario."""
    dam = [r for r in rows if r["state"] == "damaged"]
    hea = [r for r in rows if r["state"] == "healthy"]
    out = {}
    for key, label in (("gamma2", "gamma2"), ("coherence_masked_pixels", "n_masked")):
        p, delta = st.mannwhitney([r[key] for r in dam], [r[key] for r in hea])
        out[label] = {"n_damaged": len(dam), "n_healthy": len(hea),
                      "median_damaged": _med([r[key] for r in dam]),
                      "median_healthy": _med([r[key] for r in hea]),
                      "cliff_delta": delta, "mannwhitney_p": p,
                      "sds": None if delta is None else 10.0 * abs(delta)}
    sg, sn = out["gamma2"]["sds"], out["n_masked"]["sds"]
    out["ratio_sds_n_masked_over_gamma2"] = (sn / sg) if (sg and sn and sg > 0) else None
    if sg is None or sn is None:
        scen = "not_testable"
    elif sn < 0.1:
        scen = "A"
    elif sn < sg:
        scen = "B"
    else:
        scen = "C"
    out["scenario"] = scen
    # The reading is kept verbatim from the analysis project (it is its wording,
    # German, and is part of the pinned block).
    out["scenario_reading"] = {
        "A": "SDS(n_masked) ~ 0 -> kein klarer Schadensbezug der Maske",
        "B": "SDS(n_masked) < SDS(gamma2) -> Schaden -> n_masked -> gamma2 (indirekter Kanal)",
        "C": "SDS(n_masked) >= SDS(gamma2) -> Maske trennt die Zustaende mindestens so gut wie gamma2",
        "not_testable": "zu wenige Werte je Zustand",
    }[scen]
    return out


def model_block(rows):
    """The reference's ``model_block``: the nested OLS battery + inverse-N."""
    y = np.array([r["gamma2"] for r in rows], float)
    mask = np.array([r["coherence_masked_pixels"] for r in rows], float)
    dm = np.array([1.0 if r["state"] == "damaged" else 0.0 for r in rows])
    inv = 1.0 / np.maximum(mask, 1.0)
    ones = np.ones(len(rows))
    m1 = st.ols_fit(y, np.column_stack([ones, dm]), ["(intercept=healthy)", "damaged"])
    m2 = st.ols_fit(y, np.column_stack([ones, mask]), ["(intercept)", "n_masked"])
    m3 = st.ols_fit(y, np.column_stack([ones, dm, mask]),
                    ["(intercept=healthy)", "damaged", "n_masked"])
    return {
        "model1_state_only": m1, "model2_mask_only": m2,
        "model3_state_plus_mask": m3,
        "inverse_n_fit": st.ols_fit(y, np.column_stack([ones, inv]),
                                    ["(intercept)", "1/n_masked"]),
        "f_test_add_n_masked_to_state": st.f_test(m1["rss"], m3["rss"], 1, m3["df_resid"]),
        "f_test_add_state_to_n_masked": st.f_test(m2["rss"], m3["rss"], 1, m3["df_resid"]),
        "delta_r2_from_adding_n_masked": m3["r2"] - m1["r2"],
        "delta_r2_from_adding_state": m3["r2"] - m2["r2"],
    }


def channel_deltas(rows):
    """Cliff's delta + bootstrap CI + SDS of the six channels (pre vs. post)."""
    out = {}
    for k in core.FEATURES:
        a = [r[k] for r in rows if r["state"] == core.PRE and r.get(k) is not None]
        b = [r[k] for r in rows if r["state"] == core.POST and r.get(k) is not None]
        if not a or not b:
            out[k] = {"n": [len(a), len(b)]}
            continue
        out[k] = st.delta_block(a, b, f"pre_vs_post_{k}", n_boot=core.N_BOOT,
                                seed=core.RNG_SEED)
    return out


def girder_deltas(rows):
    """The same effect sizes inside each girder (the deck-mask layer)."""
    out = {}
    for g, recs in core.by_girder(rows).items():
        if not recs:
            continue
        block = {}
        for k in core.FEATURES:
            a = [r[k] for r in recs if r["state"] == core.PRE and r.get(k) is not None]
            b = [r[k] for r in recs if r["state"] == core.POST and r.get(k) is not None]
            if not a or not b:
                block[k] = {"n": [len(a), len(b)]}
                continue
            block[k] = st.delta_block(a, b, f"{g}_{k}", n_boot=core.N_BOOT,
                                      seed=core.RNG_SEED)
        out[g] = {"n_rows": len(recs),
                  "n_pre": sum(1 for r in recs if r["state"] == core.PRE),
                  "n_post": sum(1 for r in recs if r["state"] == core.POST),
                  "channels": block}
    return out


def lda_block(rows, quick=False):
    """2-class LDA of every model: LOO over the lambda grid + permutation null."""
    X6 = _feature_matrix(rows, core.FEATURES)
    y = np.array([r["state"] for r in rows], dtype=object)
    classes = list(core.STATES)
    grid = {}
    for name, keys in core.MODELS.items():
        idx = [core.FEATURES.index(k) for k in keys if k in core.FEATURES]
        if not idx:
            continue
        X = X6[:, idx]
        per_lambda = {}
        for lam in core.LAMBDA_GRID:
            _, acc, bal, conf = st.loo_cv(X, y, classes, lam)
            per_lambda[str(lam)] = {"accuracy": acc, "balanced_accuracy": bal,
                                    "confusion": conf}
        grid[name] = per_lambda
    out = {"lambda_grid": list(core.LAMBDA_GRID),
           "headline_lambda": core.HEADLINE_LAMBDA, "grid": grid,
           "headline": grid.get(core.HEADLINE_PAIR, {}).get(
               str(core.HEADLINE_LAMBDA))}
    idx = [core.FEATURES.index(k) for k in core.MODELS[core.HEADLINE_PAIR]]
    X = X6[:, idx]
    if quick:
        out["headline_permutation_test"] = {
            "recomputed": False,
            "note": "run without --quick for the label permutation"}
    else:
        out["headline_permutation_test"] = st.permutation_test(
            X, y, classes, core.HEADLINE_LAMBDA, n_perm=core.N_PERM,
            seed=core.RNG_SEED)
    return out


def oof_block(rows):
    """OOF Fisher-direction SDS per model — the pairwise CV of the package."""
    X6 = _feature_matrix(rows, core.FEATURES)
    y = np.array([r["state"] for r in rows], dtype=object)
    classes = list(core.STATES)
    q, per_model = {}, {}
    for name, keys in core.MODELS.items():
        idx = [core.FEATURES.index(k) for k in keys if k in core.FEATURES]
        if not idx:
            continue
        X = X6[:, idx]
        qt = st.oof_scores(X, y, classes[0], classes[1], core.HEADLINE_LAMBDA)
        d, sds, ref, grp = st.delta_sds(qt, y, classes[0], classes[1])
        q[name] = [float(v) for v in qt]
        per_model[name] = {"n": [len(ref), len(grp)], "delta": d, "sds": sds}
    rng = np.random.default_rng(core.RNG_SEED)
    q_base = np.asarray(q["gamma2"], float)
    for name in q:
        if name == "gamma2":
            continue
        per_model[name]["paired_vs_gamma2"] = st.paired_bootstrap_diff(
            np.asarray(q[name], float), q_base, y, classes[0], classes[1],
            core.N_BOOT, rng)
    return {"q_by_model": q, "per_model": per_model,
            "headline_model": core.HEADLINE_PAIR}


def strata_block(rows, key_fn, key="gamma2"):
    """One ``delta_block`` per stratum of one stratifier."""
    out = {}
    for skey in sorted({key_fn(r) for r in rows if key_fn(r) is not None},
                       key=lambda s: str(s)):
        sub = [r for r in rows if key_fn(r) == skey]
        a = [r[key] for r in sub if r["state"] == core.PRE and r.get(key) is not None]
        b = [r[key] for r in sub if r["state"] == core.POST and r.get(key) is not None]
        if len(a) < 3 or len(b) < 3:
            out[str(skey)] = {"n": [len(a), len(b)], "delta": None, "note": "n < 3"}
        else:
            out[str(skey)] = st.delta_block(a, b, str(skey),
                                            n_boot=st.BOOT_STRATA,
                                            seed=core.RNG_SEED)
    return out


def compute(quick=False):
    """Recompute every number of the five figures and every pinned block."""
    csv_rows = core.load_csv()
    rows = mask_rows(csv_rows)
    rrows = ref_rows(csv_rows)
    t0 = time.time()
    res = {"rows": rows, "ref_rows": rrows, "csv_rows": csv_rows,
           "meta": core.load_csv_meta(),
           "ref_states": core.load_reference(core.REF_STATES),
           "ref_mask": core.load_reference(core.REF_MASK),
           "ref_osm": core.load_reference(core.REF_OSM),
           "manifest": core.load_manifest()}
    # the reference-convention aggregates -------------------------------------
    res["states"] = state_table(rrows)
    res["by_segment"] = group_table(rrows, "segment")
    res["by_request"] = group_table(rrows, "request_short")
    res["sds"] = sds_block(rrows)
    res["models"] = model_block(rrows)
    # this package's own figures ----------------------------------------------
    res["channels"] = channel_deltas(rows)
    res["girders"] = girder_deltas(rows)
    res["strata_girder"] = strata_block(rows, lambda r: r.get("segment"))
    res["strata_orbit"] = strata_block(rows, lambda r: r.get("orbit_direction"))
    res["strata_season"] = strata_block(rows, lambda r: st.season_of(r))
    res["lda"] = lda_block(rows, quick=quick)
    res["oof"] = oof_block(rows)
    res["mask_clutter"] = {
        "pooled": st.pair_block(rrows, "coherence_masked_pixels", "gamma2"),
        "ASCENDING": st.pair_block(rrows, "coherence_masked_pixels", "gamma2",
                                   orbit="ASCENDING"),
        "DESCENDING": st.pair_block(rrows, "coherence_masked_pixels", "gamma2",
                                    orbit="DESCENDING")}
    res["summary_by_state"] = core.summary_by_state(csv_rows)
    # the per-window pin counters (reported in the Markdown, pinned below)
    by_mid = {r["id"]: r for r in rows}
    checked = matched = 0
    for rec in res["ref_states"]["records"]:
        r = by_mid.get(rec["mid"])
        if r is None:
            continue
        checked += 1
        if (abs(float(rec["gamma2"]) - float(r["gamma2"])) < 1e-12
                and int(rec["n_masked"]) == int(r["A"])):
            matched += 1
    res["pins_windows"] = {"checked": checked, "matched": matched,
                           "of": len(res["ref_states"]["records"])}
    res["n_rows"] = len(rows)
    res["windows"] = {"n_rows": len(rows),
                      "n_segments": len(res["by_segment"]),
                      "n_requests": len(res["by_request"])}
    res["runtime_s"] = time.time() - t0
    return res


# ---------------------------------------------------------------------------
# Panels
# ---------------------------------------------------------------------------
def panel(ax, letter, text):
    ax.set_title(f"({letter}) {text}", loc="left")


def jitter(n, width, seed):
    rng = np.random.default_rng(seed)
    return rng.uniform(-width, width, n)


def strip(ax, rows, key, seed0=100, span=0.16, fs=6.4):
    """Two jittered strips (pre/post) with median markers."""
    for i, state in enumerate(core.STATES):
        vals = [r[key] for r in rows if r["state"] == state and r.get(key) is not None]
        if not vals:
            continue
        x = i + jitter(len(vals), span, seed0 + i)
        ax.plot(x, vals, "o", ms=1.6, alpha=0.35, color=STATE_COLORS[state],
                markeredgewidth=0)
        ax.plot([i - 0.3, i + 0.3], [np.median(vals)] * 2, "-", lw=1.8,
                color=STATE_COLORS[state])
        ax.text(i, ax.get_ylim()[1], f"n={len(vals)}", ha="center", va="top",
                fontsize=fs, color=STATE_COLORS[state])
    ax.set_xticks(range(len(core.STATES)))
    ax.set_xticklabels([SHORT[s] for s in core.STATES])


def fig_a(res, out_dir=None):
    rows = res["rows"]
    fig, axes = plt.subplots(2, 3, figsize=(COL_WIDTH, 105 * MM))
    for ax, k in zip(axes.ravel(), core.FEATURES[:6]):
        strip(ax, rows, k, seed0=100 + 7 * core.FEATURES.index(k))
        ax.set_title(core.CHANNEL_AXIS[k], loc="left", fontsize=8.0)
    fig.suptitle(f"(A) Echo-mask channels per state — {res['n_rows']} 80x80 chips, "
                 f"{core.SITE}", fontsize=9.0, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.965))
    return save_fig(fig, "fig_carola_ocv_paper_A", out_dir)


def fig_b(res, out_dir=None):
    ch = res["channels"]
    fig, axes = plt.subplots(1, 3, figsize=(COL_WIDTH, 78 * MM))
    keys = [k for k in core.FEATURES if "delta" in ch.get(k, {})]
    ax = axes[0]
    ax.axvline(0.0, color="#888888", lw=0.8)
    for i, k in enumerate(keys):
        b = ch[k]
        ax.plot([b["delta_ci95"][0], b["delta_ci95"][1]], [i, i], "-", lw=1.2,
                color=CHANNEL_COLORS[k])
        ax.plot(b["delta"], i, "o", ms=4.0, color=CHANNEL_COLORS[k])
    ax.set_yticks(np.arange(len(keys)))
    ax.set_yticklabels(keys)
    ax.invert_yaxis()
    ax.set_xlabel("Cliff's delta  (pre-collapse -> post-collapse)")
    panel(ax, "a", "echo mask, all chips")
    ax = axes[1]
    sds = {k: 10.0 * abs(ch[k]["delta"]) for k in keys}
    ax.barh(np.arange(len(keys)), [sds[k] for k in keys],
            color=[CHANNEL_COLORS[k] for k in keys], height=0.6)
    ax.set_yticks(np.arange(len(keys)))
    ax.set_yticklabels(keys)
    ax.invert_yaxis()
    ax.set_xlabel("SDS = 10|delta|")
    for i, k in enumerate(keys):
        ax.text(sds[k], i, f" {sds[k]:.2f}", va="center", fontsize=6.6)
    panel(ax, "b", "in-sample effect size")
    ax = axes[2]
    gd = res["girders"]
    gir = [g for g in core.GIRDERS if g in gd]
    w = 0.26
    for j, k in enumerate(["gamma2", "A", "D"]):
        vals = [gd[g]["channels"].get(k, {}).get("delta") or 0.0 for g in gir]
        ax.bar(np.arange(len(gir)) + (j - 1) * w, vals, width=w,
               color=CHANNEL_COLORS[k], label=k)
    ax.axhline(0.0, color="#888888", lw=0.8)
    ax.set_xticks(range(len(gir)))
    ax.set_xticklabels([g.replace("Girder ", "") for g in gir])
    ax.set_xlabel("girder")
    ax.set_ylabel("Cliff's delta")
    ax.legend(frameon=False, ncol=3)
    panel(ax, "c", f"per-girder deck mask (n={[gd[g]['n_rows'] for g in gir]})")
    fig.suptitle("(B) Effect sizes of the echo-mask and deck-mask channels "
                 "(Cliff's delta > 0: post-collapse larger)",
                 fontsize=9.0, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    return save_fig(fig, "fig_carola_ocv_paper_B", out_dir)


def fig_c(res, out_dir=None):
    lda = res["lda"]
    fig, axes = plt.subplots(1, 3, figsize=(COL_WIDTH, 80 * MM))
    ax = axes[0]
    lams = [float(v) for v in lda["lambda_grid"]]
    for name in core.MODELS:
        if name not in lda["grid"]:
            continue
        accs = [lda["grid"][name][str(l)]["accuracy"] for l in lams]
        ax.plot(lams, accs, "-" + MODEL_MARKERS.get(name, "o"), ms=3.2, lw=1.1,
                color=MODEL_COLORS.get(name, "#333333"),
                label=SHORT_MODEL.get(name, name))
    ax.axhline(CHANCE, color="#888888", ls=":", lw=0.9)
    ax.set_xlabel("shrinkage lambda")
    ax.set_ylabel("LOO accuracy")
    ax.legend(frameon=False, ncol=2)
    panel(ax, "a", "2-class LDA, leave-one-out")

    ax = axes[1]
    head = lda.get("headline") or {}
    conf = head.get("confusion") or {}
    classes = list(core.STATES)
    M = np.array([[conf.get(a, {}).get(b, 0) for b in classes] for a in classes],
                 dtype=float)
    ax.imshow(M, cmap="Blues", vmin=0)
    for i in range(len(classes)):
        for j in range(len(classes)):
            ax.text(j, i, int(M[i, j]), ha="center", va="center", fontsize=8)
    ax.set_xticks(range(len(classes)))
    ax.set_xticklabels([core.STATE_SHORT[s] for s in classes])
    ax.set_yticks(range(len(classes)))
    ax.set_yticklabels([core.STATE_SHORT[s] for s in classes])
    ax.set_xlabel("predicted")
    ax.set_ylabel("true")
    acc, bal = head.get("accuracy"), head.get("balanced_accuracy")
    panel(ax, "b", f"{core.HEADLINE_PAIR}, lambda={lda['headline_lambda']} "
                   f"(acc={acc:.3f}, bal={bal:.3f})")

    ax = axes[2]
    pt = lda.get("headline_permutation_test") or {}
    if pt.get("recomputed", True) and "p_accuracy" in pt:
        labs = ["accuracy", "balanced"]
        obs = [pt["accuracy_observed"], pt["balanced_accuracy_observed"]]
        nul = [pt["accuracy_null_mean"], pt["balanced_accuracy_null_mean"]]
        p95 = [pt["accuracy_null_p95"], pt["balanced_accuracy_null_p95"]]
        x = np.arange(len(labs))
        ax.bar(x - 0.18, obs, width=0.34, color="#4C72B0", label="observed")
        ax.bar(x + 0.18, nul, width=0.34, color="#CCCCCC", label="null mean")
        ax.plot(x + 0.18, p95, "k_", ms=12, label="null p95")
        ax.axhline(CHANCE, color="#888888", ls=":", lw=0.9)
        ax.set_xticks(x)
        ax.set_xticklabels(labs)
        ax.legend(frameon=False)
        ax.text(0.02, 0.02, f"p_acc={pt['p_accuracy']:.4f}\n"
                            f"p_bal={pt['p_balanced_accuracy']:.4f}\n"
                            f"n_perm={pt['n_perm']}", transform=ax.transAxes,
                fontsize=6.6, va="bottom")
    else:
        ax.text(0.5, 0.5, "label permutation skipped (--quick)",
                ha="center", va="center", transform=ax.transAxes, fontsize=8)
        ax.set_axis_off()
    panel(ax, "c", "label permutation null")
    fig.suptitle("(C) Is the state recoverable from the mask vector?",
                 fontsize=9.0, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    return save_fig(fig, "fig_carola_ocv_paper_C", out_dir)


def fig_d(res, out_dir=None):
    oof = res["oof"]
    classes = list(core.STATES)
    fig, axes = plt.subplots(1, 2, figsize=(COL_WIDTH, 80 * MM))
    ax = axes[0]
    names = [n for n in core.MODELS if n in oof["per_model"]]
    for i, name in enumerate(names):
        q = np.asarray(oof["q_by_model"][name], float)
        for st in classes:
            sel = np.asarray([r["state"] == st for r in res["rows"]])
            v = q[sel]
            if v.size == 0:
                continue
            ax.plot(i + jitter(v.size, 0.18, 500 + i),
                    v, "o", ms=1.7, alpha=0.35, color=STATE_COLORS[st],
                    markeredgewidth=0)
            ax.plot([i - 0.32, i + 0.32], [np.median(v)] * 2, "-", lw=1.6,
                    color=STATE_COLORS[st])
    ax.axhline(0.0, color="#888888", lw=0.8)
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels([SHORT_MODEL.get(n, n) for n in names])
    ax.set_ylabel("out-of-fold Fisher score q")
    panel(ax, "a", "OOF scores by model")

    ax = axes[1]
    ax.axvline(0.0, color="#888888", lw=0.8)
    for i, name in enumerate(names):
        b = oof["per_model"][name]
        if b.get("sds") is None:
            continue
        lo = b.get("paired_vs_gamma2", {}).get("ci95", [0.0, 0.0])
        ax.plot(lo, [i, i], "-", lw=1.1, color=MODEL_COLORS.get(name, "#333333"),
                alpha=0.7)
        ax.plot(b["sds"], i, MODEL_MARKERS.get(name, "o"), ms=5.0,
                color=MODEL_COLORS.get(name, "#333333"))
        ax.text(b["sds"], i, f"  {b['sds']:.2f}", va="center", fontsize=6.6)
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels([SHORT_MODEL.get(n, n) for n in names])
    ax.invert_yaxis()
    ax.set_xlabel("SDS of the OOF score (bars: 95% paired bootstrap CI vs gamma2)")
    panel(ax, "b", "out-of-fold separation")
    fig.suptitle("(D) Does the 6D vector separate better than gamma2 alone?",
                 fontsize=9.0, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    return save_fig(fig, "fig_carola_ocv_paper_D", out_dir)


def fig_e(res, out_dir=None):
    rows, rrows = res["rows"], res["ref_rows"]
    fig, axes = plt.subplots(1, 3, figsize=(COL_WIDTH, 80 * MM))
    ax = axes[0]
    groups = [("girder", res["strata_girder"]),
              ("orbit", res["strata_orbit"]),
              ("season", res["strata_season"])]
    lab, val, col = [], [], []
    for gname, block in groups:
        for k in sorted(block):
            b = block[k]
            lab.append(f"{gname[:4]}:{str(k)[:9]}")
            val.append(10.0 * abs(b.get("delta")) if b.get("delta") is not None else 0.0)
            col.append("#4C72B0" if b.get("delta") is not None else "#DDDDDD")
    ax.barh(np.arange(len(lab)), val, color=col, height=0.62)
    ax.set_yticks(np.arange(len(lab)))
    ax.set_yticklabels(lab, fontsize=6.0)
    ax.invert_yaxis()
    ax.set_xlabel("SDS = 10|delta| within the stratum")
    panel(ax, "a", "strata (pre vs. post)")

    ax = axes[1]
    for sv in core.STATES:
        sel = [r for r in rrows if r["state"] == RKEY[sv]]
        ax.plot([r["coherence_masked_pixels"] for r in sel],
                [r["gamma2"] for r in sel], "o", ms=1.8, alpha=0.35,
                color=STATE_COLORS[sv], markeredgewidth=0, label=sv)
    ax.set_xlabel("A = masked pixels (n_masked)")
    ax.set_ylabel("gamma2")
    rb = res["mask_clutter"]
    ax.text(0.02, 0.98, f"rho={rb['pooled']['rho']:.3f} (pooled)\n"
                        f"  ASC {rb['ASCENDING']['rho']:.3f}\n"
                        f"  DESC {rb['DESCENDING']['rho']:.3f}",
            transform=ax.transAxes, va="top", fontsize=6.6)
    ax.legend(frameon=False, markerscale=3)
    panel(ax, "b", "the mask-size co-variate")

    ax = axes[2]
    for sv in core.STATES:
        sel = [r for r in rows if r["state"] == sv and r.get("wind_speed_ms") is not None]
        ax.plot([r["wind_speed_ms"] for r in sel], [r["gamma2"] for r in sel], "o",
                ms=1.8, alpha=0.35, color=STATE_COLORS[sv], markeredgewidth=0)
    rho_w, p_w = st.spearman([r["wind_speed_ms"] for r in rows
                              if r.get("wind_speed_ms") is not None],
                             [r["gamma2"] for r in rows
                              if r.get("wind_speed_ms") is not None])
    ax.set_xlabel("wind speed (m/s)")
    ax.set_ylabel("gamma2")
    ax.text(0.02, 0.98, f"rho(gamma2, wind)={rho_w:.3f}\np={p_w:.2e}",
            transform=ax.transAxes, va="top", fontsize=6.6)
    panel(ax, "c", "the weather co-variate")
    fig.suptitle("(E) What else moves with the collapse date?",
                 fontsize=9.0, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    return save_fig(fig, "fig_carola_ocv_paper_E", out_dir)


# ---------------------------------------------------------------------------
# Pins
# ---------------------------------------------------------------------------
# Exact doubles: a hand copy of the numbers can be off in the last bits, so the
# module tolerates 1e-12 on the few values that come out of scipy, numpy or a
# permutation/bootstrap loop. A deviation above the tolerance aborts the script.
FLOAT_TOL = 1e-12
TOL_KEYS = {"welch_t", "welch_p", "mannwhitney_p", "p", "perm_p", "F",
            "p_accuracy", "p_balanced_accuracy", "accuracy_null_mean",
            "balanced_accuracy_null_mean", "accuracy_null_p95",
            "balanced_accuracy_null_p95", "delta_r2_from_adding_n_masked",
            "delta_r2_from_adding_state",
            # the OLS report of the analysis project: the same closed-form
            # numbers, but the reference's sums ran in a different order, so the
            # last bits differ (observed <= 6e-13 absolute).
            "beta", "se_classical", "t_classical", "p_classical",
            "se_hc3", "t_hc3", "p_hc3", "r2", "adj_r2", "rss", "tss",
            "aic", "bic"}

# Structural pin of ``carola_bridge_osm.json``: the committed map extract cannot
# be re-derived from the record, so its identity is asserted instead — the file
# is the Carolabrücke corridor (B 170, Wikidata Q1044279) with the two razed
# DVB tram ways of the collapse.
OSM_EXPECTED = {
    "version": 0.6,
    "generator_prefix": "Overpass API",
    "n_elements": 4,
    "n_bridge_yes": 2,
    "n_bridge_razed": 2,
    "name": "Carolabrücke",
    "road_ref": "B 170",
    "wikidata": "Q1044279",
    "nodes_per_bridge_way": 12,
    "razed_operator": "DVB",
}


def osm_facts(d):
    """The immutable facts of the committed map extract."""
    ways = [e for e in d["elements"] if e.get("type") == "way"]
    bridge = [e for e in ways if e.get("tags", {}).get("bridge") == "yes"]
    razed = [e for e in ways if e.get("tags", {}).get("bridge") == "razed"]
    tags = bridge[0].get("tags", {}) if bridge else {}
    return {
        "version": d["version"],
        "generator_prefix": " ".join(str(d["generator"]).split()[:2]),
        "n_elements": len(d["elements"]),
        "n_bridge_yes": len(bridge),
        "n_bridge_razed": len(razed),
        "name": tags.get("name"),
        "road_ref": tags.get("ref"),
        "wikidata": tags.get("wikidata"),
        "nodes_per_bridge_way": len(bridge[0]["nodes"]) if bridge else None,
        "razed_operator": (razed[0].get("tags", {}).get("operator") if razed else None),
    }


class Pins:
    """Collects pin comparisons (recomputed vs. committed) and reports deviations."""

    def __init__(self):
        self.checks = []
        self.fails = []

    def add(self, path, got, want, tol):
        try:
            exact = bool(got == want)
        except Exception:                       # pragma: no cover - numpy arrays
            exact = False
        dev = 0.0
        if not exact:
            numeric = (isinstance(got, (int, float)) and isinstance(want, (int, float))
                       and not isinstance(got, bool) and not isinstance(want, bool))
            dev = abs(float(got) - float(want)) if numeric else None
        ok = exact or (dev is not None and dev <= tol)
        self.checks.append({"path": path, "kind": "exact" if exact else "tolerance",
                            "abs_deviation": dev, "tolerance": tol, "ok": ok})
        if not ok:
            self.fails.append({"path": path, "recomputed": repr(got)[:200],
                               "committed": repr(want)[:200], "abs_deviation": dev})

    def cmp(self, prefix, got, want, key=None):
        """Recursive comparison: ``want`` dictates the structure and expectation."""
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
            return
        if isinstance(want, (list, tuple)):
            glen = len(got) if hasattr(got, "__len__") else -1
            if not isinstance(got, (list, tuple)) or glen != len(want):
                self.add(prefix, f"<len {glen}>", f"<len {len(want)}>", 0.0)
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
            "n_within_tolerance": sum(1 for c in self.checks
                                      if c["kind"] == "tolerance"),
            "n_failed": len(self.fails),
            "max_abs_deviation": max(devs) if devs else 0.0,
            "float_tolerance": FLOAT_TOL,
            "failures": self.fails[:40],
        }


def sha256(path):
    import hashlib
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def pin_all(res, quick=False):
    """All pins against the three committed reference JSONs and the metadata."""
    meta, ref_mask, ref_states = res["meta"], res["ref_mask"], res["ref_states"]
    rows, rrows, man = res["rows"], res["ref_rows"], res["manifest"]
    P = Pins()

    # 1) the inputs themselves: vendored bytes, shapes, definitions ----------
    P.cmp("meta.generator", meta["generator"], "code/carola/carola_ocv_channels_csv.py")
    P.cmp("meta.out_csv", meta["out_csv"], "carola_ocv_channels.csv")
    P.cmp("meta.n_rows", meta["n_rows"], len(rows))
    P.cmp("meta.n_columns", meta["n_columns"], len(core.CSV_COLUMNS))
    P.cmp("meta.event", meta["event"], core.COLLAPSE_DATE)
    P.cmp("meta.ocv", meta["ocv"], list(core.OCV))
    P.cmp("meta.features", meta["features"], list(core.FEATURES))
    P.cmp("meta.girders", meta["girders"], list(core.GIRDERS))
    P.cmp("meta.mask_columns", meta["mask_columns"], list(core.MASK_COLUMNS))
    P.cmp("meta.mask_window", meta["mask_window"], "80x80")
    P.cmp("meta.inputs.mask_cache.sha256", meta["inputs"]["mask_cache"]["sha256"],
          sha256(core.CACHE_PATH))
    P.cmp("meta.inputs.mask_cache_manifest.sha256",
          meta["inputs"]["mask_cache_manifest"]["sha256"], sha256(core.MANIFEST_PATH))
    P.cmp("meta.inputs.measurements.sha256",
          meta["inputs"]["measurements"]["sha256"], sha256(core.MEAS_PATH))
    P.cmp("meta.inputs.segments.sha256", meta["inputs"]["segments"]["sha256"],
          sha256(core.SEGMENTS_PATH))
    P.cmp("meta.mask_rule", meta["mask_rule"], masks.load_manifest()["mask_rule"])
    P.cmp("meta.summary_by_state", meta["summary_by_state"], res["summary_by_state"])

    # 2) the mask cache's own guard vs. the reference's guard ----------------
    P.cmp("manifest.window_sizes", man["window_sizes"],
          ref_mask["meta"]["guard"]["window_sizes"])
    P.cmp("manifest.n_deduplicated_same_payload",
          man["n_deduplicated_same_payload"], 472)
    P.cmp("manifest.n_rows", man["n_rows"], ref_mask["meta"]["rows"])
    P.cmp("manifest.mask_rule.peak_frac", man["mask_rule"]["peak_frac"],
          ref_mask["meta"]["params"].get("peak_frac", 0.30))
    P.cmp("manifest.segment_resolved", man["segment_resolved"], True)
    P.cmp("manifest.asset_date_groups_with_several_segments",
          man["asset_date_groups_with_several_segments"], 264)

    # 3) the reference-convention aggregates, number by number ---------------
    P.cmp("states", res["states"], ref_mask["states"])
    P.cmp("by_segment", res["by_segment"], ref_mask["by_segment"])
    P.cmp("by_request", res["by_request"], ref_mask["by_request"])
    P.cmp("sds", res["sds"], ref_mask["sds"])
    P.cmp("models", res["models"], ref_mask["models"])

    # 4) the per-window records of the first analysis ------------------------
    by_mid = {r["id"]: r for r in rows}
    n_seen = n_ok = 0
    for rec in ref_states["records"]:
        r = by_mid.get(rec["mid"])
        if r is None:
            continue
        n_seen += 1
        g_ok = (abs(float(rec["gamma2"]) - float(r["gamma2"])) < 1e-12)
        a_ok = (int(rec["n_masked"]) == int(r["A"]))
        n_ok += 1 if (g_ok and a_ok) else 0
        if not (g_ok and a_ok):
            P.add(f"window[{rec['mid']}].gamma2", r["gamma2"], rec["gamma2"], FLOAT_TOL)
            P.add(f"window[{rec['mid']}].n_masked", r["A"], rec["n_masked"], 0.0)
    P.cmp("window_records.checked", n_seen, 925)
    P.cmp("window_records.matched", n_ok, 925)
    P.cmp("window_records.of", len(ref_states["records"]), 1454)

    # 5) cross-check: the two committed references must agree about the
    #    old 1,454-window basis (the second one quotes the first) ------------
    cr = ref_mask["states"]["committed_reference"]
    ps = ref_states["per_state"]
    P.cmp("committed_reference.healthy.n", cr["healthy"]["n"],
          ps[core.PRE]["n"])
    P.cmp("committed_reference.damaged.n", cr["damaged"]["n"], ps[core.POST]["n"])
    P.cmp("committed_reference.healthy.gamma2_median",
          cr["healthy"]["gamma2_median"],
          round(ps[core.PRE]["gamma2"]["median"], 4))
    P.cmp("committed_reference.damaged.gamma2_median",
          cr["damaged"]["gamma2_median"],
          round(ps[core.POST]["gamma2"]["median"], 4))
    P.cmp("committed_reference.healthy.n_masked_median",
          cr["healthy"]["n_masked_median"], ps[core.PRE]["n_masked"]["median"])
    P.cmp("committed_reference.damaged.n_masked_median",
          cr["damaged"]["n_masked_median"], ps[core.POST]["n_masked"]["median"])
    P.cmp("committed_reference.source", cr["source"],
          "carola_coherence_states.md (1454 stored windows)")

    # 6) the reference's own summary of the current export -------------------
    P.cmp("ref_mask.meta.event", ref_mask["meta"]["event"], core.COLLAPSE_DATE)
    P.cmp("ref_mask.meta.rows", ref_mask["meta"]["rows"], len(rrows))
    P.cmp("ref_osm", osm_facts(res["ref_osm"]), OSM_EXPECTED)

    # 7) this package's own numbers must be finite and self-consistent -------
    ch = res["channels"]
    for k in core.FEATURES:
        if "delta" in ch.get(k, {}):
            P.add(f"channels[{k}].sds>=0", 10.0 * abs(ch[k]["delta"]) >= 0.0, True, 0.0)
    P.cmp("n_rows", res["n_rows"], len(rows))
    P.cmp("lda.headline_lambda", res["lda"]["headline_lambda"], core.HEADLINE_LAMBDA)
    if quick:
        P.cmp("lda.perm.recomputed",
              res["lda"]["headline_permutation_test"].get("recomputed"), False)
    else:
        P.cmp("lda.perm.recomputed",
              res["lda"]["headline_permutation_test"].get("recomputed", True), True)
    return P


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------
def fmt(x, nd=3):
    if x is None:
        return "-"
    if isinstance(x, float):
        return f"{x:.{nd}f}"
    return str(x)


def md_table(header, rows):
    out = ["| " + " | ".join(header) + " |",
           "| " + " | ".join("---" for _ in header) + " |"]
    for r in rows:
        out.append("| " + " | ".join(str(c) for c in r) + " |")
    return out


def report(res, pins, quick):
    man = res["manifest"]
    ch, gd = res["channels"], res["girders"]
    md = []
    md.append(f"# {core.SITE} — observability core vector, pre/post collapse")
    md.append("")
    md.append(f"**Site** {core.SITE} · **event** {core.COLLAPSE_DATE} · "
              f"**rows** {res['n_rows']} 80x80 girder chips · "
              f"**states** `{core.PRE}` / `{core.POST}`")
    md.append("")
    md.append("```")
    md.append("OCV  = [gamma2, P, D]              observability core vector")
    md.append("x    = [gamma2, P, D, A, F, S]      echo-mask / deck-mask vector")
    md.append("mask = intensity >= 0.30*peak AND >= 5.0*np.median(intensity), A >= 2")
    md.append("layers: echo mask (whole chip)  |  deck mask (per girder A/B/C)")
    md.append("```")
    md.append("")
    md.append("## What is recomputed, and what is pinned")
    md.append("")
    md.append("| input | role |")
    md.append("| --- | --- |")
    md.append(f"| `data/carola/{os.path.basename(core.CACHE_PATH)}` | the echo mask of "
              f"every committed chip (the mask layer's input) |")
    md.append(f"| `data/carola/{os.path.basename(core.MEAS_PATH)}` | the committed "
              f"measurement extract (channels, weather, geometry) |")
    md.append(f"| `data/carola/{os.path.basename(core.SEGMENTS_PATH)}` | the girder "
              f"labels + FEM baselines |")
    md.append(f"| `data/carola/carola_ocv_channels.csv` | generated: one row per chip, "
              f"{res['meta']['n_columns']} columns |")
    md.append("")
    md.append(f"The mask cache holds **{man['n_rows']}** 80x80 chips "
              f"({man['window_sizes'].get('7x7 (excluded)', 0)} 7x7 tower chips are "
              f"excluded), after dropping "
              f"**{man['n_deduplicated_same_payload']}** byte-identical duplicates "
              f"(the per-girder re-extraction resolves the segments: "
              f"`segment_resolved = {man['segment_resolved']}`, "
              f"{man['asset_date_groups_with_differing_payloads']} acquisition groups "
              f"with differing payloads).")
    md.append("")
    md.append("### The two mask layers")
    md.append("")
    md.append("* **echo mask** — the whole 80x80 chip; the layer of the 6D vector.")
    md.append("* **deck mask** — the same rule on the girder's own chip "
              "(`segment_index` of the export index, labelled by "
              "`carola_segments.txt`).")
    md.append("")
    md.append("Neither layer's `A` is the pipeline's own `coherence_masked_pixels`: "
              "that column was computed on the pre-fix *asset-point* chip while the "
              "committed export holds the per-girder re-extraction, so the two agree "
              "only by coincidence "
              f"({res['meta']['n_coherence_masked_pixels_matched']} of "
              f"{res['meta']['n_coherence_masked_pixels_compared']} rows). The site's "
              "analysis project therefore recomputed the mask from the payloads, and "
              "this package re-derives it from the committed cache.")
    md.append("")
    md.append("## A — the channels")
    md.append("")
    md.append("Raw distributions of the six mask channels per state (fig. A).")
    md.append("")
    body = []
    for k in core.FEATURES:
        b = ch.get(k, {})
        if "delta" not in b:
            continue
        body.append([k, fmt(b["median"][0], 4), fmt(b["median"][1], 4),
                     fmt(b["delta"]), fmt(b["sds"]),
                     "yes" if b.get("ci_excludes_zero") else "no"])
    md += md_table(["channel", "pre median", "post median", "delta", "SDS",
                    "CI95 excludes 0"], body)
    md.append("")
    md.append("## B — the two layers side by side")
    md.append("")
    md.append("Cliff's delta (pre -> post) of the deck-mask channels inside each "
              "girder:")
    md.append("")
    gir = [g for g in core.GIRDERS if g in gd]
    body = []
    for g in gir:
        b = gd[g]
        body.append([g, b["n_rows"], b["n_pre"], b["n_post"]]
                    + [fmt(b["channels"].get(k, {}).get("delta")) for k in core.FEATURES])
    md += md_table(["girder", "n", "pre", "post"] + core.FEATURES, body)
    md.append("")
    md.append("The reference's own per-girder table (recomputed, fig. B/C):")
    md.append("")
    body = []
    for g, b in res["by_segment"].items():
        body.append([g, b["n"], b["n_healthy"], b["n_damaged"],
                     fmt(b["gamma2_median"], 4), fmt(b["n_masked_median"], 1),
                     fmt(b["rho"]), "yes" if b["repair_target"] else "no"])
    md += md_table(["girder", "n", "healthy", "damaged", "gamma2 median",
                    "A median", "rho(gamma2, A)", "repair target"], body)
    md.append("")
    md.append("## C — is the state recoverable?")
    md.append("")
    lda = res["lda"]
    head = lda.get("headline") or {}
    md.append(f"2-class LDA on the 6D vector, leave-one-out at "
              f"`lambda = {lda['headline_lambda']}`: accuracy "
              f"**{fmt(head.get('accuracy'))}**, balanced accuracy "
              f"**{fmt(head.get('balanced_accuracy'))}** (chance 0.5).")
    md.append("")
    pt = lda.get("headline_permutation_test") or {}
    if not quick and "p_accuracy" in pt:
        md.append(f"Label permutation ({pt['n_perm']} draws): "
                  f"p(accuracy) = {pt['p_accuracy']:.4f}, "
                  f"p(balanced) = {pt['p_balanced_accuracy']:.4f}; "
                  f"null mean accuracy {fmt(pt['accuracy_null_mean'])}.")
    else:
        md.append("Label permutation skipped (`--quick`).")
    md.append("")
    md.append("## D — does the extra structure help?")
    md.append("")
    body = []
    for name, b in res["oof"]["per_model"].items():
        pb = b.get("paired_vs_gamma2") or {}
        ci = pb.get("ci95") or [None, None]
        body.append([SHORT_MODEL.get(name, name),
                     "+".join(core.MODELS.get(name, [])),
                     fmt(b.get("sds")), fmt(pb.get("diff_mean")),
                     f"[{fmt(ci[0])}, {fmt(ci[1])}]"])
    md += md_table(["model", "features", "SDS (OOF)", "delta vs gamma2", "CI95"], body)
    md.append("")
    md.append("## E — what else moves with the collapse date?")
    md.append("")
    md.append("Stratified SDS (= 10|Cliff's delta|) of `gamma2` within each stratum:")
    md.append("")
    body = []
    for sname, block in (("girder", res["strata_girder"]),
                         ("orbit", res["strata_orbit"]),
                         ("season", res["strata_season"])):
        for k, b in block.items():
            body.append([sname, k, b["n"][0], b["n"][1],
                         fmt(10.0 * abs(b["delta"]) if b.get("delta") is not None
                             else None)])
    md += md_table(["stratifier", "stratum", "n pre", "n post", "SDS gamma2"], body)
    md.append("")
    rb = res["mask_clutter"]
    m = res["models"]
    s = res["sds"]
    md.append("The mask-size co-variate (fig. E panel b) and the nested OLS "
              "(the reference's `model_block`):")
    md.append("")
    md.append(f"* rho(gamma2, A) = {fmt(rb['pooled']['rho'])} pooled, "
              f"{fmt(rb['ASCENDING']['rho'])} ascending, "
              f"{fmt(rb['DESCENDING']['rho'])} descending.")
    md.append(f"* `gamma2 ~ damaged`: r2 = {fmt(m['model1_state_only']['r2'])}, "
              f"beta(damaged) = {fmt(m['model1_state_only']['beta'][1], 4)} "
              f"(HC3 p = {fmt(m['model1_state_only']['p_hc3'][1], 3)}).")
    md.append(f"* `gamma2 ~ A`: r2 = {fmt(m['model2_mask_only']['r2'])}, "
              f"beta(A) = {fmt(m['model2_mask_only']['beta'][1], 5)}.")
    md.append(f"* `gamma2 ~ damaged + A`: r2 = {fmt(m['model3_state_plus_mask']['r2'])}; "
              f"delta r2 = +{fmt(m['delta_r2_from_adding_n_masked'])} from A and "
              f"+{fmt(m['delta_r2_from_adding_state'])} from the state.")
    md.append(f"* inverse-N law: r2 = {fmt(m['inverse_n_fit']['r2'])}.")
    # ``scenario_reading`` of the reference is its own (German) wording and is
    # pinned as such; the report renders the reading in English.
    scen_en = {
        "A": "SDS(n_masked) ~ 0, so the mask carries no clear state effect",
        "B": "SDS(n_masked) < SDS(gamma2), so a state effect on n_masked "
             "could reach gamma2 through the mask (an indirect channel)",
    }.get(s["scenario"], s["scenario_reading"])
    md.append(f"* SDS(gamma2) = {fmt(s['gamma2']['sds'])}, "
              f"SDS(A) = {fmt(s['n_masked']['sds'])}, scenario **{s['scenario']}** — "
              f"{scen_en}.")
    md.append("")
    md.append("## Pins")
    md.append("")
    ps = pins.summary()
    md.append(f"{ps['n_checks']} checks against the committed references "
              f"({ps['n_exact']} exact, {ps['n_within_tolerance']} within "
              f"{ps['float_tolerance']:g}); **{ps['n_failed']} deviations**.")
    md.append("")
    md.append("| reference | what is pinned |")
    md.append("| --- | --- |")
    pw = res["pins_windows"]
    md.append("| `carola_coherence_states.json` | per-window `gamma2` and `n_masked` "
              f"({pw['checked']} of its {pw['of']} stored windows are in the current "
              f"export; {pw['matched']} match exactly) |")
    md.append("| `carola_echo_mask_gamma2.json` | the state table, the per-girder and "
              "per-request tables, the SDS block and the nested OLS — number by number |")
    md.append("| `carola_bridge_osm.json` | the map extract structurally "
              "(Carolabrücke B 170, Q1044279, two 12-node bridge ways, two razed DVB "
              "ways) |")
    md.append("")
    md.append("### Reading (not a claim of damage)")
    md.append("")
    md.append(f"The echo-mask coherence rises across the collapse "
              f"(SDS {fmt(s['gamma2']['sds'])}), while the mask size `A` moves much "
              f"less (SDS {fmt(s['n_masked']['sds'])}). Both states sit on the same "
              "bridge, so season, weather, traffic and pipeline generation move with "
              "the state too — figures C–E report that co-variation rather than a "
              "mechanism. What the figures support is the *observability* statement: "
              "the vector separates the two epochs of the record, and the 6D deck-mask "
              "vector does so per girder.")
    md.append("")
    md.append("```bash")
    md.append("python3 code/carola/carola_ocv_masks.py --check      # the mask rule")
    md.append("python3 code/carola/carola_ocv_channels_csv.py       # build + verify")
    md.append("python3 code/carola/fig_carola_ocv_paper.py          # figures + pins")
    md.append("```")
    return "\n".join(md) + "\n"


# ---------------------------------------------------------------------------
# JSON artifact and entry point
# ---------------------------------------------------------------------------
def json_default(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    return str(o)


def build_json(res, pins_summary, quick):
    """The result JSON — every recomputed number, without the per-row rows."""
    return {
        "site": core.SITE,
        "event": core.COLLAPSE_DATE,
        "state_labels": list(core.STATES),
        "ocv": list(core.OCV),
        "features": list(core.FEATURES),
        "girder_labels": list(core.GIRDERS),
        "headline_lambda": core.HEADLINE_LAMBDA,
        "rng_seed": core.RNG_SEED,
        "n_boot": core.N_BOOT,
        "n_perm": core.N_PERM,
        "quick": quick,
        "generated_by": "code/carola/fig_carola_ocv_paper.py",
        "windows": res["windows"],
        "meta": res["meta"],
        "pins_windows": res["pins_windows"],
        "states": res["states"],
        "by_segment": res["by_segment"],
        "by_request": res["by_request"],
        "sds": res["sds"],
        "models": res["models"],
        "mask_clutter": res["mask_clutter"],
        "channels": res["channels"],
        "girders": res["girders"],
        "strata_girder": res["strata_girder"],
        "strata_orbit": res["strata_orbit"],
        "strata_season": res["strata_season"],
        "lda": res["lda"],
        "oof_per_model": res["oof"]["per_model"],
        "summary_by_state": res["summary_by_state"],
        "pins": pins_summary,
        # no wall-clock time: with the fixed seeds two runs must produce the
        # identical JSON, so the runtime stays on stdout only.
        "figures": [f"figures/carola/fig_carola_ocv_paper_{k}.png" for k in "ABCDE"],
    }


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--quick", action="store_true",
                    help="skip the label-permutation null of panel C")
    ap.add_argument("--figdir", default=FIGDIR)
    ap.add_argument("--json", default=OUT_JSON)
    ap.add_argument("--md", default=OUT_MD)
    args = ap.parse_args(argv)

    t0 = time.time()
    print(f"computing — {core.SITE} …")
    # Developer aid: CAROLA_RES_CACHE=<path> reuses a pickled ``compute()``
    # result between runs (the bootstrap makes a cold run a few minutes).
    cache = os.environ.get("CAROLA_RES_CACHE")
    if cache and os.path.exists(cache):
        with open(cache, "rb") as fh:
            res = pickle.load(fh)
        print(f"  loaded {cache} ({res['n_rows']} chips)")
    else:
        res = compute(quick=args.quick)
        print(f"  {res['n_rows']} chips, {len(res['by_segment'])} girders, "
              f"{res['runtime_s']:.1f} s")
        if cache:
            with open(cache, "wb") as fh:
                pickle.dump(res, fh)
            print(f"  wrote {cache}")

    pins = pin_all(res, quick=args.quick)
    summary = pins.summary()
    print(f"pins: {summary['n_checks']} checks, {summary['n_exact']} exact, "
          f"{summary['n_within_tolerance']} within {FLOAT_TOL:g}, "
          f"{summary['n_failed']} failed")
    for f in summary["failures"]:
        print("PIN FAILED:", f["path"], "recomputed=", f["recomputed"],
              "committed=", f["committed"])

    for fn in (fig_a, fig_b, fig_c, fig_d, fig_e):
        fn(res, args.figdir)

    with open(args.json, "w") as fh:
        json.dump(build_json(res, summary, args.quick), fh, indent=1,
                  sort_keys=True, default=json_default)
    print(f"  {os.path.relpath(args.json)}")
    with open(args.md, "w") as fh:
        fh.write(report(res, pins, args.quick))
    print(f"  {os.path.relpath(args.md)}")

    if summary["n_failed"]:
        for f in summary["failures"][:20]:
            print("PIN FAILED:", f["path"], f["recomputed"], "!=", f["committed"])
        raise SystemExit(f"PIN VERIFICATION FAILED ({summary['n_failed']})")
    print(f"PIN VERIFICATION OK  ({time.time() - t0:.1f} s)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
