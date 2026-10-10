#!/usr/bin/env python3
"""Espoo OCV-Paper — statistics ports (verbatim from the LUMO project) plus
the two event-free layers of this site.

This module invents *nothing*: the first block is a byte-for-byte copy of
``code/lumo/lumo_ocv_stats.py``, itself a verbatim port of

  * ``analyze_lumo_tower_coherence.py``        (``_stats``, ``_test``)
  * ``analyze_lumo_gamma2_pairs.py``           (Cliff's delta, Bootstrap-CI,
                                                ``delta_block``, ``sign_test``,
                                                ``sds_from``, ``by_state``,
                                                ``vals``, ``pairs_block``,
                                                Jonckheere-Terpstra,
                                                ``season_of``)
  * ``analyze_lumo_emv_multivariate_cv.py``    (``fit_lda``, ``predict_lda``,
                                                ``loo_cv``, ``permutation_test``)
  * ``analyze_lumo_emv_pairwise_triplets_cv.py`` (``fit_linear_score``,
                                                ``score_row``, ``oof_scores``,
                                                ``delta_sds``, ``permutation_p``,
                                                ``paired_bootstrap_diff``)

The Espoo difference is *naming and scope only*, and it is a large one:

  * **the states are orbits, not damage labels.** ``STATE_ORDER`` is
    ``ASCENDING`` / ``DESCENDING`` and the state field is ``state``; the record
    of this site carries **no** ``damage_label`` / ``structural_state`` /
    ``condition_label`` on any of its 150 rows. The two states are therefore a
    *geometry* cut, not a condition cut — see ``code/espoo/README.md``.
  * **one measured mask dimension exists.** The record delivers only
    ``A = coherence_masked_pixels``; there is no echo-mask layer of this package
    (no ``peak_intensity``, no ``sub_aperture_brightness``, no ``dwell_s``), so
    ``DIMS = ["A"]`` and every mask-morphology helper of the LUMO port stays
    ported but unused.
  * **this package adds two layers the event sites do not need** — ported
    line by line from the site's own upstream scripts, so that the pin block can
    reproduce them:
      - ``group_stats`` / ``classify`` / ``reference_for`` / ``compare_orbit``
        port ``espoo_mast_observability.py`` (the per-orbit observability block
        and the orbit test);
      - ``mann_kendall`` / ``theil_sen`` / ``block_bootstrap_slope`` /
        ``lag1_autocorr`` / ``runs_test`` / ``ljung_box`` / ``annual_harmonic``
        / ``split_half_series`` / ``stationarity_flags`` / ``stationarity_of``
        / ``phase_series_block`` port ``espoo_phase_stationarity.py`` (the
        Theil-Sen / Mann-Kendall stationarity layer). The phase scan is a
        *different* input file than the channel table, but
        ``espoo_phase_stationarity.json`` commits the **full per-acquisition
        series** of all 14 of its groups under ``series_by_group``, so every
        statistical number of that file — median, IQR, Theil-Sen slope, its
        moving-block bootstrap CI, Mann-Kendall p, lag-1 autocorrelation (raw and
        detrended), runs test, Ljung-Box Q, split-half Cliff's delta with its CI,
        the annual harmonic and the monthly means — is recomputed here from the
        committed series and pinned field by field, together with the flag rule
        that turns those primitives into the committed ``verdict``. Only the keys
        the file does not carry (``dwell_s``, ``n_sub``, ``nyquist_hz``,
        ``n_observable``, ``dwell_control``) are not re-derived; the report names
        them.

Conventions (identical to the origin):
  SDS  = 10 * |Cliff's delta|
  delta(ref, grp) > 0  =>  the *group* is larger than the reference
  bootstrap seeds 7 (in-sample) resp. 11/12/13 (OOF), permutation seed 7/8
"""
from __future__ import annotations

import collections
import datetime as dt
from math import erf, sqrt

import numpy as np

try:
    from scipy import stats as sps
except ImportError:  # pragma: no cover - statistical tests are optional
    sps = None

# The two observation geometries of the record. The order is the declaration
# order of the upstream report (ASCENDING first, the geometry that carries the
# usable echo) and it is the reference side of every comparison below.
STATE_ORDER = ["ASCENDING", "DESCENDING"]
STATE_SHORT = {"ASCENDING": "ASC", "DESCENDING": "DESC"}
# ``DAMAGED`` keeps the LUMO name: it is the *group* side of ``pairs_block`` —
# with two states that is simply the second one.
DAMAGED = ["DESCENDING"]
# ``SEVERITY`` keeps the LUMO name too (the monotone axis of the strata/trend
# helpers). Here it is only an index, not a claim: 0 = ASCENDING, 1 =
# DESCENDING. No monotonicity claim follows from it and the figure script
# reports ``monotonicity`` as *not applicable* rather than reading it.
SEVERITY = {"ASCENDING": 0, "DESCENDING": 1}

BOOT_MAIN = 10000
BOOT_STRATA = 2000
WINTER_MONTHS = (10, 11, 12, 1, 2)
SUMMER_MONTHS = (3, 4, 5, 6, 7)


def _stats(vals):
    """Descriptive statistics, exactly as ``alc._stats``."""
    vals = [v for v in vals if v is not None]
    if not vals:
        return None
    a = np.asarray(vals, dtype=float)
    return {
        "n": int(len(vals)),
        "median": float(np.median(a)),
        "mean": float(np.mean(a)),
        "std": float(np.std(a, ddof=1)) if len(vals) > 1 else 0.0,
        "p5": float(np.percentile(a, 5)),
        "p25": float(np.percentile(a, 25)),
        "p75": float(np.percentile(a, 75)),
        "p95": float(np.percentile(a, 95)),
    }


def welch_mw_test(a, b):
    """Welch t + Mann-Whitney U, exactly as ``alc._test``."""
    x = np.asarray(a, dtype=float)
    y = np.asarray(b, dtype=float)
    out = {"n": [int(len(x)), int(len(y))],
           "median": [float(np.median(x)), float(np.median(y))]}
    if sps is not None and len(x) >= 3 and len(y) >= 3:
        t, p = sps.ttest_ind(x, y, equal_var=False)
        out["welch_t"] = float(t)
        out["welch_p"] = float(p)
        try:
            u, pu = sps.mannwhitneyu(x, y, alternative="two-sided")
            out["mannwhitney_u"] = float(u)
            out["mannwhitney_p"] = float(pu)
        except ValueError:
            pass
    return out


def cliffs_delta_bruteforce(ref, grp):
    """The literal definition, O(n*m) — reference implementation."""
    ref = [v for v in ref if v is not None]
    grp = [v for v in grp if v is not None]
    if not ref or not grp:
        return None
    gt = sum(1 for x in ref for y in grp if y > x)
    lt = sum(1 for x in ref for y in grp if y < x)
    return (gt - lt) / (len(ref) * len(grp))


def _u_from_ranks(a, b, rng=None):
    """``U = #(b>a) + 0.5*#(b=a)`` from mid-ranks of the pooled sample."""
    n, m = a.size, b.size
    pooled = np.concatenate([a, b]) if rng is None else np.concatenate([
        a[rng.integers(0, n, n)], b[rng.integers(0, m, m)]])
    order = pooled.argsort(kind="mergesort")
    ranks = np.empty(pooled.size, float)
    sp = pooled[order]
    i = 0
    while i < sp.size:                      # average ranks across ties
        j = i
        while j + 1 < sp.size and sp[j + 1] == sp[i]:
            j += 1
        ranks[order[i:j + 1]] = 0.5 * (i + j) + 1.0
        i = j + 1
    w_b = ranks[n:].sum()
    return w_b - m * (m + 1) / 2.0


def cliffs_delta(ref, grp):
    """Cliff's delta; positive sign => group > reference."""
    a = np.asarray([v for v in ref if v is not None], float)
    b = np.asarray([v for v in grp if v is not None], float)
    if a.size == 0 or b.size == 0:
        return None
    return float((2.0 * _u_from_ranks(a, b) - a.size * b.size) / (a.size * b.size))


def bootstrap_delta_ci(ref, grp, n=BOOT_MAIN, seed=7):
    """Percentile 95% CI of delta (both groups resampled independently)."""
    a = np.asarray([v for v in ref if v is not None], float)
    b = np.asarray([v for v in grp if v is not None], float)
    if a.size < 2 or b.size < 2:
        return None
    rng = np.random.default_rng(seed)
    out = np.empty(n)
    for i in range(n):
        u = _u_from_ranks(a, b, rng)
        out[i] = (2.0 * u - a.size * b.size) / (a.size * b.size)
    return float(np.percentile(out, 2.5)), float(np.percentile(out, 97.5))


def delta_block(ref, grp, label, n_boot=BOOT_MAIN, seed=7):
    """delta + CI + Welch/Mann-Whitney + SDS for one comparison."""
    out = {"comparison": label, "n": [len(ref), len(grp)]}
    d = cliffs_delta(ref, grp)
    if d is None:
        return out
    out["median"] = [float(np.median(ref)), float(np.median(grp))]
    out["mean"] = [float(np.mean(ref)), float(np.mean(grp))]
    out["delta"] = d
    out["sds"] = 10.0 * abs(d)
    ci = bootstrap_delta_ci(ref, grp, n=n_boot, seed=seed)
    if ci:
        out["delta_ci95"] = [ci[0], ci[1]]
        out["ci_excludes_zero"] = bool(ci[0] * ci[1] > 0)
    out.update(welch_mw_test(ref, grp))
    return out


def sign_test(diffs):
    """Two-sided exact sign test on the differences that are not zero."""
    nz = [d for d in diffs if d != 0]
    if not nz:
        return {"n": 0, "n_positive": 0, "p": None}
    pos = sum(1 for d in nz if d > 0)
    if sps is not None:
        p = float(sps.binomtest(pos, len(nz), 0.5).pvalue)
    else:  # pragma: no cover
        from math import comb
        k = min(pos, len(nz) - pos)
        p = float(min(1.0, 2 * sum(comb(len(nz), i) for i in range(k + 1)) / 2 ** len(nz)))
    return {"n": len(nz), "n_positive": pos, "p": p}


def sds_from(deltas):
    """``max_single_group`` SDS as in the backend: 10*max|delta|."""
    vals = [abs(d) for d in deltas if d is not None]
    return 10.0 * max(vals) if vals else None


def by_state(rows):
    """Rows per observation geometry (the ``state`` field of this package)."""
    return {lab: [r for r in rows if r.get("state") == lab] for lab in STATE_ORDER}


def vals(rows, label, key="gamma2"):
    """Values of one state — accepts a row list OR a ``by_state`` dict."""
    if isinstance(rows, dict):
        recs = rows.get(label) or []
    else:
        recs = [r for r in rows if r.get("state") == label]
    return [r[key] for r in recs if r.get(key) is not None]


def pairs_block(rows, n_boot=BOOT_MAIN, key="gamma2"):
    """healthy vs. each damage state + all six state pairs in one dict."""
    st = by_state(rows)
    out = {"channel": key,
           "n_by_state": {STATE_SHORT[l]: len(st[l]) for l in STATE_ORDER},
           "healthy_vs": {}, "all_pairs": {}}
    ds = []
    for lab in DAMAGED:
        if not st[lab]:
            continue
        blk = delta_block(vals(st, "healthy", key), vals(st, lab, key),
                          f"healthy_vs_{STATE_SHORT[lab]}", n_boot)
        out["healthy_vs"][STATE_SHORT[lab]] = blk
        ds.append((STATE_SHORT[lab], blk.get("delta")))
    out["sign_decomposition"] = dict(ds)
    labels = STATE_ORDER
    for i, a in enumerate(labels):
        for b in labels[i + 1:]:
            blk = delta_block(vals(st, a, key), vals(st, b, key),
                              f"{STATE_SHORT[a]}_vs_{STATE_SHORT[b]}", n_boot)
            out["all_pairs"][blk["comparison"]] = blk
    pooled = delta_block(vals(st, "healthy", key),
                         [v for lab in DAMAGED for v in vals(st, lab, key)],
                         "healthy_vs_all_damaged", n_boot)
    out["pooled_all_damaged"] = pooled
    out["max_single_group"] = {
        "group": max(ds, key=lambda kv: abs(kv[1] or 0.0))[0] if ds else None,
        "sds": sds_from([d for _l, d in ds]),
    }
    out["sds_pooled_all_damaged"] = (
        None if pooled.get("delta") is None else 10.0 * abs(pooled["delta"]))
    return out


def _jt_statistic(groups):
    """``S = sum_{i<j} [#(b>a) + 0.5*#(b=a)]`` over groups ordered by severity."""
    s = 0.0
    for i in range(len(groups)):
        for j in range(i + 1, len(groups)):
            a, b = groups[i], groups[j]
            right = np.searchsorted(b, a, side="right")
            left = np.searchsorted(b, a, side="left")
            s += float((b.size - right).sum()) + 0.5 * float((right - left).sum())
    return s


def jonckheere(rows, n_perm=20000, seed=7, key="gamma2"):
    """Jonckheere-Terpstra trend test across severity-ordered states."""
    st = by_state(rows)
    groups = [np.sort(np.asarray(vals(st, l, key), float)) for l in STATE_ORDER]
    groups = [g for g in groups if g.size]
    obs = _jt_statistic(groups)

    labels = np.concatenate([[i] * len(st[l]) for i, l in enumerate(STATE_ORDER) if st[l]])
    g2 = np.asarray([r[key] for r in rows], float)
    rng = np.random.default_rng(seed)
    lo = hi = 0
    for _ in range(n_perm):
        perm = rng.permutation(labels)
        gs = [np.sort(g2[perm == i]) for i in range(len(STATE_ORDER))]
        s = _jt_statistic([g for g in gs if g.size])
        if s <= obs:
            lo += 1
        if s >= obs:
            hi += 1
    n1 = n_perm + 1
    return {
        "statistic": obs,
        "p_decreasing": (lo + 1) / n1,
        "p_increasing": (hi + 1) / n1,
        "p_two_sided": min(1.0, 2 * min(lo + 1, hi + 1) / n1),
        "n_perm": n_perm,
        "severity_ranks": {STATE_SHORT[l]: SEVERITY[l] for l in STATE_ORDER},
    }


def monotonicity(rows, pair_block, key="gamma2", n_perm=20000):
    """The literal hypothesis of the plan plus a trend test on the severity axis."""
    d = {k: (pair_block["healthy_vs"].get(k) or {}).get("delta")
         for k in ("DAM3", "DAM4", "DAM6")}
    absd = {k: (None if v is None else abs(v)) for k, v in d.items()}
    hy = {"channel": key,
          "claim": "|delta(healthy,DAM3)| < |delta(healthy,DAM4)| < |delta(healthy,DAM6)|",
          "delta": d, "abs_delta": absd}
    if all(v is not None for v in absd.values()):
        hy["holds"] = bool(absd["DAM3"] < absd["DAM4"] < absd["DAM6"])
        hy["order_observed"] = " < ".join(
            k for k, _v in sorted(absd.items(), key=lambda kv: kv[1]))
    out = {"hypothesis": hy}
    if sps is not None:
        sev = np.asarray([SEVERITY[r.get("state")] for r in rows], float)
        v_ = np.asarray([r[key] for r in rows], float)
        rho, p = sps.spearmanr(sev, v_)
        out["spearman_severity_vs_value"] = {
            "n": int(sev.size), "rho": float(rho), "p": float(p)}
    out["jonckheere"] = jonckheere(rows, n_perm=n_perm, key=key)
    return out


def season_of(r):
    """Season from the ``month`` field (YYYY-MM) — verbatim ``gp.season_of``."""
    m = r.get("month") or ""
    m = int(m.split("-")[1]) if "-" in m else None
    if m in WINTER_MONTHS:
        return "winter_Oct-Feb"
    if m in SUMMER_MONTHS:
        return "summer_Mar-Jul"
    return "shoulder_Aug-Sep"


def strata(rows, key_fn, n_boot=BOOT_STRATA, key="gamma2"):
    """Strata controls (season/orbit direction) — port of ``gp._strata``."""
    out = {}
    for skey in sorted({key_fn(r) for r in rows if key_fn(r) is not None}):
        sub = [r for r in rows if key_fn(r) == skey]
        st = by_state(sub)
        entry = {"n_by_state": {STATE_SHORT[l]: len(st[l]) for l in STATE_ORDER}}
        for lab in DAMAGED:
            if len(st["healthy"]) < 3 or len(st[lab]) < 3:
                entry[STATE_SHORT[lab]] = {"n": [len(st["healthy"]), len(st[lab])],
                                           "delta": None, "note": "n < 3"}
                continue
            entry[STATE_SHORT[lab]] = delta_block(
                vals(st, "healthy", key), vals(st, lab, key),
                f"healthy_vs_{STATE_SHORT[lab]}_in_{skey}", n_boot)
        out[str(skey)] = entry
    return out


# ---------------------------------------------------------------------------
# 4-class LDA (shrinkage-regularised pooled covariance) — port of mcv
# ---------------------------------------------------------------------------
def fit_lda(X_train, y_train, classes, shrinkage):
    """Uniform priors (deliberate, see ``mcv.fit_lda``): otherwise the model tips
    at healthy=120 vs. DAM4/DAM6=16 into "always healthy"."""
    mean = X_train.mean(axis=0)
    std = X_train.std(axis=0)
    std[std == 0] = 1.0
    Xs = (X_train - mean) / std
    p = Xs.shape[1]

    class_means, pooled_cov, n_dof, priors = {}, np.zeros((p, p)), 0, {}
    for c in classes:
        Xc = Xs[y_train == c]
        priors[c] = 1.0 / len(classes)
        class_means[c] = Xc.mean(axis=0) if Xc.shape[0] else np.zeros(p)
        if Xc.shape[0] > 1:
            diff = Xc - class_means[c]
            pooled_cov += diff.T @ diff
            n_dof += Xc.shape[0] - 1
    pooled_cov /= max(n_dof, 1)
    diag_cov = np.diag(np.diag(pooled_cov))
    cov = (1.0 - shrinkage) * pooled_cov + shrinkage * diag_cov
    cov += 1e-9 * np.eye(p)
    return {"mean": mean, "std": std, "class_means": class_means,
            "inv_cov": np.linalg.inv(cov), "priors": priors, "classes": classes}


def predict_lda(model, x_row):
    xs = (x_row - model["mean"]) / model["std"]
    inv_cov = model["inv_cov"]
    best_c, best_score = None, -np.inf
    for c in model["classes"]:
        mu = model["class_means"][c]
        score = xs @ inv_cov @ mu - 0.5 * (mu @ inv_cov @ mu) + np.log(model["priors"][c] + 1e-12)
        if score > best_score:
            best_score, best_c = score, c
    return best_c


def loo_cv(X, y, classes, shrinkage):
    """Leave-one-out CV: each model is fitted only on the other n-1 rows."""
    n = X.shape[0]
    y_pred = np.empty(n, dtype=object)
    for i in range(n):
        mask = np.ones(n, dtype=bool)
        mask[i] = False
        model = fit_lda(X[mask], y[mask], classes, shrinkage)
        y_pred[i] = predict_lda(model, X[i])
    acc = float((y_pred == y).mean())
    recalls = []
    confusion = {a: {b: 0 for b in classes} for a in classes}
    for true_c in classes:
        idx = (y == true_c)
        if idx.sum() == 0:
            continue
        recalls.append(float((y_pred[idx] == true_c).mean()))
        for pred_c in classes:
            confusion[true_c][pred_c] = int((y_pred[idx] == pred_c).sum())
    bal_acc = float(np.mean(recalls)) if recalls else None
    return y_pred, acc, bal_acc, confusion


def permutation_test(X, y, classes, shrinkage, n_perm=1000, seed=7):
    """Label permutation: shuffle y and rerun the FULL LOO-CV pipeline."""
    _, acc_obs, bal_obs, _ = loo_cv(X, y, classes, shrinkage)
    rng = np.random.default_rng(seed)
    n = X.shape[0]
    null_acc = np.empty(n_perm)
    null_bal = np.empty(n_perm)
    for b in range(n_perm):
        y_perm = y[rng.permutation(n)]
        _, a, ba, _ = loo_cv(X, y_perm, classes, shrinkage)
        null_acc[b] = a
        null_bal[b] = ba if ba is not None else 0.0
    return {
        "accuracy_observed": acc_obs, "balanced_accuracy_observed": bal_obs,
        "accuracy_null_mean": float(null_acc.mean()),
        "accuracy_null_p95": float(np.percentile(null_acc, 95)),
        "balanced_accuracy_null_mean": float(null_bal.mean()),
        "balanced_accuracy_null_p95": float(np.percentile(null_bal, 95)),
        "p_accuracy": float((1 + (null_acc >= acc_obs).sum()) / (n_perm + 1)),
        "p_balanced_accuracy": float((1 + (null_bal >= bal_obs).sum()) / (n_perm + 1)),
        "n_perm": n_perm,
    }


# ---------------------------------------------------------------------------
# Fisher LDA direction (2 classes, pooled covariance, shrinkage) — port of the
# pairwise/triplet CV: OOF scores, SDS, permutation p, paired bootstrap
# ---------------------------------------------------------------------------
def fit_linear_score(X_train, y_train, class_a, class_b, shrinkage):
    """``w = Sigma_shrunk^-1 (mu_b - mu_a)`` on standardised features."""
    mean = X_train.mean(axis=0)
    std = X_train.std(axis=0)
    std[std == 0] = 1.0
    Xs = (X_train - mean) / std
    p = Xs.shape[1]

    Xa, Xb = Xs[y_train == class_a], Xs[y_train == class_b]
    mu_a = Xa.mean(axis=0) if Xa.shape[0] else np.zeros(p)
    mu_b = Xb.mean(axis=0) if Xb.shape[0] else np.zeros(p)

    pooled_cov, n_dof = np.zeros((p, p)), 0
    for Xc, mu in ((Xa, mu_a), (Xb, mu_b)):
        if Xc.shape[0] > 1:
            diff = Xc - mu
            pooled_cov += diff.T @ diff
            n_dof += Xc.shape[0] - 1
    pooled_cov /= max(n_dof, 1)
    diag_cov = np.diag(np.diag(pooled_cov))
    cov = (1.0 - shrinkage) * pooled_cov + shrinkage * diag_cov
    cov += 1e-9 * np.eye(p)
    return mean, std, np.linalg.solve(cov, mu_b - mu_a)


def score_row(x_row, mean, std, w):
    return float(((x_row - mean) / std) @ w)


def oof_scores(Xp, yp, class_a, class_b, shrinkage):
    """Leave-one-out: each q_t fitted only on the other n-1 rows."""
    n = Xp.shape[0]
    q = np.empty(n)
    for i in range(n):
        mask = np.ones(n, dtype=bool)
        mask[i] = False
        mean, std, w = fit_linear_score(Xp[mask], yp[mask], class_a, class_b, shrinkage)
        q[i] = score_row(Xp[i], mean, std, w)
    return q


def delta_sds(q, yp, class_a, class_b):
    ref = list(q[yp == class_a])
    grp = list(q[yp == class_b])
    d = cliffs_delta(ref, grp)
    sds = None if d is None else 10.0 * abs(d)
    return d, sds, ref, grp


def permutation_p(Xp, yp, class_a, class_b, shrinkage, observed_sds, n_perm, rng):
    n = len(yp)
    null = np.empty(n_perm)
    for b in range(n_perm):
        yp_perm = yp[rng.permutation(n)]
        q_perm = oof_scores(Xp, yp_perm, class_a, class_b, shrinkage)
        _, sds_perm, _, _ = delta_sds(q_perm, yp_perm, class_a, class_b)
        null[b] = sds_perm if sds_perm is not None else 0.0
    p = float((1 + (null >= observed_sds).sum()) / (n_perm + 1))
    return p, float(null.mean()), float(np.percentile(null, 95))


def paired_bootstrap_diff(q_x, q_base, yp, class_a, class_b, n_boot, rng):
    """Paired bootstrap of SDS(model_x) - SDS(baseline) on the same indices."""
    idx_a = np.where(yp == class_a)[0]
    idx_b = np.where(yp == class_b)[0]
    diffs = np.empty(n_boot)
    for b in range(n_boot):
        sa = rng.choice(idx_a, size=len(idx_a), replace=True)
        sb = rng.choice(idx_b, size=len(idx_b), replace=True)
        d_x = cliffs_delta(list(q_x[sa]), list(q_x[sb]))
        d_base = cliffs_delta(list(q_base[sa]), list(q_base[sb]))
        sds_x = 10.0 * abs(d_x) if d_x is not None else 0.0
        sds_base = 10.0 * abs(d_base) if d_base is not None else 0.0
        diffs[b] = sds_x - sds_base
    lo, hi = float(np.percentile(diffs, 2.5)), float(np.percentile(diffs, 97.5))
    return {"diff_mean": float(diffs.mean()), "ci95": [lo, hi],
            "ci_excludes_zero": bool(lo * hi > 0)}





# ---------------------------------------------------------------------------
# Espoo layer 1 — the per-orbit observability block
# ---------------------------------------------------------------------------
# Verbatim port of ``espoo_mast_analysis/espoo_mast_observability.py``
# (``group_stats``, ``classify``, ``CLASS_TEXT``, ``_reference_for``,
# ``per_group``, ``compare_orbit``). It lets the figure script recompute the
# committed ``per_orbit`` / ``per_pass`` / ``orbit_test`` / ``verdict`` blocks of
# ``data/espoo/reference/espoo_mast_observability.json`` from the channel table
# and abort on a single deviating number.
COHERENCE_PEAK_FRAC = 0.30
COHERENCE_MEDIAN_MULT = 5.0
COHERENCE_MIN_MASKED = 2
PRED_MIN_COHERENCE = 0.15
# LUMO reference medians — documentation constants of the upstream module, kept
# here unchanged so the classification reproduces field for field.
LUMO_HEALTHY_G2_MEDIAN = 0.0379
LUMO_HEALTHY_NMASKED_MEDIAN = 19
LUMO_ASC_G2_MEDIAN = 0.0617
LUMO_DESC_G2_MEDIAN = 0.0261
LUMO_CLUTTER_FLOOR = 0.04


def _num(v):
    """Coerce a JSON value to float, or None when absent/non-numeric."""
    try:
        return None if v is None else float(v)
    except (TypeError, ValueError):
        return None


def _pct(vals_, q):
    return float(np.percentile(vals_, q)) if vals_ else None


def group_stats(rows):
    """gamma^2 / mask-size statistics for one group of acquisitions.

    Field-for-field port of ``espoo_mast_observability.group_stats``; the
    ``reference_g2`` key is attached afterwards by ``reference_for`` exactly as
    upstream, where ``per_group`` calls ``_reference_for`` on the result.
    """
    g2 = [_num(r.get("gamma2")) for r in rows]
    g2 = [v for v in g2 if v is not None]
    px = [_num(r.get("A")) for r in rows]
    px = [v for v in px if v is not None]
    if not g2:
        return None
    above_floor = sum(1 for v in g2 if v >= PRED_MIN_COHERENCE)
    return {
        "n": len(rows),
        "n_with_gamma2": len(g2),
        "gamma2_median": float(np.median(g2)),
        "gamma2_p25": _pct(g2, 25),
        "gamma2_p75": _pct(g2, 75),
        "gamma2_min": float(min(g2)),
        "gamma2_max": float(max(g2)),
        "n_masked_median": float(np.median(px)) if px else None,
        "n_masked_max": float(max(px)) if px else None,
        "share_above_pred_floor": above_floor / len(g2),
        "share_mask_below_min": (sum(1 for v in px if v < COHERENCE_MIN_MASKED) / len(px)
                                 if px else None),
    }


CLASS_TEXT = {
    "compact_echo_strong":
        "compact coherent echo, above the LUMO reference scale and the 0.15 "
        "pipeline floor — usable",
    "compact_echo_moderate":
        "compact coherent echo, above the LUMO reference scale but below the "
        "0.15 pipeline floor",
    "coherent_but_clutter_masked":
        "coherent, but the mask is large — coherence plausibly carried by "
        "built-up clutter",
    "below_reference": "below the LUMO reference scale for this geometry",
    "insufficient_data": "insufficient data",
}



def classify(stat):
    """Observability class from the (gamma^2, mask size) pair — upstream verbatim.

    Reference-relative, because ``PRED_MIN_COHERENCE`` (0.15) is the pipeline's
    bar for *using* a sample, not a bar for whether the mast is a scatterer at
    all. Each group carries ``reference_g2`` — the matching LUMO median for its
    orbit.
    """
    g2m = stat.get("gamma2_median")
    if not stat or g2m is None:
        return "insufficient_data"
    pxm = stat.get("n_masked_median")
    ref = stat.get("reference_g2", LUMO_HEALTHY_G2_MEDIAN)
    compact = pxm is not None and pxm <= 25.0
    if g2m < ref:
        return "below_reference"
    if not compact:
        return "coherent_but_clutter_masked"
    return "compact_echo_strong" if g2m >= PRED_MIN_COHERENCE else "compact_echo_moderate"


def reference_for(key, stat):
    """Attach the matching LUMO reference median to a group's statistics."""
    if not stat:
        return stat
    k = key.upper()
    if k.startswith("DESC") or k.startswith("MORNING"):
        stat["reference_g2"] = LUMO_DESC_G2_MEDIAN
    elif k.startswith("ASC") or k.startswith("AFTERNOON"):
        stat["reference_g2"] = LUMO_ASC_G2_MEDIAN
    else:
        stat["reference_g2"] = LUMO_HEALTHY_G2_MEDIAN
    return stat


def per_group(rows):
    """Group the acquisitions by pass label and orbit direction (upstream)."""
    by_pass = collections.defaultdict(list)
    by_orbit = collections.defaultdict(list)
    for r in rows:
        by_pass[r.get("pass") or "?"].append(r)
        by_orbit[r.get("state") or "?"].append(r)
    return (
        {k: reference_for(k, group_stats(v)) for k, v in sorted(by_pass.items())},
        {k: reference_for(k, group_stats(v)) for k, v in sorted(by_orbit.items())},
    )


def compare_orbit(a, b):
    """Welch t + Mann-Whitney between two groups (gamma^2) — upstream verbatim."""
    if sps is None or not a or not b:
        return None
    if len(a) < 3 or len(b) < 3:
        return None
    try:
        t, p_t = sps.ttest_ind(a, b, equal_var=False)
        _, p_mw = sps.mannwhitneyu(a, b, alternative="two-sided")
    except Exception:  # pragma: no cover - degenerate inputs
        return None
    return {"n_a": len(a), "n_b": len(b), "welch_t": float(t),
            "welch_p": float(p_t), "mannwhitney_p": float(p_mw)}


def orbit_block(rows):
    """The whole upstream observability layer in one dict.

    Compared field by field against the committed ``per_pass``, ``per_orbit``,
    ``orbit_test`` and ``verdict.classes`` blocks.
    """
    by_pass, by_orbit = per_group(rows)
    a = vals(rows, STATE_ORDER[0], key="gamma2")
    b = vals(rows, STATE_ORDER[1], key="gamma2")
    return {
        "per_pass": by_pass,
        "per_orbit": by_orbit,
        "orbit_test": compare_orbit(a, b),
        "classes": {k: classify(s) for k, s in by_orbit.items() if s},
        "class_text": {k: CLASS_TEXT[classify(s)] for k, s in by_orbit.items() if s},
    }


# ---------------------------------------------------------------------------
# Espoo layer 2 — the stationarity layer
# ---------------------------------------------------------------------------
# Verbatim port of ``espoo_mast_analysis/espoo_phase_stationarity.py``
# (``mann_kendall``, ``theil_sen``, ``block_bootstrap_slope``, ``lag1_autocorr``,
# ``runs_test``, ``ljung_box``, ``annual_harmonic``, ``split_half``,
# ``monthly_means``) plus the flag rule of ``group_stationarity``.
#
# The upstream script computes these on ``espoo_archive_series.json`` (30 MB, a
# machine-local pipeline artefact that is **not** committed). What *is* committed
# is ``espoo_phase_stationarity.json``, and it stores the full per-acquisition
# series of every one of its 14 groups under ``series_by_group`` as
# ``[[day_ordinal, coherence], ...]``. Every statistical number of that file is
# therefore recomputed here from the committed series and pinned field by field;
# only the keys that the file does not carry (``dwell_s``, ``n_sub``,
# ``nyquist_hz``, ``n_observable``, ``dwell_control``) are not re-derived and are
# listed as such in the report.
def mann_kendall(y):
    """Mann-Kendall trend test with tie correction. Returns tau_b, S, p (two-sided)."""
    y = np.asarray(y, float)
    n = y.size
    if n < 4:
        return {"n": n}
    s = 0.0
    for k in range(n - 1):
        s += float(np.sign(y[k + 1:] - y[k]).sum())
    _, counts = np.unique(y, return_counts=True)
    n1 = n * (n - 1) / 2.0
    n2 = sum(c * (c - 1) * (2 * c + 5) for c in counts)
    var = (n * (n - 1) * (2 * n + 5) - n2) / 18.0
    if var <= 0:
        return {"n": n, "S": s, "tau": None, "p": None}
    z = (s - np.sign(s)) / np.sqrt(var)
    p = float(2 * (1 - 0.5 * (1 + erf(abs(z) / sqrt(2)))))
    pairs = sum(c * (c - 1) / 2.0 for c in counts)
    denom = np.sqrt((n1 - pairs) * (n1 - pairs))
    tau = None if denom == 0 else s / denom
    return {"n": n, "S": s, "tau": tau, "z": float(z), "p": p, "tie_pairs": float(pairs)}


def theil_sen(x, y):
    """Median of pairwise slopes — the trend the outliers cannot move."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    n = x.size
    if n < 3:
        return None
    slopes = [(y[j] - y[i]) / (x[j] - x[i])
              for i in range(n - 1) for j in range(i + 1, n) if x[j] != x[i]]
    return float(np.median(slopes)) if slopes else None


def block_bootstrap_slope(x, y, block=5, n=2000, seed=13, years=True):
    """Moving-block bootstrap CI of the Theil-Sen slope (blocks keep the autocorrelation)."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = x.size
    if m < 2 * block or m < 6:
        return None
    rng = np.random.default_rng(seed)
    nblocks = int(np.ceil(m / block))
    out = []
    for _ in range(n):
        starts = rng.integers(0, m - block + 1, nblocks)
        idx = np.concatenate([np.arange(s, s + block) for s in starts])[:m]
        s_ = theil_sen(x[idx], y[idx])
        if s_ is not None:
            out.append(s_ * (365.25 if years else 1.0))
    if not out:
        return None
    return [float(np.percentile(out, 2.5)), float(np.percentile(out, 97.5))]


def lag1_autocorr(y):
    y = np.asarray(y, float)
    if y.size < 4:
        return None
    a, b = y[:-1], y[1:]
    if a.std() == 0 or b.std() == 0:
        return None
    return float(np.corrcoef(a, b)[0, 1])


def runs_test(y):
    """Wald-Wolfowitz runs test around the median (p two-sided, normal approx)."""
    y = np.asarray(y, float)
    med = np.median(y)
    signs = np.sign(y - med)
    signs = signs[signs != 0]
    n = signs.size
    if n < 6:
        return {"n": n}
    n1 = int((signs > 0).sum())
    n2 = n - n1
    if n1 < 2 or n2 < 2:
        return {"n": n, "n_pos": n1, "n_neg": n2}
    runs = 1 + int((signs[1:] != signs[:-1]).sum())
    mu = 2.0 * n1 * n2 / n + 1
    var = (2.0 * n1 * n2 * (2.0 * n1 * n2 - n)) / (n * n * (n - 1))
    if var <= 0:
        return {"n": n, "runs": runs}
    z = (runs - mu) / np.sqrt(var)
    p = float(2 * (1 - sps.norm.cdf(abs(z)))) if sps is not None else None
    return {"n": n, "n_pos": n1, "n_neg": n2, "runs": runs, "z": float(z), "p": p}


def ljung_box(y, lags=5):
    """Ljung-Box Q on the (already detrended) series — is anything left correlated?"""
    y = np.asarray(y, float)
    m = y.size
    if m < lags + 3:
        return {"n": m}
    y = y - y.mean()
    denom = float((y * y).sum())
    if denom == 0:
        return {"n": m}
    q = 0.0
    for k in range(1, lags + 1):
        r = float((y[:-k] * y[k:]).sum() / denom)
        q += r * r / (m - k)
    q *= m * (m + 2)
    p = float(1 - sps.chi2.cdf(q, lags)) if sps is not None else None
    return {"n": m, "lags": lags, "Q": float(q), "p": p}


def annual_harmonic(x, y):
    """Least-squares annual harmonic: amplitude and day-of-year of the maximum."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    if x.size < 12:
        return None
    t = (x - x.mean()) / 365.25
    w = 2 * np.pi * t
    design = np.column_stack([np.ones_like(t), np.cos(w), np.sin(w)])
    try:
        coef, *_ = np.linalg.lstsq(design, y, rcond=None)
    except np.linalg.LinAlgError:  # pragma: no cover
        return None
    amp = float(np.hypot(coef[1], coef[2]))
    phase = float(np.arctan2(coef[2], coef[1]))
    doy = (phase / (2 * np.pi) * 365.25) % 365.25
    resid = y - design @ coef
    ss_tot = float(((y - y.mean()) ** 2).sum())
    r2 = 1 - float((resid ** 2).sum()) / ss_tot if ss_tot > 0 else None
    return {"amplitude": amp, "peak_day_of_year": doy, "r2": r2,
            "residual_std": float(resid.std(ddof=1)) if y.size > 3 else None}


def _as_list(pair):
    """JSON needs lists, the bootstrap returns tuples — same two numbers."""
    return None if pair is None else [float(v) for v in pair]


def split_half_series(series, split_day):
    """``split_half`` on a committed ``[[day, value], ...]`` series."""
    first = [v for d, v in series if d < split_day and v is not None]
    second = [v for d, v in series if d >= split_day and v is not None]
    if len(first) < 5 or len(second) < 5:
        return {"n": [len(first), len(second)], "delta": None}
    d = cliffs_delta(first, second)
    return {
        "n": [len(first), len(second)],
        "median": [float(np.median(first)), float(np.median(second))],
        "delta": d,
        "delta_ci95": _as_list(bootstrap_delta_ci(first, second, n=2000, seed=7)),
        "sds": None if d is None else 10.0 * abs(d),
    }


def monthly_means_series(series):
    """``monthly_means`` on a committed ``[[day, value], ...]`` series."""
    by_month = collections.defaultdict(list)
    for d, v in series:
        if v is not None:
            by_month[dt.date.fromordinal(int(d)).isoformat()[:7]].append(v)
    return {m: float(np.mean(v)) for m, v in sorted(by_month.items()) if len(v) >= 2}


def stationarity_flags(blk):
    """The flag rule of ``group_stationarity``, applied to one group's numbers."""
    flags = []
    if blk["slope_excludes_zero"] and blk["mann_kendall"].get("p") is not None \
            and blk["mann_kendall"]["p"] < 0.05:
        flags.append("trend")
    sh = (blk.get("split_half") or {}).get("delta_ci95")
    if sh and sh[0] * sh[1] > 0:
        flags.append("split_shift")
    l1 = blk["lag1_autocorr"]
    if l1 is not None and abs(l1) > 0.3:
        flags.append("serial_correlation")
    if blk["runs_test"].get("p") is not None and blk["runs_test"]["p"] < 0.05:
        flags.append("non_random_order")
    if blk["ljung_box"].get("p") is not None and blk["ljung_box"]["p"] < 0.05:
        flags.append("residual_autocorrelation")
    h = blk["seasonal_harmonic"]
    if h and h.get("r2") is not None and h["r2"] > 0.15:
        flags.append("seasonal")
    return flags


def stationarity_of(series, block=5, n_boot=2000):
    """Recompute the committed stationarity block of **one** phase-scan group.

    ``series`` is the committed ``series_by_group[key]`` list of
    ``[day_ordinal, coherence]`` pairs. Returns the same keys the committed file
    carries for that group (minus the source-file-only ones), so the figure
    script can pin them field by field.
    """
    x = np.asarray([d for d, _ in series], float)
    y = np.asarray([v for _, v in series], float)
    ts = theil_sen(x, y)
    slope = None if ts is None else ts * 365.25
    ci = block_bootstrap_slope(x, y, block=block, n=n_boot)
    detrended = y - (x - x.mean()) * (ts if ts is not None else 0.0)
    split_day = int(np.median(x))
    out = {
        "metric": "coherence",
        "n": int(y.size),
        "date_range": [dt.date.fromordinal(int(x.min())).isoformat(),
                       dt.date.fromordinal(int(x.max())).isoformat()],
        "median": float(np.median(y)),
        "iqr": [float(np.percentile(y, 25)), float(np.percentile(y, 75))],
        "theil_sen_per_year": slope,
        "slope_ci95_per_year": ci,
        "slope_excludes_zero": bool(ci and ci[0] * ci[1] > 0),
        "mann_kendall": mann_kendall(y),
        "lag1_autocorr": lag1_autocorr(y),
        "lag1_detrended": lag1_autocorr(detrended),
        "runs_test": runs_test(y),
        "ljung_box": ljung_box(detrended, lags=5),
        "split_half": split_half_series(series, split_day),
        "split_date": dt.date.fromordinal(split_day).isoformat(),
        "seasonal_harmonic": annual_harmonic(x, y),
        "monthly_means": monthly_means_series(series),
    }
    out["flags"] = stationarity_flags(out)
    out["verdict"] = ("stationary" if not out["flags"]
                      else "not stationary: " + ", ".join(out["flags"]))
    return out


# The memo of ``phase_series_block`` (see there). Keyed on the identity of the
# document object; the hit path re-checks that identity before trusting the key.
_PHASE_BLOCK_CACHE = {}


def phase_series_block(doc):
    """Recompute every group of a committed ``espoo_phase_stationarity.json``.

    Returns ``(groups, meta)``: ``groups`` maps the committed group key
    (``"<band> · <config> · <orbit>"``) to the recomputed stationarity block,
    ``meta`` carries the file-level quantities that the series determine
    (n rows, per-orbit counts, date range).

    The result is memoised per document object. ``stationarity_summary``,
    ``flagged_groups`` and ``primary_series`` are all defined on top of this
    function, and recomputing the 14 block-bootstrap series three times in one
    figure run is pure waste (the Theil-Sen estimator dominates the runtime).
    The cache changes no number — it is keyed on the identity of the document,
    which is re-checked on every hit.
    """
    key = id(doc)
    hit = _PHASE_BLOCK_CACHE.get(key)
    if hit is not None and hit[0] is doc:
        return hit[1]

    series_by_group = doc.get("series_by_group") or {}
    params = doc.get("params") or {}
    block = int(params.get("block", 5))
    n_boot = int(params.get("n_boot", 2000))
    groups, rows, n_by_orbit = {}, [], {"ASCENDING": 0, "DESCENDING": 0}
    for key_ in sorted(series_by_group):
        series = series_by_group[key_]
        groups[key_] = stationarity_of(series, block=block, n_boot=n_boot)
        n = groups[key_]["n"]
        rows.append(n)
        orbit = key_.rsplit(" · ", 1)[-1]
        if orbit in n_by_orbit:
            n_by_orbit[orbit] += n
    out = (groups, {
        "n_groups": len(groups),
        "n_rows": int(sum(rows)),
        "n_by_orbit": n_by_orbit,
        "block": block,
        "n_boot": n_boot,
    })
    _PHASE_BLOCK_CACHE[key] = (doc, out)
    return out


def stationarity_summary(doc):
    """``(n_stationary, n_groups)`` of a committed stationarity file, recomputed.

    The committed file states this in prose (``interpretation[0]``: *"13 of 14
    series are stationary by the tests used, 1 carry at least one flag"*); this
    recomputes it from the committed series so both can be pinned against each
    other.
    """
    groups, _meta = phase_series_block(doc)
    n_stat = sum(1 for b in groups.values() if b["verdict"] == "stationary")
    return n_stat, len(groups)


def flagged_groups(doc):
    """The recomputed group keys that carry at least one stationarity flag."""
    groups, _meta = phase_series_block(doc)
    return sorted(k for k, b in groups.items() if b["flags"])


def primary_series(doc):
    """The committed primary series key and its recomputed block.

    ``primary`` is the group the upstream file nominates as the highest-median
    series; the recomputed block is what the month companion draws.
    """
    groups, _meta = phase_series_block(doc)
    key = doc["primary"]
    return groups[key], doc.get("series_by_group", {}).get(key)

