#!/usr/bin/env python3
"""KDLO OCV-Paper — statistics ports (verbatim from the LUMO package).

This module invents *nothing*: every function is a verbatim port of
``code/lumo/lumo_ocv_stats.py``, which is itself a verbatim port of the LUMO
project analysis scripts

  * ``analyze_lumo_tower_coherence.py``            (``_stats``, ``_test``)
  * ``analyze_lumo_gamma2_pairs.py``               (Cliff's delta, Bootstrap-CI,
                                                    ``delta_block``, ``sign_test``,
                                                    ``sds_from``, ``by_state``,
                                                    ``vals``, ``season_of``)
  * ``analyze_lumo_emv_multivariate_cv.py``        (``fit_lda``, ``predict_lda``,
                                                    ``loo_cv``, ``permutation_test``)
  * ``analyze_lumo_emv_pairwise_triplets_cv.py``   (``fit_linear_score``,
                                                    ``score_row``, ``oof_scores``,
                                                    ``delta_sds``, ``permutation_p``,
                                                    ``paired_bootstrap_diff``)

The two KDLO differences are *naming only*:

  * the state field is ``state`` (not ``damage_label``) — a KDLO "state" is an
    **epoch** (a time window), not a damage label; and
  * the ordering is 2-valued and chronological: ``pre`` -> ``rebuild``.

Everything else — the exact formulas, the bootstrap seeds, the rank-based
Cliff's delta, the shrinkage LDA, the leave-one-out loops — is unchanged, so a
difference in the KDLO numbers cannot come from a difference in the statistics.

Conventions (identical to the origin):
  SDS  = 10 * |Cliff's delta|
  delta(ref, grp) > 0  =>  the *group* is larger than the reference
  bootstrap seeds 7 (in-sample) resp. 11/12/13 (OOF), permutation seed 7/8
"""
from __future__ import annotations

import numpy as np

try:
    from scipy import stats as sps
except ImportError:  # pragma: no cover - statistical tests are optional
    sps = None

# Chronological epoch order. `SEVERITY` keeps the LUMO name (it is the monotone
# axis used by the strata/trend helpers), but at KDLO it counts *time*, not
# damage: 0 = before the collapse, 1 = during the rebuild. No causal claim.
STATE_ORDER = ["pre", "rebuild"]
STATE_SHORT = {"pre": "pre", "rebuild": "rebuild"}
SEVERITY = {"pre": 0, "rebuild": 1}

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
    vals_ = [abs(d) for d in deltas if d is not None]
    return 10.0 * max(vals_) if vals_ else None


def by_state(rows):
    """Rows per epoch state (``state`` field; the ``single`` bucket is not one)."""
    return {lab: [r for r in rows if r.get("state") == lab] for lab in STATE_ORDER}


def vals(rows, label, key="gamma2"):
    """Values of one state — accepts a row list OR a ``by_state`` dict."""
    if isinstance(rows, dict):
        recs = rows.get(label) or []
    else:
        recs = [r for r in rows if r.get("state") == label]
    return [r[key] for r in recs if r.get(key) is not None]


def season_of(r):
    """Season from the ``month`` field (1..12) — verbatim ``gp.season_of``."""
    m = r.get("month")
    if m is None:
        return None
    try:
        m = int(m)
    except (TypeError, ValueError):
        return None
    if m in WINTER_MONTHS:
        return "winter_Oct-Feb"
    if m in SUMMER_MONTHS:
        return "summer_Mar-Jul"
    return "shoulder_Aug-Sep"


def strata_pair(rows, key_fn, ref="pre", grp="rebuild", n_boot=BOOT_STRATA,
                key="gamma2", min_n=3):
    """One ``delta_block`` per stratum — the KDLO form of ``gp._strata``.

    LUMO stratified ``healthy vs. one damage label``; KDLO has a single pair, so
    a stratum reports ``pre vs. rebuild`` and is marked ``n < 3`` instead of
    being silently dropped as soon as one side has fewer than ``min_n`` rows.
    """
    out = {}
    for skey in sorted({key_fn(r) for r in rows if key_fn(r) is not None},
                       key=lambda s: str(s)):
        sub = [r for r in rows if key_fn(r) == skey]
        a = vals(sub, ref, key)
        b = vals(sub, grp, key)
        entry = {"n": [len(a), len(b)]}
        if len(a) < min_n or len(b) < min_n:
            entry.update({"delta": None, "note": f"n < {min_n}"})
        else:
            entry = delta_block(a, b, f"{ref}_vs_{grp}_in_{skey}", n_boot)
        out[str(skey)] = entry
    return out


# ---------------------------------------------------------------------------
# LDA (shrinkage-regularised pooled covariance) — port of mcv
# ---------------------------------------------------------------------------
def fit_lda(X_train, y_train, classes, shrinkage):
    """Uniform priors (deliberate, see ``mcv.fit_lda``): otherwise the model tips
    towards the majority state."""
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

