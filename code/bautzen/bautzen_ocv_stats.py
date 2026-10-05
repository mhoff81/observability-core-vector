#!/usr/bin/env python3
"""OCV-Paper (Bautzen bridge site) — statistics layer.

Estimator layer: lines 47 ff. of ``lumo_ocv_stats.py`` are copied **unchanged**,
so every estimator (Cliff's delta, SDS, bootstrap CIs, LDA/LOO-CV, permutation
tests, Jonckheere-Terpstra, strata) is bit-identical to the LUMO mast.

Only two things differ, both deliberate:

  1. the state alphabet of this site is ``RS / DS1 / DS2`` (reference, damage
     state 1, damage state 2) instead of healthy / DAM3 / DAM4 / DAM6, and the
     severity axis is the ordinal position of the state;
  2. the Bautzen-specific estimators in the second half of the file
     (``wilson`` … ``env_fit``) are verbatim ports of
     ``analyze_bautzen_deck_edge_vs_state.py`` — the estimators the committed
     reference ``bautzen_deck_edge_state.json`` was computed with (same seeds:
     permutation 7, bootstrap 42, 20000 permutations, 5000 draws).

Conventions (identical to the origin):
  SDS  = 10 * |Cliff's delta|
  delta(ref, grp) > 0  =>  the *group* is larger than the reference
  bootstrap seeds 7 (in-sample) resp. 11/12/13 (OOF), permutation seed 7/8
  ``step_test`` uses the permutation p **without** the +1 correction (as in the
  Bautzen origin), the LDA tests keep the +1 correction (as in the LUMO origin)
"""
from __future__ import annotations

import datetime as dt
import math

import numpy as np

try:
    from scipy import stats as sps
except ImportError:  # pragma: no cover - statistical tests are optional
    sps = None

# ── State alphabet of this site ─────────────────────────────────────────────
STATE_ORDER = ["RS", "DS1", "DS2"]
STATE_SHORT = {"RS": "RS", "DS1": "DS1", "DS2": "DS2"}
DAMAGED = ["DS1", "DS2"]
SEVERITY = {"RS": 0, "DS1": 1, "DS2": 2}
REFERENCE_STATE = "RS"      # the pre-intervention state (= LUMO's "healthy")

# ── Seeds and draw counts (identical to the Bautzen origin) ─────────────────
N_PERM = 20000
N_BOOT = 5000
RNG_SEED = 7
BOOT_SEED = 42
PLACEBO_WEEKS = (2, 4, 6, 8, 10, 12)

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
    return {lab: [r for r in rows if r["damage_label"] == lab] for lab in STATE_ORDER}


def vals(rows, label, key="gamma2"):
    """Values of one state — accepts a row list OR a ``by_state`` dict."""
    if isinstance(rows, dict):
        recs = rows.get(label) or []
    else:
        recs = [r for r in rows if r["damage_label"] == label]
    return [r[key] for r in recs if r.get(key) is not None]


def pairs_block(rows, n_boot=BOOT_MAIN, key="gamma2"):
    """Reference state vs. each damage state + every state pair in one dict.

    The LUMO original calls this "healthy_vs"; this site's pre-intervention
    reference state is ``RS`` (see REFERENCE_STATE), so the keys are
    ``reference_vs_*``. Everything else is unchanged.
    """
    st = by_state(rows)
    out = {"channel": key,
           "n_by_state": {STATE_SHORT[l]: len(st[l]) for l in STATE_ORDER},
           "reference_vs": {}, "all_pairs": {}}
    ds = []
    for lab in DAMAGED:
        if not st[lab]:
            continue
        blk = delta_block(vals(st, REFERENCE_STATE, key), vals(st, lab, key),
                          f"reference_vs_{STATE_SHORT[lab]}", n_boot)
        out["reference_vs"][STATE_SHORT[lab]] = blk
        ds.append((STATE_SHORT[lab], blk.get("delta")))
    out["sign_decomposition"] = dict(ds)
    labels = STATE_ORDER
    for i, a in enumerate(labels):
        for b in labels[i + 1:]:
            blk = delta_block(vals(st, a, key), vals(st, b, key),
                              f"{STATE_SHORT[a]}_vs_{STATE_SHORT[b]}", n_boot)
            out["all_pairs"][blk["comparison"]] = blk
    pooled = delta_block(vals(st, REFERENCE_STATE, key),
                         [v for lab in DAMAGED for v in vals(st, lab, key)],
                         "reference_vs_all_damaged", n_boot)
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
    d = {k: (pair_block["reference_vs"].get(k) or {}).get("delta")
         for k in ("DS1", "DS2")}
    absd = {k: (None if v is None else abs(v)) for k, v in d.items()}
    hy = {"channel": key,
          "claim": "|delta(RS,DS1)| < |delta(RS,DS2)|",
          "delta": d, "abs_delta": absd}
    if all(v is not None for v in absd.values()):
        hy["holds"] = bool(absd["DS1"] < absd["DS2"])
        hy["order_observed"] = " < ".join(
            k for k, _v in sorted(absd.items(), key=lambda kv: kv[1]))
    out = {"hypothesis": hy}
    if sps is not None:
        sev = np.asarray([SEVERITY[r["damage_label"]] for r in rows], float)
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
            if len(st[REFERENCE_STATE]) < 3 or len(st[lab]) < 3:
                entry[STATE_SHORT[lab]] = {
                    "n": [len(st[REFERENCE_STATE]), len(st[lab])],
                    "delta": None, "note": "n < 3"}
                continue
            entry[STATE_SHORT[lab]] = delta_block(
                vals(st, REFERENCE_STATE, key), vals(st, lab, key),
                f"reference_vs_{STATE_SHORT[lab]}_in_{skey}", n_boot)
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
# Bautzen-specific estimators — verbatim ports of
# analyze_bautzen_deck_edge_vs_state.py (same formulas, same seeds, same
# defaults). Used by the committed reference bautzen_deck_edge_state.json.
# ---------------------------------------------------------------------------
def wilson(k, n, z=1.96):
    """Share + Wilson 95% CI (small n!)."""
    if n == 0:
        return float("nan"), float("nan"), float("nan")
    p = k / n
    den = 1.0 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return p, max(0.0, c - h), min(1.0, c + h)


def median_iqr(v):
    """Median + 25/75 percentiles; (nan, nan, nan) for an empty list."""
    if len(v) == 0:
        return float("nan"), float("nan"), float("nan")
    a = np.asarray(v, float)
    q = np.percentile(a, [25, 75])
    return float(np.median(a)), float(q[0]), float(q[1])


def thirds_drift(v):
    """Median(last third) - median(first third); nan for n < 3."""
    n = len(v)
    if n < 3:
        return float("nan")
    k = max(1, n // 3)
    return float(np.median(v[-k:]) - np.median(v[:k]))


def ols_slope(v):
    """Slope [unit per date] of a straight line through the time series."""
    n = len(v)
    if n < 3:
        return float("nan")
    return float(np.polyfit(np.arange(n, dtype=float), np.asarray(v, float), 1)[0])


def jaccard(a, b):
    """Jaccard overlap of two pixel sets (list of (row, col, freq))."""
    sa = {(r, c) for r, c, _ in a}
    sb = {(r, c) for r, c, _ in b}
    if not sa and not sb:
        return float("nan")
    return len(sa & sb) / len(sa | sb)


def boot_median_ci(v, n_boot=N_BOOT, seed=BOOT_SEED):
    """Median + percentile bootstrap 95% CI; (nan, nan, nan) for an empty list."""
    a = np.asarray(v, float)
    if a.size == 0:
        return float("nan"), float("nan"), float("nan")
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, a.size, size=(n_boot, a.size))
    boots = np.median(a[idx], axis=1)
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return float(np.median(a)), float(lo), float(hi)


def env_at(env, when):
    """Environment values at ``when`` (linear interpolation; None on a gap).

    ``env`` is the ``(timestamps, {variable: array})`` tuple built by
    ``bautzen_ocv_core.load_meteo()``.
    """
    if env is None:
        return None
    ts, cols = env
    if when < ts[0] or when > ts[-1]:
        return None
    x = np.asarray([(t - ts[0]).total_seconds() for t in ts], float)
    xq = (when - ts[0]).total_seconds()
    return {k: float(np.interp(xq, x, v)) for k, v in cols.items()}


def env_fit(x, y):
    """OLS y = a + b*x with r, r^2, p (two-sided) and se(b).

    Estimation logic identical to ``env_fit`` in
    ``analyze_bautzen_deck_edge_vs_state.py`` (and to ``linfit`` in
    ``analyze_bautzen_phase_frequency_proxy.py``).
    """
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    n = x.size
    if n < 4 or np.std(x) == 0:
        return None
    b, a = np.polyfit(x, y, 1)
    resid = y - (a + b * x)
    ss_res = float((resid ** 2).sum())
    r = float(np.corrcoef(x, y)[0, 1])
    se = float(np.sqrt(ss_res / (n - 2) / ((x - x.mean()) ** 2).sum()))
    try:
        from scipy import stats
        p = float(2 * stats.t.sf(abs(b) / se, df=n - 2)) if se > 0 else float("nan")
    except Exception:  # noqa: BLE001
        p = float("nan")
    return {"slope": float(b), "intercept": float(a), "r": r, "r2": r * r,
            "p": p, "n": int(n), "se_slope": se}


def step_test(keys, values, cut, n_perm=N_PERM, seed=RNG_SEED):
    """Median(after) - median(before) + permutation p (two-sided).

    ``keys`` are acquisition datetimes, ``cut`` a datetime. The permutation p
    carries no +1 correction (verbatim ``step_test`` of the Bautzen origin).
    """
    v = np.asarray(values, float)
    is_before = np.asarray([k < cut for k in keys])
    before, after = v[is_before], v[~is_before]
    if len(before) == 0 or len(after) == 0:
        return None
    obs = float(np.median(after) - np.median(before))
    pool = np.concatenate([before, after])
    rng = np.random.default_rng(seed)
    nb = len(before)
    cnt = 0
    for _ in range(n_perm):
        rng.shuffle(pool)
        cnt += abs(float(np.median(pool[nb:]) - np.median(pool[:nb]))) >= abs(obs)
    return {"n_before": int(nb), "n_after": int(len(after)),
            "median_before": float(np.median(before)),
            "median_after": float(np.median(after)),
            "delta": obs, "p_perm": cnt / n_perm}


def placebo_sweep(keys, values, cut):
    """The real cut against cuts shifted by +/-2..12 weeks.

    If the real cut is not the most extreme one, the effect is season.
    """
    real = step_test(keys, values, cut)
    if real is None:
        return None
    rows = []
    for w in PLACEBO_WEEKS:
        for sign in (-1, 1):
            c = cut + dt.timedelta(weeks=sign * w)
            if c <= min(keys) or c > max(keys):
                continue
            st = step_test(keys, values, c)
            if st is None:
                continue
            rows.append({"cut": c.date().isoformat(), "offset_weeks": sign * w,
                         "delta": st["delta"], "p_perm": st["p_perm"]})
    deltas = [abs(r["delta"]) for r in rows] + [abs(real["delta"])]
    return {"real": dict(cut=cut.date().isoformat(), **real),
            "placebo": rows,
            "real_is_most_extreme": bool(abs(real["delta"]) >= max(deltas) - 1e-12),
            "n_placebo": len(rows),
            "max_abs_placebo_delta": max((abs(r["delta"]) for r in rows),
                                         default=float("nan"))}
