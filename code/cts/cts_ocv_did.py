#!/usr/bin/env python3
"""CTS OCV-Paper — the difference-in-differences layer.

This is the layer the CTS site makes possible and a single-structure site cannot
have. The four Sentinel-1 footprints of the record

    cts_asc_48_iw3    the tower            role ``target``
    ctn_asc_48_iw3    Champlain Towers N   role ``control_ctn``  (165 m N)
    cte_asc_48_iw3    Champlain Towers E   role ``control_cte``  (~90 m N)
    beach_asc_48_iw3  the beach east       role ``reference_beach``

share the *same* orbit (ASC rel-48 IW3), the *same* dates and the *same* counts,
so a target and a control tower see the same atmosphere, orbit drift and any
processing effect. The difference of differences removes exactly those common
terms:

    y ~ 1 + T + P + T:P        T = is_target, P = is_post
    DiD = coefficient of T:P

The module is a verbatim port of the site's own generator
(``asset_analysis.reference_did`` — same ``cells`` / ``analyse_channel`` / ``run``
/ ``report``) onto the CTS measurement table, so the three committed references

    data/cts/reference/cts_reference_did_ctn.json               target vs. CTN
    data/cts/reference/cts_reference_did_cte.json               target vs. CTE
    data/cts/reference/cts_reference_did_placebo_ctn_vs_cte.json  placebo CTN vs. CTE

can be reproduced and pinned. The placebo (CTN vs. CTE, both controls) is the
*specificity* check: neither control saw an event, so a significant placebo DiD
would mean the contrast is not specific to the collapse.

Three things are reported beside the DiD, because a DiD without them is not
interpretable:

* **pre-period comparability** — target vs reference before the event (Welch):
  two series on different scales make the DiD an estimate of nothing;
* **the reference's own step** — if the reference *also* steps at the event date,
  the differencing assumption failed;
* **parallel pre-trends** — the pre-event drift of both series (Spearman rho).

ASCENDING and DESCENDING are never pooled (EQS-9).

One channel is degenerate: ``coherence_masked_pixels`` is the site's *fixed*
640-px quantile mask, constant for every row, so its OLS interaction is pure
floating-point noise. It is pinned structurally (``n`` / ``cells`` / ``did`` /
deltas), never on ``did_t`` / ``did_p``.

Usage
-----
  python3 code/cts/cts_ocv_did.py --verify            # pin against the three refs
  python3 code/cts/cts_ocv_did.py --verify --json      # ... machine-readable
  python3 code/cts/cts_ocv_did.py --run did_ctn        # print one DiD document
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from collections import defaultdict
from datetime import date

import numpy as np

try:
    from scipy import stats as sps
except ImportError:  # pragma: no cover - the tests are optional
    sps = None

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import cts_ocv_core as core  # noqa: E402

EVENT = core.EVENT
DID_CHANNELS = list(core.DID_CHANNELS)
DID_ORBITS = list(core.DID_ORBITS)
DID_CONTRASTS = list(core.DID_CONTRASTS)      # (target_role, reference_role, key)
DEGENERATE = set(core.DEGENERATE_CHANNELS)
REF_DID = dict(core.REF_DID)

EPOCH = "2022-01-01"     # ``io.t_days`` epoch
TOL = 1e-9               # the pin tolerance
METHOD = ("y ~ 1 + T + P + T:P over the pooled target/reference rows; the T:P "
          "coefficient is the difference of differences. Reported with the "
          "pre-period comparability, the reference's own step and both pre-trends, "
          "because a DiD without them is not interpretable.")



# ---------------------------------------------------------------------------
# io ports (asset_analysis.io) — verbatim
# ---------------------------------------------------------------------------
def _f(v):
    """Cell -> float or ``None`` (empty/``null``/unparseable -> ``None``)."""
    if v is None:
        return None
    if isinstance(v, (int, float)):
        f = float(v)
        return f if np.isfinite(f) else None
    v = str(v).strip()
    if v in ("", "null", "NULL", "None", "nan"):
        return None
    try:
        f = float(v)
    except ValueError:
        return None
    return f if np.isfinite(f) else None


def iso_date(ts):
    return (ts or "")[:10]


def day(ts):
    try:
        return date.fromisoformat(iso_date(ts))
    except ValueError:
        return None


def t_days(ts, epoch=EPOCH):
    d, e = day(ts), date.fromisoformat(epoch)
    return float((d - e).days) if d else None


def load_rows(path=core.MEAS_PATH):
    """Header-aware loader -> ``(keys, rows)``; rows as dicts of strings, sorted
    by ``acquisition_ts`` (``io.load_measurements``, unchanged)."""
    with open(path) as fh:
        lines = fh.read().splitlines()
    if not lines:
        return [], []
    keys = lines[0].split("|")
    rows = []
    for ln in lines[1:]:
        p = ln.split("|")
        if len(p) != len(keys):
            continue
        rows.append(dict(zip(keys, p)))
    rows.sort(key=lambda r: r.get("acquisition_ts", ""))
    return keys, rows


def dedup(rows, keys):
    """One observation per ``(date, pass, orbit)`` — median of the numeric keys."""
    groups = defaultdict(list)
    for r in rows:
        groups[(iso_date(r.get("acquisition_ts", "")),
                r.get("pass_label", ""),
                r.get("orbit_direction", ""))].append(r)
    out = []
    for (_d, _p, _o), grp in groups.items():
        grp.sort(key=lambda r: r.get("acquisition_ts", ""))
        rec = dict(grp[0])
        for k in keys:
            vals = [_f(r.get(k)) for r in grp]
            vals = [v for v in vals if v is not None]
            rec[k] = float(np.median(vals)) if vals else ""
        out.append(rec)
    out.sort(key=lambda r: r.get("acquisition_ts", ""))
    return out


def coverage(rows):
    """``n``, date range and month count (``io.coverage``, without a window)."""
    days = sorted({iso_date(r.get("acquisition_ts", "")) for r in rows
                   if r.get("acquisition_ts")})
    months = sorted({d[:7] for d in days})
    return {
        "n": len(rows),
        "days": len(days),
        "first_day": days[0] if days else None,
        "last_day": days[-1] if days else None,
        "months": len(months),
    }


# ---------------------------------------------------------------------------
# stats ports (asset_analysis.stats + season_temperature_models._t_p)
# ---------------------------------------------------------------------------
def welch(a, b):
    """Welch t-test -> ``(t, p)``; ``(None, None)`` when a side has <3 finite
    values or both sides are constant."""
    x = np.asarray([v for v in a if v is not None], float)
    y = np.asarray([v for v in b if v is not None], float)
    if sps is None or x.size < 3 or y.size < 3:
        return None, None
    if x.std(ddof=1) == 0 and y.std(ddof=1) == 0:
        return None, None
    res = sps.ttest_ind(x, y, equal_var=False)
    return float(res.statistic), float(res.pvalue)


def spearman(a, b):
    """Spearman ``(rho, p)``; ``(None, None)`` when <5 pairs or a side is constant."""
    x = np.asarray([v for v in a if v is not None], float)
    y = np.asarray([v for v in b if v is not None], float)
    if sps is None or x.size < 5 or np.all(x == x[0]) or np.all(y == y[0]):
        return None, None
    res = sps.spearmanr(x, y)
    return float(res.statistic), float(res.pvalue)


def ols(design, y):
    """OLS via SVD -> ``(beta, se, dof, rss, r2)``; ``None`` when rank-deficient
    (``stats.ols``, unchanged)."""
    X = np.asarray(design, dtype=float)
    y = np.asarray(y, dtype=float)
    if X.ndim != 2 or X.shape[0] != y.size or X.shape[0] <= X.shape[1]:
        return None
    if not np.all(np.isfinite(X)) or not np.all(np.isfinite(y)):
        return None
    beta, _, rank, _ = np.linalg.lstsq(X, y, rcond=None)
    if rank < X.shape[1]:
        return None
    resid = y - X @ beta
    rss = float(resid @ resid)
    dof = float(X.shape[0] - X.shape[1])
    if dof <= 0 or not np.isfinite(rss):
        return None
    s2 = rss / dof
    try:
        cov = s2 * np.linalg.inv(X.T @ X)
        se = np.sqrt(np.clip(np.diag(cov), 0.0, None))
    except np.linalg.LinAlgError:
        se = None
    tss = float(((y - y.mean()) ** 2).sum())
    r2 = 1.0 - rss / tss if tss > 0 else 0.0
    return beta, se, dof, rss, r2


def _t_p(beta, se, dof):
    """Two-sided p of one coefficient (t distribution) — ``_t_p``, verbatim."""
    if se is None or se <= 0 or dof <= 0:
        return None
    return float(2.0 * sps.t.sf(abs(beta / se), dof))



# ---------------------------------------------------------------------------
# The DiD itself (asset_analysis.reference_did) — verbatim, roles as parameters
# ---------------------------------------------------------------------------
def cells(rows, key, orbit, event, target_role="target",
          reference_role="reference"):
    """Rows of both roles as ``(y, is_target, is_post, date)``.

    The roles are parameters because the choice of reference is a *design*
    decision that has to be re-runnable: the CTS contrasts are target-vs-CTN,
    target-vs-CTE and the CTN-vs-CTE placebo.
    """
    out = []
    for r in rows:
        if r.get("orbit_direction") != orbit:
            continue
        y = _f(r.get(key))
        role = (r.get("role") or "").strip().lower()
        if y is None or role not in (target_role, reference_role):
            continue
        if not np.isfinite(y):
            continue
        out.append({
            "y": float(y),
            "target": 1.0 if role == target_role else 0.0,
            "post": 1.0 if iso_date(r.get("acquisition_ts", "")) >= event else 0.0,
            "date": iso_date(r.get("acquisition_ts", "")),
        })
    return out


def analyse_channel(rows, key, orbit, event, target_role="target",
                    reference_role="reference"):
    """The 2x2 target/post model of one channel in one orbit."""
    data = cells(rows, key, orbit, event, target_role, reference_role)
    pre_t = [d["y"] for d in data if d["target"] and not d["post"]]
    pre_r = [d["y"] for d in data if not d["target"] and not d["post"]]
    post_t = [d["y"] for d in data if d["target"] and d["post"]]
    post_r = [d["y"] for d in data if not d["target"] and d["post"]]
    if min(len(pre_t), len(pre_r), len(post_t), len(post_r)) < 5:
        return {"n": len(data), "note": "eine der vier Zellen hat < 5 Werte",
                "cells": [len(pre_t), len(pre_r), len(post_t), len(post_r)]}

    X = np.array([[1.0, d["target"], d["post"], d["target"] * d["post"]]
                  for d in data])
    y = np.array([d["y"] for d in data], dtype=float)
    fit = ols(X, y)
    if fit is None:
        return {"n": len(data), "note": "rang-defizient"}
    beta, se, dof, _rss, r2 = fit
    if se is None or se[3] <= 0:
        return {"n": len(data), "note": "kein Standardfehler fuer den Interaktionsterm"}

    pre_p = welch(pre_t, pre_r)[1]
    ref_step_p = welch(post_r, pre_r)[1]
    trend_t = spearman(
        [t_days(d["date"] + "T00:00:00Z") for d in data if d["target"] and not d["post"]],
        pre_t,
    )
    trend_r = spearman(
        [t_days(d["date"] + "T00:00:00Z") for d in data if not d["target"] and not d["post"]],
        pre_r,
    )
    return {
        "n": len(data),
        "cells": {"target_pre": len(pre_t), "ref_pre": len(pre_r),
                  "target_post": len(post_t), "ref_post": len(post_r)},
        "did": float(beta[3]),
        "did_t": float(beta[3] / se[3]),
        "did_p": _t_p(beta[3], se[3], dof),
        "r2": float(r2),
        "target_delta": float(np.median(post_t) - np.median(pre_t)),
        "ref_delta": float(np.median(post_r) - np.median(pre_r)),
        "pre_comparable_p": pre_p,
        "reference_steps_p": ref_step_p,
        "pre_trend_target_rho": trend_t[0], "pre_trend_target_p": trend_t[1],
        "pre_trend_reference_rho": trend_r[0], "pre_trend_reference_p": trend_r[1],
    }



# ---------------------------------------------------------------------------
# The document (reference_did.run / document) — one DiD file per contrast
# ---------------------------------------------------------------------------
def _asset_block():
    """The committed asset registry entry plus the DiD ``placebo`` flag.

    The site's registry entry carries ``placebo`` (whether the asset is the
    control-for of another); the committed ``cts_asset.json`` stores the rest of
    the entry, so the flag is derived from ``control_for`` exactly as
    ``asset_analysis``'s ``reference_did`` does — an asset that controls for
    another is the placebo.
    """
    a = dict(core.load_reference(core.REF_ASSET))
    a.setdefault("placebo", bool(a.get("control_for")))
    return a


def document(rows, target_role, reference_role, raw_n, event=EVENT):
    """The full DiD document for one contrast (same shape as the references)."""
    rows_dedup = dedup(rows, list(DID_CHANNELS) + ["temperature_c"])
    series = {
        orbit: {key: analyse_channel(rows_dedup, key, orbit, event,
                                     target_role, reference_role)
                for key in DID_CHANNELS}
        for orbit in DID_ORBITS
    }
    roles_present = sorted({(r.get("role") or "").strip() for r in rows})
    return {
        "asset": _asset_block(),
        "event": event,
        "roles_present": roles_present,
        "target_role": target_role,
        "reference_role": reference_role,
        "method": METHOD,
        "rows": {"rows_raw": raw_n, "rows_dedup": len(rows_dedup)},
        "coverage": coverage(rows_dedup),
        "series": series,
    }


def run(target_role, reference_role, path=core.MEAS_PATH, event=EVENT):
    """Load the measurement table and build one contrast's document."""
    _keys, rows = load_rows(path)
    return document(rows, target_role, reference_role, len(rows), event)


def run_all(path=core.MEAS_PATH, event=EVENT):
    """``{reference key: document}`` for the three committed contrasts."""
    _keys, rows = load_rows(path)
    return {key: document(rows, t, r, len(rows), event)
            for (t, r, key) in DID_CONTRASTS}


# ---------------------------------------------------------------------------
# The markdown report (reference_did.report)
# ---------------------------------------------------------------------------
def _d(x):      # Δ and DiD: four significant digits, signed
    return "%+.4g" % x if x is not None else "—"


def _f2(x):     # a t or a rho: two decimals, signed
    return "%+.2f" % x if x is not None else "—"


def _g3(x):     # a p in the table: three decimals, or an em dash
    return "%.3f" % x if x is not None else "—"


def _n3(x):     # a comparability p: three decimals, or "n/a"
    return "%.3f" % x if x is not None else "n/a"


def _s2(x):     # a rho: two decimals, signed, or "n/a"
    return "%+.2f" % x if x is not None else "n/a"


def _row(chan, a):
    if "note" in a:
        return "| `%s` | %d | — | — | — | — | — | — | — | %s |" % (
            chan, a["n"], a["note"])
    return ("| `%s` | %d | %s | %s | **%s** | %s | %s | %s | %s | %s/%s |"
            % (chan, a["n"], _d(a["target_delta"]), _d(a["ref_delta"]),
               _d(a["did"]), _f2(a["did_t"]), _g3(a["did_p"]),
               _n3(a["pre_comparable_p"]), _n3(a["reference_steps_p"]),
               _s2(a["pre_trend_target_rho"]),
               _s2(a["pre_trend_reference_rho"])))


def report(doc):
    """The committed ``*_did_*.md`` markdown."""
    a = doc["asset"]
    L = ["# Reference vs target (DiD) — %s" % a["label"], "",
         "- event: **%s**" % doc["event"],
         "- comparison: `%s` vs `%s`" % (doc["target_role"], doc["reference_role"]),
         "- roles present: %s" % ", ".join(doc["roles_present"]),
         "- rows: %d raw -> %d observations"
         % (doc["rows"]["rows_raw"], doc["rows"]["rows_dedup"])]
    head = ("| channel | n | target Δ | reference Δ | **DiD** | t | p | "
            "pre comparable p | reference steps p | pre-trend ρ target/ref |")
    rule = "|---|---:|---:|---:|---:|---:|---:|---:|---:|---|"
    for orbit in DID_ORBITS:
        L += ["", "## %s" % orbit, "", head, rule]
        L += [_row(chan, doc["series"][orbit][chan]) for chan in DID_CHANNELS]
    L += [
        "", "## Reading", "",
        "- A significant DiD in a channel whose **reference also steps** "
        "(`reference steps p` small) is not a DiD result — the differencing "
        "assumption failed.",
        "- A significant DiD with a diverging pre-trend or a small "
        "`pre comparable p` is reported as such, not as evidence.",
        "- The DiD removes *common* confounders (atmosphere, soil moisture, "
        "orbits). It cannot remove anything that affects only the bridge — that "
        "is what the temperature and season terms of the model ladder are for.",
        "",
    ]
    return "\n".join(L)



# ---------------------------------------------------------------------------
# The pin (against the three committed references)
# ---------------------------------------------------------------------------
def _cmp(path, ref, got, problems, tol=TOL, skip=()):
    if isinstance(ref, dict):
        if not isinstance(got, dict):
            problems.append((path, "type", "dict", type(got).__name__)); return
        for k, rv in ref.items():
            if k in skip:
                continue
            if k not in got:
                problems.append((path + "/" + k, "missing", rv, None)); continue
            _cmp(path + "/" + k, rv, got[k], problems, tol, skip)
    elif isinstance(ref, list):
        if not isinstance(got, list) or len(got) != len(ref):
            problems.append((path, "list", ref, got)); return
        for i, (rv, gv) in enumerate(zip(ref, got)):
            _cmp("%s[%d]" % (path, i), rv, gv, problems, tol, skip)
    elif isinstance(ref, bool):
        if bool(got) != ref:
            problems.append((path, "bool", ref, got))
    elif isinstance(ref, (int, float)):
        if isinstance(got, bool) or not isinstance(got, (int, float)):
            problems.append((path, "num", ref, got)); return
        if not (abs(float(got) - float(ref)) <= tol or
                (ref != 0 and abs((float(got) - float(ref)) / float(ref)) <= tol)):
            problems.append((path, "num", ref, got))
    else:
        if ref != got:
            problems.append((path, "val", ref, got))


def verify(references=None, tol=TOL, path=core.MEAS_PATH):
    """Pin the three documents against the committed references.

    Returns ``(problems, per_contrast)``; ``problems`` is empty iff every numeric
    value matches to ``tol`` and every structure is identical. The degenerate
    channel is pinned structurally (``did`` / ``n`` / ``cells`` / deltas / r2),
    never on ``did_t`` / ``did_p``.
    """
    references = references or REF_DID
    _keys, rows = load_rows(path)
    problems = []
    per_contrast = {}
    for target_role, reference_role, ckey in DID_CONTRASTS:
        refpath = references[ckey]
        entry = {"target_role": target_role, "reference_role": reference_role,
                 "reference": refpath, "problems": 0}
        if not os.path.exists(refpath):
            problems.append((ckey, "missing-reference", refpath, None))
            entry["problems"] = 1
            per_contrast[ckey] = entry
            continue
        ref = core.load_reference(refpath)
        got = document(rows, target_role, reference_role, len(rows))
        probs = []
        for k in ("asset", "event", "roles_present", "target_role",
                  "reference_role", "method", "rows", "coverage"):
            _cmp("%s/%s" % (ckey, k), ref[k], got[k], probs, tol)
        for orbit in DID_ORBITS:
            for chan in DID_CHANNELS:
                skip = ("did_t", "did_p") if chan in DEGENERATE else ()
                _cmp("%s/series/%s/%s" % (ckey, orbit, chan),
                     ref["series"][orbit][chan],
                     got["series"][orbit][chan], probs, tol, skip)
        entry["problems"] = len(probs)
        entry["channels"] = len(DID_ORBITS) * len(DID_CHANNELS)
        per_contrast[ckey] = entry
        problems.extend(probs)
    return problems, per_contrast



def _fl(s):
    """A table cell as a float, or ``None`` for ``n/a`` / ``—`` / a rho pair."""
    try:
        return float(s)
    except (TypeError, ValueError):
        return None


def _cell_eq(u, v, tol=TOL):
    if u == v:
        return True
    if "/" in u or "/" in v:
        us, vs = u.split("/"), v.split("/")
        if len(us) != len(vs):
            return False
        return all(_cell_eq(a, b, tol) for a, b in zip(us, vs))
    fu, fv = _fl(u), _fl(v)
    if fu is None or fv is None:
        return False
    return abs(fu - fv) <= tol or (fv != 0 and abs((fu - fv) / fv) <= tol)


def _md_equal(a, b, tol=TOL):
    """Compare two reports: prose exactly, a channel row cell-by-cell with a
    numeric tolerance (the degenerate channel's DiD is floating-point noise, so
    its *rendered* number cannot be matched by string equality)."""
    la, lb = a.rstrip("\n").split("\n"), b.rstrip("\n").split("\n")
    if len(la) != len(lb):
        return False
    for x, y in zip(la, lb):
        if not x.startswith("| `"):
            if x != y:
                return False
            continue
        cx = [c.strip().replace("**", "") for c in x.split("|")]
        cy = [c.strip().replace("**", "") for c in y.split("|")]
        if len(cx) != len(cy):
            return False
        chan = cx[1].strip("`") if len(cx) > 1 else ""
        # The t / p of the degenerate channel are noise -> not compared (as in
        # the series pin).
        skip = {6, 7} if chan in DEGENERATE else set()
        for i, (u, v) in enumerate(zip(cx, cy)):
            if i in skip:
                continue
            if not _cell_eq(u, v, tol):
                return False
    return True


def report_matches(references=None, tol=TOL, path=core.MEAS_PATH):
    """Pin the *markdown* reports too -> ``(problems, checked)``."""
    references = references or REF_DID
    _keys, rows = load_rows(path)
    problems = []
    checked = 0
    for target_role, reference_role, ckey in DID_CONTRASTS:
        md = references[ckey][:-5] + ".md"
        if not os.path.exists(md):
            continue
        checked += 1
        got = report(document(rows, target_role, reference_role, len(rows)))
        with open(md) as fh:
            want = fh.read()
        if not _md_equal(got, want, tol):
            problems.append((md, "report-mismatch", None, None))
    return problems, checked


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(description="CTS OCV-Paper DiD layer")
    ap.add_argument("--verify", action="store_true",
                    help="pin the three references (series + reports)")
    ap.add_argument("--run", metavar="KEY",
                    choices=[k for _t, _r, k in DID_CONTRASTS],
                    help="build one contrast's document")
    ap.add_argument("--report", action="store_true",
                    help="with --run, print the markdown report instead of JSON")
    ap.add_argument("--json", action="store_true",
                    help="machine-readable output")
    ap.add_argument("--path", default=core.MEAS_PATH,
                    help="the measurement table (default: the committed extract)")
    args = ap.parse_args(argv)

    if args.run:
        key = args.run
        target, reference = next((t, r) for t, r, k in DID_CONTRASTS if k == key)
        doc = run(target, reference, path=args.path)
        print(report(doc) if args.report else json.dumps(doc, indent=2))
        return 0

    if args.verify:
        problems, per = verify(path=args.path)
        rproblems, rchecked = report_matches(path=args.path)
        for ckey, entry in per.items():
            status = "OK" if entry["problems"] == 0 else "FAIL"
            print("%-24s %-4s  %d channels  (%s vs %s)"
                  % (ckey, status, entry["channels"], entry["target_role"],
                     entry["reference_role"]))
        if rchecked:
            print("markdown reports: %s (%d checked)"
                  % ("OK" if not rproblems else "FAIL", rchecked))
        if args.json:
            print(json.dumps(
                {"problems": [list(map(str, p)) for p in problems],
                 "report_problems": [list(map(str, p)) for p in rproblems],
                 "per_contrast": per}, indent=2))
        else:
            for p in problems:
                print("  MISMATCH", p[0], p[1], "ref=%r got=%r" % (p[2], p[3]))
            for p in rproblems:
                print("  REPORT", p[0], p[1])
        ok = not problems and not rproblems
        print("PIN %s" % ("PASS" if ok else "FAIL"))
        return 0 if ok else 1

    ap.print_help()
    return 0


if __name__ == "__main__":
    sys.exit(main())

