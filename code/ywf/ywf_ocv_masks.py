#!/usr/bin/env python3
"""YWF OCV — the echo mask of the committed 7x7 window-payload cache.

Yeongdeok Wind Farm (Changpo Wind Power Complex, Samgye-ri, 36.4236 N
129.4203 E) is Unit 21, a 100 m steel monopole wind-turbine tower with a 4.6 m
base diameter that **collapsed onto a public road on 2026-02-02**. The export
stores, per acquisition, one 7x7 complex window per **mast-section request**
(segment 0..4 -> ``turbine-mast-section-0`` .. ``-4``, see
``data/ywf/ywf_segments.txt``).

Unlike the bridge sites the committed payload of this site is *small* (57
windows x 49 complex samples = 38 kB), so ``data/ywf/ywf_windows_full.txt`` is
committed in full and the mask layer can be recomputed from it directly. The
cache ``data/ywf/ywf_windows_mask_cache.txt`` is still written, so the site has
the same committed mask input as Carola/Morandi/CTS and the same re-verification
step; ``ywf_ocv_masks.py --check`` re-derives every dimension from it.

The mask rule is the project-wide tower/bridge echo definition (identical to
``analyze_carola_coherence.coherence_gamma2`` and the backend's
``sarsolve::coherence::coherence_mast_echo``):

    intensity = |z|^2
    peak      = max(intensity)
    median    = np.median(intensity)
    mask      = intensity >= 0.30 * peak  and  intensity >= 5.0 * median
    A         = count(mask)                    (no echo if A < 2)
    gamma2    = |sum(z[mask])|^2 / (A * sum(|z[mask]|^2)),  clamped to [0, 1]

Derivation
----------
  A = n_masked                       (echo pixels, ``coherence_masked_pixels``)
  D = A / bbox_area                  (bounding-box fill of the echo)
  F = number of 8-neighbourhood components
  S = ||mask centroid - peak pixel|| in px
  P = mean pixel frequency inside the own mask (over all chips of the shape)

Segment resolution
------------------
The record is **not** segment-resolved: of 57 committed windows 24 are byte
identical duplicates that differ only in the request's ``segment_index``
(2 vs. 3) on a single date — i.e. both mast-section requests decoded the *same*
chip at the asset point (the pre-fix behaviour described in
``carola_ocv_masks``). The cache therefore de-duplicates on
``(asset, date, md5(payload))``, first occurrence wins, and the manifest records
``segment_resolved: false``. Consequently this package carries **one** mask
layer, the echo mask: a per-mast-section layer would be dishonest.

Usage
-----
  python3 code/ywf/ywf_ocv_masks.py --extract   # recompute cache + manifest
  python3 code/ywf/ywf_ocv_masks.py --check     # acceptance test of the cache
  python3 code/ywf/ywf_ocv_masks.py --json      # print the derived 6D vectors
Command lines are always run from the **repository root**.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
from collections import Counter

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, os.pardir, os.pardir, "data", "ywf")

CACHE_PATH = os.path.join(DATA, "ywf_windows_mask_cache.txt")
MANIFEST_PATH = os.path.join(DATA, "ywf_windows_mask_cache.manifest.json")
MEAS_PATH = os.path.join(DATA, "ywf_measurements_full.txt")
INDEX_PATH = os.path.join(DATA, "ywf_windows_index.txt")
SEGMENTS_PATH = os.path.join(DATA, "ywf_segments.txt")
DEFAULT_WINDOWS = os.path.join(DATA, "ywf_windows_full.txt")

PEAK_FRAC = 0.30
MEDIAN_MULT = 5.0
MIN_N_MASKED = 2
WINDOW_SIZE = 7
RULE_TOL = 1e-9

CACHE_HEADER = (
    "# ywf windows mask cache v1 — the echo mask of every committed 7x7 window\n"
    "# rule: intensity >= 0.3 * max(intensity) AND intensity >= 5.0 * np.median(intensity); min_n_masked=2\n"
    "# dedup: (asset, date, md5(payload)), first occurrence wins\n"
    "# extractor: code/ywf/ywf_ocv_masks.py --extract\n"
)
CACHE_COLUMNS = ["mid", "asset", "segment", "date", "orbit", "w", "h",
                 "peak", "median", "max_masked", "min_masked",
                 "ss_re", "ss_im", "ss_abs2", "n_masked",
                 "peak_row", "peak_col", "n_components", "largest_component", "px"]

# The mask dimensions the cache carries (gamma2 is derived, not stored).
MASK_COLUMNS = ["A", "D", "F", "S", "P", "gamma2"]


# ---------------------------------------------------------------------------
# The echo mask itself (port of carola_ocv_masks.echo_mask / largest_component)
# ---------------------------------------------------------------------------
def largest_component(mask):
    """``(n_components, largest_size)`` of a 2-D boolean mask, 8-connectivity.

    Dependency-free BFS, verbatim from ``carola_ocv_masks.largest_component``.
    """
    h, w = mask.shape
    seen = np.zeros_like(mask, dtype=bool)
    sizes = []
    for r0 in range(h):
        for c0 in range(w):
            if not mask[r0, c0] or seen[r0, c0]:
                continue
            stack = [(r0, c0)]
            seen[r0, c0] = True
            size = 0
            while stack:
                r, c = stack.pop()
                size += 1
                for dr in (-1, 0, 1):
                    for dc in (-1, 0, 1):
                        rr, cc = r + dr, c + dc
                        if (0 <= rr < h and 0 <= cc < w and mask[rr, cc]
                                and not seen[rr, cc]):
                            seen[rr, cc] = True
                            stack.append((rr, cc))
            sizes.append(size)
    return len(sizes), (max(sizes) if sizes else 0)


def echo_mask(z, peak_frac=PEAK_FRAC, median_mult=MEDIAN_MULT,
              min_n_masked=MIN_N_MASKED):
    """The project-wide echo mask of one complex window.

    ``z`` is the ``(h, w)`` complex chip. Returns ``None`` when the window has
    no echo (``peak <= 0`` or fewer than ``min_n_masked`` masked pixels), else a
    dict with the mask and the sufficient statistics of the six dimensions.
    """
    intensity = np.abs(z) ** 2
    peak = float(intensity.max())
    if peak <= 0.0:
        return None
    median = float(np.median(intensity))
    mask = (intensity >= peak_frac * peak) & (intensity >= median_mult * median)
    n_masked = int(mask.sum())
    if n_masked < min_n_masked:
        return None
    rows, cols = np.nonzero(mask)
    sel = z[mask]
    masked_int = intensity[mask]
    peak_row, peak_col = np.unravel_index(int(np.argmax(intensity)), intensity.shape)
    n_comp, largest = largest_component(mask)
    return {
        "mask": mask,
        "shape": (int(z.shape[0]), int(z.shape[1])),
        "A": n_masked,
        "n_masked": n_masked,
        "peak": peak,
        "median": median,
        "max_masked": float(masked_int.max()),
        "min_masked": float(masked_int.min()),
        "ss_re": float(sel.real.sum()),
        "ss_im": float(sel.imag.sum()),
        "ss_abs2": float((sel.real ** 2 + sel.imag ** 2).sum()),
        "peak_row": int(peak_row),
        "peak_col": int(peak_col),
        "px": sorted((int(r), int(c)) for r, c in zip(rows, cols)),
        "n_components": n_comp,
        "largest_component": largest,
    }


def gamma2_of(ss_re, ss_im, ss_abs2, n_masked):
    """``|sum z|^2 / (A * sum |z|^2)`` from the committed sums (clamped)."""
    den = float(n_masked) * float(ss_abs2)
    if den <= 0.0:
        return None
    return float(min(max((ss_re * ss_re + ss_im * ss_im) / den, 0.0), 1.0))


def vector_from_row(row):
    """The 6D echo-mask vector of one cache row (mask rule re-verified)."""
    px = row["px"]
    A = len(px)
    if A != row["n_masked"]:
        raise ValueError(f"{row['mid']}: n_masked {row['n_masked']} != |px| {A}")
    if A < MIN_N_MASKED:
        raise ValueError(f"{row['mid']}: A={A} < {MIN_N_MASKED}")
    peak, median = row["peak"], row["median"]
    lo_thr = max(PEAK_FRAC * peak, MEDIAN_MULT * median)
    if row["min_masked"] * (1.0 + RULE_TOL) < lo_thr:
        raise ValueError(f"{row['mid']}: committed mask violates the rule "
                         f"(min|z|^2={row['min_masked']!r} < {lo_thr!r})")
    if abs(row["max_masked"] - peak) > RULE_TOL * max(peak, 1e-30):
        raise ValueError(f"{row['mid']}: peak {peak!r} is not a masked pixel "
                         f"(max_masked={row['max_masked']!r})")
    rows = np.array([r for r, _ in px], dtype=float)
    cols = np.array([c for _, c in px], dtype=float)
    row_span = int(rows.max() - rows.min() + 1)
    col_span = int(cols.max() - cols.min() + 1)
    bbox_area = int(row_span * col_span)
    cen_r, cen_c = float(rows.mean()), float(cols.mean())
    return {
        "A": A,
        "D": (A / bbox_area) if bbox_area > 0 else None,
        "F": row["n_components"],
        "S": float(np.hypot(cen_r - row["peak_row"], cen_c - row["peak_col"])),
        "gamma2": gamma2_of(row["ss_re"], row["ss_im"], row["ss_abs2"], A),
        "bbox_area": bbox_area,
        "row_span": row_span,
        "col_span": col_span,
        "centroid_row": cen_r,
        "centroid_col": cen_c,
        "peak_row": row["peak_row"],
        "peak_col": row["peak_col"],
        "n_components": row["n_components"],
        "largest_component": row["largest_component"],
        "fragmentation": (1.0 - row["largest_component"] / float(max(A, 1))),
    }


# ---------------------------------------------------------------------------
# P — the mask persistence over all chips of the shape
# ---------------------------------------------------------------------------
def pixel_frequencies(rows):
    """``{(r, c): frequency}`` — share of chips whose echo mask holds the pixel.

    The YWF analogue of the LUMO/Carola persistence reference: every chip here
    has the same 7x7 shape, so the shape filter is trivially satisfied and the
    chip count is the deduplicated one.
    """
    n_chips = len(rows)
    counts = Counter()
    for r in rows:
        for px in r["px"]:
            counts[px] += 1
    if not n_chips:
        return {}
    return {px: n / n_chips for px, n in counts.items()}


def add_persistence(rows):
    """Append ``P`` (and its ingredient) to every cache row, in place."""
    freq = pixel_frequencies(rows)
    for r in rows:
        vals = [freq[px] for px in r["px"]]
        r["P"] = float(np.mean(vals)) if vals else None
        r["P_reference_chips"] = len(rows)
    return rows


# ---------------------------------------------------------------------------
# I/O of the committed cache
# ---------------------------------------------------------------------------
def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _row_to_line(r):
    px = ";".join(f"{r},{c}" for r, c in r["px"])
    cols = {
        "mid": r["mid"], "asset": r["asset"], "segment": r["segment"],
        "date": r["date"], "orbit": r["orbit"], "w": r["shape"][1],
        "h": r["shape"][0], "peak": repr(r["peak"]), "median": repr(r["median"]),
        "max_masked": repr(r["max_masked"]), "min_masked": repr(r["min_masked"]),
        "ss_re": repr(r["ss_re"]), "ss_im": repr(r["ss_im"]),
        "ss_abs2": repr(r["ss_abs2"]), "n_masked": r["n_masked"],
        "peak_row": r["peak_row"], "peak_col": r["peak_col"],
        "n_components": r["n_components"],
        "largest_component": r["largest_component"], "px": px,
    }
    return "|".join(str(cols[c]) for c in CACHE_COLUMNS)


def load_cache(path=CACHE_PATH):
    """The committed mask cache as a list of row dicts (numbers parsed)."""
    rows = []
    if not os.path.isfile(path):
        return rows
    with open(path) as fh:
        lines = fh.read().splitlines()
    body = [ln for ln in lines if ln and not ln.startswith("#")]
    if not body:
        return rows
    keys = body[0].split("|")
    for line in body[1:]:
        p = line.split("|")
        if len(p) != len(keys):
            continue
        d = dict(zip(keys, p))
        for k in ("w", "h", "n_masked", "peak_row", "peak_col",
                  "n_components", "largest_component"):
            d[k] = int(d[k])
        for k in ("peak", "median", "max_masked", "min_masked",
                  "ss_re", "ss_im", "ss_abs2"):
            d[k] = float(d[k])
        d["shape"] = (d["h"], d["w"])
        d["px"] = [tuple(int(v) for v in ch.split(","))
                   for ch in d["px"].split(";") if ch]
        rows.append(d)
    add_persistence(rows)
    return rows


def load_manifest(path=MANIFEST_PATH):
    if not os.path.isfile(path):
        return {}
    with open(path) as fh:
        return json.load(fh)


# ---------------------------------------------------------------------------
# Extraction: windows payload -> cache + manifest
# ---------------------------------------------------------------------------
def read_index(path=INDEX_PATH):
    """``measurement_id -> index row`` of the committed window index.

    The committed index has **no header line**: every line is

        measurement_id|asset_id|asset_name|segment_index|acquisition_ts|
        request_id|window_width|window_height
    """
    out = {}
    if not os.path.isfile(path):
        return out
    for line in open(path).read().splitlines():
        p = line.split("|")
        if len(p) < 8:
            continue
        out[p[0]] = {"asset_id": p[1], "asset": p[2], "segment": p[3],
                     "ts": p[4], "date": p[4][:10], "request_id": p[5],
                     "w": int(p[6]), "h": int(p[7])}
    return out


def read_measurements(path=MEAS_PATH):
    """``measurement_id -> column dict`` of the committed measurement extract."""
    out = {}
    if not os.path.isfile(path):
        return out
    lines = open(path).read().splitlines()
    if not lines:
        return out
    keys = lines[0].split("|")
    for line in lines[1:]:
        p = line.split("|")
        if len(p) != len(keys):
            continue
        d = dict(zip(keys, p))
        if d.get("id"):
            out[d["id"]] = d
    return out


def extract(windows_path=DEFAULT_WINDOWS, out_path=CACHE_PATH,
            manifest_path=MANIFEST_PATH):
    """Recompute the mask cache from the committed 7x7 window payloads."""
    index = read_index()
    meas = read_measurements()
    rows = []
    n_dup = 0
    n_no_index = 0
    n_no_echo = 0
    sizes = Counter()
    seen = {}
    groups = {}
    for line in open(windows_path).read().splitlines():
        if not line.strip():
            continue
        mid, w, h, payload = line.split("|", 3)
        meta = index.get(mid)
        if meta is None:
            n_no_index += 1
            continue
        digest = hashlib.md5(payload.encode()).hexdigest()[:12]
        # every committed line enters the segment diagnostic, duplicates too
        groups.setdefault((meta["asset"], meta["date"]), []).append(
            (meta["segment"], digest))
        key = (meta["asset"], meta["date"], digest)
        if key in seen:
            n_dup += 1
            continue
        samples = np.array([complex(p[0], p[1])
                            for p in json.loads(payload)], dtype=complex)
        hh, ww = int(h), int(w)
        if samples.size != hh * ww:
            raise SystemExit(f"{mid}: payload {samples.size} != {hh}x{ww}")
        z = samples.reshape(hh, ww)
        sizes[(hh, ww)] += 1
        res = echo_mask(z)
        if res is None:
            n_no_echo += 1
        else:
            m = meas.get(mid, {})
            res.update({"mid": mid, "asset": meta["asset"],
                        "segment": meta["segment"], "date": meta["date"],
                        "orbit": m.get("orbit_direction"),
                        "pass_label": m.get("pass_label"),
                        "request_id": meta["request_id"]})
            rows.append(res)
        seen[key] = mid
    add_persistence(rows)

    size = max(sizes, key=sizes.get)[0] if sizes else WINDOW_SIZE
    multi = 0
    multi_differing = 0
    repeated_segment = 0
    for v in groups.values():
        segs = {s for s, _ in v}
        digs = {p for _, p in v}
        if len(segs) > 1:
            multi += 1
            if len(digs) > 1:
                multi_differing += 1
        if len(v) > len(segs):
            repeated_segment += 1
    resolved = multi_differing
    manifest = {
        "generator": "code/ywf/ywf_ocv_masks.py --extract",
        "source_windows": os.path.basename(windows_path),
        "source_windows_sha256": _sha256(windows_path),
        "source_index": os.path.basename(INDEX_PATH),
        "source_segments": os.path.basename(SEGMENTS_PATH),
        "mask_rule": {"peak_frac": PEAK_FRAC, "median_mult": MEDIAN_MULT,
                      "median_kind": "np.median", "min_n_masked": MIN_N_MASKED},
        "window_size": size,
        "n_rows": len(rows),
        "n_deduplicated_same_payload": n_dup,
        "n_without_index": n_no_index,
        "n_without_echo": n_no_echo,
        "window_sizes": {f"{a}x{b}": n for (a, b), n in sizes.items()},
        "asset_date_groups": len(groups),
        "asset_date_groups_with_several_segments": multi,
        "asset_date_groups_with_several_segments_with_differing_payloads":
            multi_differing,
        "asset_date_groups_with_differing_payloads": resolved,
        "asset_date_groups_with_repeated_segment": repeated_segment,
        "segment_resolved": bool(multi and multi_differing == multi),
    }
    with open(out_path, "w") as fh:
        fh.write(CACHE_HEADER)
        fh.write("|".join(CACHE_COLUMNS) + "\n")
        for r in rows:
            fh.write(_row_to_line(r) + "\n")
    with open(manifest_path, "w") as fh:
        json.dump(manifest, fh, indent=1, sort_keys=True)
    print(f"written: {os.path.relpath(out_path, DATA)} ({len(rows)} rows)")
    print(json.dumps(manifest, indent=1, sort_keys=True))
    return rows, manifest


# ---------------------------------------------------------------------------
# Acceptance test
# ---------------------------------------------------------------------------
def check(rows, manifest=None):
    """Re-verify the mask rule and the derived geometry of every cache row."""
    problems = []
    for r in rows:
        try:
            vector_from_row(r)
        except ValueError as exc:
            problems.append(str(exc))
        if tuple(r["shape"]) != (WINDOW_SIZE, WINDOW_SIZE):
            problems.append(f"{r['mid']}: shape {r['shape']} != "
                            f"({WINDOW_SIZE}, {WINDOW_SIZE})")
        if len(set(r["px"])) != len(r["px"]):
            problems.append(f"{r['mid']}: duplicate mask pixels")
    shapes = Counter(f"{h}x{w}" for (h, w) in (r["shape"] for r in rows))
    return {"n_checked": len(rows), "shapes": dict(shapes),
            "problems": problems[:40], "n_problems": len(problems)}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--extract", action="store_true",
                    help="recompute the cache from the committed window payloads")
    ap.add_argument("--windows", default=DEFAULT_WINDOWS)
    ap.add_argument("--out", default=CACHE_PATH)
    ap.add_argument("--manifest", default=MANIFEST_PATH)
    ap.add_argument("--check", action="store_true",
                    help="acceptance test of the committed cache")
    ap.add_argument("--json", action="store_true",
                    help="print the derived 6D vectors")
    args = ap.parse_args()

    if args.extract:
        extract(args.windows, args.out, args.manifest)

    rows = load_cache(args.out)
    rep = check(rows, load_manifest(args.manifest))
    if args.json:
        out = []
        for r in rows:
            v = vector_from_row(r)
            v["mid"] = r["mid"]
            v["segment"] = r["segment"]
            v["date"] = r["date"]
            v["P"] = r.get("P")
            out.append(v)
        print(json.dumps(out, indent=1, sort_keys=True))
    print(json.dumps({k: v for k, v in rep.items() if k != "problems"},
                     indent=1, sort_keys=True))
    if rep["problems"]:
        for p in rep["problems"][:40]:
            print("PROBLEM:", p)
        raise SystemExit(f"VERIFICATION FAILED ({rep['n_problems']} problems)")
    print("VERIFICATION OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
