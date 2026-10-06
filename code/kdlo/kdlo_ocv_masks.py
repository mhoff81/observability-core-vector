#!/usr/bin/env python3
"""KDLO OCV-Paper — echo-mask geometry from the committed strip cache.

Port of ``lumo_ocv_masks.py`` (itself a verbatim port of the LUMO project's
``analyze_lumo_tower_coherence.read_complex_bin`` and
``analyze_lumo_echo_mask_vector.connected_components_8`` / ``vector_from_strip``
/ ``add_persistence``), with **one deliberate difference**: the background
median of the mask rule is the *upper* median ``sorted[len // 2]``, because that
is what the tower pipeline uses.

Why the upper median
--------------------
``onboarder.insar_measurements.coherence_masked_pixels`` / ``coherence_gamma2``
are written by the Rust function ``sarsolve::coherence::coherence_mast_echo``:

    intensity  = |z|^2                       (strictly the strip's 11 x 400 samples)
    peak       = max(intensity)
    median     = sorted(intensity)[len // 2]               <- upper median
    mask       = intensity >= 0.30 * peak  and  intensity >= 5.0 * median
    A          = count(mask)                 (error if A < 2)
    gamma2     = |sum(z[mask])|^2 / (A * sum(|z[mask]|^2)),  clamped to [0, 1]

with ``COHERENCE_PEAK_FRAC_DEFAULT = 0.30``, ``COHERENCE_MEDIAN_MULT_DEFAULT =
5.0``, ``COHERENCE_MIN_MASKED_DEFAULT = 2``. Using ``np.median`` (the LUMO
convention, mean of the two central values for even counts) would shift the
threshold and change ``A`` on this record — so the reported value uses the
upper median and ``np.median`` is kept only as the diagnostic
``A_np_median`` / ``gamma2_np_median`` pair.

Input
-----
The committed strips ``data/kdlo/strips/strip_<date>.bin``, written by the Rust
capture binary ``tower/backend/src/cli/analysis/kdlo_strip_capture.rs``, which
decodes exactly what the pipeline decoded:
``decode_slc_rect_at_location(annotation, measurement.tiff, lat, lon, 11, 400)``
— the full-dwell tower strip (``PHASE_STRIP_LINES`` x ``TOWER_PHASE_STRIP_WIDTH``)
of the burst recorded in ``insar_measurements.burst_id``.
Format: ``u32 width, u32 height`` then ``2*width*height`` little-endian f64
(re, im), row-major — the same format LUMO's ``read_complex_bin`` reads.

The masks themselves are *not* committed: they are a deterministic function of
the strips, so ``kdlo_ocv_channels_csv.py`` recomputes them on every build and
verifies ``A == coherence_masked_pixels`` against the committed DB extraction.
``python3 code/kdlo/kdlo_ocv_masks.py --check`` runs that acceptance test
standalone against ``strips/manifest.json`` (no database needed).

Derivation (identical to LUMO, see ``lumo_ocv_masks.py``):
  A = n_masked
  D = n_masked / bbox_area
  F = number of 8-neighbourhood components
  S = ||centroid - peak|| in px
  P = mean pixel frequency across the dates of the majority chip shape
"""
from __future__ import annotations

import argparse
import json
import os
import struct
import sys
from collections import Counter

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, os.pardir, os.pardir, "data", "kdlo")
STRIPS_DIR = os.path.join(DATA, "strips")
MANIFEST_PATH = os.path.join(STRIPS_DIR, "manifest.json")

PEAK_FRAC = 0.30
MEDIAN_MULT = 5.0
MIN_N_MASKED = 2
TOL_GAMMA2 = 1e-12


def read_complex_bin(path):
    """Format: (u32 w, u32 h) + 2*w*h interleaved f64 (real, imag).

    Returns ``(z, width, height)`` with ``z.shape == (height, width)``.
    """
    with open(path, "rb") as fh:
        width, height = struct.unpack("<II", fh.read(8))
        n = width * height
        raw = fh.read(16 * n)
        if len(raw) < 16 * n:
            raise ValueError(f"{path}: truncated complex payload")
        vals = np.frombuffer(raw, dtype="<f8").reshape(n, 2)
        z = (vals[:, 0] + 1j * vals[:, 1]).reshape(height, width)
    return z, width, height


def connected_components_8(mask):
    """Number of connected components (8-neighbourhood, iterative BFS)."""
    h, w = mask.shape
    visited = np.zeros_like(mask, dtype=bool)
    n_components = 0
    for r0 in range(h):
        for c0 in range(w):
            if not mask[r0, c0] or visited[r0, c0]:
                continue
            n_components += 1
            stack = [(r0, c0)]
            visited[r0, c0] = True
            while stack:
                r, c = stack.pop()
                for dr in (-1, 0, 1):
                    for dc in (-1, 0, 1):
                        if dr == 0 and dc == 0:
                            continue
                        nr, nc = r + dr, c + dc
                        if (0 <= nr < h and 0 <= nc < w
                                and mask[nr, nc] and not visited[nr, nc]):
                            visited[nr, nc] = True
                            stack.append((nr, nc))
    return n_components

def vector_from_strip(z, peak_frac=PEAK_FRAC, median_mult=MEDIAN_MULT,
                      min_n_masked=MIN_N_MASKED):
    """Mask + A/D/F/S + gamma2 from one strip (P needs the frequency map, see
    ``add_persistence``). Returns ``None`` if the strip carries no echo mask.

    ``gamma2`` is the pipeline's masked-coherence estimate and must equal
    ``insar_measurements.coherence_gamma2``; the ``*_np_median`` twins are
    diagnostics (LUMO's mask rule) and are never reported as the record's value.
    """
    intensity = z.real ** 2 + z.imag ** 2
    flat = np.asarray(intensity).ravel()          # row-major, as the Rust loop
    peak = float(flat.max())
    if peak <= 0.0:
        return None
    med = float(np.sort(flat)[flat.size // 2])    # upper median (pipeline rule)
    med_np = float(np.median(flat))               # LUMO rule, diagnostic only

    mask = (intensity >= peak_frac * peak) & (intensity >= median_mult * med)
    n_masked = int(mask.sum())
    if n_masked < min_n_masked:
        return None

    rows, cols = np.nonzero(mask)
    row_span = int(rows.max() - rows.min() + 1)
    col_span = int(cols.max() - cols.min() + 1)
    bbox_area = row_span * col_span

    peak_row, peak_col = np.unravel_index(np.argmax(intensity), intensity.shape)
    centroid_row, centroid_col = float(rows.mean()), float(cols.mean())

    sel = z[mask]
    den = n_masked * float(np.sum(sel.real ** 2 + sel.imag ** 2))
    gamma2 = None if den <= 0.0 else float(min(max(abs(sel.sum()) ** 2 / den, 0.0), 1.0))

    mask_np = (intensity >= peak_frac * peak) & (intensity >= median_mult * med_np)
    sel_np = z[mask_np]
    den_np = int(mask_np.sum()) * float(np.sum(sel_np.real ** 2 + sel_np.imag ** 2))
    gamma2_np = (None if den_np <= 0.0
                 else float(min(max(abs(sel_np.sum()) ** 2 / den_np, 0.0), 1.0)))

    return {
        "mask": mask,
        "shape": mask.shape,
        "A": n_masked,
        "D": (n_masked / bbox_area) if bbox_area > 0 else None,
        "F": connected_components_8(mask),
        "S": float(np.hypot(centroid_row - peak_row, centroid_col - peak_col)),
        "bbox_area": bbox_area,
        "row_span": row_span,
        "col_span": col_span,
        "centroid_row": centroid_row,
        "centroid_col": centroid_col,
        "peak_row": int(peak_row),
        "peak_col": int(peak_col),
        "peak": peak,
        "med": med,
        "med_np_median": med_np,
        "gamma2": gamma2,
        "gamma2_np_median": gamma2_np,
        "A_np_median": int(mask_np.sum()),
    }


def add_persistence(rows):
    """P = mean pixel frequency inside the own mask.

    A pixel-exact frequency map is only comparable across dates with the same
    chip shape; dates with a different shape are excluded (P = None) instead of
    being silently misaligned. Every KDLO strip is 400x11, so the majority shape
    is the whole record.
    """
    if not rows:
        return {"majority_shape": None, "n_majority_shape": 0,
                "n_excluded_other_shape": 0}
    shape_counts = Counter(tuple(r["shape"]) for r in rows)
    majority_shape, n_majority = shape_counts.most_common(1)[0]

    freq = np.zeros(majority_shape, dtype=float)
    n_used = 0
    for r in rows:
        if tuple(r["shape"]) == majority_shape:
            freq += r["mask"]
            n_used += 1
    freq /= max(n_used, 1)

    n_excluded = 0
    for r in rows:
        if tuple(r["shape"]) != majority_shape:
            r["P"] = None
            n_excluded += 1
            continue
        m = r["mask"]
        r["P"] = float(freq[m].mean())
        r["P_frac_ge_50pct"] = float((freq[m] >= 0.5).mean())
        r["P_frac_ge_30pct"] = float((freq[m] >= 0.3).mean())
    return {"majority_shape": list(majority_shape), "n_majority_shape": n_majority,
            "n_excluded_other_shape": n_excluded}



def load_manifest(manifest_path=MANIFEST_PATH):
    with open(manifest_path) as fh:
        return json.load(fh)


def load_vectors(strips_dir=STRIPS_DIR, manifest_path=MANIFEST_PATH):
    """All strips that carry an echo mask, in manifest order, with P added.

    Returns ``(rows, meta)``; each row keeps ``date``, ``burst_id``, the
    manifest's ``rust_*`` values and the geometry of ``vector_from_strip``.
    """
    manifest = load_manifest(manifest_path)
    rows, skipped = [], []
    for e in manifest["entries"]:
        path = os.path.join(strips_dir, e["strip_file"])
        z, width, height = read_complex_bin(path)
        v = vector_from_strip(z)
        if v is None:
            skipped.append(e["date"])
            continue
        v["date"] = e["date"]
        v["burst_id"] = e["burst_id"]
        v["subswath"] = e["subswath"]
        v["polarisation"] = e["polarisation"]
        v["cols"], v["rows"] = width, height
        v["strip_sha256"] = e.get("strip_sha256")
        v["rust_masked_pixels"] = e.get("rust_masked_pixels")
        v["rust_gamma2"] = e.get("rust_gamma2")
        rows.append(v)
    meta = add_persistence(rows)
    meta["n_strips"] = len(manifest["entries"])
    meta["n_with_mask"] = len(rows)
    meta["n_without_mask"] = len(skipped)
    meta["dates_without_mask"] = skipped
    return rows, meta


def check(rows, meta, tol_gamma2=TOL_GAMMA2):
    """Acceptance test: recomputed A/gamma2 == the pipeline's recorded values."""
    problems = []
    shape_counts = Counter(tuple(r["shape"]) for r in rows)
    for r in rows:
        tag = f"{r['date']} {r['subswath']}/{r['polarisation']}"
        if r["rust_masked_pixels"] is None:
            problems.append(f"{tag}: manifest has no rust_masked_pixels")
            continue
        if r["A"] != int(r["rust_masked_pixels"]):
            problems.append(f"{tag}: A={r['A']} != coherence_masked_pixels="
                            f"{r['rust_masked_pixels']}")
        if r["rust_gamma2"] is None:
            problems.append(f"{tag}: manifest has no rust_gamma2")
        elif abs(r["gamma2"] - float(r["rust_gamma2"])) > tol_gamma2:
            problems.append(f"{tag}: gamma2={r['gamma2']!r} != recorded "
                            f"{r['rust_gamma2']!r}")
        if tuple(r["shape"]) != (400, 11):
            problems.append(f"{tag}: strip shape {r['shape']} != (400, 11)")
    return {"n_checked": len(rows), "problems": problems,
            "shapes": {f"{a}x{b}": n for (a, b), n in shape_counts.items()},
            "n_P_set": sum(1 for r in rows if r.get("P") is not None),
            "n_P_excluded_other_shape": meta["n_excluded_other_shape"]}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--strips-dir", default=STRIPS_DIR)
    ap.add_argument("--manifest", default=MANIFEST_PATH)
    ap.add_argument("--json", action="store_true",
                    help="print the derived vectors as JSON (masks omitted)")
    args = ap.parse_args()

    rows, meta = load_vectors(args.strips_dir, args.manifest)
    rep = check(rows, meta)
    if args.json:
        print(json.dumps([{k: v for k, v in r.items() if k != "mask"}
                          for r in rows], indent=1, sort_keys=True))
    print(json.dumps({k: v for k, v in rep.items() if k != "problems"},
                     indent=1, sort_keys=True))
    print(json.dumps(meta, indent=1, sort_keys=True))
    if rep["problems"]:
        for p in rep["problems"][:40]:
            print("PROBLEM:", p)
        raise SystemExit(f"VERIFICATION FAILED ({len(rep['problems'])} problems)")
    print("VERIFICATION OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())

