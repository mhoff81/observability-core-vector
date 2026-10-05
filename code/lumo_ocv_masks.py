#!/usr/bin/env python3
"""OCV-Paper — mask-geometry ports (verbatim from the LUMO project).

Port of ``analyze_lumo_tower_coherence.read_complex_bin``,
``analyze_lumo_echo_mask_vector.connected_components_8`` /
``vector_from_strip`` / ``add_persistence``.

This file is needed *only* by the one-shot generator ``lumo_ocv_channels_csv.py``
(masks are not committed because they require the 1394-directory burst cache).
The figure pipeline itself reads only the derived ``data/lumo_ocv_channels.csv``
and needs neither cache nor mask.
"""
from __future__ import annotations

import struct
from collections import Counter

import numpy as np

PEAK_FRAC = 0.30
MEDIAN_MULT = 5.0
MIN_N_MASKED = 2


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
    """Mask + A/D/F/S (P needs the global frequency map, see add_persistence)."""
    intensity = np.abs(z) ** 2
    peak = float(intensity.max())
    if peak <= 0.0:
        return None
    med = float(np.median(intensity))
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
    }


def add_persistence(rows):
    """P = mean pixel frequency inside the own mask.

    A pixel-exact frequency map is only comparable across dates with the same
    chip shape; dates with a different shape are excluded (P = None) instead
    of being silently misaligned.
    """
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
