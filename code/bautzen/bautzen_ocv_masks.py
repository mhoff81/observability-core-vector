#!/usr/bin/env python3
"""OCV-Paper (Bautzen bridge site) — mask layer (the *swapped* layer).

This is the one layer that differs from the LUMO package; everything else of
this site mirrors it. The *echo mask* of the LUMO site
(``|z|^2 >= 0.30 * peak`` and ``>= 5.0 * median``, at least 2 px) is replaced by
the **deck-edge mask** of the openBridgeLAB bridge at Bautzen.

Why the echo mask cannot be used here (recomputed, see the mask reference JSON):
the rect chips are clutter-dominated, so the global ``0.30*peak and 5*median``
rule fires far away from the bridge and yields ``n ~ 0`` at the deck; no
echo-mask vector of the deck exists. The deck line is only reachable through the
*local, geometrically anchored* band:

  anchor     row 40 (rect centre = deck centre, deliberately NOT searched: a
             data-driven anchor self-selects, because the brightest row band
             would be declared "the" band and the control area would answer too)
  length     columns 36..44 (deck axis 28.90 m / 3.37 m per px -> 8.6 px)
  band       anchor row +/-1 (the two long deck edges are < 1 px apart and not
             separable at S1-SLC resolution, so the mask captures the *fused*
             deck dihedral line, not both edges)
  threshold  |z|^2 >= 1.5 * median(neighbouring azimuth rows) (rest of the ROI)
  gate       n_masked >= 2 and band contrast >= 1.5 => "deck edge present"

Ports (verbatim, all thresholds included) from
``analyze_bautzen_deck_edge_mask.py``: ``load_intensity``, ``load_complex``,
``largest_component``, ``date_mask``, ``analyse_series``; and from
``analyze_bautzen_deck_edge_vs_state.py``: ``gamma2``, ``house_mask``.

One documented deviation from the LUMO vector definition (and only one):
  the peak pixel behind ``S`` and ``peak_r``/``peak_c`` is the **deck-local peak
  inside the band**, not the global chip maximum. The global maximum is clutter
  (see above); an ``S`` measured against it would describe the distance to the
  chip's brightest clutter cell, not the geometry of the deck line. The
  deck-local peak is stored as ``band_peak_r`` / ``band_peak_c`` and is pinned
  explicitly against the reference.
"""
from __future__ import annotations

import json
import math
from collections import Counter

import numpy as np

import bautzen_ocv_stats as st

# ── Geometry (verified at the geolocation grid) ─────────────────────────────
RANGE_M_PER_PX = 3.37      # column spacing (ground range)
AZIMUTH_M_PER_PX = 12.38   # row spacing (azimuth)
DECK_SPAN_M = 28.90        # deck axis (2 x 14.45 m)
DECK_WIDTH_M = 5.90        # deck width (across)
CENTER_ROW = 40            # deck centre = rect centre (geodetic anchor)
CENTER_COL = 40
ANCHOR_HALF = 4            # ROI half width (azimuth rows) around the deck band
COL_HALF = 4               # deck length +/-4 columns (8.6 px -> 9 columns)
BAND_HALF = 1              # deck band = anchor row +/-1
K_LOCAL = 1.5              # local threshold relative to the neighbour azimuth
MIN_N_MASKED = 2           # house gate (as in the LUMO site)
CONTRAST_DETECT = 1.5      # band/neighbour from here on = deck edge visible

# ── Reference annulus (the v2 rule) ─────────────────────────────────────────
# The committed rule takes its local threshold reference from the *neighbouring
# azimuth rows* inside a 9 x 9 ROI. That reference sits right beside the deck:
# measured on the committed chips it is 23 % (ASC) / 7 % (DESC) *brighter* than
# the surrounding land, so it raises the threshold and hides deck edges. The ring
# below takes the reference from an annulus around the support instead, and is
# therefore rule-independent in the same sense as the OCV's background rim.
#
# The chip must carry the whole annulus. The committed 80 x 80 chip cannot (a
# 200 m ring needs 119 x 33 at 3.70 m/bin and 129 x 33 at 3.38 m/bin), which is
# why `reference_mask` reports `complete=False` for it instead of clipping.
RING_GUARD_M = 30.0      # the annulus starts here, outside the deck edge
RING_M = 200.0           # and ends here
REFERENCE_MODES = ("neighbour", "ring")

RING_CONFIG = {
    "reference": "ring",
    "guard_m": RING_GUARD_M,
    "ring_m": RING_M,
    "band_lines": 2 * BAND_HALF + 1,
    "band_bins": 2 * COL_HALF + 1,
}


def anchor_of(shape):
    """``(row, col)`` of the deck anchor inside a chip: ``rows // 2``, ``cols // 2``.

    This is the committed convention *and* the new one. The published 80 x 80
    chips were cut as ``[line_c - 40, line_c + 40)``, so the geodetic anchor sits
    at index 40 of 0..79 — which is ``80 // 2``. An odd chip of 33 x 119 puts it
    at 16 x 59 — which is also ``n // 2``. One rule, two chip shapes.
    """
    return int(shape[0]) // 2, int(shape[1]) // 2


def _block(shape, row0, row1, col0, col1):
    """A boolean block, clipped to the chip (clipping is *reported*, not hidden)."""
    mask = np.zeros(shape, dtype=bool)
    mask[max(row0, 0):max(row1 + 1, 0), max(col0, 0):max(col1 + 1, 0)] = True
    return mask


def reference_mask(shape, spacing, guard_m=RING_GUARD_M, ring_m=RING_M):
    """The support band, its guard and the reference annulus of one chip.

    ``spacing`` is ``(ground range, azimuth)`` metres per sample — the channel's
    own, so the annulus is metric in both axes. Returns the three masks, their
    cell counts, the chip size the annulus *needs*, and ``complete``: a chip too
    small for the annulus is refused rather than measured against a truncated
    reference.
    """
    rows, cols = int(shape[0]), int(shape[1])
    line_c, col_c = anchor_of(shape)
    half_l, half_b = BAND_HALF, COL_HALF
    across_m, along_m = float(spacing[0]), float(spacing[1])
    guard_b = int(math.ceil(guard_m / across_m))
    guard_l = int(math.ceil(guard_m / along_m))
    ring_b = int(math.ceil(ring_m / across_m))
    ring_l = int(math.ceil(ring_m / along_m))

    support = _block(shape, line_c - half_l, line_c + half_l,
                     col_c - half_b, col_c + half_b)
    inner = _block(shape, line_c - half_l - guard_l, line_c + half_l + guard_l,
                   col_c - half_b - guard_b, col_c + half_b + guard_b)
    outer = _block(shape, line_c - half_l - ring_l, line_c + half_l + ring_l,
                   col_c - half_b - ring_b, col_c + half_b + ring_b)
    needed = (2 * (half_l + ring_l) + 1, 2 * (half_b + ring_b) + 1)
    return {
        "support": support,
        "guard": inner & ~support,
        "ring": outer & ~inner,
        "counts": {"support": int(support.sum()),
                   "guard": int((inner & ~support).sum()),
                   "ring": int((outer & ~inner).sum()),
                   "chip": int(rows * cols)},
        "needed": {"lines": needed[0], "bins": needed[1]},
        "shape": (rows, cols),
        "spacing": (across_m, along_m),
        "guard_m": float(guard_m), "ring_m": float(ring_m),
        # A chip smaller than `needed` clips the annulus at its border: the
        # reference would then be a strip, not an annulus. Refuse it.
        "complete": bool(needed[0] <= rows and needed[1] <= cols),
    }

GAMMA2_PEAK_FRAC = 0.30
GAMMA2_MEDIAN_MULT = 5.0
GAMMA2_MIN_MASKED = 2
LOW_N_GATE = 5             # floor correction is trustworthy from here on

# The config block of the committed mask reference, key for key.
CONFIG = {
    "range_m_per_px": RANGE_M_PER_PX,
    "azimuth_m_per_px": AZIMUTH_M_PER_PX,
    "center_row": CENTER_ROW,
    "center_col": CENTER_COL,
    "anchor_half": ANCHOR_HALF,
    "col_half": COL_HALF,
    "band_half": BAND_HALF,
    "k_local": K_LOCAL,
    "min_n_masked": MIN_N_MASKED,
    "contrast_detect": CONTRAST_DETECT,
}


def load_intensity(path):
    """Rect file (``date|orbit|rows|cols|[[re,im],...]``) -> list (date, |z|^2)."""
    recs = []
    with open(path) as fh:
        for line in fh:
            parts = line.strip().split("|")
            if len(parts) < 5:
                continue
            rows, cols = int(parts[2]), int(parts[3])
            flat = np.asarray(json.loads(parts[4]), dtype=float)
            z = (flat[:, 0] + 1j * flat[:, 1]).reshape(rows, cols)
            recs.append((parts[0], np.abs(z) ** 2))
    return recs


def load_complex(path):
    """Rect file -> list ``(date, z)`` with the **complex** chip per date.

    Parsing deliberately identical to ``load_intensity`` (same JSON floats in the
    same order), so ``date_mask(np.abs(z) ** 2)`` is bit-identical to the mask of
    the inventory script: no second threshold set, no method drift.
    """
    recs = []
    with open(path) as fh:
        for line in fh:
            parts = line.strip().split("|")
            if len(parts) < 5:
                continue
            rows, cols = int(parts[2]), int(parts[3])
            flat = np.asarray(json.loads(parts[4]), dtype=float)
            z = (flat[:, 0] + 1j * flat[:, 1]).reshape(rows, cols)
            recs.append((parts[0], z))
    return recs


def largest_component(mask):
    """``(n_components, largest_size)`` of a 2-D bool mask, 8-neighbourhood."""
    h, w = mask.shape
    seen = np.zeros_like(mask, dtype=bool)
    comps = []
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
            comps.append(size)
    if not comps:
        return 0, 0
    return len(comps), max(comps)


def gamma2(z_sel):
    """gamma^2 = |sum z|^2 / (N * sum |z|^2) plus the floor correction.

    Returns ``(raw, floor_corrected, N)``; ``(None, None, N)`` for N < 2 or
    vanishing energy. ``floor_corrected = (N*raw - 1)/(N - 1)`` **without
    clipping**: negative values stay, because they mean "below the noise floor
    assumed for this N", not negative coherence.
    """
    n = int(np.size(z_sel))
    if n < 2:
        return None, None, n
    s = float(np.sum(np.abs(z_sel) ** 2))
    if not s > 0.0:
        return None, None, n
    raw = float(abs(np.sum(z_sel)) ** 2 / (n * s))
    return raw, (n * raw - 1.0) / (n - 1.0), n


def house_mask(inten):
    """Signal-based house/peak mask — **not** structurally attributed.

    Rule as in ``sarsolve::coherence``: ``|z|^2 >= 0.30*peak`` and
    ``>= 5.0*median``, at least 2 pixels. The result must **not** be read as a
    deck-edge gamma^2: the mask picks the brightest pixels of the chip (chip
    edge, buildings, car park), not the deck.
    """
    peak = float(inten.max()) if inten.size else 0.0
    if not peak > 0.0:
        return np.zeros(inten.shape, dtype=bool)
    m = ((inten >= GAMMA2_PEAK_FRAC * peak)
         & (inten >= GAMMA2_MEDIAN_MULT * float(np.median(inten))))
    return m if int(m.sum()) >= GAMMA2_MIN_MASKED else np.zeros(inten.shape, dtype=bool)


def date_mask(inten, ref_mask=None):
    """Deck-edge mask of one date: anchor row, band, local threshold, vector.

    The mask itself is a verbatim port (same thresholds, same order of
    operations) of ``analyze_bautzen_deck_edge_mask.date_mask``, generalised in
    exactly two places:

    * the anchor is the chip's own centre (``anchor_of``), which is 40 x 40 for
      the committed 80 x 80 chips and 16 x 59 for the new 33 x 119 ones — the
      committed numbers are therefore reproduced bit for bit;
    * ``ref_mask`` selects the **reference cells**. ``None`` keeps the committed
      neighbour rows (the v1 rule); a ring mask gives the background annulus of
      ``reference_mask``. Everything the rule then computes is identical.

    The returned ``info`` carries the OCV vector of the mask (``A``, ``D``,
    ``F``, ``S``), the deck-local peak behind ``S``, and the provenance of the
    reference (``ref_source``, ``ref_cells``, ``ref_median``).

    Returns ``(info, full)`` with the full-chip bool mask.
    """
    rows, cols = inten.shape
    anchor_row, anchor_col = anchor_of(inten.shape)
    c_lo, c_hi = anchor_col - COL_HALF, anchor_col + COL_HALF
    anchor = anchor_row
    r_lo0, r_hi0 = anchor - ANCHOR_HALF, anchor + ANCHOR_HALF
    roi = inten[r_lo0:r_hi0 + 1, c_lo:c_hi + 1]
    band = np.zeros_like(roi, dtype=bool)
    band[ANCHOR_HALF - BAND_HALF:ANCHOR_HALF + BAND_HALF + 1, :] = True
    # Diagnostic: row of the row-profile maximum (report only, NOT the anchor).
    prof = roi.mean(axis=1)
    peak_row = r_lo0 + int(np.argmax(prof))
    if ref_mask is None:
        ref = float(np.median(roi[~band])) or 1.0
        ref_source = "neighbour_rows"
        ref_cells = int((~band).sum())
    else:
        ref_mask = np.asarray(ref_mask, dtype=bool)
        if ref_mask.shape != inten.shape:
            raise ValueError(
                f"the reference mask is {ref_mask.shape}, the chip is "
                f"{inten.shape} -- they must be the same chip")
        ref = float(np.median(inten[ref_mask])) or 1.0
        ref_source = "ring"
        ref_cells = int(ref_mask.sum())
    mask_roi = band & (roi >= K_LOCAL * ref)
    n_masked = int(mask_roi.sum())
    # Reference filter inside the band (house definition, cross-check only)
    band_vals = roi[band]
    peak_b = float(band_vals.max())
    med_b = float(np.median(band_vals))
    house = band & (roi >= 0.30 * peak_b) & (roi >= 5.0 * med_b)
    band_contrast = float(np.median(band_vals) / ref) if ref > 0 else float("nan")
    full = np.zeros((rows, cols), dtype=bool)
    full[r_lo0:r_hi0 + 1, c_lo:c_hi + 1] = mask_roi
    n_comp, largest = largest_component(mask_roi)
    info = {
        "anchor_row": int(anchor),
        "anchor_col": int(anchor_col),
        "peak_row": int(peak_row),
        "peak_row_offset": int(peak_row - anchor),
        "n_masked": n_masked,
        "band_contrast": band_contrast,
        "ref_source": ref_source,
        "ref_cells": ref_cells,
        "ref_median": float(ref),
        "n_masked_house": int(house.sum()),
        "n_components": n_comp,
        "largest_component": largest,
        "has_mask": (n_masked >= MIN_N_MASKED and band_contrast >= CONTRAST_DETECT),
    }
    if n_masked > 0:
        rr, cc = np.nonzero(mask_roi)
        info["centroid_row"] = float(rr.mean() + r_lo0)
        info["centroid_col"] = float(cc.mean() + c_lo)
    # ---- OCV vector of this mask (the swapped-in mask layer) ----------------
    # Deck-local peak = brightest pixel *inside the band* (see module docstring).
    br, bc = np.unravel_index(np.argmax(np.where(band, roi, -np.inf)), roi.shape)
    info["band_peak_r"] = int(r_lo0 + br)
    info["band_peak_c"] = int(c_lo + bc)
    info["A"] = n_masked
    rr, cc = np.nonzero(full)
    if n_masked > 0:
        row_span = int(rr.max() - rr.min() + 1)
        col_span = int(cc.max() - cc.min() + 1)
        bbox_area = row_span * col_span
        info["row_span"] = row_span
        info["col_span"] = col_span
        info["bbox_area"] = bbox_area
        info["D"] = (n_masked / bbox_area) if bbox_area > 0 else None
        info["F"] = n_comp
        info["S"] = float(np.hypot(info["centroid_row"] - info["band_peak_r"],
                                   info["centroid_col"] - info["band_peak_c"]))
    else:
        info["row_span"] = 0
        info["col_span"] = 0
        info["bbox_area"] = 0
        info["D"] = None
        info["F"] = 0
        info["S"] = None
    return info, full


def analyse_series(recs, ref_mask=None):
    """One mask per date + series summary + union map (verbatim port).

    ``ref_mask`` is forwarded to :func:`date_mask`: ``None`` is the committed
    neighbour-row reference, a ring mask is the v2 reference. The masks, and
    therefore every summary number, follow from it.
    """
    per_date, counts = [], None
    for date, inten in recs:
        info, full = date_mask(inten, ref_mask=ref_mask)
        info["date"] = date
        per_date.append(info)
        counts = full.astype(int) if counts is None else counts + full.astype(int)
    n = len(per_date)
    contrast = np.asarray([d["band_contrast"] for d in per_date], float)
    offsets = np.asarray([abs(d["peak_row_offset"]) for d in per_date], float)
    nm = np.asarray([d["n_masked"] for d in per_date], float)
    with_mask = int(sum(1 for d in per_date if d["has_mask"]))
    # Persistence over all dates: pixels masked in >= 50 % / >= 30 % of them.
    union50 = [(int(r), int(c), float(counts[r, c] / n))
               for r, c in zip(*np.nonzero(counts >= np.ceil(0.5 * n)))] if n else []
    union30 = [(int(r), int(c), float(counts[r, c] / n))
               for r, c in zip(*np.nonzero(counts >= np.ceil(0.3 * n)))] if n else []
    summary = {
        "n_dates": n,
        "n_dates_with_mask": with_mask,
        "frac_dates_with_mask": with_mask / n if n else None,
        "band_contrast_median": float(np.median(contrast)) if n else None,
        "band_contrast_max": float(contrast.max()) if n else None,
        "frac_dates_contrast_gt_thr": float((contrast >= CONTRAST_DETECT).mean()) if n else None,
        "peak_row_median": float(np.median([d["peak_row"] for d in per_date])) if n else None,
        "peak_row_offset_median": float(np.median(offsets)) if n else None,
        "frac_dates_peak_within_1": float((offsets <= 1).mean()) if n else None,
        "n_masked_median": float(np.median(nm)) if n else None,
        "union50_n": len(union50),
        "union30_n": len(union30),
        "n_union_dates": n,
    }
    detect = (n > 0 and summary["band_contrast_median"] >= CONTRAST_DETECT
              and summary["frac_dates_contrast_gt_thr"] >= 0.5
              and with_mask >= 0.5 * n)
    summary["deck_edge_detectable"] = bool(detect)
    return per_date, summary, union50, union30, counts


def add_persistence(rows):
    """P = mean pixel frequency inside the own mask (deck band, per series).

    Verbatim definition of ``lumo_ocv_masks.add_persistence``: a pixel-exact
    frequency map is only comparable across dates with the same chip shape, so
    dates with a deviating shape get ``P = None`` instead of a silent
    misalignment. Unlike the LUMO site, the frequency map is built **per
    series** (bridge/control, ASC/DESC), because the deck mask differs per
    series; ``rows`` must therefore be the dates of one series and must carry
    the keys ``mask`` (full-chip bool array) and ``shape``.
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
        if not m.any():
            # No deck-edge pixel on this date: there is no own mask to average
            # over, so P is undefined (empty in the CSV) -- not 0.0 and not a
            # NaN produced by mean() of an empty slice.
            r["P"] = None
            r["P_frac_ge_50pct"] = None
            r["P_frac_ge_30pct"] = None
            continue
        r["P"] = float(freq[m].mean())
        r["P_frac_ge_50pct"] = float((freq[m] >= 0.5).mean())
        r["P_frac_ge_30pct"] = float((freq[m] >= 0.3).mean())
    return {"majority_shape": list(majority_shape), "n_majority_shape": n_majority,
            "n_excluded_other_shape": n_excluded}


# The two documentation strings of the gamma^2 block are stated once, because
# ``bautzen_ocv_core.LABEL_EN`` maps the committed German texts of the reference
# JSON onto exactly these English ones (the pin block compares them).
GAMMA2_NOTE = ("gamma2_band_* = structurally attributed deck band mask "
               "(deck-edge statement); gamma2_house_* = signal-based "
               "house/peak mask, not structurally attributed "
               "(no deck-edge statement)")
GAMMA2_FLOOR_CORRECTION = ("(N*raw - 1)/(N - 1) under the incoherent null "
                           "model; negative values are kept and mean "
                           "'below the finite-N noise floor'")


# ---------------------------------------------------------------------------
# Per-state aggregation of the deck-edge features — verbatim port of
# analyze_bautzen_deck_edge_vs_state.state_metrics / .discrimination
# ---------------------------------------------------------------------------
def state_metrics(infos, counts):
    """The four deck-edge features of one state (or of all dates of a series)."""
    if not infos:
        return None
    n = len(infos)
    contrast = [d["band_contrast"] for d in infos]
    peak = [d["peak_row"] for d in infos]
    offs = [abs(d["peak_row_offset"]) for d in infos]
    nmask = [d["n_masked"] for d in infos]
    cen = [d["centroid_row"] for d in infos if d.get("centroid_row") is not None]
    n_has = sum(1 for d in infos if d["has_mask"])
    p, p_lo, p_hi = st.wilson(n_has, n)
    c_med, c_q1, c_q3 = st.median_iqr(contrast)
    p_med, p_q1, p_q3 = st.median_iqr(peak)
    o_med, o_q1, o_q3 = st.median_iqr(offs)
    e_med, e_q1, e_q3 = st.median_iqr(cen)
    # Persistence **within** this state
    u50, u30 = [], []
    if counts is not None and n:
        u50 = [(int(r), int(c), float(counts[r, c] / n))
               for r, c in zip(*np.nonzero(counts >= np.ceil(0.5 * n)))]
        u30 = [(int(r), int(c), float(counts[r, c] / n))
               for r, c in zip(*np.nonzero(counts >= np.ceil(0.3 * n)))]
    frac_thr = float(np.mean(np.asarray(contrast) >= CONTRAST_DETECT))
    # 5 -- gamma^2 phase consistency (complementary channel, see module docstring).
    #      `_band_*` = structurally attributed deck band mask (deck-edge statement),
    #      `_house_*` = signal-based house/peak mask (NO deck-edge statement).
    g_block = {
        "note": GAMMA2_NOTE,
        "floor_correction": GAMMA2_FLOOR_CORRECTION,
        "low_n_gate": LOW_N_GATE,
    }
    for tag, key in (("band", "gamma2_band"), ("house", "gamma2_house")):
        ns = [d[f"n_{tag}"] for d in infos if d.get(f"n_{tag}") is not None]
        for kind in ("raw", "floor_corrected"):
            kk = f"{key}_{kind}"
            v = [d[kk] for d in infos if d.get(kk) is not None]
            med, lo, hi = st.boot_median_ci(v)
            g_block[f"{kk}_median"] = med
            g_block[f"{kk}_ci95"] = [lo, hi]
            g_block[f"{kk}_n_dates"] = len(v)
        vc = [d[f"{key}_floor_corrected"] for d in infos
              if d.get(f"{key}_floor_corrected") is not None]
        g_block[f"n_{tag}_median"] = float(np.median(ns)) if ns else float("nan")
        g_block[f"n_{tag}_range"] = [int(min(ns)), int(max(ns))] if ns else None
        g_block[f"n_{tag}_below_{LOW_N_GATE}"] = int(sum(1 for x in ns if x < LOW_N_GATE))
        g_block[f"n_negative_floor_corrected_{tag}"] = int(sum(1 for x in vc if x < 0))
        g_block[f"{key}_floor_corrected_slope_per_date"] = st.ols_slope(vc)
        g_block[f"{key}_floor_corrected_drift"] = st.thirds_drift(vc)
    return {
        "n_dates": n,
        "date_min": min(d["date"] for d in infos),
        "date_max": max(d["date"] for d in infos),
        # 1 -- Visibility
        "visibility_n_with_mask": n_has,
        "visibility": p,
        "visibility_wilson_lo": p_lo,
        "visibility_wilson_hi": p_hi,
        # 2 -- Contrast
        "contrast_median": c_med,
        "contrast_iqr": [c_q1, c_q3],
        "contrast_min": float(np.min(contrast)),
        "contrast_max": float(np.max(contrast)),
        "contrast_frac_gt_thr": frac_thr,
        "contrast_detectable": bool(c_med >= CONTRAST_DETECT and frac_thr >= 0.5
                                    and n_has >= 0.5 * n),
        # 3 -- Position
        "position_peak_median": p_med,
        "position_peak_iqr": [p_q1, p_q3],
        "position_centroid_median": e_med,
        "position_centroid_iqr": [e_q1, e_q3],
        "position_offset_median": o_med,
        "position_offset_iqr": [o_q1, o_q3],
        "position_frac_within_1": float(np.mean(np.asarray(offs) <= 1)),
        "position_drift_px": st.thirds_drift(peak),
        "position_slope_px_per_date": st.ols_slope(peak),
        # 4 -- Persistence
        "persistence_n50": len(u50),
        "persistence_n30": len(u30),
        "persistence_n_masked_median": float(np.median(nmask)),
        # 5 -- gamma^2 (complementary channel)
        "gamma2": g_block,
        "u50": u50,
        "u30": u30,
    }


def discrimination(by_bridge, by_ctrl):
    """Median contrast bridge/control per state, date-paired + bootstrap CI."""
    out = {}
    for s in st.STATE_ORDER:
        mb = {r["date"]: r["band_contrast"] for r in by_bridge[s]}
        mk = {r["date"]: r["band_contrast"] for r in by_ctrl[s]}
        dates = sorted(set(mb) & set(mk))
        if not dates:
            out[s] = None
            continue
        b = np.asarray([mb[d] for d in dates], float)
        k = np.asarray([mk[d] for d in dates], float)
        med_k = float(np.median(k))
        ratio = float(np.median(b) / med_k) if med_k else float("nan")
        rng = np.random.default_rng(st.BOOT_SEED)
        idx = rng.integers(0, len(dates), size=(st.N_BOOT, len(dates)))
        boots = np.median(b[idx], axis=1) / np.median(k[idx], axis=1)
        lo, hi = np.percentile(boots, [2.5, 97.5])
        out[s] = {"n_pairs": len(dates), "bridge_median": float(np.median(b)),
                  "ctrl_median": med_k, "ratio": ratio,
                  "ratio_ci95": [float(lo), float(hi)]}
    return out


def counts_of(infos):
    """Pixel count map (full chip) over the given records.

    ``infos`` must carry the full-chip mask under the key ``mask`` (the records
    of ``date_mask``/``analyse_series`` do). Used for the per-state persistence
    maps, in the same way as ``build_series`` of the origin.
    """
    counts = None
    for d in infos:
        m = d["mask"].astype(int)
        counts = m if counts is None else counts + m
    return counts


def unions(counts, n):
    """``(>= 50 %, >= 30 %)`` persistence pixel lists of a count map.

    Same formulas as ``analyse_series``/``state_metrics`` of this module (and of
    the origin): a pixel is listed when it is masked in at least ceil(f * n)
    dates; the third value is its frequency.
    """
    if counts is None or not n:
        return [], []
    u50 = [(int(r), int(c), float(counts[r, c] / n))
           for r, c in zip(*np.nonzero(counts >= np.ceil(0.5 * n)))]
    u30 = [(int(r), int(c), float(counts[r, c] / n))
           for r, c in zip(*np.nonzero(counts >= np.ceil(0.3 * n)))]
    return u50, u30



