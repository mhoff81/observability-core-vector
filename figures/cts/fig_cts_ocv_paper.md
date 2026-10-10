# Champlain Towers South, Surfside (FL) — the observability core vector around the collapse

**Site** Champlain Towers South, Surfside (FL) · **event** 2021-06-24 · **chips** 711 (543 pre-collapse (healthy) / 168 post-collapse)

```
state = pre-collapse (healthy)  day <  2021-06-24
        post-collapse            day >= 2021-06-24  (bare sand)
x     = [gamma2, P, D, A, F, S]   echo-mask vector
OCV   = [gamma2, P, D]            observability core vector
tracks = [cts (target), ctn, cte (controls), beach (reference)]
mask  = intensity >= 0.3*peak AND >= 5.0*np.median(intensity), A >= 2
DiD   = y ~ 1 + T + P + T:P over the pooled target/reference rows
```

## Honest null — read this first

**168** of **711** chips follow the collapse against **543** before it, so every post-state interval is wide. CTS is the **difference-in-differences** site: the post state is *bare sand* whose intensity is moisture-driven, so a change in the footprint is not evidence about the collapse — panel F differences it against the two control towers and the CTN-vs-CTE placebo. No channel is significant at the 5 % level in either contrast, and several channels' own reference steps; the honest result is a null, and it must never be read as a damage verdict.

## The difference-in-differences layer (ASCENDING)

The three committed references are re-derived by this script and pinned to 1e-9 (see `cts_ocv_did.py --verify`). A significant DiD in a channel whose **reference also steps** (`reference steps p` small) is not a DiD result.

### target − CTN  (`did_ctn`)

| channel | n | target Δ | reference Δ | DiD | t | p | pre comparable p | reference steps p | pre-trend ρ tgt/ref |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `phase_coherence` | 356 | -0.001187 | -0.002953 | **+0.002878** | +1.45 | 0.148 | 0.957 | 0.004 | +0.16 / -0.07 |
| `phase_snr_db` | 0 | — | — | — | — | — | — | — | eine der vier Zellen hat < 5 Werte |
| `phase_scatterer_count` | 356 | -12.5 | -16 | **-9.69** | -0.71 | 0.480 | 1.000 | 0.574 | +0.05 / -0.05 |
| `displacement_los_m` | 0 | — | — | — | — | — | — | — | eine der vier Zellen hat < 5 Werte |
| `coherence_gamma2` | 356 | -3.625e-05 | -4.79e-05 | **+4.164e-05** | +0.53 | 0.596 | 0.191 | 0.017 | +0.10 / -0.00 |
| `coherence` | 350 | -0.02372 | -0.004949 | **-0.02116** | -0.65 | 0.518 | 0.667 | 0.090 | +0.25 / -0.00 |
| `coherence_masked_pixels` | 356 | +0 | +0 | -5.07e-14* | — | — | n/a | n/a | n/a / n/a |
| `intensity` | 356 | +643.3 | +3499 | **-2830** | -1.88 | 0.060 | 0.000 | 0.004 | -0.08 / +0.07 |
| `brightness_ratio` | 356 | +0.016 | +0.07695 | **-0.06023** | -1.71 | 0.088 | 0.755 | 0.004 | -0.08 / +0.07 |

### target − CTE  (`did_cte`)

| channel | n | target Δ | reference Δ | DiD | t | p | pre comparable p | reference steps p | pre-trend ρ tgt/ref |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `phase_coherence` | 356 | -0.001187 | -0.003351 | **+0.001819** | +0.90 | 0.369 | 0.875 | 0.073 | +0.16 / -0.05 |
| `phase_snr_db` | 0 | — | — | — | — | — | — | — | eine der vier Zellen hat < 5 Werte |
| `phase_scatterer_count` | 356 | -12.5 | -29 | **-2.5** | -0.19 | 0.850 | 1.000 | 0.167 | +0.05 / -0.09 |
| `displacement_los_m` | 0 | — | — | — | — | — | — | — | eine der vier Zellen hat < 5 Werte |
| `coherence_gamma2` | 356 | -3.625e-05 | -2.836e-05 | **-1.14e-05** | -0.14 | 0.890 | 0.671 | 0.394 | +0.10 / -0.00 |
| `coherence` | 350 | -0.02372 | -0.009867 | **-0.008225** | -0.24 | 0.807 | 0.904 | 0.029 | +0.25 / +0.09 |
| `coherence_masked_pixels` | 356 | +0 | +0 | -5.07e-14* | — | — | n/a | n/a | n/a / n/a |
| `intensity` | 356 | +643.3 | +145.2 | **-118.6** | -0.08 | 0.934 | 0.838 | 0.469 | -0.08 / -0.13 |
| `brightness_ratio` | 356 | +0.016 | +0.003615 | **-0.002971** | -0.08 | 0.934 | 0.886 | 0.469 | -0.08 / -0.13 |

### placebo CTN − CTE  (`did_placebo_ctn_vs_cte`)

| channel | n | target Δ | reference Δ | DiD | t | p | pre comparable p | reference steps p | pre-trend ρ tgt/ref |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `phase_coherence` | 356 | -0.002953 | -0.003351 | **-0.001059** | -0.54 | 0.587 | 0.828 | 0.073 | -0.07 / -0.05 |
| `phase_snr_db` | 0 | — | — | — | — | — | — | — | eine der vier Zellen hat < 5 Werte |
| `phase_scatterer_count` | 356 | -16 | -29 | **+7.19** | +0.52 | 0.603 | 1.000 | 0.167 | -0.05 / -0.09 |
| `displacement_los_m` | 0 | — | — | — | — | — | — | — | eine der vier Zellen hat < 5 Werte |
| `coherence_gamma2` | 356 | -4.79e-05 | -2.836e-05 | **-5.304e-05** | -0.71 | 0.478 | 0.352 | 0.394 | -0.00 / -0.00 |
| `coherence` | 350 | -0.004949 | -0.009867 | **+0.01294** | +0.38 | 0.703 | 0.591 | 0.029 | -0.00 / +0.09 |
| `coherence_masked_pixels` | 356 | +0 | +0 | -5.07e-14* | — | — | n/a | n/a | n/a / n/a |
| `intensity` | 356 | +3499 | +145.2 | **+2712** | +1.84 | 0.067 | 0.000 | 0.469 | +0.07 / -0.13 |
| `brightness_ratio` | 356 | +0.07695 | +0.003615 | **+0.05726** | +1.66 | 0.097 | 0.866 | 0.469 | +0.07 / -0.13 |

`*` `coherence_masked_pixels` is the site's **fixed 640-px quantile mask** (constant for every row), so its OLS interaction is floating-point noise: its t/p are not reported and only its structure is pinned.

## Effect sizes of the six channels (post − pre)

| channel | Cliff's delta | 95 % CI | CI≠0 | SDS |
| --- | --- | --- | --- | --- |
| gamma2 | -0.09298 | [-0.1901, 0.007696] | no | 0.9298 |
| P | -0.1678 | [-0.2708, -0.06259] | yes | 1.678 |
| D | 0.2756 | [0.1814, 0.3672] | yes | 2.756 |
| A | 0.1672 | [0.05651, 0.2751] | yes | 1.672 |
| F | 0.1532 | [0.04419, 0.2619] | yes | 1.532 |
| S | -0.1184 | [-0.2169, -0.01778] | yes | 1.184 |

## Per-footprint effect sizes

**cts_asc_48_iw3** — target

| channel | Cliff's delta | 95 % CI | CI≠0 |
| --- | --- | --- | --- |
| gamma2 | -0.1033 | [-0.3074, 0.124] | no |
| P | -0.3428 | [-0.5035, -0.1346] | yes |
| D | 0.3732 | [0.1645, 0.529] | yes |
| A | 0.3643 | [0.1691, 0.5447] | yes |
| F | 0.3704 | [0.1764, 0.5495] | yes |
| S | -0.001401 | [-0.2106, 0.1847] | no |

**ctn_asc_48_iw3** — control_ctn

| channel | Cliff's delta | 95 % CI | CI≠0 |
| --- | --- | --- | --- |
| gamma2 | -0.09314 | [-0.2982, 0.1024] | no |
| P | 0.1766 | [-0.04321, 0.3901] | no |
| D | 0.3382 | [0.1268, 0.509] | yes |
| A | -0.1991 | [-0.3743, -0.00822] | yes |
| F | -0.2705 | [-0.4264, -0.08129] | yes |
| S | -0.2768 | [-0.4236, -0.08451] | yes |

**cte_asc_48_iw3** — control_cte

| channel | Cliff's delta | 95 % CI | CI≠0 |
| --- | --- | --- | --- |
| gamma2 | -0.1326 | [-0.3564, 0.07972] | no |
| P | -0.3159 | [-0.4835, -0.1339] | yes |
| D | 0.1402 | [-0.06667, 0.3418] | no |
| A | 0.3453 | [0.1286, 0.5493] | yes |
| F | 0.2963 | [0.09194, 0.4864] | yes |
| S | -0.1416 | [-0.3356, 0.05441] | no |

**beach_asc_48_iw3** — reference_beach

| channel | Cliff's delta | 95 % CI | CI≠0 |
| --- | --- | --- | --- |
| gamma2 | -0.04237 | [-0.2493, 0.1633] | no |
| P | -0.2794 | [-0.464, -0.07466] | yes |
| D | 0.2349 | [0.006167, 0.4228] | yes |
| A | 0.2654 | [0.06766, 0.4556] | yes |
| F | 0.2516 | [0.05808, 0.4429] | yes |
| S | 0.1005 | [-0.09206, 0.2967] | no |

## Is the state learnable? (2-class LDA, leave-one-out)

| model | dims | LOO acc. | balanced acc. | p (perm.) |
| --- | --- | --- | --- | --- |
| gamma2 | 1 | 0.584 | 0.436 | – |
| PD | 2 | 0.678 | 0.602 | – |
| gamma2PD | 3 | 0.677 | 0.605 | – |
| x_6d | 6 | 0.627 | 0.626 | – |
| weather | 2 | 0.696 | 0.579 | – |

## What this does and does not say

Both sides of the event are the *same* footprint, so season, weather and pipeline generation move with the state — and the post side is bare sand. Panel F is the only reading that removes the *common* confounders (atmosphere, soil moisture, orbits): it is a null, and no contrast here is a damage verdict. The `A` column of the site's export is a *fixed 640-px quantile mask*, not the echo mask this package recomputes — the two never coincide (panel E).

```bash
python3 code/cts/cts_ocv_channels_csv.py   # build + verify the table
python3 code/cts/cts_ocv_did.py --verify    # pin the three DiD references
python3 code/cts/fig_cts_ocv_paper.py      # figures A–F + this report
```
