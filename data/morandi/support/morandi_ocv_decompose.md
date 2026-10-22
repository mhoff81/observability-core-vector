# Ponte Morandi — decomposition of the post-collapse signature

* table `data/morandi/support/morandi_ocv_channels.csv` (sha256 `249c9ee6630c21df…`)
* event **2018-08-14**
* strata: ALL
* targets: gamma2, P, D, F, S — `A` is the explanatory variable, never a target
* moving-block bootstrap: block 5 (ladder 3, 5, 7), 1000 draws, seed 7

## How to read the two ladders

`+A` and `A-free` are two parallel ladders over the same data:

```
A-adjusted   M0 (state) -> M1 (+log A) -> M2 (+T) -> M3 (+temperature and season)
A-free       M0 (state) -> M2- (+T)   -> M3- (+temperature and season)
```

A post effect that **attenuates under A adjustment while staying put in the A-free ladder** *co-moves with the observed area*. That is a statement about co-movement, not about cause: the area is itself moved by the collapse (DESC area falls by 77 %), so conditioning on it can remove the very change being measured — and, for a mediator, can also *create* association where the raw data show none. Both directions occur in this record and are flagged in the verdict column.

## The minimal table

| stratum | target | raw | +A | A-free | +T/season | final 95% CI | Holm |
|---|---|---|---|---|---|---|---|
| ALL | `gamma2` | 1.163 | 1.062 | 0.959 | 1.062 | [0.627, 1.445] | 0.756 |
| ALL | `P` | 1.058 | 1.090 | 1.071 | 1.090 | [0.982, 1.167] | 0.288 |
| ALL | `D` | 0.389 | 0.305 | 0.304 | 0.305 | [0.220, 0.515] | <0.001 |
| ALL | `F` | 1.853 | 1.266 | 1.374 | 1.266 | [1.139, 1.380] | <0.001 |
| ALL | `S` | +1.75 | +4.595 | +4.757 | +4.595 | [-0.459, +8.739] | 0.288 |

Ratio columns are `exp(beta_post)` for `gamma2, P, D, F` (adjusted post level as a fraction of pre); `S` is in pixels and shows the signed shift. `+T/season` is the richest A-adjusted model this site could fit (see the tier and the estimability rule above); the Holm column is family A: this site's targets **within** a stratum.

## ALL

*n* = 198 pre / 9 post. Missing values (never imputed):  — none

### Common support (pre 5th-95th percentile)

| variable | pre p5 | pre p95 | post range | post outside | overlap |
|---|---|---|---|---|---|
| `logA` | 2.08 | 6.51 | [2.40, 5.77] | **0.0 %** | 0.47 |
| `T` | 2.56 | 25.03 | [17.90, 27.80] | **33.3 %** | 0.28 |
| `T_smooth30` | 4.14 | 24.35 | [20.17, 25.32] | **22.2 %** | 0.19 |

Common-support subset (the pre p5-p95 box of `log A` **and** `T`): 167 of 207 rows kept, 6 of 9 post rows (67 %).

### Raw, unadjusted (reference only)

| target | n pre/post | median shift (g) | ratio | Cliff's delta [CI] | exp(beta) M0 [CI] |
|---|---|---|---|---|---|
| `gamma2` | 198/9 | +0.151 | 1.163 | -0.08 [-0.25, +0.09] | 0.959 [0.468, 1.734] |
| `P` | 198/9 | +0.057 | 1.058 | +0.08 [-0.08, +0.22] | 1.071 [0.925, 1.223] |
| `D` | 198/9 | -0.945 | 0.389 | -0.43 [-0.52, -0.21] | 0.304 [0.247, 0.509] |
| `F` | 198/9 | +0.617 | 1.853 | +0.15 [+0.00, +0.28] | 1.374 [0.916, 2.048] |
| `S` | 198/9 | +1.748 | — | +0.20 [-0.03, +0.33] | +4.757 [-0.630, +7.081] |

### A adjusted ladder

| target | model | post effect | 95% CI | p (bootstrap) | n pre/post | n |
|---|---|---|---|---|---|---|
| `gamma2` | M0 | 0.959 | [0.389, 1.511] | 0.504 | 198/9 | 207 |
| `gamma2` | M1 | 1.062 | [0.627, 1.445] | 0.756 | 198/9 | 207 |
| `P` | M0 | 1.071 | [0.960, 1.124] | 0.366 | 198/9 | 207 |
| `P` | M1 | 1.090 | [0.982, 1.167] | 0.108 | 198/9 | 207 |
| `D` | M0 | 0.304 | [0.221, 0.513] | <0.001 | 198/9 | 207 |
| `D` | M1 | 0.305 | [0.220, 0.515] | <0.001 | 198/9 | 207 |
| `F` | M0 | 1.374 | [1.075, 1.989] | 0.024 | 198/9 | 207 |
| `F` | M1 | 1.266 | [1.139, 1.380] | <0.001 | 198/9 | 207 |
| `S` | M0 | +4.757 | [+0.395, +8.624] | 0.010 | 198/9 | 207 |
| `S` | M1 | +4.595 | [-0.459, +8.739] | 0.096 | 198/9 | 207 |

### A free ladder

| target | model | post effect | 95% CI | p (bootstrap) | n pre/post | n |
|---|---|---|---|---|---|---|
| `gamma2` | M0 | 0.959 | [0.389, 1.511] | 0.504 | 198/9 | 207 |
| `P` | M0 | 1.071 | [0.960, 1.124] | 0.366 | 198/9 | 207 |
| `D` | M0 | 0.304 | [0.221, 0.513] | <0.001 | 198/9 | 207 |
| `F` | M0 | 1.374 | [1.075, 1.989] | 0.024 | 198/9 | 207 |
| `S` | M0 | +4.757 | [+0.395, +8.624] | 0.010 | 198/9 | 207 |

### Attenuation (descriptive, no threshold)

| target | raw -> +A | raw -> +A+T+season | raw -> A-free final | verdict |
|---|---|---|---|---|
| `gamma2` | -0.41 | -0.41 | 0.00 | not distinguishable from no effect after the covariates |
| `P` | -0.25 | -0.25 | 0.00 | not distinguishable from no effect after the covariates |
| `D` | 0.00 | 0.00 | 0.00 | survives area, temperature and season (both ladders exclude no effect, same sign) |
| `F` | 0.26 | 0.26 | 0.00 | survives area, temperature and season (both ladders exclude no effect, same sign) |
| `S` | 0.03 | 0.03 | 0.00 | co-moves with the observed area (A-free survives, A-adjusted does not) |

### Nested tests (partial F)

| target | +temperature | +season (A-adjusted) | +season (A-free) |
|---|---|---|---|
| `gamma2` | F 0.12 (df 1), p 0.734 | F 0.90 (df 2), p 0.409 | F 6.57 (df 2), p 0.002 |
| `P` | F 2.68 (df 1), p 0.103 | F 11.01 (df 2), p <0.001 | F 19.13 (df 2), p <0.001 |
| `D` | F 17.54 (df 1), p <0.001 | F 11.01 (df 2), p <0.001 | F 11.81 (df 2), p <0.001 |
| `F` | F 6.39 (df 1), p 0.012 | F 11.62 (df 2), p <0.001 | F 17.29 (df 2), p <0.001 |
| `S` | F 0.77 (df 1), p 0.381 | F 9.04 (df 2), p <0.001 | F 15.79 (df 2), p <0.001 |

### Block-length stability of the final post effect

| target | block 3 | block 5 | block 7 | same verdict |
|---|---|---|---|---|
| `gamma2` | 1.010 [0.550, 1.412] | 1.010 [0.622, 1.373] | 1.010 [0.574, 1.255] | yes |
| `P` | 1.172 [1.019, 1.283] | 1.172 [0.961, 1.375] | 1.172 [0.926, 1.370] | **no — unstable** |
| `D` | 0.507 [0.360, 0.866] | 0.507 [0.371, 0.924] | 0.507 [0.349, 0.833] | yes |
| `F` | 1.128 [1.025, 1.227] | 1.128 [1.018, 1.248] | 1.128 [1.021, 1.211] | yes |
| `S` | +2.712 [-2.476, +7.029] | +2.712 [-2.763, +5.285] | +2.712 [-1.192, +6.019] | yes |

### Residual check (model fitted on pre only)

| target | pre median resid | post median resid | delta [CI] |
|---|---|---|---|
| `gamma2` | +0.218 | +0.302 | +0.03 [-0.21, +0.16] |
| `P` | +0.170 | +0.295 | +0.14 [-0.06, +0.33] |
| `D` | -0.152 | -1.087 | -0.42 [-0.52, -0.20] |
| `F` | -0.044 | +0.240 | +0.46 [+0.27, +0.59] |
| `S` | -2.893 | -3.927 | +0.13 [-0.15, +0.24] |

### Residual diagnostics (final model)

| target | lag-1 | runs p | acf 1..4 |
|---|---|---|---|
| `gamma2` | +0.00 | 0.124 | +0.00, -0.04, +0.05, +0.02 |
| `P` | -0.04 | 0.889 | -0.04, +0.15, +0.19, +0.44 |
| `D` | -0.30 | <0.001 | -0.30, +0.08, -0.08, +0.18 |
| `F` | -0.21 | 0.124 | -0.21, -0.08, +0.08, +0.13 |
| `S` | -0.13 | 0.069 | -0.13, +0.04, -0.11, +0.06 |

### Same calendar month, post versus pre

| target | months used | median ratio (log targets) / shift (S) |
|---|---|---|
| `gamma2` | 2 | 1.951 |
| `P` | 2 | 1.530 |
| `D` | 2 | 0.850 |
| `F` | 2 | 1.286 |
| `S` | 2 | +1.30 |

## The whole vector (secondary)

Standardised residuals of `M3-` (temperature and season, no area, no post) fitted on the pre data, with a circular-shift null that keeps the serial dependence and the vector's internal correlation.

| stratum | n | n pre/post | statistic | p (shift) | dimensions |
|---|---|---|---|---|---|
| ALL | 207 | 198/9 | 0.470 | 0.500 | gamma2, P, D, F, S |

## Placebo events (exploratory, family C)

Pseudo events inside the pre period, run through the same estimator. The pseudo-post group is smaller than the real one, so the size-matched rank is the fair comparison; both are shown.

| stratum | spec | real | rank (raw) | matched size | rank (matched) | of |
|---|---|---|---|---|---|---|
| ALL | M0 | 0.987 | 4 | 9 | 4 | 8 |
| ALL | M3_free | 0.389 | 9 | 9 | 9 | 8 |

A rank of 1 means the real event produced a stronger signature than every pseudo event; a rank of n+1 means all pseudo events were stronger.

## Holm across all target x stratum post effects (sensitivity)

Family size 5; values in `report['holm_all_strata']`.

| stratum | target | Holm p (all 15) |
|---|---|---|
| ALL | `D` | <0.001 |
| ALL | `F` | <0.001 |
| ALL | `P` | 0.288 |
| ALL | `S` | 0.288 |
| ALL | `gamma2` | 0.756 |

## Decision logic, applied

| stratum | survives | co-moves with A | sign flip | no effect | unstable |
|---|---|---|---|---|---|
| ALL | 2 (D, F) | 1 (S) | 0 | 2 (gamma2, P) | 0 |

The verdict is about the **covariates**, not about multiplicity: a cell can survive the area/temperature/season control and still miss its Holm threshold (see the last column of the minimal table). Read the two together.

## What this does not say

* `A` is not a neutral confounder here. The area is part of what the collapse changes, so the A-adjusted column measures **co-movement with the observed area**, not a correction of an acquisition artefact. Where the two ladders disagree in sign, the A-adjusted estimate is the one to distrust: conditioning on a variable that the event moves can create association as easily as it removes it.
* `T24` and `T72` do not exist in this record (one temperature per overpass). `T0` is the primary temperature and `T_smooth30` the sensitivity; neither is a measurement of the deck temperature.
* `Post` is a step at a single date. These numbers describe this one collapse; they are not a causal effect of collapsing and not a forecast.
* The placebos live inside the pre period and cannot reproduce the real group sizes; the size-matched rank is reported beside the raw one.
* Descending geometry has 56 % of its post area below the pre 5th percentile: its A-adjusted estimate is an **extrapolation** along the area relation, and the common-support subset is the number to read.

