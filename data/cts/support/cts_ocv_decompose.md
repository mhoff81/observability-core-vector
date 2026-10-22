# CTS deck — decomposition of the post-collapse signature

* table `data/cts/support/cts_ocv_channels.csv` (sha256 `3965a15b2053dc0b…`)
* event **2021-06-24**
* strata: ASC
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
| ASC | `gamma2` | 0.858 | 1.250 | 1.147 | 1.255 | [0.736, 1.899] | 1.000 |
| ASC | `P` | 0.880 | 0.925 | 0.904 | 0.933 | [0.880, 0.997] | 0.210 |
| ASC | `D` | 1.185 | 1.010 | 1.085 | 1.006 | [0.975, 1.044] | 1.000 |
| ASC | `F` | 1.223 | 1.004 | 1.073 | 1.006 | [0.979, 1.044] | 1.000 |
| ASC | `S` | -3.25 | -1.734 | -2.308 | -1.908 | [-5.384, +1.709] | 1.000 |

Ratio columns are `exp(beta_post)` for `gamma2, P, D, F` (adjusted post level as a fraction of pre); `S` is in pixels and shows the signed shift. `+T/season` is the richest A-adjusted model this site could fit (see the tier and the estimability rule above); the Holm column is family A: this site's targets **within** a stratum.

## ASC

*n* = 136 pre / 42 post. Missing values (never imputed):  — none

### Common support (pre 5th-95th percentile)

| variable | pre p5 | pre p95 | post range | post outside | overlap |
|---|---|---|---|---|---|
| `logA` | 4.14 | 6.07 | [3.81, 6.29] | **11.9 %** | 0.69 |
| `T` | 18.23 | 28.20 | [7.90, 28.00] | **4.8 %** | 0.82 |
| `T_smooth30` | 19.57 | 28.03 | [14.30, 27.70] | **4.8 %** | 0.75 |

Common-support subset (the pre p5-p95 box of `log A` **and** `T`): 151 of 178 rows kept, 36 of 42 post rows (86 %).

### Raw, unadjusted (reference only)

| target | n pre/post | median shift (g) | ratio | Cliff's delta [CI] | exp(beta) M0 [CI] |
|---|---|---|---|---|---|
| `gamma2` | 136/42 | -0.154 | 0.858 | +0.02 [-0.15, +0.20] | 1.141 [0.720, 1.803] |
| `P` | 136/42 | -0.128 | 0.880 | -0.23 [-0.36, -0.06] | 0.898 [0.842, 0.974] |
| `D` | 136/42 | +0.170 | 1.185 | +0.14 [-0.04, +0.29] | 1.087 [0.936, 1.242] |
| `F` | 136/42 | +0.201 | 1.223 | +0.12 [-0.07, +0.29] | 1.068 [0.938, 1.204] |
| `S` | 136/42 | -3.253 | — | -0.17 [-0.43, +0.11] | -2.110 [-5.690, +1.748] |

### A adjusted ladder

| target | model | post effect | 95% CI | p (bootstrap) | n pre/post | n |
|---|---|---|---|---|---|---|
| `gamma2` | M0 | 1.141 | [0.720, 1.803] | 0.630 | 136/42 | 178 |
| `gamma2` | M1 | 1.250 | [0.768, 1.860] | 0.430 | 136/42 | 178 |
| `gamma2` | M2 | 1.249 | [0.746, 1.828] | 0.510 | 136/42 | 178 |
| `gamma2` | M3 | 1.255 | [0.736, 1.899] | 0.490 | 136/42 | 178 |
| `P` | M0 | 0.898 | [0.842, 0.974] | 0.004 | 136/42 | 178 |
| `P` | M1 | 0.925 | [0.875, 0.979] | 0.006 | 136/42 | 178 |
| `P` | M2 | 0.925 | [0.878, 0.980] | 0.004 | 136/42 | 178 |
| `P` | M3 | 0.933 | [0.880, 0.997] | 0.042 | 136/42 | 178 |
| `D` | M0 | 1.087 | [0.936, 1.242] | 0.336 | 136/42 | 178 |
| `D` | M1 | 1.010 | [0.985, 1.052] | 0.344 | 136/42 | 178 |
| `D` | M2 | 1.010 | [0.984, 1.049] | 0.360 | 136/42 | 178 |
| `D` | M3 | 1.006 | [0.975, 1.044] | 0.502 | 136/42 | 178 |
| `F` | M0 | 1.068 | [0.938, 1.204] | 0.392 | 136/42 | 178 |
| `F` | M1 | 1.004 | [0.975, 1.039] | 0.568 | 136/42 | 178 |
| `F` | M2 | 1.004 | [0.977, 1.039] | 0.566 | 136/42 | 178 |
| `F` | M3 | 1.006 | [0.979, 1.044] | 0.476 | 136/42 | 178 |
| `S` | M0 | -2.110 | [-5.690, +1.748] | 0.252 | 136/42 | 178 |
| `S` | M1 | -1.734 | [-5.169, +1.792] | 0.320 | 136/42 | 178 |
| `S` | M2 | -1.732 | [-5.127, +1.934] | 0.336 | 136/42 | 178 |
| `S` | M3 | -1.908 | [-5.384, +1.709] | 0.290 | 136/42 | 178 |

### A free ladder

| target | model | post effect | 95% CI | p (bootstrap) | n pre/post | n |
|---|---|---|---|---|---|---|
| `gamma2` | M0 | 1.141 | [0.720, 1.803] | 0.630 | 136/42 | 178 |
| `gamma2` | M2 | 1.143 | [0.710, 1.770] | 0.702 | 136/42 | 178 |
| `gamma2` | M3 | 1.147 | [0.701, 1.796] | 0.684 | 136/42 | 178 |
| `P` | M0 | 0.898 | [0.842, 0.974] | 0.004 | 136/42 | 178 |
| `P` | M2 | 0.897 | [0.843, 0.974] | 0.004 | 136/42 | 178 |
| `P` | M3 | 0.904 | [0.843, 0.995] | 0.046 | 136/42 | 178 |
| `D` | M0 | 1.087 | [0.936, 1.242] | 0.336 | 136/42 | 178 |
| `D` | M2 | 1.086 | [0.935, 1.244] | 0.334 | 136/42 | 178 |
| `D` | M3 | 1.085 | [0.922, 1.242] | 0.360 | 136/42 | 178 |
| `F` | M0 | 1.068 | [0.938, 1.204] | 0.392 | 136/42 | 178 |
| `F` | M2 | 1.068 | [0.941, 1.205] | 0.368 | 136/42 | 178 |
| `F` | M3 | 1.073 | [0.930, 1.201] | 0.340 | 136/42 | 178 |
| `S` | M0 | -2.110 | [-5.690, +1.748] | 0.252 | 136/42 | 178 |
| `S` | M2 | -2.114 | [-5.668, +1.856] | 0.256 | 136/42 | 178 |
| `S` | M3 | -2.308 | [-5.820, +1.807] | 0.256 | 136/42 | 178 |

### Attenuation (descriptive, no threshold)

| target | raw -> +A | raw -> +A+T+season | raw -> A-free final | verdict |
|---|---|---|---|---|
| `gamma2` | -0.69 | -0.72 | -0.04 | not distinguishable from no effect after the covariates |
| `P` | 0.27 | 0.36 | 0.07 | survives area, temperature and season (both ladders exclude no effect, same sign) |
| `D` | 0.88 | 0.92 | 0.02 | not distinguishable from no effect after the covariates |
| `F` | 0.94 | 0.92 | -0.06 | not distinguishable from no effect after the covariates |
| `S` | 0.18 | 0.10 | -0.09 | not distinguishable from no effect after the covariates |

### Nested tests (partial F)

| target | +temperature | +season (A-adjusted) | +season (A-free) |
|---|---|---|---|
| `gamma2` | F 1.81 (df 1), p 0.180 | F 0.78 (df 2), p 0.459 | F 1.60 (df 2), p 0.204 |
| `P` | F 1.45 (df 1), p 0.229 | F 2.43 (df 2), p 0.091 | F 0.53 (df 2), p 0.591 |
| `D` | F 0.75 (df 1), p 0.387 | F 1.36 (df 2), p 0.260 | F 1.13 (df 2), p 0.325 |
| `F` | F 0.00 (df 1), p 0.989 | F 0.47 (df 2), p 0.627 | F 1.19 (df 2), p 0.306 |
| `S` | F 1.48 (df 1), p 0.225 | F 0.64 (df 2), p 0.530 | F 0.49 (df 2), p 0.611 |

### Block-length stability of the final post effect

| target | block 3 | block 5 | block 7 | same verdict |
|---|---|---|---|---|
| `gamma2` | 1.255 [0.768, 1.966] | 1.255 [0.736, 1.899] | 1.255 [0.715, 1.985] | yes |
| `P` | 0.933 [0.875, 0.997] | 0.933 [0.880, 0.997] | 0.933 [0.882, 0.991] | yes |
| `D` | 1.006 [0.975, 1.047] | 1.006 [0.975, 1.044] | 1.006 [0.976, 1.047] | yes |
| `F` | 1.006 [0.981, 1.040] | 1.006 [0.979, 1.044] | 1.006 [0.976, 1.045] | yes |
| `S` | -1.908 [-4.818, +1.307] | -1.908 [-5.384, +1.709] | -1.908 [-5.460, +1.619] | yes |

### Residual check (model fitted on pre only)

| target | pre median resid | post median resid | delta [CI] |
|---|---|---|---|
| `gamma2` | +0.383 | +0.613 | +0.08 [-0.11, +0.24] |
| `P` | +0.062 | -0.029 | -0.35 [-0.54, -0.16] |
| `D` | -0.008 | +0.037 | +0.17 [+0.02, +0.38] |
| `F` | +0.000 | +0.013 | +0.01 [-0.17, +0.25] |
| `S` | +1.901 | -0.049 | -0.13 [-0.40, +0.14] |

### Residual diagnostics (final model)

| target | lag-1 | runs p | acf 1..4 |
|---|---|---|---|
| `gamma2` | +0.04 | 0.881 | +0.04, +0.01, +0.08, -0.04 |
| `P` | +0.17 | 0.007 | +0.17, -0.04, -0.03, +0.05 |
| `D` | +0.03 | 1.000 | +0.03, +0.07, +0.01, -0.07 |
| `F` | -0.04 | 0.548 | -0.04, +0.08, +0.04, +0.20 |
| `S` | +0.13 | 0.016 | +0.13, -0.04, +0.10, +0.04 |

### Same calendar month, post versus pre

| target | months used | median ratio (log targets) / shift (S) |
|---|---|---|
| `gamma2` | 9 | 1.041 |
| `P` | 9 | 0.863 |
| `D` | 9 | 1.056 |
| `F` | 9 | 1.034 |
| `S` | 9 | -2.17 |

## The whole vector (secondary)

Standardised residuals of `M3-` (temperature and season, no area, no post) fitted on the pre data, with a circular-shift null that keeps the serial dependence and the vector's internal correlation.

| stratum | n | n pre/post | statistic | p (shift) | dimensions |
|---|---|---|---|---|---|
| ASC | 178 | 136/42 | 0.267 | 0.481 | gamma2, P, D, F, S |

## Placebo events (exploratory, family C)

Pseudo events inside the pre period, run through the same estimator. The pseudo-post group is smaller than the real one, so the size-matched rank is the fair comparison; both are shown.

| stratum | spec | real | rank (raw) | matched size | rank (matched) | of |
|---|---|---|---|---|---|---|
| ASC | M0 | 0.292 | 2 | 42 | 2 | 8 |
| ASC | M3_free | 0.290 | 1 | 42 | 1 | 8 |

A rank of 1 means the real event produced a stronger signature than every pseudo event; a rank of n+1 means all pseudo events were stronger.

## Holm across all target x stratum post effects (sensitivity)

Family size 5; values in `report['holm_all_strata']`.

| stratum | target | Holm p (all 15) |
|---|---|---|
| ASC | `D` | 1.000 |
| ASC | `F` | 1.000 |
| ASC | `P` | 0.210 |
| ASC | `S` | 1.000 |
| ASC | `gamma2` | 1.000 |

## Decision logic, applied

| stratum | survives | co-moves with A | sign flip | no effect | unstable |
|---|---|---|---|---|---|
| ASC | 1 (P) | 0 | 0 | 4 (gamma2, D, F, S) | 0 |

The verdict is about the **covariates**, not about multiplicity: a cell can survive the area/temperature/season control and still miss its Holm threshold (see the last column of the minimal table). Read the two together.

## What this does not say

* `A` is not a neutral confounder here. The area is part of what the collapse changes, so the A-adjusted column measures **co-movement with the observed area**, not a correction of an acquisition artefact. Where the two ladders disagree in sign, the A-adjusted estimate is the one to distrust: conditioning on a variable that the event moves can create association as easily as it removes it.
* `T24` and `T72` do not exist in this record (one temperature per overpass). `T0` is the primary temperature and `T_smooth30` the sensitivity; neither is a measurement of the deck temperature.
* `Post` is a step at a single date. These numbers describe this one collapse; they are not a causal effect of collapsing and not a forecast.
* The placebos live inside the pre period and cannot reproduce the real group sizes; the size-matched rank is reported beside the raw one.
* Descending geometry has 56 % of its post area below the pre 5th percentile: its A-adjusted estimate is an **extrapolation** along the area relation, and the common-support subset is the number to read.

