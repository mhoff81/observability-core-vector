# Haltern canal bridge — decomposition of the post-collapse signature

* table `data/haltern/support/haltern_ocv_channels.csv` (sha256 `cb39c7473b291540…`)
* event **2021-09-29**
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
| ALL | `gamma2` | 0.742 | 1.397 | 1.208 | 1.397 | [1.091, 2.516] | 0.090 |
| ALL | `P` | 1.015 | 0.997 | 1.001 | 0.997 | [0.971, 1.029] | 0.974 |
| ALL | `D` | 1.139 | 0.961 | 1.142 | 0.961 | [0.930, 1.000] | 0.200 |
| ALL | `F` | 1.103 | 1.071 | 1.168 | 1.071 | [0.995, 1.099] | 0.210 |
| ALL | `S` | -2.75 | +4.221 | +3.800 | +4.221 | [-5.058, +13.454] | 0.592 |

Ratio columns are `exp(beta_post)` for `gamma2, P, D, F` (adjusted post level as a fraction of pre); `S` is in pixels and shows the signed shift. `+T/season` is the richest A-adjusted model this site could fit (see the tier and the estimability rule above); the Holm column is family A: this site's targets **within** a stratum.

## ALL

*n* = 96 pre / 11 post. Missing values (never imputed):  — none

### Common support (pre 5th-95th percentile)

| variable | pre p5 | pre p95 | post range | post outside | overlap |
|---|---|---|---|---|---|
| `logA` | 5.31 | 6.64 | [5.65, 6.70] | **27.3 %** | 0.55 |
| `T` | -1.00 | 23.05 | [6.50, 15.90] | **0.0 %** | 0.46 |
| `T_smooth30` | 1.15 | 20.09 | [10.68, 11.66] | **0.0 %** | 0.01 |

Common-support subset (the pre p5-p95 box of `log A` **and** `T`): 86 of 107 rows kept, 8 of 11 post rows (73 %).

### Raw, unadjusted (reference only)

| target | n pre/post | median shift (g) | ratio | Cliff's delta [CI] | exp(beta) M0 [CI] |
|---|---|---|---|---|---|
| `gamma2` | 96/11 | -0.298 | 0.742 | +0.03 [-0.19, +0.35] | 1.207 [0.840, 2.296] |
| `P` | 96/11 | +0.015 | 1.015 | +0.12 [-0.08, +0.29] | 1.000 [0.970, 1.029] |
| `D` | 96/11 | +0.130 | 1.139 | +0.08 [-0.28, +0.31] | 1.137 [0.945, 1.354] |
| `F` | 96/11 | +0.098 | 1.103 | +0.37 [+0.01, +0.54] | 1.169 [1.057, 1.241] |
| `S` | 96/11 | -2.748 | — | +0.12 [-0.16, +0.41] | +3.851 [-4.831, +14.333] |

### A adjusted ladder

| target | model | post effect | 95% CI | p (bootstrap) | n pre/post | n |
|---|---|---|---|---|---|---|
| `gamma2` | M0 | 1.207 | [0.912, 2.296] | 0.118 | 96/11 | 107 |
| `gamma2` | M1 | 1.397 | [1.091, 2.516] | 0.018 | 96/11 | 107 |
| `P` | M0 | 1.000 | [0.970, 1.031] | 0.804 | 96/11 | 107 |
| `P` | M1 | 0.997 | [0.971, 1.029] | 0.974 | 96/11 | 107 |
| `D` | M0 | 1.137 | [0.934, 1.367] | 0.172 | 96/11 | 107 |
| `D` | M1 | 0.961 | [0.930, 1.000] | 0.050 | 96/11 | 107 |
| `F` | M0 | 1.169 | [1.056, 1.258] | 0.002 | 96/11 | 107 |
| `F` | M1 | 1.071 | [0.995, 1.099] | 0.070 | 96/11 | 107 |
| `S` | M0 | +3.851 | [-5.102, +13.299] | 0.332 | 96/11 | 107 |
| `S` | M1 | +4.221 | [-5.058, +13.454] | 0.296 | 96/11 | 107 |

### A free ladder

| target | model | post effect | 95% CI | p (bootstrap) | n pre/post | n |
|---|---|---|---|---|---|---|
| `gamma2` | M0 | 1.207 | [0.912, 2.296] | 0.118 | 96/11 | 107 |
| `gamma2` | M2 | 1.208 | [0.914, 2.283] | 0.130 | 96/11 | 107 |
| `P` | M0 | 1.000 | [0.970, 1.031] | 0.804 | 96/11 | 107 |
| `P` | M2 | 1.001 | [0.969, 1.033] | 0.816 | 96/11 | 107 |
| `D` | M0 | 1.137 | [0.934, 1.367] | 0.172 | 96/11 | 107 |
| `D` | M2 | 1.142 | [0.930, 1.370] | 0.164 | 96/11 | 107 |
| `F` | M0 | 1.169 | [1.056, 1.258] | 0.002 | 96/11 | 107 |
| `F` | M2 | 1.168 | [1.054, 1.259] | 0.002 | 96/11 | 107 |
| `S` | M0 | +3.851 | [-5.102, +13.299] | 0.332 | 96/11 | 107 |
| `S` | M2 | +3.800 | [-4.956, +13.248] | 0.342 | 96/11 | 107 |

### Attenuation (descriptive, no threshold)

| target | raw -> +A | raw -> +A+T+season | raw -> A-free final | verdict |
|---|---|---|---|---|
| `gamma2` | -0.78 | -0.78 | -0.00 | unstable (A-adjusted survives, A-free does not) |
| `P` | -9.89 | -9.89 | -3.60 | not distinguishable from no effect after the covariates |
| `D` | 0.69 | 0.69 | -0.03 | unstable (A-adjusted survives, A-free does not) |
| `F` | 0.56 | 0.56 | 0.00 | co-moves with the observed area (A-free survives, A-adjusted does not) |
| `S` | -0.10 | -0.10 | 0.01 | not distinguishable from no effect after the covariates |

### Nested tests (partial F)

| target | +temperature | +season (A-adjusted) | +season (A-free) |
|---|---|---|---|
| `gamma2` | F 0.07 (df 1), p 0.791 | F 0.18 (df 2), p 0.837 | F 0.01 (df 2), p 0.987 |
| `P` | F 1.46 (df 1), p 0.230 | F 4.90 (df 2), p 0.009 | F 5.02 (df 2), p 0.008 |
| `D` | F 5.18 (df 1), p 0.025 | F 4.75 (df 2), p 0.011 | F 7.33 (df 2), p 0.001 |
| `F` | F 1.71 (df 1), p 0.195 | F 2.70 (df 2), p 0.072 | F 0.61 (df 2), p 0.544 |
| `S` | F 0.21 (df 1), p 0.652 | F 0.45 (df 2), p 0.637 | F 0.37 (df 2), p 0.692 |

### Block-length stability of the final post effect

| target | block 3 | block 5 | block 7 | same verdict |
|---|---|---|---|---|
| `gamma2` | 1.290 [0.817, 3.045] | 1.290 [0.659, 3.718] | 1.290 [0.674, 3.860] | yes |
| `P` | 1.002 [0.956, 1.050] | 1.002 [0.944, 1.039] | 1.002 [0.933, 1.023] | yes |
| `D` | 0.974 [0.922, 1.031] | 0.974 [0.905, 1.032] | 0.974 [0.897, 1.029] | yes |
| `F` | 1.039 [0.945, 1.102] | 1.039 [0.952, 1.143] | 1.039 [0.978, 1.134] | yes |
| `S` | +6.050 [-6.438, +18.838] | +6.050 [-4.666, +20.619] | +6.050 [-0.012, +22.356] | yes |

### Residual check (model fitted on pre only)

| target | pre median resid | post median resid | delta [CI] |
|---|---|---|---|
| `gamma2` | +0.243 | +0.563 | +0.12 [-0.04, +0.44] |
| `P` | -0.042 | -0.026 | +0.12 [-0.09, +0.28] |
| `D` | -0.075 | -0.079 | -0.08 [-0.35, +0.06] |
| `F` | +0.099 | +0.116 | +0.20 [-0.05, +0.38] |
| `S` | -2.121 | -4.148 | +0.14 [-0.14, +0.42] |

### Residual diagnostics (final model)

| target | lag-1 | runs p | acf 1..4 |
|---|---|---|---|
| `gamma2` | -0.03 | 0.118 | -0.03, -0.09, -0.03, -0.08 |
| `P` | -0.19 | 0.329 | -0.19, -0.27, +0.12, +0.09 |
| `D` | -0.26 | 0.435 | -0.26, -0.20, +0.09, +0.14 |
| `F` | -0.23 | 0.435 | -0.23, -0.40, +0.29, +0.17 |
| `S` | -0.07 | 0.696 | -0.07, -0.15, -0.06, -0.11 |

### Same calendar month, post versus pre

| target | months used | median ratio (log targets) / shift (S) |
|---|---|---|
| `gamma2` | — | — |
| `P` | — | — |
| `D` | — | — |
| `F` | — | — |
| `S` | — | — |

## The whole vector (secondary)

Standardised residuals of `M3-` (temperature and season, no area, no post) fitted on the pre data, with a circular-shift null that keeps the serial dependence and the vector's internal correlation.

| stratum | n | n pre/post | statistic | p (shift) | dimensions |
|---|---|---|---|---|---|
| ALL | 107 | 96/11 | 0.262 | 0.468 | gamma2, P, D, F, S |

## Placebo events (exploratory, family C)

Pseudo events inside the pre period, run through the same estimator. The pseudo-post group is smaller than the real one, so the size-matched rank is the fair comparison; both are shown.

| stratum | spec | real | rank (raw) | matched size | rank (matched) | of |
|---|---|---|---|---|---|---|
| ALL | M0 | 0.475 | 6 | 11 | 6 | 8 |
| ALL | M3_free | 0.258 | 8 | 11 | 8 | 8 |

A rank of 1 means the real event produced a stronger signature than every pseudo event; a rank of n+1 means all pseudo events were stronger.

## Holm across all target x stratum post effects (sensitivity)

Family size 5; values in `report['holm_all_strata']`.

| stratum | target | Holm p (all 15) |
|---|---|---|
| ALL | `D` | 0.200 |
| ALL | `F` | 0.210 |
| ALL | `P` | 0.974 |
| ALL | `S` | 0.592 |
| ALL | `gamma2` | 0.090 |

## Decision logic, applied

| stratum | survives | co-moves with A | sign flip | no effect | unstable |
|---|---|---|---|---|---|
| ALL | 0 | 1 (F) | 0 | 2 (P, S) | 2 (gamma2, D) |

The verdict is about the **covariates**, not about multiplicity: a cell can survive the area/temperature/season control and still miss its Holm threshold (see the last column of the minimal table). Read the two together.

## What this does not say

* `A` is not a neutral confounder here. The area is part of what the collapse changes, so the A-adjusted column measures **co-movement with the observed area**, not a correction of an acquisition artefact. Where the two ladders disagree in sign, the A-adjusted estimate is the one to distrust: conditioning on a variable that the event moves can create association as easily as it removes it.
* `T24` and `T72` do not exist in this record (one temperature per overpass). `T0` is the primary temperature and `T_smooth30` the sensitivity; neither is a measurement of the deck temperature.
* `Post` is a step at a single date. These numbers describe this one collapse; they are not a causal effect of collapsing and not a forecast.
* The placebos live inside the pre period and cannot reproduce the real group sizes; the size-matched rank is reported beside the raw one.
* Descending geometry has 56 % of its post area below the pre 5th percentile: its A-adjusted estimate is an **extrapolation** along the area relation, and the common-support subset is the number to read.

