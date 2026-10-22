# Carolabrücke deck — decomposition of the post-collapse signature

* table `data/carola/support/carola_ocv_channels.csv` (sha256 `3137d0b31fef80b5…`)
* event **2024-09-11**
* strata: ASC 36.9, ASC 44.8, DESC
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
| ASC 36.9 | `gamma2` | 1.714 | 1.340 | 1.716 | 1.315 | [0.816, 1.847] | 0.972 |
| ASC 36.9 | `P` | 1.029 | 1.044 | 1.064 | 1.053 | [0.995, 1.138] | 0.350 |
| ASC 36.9 | `D` | 0.829 | 1.017 | 0.833 | 1.016 | [0.988, 1.053] | 0.888 |
| ASC 36.9 | `F` | 0.850 | 1.009 | 0.851 | 1.011 | [0.985, 1.043] | 0.972 |
| ASC 36.9 | `S` | +0.67 | +0.111 | +0.917 | -0.030 | [-1.820, +1.146] | 0.972 |
| ASC 44.8 | `gamma2` | 0.534 | 0.872 | 0.519 | 0.950 | [0.649, 1.375] | 1.000 |
| ASC 44.8 | `P` | 0.797 | 0.972 | 0.802 | 0.970 | [0.913, 1.040] | 1.000 |
| ASC 44.8 | `D` | 1.495 | 1.004 | 1.562 | 1.005 | [0.971, 1.042] | 1.000 |
| ASC 44.8 | `F` | 1.488 | 0.962 | 1.569 | 0.965 | [0.927, 0.992] | 0.110 |
| ASC 44.8 | `S` | -3.75 | -0.830 | -3.997 | -0.724 | [-2.686, +1.278] | 1.000 |
| DESC | `gamma2` | 4.824 | 1.329 | 4.980 | 1.322 | [0.890, 1.766] | 0.432 |
| DESC | `P` | 0.808 | 0.815 | 0.642 | 0.797 | [0.678, 0.956] | 0.100 |
| DESC | `D` | 0.316 | 0.987 | 0.443 | 0.978 | [0.880, 1.067] | 0.644 |
| DESC | `F` | 0.311 | 1.050 | 0.367 | 1.048 | [0.994, 1.098] | 0.432 |
| DESC | `S` | +9.29 | -0.989 | +9.014 | -1.194 | [-3.070, +0.989] | 0.644 |

Ratio columns are `exp(beta_post)` for `gamma2, P, D, F` (adjusted post level as a fraction of pre); `S` is in pixels and shows the signed shift. `+T/season` is the richest A-adjusted model this site could fit (see the tier and the estimability rule above); the Holm column is family A: this site's targets **within** a stratum.

## ASC 36.9

*n* = 62 pre / 85 post. Missing values (never imputed):  — none

### Common support (pre 5th-95th percentile)

| variable | pre p5 | pre p95 | post range | post outside | overlap |
|---|---|---|---|---|---|
| `logA` | 4.83 | 6.50 | [3.93, 6.85] | **22.4 %** | 0.71 |
| `T` | -0.82 | 27.30 | [-3.70, 30.90] | **11.8 %** | 0.74 |
| `T_smooth30` | 1.53 | 23.38 | [-2.11, 22.00] | **10.6 %** | 0.69 |

Common-support subset (the pre p5-p95 box of `log A` **and** `T`): 108 of 147 rows kept, 62 of 85 post rows (73 %).

### Raw, unadjusted (reference only)

| target | n pre/post | median shift (g) | ratio | Cliff's delta [CI] | exp(beta) M0 [CI] |
|---|---|---|---|---|---|
| `gamma2` | 62/85 | +0.539 | 1.714 | +0.23 [+0.03, +0.44] | 1.741 [1.086, 2.798] |
| `P` | 62/85 | +0.029 | 1.029 | +0.15 [-0.07, +0.39] | 1.059 [0.975, 1.177] |
| `D` | 62/85 | -0.188 | 0.829 | -0.16 [-0.45, +0.09] | 0.820 [0.612, 1.049] |
| `F` | 62/85 | -0.162 | 0.850 | -0.15 [-0.44, +0.08] | 0.841 [0.664, 1.020] |
| `S` | 62/85 | +0.675 | — | +0.12 [-0.14, +0.36] | +1.021 [-1.034, +2.996] |

### A adjusted ladder

| target | model | post effect | 95% CI | p (bootstrap) | n pre/post | n |
|---|---|---|---|---|---|---|
| `gamma2` | M0 | 1.741 | [1.086, 2.798] | 0.012 | 62/85 | 147 |
| `gamma2` | M1 | 1.340 | [0.831, 1.914] | 0.288 | 62/85 | 147 |
| `gamma2` | M2 | 1.326 | [0.814, 1.831] | 0.306 | 62/85 | 147 |
| `gamma2` | M3 | 1.315 | [0.816, 1.847] | 0.324 | 62/85 | 147 |
| `P` | M0 | 1.059 | [0.975, 1.177] | 0.156 | 62/85 | 147 |
| `P` | M1 | 1.044 | [0.962, 1.154] | 0.272 | 62/85 | 147 |
| `P` | M2 | 1.048 | [0.987, 1.140] | 0.124 | 62/85 | 147 |
| `P` | M3 | 1.053 | [0.995, 1.138] | 0.070 | 62/85 | 147 |
| `D` | M0 | 0.820 | [0.612, 1.049] | 0.098 | 62/85 | 147 |
| `D` | M1 | 1.017 | [0.999, 1.051] | 0.064 | 62/85 | 147 |
| `D` | M2 | 1.016 | [0.995, 1.051] | 0.142 | 62/85 | 147 |
| `D` | M3 | 1.016 | [0.988, 1.053] | 0.222 | 62/85 | 147 |
| `F` | M0 | 0.841 | [0.664, 1.020] | 0.072 | 62/85 | 147 |
| `F` | M1 | 1.009 | [0.984, 1.039] | 0.498 | 62/85 | 147 |
| `F` | M2 | 1.010 | [0.987, 1.041] | 0.384 | 62/85 | 147 |
| `F` | M3 | 1.011 | [0.985, 1.043] | 0.402 | 62/85 | 147 |
| `S` | M0 | +1.021 | [-1.034, +2.996] | 0.306 | 62/85 | 147 |
| `S` | M1 | +0.111 | [-1.985, +1.845] | 0.986 | 62/85 | 147 |
| `S` | M2 | +0.033 | [-1.753, +1.306] | 0.748 | 62/85 | 147 |
| `S` | M3 | -0.030 | [-1.820, +1.146] | 0.722 | 62/85 | 147 |

### A free ladder

| target | model | post effect | 95% CI | p (bootstrap) | n pre/post | n |
|---|---|---|---|---|---|---|
| `gamma2` | M0 | 1.741 | [1.086, 2.798] | 0.012 | 62/85 | 147 |
| `gamma2` | M2 | 1.745 | [1.081, 2.767] | 0.014 | 62/85 | 147 |
| `gamma2` | M3 | 1.716 | [1.043, 2.766] | 0.030 | 62/85 | 147 |
| `P` | M0 | 1.059 | [0.975, 1.177] | 0.156 | 62/85 | 147 |
| `P` | M2 | 1.056 | [0.995, 1.152] | 0.076 | 62/85 | 147 |
| `P` | M3 | 1.064 | [1.003, 1.149] | 0.042 | 62/85 | 147 |
| `D` | M0 | 0.820 | [0.612, 1.049] | 0.098 | 62/85 | 147 |
| `D` | M2 | 0.823 | [0.615, 1.035] | 0.100 | 62/85 | 147 |
| `D` | M3 | 0.833 | [0.608, 1.052] | 0.118 | 62/85 | 147 |
| `F` | M0 | 0.841 | [0.664, 1.020] | 0.072 | 62/85 | 147 |
| `F` | M2 | 0.843 | [0.662, 1.013] | 0.066 | 62/85 | 147 |
| `F` | M3 | 0.851 | [0.655, 1.028] | 0.096 | 62/85 | 147 |
| `S` | M0 | +1.021 | [-1.034, +2.996] | 0.306 | 62/85 | 147 |
| `S` | M2 | +1.054 | [-0.845, +2.839] | 0.274 | 62/85 | 147 |
| `S` | M3 | +0.917 | [-0.852, +2.631] | 0.326 | 62/85 | 147 |

### Attenuation (descriptive, no threshold)

| target | raw -> +A | raw -> +A+T+season | raw -> A-free final | verdict |
|---|---|---|---|---|
| `gamma2` | 0.47 | 0.51 | 0.03 | co-moves with the observed area (A-free survives, A-adjusted does not) |
| `P` | 0.24 | 0.09 | -0.09 | co-moves with the observed area (A-free survives, A-adjusted does not) |
| `D` | 0.92 | 0.92 | 0.08 | not distinguishable from no effect after the covariates |
| `F` | 0.95 | 0.94 | 0.07 | not distinguishable from no effect after the covariates |
| `S` | 0.89 | 0.97 | 0.10 | not distinguishable from no effect after the covariates |

### Nested tests (partial F)

| target | +temperature | +season (A-adjusted) | +season (A-free) |
|---|---|---|---|
| `gamma2` | F 5.32 (df 1), p 0.023 | F 0.82 (df 2), p 0.442 | F 0.51 (df 2), p 0.601 |
| `P` | F 13.18 (df 1), p <0.001 | F 4.35 (df 2), p 0.015 | F 3.86 (df 2), p 0.023 |
| `D` | F 2.91 (df 1), p 0.090 | F 1.47 (df 2), p 0.235 | F 2.87 (df 2), p 0.060 |
| `F` | F 7.27 (df 1), p 0.008 | F 0.88 (df 2), p 0.415 | F 1.91 (df 2), p 0.152 |
| `S` | F 25.96 (df 1), p <0.001 | F 1.96 (df 2), p 0.145 | F 3.46 (df 2), p 0.034 |

### Block-length stability of the final post effect

| target | block 3 | block 5 | block 7 | same verdict |
|---|---|---|---|---|
| `gamma2` | 1.315 [0.833, 1.881] | 1.315 [0.816, 1.847] | 1.315 [0.785, 1.866] | yes |
| `P` | 1.053 [0.988, 1.132] | 1.053 [0.995, 1.138] | 1.053 [0.996, 1.148] | yes |
| `D` | 1.016 [0.993, 1.056] | 1.016 [0.988, 1.053] | 1.016 [0.980, 1.056] | yes |
| `F` | 1.011 [0.984, 1.041] | 1.011 [0.985, 1.043] | 1.011 [0.986, 1.046] | yes |
| `S` | -0.030 [-1.440, +1.145] | -0.030 [-1.820, +1.146] | -0.030 [-2.074, +1.375] | yes |

### Residual check (model fitted on pre only)

| target | pre median resid | post median resid | delta [CI] |
|---|---|---|---|
| `gamma2` | +0.169 | +0.341 | +0.10 [-0.16, +0.29] |
| `P` | +0.059 | +0.071 | +0.09 [-0.16, +0.36] |
| `D` | -0.001 | +0.000 | +0.13 [-0.08, +0.36] |
| `F` | +0.005 | +0.001 | +0.03 [-0.19, +0.28] |
| `S` | -0.356 | +0.612 | +0.07 [-0.23, +0.33] |

### Residual diagnostics (final model)

| target | lag-1 | runs p | acf 1..4 |
|---|---|---|---|
| `gamma2` | +0.06 | 0.506 | +0.06, +0.11, -0.12, +0.01 |
| `P` | -0.01 | 0.618 | -0.01, +0.01, -0.11, +0.02 |
| `D` | -0.03 | 0.008 | -0.03, +0.03, -0.03, +0.01 |
| `F` | +0.20 | 0.406 | +0.20, +0.18, -0.02, -0.11 |
| `S` | +0.03 | 0.245 | +0.03, +0.12, +0.10, +0.19 |

### Same calendar month, post versus pre

| target | months used | median ratio (log targets) / shift (S) |
|---|---|---|
| `gamma2` | 12 | 1.243 |
| `P` | 12 | 1.116 |
| `D` | 12 | 0.929 |
| `F` | 12 | 0.941 |
| `S` | 12 | +0.36 |

## ASC 44.8

*n* = 62 pre / 81 post. Missing values (never imputed):  — none

### Common support (pre 5th-95th percentile)

| variable | pre p5 | pre p95 | post range | post outside | overlap |
|---|---|---|---|---|---|
| `logA` | 3.69 | 5.60 | [4.13, 6.33] | **24.7 %** | 0.59 |
| `T` | 0.02 | 27.29 | [-3.40, 27.80] | **8.6 %** | 0.76 |
| `T_smooth30` | 1.85 | 23.44 | [-2.52, 22.19] | **11.1 %** | 0.73 |

Common-support subset (the pre p5-p95 box of `log A` **and** `T`): 104 of 143 rows kept, 57 of 81 post rows (70 %).

### Raw, unadjusted (reference only)

| target | n pre/post | median shift (g) | ratio | Cliff's delta [CI] | exp(beta) M0 [CI] |
|---|---|---|---|---|---|
| `gamma2` | 62/81 | -0.627 | 0.534 | -0.30 [-0.49, -0.12] | 0.483 [0.301, 0.719] |
| `P` | 62/81 | -0.227 | 0.797 | -0.44 [-0.64, -0.27] | 0.813 [0.728, 0.886] |
| `D` | 62/81 | +0.402 | 1.495 | +0.48 [+0.29, +0.70] | 1.538 [1.285, 1.933] |
| `F` | 62/81 | +0.398 | 1.488 | +0.43 [+0.20, +0.66] | 1.553 [1.239, 2.075] |
| `S` | 62/81 | -3.749 | — | -0.39 [-0.62, -0.15] | -4.039 [-6.422, -1.733] |

### A adjusted ladder

| target | model | post effect | 95% CI | p (bootstrap) | n pre/post | n |
|---|---|---|---|---|---|---|
| `gamma2` | M0 | 0.483 | [0.301, 0.719] | 0.004 | 62/81 | 143 |
| `gamma2` | M1 | 0.872 | [0.587, 1.259] | 0.436 | 62/81 | 143 |
| `gamma2` | M2 | 0.868 | [0.600, 1.211] | 0.396 | 62/81 | 143 |
| `gamma2` | M3 | 0.950 | [0.649, 1.375] | 0.718 | 62/81 | 143 |
| `P` | M0 | 0.813 | [0.728, 0.886] | <0.001 | 62/81 | 143 |
| `P` | M1 | 0.972 | [0.908, 1.048] | 0.372 | 62/81 | 143 |
| `P` | M2 | 0.973 | [0.912, 1.047] | 0.394 | 62/81 | 143 |
| `P` | M3 | 0.970 | [0.913, 1.040] | 0.376 | 62/81 | 143 |
| `D` | M0 | 1.538 | [1.285, 1.933] | <0.001 | 62/81 | 143 |
| `D` | M1 | 1.004 | [0.979, 1.036] | 0.630 | 62/81 | 143 |
| `D` | M2 | 1.004 | [0.975, 1.036] | 0.678 | 62/81 | 143 |
| `D` | M3 | 1.005 | [0.971, 1.042] | 0.596 | 62/81 | 143 |
| `F` | M0 | 1.553 | [1.239, 2.075] | <0.001 | 62/81 | 143 |
| `F` | M1 | 0.962 | [0.923, 0.992] | 0.014 | 62/81 | 143 |
| `F` | M2 | 0.963 | [0.926, 0.990] | 0.012 | 62/81 | 143 |
| `F` | M3 | 0.965 | [0.927, 0.992] | 0.022 | 62/81 | 143 |
| `S` | M0 | -4.039 | [-6.422, -1.733] | 0.002 | 62/81 | 143 |
| `S` | M1 | -0.830 | [-2.548, +0.975] | 0.392 | 62/81 | 143 |
| `S` | M2 | -0.832 | [-2.576, +0.943] | 0.392 | 62/81 | 143 |
| `S` | M3 | -0.724 | [-2.686, +1.278] | 0.588 | 62/81 | 143 |

### A free ladder

| target | model | post effect | 95% CI | p (bootstrap) | n pre/post | n |
|---|---|---|---|---|---|---|
| `gamma2` | M0 | 0.483 | [0.301, 0.719] | 0.004 | 62/81 | 143 |
| `gamma2` | M2 | 0.515 | [0.342, 0.726] | <0.001 | 62/81 | 143 |
| `gamma2` | M3 | 0.519 | [0.335, 0.746] | <0.001 | 62/81 | 143 |
| `P` | M0 | 0.813 | [0.728, 0.886] | <0.001 | 62/81 | 143 |
| `P` | M2 | 0.814 | [0.725, 0.886] | <0.001 | 62/81 | 143 |
| `P` | M3 | 0.802 | [0.717, 0.876] | <0.001 | 62/81 | 143 |
| `D` | M0 | 1.538 | [1.285, 1.933] | <0.001 | 62/81 | 143 |
| `D` | M2 | 1.509 | [1.285, 1.859] | <0.001 | 62/81 | 143 |
| `D` | M3 | 1.562 | [1.345, 1.928] | <0.001 | 62/81 | 143 |
| `F` | M0 | 1.553 | [1.239, 2.075] | <0.001 | 62/81 | 143 |
| `F` | M2 | 1.512 | [1.240, 1.949] | <0.001 | 62/81 | 143 |
| `F` | M3 | 1.569 | [1.287, 2.017] | <0.001 | 62/81 | 143 |
| `S` | M0 | -4.039 | [-6.422, -1.733] | 0.002 | 62/81 | 143 |
| `S` | M2 | -3.873 | [-6.348, -1.650] | 0.002 | 62/81 | 143 |
| `S` | M3 | -3.997 | [-6.464, -1.700] | 0.002 | 62/81 | 143 |

### Attenuation (descriptive, no threshold)

| target | raw -> +A | raw -> +A+T+season | raw -> A-free final | verdict |
|---|---|---|---|---|
| `gamma2` | 0.81 | 0.93 | 0.10 | co-moves with the observed area (A-free survives, A-adjusted does not) |
| `P` | 0.86 | 0.85 | -0.06 | co-moves with the observed area (A-free survives, A-adjusted does not) |
| `D` | 0.99 | 0.99 | -0.04 | co-moves with the observed area (A-free survives, A-adjusted does not) |
| `F` | 0.91 | 0.92 | -0.02 | sign flips between the ladders -- co-movement with the observed area, not a stable direction |
| `S` | 0.79 | 0.82 | 0.01 | co-moves with the observed area (A-free survives, A-adjusted does not) |

### Nested tests (partial F)

| target | +temperature | +season (A-adjusted) | +season (A-free) |
|---|---|---|---|
| `gamma2` | F 3.96 (df 1), p 0.049 | F 1.32 (df 2), p 0.271 | F 0.02 (df 2), p 0.981 |
| `P` | F 6.67 (df 1), p 0.011 | F 1.14 (df 2), p 0.323 | F 4.18 (df 2), p 0.017 |
| `D` | F 0.35 (df 1), p 0.555 | F 0.22 (df 2), p 0.802 | F 4.07 (df 2), p 0.019 |
| `F` | F 8.36 (df 1), p 0.004 | F 1.22 (df 2), p 0.297 | F 4.70 (df 2), p 0.011 |
| `S` | F 0.04 (df 1), p 0.847 | F 1.11 (df 2), p 0.332 | F 2.35 (df 2), p 0.099 |

### Block-length stability of the final post effect

| target | block 3 | block 5 | block 7 | same verdict |
|---|---|---|---|---|
| `gamma2` | 0.950 [0.634, 1.456] | 0.950 [0.649, 1.375] | 0.950 [0.654, 1.411] | yes |
| `P` | 0.970 [0.905, 1.044] | 0.970 [0.913, 1.040] | 0.970 [0.915, 1.050] | yes |
| `D` | 1.005 [0.972, 1.040] | 1.005 [0.971, 1.042] | 1.005 [0.970, 1.040] | yes |
| `F` | 0.965 [0.930, 0.996] | 0.965 [0.927, 0.992] | 0.965 [0.927, 0.991] | yes |
| `S` | -0.724 [-2.543, +1.282] | -0.724 [-2.686, +1.278] | -0.724 [-2.658, +1.436] | yes |

### Residual check (model fitted on pre only)

| target | pre median resid | post median resid | delta [CI] |
|---|---|---|---|
| `gamma2` | +0.183 | +0.044 | -0.11 [-0.28, +0.06] |
| `P` | +0.061 | -0.022 | -0.23 [-0.44, -0.01] |
| `D` | -0.007 | +0.032 | +0.27 [+0.10, +0.48] |
| `F` | +0.011 | -0.049 | -0.38 [-0.60, -0.19] |
| `S` | +0.371 | +0.906 | -0.01 [-0.17, +0.21] |

### Residual diagnostics (final model)

| target | lag-1 | runs p | acf 1..4 |
|---|---|---|---|
| `gamma2` | -0.09 | 0.866 | -0.09, +0.08, -0.13, -0.03 |
| `P` | +0.18 | 0.130 | +0.18, -0.03, -0.05, -0.02 |
| `D` | +0.05 | 0.312 | +0.05, +0.09, +0.06, -0.14 |
| `F` | +0.01 | 0.130 | +0.01, +0.05, +0.06, +0.02 |
| `S` | +0.04 | 0.613 | +0.04, -0.01, -0.05, -0.07 |

### Same calendar month, post versus pre

| target | months used | median ratio (log targets) / shift (S) |
|---|---|---|
| `gamma2` | 12 | 0.452 |
| `P` | 12 | 0.870 |
| `D` | 12 | 1.603 |
| `F` | 12 | 1.613 |
| `S` | 12 | -3.90 |

## DESC

*n* = 60 pre / 86 post. Missing values (never imputed):  — none

### Common support (pre 5th-95th percentile)

| variable | pre p5 | pre p95 | post range | post outside | overlap |
|---|---|---|---|---|---|
| `logA` | 4.51 | 6.58 | [2.08, 6.36] | **55.8 %** | 0.47 |
| `T` | -2.08 | 18.02 | [-12.60, 24.10] | **10.5 %** | 0.77 |
| `T_smooth30` | 1.32 | 21.80 | [-2.37, 22.27] | **11.6 %** | 0.66 |

Common-support subset (the pre p5-p95 box of `log A` **and** `T`): 83 of 146 rows kept, 35 of 86 post rows (41 %).

### Raw, unadjusted (reference only)

| target | n pre/post | median shift (g) | ratio | Cliff's delta [CI] | exp(beta) M0 [CI] |
|---|---|---|---|---|---|
| `gamma2` | 60/86 | +1.574 | 4.824 | +0.55 [+0.40, +0.67] | 4.972 [3.185, 7.809] |
| `P` | 60/86 | -0.214 | 0.808 | -0.29 [-0.47, -0.09] | 0.643 [0.546, 0.764] |
| `D` | 60/86 | -1.151 | 0.316 | -0.57 [-0.74, -0.41] | 0.426 [0.318, 0.554] |
| `F` | 60/86 | -1.167 | 0.311 | -0.58 [-0.75, -0.41] | 0.359 [0.262, 0.479] |
| `S` | 60/86 | +9.293 | — | +0.42 [+0.25, +0.63] | +9.232 [+6.096, +12.664] |

### A adjusted ladder

| target | model | post effect | 95% CI | p (bootstrap) | n pre/post | n |
|---|---|---|---|---|---|---|
| `gamma2` | M0 | 4.972 | [3.185, 7.809] | <0.001 | 60/86 | 146 |
| `gamma2` | M1 | 1.329 | [0.930, 1.748] | 0.108 | 60/86 | 146 |
| `gamma2` | M2 | 1.310 | [0.913, 1.700] | 0.156 | 60/86 | 146 |
| `gamma2` | M3 | 1.322 | [0.890, 1.766] | 0.144 | 60/86 | 146 |
| `P` | M0 | 0.643 | [0.546, 0.764] | <0.001 | 60/86 | 146 |
| `P` | M1 | 0.815 | [0.694, 0.952] | 0.010 | 60/86 | 146 |
| `P` | M2 | 0.803 | [0.683, 0.952] | 0.010 | 60/86 | 146 |
| `P` | M3 | 0.797 | [0.678, 0.956] | 0.020 | 60/86 | 146 |
| `D` | M0 | 0.426 | [0.318, 0.554] | <0.001 | 60/86 | 146 |
| `D` | M1 | 0.987 | [0.888, 1.058] | 0.650 | 60/86 | 146 |
| `D` | M2 | 0.969 | [0.878, 1.046] | 0.398 | 60/86 | 146 |
| `D` | M3 | 0.978 | [0.880, 1.067] | 0.620 | 60/86 | 146 |
| `F` | M0 | 0.359 | [0.262, 0.479] | <0.001 | 60/86 | 146 |
| `F` | M1 | 1.050 | [1.004, 1.095] | 0.034 | 60/86 | 146 |
| `F` | M2 | 1.054 | [1.005, 1.104] | 0.032 | 60/86 | 146 |
| `F` | M3 | 1.048 | [0.994, 1.098] | 0.108 | 60/86 | 146 |
| `S` | M0 | +9.232 | [+6.096, +12.664] | <0.001 | 60/86 | 146 |
| `S` | M1 | -0.989 | [-2.579, +0.795] | 0.336 | 60/86 | 146 |
| `S` | M2 | -1.273 | [-3.264, +0.731] | 0.244 | 60/86 | 146 |
| `S` | M3 | -1.194 | [-3.070, +0.989] | 0.322 | 60/86 | 146 |

### A free ladder

| target | model | post effect | 95% CI | p (bootstrap) | n pre/post | n |
|---|---|---|---|---|---|---|
| `gamma2` | M0 | 4.972 | [3.185, 7.809] | <0.001 | 60/86 | 146 |
| `gamma2` | M2 | 5.076 | [3.281, 7.638] | <0.001 | 60/86 | 146 |
| `gamma2` | M3 | 4.980 | [3.178, 7.713] | <0.001 | 60/86 | 146 |
| `P` | M0 | 0.643 | [0.546, 0.764] | <0.001 | 60/86 | 146 |
| `P` | M2 | 0.636 | [0.543, 0.759] | <0.001 | 60/86 | 146 |
| `P` | M3 | 0.642 | [0.543, 0.775] | <0.001 | 60/86 | 146 |
| `D` | M0 | 0.426 | [0.318, 0.554] | <0.001 | 60/86 | 146 |
| `D` | M2 | 0.416 | [0.325, 0.520] | <0.001 | 60/86 | 146 |
| `D` | M3 | 0.443 | [0.355, 0.550] | <0.001 | 60/86 | 146 |
| `F` | M0 | 0.359 | [0.262, 0.479] | <0.001 | 60/86 | 146 |
| `F` | M2 | 0.352 | [0.262, 0.462] | <0.001 | 60/86 | 146 |
| `F` | M3 | 0.367 | [0.272, 0.488] | <0.001 | 60/86 | 146 |
| `S` | M0 | +9.232 | [+6.096, +12.664] | <0.001 | 60/86 | 146 |
| `S` | M2 | +9.324 | [+6.219, +12.721] | <0.001 | 60/86 | 146 |
| `S` | M3 | +9.014 | [+6.014, +12.468] | <0.001 | 60/86 | 146 |

### Attenuation (descriptive, no threshold)

| target | raw -> +A | raw -> +A+T+season | raw -> A-free final | verdict |
|---|---|---|---|---|
| `gamma2` | 0.82 | 0.83 | -0.00 | co-moves with the observed area (A-free survives, A-adjusted does not) |
| `P` | 0.54 | 0.49 | -0.00 | survives area, temperature and season (both ladders exclude no effect, same sign) |
| `D` | 0.98 | 0.97 | 0.05 | co-moves with the observed area (A-free survives, A-adjusted does not) |
| `F` | 0.95 | 0.95 | 0.02 | co-moves with the observed area (A-free survives, A-adjusted does not) |
| `S` | 0.89 | 0.87 | 0.02 | co-moves with the observed area (A-free survives, A-adjusted does not) |

### Nested tests (partial F)

| target | +temperature | +season (A-adjusted) | +season (A-free) |
|---|---|---|---|
| `gamma2` | F 0.20 (df 1), p 0.652 | F 1.54 (df 2), p 0.218 | F 0.10 (df 2), p 0.905 |
| `P` | F 1.05 (df 1), p 0.306 | F 0.79 (df 2), p 0.454 | F 1.37 (df 2), p 0.256 |
| `D` | F 3.37 (df 1), p 0.069 | F 3.01 (df 2), p 0.053 | F 4.93 (df 2), p 0.009 |
| `F` | F 0.73 (df 1), p 0.394 | F 3.30 (df 2), p 0.040 | F 2.14 (df 2), p 0.122 |
| `S` | F 2.62 (df 1), p 0.108 | F 0.71 (df 2), p 0.494 | F 0.91 (df 2), p 0.407 |

### Block-length stability of the final post effect

| target | block 3 | block 5 | block 7 | same verdict |
|---|---|---|---|---|
| `gamma2` | 1.322 [0.892, 1.853] | 1.322 [0.890, 1.766] | 1.322 [0.903, 1.772] | yes |
| `P` | 0.797 [0.677, 0.949] | 0.797 [0.678, 0.956] | 0.797 [0.659, 0.947] | yes |
| `D` | 0.978 [0.870, 1.067] | 0.978 [0.880, 1.067] | 0.978 [0.869, 1.071] | yes |
| `F` | 1.048 [0.998, 1.104] | 1.048 [0.994, 1.098] | 1.048 [0.993, 1.103] | yes |
| `S` | -1.194 [-2.958, +0.736] | -1.194 [-3.070, +0.989] | -1.194 [-2.937, +1.121] | yes |

### Residual check (model fitted on pre only)

| target | pre median resid | post median resid | delta [CI] |
|---|---|---|---|
| `gamma2` | +0.140 | +0.637 | +0.24 [+0.10, +0.38] |
| `P` | +0.019 | -0.119 | -0.30 [-0.49, -0.09] |
| `D` | -0.001 | +0.019 | +0.41 [+0.23, +0.59] |
| `F` | +0.012 | -0.017 | -0.16 [-0.38, +0.07] |
| `S` | -0.870 | +2.021 | +0.19 [+0.00, +0.39] |

### Residual diagnostics (final model)

| target | lag-1 | runs p | acf 1..4 |
|---|---|---|---|
| `gamma2` | -0.14 | 0.319 | -0.14, -0.13, -0.13, +0.14 |
| `P` | -0.12 | 0.068 | -0.12, +0.19, +0.09, -0.20 |
| `D` | -0.12 | 0.868 | -0.12, +0.12, +0.05, +0.04 |
| `F` | -0.04 | 0.406 | -0.04, +0.16, -0.01, +0.03 |
| `S` | -0.17 | 0.245 | -0.17, +0.17, -0.10, -0.11 |

### Same calendar month, post versus pre

| target | months used | median ratio (log targets) / shift (S) |
|---|---|---|
| `gamma2` | 12 | 4.175 |
| `P` | 12 | 0.898 |
| `D` | 12 | 0.471 |
| `F` | 12 | 0.438 |
| `S` | 12 | +5.61 |

## Geometry — the joint ASC model (family B)

Both ASC tracks in one model: `post`, `angle` and `post x angle`, with `log A`, temperature and season. With two angle values the interaction is the stratum contrast, so this adds no information beyond the two stratified fits — it is reported as its own family, not as a headline.

| target | post | angle | post x angle | F (interaction) | p |
|---|---|---|---|---|---|
| `gamma2` | +0.114 [-0.203, +0.359] | -0.165 [-0.232, -0.103] | -0.046 [-0.108, +0.031] | 1.65 | 0.200 |
| `P` | -0.045 [-0.111, +0.017] | +0.042 [+0.030, +0.056] | -0.017 [-0.029, -0.008] | 6.29 | 0.013 |
| `D` | -0.015 [-0.048, +0.024] | +0.003 [-0.007, +0.010] | -0.005 [-0.010, +0.000] | 1.79 | 0.182 |
| `F` | +0.005 [-0.018, +0.026] | -0.019 [-0.024, -0.013] | -0.003 [-0.010, +0.002] | 1.40 | 0.237 |
| `S` | -0.675 [-1.647, +0.376] | +0.667 [+0.367, +0.904] | -0.160 [-0.453, +0.221] | 1.36 | 0.245 |

## The whole vector (secondary)

Standardised residuals of `M3-` (temperature and season, no area, no post) fitted on the pre data, with a circular-shift null that keeps the serial dependence and the vector's internal correlation.

| stratum | n | n pre/post | statistic | p (shift) | dimensions |
|---|---|---|---|---|---|
| ASC 36.9 | 147 | 62/85 | 0.614 | 0.316 | gamma2, P, D, F, S |
| ASC 44.8 | 143 | 62/81 | 3.675 | 0.091 | gamma2, P, D, F, S |
| DESC | 146 | 60/86 | 12.869 | 0.175 | gamma2, P, D, F, S |

## Placebo events (exploratory, family C)

Pseudo events inside the pre period, run through the same estimator. The pseudo-post group is smaller than the real one, so the size-matched rank is the fair comparison; both are shown.

| stratum | spec | real | rank (raw) | matched size | rank (matched) | of |
|---|---|---|---|---|---|---|
| ASC 36.9 | M0 | 0.472 | 9 | 37 | 9 | 8 |
| ASC 36.9 | M3_free | 0.436 | 9 | 37 | 9 | 8 |
| ASC 44.8 | M0 | 2.743 | 1 | 37 | 1 | 8 |
| ASC 44.8 | M3_free | 2.839 | 3 | 37 | 3 | 8 |
| DESC | M0 | 4.090 | 1 | 36 | 1 | 8 |
| DESC | M3_free | 3.934 | 1 | 36 | 1 | 8 |

A rank of 1 means the real event produced a stronger signature than every pseudo event; a rank of n+1 means all pseudo events were stronger.

## Holm across all target x stratum post effects (sensitivity)

Family size 15; values in `report['holm_all_strata']`.

| stratum | target | Holm p (all 15) |
|---|---|---|
| ASC 36.9 | `D` | 1.000 |
| ASC 36.9 | `F` | 1.000 |
| ASC 36.9 | `P` | 0.910 |
| ASC 36.9 | `S` | 1.000 |
| ASC 36.9 | `gamma2` | 1.000 |
| ASC 44.8 | `D` | 1.000 |
| ASC 44.8 | `F` | 0.308 |
| ASC 44.8 | `P` | 1.000 |
| ASC 44.8 | `S` | 1.000 |
| ASC 44.8 | `gamma2` | 1.000 |
| DESC | `D` | 1.000 |
| DESC | `F` | 1.000 |
| DESC | `P` | 0.300 |
| DESC | `S` | 1.000 |
| DESC | `gamma2` | 1.000 |

## Decision logic, applied

| stratum | survives | co-moves with A | sign flip | no effect | unstable |
|---|---|---|---|---|---|
| ASC 36.9 | 0 | 2 (gamma2, P) | 0 | 3 (D, F, S) | 0 |
| ASC 44.8 | 0 | 4 (gamma2, P, D, S) | 1 (F) | 0 | 0 |
| DESC | 1 (P) | 4 (gamma2, D, F, S) | 0 | 0 | 0 |

The verdict is about the **covariates**, not about multiplicity: a cell can survive the area/temperature/season control and still miss its Holm threshold (see the last column of the minimal table). Read the two together.

## What this does not say

* `A` is not a neutral confounder here. The area is part of what the collapse changes, so the A-adjusted column measures **co-movement with the observed area**, not a correction of an acquisition artefact. Where the two ladders disagree in sign, the A-adjusted estimate is the one to distrust: conditioning on a variable that the event moves can create association as easily as it removes it.
* `T24` and `T72` do not exist in this record (one temperature per overpass). `T0` is the primary temperature and `T_smooth30` the sensitivity; neither is a measurement of the deck temperature.
* `Post` is a step at a single date. These numbers describe this one collapse; they are not a causal effect of collapsing and not a forecast.
* The placebos live inside the pre period and cannot reproduce the real group sizes; the size-matched rank is reported beside the raw one.
* Descending geometry has 56 % of its post area below the pre 5th percentile: its A-adjusted estimate is an **extrapolation** along the area relation, and the common-support subset is the number to read.

