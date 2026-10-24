# LUMO lattice tower — correspondence between the physical and the OCV plane

**Unit of inference: the campaign / visit, never the single overpass.** Several SAR acquisitions inside one structural campaign improve the estimate of (gamma2, P); they do not create additional independent measurements of (f, zeta). Every *inference* interval (`ci_inference`, the permutation p-values, the Holm adjustment) therefore resamples campaigns; the date-level bootstrap appears only as `ci_measurement` inside a campaign and never increases n.

| property | value |
| --- | --- |
| asset | LUMO lattice tower |
| unit | campaign (visit) |
| reference | healthy |
| channel transform | log |
| x | H1 mode: signed df_H1/f_H1,0 and dzeta_H1/zeta_H1,0 against the healthy block of the same campaign |
| y | z of the campaign median log(gamma2) and log(P) against the healthy median/MAD of the same orbit |
| data | `data/lumo/support/lumo_ocv_channels.csv` sha256 `4874a2c8c9…` |
| physics | `data/lumo/reference/lumo_state_dynamics.json` sha256 `bdcd8fb071…` |

## 1. Campaigns

| campaign | state | visit | x (df/f0, dz/z0) | y ASC | y DESC | n dates ASC/DESC | reference |
| --- | --- | --- | --- | --- | --- | --- | --- |
| DAM6 v1 | DAM 6 | 1 | (-0.1296, -0.2802) | (+1.33, +0.95) | (-1.08, +0.05) | 2/6 |  |
| DAM4 v1 | DAM 4 | 1 | (-0.0741, -0.2522) | (-1.79, +1.16) | (+1.03, +0.76) | 3/5 |  |
| DAM3 v1 | DAM 3 | 1 | (-0.0741, -0.2368) | (-0.46, -0.51) | (-1.88, -1.03) | 5/12 |  |
| DAM6 v2 | DAM 6 | 2 | (-0.1296, -0.2802) | (+0.91, +0.10) | (-1.22, -0.28) | 2/6 |  |
| DAM4 v2 | DAM 4 | 2 | (-0.0741, -0.2522) | (+0.35, -0.58) | (+0.10, +0.09) | 2/6 |  |
| DAM3 v2 | DAM 3 | 2 | (-0.0741, -0.2368) | (+0.34, -1.54) | (+0.10, -0.92) | 3/5 |  |

The healthy reference is not a campaign of its own: it is the healthy block of each campaign.

### OCV vectors with measurement intervals

`ci_measurement` resamples the dates **inside** a campaign; it is the uncertainty of that campaign's median, not an increase in the number of campaigns.

| campaign | orbit | y | ci gamma2 | ci P | n dates | x |
| --- | --- | --- | --- | --- | --- | --- |
| DAM6 v1 | ASC | (+1.33, +0.95) | [0.81; 1.85] | [0.01; 1.89] | 2 | (-0.1296, -0.2802) |
| DAM6 v1 | DESC | (-1.08, +0.05) | [-3.87; 1.19] | [-1.74; 1.77] | 6 | (-0.1296, -0.2802) |
| DAM4 v1 | ASC | (-1.79, +1.16) | [-4.38; 0.66] | [-0.39; 1.67] | 3 | (-0.0741, -0.2522) |
| DAM4 v1 | DESC | (+1.03, +0.76) | [-2.01; 1.18] | [-1.00; 1.98] | 5 | (-0.0741, -0.2522) |
| DAM3 v1 | ASC | (-0.46, -0.51) | [-5.68; 1.14] | [-1.07; -0.01] | 5 | (-0.0741, -0.2368) |
| DAM3 v1 | DESC | (-1.88, -1.03) | [-5.19; 0.48] | [-2.17; 0.22] | 12 | (-0.0741, -0.2368) |
| DAM6 v2 | ASC | (+0.91, +0.10) | [-0.34; 2.16] | [-0.09; 0.29] | 2 | (-0.1296, -0.2802) |
| DAM6 v2 | DESC | (-1.22, -0.28) | [-2.29; 0.03] | [-1.55; 0.99] | 6 | (-0.1296, -0.2802) |
| DAM4 v2 | ASC | (+0.35, -0.58) | [-0.77; 1.47] | [-0.75; -0.42] | 2 | (-0.0741, -0.2522) |
| DAM4 v2 | DESC | (+0.10, +0.09) | [-1.49; 1.09] | [-1.87; 1.08] | 6 | (-0.0741, -0.2522) |
| DAM3 v2 | ASC | (+0.34, -1.54) | [-0.33; 1.61] | [-2.67; -1.44] | 3 | (-0.0741, -0.2368) |
| DAM3 v2 | DESC | (+0.10, -0.92) | [-7.21; 1.21] | [-0.97; 0.56] | 5 | (-0.0741, -0.2368) |

## 2. A — is the OCV response repeatable?

| orbit | state | n dates | d_visit | cos | angle [deg] | ci_measurement | d_visit / d_between |
| --- | --- | --- | --- | --- | --- | --- | --- |
| ASC | DAM 3 | 5/3 | 1.31 | 0.577 | 54.7 | [0.94; 6.07] | 0.71 |
| ASC | DAM 4 | 3/2 | 2.76 | -0.897 | 153.8 | [0.37; 6.33] | 1.49 |
| ASC | DAM 6 | 2/2 | 0.95 | 0.871 | 29.4 | [0.13; 2.96] | 0.51 |
| DESC | DAM 3 | 12/5 | 1.98 | 0.387 | 67.2 | [0.25; 6.56] | 1.10 |
| DESC | DAM 4 | 5/6 | 1.15 | 0.993 | 6.8 | [0.22; 3.45] | 0.64 |
| DESC | DAM 6 | 6/6 | 0.36 | 0.963 | 15.7 | [0.34; 3.56] | 0.20 |

Median ratio **R = 0.672** against the pre-declared rule `R < 1`: the mapping is **identifiable in 4 of 6 state-orbit cells**.

Visits with a cosine below 0.9: DAM 3 ASC, DAM 4 ASC, DAM 6 ASC, DAM 3 DESC.

## 3. B — does the acquisition geometry change the mapping?

| campaign | state | dy = y_ASC - y_DESC | \|dy\| | angle [deg] | ci_measurement \|dy\| | n dates |
| --- | --- | --- | --- | --- | --- | --- |
| DAM6 v1 | DAM 6 | (+2.41, +0.90) | 2.58 | 20.5 | [0.37; 6.10] | 2/6 |
| DAM4 v1 | DAM 4 | (-2.82, +0.40) | 2.85 | 172.0 | [0.49; 5.63] | 3/5 |
| DAM3 v1 | DAM 3 | (+1.42, +0.52) | 1.51 | 20.2 | [0.34; 5.84] | 5/12 |
| DAM6 v2 | DAM 6 | (+2.13, +0.38) | 2.17 | 10.1 | [0.56; 3.95] | 2/6 |
| DAM4 v2 | DAM 4 | (+0.25, -0.68) | 0.72 | -69.8 | [0.21; 2.40] | 2/6 |
| DAM3 v2 | DAM 3 | (+0.25, -0.62) | 0.67 | -68.3 | [0.57; 7.57] | 3/5 |

Between-state share of the `dy` variance: **0.666**, exact permutation over 720 relabellings, p = **0.133** (null median 0.369).

A large share with a small p means the ASC-DESC contrast changes with the state, i.e. the physical response interacts with the geometry. With two visits per state it cannot be separated from the visit scatter itself, which is why the per-contrast intervals are printed next to it.

## 4. C — do physical and OCV distances correspond?

| orbit | campaigns | rho (raw) | p (Mantel) | ci_inference (raw) | rho (standardised) | p (standardised) |
| --- | --- | --- | --- | --- | --- | --- |
| ASC | 6 | 0.184 | 0.244 | [0.02; 1.00] | 0.184 | 0.244 |
| DESC | 6 | -0.015 | 0.533 | [-0.18; 1.00] | -0.015 | 0.533 |

`d_phys` is the Euclidean distance in (df/f0, dz/z0): the physical meaning is kept, so the damping axis dominates by its own scale. The standardised column is the **sensitivity** that asks whether the conclusion depends on that scale; it is never used to make the planes look alike. The Mantel test permutes which physical vector belongs to which campaign, exhaustively where the campaign count allows.

Here the raw and the standardised rho are **identical**: with three distinct physical vectors the rank order of the pair distances does not change when the axes are rescaled, so the conclusion does not depend on axis scaling.

All campaign pairs (ASC, raw scale):

| pair | d_phys | d_OCV |
| --- | --- | --- |
| DAM6 v1 - DAM4 v1 | 0.0622 | 3.12 |
| DAM6 v1 - DAM3 v1 | 0.0705 | 2.31 |
| DAM6 v1 - DAM6 v2 | 0.0000 | 0.95 |
| DAM6 v1 - DAM4 v2 | 0.0622 | 1.82 |
| DAM6 v1 - DAM3 v2 | 0.0705 | 2.69 |
| DAM4 v1 - DAM3 v1 | 0.0154 | 2.13 |
| DAM4 v1 - DAM6 v2 | 0.0622 | 2.90 |
| DAM4 v1 - DAM4 v2 | 0.0000 | 2.76 |
| DAM4 v1 - DAM3 v2 | 0.0154 | 3.44 |
| DAM3 v1 - DAM6 v2 | 0.0705 | 1.50 |
| DAM3 v1 - DAM4 v2 | 0.0154 | 0.81 |
| DAM3 v1 - DAM3 v2 | 0.0000 | 1.31 |
| DAM6 v2 - DAM4 v2 | 0.0622 | 0.88 |
| DAM6 v2 - DAM3 v2 | 0.0705 | 1.74 |
| DAM4 v2 - DAM3 v2 | 0.0154 | 0.96 |

## 5. D — the closest physical pair

Contrast pair: **DAM 3 vs DAM 4** (physical distance 0.0154). `x` takes only **3 distinct values** while the model has **3 parameters** per channel, so `y = a + B x` is saturated at the state level and the residual contrast would be zero for any data. It is reported as **not identifiable**; the raw contrast is measured instead.

| orbit | raw contrast | \|contrast\| | ci_measurement | median visit noise | ratio |
| --- | --- | --- | --- | --- | --- |
| ASC | (+0.66, -1.31) | 1.47 | [0.43; 3.72] | 1.31 | 1.12 |
| DESC | (-1.46, -1.40) | 2.02 | [0.32; 5.33] | 1.15 | 1.76 |

A ratio near or below one means the OCV contrast between the two physically near-identical states is no larger than the scatter between two visits of the *same* state.

## 6. Nested models

| model | formula | parameters / channel | status |
| --- | --- | --- | --- |
| M0 | y ~ orbit | 2 | fitted |
| M1 | y ~ orbit + x | 4 | fitted |
| M2 | y ~ M1 + orbit x x | 6 | **not identifiable** |

M2: six independent campaigns and six scalar parameters per channel: the two orbits of a campaign share x, so the twelve campaign-orbit cells are not twelve independent units and M2 has no residual degrees of freedom. It is reported, not fitted -- a ridge fit would produce numbers without producing information.

| orbit | LOCO RMSE M0 | LOCO RMSE M1 | M0 - M1 | p (x permuted) | folds |
| --- | --- | --- | --- | --- | --- |
| ASC | 1.168 | 1.305 | -0.137 | 0.078 | 6 folds |
| DESC | 0.987 | 0.944 | 0.042 | 0.178 | 6 folds |

`M0 - M1 > 0` means the physical vector improves out-of-campaign prediction; the p-value permutes the assignment of physical vectors to campaigns. the table above is leave-one-campaign-out only: in-sample R2 would be meaningless because x takes one value per state, so M1 reproduces the three state means exactly and only the within-state variation is left as residual.

## 7. Multiplicity

| test | p | p (Holm) |
| --- | --- | --- |
| B.geometry | 0.133 | 0.533 |
| C.rho.ASC | 0.244 | 0.533 |
| C.rho.DESC | 0.533 | 0.533 |
| M1.loco.ASC | 0.078 | 0.389 |
| M1.loco.DESC | 0.178 | 0.533 |

Holm inside this five-test family; the standardised distance variant is a sensitivity and stays outside it.

## 8. Cross-checks

**17 checks, 0 failures.**

| group | checks |
| --- | --- |
| lumo | 17 |

## 9. What this does and does not show

- **No p-value in A or D.** Three repeat pairs and one contrast cannot support a test; both are reported as effect sizes with measurement intervals.
- **Few effective physical points.** The campaign pairs are not independent observations: the visits of a state share `x`, so the information about the mapping comes from the number of distinct configurations, not from the number of pairs.
- **A permutation p answers a narrow question**: is the observed assignment of physical vectors to campaigns unusual, not whether a mechanism exists. The permutation is exact for the campaign count used.
- **A matching direction would not show that (gamma2, P) measures stiffness or damping.** It would show that two measurement channels agree about which state is which.
- **Assets are not pooled.** LUMO and Bautzen are modelled separately; only the pattern (sign, order of magnitude) is comparable, never coefficients or vector lengths.

---

Generated by `code/cases/ocv_correspondence.py --asset lumo` (17 cross-checks, 0 failures).
