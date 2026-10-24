# Bautzen OpenLabs bridge — correspondence between the physical and the OCV plane

**Unit of inference: the campaign / visit, never the single overpass.** Several SAR acquisitions inside one structural campaign improve the estimate of (gamma2, P); they do not create additional independent measurements of (f, zeta). Every *inference* interval (`ci_inference`, the permutation p-values, the Holm adjustment) therefore resamples campaigns; the date-level bootstrap appears only as `ci_measurement` inside a campaign and never increases n.

| property | value |
| --- | --- |
| asset | Bautzen OpenLabs bridge |
| unit | campaign (visit) |
| reference | RS |
| channel transform | identity |
| x | primary matched mode pair: signed df/f0 and dzeta/zeta0 against the RS night; DS2 has no independent physics |
| y | z of the median paired difference (bridge - control) of gamma2 and P against the RS distribution of the same paired difference, same orbit; identity scale because the floor-corrected gamma2 is not strictly positive |
| data | `data/bautzen/bautzen_ocv_channels.csv` sha256 `32b2180cc2…` |
| physics | `data/bautzen/reference/bautzen_state_dynamics.json` sha256 `c1c61eb935…` |

`gamma2` uses the **identity** transform, not `log`: the floor-corrected Bautzen coherence takes non-positive values, so a logarithm is undefined. The statistic is identical to LUMO; only this declared per-asset transform differs.

## 1. Campaigns

| campaign | state | visit | x (df/f0, dz/z0) | y ASC | y DESC | n dates ASC/DESC | reference |
| --- | --- | --- | --- | --- | --- | --- | --- |
| RS | RS | 1 | (+0.0000, +0.0000) | (+0.00, +0.00) | (+0.00, +0.00) | 6/7 | yes |
| DS1 | DS1 | 1 | (-0.0100, -0.1234) | (+1.68, +0.48) | (-0.44, +2.58) | 10/23 |  |
| DS2 | DS2 | 1 | not available | (+2.04, -0.14) | (-0.41, +4.09) | 6/5 |  |

The reference campaign sits at the origin of both planes by construction and is not a subject of the tests.

**Pairing.** The OCV vector is the per-date difference bridge − control; a date enters only when **both channels exist in both series**. Pairing is exact (each bridge acquisition has a control acquisition on the same date), but the control chips are missing on some ASC dates, which is why the ASC counts are lower:

| cell | dates | complete pairs | incomplete |
| --- | --- | --- | --- |
| DS1.ASC | 23 | 10 | 13 |
| DS1.DESC | 23 | 23 | 0 |
| DS2.ASC | 6 | 6 | 0 |
| DS2.DESC | 5 | 5 | 0 |
| RS.ASC | 7 | 6 | 1 |
| RS.DESC | 7 | 7 | 0 |

### OCV vectors with measurement intervals

`ci_measurement` resamples the dates **inside** a campaign; it is the uncertainty of that campaign's median, not an increase in the number of campaigns.

| campaign | orbit | y | ci gamma2 | ci P | n dates | x |
| --- | --- | --- | --- | --- | --- | --- |
| RS | ASC | (+0.00, +0.00) | [-1.00; 2.36] | [-1.37; 1.08] | 6 | (+0.0000, +0.0000) |
| RS | DESC | (+0.00, +0.00) | [-1.00; 1.01] | [-1.00; 1.99] | 7 | (+0.0000, +0.0000) |
| DS1 | ASC | (+1.68, +0.48) | [-0.21; 2.27] | [-0.46; 1.20] | 10 | (-0.0100, -0.1234) |
| DS1 | DESC | (-0.44, +2.58) | [-1.18; 0.33] | [-0.20; 3.45] | 23 | (-0.0100, -0.1234) |
| DS2 | ASC | (+2.04, -0.14) | [-0.04; 2.32] | [-2.68; 0.54] | 6 | not available |
| DS2 | DESC | (-0.41, +4.09) | [-3.23; 5.30] | [-6.64; 8.12] | 5 | not available |

## 2. A — is the OCV response repeatable?

**Not defined.** no state was visited twice: every configuration has exactly one campaign, so a repeatability estimate would have to come from dates inside one campaign, which are not independent visits of the tower

## 3. B — does the acquisition geometry change the mapping?

**Not defined.** 2 non-reference campaigns carry both geometries; a between-state share needs at least three

## 4. C — do physical and OCV distances correspond?

**ASC: not defined.** fewer than three campaigns carry both a physical and an OCV vector, so there are no distance matrices to compare
**DESC: not defined.** fewer than three campaigns carry both a physical and an OCV vector, so there are no distance matrices to compare

## 5. D — the closest physical pair

Contrast pair: **DS1 vs RS** (physical distance 0.1238). `x` takes only **2 distinct values** while the model has **3 parameters** per channel, so `y = a + B x` is saturated at the state level and the residual contrast would be zero for any data. It is reported as **not identifiable**; the raw contrast is measured instead.

| orbit | raw contrast | \|contrast\| | ci_measurement | median visit noise | ratio |
| --- | --- | --- | --- | --- | --- |
| ASC | (+1.68, +0.48) | 1.75 | [0.35; 2.97] | -- | -- |
| DESC | (-0.44, +2.58) | 2.61 | [0.23; 4.71] | -- | -- |

This asset has no repeated visit, so there is no within-state yardstick for the contrast: its own date-level interval is the only uncertainty it carries, and no repeatability claim can be attached to it.

## 6. Nested models

| model | formula | parameters / channel | status |
| --- | --- | --- | --- |
| M0 | y ~ orbit | 2 | fitted |
| M1 | y ~ orbit + x | 4 | not identifiable |
| M2 | y ~ M1 + orbit x x | 6 | **not identifiable** |

M2: six independent campaigns and six scalar parameters per channel: the two orbits of a campaign share x, so the twelve campaign-orbit cells are not twelve independent units and M2 has no residual degrees of freedom. It is reported, not fitted -- a ridge fit would produce numbers without producing information.

No orbit carries enough campaigns for the M0/M1 comparison either.

## 7. Multiplicity

**No test in this family is defined for this asset**, so there is nothing to adjust; the report above states for each test why.

## 8. Cross-checks

**16 checks, 0 failures.**

| group | checks |
| --- | --- |
| bautzen | 16 |

## 9. What this does and does not show

- **No p-value in A or D.** Three repeat pairs and one contrast cannot support a test; both are reported as effect sizes with measurement intervals.
- **Few effective physical points.** The campaign pairs are not independent observations: the visits of a state share `x`, so the information about the mapping comes from the number of distinct configurations, not from the number of pairs.
- **A permutation p answers a narrow question**: is the observed assignment of physical vectors to campaigns unusual, not whether a mechanism exists. The permutation is exact for the campaign count used.
- **A matching direction would not show that (gamma2, P) measures stiffness or damping.** It would show that two measurement channels agree about which state is which.
- **Assets are not pooled.** LUMO and Bautzen are modelled separately; only the pattern (sign, order of magnitude) is comparable, never coefficients or vector lengths.

---

Generated by `code/cases/ocv_correspondence.py --asset bautzen` (16 cross-checks, 0 failures).
