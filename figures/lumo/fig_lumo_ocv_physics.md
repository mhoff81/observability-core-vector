# LUMO — physical dynamics against OCV behaviour

Two planes, same three states.  Panel A is the **independently measured physical** response of the tower, panel B the **OCV behaviour** of the same states.  The question is whether states that sit close together in A also sit close together in B.

```
v_phys = [ df_H1 / f_H1,0 ,  dzeta_H1 / zeta_H1,0 ]
v_OCV  = [ z(gamma2)      ,  z(P)               ]
```

`f_H1` is the first higher mode of the Welch mode consensus, `zeta_H1` its half-power bandwidth damping, both measured in the controlled campaign of the state (five healthy recordings, then five damaged).  `z` is the median log channel of the state against the healthy median of the **same orbit**, in healthy MADs of the log channel.  The origin of both panels is the corresponding healthy state.

Inputs: `data/lumo/support/lumo_ocv_channels.csv` (support chip, 178 overpasses, sha256 `4874a2c8c938…`) and `data/lumo/reference/lumo_state_dynamics.json` (107 reference checks, 0 failures).

## 1. The asymmetry that shapes the figure

| plane | objects | replicate |
| --- | --- | --- |
| A physical | 3 states | **none** — one campaign per state (5 + 5 recordings) |
| B OCV | 3 states x 2 orbits x 2 visits | 2 SAR visits (13 consecutive state runs in the record) |

So visit-to-visit repeatability can only be assessed on the OCV side.  Panel A carries a 95 % bootstrap interval over the ten recordings instead, and its arrows are drawn only for the unambiguous healthy → state pairing of one campaign.

## 2. Panel A — physical response

| state | healthy | damaged | df/f0 | dzeta/zeta0 | 95 % CI (damping) | \|v\| | angle [deg] | fund. dzeta |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| DAM 3 | 05_Healthy | 06_DAM3_111 | -7.41 % | -23.68 % | [-40.0 %; -5.4 %] | 0.248 | -107.4 | -2.3 % |
| DAM 4 | 03_Healthy | 04_DAM4_111 | -7.41 % | -25.22 % | [-44.1 %; -11.7 %] | 0.263 | -106.4 | +11.0 % |
| DAM 6 | 01_Healthy | 02_DAM6_111 | -12.96 % | -28.02 % | [-40.7 %; -15.1 %] | 0.309 | -114.8 | +9.6 % |

`angle` is `atan2(dzeta/zeta0, df/f0)` in the plotted plane; the two axes are both relative changes, so the direction is scale free.

## 3. Panel B — OCV behaviour per orbit (never pooled)

| orbit | state | n | z(gamma2) | z(P) | \|v\| | angle [deg] | n gamma2/P |
| --- | --- | --- | --- | --- | --- | --- | --- |
| ASC | DAM 3 | 8 | +0.01 | -0.94 | 0.94 | -89.6 | 8/8 |
| ASC | DAM 4 | 5 | -0.77 | -0.39 | 0.87 | -153.5 | 5/5 |
| ASC | DAM 6 | 4 | +1.33 | +0.15 | 1.34 | 6.5 | 4/4 |
| DESC | DAM 3 | 17 | -0.23 | -0.92 | 0.95 | -103.8 | 17/17 |
| DESC | DAM 4 | 11 | +0.28 | +0.57 | 0.64 | 64.2 | 11/11 |
| DESC | DAM 6 | 12 | -1.22 | -0.11 | 1.23 | -175.0 | 12/12 |

Healthy reference (median and MAD of the log channel):

| orbit | n gamma2 | median log gamma2 | MAD log gamma2 | n P | median log P | MAD log P |
| --- | --- | --- | --- | --- | --- | --- |
| ASC | 40 | -2.0872 | 0.7706 | 40 | -2.3695 | 0.4327 |
| DESC | 81 | -0.4294 | 0.2960 | 81 | -1.6520 | 0.3928 |

### Visit-level points

| orbit | state | visit | first | last | n | z(gamma2) | z(P) | n < 3 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| ASC | DAM 3 | 1 | 2021-03-18 | 2021-04-18 | 5 | -0.46 | -0.51 |  |
| ASC | DAM 3 | 2 | 2021-06-26 | 2021-07-11 | 3 | +0.34 | -1.54 |  |
| ASC | DAM 4 | 1 | 2020-11-10 | 2020-11-24 | 3 | -1.79 | +1.16 |  |
| ASC | DAM 4 | 2 | 2021-05-29 | 2021-06-14 | 2 | +0.35 | -0.58 | yes |
| ASC | DAM 6 | 1 | 2020-10-13 | 2020-10-26 | 2 | +1.33 | +0.95 | yes |
| ASC | DAM 6 | 2 | 2021-05-05 | 2021-05-18 | 2 | +0.91 | +0.10 | yes |
| DESC | DAM 3 | 1 | 2021-03-18 | 2021-04-18 | 12 | -1.88 | -1.03 |  |
| DESC | DAM 3 | 2 | 2021-06-26 | 2021-07-11 | 5 | +0.10 | -0.92 |  |
| DESC | DAM 4 | 1 | 2020-11-10 | 2020-11-24 | 5 | +1.03 | +0.76 |  |
| DESC | DAM 4 | 2 | 2021-05-29 | 2021-06-14 | 6 | +0.10 | +0.09 |  |
| DESC | DAM 6 | 1 | 2020-10-13 | 2020-10-26 | 6 | -1.08 | +0.05 |  |
| DESC | DAM 6 | 2 | 2021-05-05 | 2021-05-18 | 6 | -1.22 | -0.28 |  |

## 4. Visit-to-visit repeatability (OCV side only)

| orbit | state | n visit 1 | n visit 2 | cos(v1, v2) | angle [deg] | thin |
| --- | --- | --- | --- | --- | --- | --- |
| ASC | DAM 3 | 5 | 3 | +0.58 | 54.7 |  |
| ASC | DAM 4 | 3 | 2 | -0.90 | 153.8 | yes |
| ASC | DAM 6 | 2 | 2 | +0.87 | 29.4 | yes |
| DESC | DAM 3 | 12 | 5 | +0.39 | 67.2 |  |
| DESC | DAM 4 | 5 | 6 | +0.99 | 6.8 |  |
| DESC | DAM 6 | 6 | 6 | +0.96 | 15.7 |  |

Cosine above 0.9: DAM 4 DESC, DAM 6 DESC.  Cosine below 0.5: DAM 4 ASC, DAM 3 DESC.

The two visits are months to a year apart, so a low cosine is a statement about the **record**, not about a laboratory repeat: season, wind climate, chip shape and sub-swath all change between the visits.  Where n is 2 or 3 per visit the point is drawn with a grey ring and its cosine carries almost no information.

## 5. Ordering and cross-plane agreement (descriptive)

| orbit | order by \|v_phys\| | order by \|v_OCV\| | same order |
| --- | --- | --- | --- |
| ASC | DAM3 < DAM4 < DAM6 | DAM4 < DAM3 < DAM6 | **no** |
| DESC | DAM3 < DAM4 < DAM6 | DAM4 < DAM3 < DAM6 | **no** |

Per state, the cosine between the physical and the OCV vector:

| orbit | state | cos(v_phys, v_OCV) | angle [deg] | \|v_phys\| | \|v_OCV\| | n OCV |
| --- | --- | --- | --- | --- | --- | --- |
| ASC | DAM 3 | +0.95 | 17.8 | 0.248 | 0.94 | 8 |
| ASC | DAM 4 | +0.68 | 47.2 | 0.263 | 0.87 | 5 |
| ASC | DAM 6 | -0.52 | 121.4 | 0.309 | 1.34 | 4 |
| DESC | DAM 3 | +1.00 | 3.6 | 0.248 | 0.95 | 17 |
| DESC | DAM 4 | -0.99 | 170.5 | 0.263 | 0.64 | 11 |
| DESC | DAM 6 | +0.50 | 60.2 | 0.309 | 1.23 | 12 |

The informative comparison is not this cosine but the **pairwise angle pattern**, because the three physical vectors are nearly parallel:

| orbit | pair | angle physical [deg] | angle OCV [deg] | \|difference\| |
| --- | --- | --- | --- | --- |
| ASC | DAM3 vs DAM4 | 1.0 | 63.9 | 62.9 |
| ASC | DAM3 vs DAM6 | 7.5 | 96.1 | 88.7 |
| ASC | DAM4 vs DAM6 | 8.5 | 160.1 | 151.6 |
| DESC | DAM3 vs DAM4 | 1.0 | 168.0 | 167.0 |
| DESC | DAM3 vs DAM6 | 7.5 | 71.2 | 63.7 |
| DESC | DAM4 vs DAM6 | 8.5 | 120.8 | 112.4 |

**Answer to the figure question, as far as three states can carry it:** **no.**  The three physical vectors are nearly parallel (1.0°, 7.5° and 8.5° apart), so panel A orders the states by distance along one direction; panel B orders them differently (DAM4 closest to healthy in both orbits), its visit-to-visit repeatability is good only where n is largest, and its vectors change sign with the chip definition.  What the data support is a monotone **dose** response in the physical plane and a state- and orbit-dependent OCV signature — not a shared geometry.

`cos(v_phys, v_OCV)` is reported for completeness.  It is a **descriptive** statistic on three points per orbit: the physical vectors point almost the same way, so the cosine is essentially the direction of the OCV vector, and a small angle between the two planes is therefore not evidence of a shared mechanism.  A matching direction would support *consistency between two measurement channels*; it would **not** show that (gamma2, P) measures stiffness or damping.

## 6. Sensitivity — chip definition and damping mode

Panel B is drawn on the **support** chip.  The two chips define gamma2 differently: the support re-cut measures the echo cells inside the dwell support against a disjoint ring reference and raises gamma2 by an order of magnitude (see `data/lumo/support/README.md`).  Repeating the whole panel on the **committed** chip:

| orbit | state | z (support) | z (committed) | angle between [deg] |
| --- | --- | --- | --- | --- |
| ASC | DAM 3 | +0.01 / -0.94 | -0.08 / -0.87 | 5.8 |
| ASC | DAM 4 | -0.77 / -0.39 | +0.00 / -1.17 | 63.5 |
| ASC | DAM 6 | +1.33 / +0.15 | -0.36 / -1.19 | 113.4 |
| DESC | DAM 3 | -0.23 / -0.92 | +0.12 / -0.43 | 29.0 |
| DESC | DAM 4 | +0.28 / +0.57 | -0.14 / +0.10 | 81.5 |
| DESC | DAM 6 | -1.22 / -0.11 | -0.17 / +1.07 | 86.1 |

The chip definition moves the OCV points by tens of degrees and in places reverses them, so the OCV plane is **definition dependent** at this sample size, and any statement that joins the two planes inherits that dependence.

| state | df/f0 (H1) | dzeta H1 | angle H1 | df/f0 (fund.) | dzeta fund. | angle fund. |
| --- | --- | --- | --- | --- | --- | --- |
| DAM 3 | -7.41 % | -23.68 % | -107.4 | +0.00 % | -2.28 % | -90.0 |
| DAM 4 | -7.41 % | -25.22 % | -106.4 | +0.00 % | +11.04 % | 90.0 |
| DAM 6 | -12.96 % | -28.02 % | -114.8 | +0.00 % | +9.61 % | 90.0 |

The fundamental mode barely moves in frequency but changes its damping in the **opposite direction** for DAM4 and DAM6.  Because its frequency shift is 0.00 % on the 0.25 Hz grid, its direction is ±90°: the damping axis alone, which is exactly why it is a sensitivity and not the primary pair.  H1 is the primary axis only because it is the mode that supplies the frequency shift.

## 7. Cross-checks

The script recomputes every plotted number and checks it: **121 checks, 0 failures** (107 of them verify the dynamics reference itself: medians, spreads, relative changes, campaign names and the agreement with the committed frequency reference).

| check group | checks |
| --- | --- |
| data | 7 |
| panel_a | 2 |
| panel_b | 3 |
| reference | 107 |
| sensitivity | 2 |

## 8. Caveats

- **No physical replicate.** Each state was measured once (5 healthy + 5 damaged recordings in one campaign).  The damping interval is a bootstrap over those ten recordings, not a repeat of the experiment.
- **The physical plane is nearly one-dimensional.** DAM3 and DAM4 sit almost on top of each other (identical 1 Hz H1 shift on the 0.25 Hz grid, damping differing by 1.5 percentage points), so the plane cannot separate them; the ordering by |v| rests on DAM6 being about twice as far out.
- **The H1 frequency shift is quantised.** The mode consensus resolves 0.25 Hz; the drawn x error is ±0.125 Hz (±0.9 %), and the identical −7.41 % of DAM3 and DAM4 is the same bin difference.
- **Thin visit points.** ASC visits carry 2–5 dates each; only the DESC visits reach n = 5–12.  Points with n < 3 are ringed.
- **The OCV plane depends on the chip definition** (section 6), and the two orbits are not comparable: they differ in viewing geometry and sub-swath, and are therefore never pooled.
- **Damping mode matters.** The H1 damping falls for all three states, the fundamental-mode damping rises for DAM4 and DAM6 (section 6).  The second-higher mode is carried in the reference as context only; its DAM3 value (a factor of 13) is implausible and shows how unstable that estimate is.
- **State and season are confounded**, as in the rest of the LUMO record: the visits are separated by months to a year.

---

Generated by `code/lumo/fig_lumo_ocv_physics.py` (121 cross-checks, 0 failures).
