# Carolabrücke, Dresden — OCV trajectory on the support-sized chip

**Site** Carolabrücke, Dresden · **event** 2024-09-11 · **reference state** `pre-collapse (healthy)` (pooled over the orbits) · **table** `data/carola/support/carola_ocv_channels.csv` (sha256 `3137d0b31fef80b5`)

```
standardisation  z = (t(x) - median_ref) / (1.4826 * MAD_ref),  |z| <= 6
                 t = log for gamma2, P, A, D, F; identity for S
blocks           behaviour = ['gamma2', 'P']   geometry = ['A', 'D', 'F', 'S']   (ocv_reconfigure.BLOCKS)
reference per block  behaviour: echo-bearing reference chips; geometry: all reference chips (R3)
bootstrap        2000 circular moving-block draws, block 3, seed 7 (R5)
```

## The window

| part | state | window | chips | ASC | DESC | echo-bearing | months |
| --- | --- | --- | --- | --- | --- | --- | --- |
| frozen pre-event reference (chosen) | `pre-collapse (healthy)` | 2022-09-03 → 2023-12-27 | **121** | 81 | 40 | 121 | 16 |
| held-out pre-event | `pre-collapse (healthy)` | 2024-01-01 → 2024-09-09 | 63 | 43 | 20 | 63 | 9 |
| frozen pre-event state | `pre-collapse (healthy)` | 2022-09-03 → 2024-09-09 | 184 | 124 | 60 | 184 | 25 |
| tail | `post-collapse` | 2024-09-13 → 2026-08-23 | **252** | 166 | 86 | 252 | 24 |

The reference window is a **stated choice per asset** (R2): `2022-09` → `2023-12`, inclusive calendar months. The `pre-collapse (healthy)` state holds **184** chips in all (2022-09-03 → 2024-09-09); the choice keeps **121** and leaves **63** chips outside it on the pre-event side (2024-01-01 → 2024-09-09) — drawn in grey, standardising nothing.

The cut holds **436** of the table's **436** rows (0 refused, none excluded); the plate draws 49 monthly bins from 2022-09 to 2026-08.

### The three frozen classes

| class | label in this report | chips | ASC | DESC | echo-bearing | monthly intervals |
| --- | --- | --- | --- | --- | --- | --- |
| `reference` | frozen pre-event reference | **121** | 81 | 40 | 121 | 16 |
| `outside` | held-out pre-event | **63** | 43 | 20 | 63 | 9 |
| `post` | post-event | **252** | 166 | 86 | 252 | 24 |

OUTSIDE is strictly before 2024-09-11; POST is the registry's post state(s); a chip on the event day is assigned by its registry state, and a reference-state chip on or after the event is reported as unassigned.

The classes reconcile with the cut: **436** kept = **121** reference + **63** held-out pre-event + **252** post-event + **0** unassigned. The event month is split into a pre and a post bin wherever the record holds chips on both sides, so no monthly aggregation combines the two.

## The reference (R3 / R4)

**pooled**

| dimension | block | standardised against | n | t | median | scale | from | z range |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `gamma2` | behaviour | echo-bearing reference chips | 121 | log | -6.105 | 1.374 | mad | -3.63 … +2.33 |
| `P` | behaviour | echo-bearing reference chips | 121 | log | -2.203 | 0.3527 | mad | -2.47 … +1.91 |
| `A` | geometry | all reference chips | 121 | log | 5.421 | 0.981 | mad | -2.09 … +1.57 |
| `D` | geometry | all reference chips | 121 | log | -3.642 | 0.9082 | mad | -1.74 … +1.69 |
| `F` | geometry | all reference chips | 121 | log | 4.304 | 0.7773 | mad | -2.34 … +1.37 |
| `S` | geometry | all reference chips | 121 | identity | 14.67 | 6.643 | mad | -2.12 … +3.29 |

### Environmental context (panel c1)

Standardised on the same reference, on their own (identity) scale.

| column | group | n | distinct | median | scale | from |
| --- | --- | --- | --- | --- | --- | --- |
| `temperature_c` | pooled | 121 | 98 | 11.4 | 7.858 | mad |
| `wind_speed_ms` | pooled | 121 | 110 | 2.51 | 1.394 | mad |

## The reference distribution and the monthly envelope

The decision distribution is built at **monthly** aggregation: each of the **16** reference months contributes one state (the monthly median of its chips), and the envelope is the empirical q05/q50/q95 of those states. A held-out month is compared only with this monthly distribution -- never with the chip-level reference and never with the reference maximum, which is reported below for transparency and is **not** a decision boundary. With 16 monthly states the quantiles are empirical order statistics, not confidence bounds.

| quantity | monthly ref states | q05 | q50 | q95 | min | max | used for |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `gamma2` | 16 | -0.51 | +0.16 | +0.70 | -1.15 | +0.99 | outside [q05, q95] = signed departure |
| `P` | 16 | -0.33 | +0.11 | +0.56 | -0.36 | +0.83 | outside [q05, q95] = signed departure |
| `A` | 16 | -0.69 | +0.04 | +0.58 | -0.93 | +0.71 | outside [q05, q95] = signed departure |
| `D` | 16 | -0.68 | +0.05 | +0.63 | -0.80 | +0.77 | outside [q05, q95] = signed departure |
| `F` | 16 | -0.81 | -0.07 | +0.35 | -1.25 | +0.39 | outside [q05, q95] = signed departure |
| `S` | 16 | -0.52 | +0.07 | +0.64 | -0.79 | +0.82 | outside [q05, q95] = signed departure |
| `R_B` | 16 | +0.42 | +0.82 | +1.20 | +0.36 | +1.51 | above q95 = unusual magnitude |
| `R_G` | 16 | +0.47 | +0.83 | +0.99 | +0.39 | +1.01 | above q95 = unusual magnitude |

The chip-level distributions are stored beside them in the JSON (`reference_envelope.chip`) and are **descriptive only**: nothing in the held-out assessment compares a monthly state with a single chip.

## The monthly trajectory

| month | n | ASC | DESC | echo | z(γ²) | z(P) | z(A) | z(D) | z(F) | z(S) | R_B | R_G |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 2022-09 | 7 | 5 | 2 | 7 | +0.09 | +0.09 | -0.61 | -0.64 | -0.38 | +0.58 | +1.00 | +0.87 |
| 2022-10 | 8 | 5 | 3 | 8 | +0.23 | +0.28 | -0.25 | -0.19 | -0.25 | +0.82 | +0.44 | +0.67 |
| 2022-11 | 7 | 5 | 2 | 7 | +0.29 | +0.45 | -0.13 | -0.12 | -0.67 | +0.28 | +0.61 | +0.72 |
| 2022-12 | 8 | 5 | 3 | 8 | -0.25 | -0.22 | +0.09 | +0.10 | -0.06 | +0.52 | +0.57 | +1.01 |
| 2023-01 | 8 | 6 | 2 | 8 | +0.24 | -0.33 | +0.46 | +0.50 | +0.39 | +0.13 | +0.65 | +0.82 |
| 2023-02 | 7 | 4 | 3 | 7 | -0.09 | +0.37 | +0.30 | +0.33 | +0.08 | -0.07 | +0.69 | +0.49 |
| 2023-03 | 8 | 6 | 2 | 8 | -1.15 | -0.36 | +0.71 | +0.77 | +0.33 | -0.34 | +0.91 | +0.70 |
| 2023-04 | 7 | 4 | 3 | 7 | +0.56 | +0.47 | -0.18 | -0.19 | -0.27 | +0.38 | +1.10 | +0.87 |
| 2023-05 | 8 | 6 | 2 | 8 | +0.03 | +0.29 | -0.41 | -0.37 | -0.67 | +0.17 | +0.53 | +0.91 |
| 2023-06 | 8 | 5 | 3 | 8 | +0.12 | +0.17 | +0.22 | +0.24 | -0.08 | -0.06 | +0.76 | +0.85 |
| 2023-07 | 7 | 5 | 2 | 7 | +0.36 | +0.83 | -0.93 | -0.80 | -1.25 | +0.52 | +1.08 | +0.97 |
| 2023-08 | 8 | 5 | 3 | 8 | +0.02 | -0.30 | +0.48 | +0.52 | +0.33 | -0.30 | +1.06 | +0.86 |
| 2023-09 | 7 | 5 | 2 | 7 | +0.61 | -0.16 | +0.54 | +0.59 | +0.24 | -0.79 | +1.51 | +0.98 |
| 2023-10 | 8 | 5 | 3 | 8 | +0.20 | +0.14 | -0.43 | -0.40 | -0.39 | -0.17 | +0.95 | +0.61 |
| 2023-11 | 8 | 5 | 3 | 8 | -0.29 | +0.09 | +0.02 | +0.03 | +0.14 | -0.43 | +0.36 | +0.52 |
| 2023-12 | 7 | 5 | 2 | 7 | +0.99 | -0.09 | +0.05 | +0.06 | +0.08 | +0.01 | +0.88 | +0.39 |
| 2024-01 | 8 | 5 | 3 | 8 | -0.57 | -0.37 | +0.16 | +0.21 | +0.21 | -0.58 | +0.77 | +0.44 |
| 2024-02 | 7 | 5 | 2 | 7 | -1.26 | -0.73 | +0.58 | +0.63 | +0.62 | -0.33 | +1.03 | +0.76 |
| 2024-03 | 8 | 5 | 3 | 8 | +0.26 | +0.12 | -0.28 | -0.27 | -0.34 | -0.03 | +1.19 | +0.63 |
| 2024-04 | 8 | 6 | 2 | 8 | +0.37 | -0.22 | +0.17 | +0.20 | +0.14 | +0.06 | +0.71 | +0.58 |
| 2024-05 | 8 | 5 | 3 | 8 | +0.21 | -0.42 | -0.57 | -0.57 | -0.55 | +0.91 | +0.81 | +0.74 |
| 2024-06 | 7 | 5 | 2 | 7 | +0.73 | -0.71 | -0.40 | -0.41 | -0.58 | +0.79 | +1.24 | +0.47 |
| 2024-07 | 8 | 5 | 3 | 8 | +0.57 | -0.49 | +0.05 | +0.05 | +0.05 | +0.42 | +0.77 | +0.66 |
| 2024-08 | 7 | 5 | 2 | 7 | +0.78 | +0.29 | +0.02 | +0.02 | -0.27 | -0.05 | +0.82 | +0.72 |
| 2024-09 (pre) | 2 | 2 | 0 | 2 | -0.17 | -0.02 | +0.14 | +0.19 | +0.21 | +0.34 | +0.38 | +0.39 |
| 2024-09 (post) | 5 | 3 | 2 | 5 | +0.01 | -0.48 | +0.91 | +0.99 | +0.53 | -0.34 | +0.55 | +0.90 |
| 2024-10 | 8 | 5 | 3 | 8 | +0.48 | +1.05 | -0.71 | -0.65 | -0.91 | +1.80 | +1.02 | +1.18 |
| 2024-11 | 7 | 5 | 2 | 7 | -0.40 | -0.16 | +0.56 | +0.60 | +0.51 | -0.38 | +0.53 | +0.63 |
| 2024-12 | 8 | 5 | 3 | 8 | +0.00 | +0.32 | -0.06 | -0.06 | -0.09 | -0.14 | +0.45 | +0.59 |
| 2025-01 | 8 | 6 | 2 | 8 | +0.31 | +0.02 | +0.21 | +0.23 | +0.16 | +0.20 | +0.53 | +0.46 |
| 2025-02 | 7 | 4 | 3 | 7 | +0.05 | +0.26 | -0.11 | -0.12 | -0.05 | -1.05 | +0.76 | +0.87 |
| 2025-03 | 8 | 5 | 3 | 8 | +0.16 | -0.13 | +0.09 | +0.10 | +0.20 | -0.16 | +0.67 | +0.51 |
| 2025-04 | 12 | 8 | 4 | 12 | +0.61 | -0.40 | -0.02 | +0.01 | +0.07 | +0.31 | +0.98 | +0.79 |
| 2025-05 | 12 | 8 | 4 | 12 | +0.74 | -0.34 | -0.32 | -0.30 | -0.37 | +0.92 | +0.92 | +0.81 |
| 2025-06 | 12 | 8 | 4 | 12 | +0.20 | -0.65 | -0.82 | -0.82 | -0.79 | +0.89 | +1.02 | +1.29 |
| 2025-07 | 12 | 8 | 4 | 12 | +0.52 | -0.46 | -0.08 | -0.05 | -0.12 | +0.85 | +1.07 | +0.81 |
| 2025-08 | 12 | 8 | 4 | 12 | +0.37 | -0.48 | -0.33 | -0.27 | -0.45 | +1.01 | +0.69 | +0.93 |
| 2025-09 | 12 | 8 | 4 | 12 | -0.70 | -0.75 | +0.04 | +0.28 | -0.30 | +0.04 | +1.18 | +0.82 |
| 2025-10 | 12 | 8 | 4 | 12 | +0.31 | -0.32 | -0.30 | -0.23 | -0.47 | +0.87 | +0.93 | +0.79 |
| 2025-11 | 12 | 8 | 4 | 12 | -0.15 | -0.22 | +0.08 | +0.09 | +0.03 | -0.00 | +1.05 | +0.65 |
| 2025-12 | 12 | 8 | 4 | 12 | -0.00 | -0.48 | -0.03 | -0.03 | -0.13 | +1.07 | +0.81 | +0.58 |
| 2026-01 | 12 | 8 | 4 | 12 | +0.52 | -0.32 | -0.16 | -0.14 | -0.29 | -0.21 | +0.73 | +0.55 |
| 2026-02 | 12 | 8 | 4 | 12 | +0.86 | -0.04 | -0.51 | -0.48 | -0.53 | +0.50 | +1.18 | +0.71 |
| 2026-03 | 12 | 8 | 4 | 12 | +0.25 | -0.41 | -0.25 | -0.17 | -0.24 | +0.28 | +1.20 | +0.49 |
| 2026-04 | 12 | 8 | 4 | 12 | +0.88 | -0.59 | -0.55 | -0.50 | -0.61 | +1.13 | +1.07 | +0.87 |
| 2026-05 | 11 | 7 | 4 | 11 | +0.72 | -1.08 | +0.18 | +0.21 | -0.05 | +0.05 | +0.92 | +0.67 |
| 2026-06 | 10 | 6 | 4 | 10 | +0.11 | -1.53 | -0.46 | -0.46 | -0.49 | +0.64 | +1.54 | +0.97 |
| 2026-07 | 12 | 8 | 4 | 12 | +1.31 | -0.79 | -0.28 | -0.17 | -0.45 | +0.61 | +1.04 | +0.73 |
| 2026-08 | 12 | 8 | 4 | 12 | +0.89 | -0.98 | -0.86 | -0.81 | -1.01 | +1.08 | +1.06 | +1.05 |

A block magnitude is computed **per chip** (the RMS of that chip's standardised dimensions, R4) and the panel and the table show the **median of those chip values** — so `R_B` is the typical joint displacement of a single overpass. It is therefore *not* the RMS of the two columns beside it, and it does not vanish when those two medians are near zero: a month whose chips scatter widely in opposite directions has a zero median in each dimension and a large magnitude per chip.

## Frozen reference / held-out pre-event / post-event

Each cell is the median over the chips of that class, in the reference's own scales. The reference column is zero by construction (the standardisation rests on its median); it is printed so the other two columns are read against the frame they are measured in.

| quantity | frozen pre-event reference | held-out pre-event | post-event |
| --- | --- | --- | --- |
| `gamma2` | +0.00 | +0.20 | +0.36 |
| `P` | +0.00 | -0.27 | -0.44 |
| `R_B` | +0.83 | +0.84 | +0.94 |
| `A` | +0.00 | +0.02 | -0.13 |
| `D` | +0.00 | +0.02 | -0.11 |
| `F` | +0.00 | -0.07 | -0.20 |
| `S` | +0.00 | +0.04 | +0.56 |
| `R_G` | +0.76 | +0.71 | +0.75 |

n = **121** / **63** / **252** chips over **16** / **9** / **24** monthly intervals.

## The contrast (post minus reference, in reference scales)

Circular moving-block bootstrap, **2000** draws, block **3**, seed **7**, one resample per group; the reference side holds **121** chips and the tail **252**.

| quantity | n_ref | n_post | delta z | 95 % CI | verdict |
| --- | --- | --- | --- | --- | --- |
| `gamma2` | 121 | 252 | +0.36 | +0.07 … +0.60 | **excludes zero** |
| `P` | 121 | 252 | -0.44 | -0.67 … -0.19 | **excludes zero** |
| `A` | 121 | 252 | -0.13 | -0.51 … +0.02 | spans zero |
| `D` | 121 | 252 | -0.11 | -0.51 … +0.03 | spans zero |
| `F` | 121 | 252 | -0.20 | -0.40 … +0.08 | spans zero |
| `S` | 121 | 252 | +0.56 | +0.08 … +0.84 | **excludes zero** |
| `R_B` | 121 | 252 | +0.12 | -0.04 … +0.30 | spans zero |
| `R_G` | 121 | 252 | -0.02 | -0.11 … +0.07 | spans zero |

The interval is on the **difference of medians**, not on the curve: a month that moves inside it is not a finding, and a month that leaves it is not yet one either — this record carries one reference and one tail, not a series of states.

## The held-out pre-event class (OUTSIDE)

OUTSIDE is the third frozen class: **63** chips of `pre-collapse (healthy)` that the reference window leaves out on the pre-event side (2024-01-01 → 2024-09-09), over **9** monthly intervals (2024-01, 2024-02, 2024-03, 2024-04, 2024-05, 2024-06, 2024-07, 2024-08, 2024-09 (pre)). They standardise nothing -- the centre and the scale are the reference's -- so every number below reads them in the frozen reference's own frame. The interval is the plate's own bootstrap (**2000** draws, block **3**, seed **7**, one resample per side).

**Signed channels** -- `r` is the signed median departure in transform units, `z = r / scale_ref` the same departure in reference scales, and the interval is on the difference of medians.

| quantity | n chips | n int | r | z | 95 % CI (z) | q05 | q95 | monthly level | position | months outside envelope |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `gamma2` | 63 | 9 | +0.280 | +0.20 | -0.09 … +0.60 | -0.51 | +0.70 | +0.26 | inside q05..q95 | 2024-01, 2024-02, 2024-06, 2024-08 |
| `P` | 63 | 9 | -0.095 | -0.27 | -0.74 … +0.05 | -0.33 | +0.56 | -0.37 | below q05 | 2024-01, 2024-02, 2024-05, 2024-06, 2024-07 |
| `A` | 63 | 9 | +0.018 | +0.02 | -0.42 … +0.22 | -0.69 | +0.58 | +0.05 | inside q05..q95 | none |
| `D` | 63 | 9 | +0.018 | +0.02 | -0.44 … +0.23 | -0.68 | +0.63 | +0.05 | inside q05..q95 | none |
| `F` | 63 | 9 | -0.056 | -0.07 | -0.34 … +0.29 | -0.81 | +0.35 | +0.05 | inside q05..q95 | 2024-02 |
| `S` | 63 | 9 | +0.288 | +0.04 | -0.28 … +0.67 | -0.52 | +0.64 | +0.06 | inside q05..q95 | 2024-01, 2024-05, 2024-06 |

**Block magnitudes** -- a magnitude has no sign, so the envelope test is one-sided: only a monthly value above q95 counts as an exceedance. A value below q05 is reported as a *lower* magnitude, never as an anomaly. `position` is where the class's median monthly state sits; the last-but-one column lists the individual months.

| quantity | n chips | n int | median R | delta | 95 % CI | q95 | monthly level | position | months above q95 | longest run | sustained |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `R_B` | 63 | 9 | 0.837 | +0.01 | -0.15 … +0.25 | +1.20 | +0.81 | inside q05..q95 | 2024-06 | 1 (2024-06 → 2024-06) | no |
| `R_G` | 63 | 9 | 0.711 | -0.05 | -0.25 … +0.01 | +0.99 | +0.63 | inside q05..q95 | none | none | no |

**Criterion (fixed before the numbers were read).** a held-out month is an unusual departure only when its monthly value exceeds the monthly reference q95 of the same quantity (magnitudes, one-sided) or leaves [q05, q95] (signed z); a run of at least 2 consecutive such months is sustained.

**Directional pattern (descriptive only).** The joint pattern `z(gamma2) > 0 and z(P) < 0` holds in **4** held-out months; its longest run is **4** month(s) (2024-04 → 2024-07). No statistic is derived from these runs.

## Held-out pre-event assessment

The three statements are kept apart throughout:

* **Directional reconfiguration** -- a signed monthly median leaves the monthly reference envelope; a statement about direction, not size.
* **Reconfiguration magnitude** -- R_B or R_G above the monthly reference q95; a statement about the joint size of the displacement, with no direction.
* **Unusual reference departure** -- only the magnitude statement, and only when the pre-specified criterion is met; a run of at least 2 consecutive months makes it sustained.

**Q. Is there directional behaviour reconfiguration in the held-out period?**

Yes - directional reconfiguration is visible: z(P) lies below its monthly reference envelope in 5 of 9 held-out months (2024-01, 2024-02, 2024-05, 2024-06, 2024-07); z(gamma2) lies above it in 2024-06, 2024-08; z(gamma2) lies below it in 2024-01, 2024-02; the joint pattern z(gamma2) > 0 and z(P) < 0 holds for 4 consecutive months (2024-04 to 2024-07). This is a statement about direction, not about the size of the displacement, and it is not an anomaly verdict by itself.

**Q. Does the held-out behaviour magnitude (R_B) exceed the monthly reference envelope?**

Yes - 1 of 9 held-out months lies above the monthly reference q95 of `R_B` (2024-06); the longest consecutive run is 1 month(s) (2024-06 to 2024-06), so the departure is not sustained under the pre-specified rule.

**Q. Does the held-out geometry magnitude (R_G) exceed the monthly reference envelope?**

No - no held-out month exceeds the monthly reference q95 of `R_G`: the highest monthly value is 0.756 against a threshold of 0.986.

**Q. Is any held-out departure temporally persistent?**

`R_B`: 1 of 9 months above the q95 threshold, longest run 1 month(s) (2024-06 to 2024-06); `R_G`: 0 of 9 months above the q95 threshold, longest run 0 month(s) (no run). Persistence is decided by the fixed rule (at least 2 consecutive months), not by the shape of the curve.

**Q. Does the last summer before the event differ from the rest of the held-out period?**

Post-hoc, descriptive only: the Summer-2024 subset (22 chips, 3 monthly intervals, 2024-06-01 to 2024-08-31) has median z(gamma2) = +0.75, z(P) = -0.30, R_B = 0.836, R_G = 0.661; the remaining held-out chips (41) have +0.00, -0.27, 0.837, 0.711. The difference is not tested, no interval is attached to it, and it is not an input to the verdict.

**Q. Do the pronounced held-out departures coincide with environmental variation?**

The month(s) above the magnitude envelope - 2024-06 - coincide with environmental variation (2024-06: temperature_c z = +1.07, wind_speed_ms z = +0.19). This is a temporal coincidence only: no weather correction is applied and no causal attribution is made.

**Verdict — isolated unusual departure.** It is derived from the R_B and R_G persistence rules alone (`verdict_inputs`: `R_B`, `R_G`); the Summer-2024 subset and the environment are not inputs. It does not claim damage, a collapse prediction, a failure precursor or causation.

### Summer 2024 — post-hoc descriptive subset

**Post-hoc and descriptive only -- no interval, no p-value, no new statistic, no verdict input.** The window `2024-06-01` → `2024-08-31` was fixed after OUTSIDE was defined, and it is a strict subset of it: **22** chips over **3** monthly intervals (2024-06, 2024-07, 2024-08). Medians: z(γ²) +0.75, z(P) -0.30, R_B 0.836, R_G 0.661; months above the monthly q95: R_B 2024-06, R_G none.

### Environmental context of the held-out departures

The co-variates are reported beside the OCV quantities, not inside them: no weather correction is applied and no causal attribution is made.

| month | R_B | R_G | z(γ²) | z(P) | T [C] | wind [m/s] | above q95 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 2024-01 | 0.766 | 0.438 | -0.57 | -0.37 | -1.36 | +0.34 | none |
| 2024-02 | 1.03 | 0.756 | -1.26 | -0.73 | -0.43 | +0.42 | none |
| 2024-03 | 1.19 | 0.629 | +0.26 | +0.12 | -0.41 | -0.10 | none |
| 2024-04 | 0.714 | 0.577 | +0.37 | -0.22 | +0.04 | +0.65 | none |
| 2024-05 | 0.807 | 0.742 | +0.21 | -0.42 | +0.73 | +0.28 | none |
| 2024-06 | 1.24 | 0.473 | +0.73 | -0.71 | +1.07 | +0.19 | R_B |
| 2024-07 | 0.772 | 0.661 | +0.57 | -0.49 | +1.67 | -0.62 | none |
| 2024-08 | 0.817 | 0.725 | +0.78 | +0.29 | +1.77 | +0.31 | none |
| 2024-09 (pre) | 0.381 | 0.387 | -0.17 | -0.02 | +1.81 | +0.53 | none |

The month(s) above the magnitude envelope -- 2024-06 -- coincide with the co-variation recorded above. This is a temporal coincidence, not an explanation: OCV is not weather-corrected in this repository.

## The panels

| panel | block | what it draws |
| --- | --- | --- |
| **a1** | behaviour | z(gamma2) -- monthly median, with the monthly reference envelope |
| **a2** | behaviour | z(P) -- monthly median, with the monthly reference envelope |
| **a3** | behaviour | the behaviour plane (z(gamma2), z(P)) |
| **a4** | behaviour | R_B -- the behaviour block magnitude, with its monthly reference envelope |
| **b1** | geometry | z(A) -- monthly median, with the monthly reference envelope |
| **b2** | geometry | z(F) -- monthly median, with the monthly reference envelope |
| **b3** | geometry | the geometry plane (z(A), z(F)) |
| **b4** | geometry | R_G -- the geometry block magnitude, with its monthly reference envelope |
| **c1** | context | environmental context (weather, descriptive only) |
| **c2** | context | the observation cadence per orbit |
| **c3** | context | observation availability: solid = an echo, hollow = a mask without one |
| **c4** | context | the trajectory matrix: six dimensions x every month |

## What the plate does not claim

* **No damage verdict.** Every ordinate is a deviation from the case's *own* reference, in that reference's robust scales — two plates are not comparable in magnitude, which is why the comparison they support is a comparison of *shape* and of *where the movement sits*.
* **No significance on the curve.** Both intervals of the plate — the contrast's and the held-out class's — are differences of medians, not bands on the curve: the monthly medians are drawn with their `n` and without one.
* **No directional score.** `R_B` and `R_G` keep their frozen definitions; no directional or persistence score is added, and the verdict never reads one.
* **The reference maximum is not a boundary.** The envelope is the q05/q50/q95 of the monthly reference states; the maximum is printed only so the reader can see how far q95 sits from it.
* **The Summer-2024 subset is post-hoc.** It is descriptive, untested and not an input to the verdict.
* **The environment is beside the OCV numbers, never inside them.** No weather correction is applied and no causal attribution is made.
* **No interpolation.** A month without a chip, the split event month and the event itself break the line: nothing is drawn across a hole.
* **A month is a bin, not a date.** Its median can mix both orbits; panel c2 says how many chips of each the bin holds, and on the sites with two orbits the two directions alternate in blocks of days.
* **No pooling across sites.** The plate is per site by construction; `figures/cases/fig_ocv_trajectory.sh` runs the six and nothing else.

## Pins

`--verify` re-derives every number of this report from `data/carola/support/carola_ocv_channels.csv` and fails if the JSON, the report and the table disagree; the table's sha256 is `3137d0b31fef80b5075378cb0f1bc71aff655d0f9c8554bc23fede56b7f4d135`.

The assessment's pins are constants, asserted by `--selftest`: seed **7**, block **3**, draws **2000**, decision quantiles **q05/q50/q95**, the magnitude boundary **q95**, and a persistence rule of **2** consecutive months. The bootstrap is seeded and no wall-clock value enters any output, so the JSON, the report and the plate reproduce deterministically.
