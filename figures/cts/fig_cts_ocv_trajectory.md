# Champlain Towers South, Surfside — OCV trajectory on the support-sized chip

**Site** Champlain Towers South, Surfside · **event** 2021-06-24 · **reference state** `pre-collapse (healthy)` (pooled over the orbits) · **table** `data/cts/support/cts_ocv_channels.csv` (sha256 `3965a15b2053dc0b`)

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
| frozen pre-event reference (chosen) | `pre-collapse (healthy)` | 2015-09-21 → 2021-06-21 | **136** | 136 | 0 | 136 | 62 |
| tail | `post-collapse` | 2021-07-03 → 2022-12-25 | **42** | 42 | 0 | 42 | 18 |

The reference window is the state's own span (R2, `None`): no sub-window is chosen, so all **136** chips of `pre-collapse (healthy)` (2015-09-21 → 2021-06-21) standardise, and the plate draws nothing outside it.

The cut holds **178** of the table's **178** rows (0 refused, none excluded); the plate draws 80 monthly bins from 2015-09 to 2022-12.

### The three frozen classes

| class | label in this report | chips | ASC | DESC | echo-bearing | monthly intervals |
| --- | --- | --- | --- | --- | --- | --- |
| `reference` | frozen pre-event reference | **136** | 136 | 0 | 136 | 62 |
| `outside` | held-out pre-event | 0 | – | – | – | 0 |
| `post` | post-event | **42** | 42 | 0 | 42 | 18 |

OUTSIDE is strictly before 2021-06-24; POST is the registry's post state(s); a chip on the event day is assigned by its registry state, and a reference-state chip on or after the event is reported as unassigned.

The classes reconcile with the cut: **178** kept = **136** reference + **0** held-out pre-event + **42** post-event + **0** unassigned. The event month is split into a pre and a post bin wherever the record holds chips on both sides, so no monthly aggregation combines the two.

## The reference (R3 / R4)

**pooled**

| dimension | block | standardised against | n | t | median | scale | from | z range |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `gamma2` | behaviour | echo-bearing reference chips | 136 | log | -4.697 | 1.411 | mad | -3.57 … +1.41 |
| `P` | behaviour | echo-bearing reference chips | 136 | log | -1.283 | 0.2751 | mad | -2.98 … +1.86 |
| `A` | geometry | all reference chips | 136 | log | 5.263 | 0.5876 | mad | -2.81 … +1.85 |
| `D` | geometry | all reference chips | 136 | log | -3.172 | 0.4803 | mad | -2.38 … +2.03 |
| `F` | geometry | all reference chips | 136 | log | 4.025 | 0.3923 | mad | -3.19 … +1.70 |
| `S` | geometry | all reference chips | 136 | identity | 19.95 | 6.248 | mad | -3.02 … +1.98 |

### Environmental context (panel c1)

Standardised on the same reference, on their own (identity) scale.

| column | group | n | distinct | median | scale | from |
| --- | --- | --- | --- | --- | --- | --- |
| `precipitation_mm` | pooled | 136 | 61 | 0.85 | 1.26 | mad |
| `temperature_c` | pooled | 136 | 79 | 24.8 | 3.039 | mad |

## The reference distribution and the monthly envelope

The decision distribution is built at **monthly** aggregation: each of the **62** reference months contributes one state (the monthly median of its chips), and the envelope is the empirical q05/q50/q95 of those states. A held-out month is compared only with this monthly distribution -- never with the chip-level reference and never with the reference maximum, which is reported below for transparency and is **not** a decision boundary. With 62 monthly states the quantiles are empirical order statistics, not confidence bounds.

| quantity | monthly ref states | q05 | q50 | q95 | min | max | used for |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `gamma2` | 62 | -1.34 | -0.12 | +0.89 | -2.74 | +1.41 | outside [q05, q95] = signed departure |
| `P` | 62 | -1.62 | -0.02 | +0.90 | -2.78 | +1.62 | outside [q05, q95] = signed departure |
| `A` | 62 | -1.17 | -0.17 | +1.09 | -1.80 | +1.53 | outside [q05, q95] = signed departure |
| `D` | 62 | -1.16 | -0.09 | +1.16 | -1.61 | +1.64 | outside [q05, q95] = signed departure |
| `F` | 62 | -1.34 | -0.05 | +1.15 | -2.06 | +1.54 | outside [q05, q95] = signed departure |
| `S` | 62 | -1.98 | -0.21 | +0.79 | -2.86 | +1.21 | outside [q05, q95] = signed departure |
| `R_B` | 62 | +0.31 | +0.76 | +1.69 | +0.18 | +2.18 | above q95 = unusual magnitude |
| `R_G` | 62 | +0.34 | +0.86 | +1.58 | +0.26 | +1.96 | above q95 = unusual magnitude |

The chip-level distributions are stored beside them in the JSON (`reference_envelope.chip`) and are **descriptive only**: nothing in the held-out assessment compares a monthly state with a single chip.

## The monthly trajectory

| month | n | ASC | DESC | echo | z(γ²) | z(P) | z(A) | z(D) | z(F) | z(S) | R_B | R_G |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 2015-09 | 1 | 1 | 0 | 1 | +1.41 | +0.55 | -1.80 | -1.61 | -1.35 | +0.52 | +1.07 | +1.41 |
| 2015-10 | 1 | 1 | 0 | 1 | +0.84 | -2.13 | -0.32 | -0.44 | -0.09 | -0.02 | +1.62 | +0.28 |
| 2016-04 | 1 | 1 | 0 | 1 | -1.32 | -2.78 | +1.53 | +1.64 | +1.12 | -1.86 | +2.18 | +1.56 |
| 2016-07 | 1 | 1 | 0 | 1 | -0.15 | -1.15 | +1.22 | +1.26 | +1.50 | -1.52 | +0.82 | +1.38 |
| 2016-09 | 1 | 1 | 0 | 1 | +1.03 | +1.62 | -1.80 | -1.57 | -2.06 | -0.01 | +1.36 | +1.58 |
| 2016-10 | 3 | 3 | 0 | 3 | +0.14 | -0.60 | +0.04 | -0.05 | +0.26 | -1.59 | +0.76 | +1.12 |
| 2016-11 | 3 | 3 | 0 | 3 | -2.04 | -1.63 | +0.42 | +0.41 | +0.22 | +0.05 | +1.85 | +0.60 |
| 2016-12 | 2 | 2 | 0 | 2 | -0.27 | -0.23 | -0.43 | -0.19 | -0.42 | -1.88 | +0.74 | +1.96 |
| 2017-01 | 3 | 3 | 0 | 3 | -0.53 | +0.32 | -0.38 | -0.62 | +0.05 | -0.86 | +0.61 | +0.75 |
| 2017-02 | 1 | 1 | 0 | 1 | +1.01 | +0.13 | -0.52 | -0.79 | -0.39 | +0.94 | +0.72 | +0.69 |
| 2017-03 | 3 | 3 | 0 | 3 | -0.83 | +0.71 | -0.82 | -0.90 | -0.50 | -0.71 | +0.89 | +0.92 |
| 2017-04 | 1 | 1 | 0 | 1 | +0.04 | -0.66 | +1.09 | +1.17 | +1.24 | -2.86 | +0.47 | +1.75 |
| 2017-05 | 2 | 2 | 0 | 2 | +0.77 | -0.15 | -0.20 | -0.25 | +0.03 | -1.43 | +0.70 | +1.17 |
| 2017-06 | 1 | 1 | 0 | 1 | +0.19 | +0.17 | -0.51 | -0.14 | -0.14 | -0.41 | +0.18 | +0.34 |
| 2017-07 | 1 | 1 | 0 | 1 | +0.25 | -0.76 | +0.05 | -0.12 | +0.38 | -2.63 | +0.57 | +1.33 |
| 2017-08 | 2 | 2 | 0 | 2 | -0.55 | -2.12 | +1.46 | +1.61 | +1.54 | -1.08 | +1.56 | +1.44 |
| 2017-09 | 1 | 1 | 0 | 1 | +0.89 | -1.24 | +0.57 | +0.36 | +0.68 | -0.57 | +1.08 | +0.56 |
| 2017-10 | 3 | 3 | 0 | 3 | -2.74 | -1.18 | +0.56 | +0.50 | +0.68 | -1.75 | +2.14 | +1.01 |
| 2017-11 | 2 | 2 | 0 | 2 | -0.47 | +0.34 | -0.22 | -0.30 | -0.13 | -0.08 | +0.98 | +0.55 |
| 2017-12 | 2 | 2 | 0 | 2 | +0.89 | +0.53 | -0.57 | -0.45 | -0.43 | +0.44 | +0.75 | +0.66 |
| 2018-01 | 2 | 2 | 0 | 2 | -0.98 | -0.58 | +0.44 | +0.35 | +0.40 | +0.29 | +1.06 | +0.39 |
| 2018-02 | 3 | 3 | 0 | 3 | -0.25 | -0.67 | +1.07 | +1.01 | +1.15 | -0.00 | +0.62 | +0.91 |
| 2018-03 | 2 | 2 | 0 | 2 | -0.01 | -0.54 | +0.37 | +0.28 | +0.30 | -0.30 | +0.95 | +0.95 |
| 2018-04 | 3 | 3 | 0 | 3 | -0.53 | +0.30 | -0.46 | -0.60 | -0.39 | +0.41 | +0.76 | +0.80 |
| 2018-05 | 2 | 2 | 0 | 2 | -0.11 | +0.90 | -0.56 | -0.66 | -0.56 | -0.45 | +0.95 | +0.71 |
| 2018-06 | 2 | 2 | 0 | 2 | -0.68 | +0.19 | +0.01 | -0.05 | -0.05 | +0.35 | +0.54 | +0.27 |
| 2018-07 | 3 | 3 | 0 | 3 | -0.15 | -0.34 | +0.11 | +0.05 | +0.09 | +0.32 | +0.37 | +0.83 |
| 2018-08 | 1 | 1 | 0 | 1 | +0.14 | +0.99 | -1.07 | -1.12 | -1.43 | +0.04 | +0.70 | +1.05 |
| 2018-09 | 3 | 3 | 0 | 3 | +0.71 | +0.27 | -0.14 | -0.16 | -0.29 | +0.72 | +0.64 | +1.03 |
| 2018-10 | 2 | 2 | 0 | 2 | -0.33 | +0.51 | -0.39 | -0.19 | -0.81 | -1.98 | +1.03 | +1.72 |
| 2018-11 | 3 | 3 | 0 | 3 | -1.02 | -0.05 | +0.19 | +0.18 | +0.18 | +0.37 | +0.75 | +0.56 |
| 2018-12 | 1 | 1 | 0 | 1 | +0.14 | +0.34 | +0.21 | +0.14 | +0.13 | +0.76 | +0.26 | +0.41 |
| 2019-01 | 3 | 3 | 0 | 3 | +0.87 | +0.42 | -0.81 | -0.92 | -0.45 | +0.79 | +0.69 | +0.80 |
| 2019-02 | 2 | 2 | 0 | 2 | +0.67 | +0.33 | +0.08 | +0.11 | +0.16 | -0.99 | +0.61 | +0.94 |
| 2019-03 | 3 | 3 | 0 | 3 | -0.13 | -0.07 | +0.02 | -0.00 | -0.05 | +0.51 | +0.31 | +0.26 |
| 2019-04 | 2 | 2 | 0 | 2 | -0.03 | +0.74 | -0.69 | -0.61 | -0.61 | -0.50 | +0.60 | +0.63 |
| 2019-05 | 3 | 3 | 0 | 3 | -0.24 | +0.14 | -0.26 | +0.00 | -0.14 | -0.13 | +0.58 | +1.44 |
| 2019-06 | 2 | 2 | 0 | 2 | -0.49 | -0.88 | +0.94 | +0.88 | +0.69 | -0.76 | +0.74 | +0.83 |
| 2019-07 | 3 | 3 | 0 | 3 | +0.85 | -0.37 | -0.38 | -0.16 | -0.19 | +0.72 | +1.36 | +0.64 |
| 2019-08 | 3 | 3 | 0 | 3 | +0.13 | +1.06 | -1.00 | -1.05 | -1.13 | -0.57 | +0.88 | +0.96 |
| 2019-09 | 3 | 3 | 0 | 3 | +0.22 | +0.25 | +0.26 | +0.14 | +0.09 | +0.03 | +0.24 | +0.61 |
| 2019-10 | 3 | 3 | 0 | 3 | -0.69 | -1.30 | +0.81 | +0.82 | +0.68 | +0.02 | +0.94 | +0.68 |
| 2019-11 | 2 | 2 | 0 | 2 | -1.34 | +0.21 | +0.01 | -0.01 | -0.04 | +0.46 | +1.50 | +0.38 |
| 2019-12 | 3 | 3 | 0 | 3 | -0.41 | +0.50 | -0.37 | +0.12 | -0.67 | -0.38 | +1.11 | +1.04 |
| 2020-01 | 2 | 2 | 0 | 2 | -1.94 | +0.70 | -0.40 | -0.49 | -0.32 | -1.07 | +1.69 | +1.23 |
| 2020-02 | 3 | 3 | 0 | 3 | -0.66 | -0.19 | +0.45 | +0.40 | +0.13 | +0.08 | +0.69 | +0.41 |
| 2020-03 | 2 | 2 | 0 | 2 | +0.49 | +0.08 | -0.23 | -0.35 | -0.16 | +0.67 | +0.42 | +0.45 |
| 2020-04 | 3 | 3 | 0 | 3 | +0.47 | +0.01 | -0.40 | -0.46 | -0.34 | +0.51 | +1.01 | +0.49 |
| 2020-05 | 2 | 2 | 0 | 2 | -1.15 | -0.66 | -0.56 | -0.28 | -0.76 | -1.31 | +0.96 | +1.26 |
| 2020-06 | 3 | 3 | 0 | 3 | -0.67 | +0.72 | -0.56 | -0.70 | -0.24 | +0.82 | +0.82 | +0.52 |
| 2020-07 | 2 | 2 | 0 | 2 | +0.00 | -0.08 | -1.18 | -1.16 | -1.01 | -0.62 | +0.92 | +1.07 |
| 2020-08 | 3 | 3 | 0 | 3 | +0.64 | +0.46 | -0.12 | +0.04 | +0.09 | -1.22 | +0.87 | +1.45 |
| 2020-09 | 2 | 2 | 0 | 2 | +0.66 | -0.30 | +0.01 | +0.09 | +0.14 | -1.31 | +0.95 | +0.93 |
| 2020-10 | 2 | 2 | 0 | 2 | +0.76 | +0.54 | -0.84 | -0.64 | -0.82 | +1.21 | +0.67 | +0.92 |
| 2020-11 | 2 | 2 | 0 | 2 | -0.77 | -0.31 | +0.80 | +0.85 | +0.86 | -0.92 | +0.76 | +0.89 |
| 2020-12 | 2 | 2 | 0 | 2 | -0.49 | -0.54 | +0.02 | -0.01 | +0.30 | -0.06 | +0.62 | +0.38 |
| 2021-01 | 2 | 2 | 0 | 2 | -1.24 | -0.35 | +0.73 | +0.76 | +0.90 | -1.67 | +0.95 | +1.26 |
| 2021-02 | 2 | 2 | 0 | 2 | +0.72 | +0.74 | -0.80 | -0.89 | -0.76 | +0.16 | +0.73 | +0.82 |
| 2021-03 | 3 | 3 | 0 | 3 | -0.42 | -0.25 | +0.66 | +0.68 | +0.84 | -2.78 | +0.60 | +1.53 |
| 2021-04 | 2 | 2 | 0 | 2 | +0.77 | +0.21 | -1.48 | -1.41 | -1.44 | -0.10 | +0.69 | +1.44 |
| 2021-05 | 3 | 3 | 0 | 3 | -0.95 | -1.04 | +0.53 | +0.52 | +0.84 | -0.33 | +1.01 | +0.59 |
| 2021-06 | 2 | 2 | 0 | 2 | +0.27 | -0.75 | -0.40 | -0.55 | -0.14 | -0.29 | +0.62 | +0.41 |
| 2021-07 | 3 | 3 | 0 | 3 | -0.01 | -0.47 | +0.01 | -0.09 | -0.29 | -1.27 | +0.81 | +0.65 |
| 2021-08 | 2 | 2 | 0 | 2 | +0.54 | -0.88 | +0.01 | +0.01 | +0.10 | -2.08 | +0.93 | +1.33 |
| 2021-09 | 3 | 3 | 0 | 3 | -0.28 | +0.54 | -0.91 | -1.11 | -0.73 | -0.54 | +0.46 | +0.87 |
| 2021-10 | 2 | 2 | 0 | 2 | -0.79 | -0.23 | +0.74 | +0.79 | +0.79 | -2.65 | +0.64 | +1.49 |
| 2021-11 | 2 | 2 | 0 | 2 | -0.40 | -1.59 | +0.90 | +0.90 | +0.90 | -1.91 | +1.37 | +1.41 |
| 2021-12 | 3 | 3 | 0 | 3 | +0.13 | +0.33 | -0.88 | -0.61 | -0.50 | -1.11 | +0.67 | +0.97 |
| 2022-01 | 2 | 2 | 0 | 2 | +0.73 | -0.41 | -0.73 | -0.64 | -0.44 | +0.20 | +1.27 | +1.07 |
| 2022-02 | 3 | 3 | 0 | 3 | -0.42 | +0.04 | +0.24 | +0.33 | +0.22 | -0.51 | +0.82 | +0.76 |
| 2022-03 | 2 | 2 | 0 | 2 | +0.30 | -1.19 | +0.96 | +1.11 | +0.92 | -2.87 | +0.92 | +1.69 |
| 2022-04 | 3 | 3 | 0 | 3 | +0.24 | -0.71 | +0.51 | +0.44 | +0.46 | -0.47 | +0.51 | +0.56 |
| 2022-05 | 1 | 1 | 0 | 1 | -0.13 | -0.46 | +0.50 | +0.69 | +0.64 | -2.79 | +0.34 | +1.49 |
| 2022-06 | 3 | 3 | 0 | 3 | +0.22 | -0.43 | +0.22 | +0.38 | +0.60 | +0.02 | +1.14 | +1.22 |
| 2022-07 | 2 | 2 | 0 | 2 | -0.03 | -0.48 | -0.32 | +0.15 | -0.57 | +0.97 | +0.39 | +1.31 |
| 2022-08 | 3 | 3 | 0 | 3 | -0.67 | -1.45 | +1.07 | +1.14 | +0.84 | -0.48 | +1.90 | +1.01 |
| 2022-09 | 1 | 1 | 0 | 1 | -1.11 | +0.42 | -1.52 | -1.29 | -1.86 | +0.79 | +0.84 | +1.42 |
| 2022-10 | 2 | 2 | 0 | 2 | +0.30 | -0.43 | -0.41 | -0.33 | -0.43 | +0.22 | +0.62 | +0.76 |
| 2022-11 | 2 | 2 | 0 | 2 | -1.80 | +0.06 | -0.02 | -0.00 | -0.03 | +0.09 | +1.30 | +0.32 |
| 2022-12 | 3 | 3 | 0 | 3 | -0.12 | -1.49 | +0.97 | +0.95 | +1.00 | +0.27 | +1.06 | +0.83 |

A block magnitude is computed **per chip** (the RMS of that chip's standardised dimensions, R4) and the panel and the table show the **median of those chip values** — so `R_B` is the typical joint displacement of a single overpass. It is therefore *not* the RMS of the two columns beside it, and it does not vanish when those two medians are near zero: a month whose chips scatter widely in opposite directions has a zero median in each dimension and a large magnitude per chip.

## Frozen reference / held-out pre-event / post-event

Each cell is the median over the chips of that class, in the reference's own scales. The reference column is zero by construction (the standardisation rests on its median); it is printed so the other two columns are read against the frame they are measured in.

| quantity | frozen pre-event reference | held-out pre-event | post-event |
| --- | --- | --- | --- |
| `gamma2` | +0.00 | – | -0.11 |
| `P` | +0.00 | – | -0.47 |
| `R_B` | +0.80 | – | +0.82 |
| `A` | +0.00 | – | +0.36 |
| `D` | +0.00 | – | +0.35 |
| `F` | +0.00 | – | +0.51 |
| `S` | +0.00 | – | -0.52 |
| `R_G` | +0.83 | – | +0.99 |

n = **136** / **0** / **42** chips over **62** / **0** / **18** monthly intervals.

## The contrast (post minus reference, in reference scales)

Circular moving-block bootstrap, **2000** draws, block **3**, seed **7**, one resample per group; the reference side holds **136** chips and the tail **42**.

| quantity | n_ref | n_post | delta z | 95 % CI | verdict |
| --- | --- | --- | --- | --- | --- |
| `gamma2` | 136 | 42 | -0.11 | -0.37 … +0.21 | spans zero |
| `P` | 136 | 42 | -0.47 | -0.87 … -0.12 | **excludes zero** |
| `A` | 136 | 42 | +0.36 | +0.07 … +0.77 | **excludes zero** |
| `D` | 136 | 42 | +0.35 | +0.03 … +0.68 | **excludes zero** |
| `F` | 136 | 42 | +0.51 | -0.00 … +0.81 | spans zero |
| `S` | 136 | 42 | -0.52 | -1.21 … +0.00 | spans zero |
| `R_B` | 136 | 42 | +0.02 | -0.17 … +0.21 | spans zero |
| `R_G` | 136 | 42 | +0.16 | -0.05 … +0.54 | spans zero |

The interval is on the **difference of medians**, not on the curve: a month that moves inside it is not a finding, and a month that leaves it is not yet one either — this record carries one reference and one tail, not a series of states.

## The held-out pre-event class (OUTSIDE)

the reference window is the state's own span (R2), so no healthy chip sits outside it and there is nothing for the plate to test.

## Held-out pre-event assessment

there is no held-out pre-event class on this site: the reference window is the reference state's own span, so there is nothing to assess against it.

### Summer 2024 — post-hoc descriptive subset

no held-out pre-event class on this site.

### Environmental context of the held-out departures

No held-out class, so there is no held-out departure to compare with the environment.

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

`--verify` re-derives every number of this report from `data/cts/support/cts_ocv_channels.csv` and fails if the JSON, the report and the table disagree; the table's sha256 is `3965a15b2053dc0b28a0cdd741092a415e7acd2bd743a788f2462d04bda28826`.

The assessment's pins are constants, asserted by `--selftest`: seed **7**, block **3**, draws **2000**, decision quantiles **q05/q50/q95**, the magnitude boundary **q95**, and a persistence rule of **2** consecutive months. The bootstrap is seeded and no wall-clock value enters any output, so the JSON, the report and the plate reproduce deterministically.
