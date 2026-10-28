# KDLO media tower, Garden City SD — OCV trajectory on the support-sized chip

**Site** KDLO media tower, Garden City SD · **event** 2022-12-11 · **reference state** `pre` (pooled over the orbits) · **table** `data/kdlo/support/kdlo_ocv_channels.csv` (sha256 `9c6bb9197197b382`)

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
| frozen pre-event reference (chosen) | `pre` | 2022-06-02 → 2022-12-11 | **16** | 0 | 16 | 16 | 7 |
| tail | `rebuild` | 2023-05-04 → 2024-09-19 | **43** | 0 | 43 | 43 | 17 |

The reference window is the state's own span (R2, `None`): no sub-window is chosen, so all **16** chips of `pre` (2022-06-02 → 2022-12-11) standardise, and the plate draws nothing outside it.

The cut holds **59** of the table's **60** rows (0 refused, `single` × 1 excluded); the plate draws 24 monthly bins from 2022-06 to 2024-09.

* `single` — one acquisition at the collapse itself (2022-12-23): one unit is not a state.

### The three frozen classes

| class | label in this report | chips | ASC | DESC | echo-bearing | monthly intervals |
| --- | --- | --- | --- | --- | --- | --- |
| `reference` | frozen pre-event reference | **16** | 0 | 16 | 16 | 7 |
| `outside` | held-out pre-event | 0 | – | – | – | 0 |
| `post` | post-event | **43** | 0 | 43 | 43 | 17 |

OUTSIDE is strictly before 2022-12-11; POST is the registry's post state(s); a chip on the event day is assigned by its registry state, and a reference-state chip on or after the event is reported as unassigned.

The classes reconcile with the cut: **59** kept = **16** reference + **0** held-out pre-event + **43** post-event + **0** unassigned. The event month is split into a pre and a post bin wherever the record holds chips on both sides, so no monthly aggregation combines the two.

## The reference (R3 / R4)

**pooled**

| dimension | block | standardised against | n | t | median | scale | from | z range |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `gamma2` | behaviour | echo-bearing reference chips | 16 | log | -3.15 | 0.8872 | mad | -4.26 … +1.19 |
| `P` | behaviour | echo-bearing reference chips | 16 | log | -0.553 | 0.1011 | mad | -6.78 … +1.05 |
| `A` | geometry | all reference chips | 16 | log | 3.989 | 0.2058 | mad | -1.84 … +6.01 |
| `D` | geometry | all reference chips | 16 | log | -2.657 | 0.6214 | mad | -2.73 … +0.88 |
| `F` | geometry | all reference chips | 16 | log | 2.303 | 0.2703 | mad | -1.32 … +4.19 |
| `S` | geometry | all reference chips | 16 | identity | 22.93 | 0.5429 | mad | -32.86 … +0.93 |

### Environmental context (panel c1)

Standardised on the same reference, on their own (identity) scale.

| column | group | n | distinct | median | scale | from |
| --- | --- | --- | --- | --- | --- | --- |
| `temperature_c` | pooled | 16 | 16 | 18.15 | 9.192 | mad |
| `wind_speed_ms` | pooled | 16 | 16 | 4.69 | 2.565 | mad |

**40** of the standardised values were clipped at ±6; they are counted here and drawn at the clip, never dropped.

## The reference distribution and the monthly envelope

The decision distribution is built at **monthly** aggregation: each of the **7** reference months contributes one state (the monthly median of its chips), and the envelope is the empirical q05/q50/q95 of those states. A held-out month is compared only with this monthly distribution -- never with the chip-level reference and never with the reference maximum, which is reported below for transparency and is **not** a decision boundary. With 7 monthly states the quantiles are empirical order statistics, not confidence bounds.

| quantity | monthly ref states | q05 | q50 | q95 | min | max | used for |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `gamma2` | 7 | -1.63 | -0.19 | +0.50 | -1.72 | +0.53 | outside [q05, q95] = signed departure |
| `P` | 7 | -4.53 | -0.16 | +0.45 | -6.00 | +0.46 | outside [q05, q95] = signed departure |
| `A` | 7 | -0.61 | -0.02 | +4.52 | -0.79 | +6.00 | outside [q05, q95] = signed departure |
| `D` | 7 | -2.18 | -0.45 | +0.56 | -2.73 | +0.65 | outside [q05, q95] = signed departure |
| `F` | 7 | -0.95 | +0.00 | +3.13 | -1.07 | +4.19 | outside [q05, q95] = signed departure |
| `S` | 7 | -5.06 | -1.34 | +0.65 | -6.00 | +0.75 | outside [q05, q95] = signed departure |
| `R_B` | 7 | +0.36 | +0.82 | +3.53 | +0.32 | +4.24 | above q95 = unusual magnitude |
| `R_G` | 7 | +0.50 | +0.94 | +4.07 | +0.45 | +4.92 | above q95 = unusual magnitude |

The chip-level distributions are stored beside them in the JSON (`reference_envelope.chip`) and are **descriptive only**: nothing in the held-out assessment compares a monthly state with a single chip.

## The monthly trajectory

| month | n | ASC | DESC | echo | z(γ²) | z(P) | z(A) | z(D) | z(F) | z(S) | R_B | R_G |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 2022-06 | 2 | 0 | 2 | 2 | -0.42 | -0.34 | +1.06 | -0.45 | -0.66 | -1.36 | +0.46 | +1.25 |
| 2022-07 | 2 | 0 | 2 | 2 | +0.45 | -0.01 | -0.02 | -0.89 | +0.14 | -2.88 | +0.82 | +2.07 |
| 2022-08 | 3 | 0 | 3 | 3 | +0.53 | -0.16 | -0.09 | -0.84 | +0.00 | +0.30 | +0.81 | +0.45 |
| 2022-09 | 3 | 0 | 3 | 3 | +0.13 | +0.44 | -0.18 | -0.00 | -0.39 | +0.75 | +0.32 | +0.64 |
| 2022-10 | 2 | 0 | 2 | 2 | -1.72 | +0.46 | -0.79 | +0.35 | -1.07 | +0.44 | +1.86 | +0.94 |
| 2022-11 | 3 | 0 | 3 | 3 | -1.44 | -1.10 | +1.05 | +0.65 | +0.67 | -1.34 | +1.35 | +0.79 |
| 2022-12 | 1 | 0 | 1 | 1 | -0.19 | -6.00 | +6.00 | -2.73 | +4.19 | -6.00 | +4.24 | +4.92 |
| 2023-05 | 3 | 0 | 3 | 3 | -1.67 | -2.36 | +2.99 | -0.10 | +2.17 | -4.83 | +3.46 | +3.44 |
| 2023-06 | 2 | 0 | 2 | 2 | -0.52 | -0.74 | +1.42 | +0.30 | +0.29 | -4.44 | +0.74 | +2.44 |
| 2023-07 | 3 | 0 | 3 | 3 | -0.41 | +0.31 | +0.26 | -0.39 | -1.89 | -6.00 | +0.61 | +3.17 |
| 2023-08 | 2 | 0 | 2 | 2 | -0.34 | -3.29 | +3.89 | -0.90 | +2.59 | +0.00 | +2.34 | +4.29 |
| 2023-09 | 3 | 0 | 3 | 3 | -0.62 | -1.24 | +2.26 | +0.02 | -0.83 | -4.97 | +1.13 | +2.78 |
| 2023-10 | 3 | 0 | 3 | 3 | +0.61 | -0.57 | +1.66 | +0.04 | +0.00 | -1.82 | +0.87 | +1.54 |
| 2023-11 | 2 | 0 | 2 | 2 | -1.90 | -1.41 | +2.19 | -0.77 | +0.66 | -5.79 | +1.98 | +3.18 |
| 2023-12 | 3 | 0 | 3 | 3 | -0.23 | -2.43 | +2.09 | -0.58 | +1.24 | -6.00 | +1.79 | +3.28 |
| 2024-01 | 2 | 0 | 2 | 2 | +0.03 | -2.28 | +2.27 | -0.87 | +0.43 | -6.00 | +1.91 | +3.29 |
| 2024-02 | 3 | 0 | 3 | 3 | -0.05 | -1.87 | +3.09 | -0.29 | +0.67 | -6.00 | +1.36 | +3.38 |
| 2024-03 | 2 | 0 | 2 | 2 | +0.08 | -2.84 | +2.78 | -0.41 | +1.21 | -6.00 | +2.06 | +3.40 |
| 2024-04 | 3 | 0 | 3 | 3 | +0.19 | -1.44 | +2.21 | +0.18 | +0.00 | -6.00 | +1.70 | +3.20 |
| 2024-05 | 2 | 0 | 2 | 2 | +1.43 | -3.32 | +1.41 | -0.09 | +0.56 | -6.00 | +2.56 | +3.18 |
| 2024-06 | 3 | 0 | 3 | 3 | +0.42 | -2.27 | +0.67 | +0.51 | -0.39 | -6.00 | +1.63 | +3.06 |
| 2024-07 | 2 | 0 | 2 | 2 | -1.44 | -2.40 | +0.79 | -0.72 | +1.49 | -6.00 | +2.27 | +3.14 |
| 2024-08 | 3 | 0 | 3 | 3 | +0.33 | -1.70 | +1.60 | -0.49 | -0.39 | -6.00 | +1.33 | +3.13 |
| 2024-09 | 2 | 0 | 2 | 2 | -0.48 | -2.32 | +3.54 | -2.61 | +2.10 | +0.00 | +1.67 | +3.91 |

A block magnitude is computed **per chip** (the RMS of that chip's standardised dimensions, R4) and the panel and the table show the **median of those chip values** — so `R_B` is the typical joint displacement of a single overpass. It is therefore *not* the RMS of the two columns beside it, and it does not vanish when those two medians are near zero: a month whose chips scatter widely in opposite directions has a zero median in each dimension and a large magnitude per chip.

## Frozen reference / held-out pre-event / post-event

Each cell is the median over the chips of that class, in the reference's own scales. The reference column is zero by construction (the standardisation rests on its median); it is printed so the other two columns are read against the frame they are measured in.

| quantity | frozen pre-event reference | held-out pre-event | post-event |
| --- | --- | --- | --- |
| `gamma2` | +0.00 | – | -0.15 |
| `P` | +0.00 | – | -1.60 |
| `R_B` | +0.80 | – | +1.33 |
| `A` | +0.00 | – | +1.85 |
| `D` | +0.00 | – | -0.49 |
| `F` | +0.00 | – | +0.35 |
| `S` | -0.00 | – | -6.00 |
| `R_G` | +0.81 | – | +3.16 |

n = **16** / **0** / **43** chips over **7** / **0** / **17** monthly intervals.

## The contrast (post minus reference, in reference scales)

Circular moving-block bootstrap, **2000** draws, block **3**, seed **7**, one resample per group; the reference side holds **16** chips and the tail **43**.

| quantity | n_ref | n_post | delta z | 95 % CI | verdict |
| --- | --- | --- | --- | --- | --- |
| `gamma2` | 16 | 43 | -0.15 | -0.75 … +0.91 | spans zero |
| `P` | 16 | 43 | -1.60 | -2.36 … -0.53 | **excludes zero** |
| `A` | 16 | 43 | +1.85 | +0.06 … +2.57 | **excludes zero** |
| `D` | 16 | 43 | -0.49 | -0.86 … +0.34 | spans zero |
| `F` | 16 | 43 | +0.35 | -0.73 … +1.36 | spans zero |
| `S` | 16 | 43 | -6.00 | -6.30 … -3.12 | **excludes zero** |
| `R_B` | 16 | 43 | +0.53 | +0.04 … +1.21 | **excludes zero** |
| `R_G` | 16 | 43 | +2.35 | +1.33 … +2.51 | **excludes zero** |

The interval is on the **difference of medians**, not on the curve: a month that moves inside it is not a finding, and a month that leaves it is not yet one either — the tail is the rebuild, so the contrast answers a different question than the collapse plates of the other sites.

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

`--verify` re-derives every number of this report from `data/kdlo/support/kdlo_ocv_channels.csv` and fails if the JSON, the report and the table disagree; the table's sha256 is `9c6bb9197197b38210f543702c425c5c7b1d7facb25c700907e89cdb2dec6381`.

The assessment's pins are constants, asserted by `--selftest`: seed **7**, block **3**, draws **2000**, decision quantiles **q05/q50/q95**, the magnitude boundary **q95**, and a persistence rule of **2** consecutive months. The bootstrap is seeded and no wall-clock value enters any output, so the JSON, the report and the plate reproduce deterministically.
