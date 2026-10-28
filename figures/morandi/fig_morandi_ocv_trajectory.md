# Ponte Morandi, Genova — OCV trajectory on the support-sized chip

**Site** Ponte Morandi, Genova · **event** 2018-08-14 · **reference state** `pre-collapse (healthy)` (pooled over the orbits) · **table** `data/morandi/support/morandi_ocv_channels.csv` (sha256 `249c9ee6630c21df`)

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
| frozen pre-event reference (chosen) | `pre-collapse (healthy)` | 2016-08-24 → 2018-08-08 | **198** | 91 | 107 | 198 | 25 |
| tail | `post-collapse` | 2018-08-14 → 2018-09-13 | **9** | 3 | 6 | 9 | 2 |

The reference window is the state's own span (R2, `None`): no sub-window is chosen, so all **198** chips of `pre-collapse (healthy)` (2016-08-24 → 2018-08-08) standardise, and the plate draws nothing outside it.

The cut holds **207** of the table's **207** rows (0 refused, none excluded); the plate draws 27 monthly bins from 2016-08 to 2018-09.

### The three frozen classes

| class | label in this report | chips | ASC | DESC | echo-bearing | monthly intervals |
| --- | --- | --- | --- | --- | --- | --- |
| `reference` | frozen pre-event reference | **198** | 91 | 107 | 198 | 25 |
| `outside` | held-out pre-event | 0 | – | – | – | 0 |
| `post` | post-event | **9** | 3 | 6 | 9 | 2 |

OUTSIDE is strictly before 2018-08-14; POST is the registry's post state(s); a chip on the event day is assigned by its registry state, and a reference-state chip on or after the event is reported as unassigned.

The classes reconcile with the cut: **207** kept = **198** reference + **0** held-out pre-event + **9** post-event + **0** unassigned. The event month is split into a pre and a post bin wherever the record holds chips on both sides, so no monthly aggregation combines the two.

## The reference (R3 / R4)

**pooled**

| dimension | block | standardised against | n | t | median | scale | from | z range |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `gamma2` | behaviour | echo-bearing reference chips | 198 | log | -3.9 | 2.491 | mad | -2.94 … +1.17 |
| `P` | behaviour | echo-bearing reference chips | 198 | log | -2.366 | 0.7058 | mad | -2.35 … +1.83 |
| `A` | geometry | all reference chips | 198 | log | 3.238 | 1.432 | mad | -1.49 … +2.38 |
| `D` | geometry | all reference chips | 198 | log | -3.48 | 1.77 | mad | -1.25 … +1.97 |
| `F` | geometry | all reference chips | 198 | log | 2.518 | 2.066 | mad | -1.22 … +1.40 |
| `S` | geometry | all reference chips | 198 | identity | 6.616 | 6.541 | mad | -1.01 … +5.25 |

### Environmental context (panel c1)

Standardised on the same reference, on their own (identity) scale.

| column | group | n | distinct | median | scale | from |
| --- | --- | --- | --- | --- | --- | --- |
| `temperature_c` | pooled | 198 | 133 | 14.7 | 8.673 | mad |
| `wind_speed_ms` | pooled | 198 | 155 | 2.42 | 1.994 | mad |

## The reference distribution and the monthly envelope

The decision distribution is built at **monthly** aggregation: each of the **25** reference months contributes one state (the monthly median of its chips), and the envelope is the empirical q05/q50/q95 of those states. A held-out month is compared only with this monthly distribution -- never with the chip-level reference and never with the reference maximum, which is reported below for transparency and is **not** a decision boundary. With 25 monthly states the quantiles are empirical order statistics, not confidence bounds.

| quantity | monthly ref states | q05 | q50 | q95 | min | max | used for |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `gamma2` | 25 | -0.88 | -0.04 | +0.50 | -0.99 | +0.54 | outside [q05, q95] = signed departure |
| `P` | 25 | -0.63 | +0.20 | +0.48 | -1.19 | +0.50 | outside [q05, q95] = signed departure |
| `A` | 25 | -0.52 | +0.25 | +1.66 | -0.59 | +1.84 | outside [q05, q95] = signed departure |
| `D` | 25 | -0.60 | +0.10 | +0.85 | -0.63 | +1.32 | outside [q05, q95] = signed departure |
| `F` | 25 | -0.56 | +0.06 | +0.98 | -0.88 | +1.12 | outside [q05, q95] = signed departure |
| `S` | 25 | -0.64 | +0.01 | +0.85 | -0.81 | +0.86 | outside [q05, q95] = signed departure |
| `R_B` | 25 | +0.55 | +0.79 | +0.98 | +0.45 | +1.10 | above q95 = unusual magnitude |
| `R_G` | 25 | +0.61 | +0.87 | +1.21 | +0.57 | +1.23 | above q95 = unusual magnitude |

The chip-level distributions are stored beside them in the JSON (`reference_envelope.chip`) and are **descriptive only**: nothing in the held-out assessment compares a monthly state with a single chip.

## The monthly trajectory

| month | n | ASC | DESC | echo | z(γ²) | z(P) | z(A) | z(D) | z(F) | z(S) | R_B | R_G |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 2016-08 | 1 | 0 | 1 | 1 | -0.41 | -1.19 | +0.90 | -0.45 | +0.58 | +0.86 | +0.89 | +0.72 |
| 2016-09 | 6 | 3 | 3 | 6 | +0.36 | -0.06 | +0.06 | -0.63 | -0.03 | +0.81 | +0.71 | +0.67 |
| 2016-10 | 10 | 5 | 5 | 10 | +0.47 | -0.63 | -0.59 | +0.30 | -0.56 | -0.36 | +0.76 | +0.94 |
| 2016-11 | 9 | 5 | 4 | 9 | +0.54 | -0.01 | -0.53 | +0.89 | -0.55 | -0.50 | +0.93 | +0.87 |
| 2016-12 | 10 | 5 | 5 | 10 | +0.11 | -0.08 | -0.26 | +0.58 | -0.30 | -0.48 | +0.89 | +1.02 |
| 2017-01 | 9 | 5 | 4 | 9 | +0.32 | +0.21 | -0.47 | +1.32 | -0.88 | -0.81 | +0.45 | +1.00 |
| 2017-02 | 7 | 5 | 2 | 7 | +0.40 | +0.32 | -0.47 | +0.00 | -0.55 | -0.39 | +0.80 | +0.68 |
| 2017-03 | 9 | 5 | 4 | 9 | +0.51 | +0.50 | -0.42 | +0.54 | -0.55 | -0.68 | +0.71 | +0.74 |
| 2017-04 | 10 | 5 | 5 | 10 | +0.24 | +0.25 | -0.12 | -0.58 | -0.14 | +0.37 | +0.67 | +1.02 |
| 2017-05 | 9 | 5 | 4 | 9 | +0.25 | +0.41 | -0.01 | -0.57 | -0.35 | -0.21 | +0.79 | +0.60 |
| 2017-06 | 10 | 5 | 5 | 10 | +0.13 | +0.05 | -0.16 | -0.60 | -0.16 | +0.16 | +0.79 | +0.67 |
| 2017-07 | 9 | 5 | 4 | 9 | +0.00 | +0.49 | -0.04 | -0.59 | -0.16 | +0.38 | +0.73 | +0.57 |
| 2017-08 | 11 | 5 | 6 | 11 | -0.04 | -0.14 | +0.64 | -0.52 | +0.50 | +0.19 | +0.53 | +0.63 |
| 2017-09 | 10 | 5 | 5 | 10 | -0.07 | +0.30 | +0.02 | -0.51 | -0.10 | -0.07 | +0.72 | +0.67 |
| 2017-10 | 10 | 5 | 5 | 10 | +0.10 | +0.28 | +0.25 | +0.71 | +0.10 | -0.02 | +0.97 | +1.22 |
| 2017-11 | 9 | 4 | 5 | 9 | -0.47 | +0.22 | +0.46 | +0.10 | +0.06 | -0.09 | +0.73 | +1.17 |
| 2017-12 | 8 | 3 | 5 | 8 | -0.72 | +0.26 | +1.28 | +0.46 | +0.70 | -0.39 | +0.89 | +0.85 |
| 2018-01 | 7 | 2 | 5 | 7 | -0.71 | -0.25 | +1.84 | +0.61 | +1.12 | +0.01 | +0.87 | +1.23 |
| 2018-02 | 7 | 2 | 5 | 7 | -0.89 | -0.09 | +1.75 | +0.33 | +1.05 | -0.11 | +0.69 | +1.12 |
| 2018-03 | 8 | 3 | 5 | 8 | -0.99 | +0.20 | +1.12 | +0.47 | +0.57 | +0.13 | +0.98 | +0.93 |
| 2018-04 | 6 | 2 | 4 | 6 | -0.84 | +0.46 | +1.13 | +0.24 | +0.62 | +0.03 | +0.86 | +1.04 |
| 2018-05 | 7 | 2 | 5 | 7 | -0.46 | +0.27 | +0.43 | -0.30 | +0.30 | +0.11 | +0.63 | +0.90 |
| 2018-06 | 6 | 2 | 4 | 6 | -0.23 | -0.24 | +0.97 | -0.05 | +0.63 | +0.29 | +0.69 | +1.08 |
| 2018-07 | 7 | 2 | 5 | 7 | -0.41 | -0.05 | +0.47 | -0.45 | +0.36 | +0.50 | +0.86 | +0.75 |
| 2018-08 (pre) | 3 | 1 | 2 | 3 | -0.67 | -0.61 | +0.54 | +0.72 | +0.28 | +0.86 | +1.10 | +0.75 |
| 2018-08 (post) | 5 | 2 | 3 | 5 | -0.00 | +0.69 | +0.56 | -0.46 | +0.30 | +0.27 | +0.75 | +1.02 |
| 2018-09 | 4 | 1 | 3 | 4 | +0.08 | -0.19 | +0.26 | -0.88 | +0.21 | +0.92 | +0.52 | +1.03 |

A block magnitude is computed **per chip** (the RMS of that chip's standardised dimensions, R4) and the panel and the table show the **median of those chip values** — so `R_B` is the typical joint displacement of a single overpass. It is therefore *not* the RMS of the two columns beside it, and it does not vanish when those two medians are near zero: a month whose chips scatter widely in opposite directions has a zero median in each dimension and a large magnitude per chip.

## Frozen reference / held-out pre-event / post-event

Each cell is the median over the chips of that class, in the reference's own scales. The reference column is zero by construction (the standardisation rests on its median); it is printed so the other two columns are read against the frame they are measured in.

| quantity | frozen pre-event reference | held-out pre-event | post-event |
| --- | --- | --- | --- |
| `gamma2` | +0.00 | – | +0.06 |
| `P` | +0.00 | – | +0.08 |
| `R_B` | +0.78 | – | +0.63 |
| `A` | +0.00 | – | +0.56 |
| `D` | +0.00 | – | -0.53 |
| `F` | -0.00 | – | +0.30 |
| `S` | +0.00 | – | +0.27 |
| `R_G` | +0.91 | – | +1.03 |

n = **198** / **0** / **9** chips over **25** / **0** / **2** monthly intervals.

## The contrast (post minus reference, in reference scales)

Circular moving-block bootstrap, **2000** draws, block **3**, seed **7**, one resample per group; the reference side holds **198** chips and the tail **9**.

| quantity | n_ref | n_post | delta z | 95 % CI | verdict |
| --- | --- | --- | --- | --- | --- |
| `gamma2` | 198 | 9 | +0.06 | -1.06 … +0.30 | spans zero |
| `P` | 198 | 9 | +0.08 | -0.24 … +0.74 | spans zero |
| `A` | 198 | 9 | +0.56 | -0.45 … +0.97 | spans zero |
| `D` | 198 | 9 | -0.53 | -0.95 … -0.12 | **excludes zero** |
| `F` | 198 | 9 | +0.30 | -0.18 … +0.70 | spans zero |
| `S` | 198 | 9 | +0.27 | -0.04 … +2.21 | spans zero |
| `R_B` | 198 | 9 | -0.15 | -0.31 … +0.02 | spans zero |
| `R_G` | 198 | 9 | +0.12 | +0.03 … +0.26 | **excludes zero** |

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

`--verify` re-derives every number of this report from `data/morandi/support/morandi_ocv_channels.csv` and fails if the JSON, the report and the table disagree; the table's sha256 is `249c9ee6630c21dfbe3e01366c9863bf1693bc4376d08e93da7b913a0b7a2745`.

The assessment's pins are constants, asserted by `--selftest`: seed **7**, block **3**, draws **2000**, decision quantiles **q05/q50/q95**, the magnitude boundary **q95**, and a persistence rule of **2** consecutive months. The bootstrap is seeded and no wall-clock value enters any output, so the JSON, the report and the plate reproduce deterministically.
