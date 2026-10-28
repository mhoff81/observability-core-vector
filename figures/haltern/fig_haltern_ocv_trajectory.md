# Haltern AV9 wind turbine — OCV trajectory on the support-sized chip

**Site** Haltern AV9 wind turbine · **event** 2021-09-29 · **reference state** `pre-collapse (healthy)` (pooled over the orbits) · **table** `data/haltern/support/haltern_ocv_channels.csv` (sha256 `cb39c7473b291540`)

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
| frozen pre-event reference (chosen) | `pre-collapse (healthy)` | 2021-01-01 → 2021-09-17 | **96** | 44 | 52 | 96 | 9 |
| tail | `post-collapse` | 2021-10-01 → 2021-10-17 | **11** | 5 | 6 | 11 | 1 |

The reference window is the state's own span (R2, `None`): no sub-window is chosen, so all **96** chips of `pre-collapse (healthy)` (2021-01-01 → 2021-09-17) standardise, and the plate draws nothing outside it.

The cut holds **107** of the table's **107** rows (0 refused, none excluded); the plate draws 10 monthly bins from 2021-01 to 2021-10.

### The three frozen classes

| class | label in this report | chips | ASC | DESC | echo-bearing | monthly intervals |
| --- | --- | --- | --- | --- | --- | --- |
| `reference` | frozen pre-event reference | **96** | 44 | 52 | 96 | 9 |
| `outside` | held-out pre-event | 0 | – | – | – | 0 |
| `post` | post-event | **11** | 5 | 6 | 11 | 1 |

OUTSIDE is strictly before 2021-09-29; POST is the registry's post state(s); a chip on the event day is assigned by its registry state, and a reference-state chip on or after the event is reported as unassigned.

The classes reconcile with the cut: **107** kept = **96** reference + **0** held-out pre-event + **11** post-event + **0** unassigned. The event month is split into a pre and a post bin wherever the record holds chips on both sides, so no monthly aggregation combines the two.

## The reference (R3 / R4)

**pooled**

| dimension | block | standardised against | n | t | median | scale | from | z range |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `gamma2` | behaviour | echo-bearing reference chips | 96 | log | -5.536 | 1.086 | mad | -4.31 … +2.22 |
| `P` | behaviour | echo-bearing reference chips | 96 | log | -2.68 | 0.0637 | mad | -2.06 … +5.83 |
| `A` | geometry | all reference chips | 96 | log | 5.924 | 0.227 | mad | -5.31 … +4.40 |
| `D` | geometry | all reference chips | 96 | log | -3.046 | 0.2381 | mad | -4.89 … +4.68 |
| `F` | geometry | all reference chips | 96 | log | 5.239 | 0.1566 | mad | -6.80 … +1.98 |
| `S` | geometry | all reference chips | 96 | identity | 19.13 | 17.82 | mad | -1.02 … +3.46 |

### Environmental context (panel c1)

Standardised on the same reference, on their own (identity) scale.

| column | group | n | distinct | median | scale | from |
| --- | --- | --- | --- | --- | --- | --- |
| `temperature_c` | pooled | 96 | 79 | 10.15 | 9.044 | mad |
| `wind_speed_ms` | pooled | 96 | 87 | 2.95 | 1.846 | mad |

**3** of the standardised values were clipped at ±6; they are counted here and drawn at the clip, never dropped.

## The reference distribution and the monthly envelope

The decision distribution is built at **monthly** aggregation: each of the **9** reference months contributes one state (the monthly median of its chips), and the envelope is the empirical q05/q50/q95 of those states. A held-out month is compared only with this monthly distribution -- never with the chip-level reference and never with the reference maximum, which is reported below for transparency and is **not** a decision boundary. With 9 monthly states the quantiles are empirical order statistics, not confidence bounds.

| quantity | monthly ref states | q05 | q50 | q95 | min | max | used for |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `gamma2` | 9 | -0.25 | +0.03 | +0.41 | -0.31 | +0.55 | outside [q05, q95] = signed departure |
| `P` | 9 | -0.45 | +0.00 | +0.27 | -0.53 | +0.31 | outside [q05, q95] = signed departure |
| `A` | 9 | -0.28 | +0.12 | +0.75 | -0.32 | +1.00 | outside [q05, q95] = signed departure |
| `D` | 9 | -0.28 | +0.03 | +0.87 | -0.36 | +0.95 | outside [q05, q95] = signed departure |
| `F` | 9 | -0.41 | -0.15 | +0.56 | -0.44 | +0.61 | outside [q05, q95] = signed departure |
| `S` | 9 | -0.40 | -0.01 | +0.26 | -0.41 | +0.28 | outside [q05, q95] = signed departure |
| `R_B` | 9 | +0.80 | +1.05 | +1.48 | +0.80 | +1.49 | above q95 = unusual magnitude |
| `R_G` | 9 | +0.63 | +0.90 | +1.48 | +0.63 | +1.56 | above q95 = unusual magnitude |

The chip-level distributions are stored beside them in the JSON (`reference_envelope.chip`) and are **descriptive only**: nothing in the held-out assessment compares a monthly state with a single chip.

## The monthly trajectory

| month | n | ASC | DESC | echo | z(γ²) | z(P) | z(A) | z(D) | z(F) | z(S) | R_B | R_G |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 2021-01 | 11 | 5 | 6 | 11 | +0.55 | -0.23 | -0.23 | -0.12 | -0.15 | -0.41 | +1.24 | +0.68 |
| 2021-02 | 11 | 5 | 6 | 11 | -0.07 | -0.34 | -0.10 | -0.17 | -0.19 | -0.40 | +0.84 | +0.64 |
| 2021-03 | 11 | 5 | 6 | 11 | -0.31 | +0.11 | -0.32 | -0.36 | -0.44 | +0.01 | +1.49 | +0.63 |
| 2021-04 | 10 | 4 | 6 | 10 | +0.06 | -0.12 | +0.37 | +0.74 | +0.48 | +0.28 | +0.80 | +0.90 |
| 2021-05 | 11 | 5 | 6 | 11 | +0.18 | -0.53 | +0.26 | +0.23 | -0.37 | -0.07 | +1.29 | +0.69 |
| 2021-06 | 11 | 6 | 5 | 11 | +0.15 | +0.31 | +0.12 | +0.03 | -0.22 | +0.04 | +1.47 | +1.56 |
| 2021-07 | 10 | 5 | 5 | 10 | +0.01 | +0.15 | +0.24 | +0.16 | +0.34 | -0.23 | +0.80 | +1.37 |
| 2021-08 | 11 | 5 | 6 | 11 | -0.16 | +0.22 | -0.10 | -0.11 | +0.25 | -0.01 | +1.05 | +1.02 |
| 2021-09 | 10 | 4 | 6 | 10 | +0.03 | +0.00 | +1.00 | +0.95 | +0.61 | +0.24 | +0.82 | +1.01 |
| 2021-10 | 11 | 5 | 6 | 11 | -0.27 | +0.24 | +0.64 | +0.55 | +0.63 | -0.15 | +0.88 | +1.13 |

A block magnitude is computed **per chip** (the RMS of that chip's standardised dimensions, R4) and the panel and the table show the **median of those chip values** — so `R_B` is the typical joint displacement of a single overpass. It is therefore *not* the RMS of the two columns beside it, and it does not vanish when those two medians are near zero: a month whose chips scatter widely in opposite directions has a zero median in each dimension and a large magnitude per chip.

## Frozen reference / held-out pre-event / post-event

Each cell is the median over the chips of that class, in the reference's own scales. The reference column is zero by construction (the standardisation rests on its median); it is printed so the other two columns are read against the frame they are measured in.

| quantity | frozen pre-event reference | held-out pre-event | post-event |
| --- | --- | --- | --- |
| `gamma2` | +0.00 | – | -0.27 |
| `P` | -0.00 | – | +0.24 |
| `R_B` | +1.02 | – | +0.88 |
| `A` | +0.00 | – | +0.64 |
| `D` | +0.00 | – | +0.55 |
| `F` | +0.00 | – | +0.63 |
| `S` | +0.00 | – | -0.15 |
| `R_G` | +0.83 | – | +1.13 |

n = **96** / **0** / **11** chips over **9** / **0** / **1** monthly intervals.

## The contrast (post minus reference, in reference scales)

Circular moving-block bootstrap, **2000** draws, block **3**, seed **7**, one resample per group; the reference side holds **96** chips and the tail **11**.

| quantity | n_ref | n_post | delta z | 95 % CI | verdict |
| --- | --- | --- | --- | --- | --- |
| `gamma2` | 96 | 11 | -0.27 | -0.56 … +0.96 | spans zero |
| `P` | 96 | 11 | +0.24 | -0.27 … +0.55 | spans zero |
| `A` | 96 | 11 | +0.64 | -0.44 … +1.96 | spans zero |
| `D` | 96 | 11 | +0.55 | -0.52 … +1.79 | spans zero |
| `F` | 96 | 11 | +0.63 | +0.00 … +1.19 | spans zero |
| `S` | 96 | 11 | -0.15 | -0.42 … +1.13 | spans zero |
| `R_B` | 96 | 11 | -0.14 | -0.55 … +0.31 | spans zero |
| `R_G` | 96 | 11 | +0.30 | -0.35 … +0.86 | spans zero |

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

`--verify` re-derives every number of this report from `data/haltern/support/haltern_ocv_channels.csv` and fails if the JSON, the report and the table disagree; the table's sha256 is `cb39c7473b291540432f65f2428a61b31290bbcf3c8e8021044ff0c0e3b36ee1`.

The assessment's pins are constants, asserted by `--selftest`: seed **7**, block **3**, draws **2000**, decision quantiles **q05/q50/q95**, the magnitude boundary **q95**, and a persistence rule of **2** consecutive months. The bootstrap is seeded and no wall-clock value enters any output, so the JSON, the report and the plate reproduce deterministically.
