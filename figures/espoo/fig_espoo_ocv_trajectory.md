# Espoo Kurttila mast (event-free control) — OCV trajectory on the support-sized chip

**Site** Espoo Kurttila mast (event-free control) · **event** none · **reference state** `all` (per orbit) · **table** `data/espoo/support/espoo_ocv_channels.csv` (sha256 `82a8cc04208e5105`)

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
| frozen pre-event reference (chosen) | `all` | 2024-09-05 → 2026-09-07 | **75** | 31 | 44 | 75 | 25 |

The reference window is the state's own span (R2, `None`): no sub-window is chosen, so all **75** chips of `all` (2024-09-05 → 2026-09-07) standardise, and the plate draws nothing outside it.

The cut holds **75** of the table's **75** rows (0 refused, none excluded); the plate draws 25 monthly bins from 2024-09 to 2026-09.

### The three frozen classes

| class | label in this report | chips | ASC | DESC | echo-bearing | monthly intervals |
| --- | --- | --- | --- | --- | --- | --- |
| `reference` | frozen pre-event reference | **75** | 31 | 44 | 75 | 25 |
| `outside` | held-out pre-event | 0 | – | – | – | 0 |
| `post` | post-event | 0 | – | – | – | 0 |

OUTSIDE is strictly before the (absent) event; POST is the registry's post state(s); a chip on the event day is assigned by its registry state, and a reference-state chip on or after the event is reported as unassigned.

The classes reconcile with the cut: **75** kept = **75** reference + **0** held-out pre-event + **0** post-event + **0** unassigned. The event month is split into a pre and a post bin wherever the record holds chips on both sides, so no monthly aggregation combines the two.

## The reference (R3 / R4)

**ASCENDING** — per orbit, never pooled

| dimension | block | standardised against | n | t | median | scale | from | z range |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `gamma2` | behaviour | echo-bearing reference chips | 31 | log | -4.111 | 1.233 | mad | -4.08 … +1.11 |
| `P` | behaviour | echo-bearing reference chips | 31 | log | -1.453 | 0.08149 | mad | -5.75 … +2.11 |
| `A` | geometry | all reference chips | 31 | log | 4.331 | 0.2781 | mad | -2.05 … +3.65 |
| `D` | geometry | all reference chips | 31 | log | -3.863 | 0.1635 | mad | -2.82 … +5.26 |
| `F` | geometry | all reference chips | 31 | log | 2.944 | 0.2174 | mad | -1.40 … +3.65 |
| `S` | geometry | all reference chips | 31 | identity | 25.72 | 3.546 | mad | -4.93 … +1.74 |

**DESCENDING** — per orbit, never pooled

| dimension | block | standardised against | n | t | median | scale | from | z range |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `gamma2` | behaviour | echo-bearing reference chips | 44 | log | -7.26 | 1.218 | mad | -4.21 … +1.68 |
| `P` | behaviour | echo-bearing reference chips | 44 | log | -1.916 | 0.07715 | mad | -3.27 … +3.87 |
| `A` | geometry | all reference chips | 44 | log | 6.268 | 0.3154 | mad | -3.65 … +1.27 |
| `D` | geometry | all reference chips | 44 | log | -2.147 | 0.3154 | mad | -3.65 … +1.27 |
| `F` | geometry | all reference chips | 44 | log | 4.89 | 0.1117 | mad | -4.66 … +1.60 |
| `S` | geometry | all reference chips | 44 | identity | 14.39 | 8.674 | mad | -1.41 … +1.26 |

### Environmental context (panel c1)

Standardised on the same reference, on their own (identity) scale.

| column | group | n | distinct | median | scale | from |
| --- | --- | --- | --- | --- | --- | --- |
| `temperature_c` | ASCENDING | 31 | 5 | 10 | 2.086 | sd |
| `temperature_c` | DESCENDING | 44 | 7 | 10 | 2.765 | sd |
| `wind_speed_ms` | ASCENDING | 31 | 5 | 5 | 0.5763 | sd |
| `wind_speed_ms` | DESCENDING | 44 | 7 | 5 | 0.9466 | sd |

A co-variate that carries only a handful of distinct values is a climatology, not a per-acquisition measurement — `temperature_c` (5 values over up to 44 chips); `wind_speed_ms` (5 values over up to 44 chips) — so panel c1 steps between those values and cannot resolve anything finer. It is drawn for exactly that reason: to keep the reader from reading the other panels' wiggles as weather.

## The reference distribution and the monthly envelope

The decision distribution is built at **monthly** aggregation: each of the **25** reference months contributes one state (the monthly median of its chips), and the envelope is the empirical q05/q50/q95 of those states. A held-out month is compared only with this monthly distribution -- never with the chip-level reference and never with the reference maximum, which is reported below for transparency and is **not** a decision boundary. With 25 monthly states the quantiles are empirical order statistics, not confidence bounds.

| quantity | monthly ref states | q05 | q50 | q95 | min | max | used for |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `gamma2` | 25 | -1.26 | +0.07 | +0.66 | -1.93 | +0.80 | outside [q05, q95] = signed departure |
| `P` | 25 | -1.72 | +0.00 | +1.75 | -2.33 | +1.94 | outside [q05, q95] = signed departure |
| `A` | 25 | -1.26 | +0.06 | +0.99 | -1.62 | +1.16 | outside [q05, q95] = signed departure |
| `D` | 25 | -1.36 | +0.00 | +0.91 | -1.62 | +1.00 | outside [q05, q95] = signed departure |
| `F` | 25 | -1.04 | +0.00 | +0.84 | -2.03 | +1.25 | outside [q05, q95] = signed departure |
| `S` | 25 | -0.67 | +0.13 | +0.86 | -1.09 | +1.08 | outside [q05, q95] = signed departure |
| `R_B` | 25 | +0.44 | +0.90 | +1.84 | +0.23 | +1.97 | above q95 = unusual magnitude |
| `R_G` | 25 | +0.55 | +0.83 | +2.25 | +0.43 | +2.91 | above q95 = unusual magnitude |

The chip-level distributions are stored beside them in the JSON (`reference_envelope.chip`) and are **descriptive only**: nothing in the held-out assessment compares a monthly state with a single chip.

## The monthly trajectory

| month | n | ASC | DESC | echo | z(γ²) | z(P) | z(A) | z(D) | z(F) | z(S) | R_B | R_G |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 2024-09 | 3 | 2 | 1 | 3 | +0.80 | +1.81 | -0.79 | -0.94 | -0.25 | -0.56 | +1.32 | +0.72 |
| 2024-10 | 3 | 1 | 2 | 3 | +0.26 | -0.34 | +0.93 | +0.93 | +0.46 | -0.64 | +1.85 | +0.91 |
| 2024-11 | 3 | 2 | 1 | 3 | -0.30 | +0.24 | +0.14 | +0.33 | +0.88 | +0.58 | +1.10 | +0.56 |
| 2024-12 | 3 | 1 | 2 | 3 | +0.29 | +0.00 | +0.38 | +0.38 | +1.25 | +0.73 | +0.68 | +0.81 |
| 2025-01 | 3 | 2 | 1 | 3 | -0.71 | -0.00 | -0.10 | -0.21 | +0.46 | -0.39 | +0.78 | +0.71 |
| 2025-02 | 3 | 1 | 2 | 3 | -0.05 | -0.45 | +0.29 | +0.29 | -0.14 | +0.78 | +0.90 | +0.68 |
| 2025-03 | 3 | 2 | 1 | 3 | +0.39 | +1.51 | -0.93 | -1.37 | -0.79 | +1.08 | +1.10 | +1.10 |
| 2025-04 | 3 | 1 | 2 | 3 | +0.37 | -1.83 | +1.16 | +0.82 | +0.65 | -0.67 | +1.54 | +0.95 |
| 2025-05 | 3 | 1 | 2 | 3 | +0.47 | +1.94 | -0.73 | -1.19 | -1.09 | +0.13 | +1.41 | +1.16 |
| 2025-06 | 3 | 1 | 2 | 3 | +0.48 | -0.63 | +0.18 | +0.18 | +0.59 | -0.10 | +0.77 | +1.06 |
| 2025-07 | 3 | 1 | 2 | 3 | +0.53 | -0.10 | -0.85 | -0.67 | -0.51 | +0.81 | +0.67 | +2.53 |
| 2025-08 | 3 | 2 | 1 | 3 | -0.54 | +0.74 | +0.06 | +0.59 | +0.40 | -0.56 | +0.65 | +0.83 |
| 2025-09 | 3 | 1 | 2 | 3 | -0.70 | +0.35 | +0.27 | +0.40 | +0.00 | +0.42 | +0.90 | +0.61 |
| 2025-10 | 3 | 1 | 2 | 3 | +0.52 | -0.01 | -0.44 | -0.44 | +0.20 | +0.23 | +0.48 | +0.67 |
| 2025-11 | 3 | 1 | 2 | 3 | +0.07 | +0.24 | +0.70 | +0.70 | +0.00 | +0.47 | +0.44 | +0.60 |
| 2025-12 | 3 | 2 | 1 | 3 | -0.23 | +0.00 | +0.44 | +0.00 | +0.24 | +0.23 | +0.23 | +0.56 |
| 2026-01 | 4 | 1 | 3 | 4 | -1.93 | -1.26 | -0.09 | -0.09 | +0.13 | +0.21 | +1.83 | +0.43 |
| 2026-02 | 3 | 1 | 2 | 3 | -0.80 | +0.06 | -1.62 | -1.62 | -2.03 | -1.09 | +1.97 | +2.91 |
| 2026-03 | 3 | 1 | 2 | 3 | -0.35 | +1.13 | -1.34 | -1.34 | -0.27 | -0.13 | +1.02 | +1.16 |
| 2026-04 | 3 | 2 | 1 | 3 | +0.07 | -2.33 | -0.79 | -0.79 | -0.85 | -0.43 | +1.65 | +0.90 |
| 2026-05 | 3 | 1 | 2 | 3 | +0.14 | -0.05 | -0.06 | -0.06 | +0.67 | +0.87 | +0.82 | +0.92 |
| 2026-06 | 3 | 1 | 2 | 3 | -0.58 | -0.68 | +0.71 | +0.71 | +0.00 | -0.53 | +0.73 | +0.74 |
| 2026-07 | 2 | 0 | 2 | 2 | -1.37 | +0.16 | +0.49 | +0.49 | -0.24 | -0.25 | +1.00 | +0.55 |
| 2026-08 | 3 | 1 | 2 | 3 | +0.69 | +0.13 | -0.35 | -0.67 | -0.20 | -0.02 | +0.73 | +0.85 |
| 2026-09 | 3 | 1 | 2 | 3 | -0.41 | -0.21 | +1.00 | +1.00 | -0.25 | +0.16 | +0.69 | +0.83 |

A block magnitude is computed **per chip** (the RMS of that chip's standardised dimensions, R4) and the panel and the table show the **median of those chip values** — so `R_B` is the typical joint displacement of a single overpass. It is therefore *not* the RMS of the two columns beside it, and it does not vanish when those two medians are near zero: a month whose chips scatter widely in opposite directions has a zero median in each dimension and a large magnitude per chip.

## Frozen reference / held-out pre-event / post-event

Each cell is the median over the chips of that class, in the reference's own scales. The reference column is zero by construction (the standardisation rests on its median); it is printed so the other two columns are read against the frame they are measured in.

| quantity | frozen pre-event reference | held-out pre-event | post-event |
| --- | --- | --- | --- |
| `gamma2` | +0.00 | – | – |
| `P` | +0.00 | – | – |
| `R_B` | +0.90 | – | – |
| `A` | +0.00 | – | – |
| `D` | +0.00 | – | – |
| `F` | +0.00 | – | – |
| `S` | +0.00 | – | – |
| `R_G` | +0.80 | – | – |

n = **75** / **0** / **0** chips over **25** / **0** / **0** monthly intervals.

## The contrast (post minus reference, in reference scales)

the record labels no event: there is no post side to contrast the reference with, and the plate draws the trajectory alone.

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
* **No contrast at all on this site.** The event-free control has no tail, so the plate is the trajectory alone — which is the point of a control.
* **No event line.** There is no event in this record to draw, so the plate carries none; the reference window is the whole record.

## Pins

`--verify` re-derives every number of this report from `data/espoo/support/espoo_ocv_channels.csv` and fails if the JSON, the report and the table disagree; the table's sha256 is `82a8cc04208e5105d04c3780ea418f6c7be80d04b8cf0bc0aad5cb87dd0dd1f5`.

The assessment's pins are constants, asserted by `--selftest`: seed **7**, block **3**, draws **2000**, decision quantiles **q05/q50/q95**, the magnitude boundary **q95**, and a persistence rule of **2** consecutive months. The bootstrap is seeded and no wall-clock value enters any output, so the JSON, the report and the plate reproduce deterministically.
