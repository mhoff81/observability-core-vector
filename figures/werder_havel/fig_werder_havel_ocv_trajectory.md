# Havelbrücke Werder, A 10 — OCV trajectory on the support-sized chip

**Asset** `havelbruecke_werder_a10` · **documented interval** 2024-07-01 → 2025-09-05 (431 d) · **reference state** `W2_bestand` (pooled over the orbits) · **table** `data/werder_havel/support/werder_havel_ocv_channels.csv` (sha256 `5915a80fdf9841a2`) · **plate** `figures/werder_havel/fig_werder_havel_ocv_trajectory.jpg`

```
standardisation  z = (t(x) - centre) / scale; centre = median_ref, scale = 1.4826 * MAD_ref, |z| <= 6
                 t = log for gamma2, P, A, D, F; identity for S (R4)
blocks           behaviour = ['gamma2', 'P']   geometry = ['A', 'D', 'F', 'S']   (ocv_reconfigure.BLOCKS)
reference per block  behaviour: echo-bearing reference chips; geometry: every reference chip (R3)
bootstrap        2000 circular moving-block draws, block 3, seed 7 (R5)
bins             calendar months, split at a documented boundary where the record holds chips on both sides
```

The crossing's three strata are the export's own (`W2_bestand` / `W3_sanierung` / `W4_post`), cut on the two documented dates of the 2024/2025 renewal.

## The documented interval

The interval **2024-07-01 → 2025-09-05** (431 d) is read from `data/werder_havel/reference/werder_havel_aoi.json` and never typed into the plate. Its stated kind is *documented change (Autobahn GmbH, A 10 renewal)*, and it does **not** come from the map: none of the way's five versions names the 2024/2025 renewal: v4/v5 (2026) add the structural tags, and the `description` names the *construction* of the two superstructures (1997 / 2004 — the 'Richtungsfahrbahn Nauen' and 'Richtungsfahrbahn Werder' decks). The event is the external chain's, corroborated by the record's own strata only (W3_sanierung / W4_post are cut on exactly these two dates)

| source | role | outcome |
| --- | --- | --- |
| `data/werder_havel/reference/werder_havel_aoi.json` | committed reference (read) | read -- the interval is taken from here |
| `data/werder_havel/support/werder_havel_ocv_channels_meta.json` | generator's own provenance | agrees with the block above |
| `/home/projects/havel_werder_analysis/havel_werder_events.json` | operator record (agrees) | agrees on both dates -- cross-check passed |

The external operator record is not taken on trust: it is read as data (`/home/projects/havel_werder_analysis/havel_werder_events.json`, sha256 `cad5d6aa9c74a870`), its two entries are matched by label to the two documented dates, and any disagreement is a failure of the run, not a note. It names `sanierung_beginn` 2024-07-01 and `sanierung_ende` 2025-09-05 -- the same two dates, to the day.

The metadata is also checked against the data it labels, not only against another document: the export's strata must be cut on these dates, within a tolerance of **7 d** -- the record's own longest acquisition gap, because a boundary can only be located to within one gap. The observed gaps are W2_bestand -> W3_sanierung 4 d, W3_sanierung -> W4_post 1 d.

| state | first chip | last chip | chips |
| --- | --- | --- | --- |
| `W2_bestand` | 2022-07-02 | 2024-06-29 | 241 |
| `W3_sanierung` | 2024-07-03 | 2025-09-04 | 189 |
| `W4_post` | 2025-09-05 | 2025-12-04 | 61 |

The export's strata are cut on the documented dates: W2_bestand 2022-07-02 -> 2024-06-29, W3_sanierung 2024-07-03 -> 2025-09-04, W4_post 2025-09-05 -> 2025-12-04 (W2_bestand -> W3_sanierung 4 d, W3_sanierung -> W4_post 1 d).

## The record and its cut

The support cut holds **491** of the table's **991** rows (0 refused, 500 not a target chip), over **491** distinct acquisition dates from 2022-07-02 to 2025-12-04. The plate draws **43** monthly bins.

| class | label in this report | state | chips | ASC | DESC | echo-bearing | months |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `reference` | frozen pre-intervention reference | `W2_bestand` | **119** | 61 | 58 | 119 | 2022-07 → 2023-06 (12) |
| `outside` | held-out pre-intervention | `W2_bestand` | **122** | 61 | 61 | 122 | 2023-07 → 2024-06 (12) |
| `during` | intervention period | `W3_sanierung` | **189** | 91 | 98 | 189 | 2024-07 → 2025-09 (15) |
| `post` | post-intervention | `W4_post` | **61** | 31 | 30 | 61 | 2025-09 → 2025-12 (4) |

The classes are disjoint by construction -- a chip is placed by the export's state window and, for the reference state only, by the month its acquisition falls in -- so **491** kept = **119** reference + **122** held-out pre-intervention + **189** intervention + **61** post, with **0** unassigned.

The frozen reference is **2022-07 → 2023-06** (inclusive calendar months, 12 bins, 119 chips): the record's first full annual cycle, fixed before any result of this plate was read. The remaining 122 chips of the same state (2023-07 → 2024-06) are drawn hollow and grey as the held-out pre-intervention year; they standardise nothing.

A month that holds chips on both sides of a documented boundary is split, so no bin mixes two sides of the interval: `2025-09` split at 2025-09-05 (end). A month is otherwise one bin, and no line is drawn across a boundary.

## The reference (R3 / R4)

**pooled over the orbits** — 119 chips, 2022-07 → 2023-06. Each channel is standardised on the reference alone.

| dimension | block | standardised against | n | t | centre | scale | from | z range (before the clip) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `gamma2` | behaviour | echo-bearing reference chips | 119 | log | -2.726 | 1.357 | mad | -5.00 … +1.95 |
| `P` | behaviour | echo-bearing reference chips | 119 | log | -3.207 | 0.5157 | mad | -5.80 … +1.68 |
| `A` | geometry | all reference chips | 119 | log | 2.485 | 0.6011 | mad | -2.31 … +2.20 |
| `D` | geometry | all reference chips | 119 | log | -4.679 | 0.6602 | mad | -1.90 … +7.09 |
| `F` | geometry | all reference chips | 119 | log | 1.609 | 0.4989 | mad | -3.23 … +2.33 |
| `S` | geometry | all reference chips | 119 | identity | 14.25 | 9.37 | mad | -1.45 … +2.44 |

A scale of zero is refused by `scale_of`, so no channel is standardised by a spread that is not its own. The geometry block is standardised on every reference chip; the behaviour block only on the reference chips that carry an echo, because on a chip without one a squared coherence and a scatterer count have no population to be drawn from (R3).

**27** value(s) lie beyond |z| = 6 and are drawn at the clip, never dropped; they are listed in `standardisation.clipping`, and the range in the table above is the raw one, before the clip.

### Environmental context (panel c1)

Standardised on the same reference chips, on their own (identity) scale, and reported beside the OCV numbers -- never inside them.

| column | n | centre | scale | from |
| --- | --- | --- | --- | --- |
| `temperature_c` | 119 | 10.5 | 7.71 | mad |
| `precipitation_mm` | 119 | 0.4 | 0.593 | mad |

## The reference distribution and the monthly envelope (R5)

The decision distribution is built at **monthly** aggregation: each of the **12** reference months contributes one state (the median of its chips), and the envelope is the empirical q05/q50/q95 of those states. A held-out month is compared only with this monthly distribution -- never with a single reference chip and never with the reference maximum, which is printed below for transparency and is **not** a decision boundary. With 12 monthly states the quantiles are empirical order statistics, not confidence bounds.

| quantity | monthly reference states | q05 | q50 | q95 | min | max | read as |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `gamma2` | 12 | -0.71 | -0.04 | +0.72 | -1.26 | +1.15 | outside [q05, q95] = signed departure |
| `P` | 12 | -0.78 | +0.08 | +0.48 | -0.80 | +0.50 | outside [q05, q95] = signed departure |
| `A` | 12 | -0.77 | +0.16 | +0.56 | -0.90 | +0.66 | outside [q05, q95] = signed departure |
| `D` | 12 | -0.43 | +0.02 | +0.21 | -0.50 | +0.21 | outside [q05, q95] = signed departure |
| `F` | 12 | -0.45 | +0.09 | +0.79 | -0.45 | +0.94 | outside [q05, q95] = signed departure |
| `S` | 12 | -0.60 | -0.03 | +0.64 | -0.71 | +0.86 | outside [q05, q95] = signed departure |
| `R_B` | 12 | +0.73 | +0.81 | +1.16 | +0.72 | +1.34 | above q95 = unusual magnitude |
| `R_G` | 12 | +0.67 | +0.84 | +1.07 | +0.65 | +1.21 | above q95 = unusual magnitude |

The chip-level distributions are stored beside them in the JSON (`envelope.chip`) and are **descriptive only**: nothing in the held-out reading compares a monthly state with a single chip, and the chip band is the light band behind the plate's panels only.

## The monthly trajectory

| bin | class | n | ASC | DESC | echo | z(γ²) | z(P) | z(A) | z(D) | z(F) | z(S) | R_B | R_G |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 2022-07 | reference | 11 | 5 | 6 | 11 | +0.12 | +0.25 | +0.13 | -0.37 | +0.00 | -0.23 | +0.82 | +0.76 |
| 2022-08 | reference | 7 | 5 | 2 | 7 | +0.08 | -0.55 | -0.90 | +0.15 | -0.45 | -0.71 | +0.79 | +0.92 |
| 2022-09 | reference | 10 | 5 | 5 | 10 | -1.26 | -0.59 | +0.45 | -0.50 | +0.94 | +0.21 | +1.01 | +0.97 |
| 2022-10 | reference | 10 | 5 | 5 | 10 | -0.27 | +0.12 | +0.48 | +0.21 | +0.00 | -0.38 | +0.72 | +0.65 |
| 2022-11 | reference | 10 | 5 | 5 | 10 | -0.10 | +0.47 | +0.19 | -0.08 | +0.18 | +0.46 | +0.79 | +0.72 |
| 2022-12 | reference | 10 | 5 | 5 | 10 | -0.15 | -0.19 | +0.37 | -0.18 | +0.37 | -0.50 | +0.93 | +0.91 |
| 2023-01 | reference | 11 | 6 | 5 | 11 | -0.05 | +0.18 | -0.14 | -0.28 | +0.00 | -0.09 | +0.75 | +0.68 |
| 2023-02 | reference | 9 | 4 | 5 | 9 | +1.15 | +0.50 | -0.48 | -0.38 | -0.45 | +0.25 | +0.89 | +0.73 |
| 2023-03 | reference | 11 | 6 | 5 | 11 | -0.03 | +0.10 | -0.67 | +0.21 | -0.45 | +0.40 | +0.93 | +0.89 |
| 2023-04 | reference | 9 | 4 | 5 | 9 | +0.36 | -0.77 | +0.00 | +0.13 | +0.37 | +0.02 | +1.34 | +0.86 |
| 2023-05 | reference | 11 | 6 | 5 | 11 | +0.00 | -0.80 | +0.26 | +0.20 | +0.37 | -0.09 | +0.76 | +0.82 |
| 2023-06 | reference | 10 | 5 | 5 | 10 | -0.22 | +0.05 | +0.66 | +0.12 | +0.67 | +0.86 | +0.74 | +1.21 |
| 2023-07 | outside | 10 | 5 | 5 | 10 | -0.66 | +0.13 | +0.84 | +0.24 | +1.18 | -0.36 | +0.68 | +1.06 |
| 2023-08 | outside | 11 | 5 | 6 | 11 | -0.30 | +0.69 | +0.13 | -0.37 | +0.37 | -0.16 | +0.97 | +0.82 |
| 2023-09 | outside | 9 | 5 | 4 | 9 | -0.61 | +0.11 | +0.00 | -0.16 | +0.00 | -0.44 | +0.84 | +0.60 |
| 2023-10 | outside | 11 | 5 | 6 | 11 | -0.10 | -0.73 | +0.00 | +0.46 | +0.37 | -0.54 | +0.86 | +0.87 |
| 2023-11 | outside | 10 | 5 | 5 | 10 | +0.46 | -0.43 | -0.22 | +0.32 | -0.45 | -0.83 | +1.03 | +1.17 |
| 2023-12 | outside | 10 | 5 | 5 | 10 | -0.21 | -0.93 | -0.14 | -0.47 | +0.00 | +0.41 | +1.06 | +0.84 |
| 2024-01 | outside | 10 | 5 | 5 | 10 | +0.44 | -0.89 | -0.05 | -0.03 | +0.00 | -0.17 | +1.02 | +1.00 |
| 2024-02 | outside | 10 | 5 | 5 | 10 | -0.21 | -0.39 | -0.15 | +0.35 | +0.00 | -0.36 | +1.39 | +1.21 |
| 2024-03 | outside | 10 | 5 | 5 | 10 | -0.69 | -0.87 | +0.40 | +0.25 | +0.52 | +0.36 | +1.09 | +0.60 |
| 2024-04 | outside | 11 | 6 | 5 | 11 | -0.16 | -0.30 | -0.14 | -0.46 | +0.00 | +0.26 | +0.77 | +0.82 |
| 2024-05 | outside | 10 | 5 | 5 | 10 | +0.16 | -0.21 | +0.19 | -0.03 | +0.18 | +0.59 | +0.68 | +0.86 |
| 2024-06 | outside | 10 | 5 | 5 | 10 | +0.18 | -0.00 | +0.24 | -0.35 | -0.04 | +0.12 | +0.59 | +0.80 |
| 2024-07 | during | 10 | 5 | 5 | 10 | -0.17 | +0.13 | +0.13 | -0.14 | -0.04 | +0.21 | +0.55 | +1.15 |
| 2024-08 | during | 10 | 5 | 5 | 10 | +0.18 | -0.19 | +0.00 | -0.12 | +0.52 | +0.34 | +1.09 | +1.51 |
| 2024-09 | during | 10 | 5 | 5 | 10 | -0.15 | +0.29 | -0.11 | -0.44 | -0.22 | +0.18 | +1.01 | +0.94 |
| 2024-10 | during | 10 | 5 | 5 | 10 | -0.06 | -0.18 | +0.37 | -0.58 | +0.37 | +0.33 | +0.73 | +1.13 |
| 2024-11 | during | 10 | 5 | 5 | 10 | +0.29 | -0.01 | -0.24 | -0.27 | +0.18 | -0.31 | +0.75 | +0.74 |
| 2024-12 | during | 10 | 5 | 5 | 10 | -0.52 | -0.32 | +0.24 | +0.06 | +0.65 | -0.26 | +1.31 | +1.31 |
| 2025-01 | during | 11 | 6 | 5 | 11 | +0.37 | +0.42 | +0.26 | +0.16 | +0.37 | +0.27 | +0.73 | +0.88 |
| 2025-02 | during | 9 | 4 | 5 | 9 | +0.35 | -0.10 | +0.00 | -0.87 | +0.67 | +0.75 | +0.95 | +1.22 |
| 2025-03 | during | 11 | 5 | 6 | 11 | +0.24 | -0.74 | +0.13 | -0.46 | +0.37 | +0.72 | +1.26 | +1.24 |
| 2025-04 | during | 18 | 8 | 10 | 18 | -0.32 | -0.57 | +0.31 | -0.84 | +0.67 | +1.82 | +0.83 | +1.46 |
| 2025-05 | during | 19 | 9 | 10 | 19 | +0.16 | -0.12 | +0.26 | -0.91 | +0.67 | +0.51 | +1.08 | +0.91 |
| 2025-06 | during | 19 | 9 | 10 | 19 | -0.26 | -0.58 | +0.58 | -0.65 | +0.94 | +1.55 | +0.90 | +1.10 |
| 2025-07 | during | 19 | 9 | 10 | 19 | -0.21 | -0.66 | +0.37 | -0.16 | +0.67 | +0.34 | +1.03 | +1.15 |
| 2025-08 | during | 20 | 10 | 10 | 20 | -0.41 | -0.85 | +1.00 | -0.36 | +1.03 | +0.45 | +0.90 | +1.15 |
| 2025-09 (pre) | during | 3 | 1 | 2 | 3 | -0.19 | +0.19 | -0.48 | -0.49 | +0.00 | +1.02 | +0.42 | +0.89 |
| 2025-09 (post) | post | 17 | 9 | 8 | 17 | +0.13 | -0.23 | +0.58 | -0.35 | +0.67 | +1.01 | +0.72 | +1.12 |
| 2025-10 | post | 21 | 10 | 11 | 21 | +0.03 | +0.20 | +0.26 | -0.48 | +0.37 | +0.46 | +1.12 | +0.88 |
| 2025-11 | post | 20 | 10 | 10 | 20 | -0.09 | -0.15 | +0.13 | -0.16 | +0.00 | +0.02 | +0.86 | +0.99 |
| 2025-12 | post | 3 | 2 | 1 | 3 | +0.19 | +0.70 | +0.67 | +0.06 | +0.37 | +1.04 | +1.74 | +1.27 |

A block magnitude is computed **per chip** (the RMS of that chip's standardised dimensions) and the table and the plate show the **median of those chip values** -- so `R_B` is the typical joint displacement of a single overpass, and it is *not* the RMS of the two columns beside it: a bin whose chips scatter widely in opposite directions has a near-zero median in each dimension and a large magnitude per chip.

## The frozen reference, the held-out year, the interval and the post period

Each cell is the median over the chips of that class, in the reference's own scales. The reference column is zero by construction (the standardisation rests on its median); it is printed so the other columns are read against the frame they are measured in.

| quantity | reference | held-out pre-intervention | intervention | post |
| --- | --- | --- | --- | --- |
| `gamma2` | +0.00 | -0.10 | -0.01 | -0.03 |
| `P` | +0.00 | -0.30 | -0.27 | -0.22 |
| `R_B` | +0.84 | +0.88 | +0.91 | +1.04 |
| `A` | +0.00 | +0.00 | +0.26 | +0.37 |
| `D` | +0.00 | -0.13 | -0.44 | -0.28 |
| `F` | +0.00 | +0.00 | +0.37 | +0.37 |
| `S` | +0.00 | -0.05 | +0.34 | +0.51 |
| `R_G` | +0.85 | +0.86 | +1.10 | +1.02 |

n = reference **119** chips / **12** bins, outside **122** chips / **12** bins, during **189** chips / **15** bins, post **61** chips / **4** bins.

### Each held-out class against the reference

The interval is R5's seeded circular block bootstrap on `median(class) - median(reference)`; `outside` counts the class's monthly states beyond the reference envelope (above q95 for the two magnitudes, outside [q05, q95] for a signed z).

A bound printed as `+0.000` is not a small number rounded to nothing: it is a draw in which the two medians coincide exactly, which is what a quantity with a mass at zero does (`A`, `F` sit at zero in at least one class here). Such an interval reads `excludes 0 = no`, and that is the honest reading of it.

**held-out pre-intervention** — 122 chips, 12 bins, 2023-07 → 2024-06.

| quantity | median | monthly state | delta | interval | excludes 0 | months outside | direction |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `gamma2` | -0.102 | -0.183 | -0.102 | [-0.430, +0.215] | no | 0 | down |
| `P` | -0.299 | -0.345 | -0.299 | [-0.721, +0.119] | no | 4 | down |
| `R_B` | +0.885 | +0.918 | +0.049 | [-0.062, +0.215] | no | 1 | up |
| `A` | +0.000 | +0.000 | +0.000 | [-0.278, +0.278] | no | 1 | flat |
| `D` | -0.131 | -0.032 | -0.131 | [-0.397, +0.361] | no | 7 | down |
| `F` | +0.000 | +0.000 | +0.000 | [-0.365, +0.365] | no | 1 | flat |
| `S` | -0.049 | -0.166 | -0.049 | [-0.391, +0.301] | no | 1 | down |
| `R_G` | +0.859 | +0.851 | +0.006 | [-0.096, +0.189] | no | 2 | up |

**intervention period** — 189 chips, 15 bins, 2024-07 → 2025-09.

| quantity | median | monthly state | delta | interval | excludes 0 | months outside | direction |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `gamma2` | -0.008 | -0.154 | -0.008 | [-0.365, +0.309] | no | 0 | down |
| `P` | -0.271 | -0.176 | -0.271 | [-0.640, +0.097] | no | 1 | down |
| `R_B` | +0.914 | +0.898 | +0.078 | [-0.062, +0.233] | no | 2 | up |
| `A` | +0.256 | +0.239 | +0.256 | [+0.000, +0.516] | no | 2 | up |
| `D` | -0.440 | -0.441 | -0.440 | [-0.654, -0.094] | yes | 8 | down |
| `F` | +0.365 | +0.520 | +0.365 | [+0.000, +0.674] | no | 2 | up |
| `S` | +0.343 | +0.339 | +0.343 | [+0.092, +0.608] | yes | 5 | up |
| `R_G` | +1.096 | +1.148 | +0.243 | [+0.092, +0.416] | yes | 10 | up |

**post-intervention** — 61 chips, 4 bins, 2025-09 → 2025-12.

| quantity | median | monthly state | delta | interval | excludes 0 | months outside | direction |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `gamma2` | -0.033 | +0.080 | -0.033 | [-0.475, +0.318] | no | 0 | down |
| `P` | -0.220 | +0.028 | -0.220 | [-0.452, +0.325] | no | 1 | down |
| `R_B` | +1.044 | +0.993 | +0.209 | [-0.073, +0.330] | no | 1 | up |
| `A` | +0.371 | +0.418 | +0.371 | [-0.133, +0.579] | no | 2 | up |
| `D` | -0.276 | -0.257 | -0.276 | [-0.521, +0.086] | no | 1 | down |
| `F` | +0.365 | +0.365 | +0.365 | [-0.365, +0.674] | no | 0 | up |
| `S` | +0.507 | +0.734 | +0.507 | [+0.158, +0.928] | yes | 2 | up |
| `R_G` | +1.024 | +1.053 | +0.171 | [+0.016, +0.323] | yes | 2 | up |

The bootstrap resamples the reference first and then the class, from a fresh `random.Random(7)` per class, so a class's interval does not depend on which other classes were computed.

## Repair/intervention assessment

The nine questions below are the registered ones (`QUESTIONS` in the generator), answered in that order and only from the frozen classes and the registered rules; the word in parentheses under each heading is the question's key in the JSON.

The three statements are kept apart throughout:

* **Directional departure** -- a signed monthly median leaves the monthly reference envelope; a statement about direction, not about size.
* **Departure magnitude** -- `R_B` or `R_G` above the monthly reference q95; a statement about the joint size of the displacement, with no direction.
* **Sustained departure** -- a run of at least **2** consecutive months that meet the criterion, decided by the fixed rule and not by the shape of the curve.

**Q1 · Does the vector depart from the frozen reference while the documented works run?** (`onset`)

'sustained unusual geometry departure' -- over the 15 monthly bins of the documented interval, R_B exceeds the reference's monthly q95 in 2 of them and R_G in 10. The word is a function of those two persistence rules and of nothing else.

**Q2 · Which response block departs - behaviour (γ², P), geometry (A, D, F, S), or both?** (`block`)

The sustained departure sits in the geometry block (A, D, F, S): R_G is above its monthly q95 for 3 consecutive months (2025-02 → 2025-04), R_B for a single month (2024-12), against a threshold of 2 consecutive months.

**Q3 · Which channel carries the departure, and in which direction does it sit in the reference's frame?** (`channel`)

3 of the 8 registered quantities leave the reference's monthly envelope in at least one interval month while their R5 interval excludes zero: R_G (up, 10 months), z(D) (down, 8 months), z(S) (up, 5 months).

**Q4 · Had the vector already departed before the documented start?** (`precursor`)

In the held-out pre-intervention year, R_B exceeds the reference's monthly q95 in 1 of 12 months and R_G in 2 of 12; that class's own verdict is 'isolated unusual departure'.

**Q5 · Is a departure sustained across consecutive months of the interval?** (`persistence`)

R_B: a single month (2024-12); R_G: 3 consecutive months (2025-02 → 2025-04) above the reference's monthly q95. Sustained means at least 2 consecutive such months (R1).

**Q6 · Does the vector move back towards the reference after the documented end, or does it stay displaced?** (`recovery`)

After the documented end the post class's monthly median sits inside q05..q95 for both R_B and R_G; counted month by month, 1 of 4 post bins above the R_B q95 and 2 of 4 above the R_G q95. The extent word, from the post block's own persistence rules, is 'confined to the intervention interval'.

**Q7 · Is the post-intervention state distinguishable from the held-out pre-intervention year?** (`distinction`)

Post minus held-out pre-intervention: the R5 interval excludes zero for 1 of the 8 registered quantities (z(S)). This is a descriptive difference between two held-out classes, not a test.

**Q8 · Do the pronounced months coincide with environmental variation?** (`environment`)

14 pronounced months outside the reference (2023-11, 2024-02, 2024-07, 2024-08, 2024-10, 2024-12, 2025-02, 2025-03, 2025-04, 2025-06, 2025-07, 2025-08, 2025-09 (post), 2025-12); 6 of them also carry a co-variate departure larger than one robust scale (2024-02, 2024-07, 2024-08, 2025-02, 2025-06, 2025-07). The co-variates are reported beside the OCV quantities: no correction, no attribution.

**Q9 · Does the pattern hold in both orbit directions, or is it one pass geometry?** (`orbit`)

Splitting the intervention interval by pass direction, the departure from the reference has z(gamma2): opposite sign in the two orbits; z(P): same sign in both orbits; R_B: opposite sign in the two orbits; R_G: same sign in both orbits. The reference is pooled, so this split is the check that the pattern is not one pass geometry's.

**Verdict — sustained unusual geometry departure.** The extent word is **confined to the intervention interval**. The verdict is derived from the intervention interval's two block magnitudes' persistence rules alone (`verdict_word` on `blocks.during.persistence.R_B` / `R_G`) and the extent word from the post block's persistence, so neither reads the co-variates, the orbit split or any interval. It names no damage, no failure precursor, no repair effect and no cause; the plate measures the observability of the crossing.

## The orbit split (Q9)

The reference is pooled, so the pooled median of a class says nothing about whether both pass geometries moved. This is a chip-level split in the pooled frame, and the verdict never reads it: a departure that only one orbit carries is reported as one pass geometry's departure, not as a site-wide one.

| class | orbit | n | chips' medians: `gamma2` | `P` | `R_B` | `R_G` |
| --- | --- | --- | --- | --- | --- | --- |
| outside | ASC | 61 | -0.211 | -0.304 | +0.715 | +0.714 |
| outside | DESC | 61 | +0.011 | -0.257 | +1.064 | +1.116 |
| during | ASC | 91 | +0.090 | -0.304 | +0.798 | +0.868 |
| during | DESC | 98 | -0.173 | -0.248 | +1.206 | +1.339 |
| post | ASC | 31 | -0.131 | -0.394 | +0.830 | +1.045 |
| post | DESC | 30 | +0.592 | +0.173 | +1.176 | +0.962 |

Sign of the departure from the reference, per orbit: outside -- `gamma2` opposite sign in the two orbits, `P` same sign in both orbits, `R_B` opposite sign in the two orbits, `R_G` opposite sign in the two orbits; during -- `gamma2` opposite sign in the two orbits, `P` same sign in both orbits, `R_B` opposite sign in the two orbits, `R_G` same sign in both orbits; post -- `gamma2` opposite sign in the two orbits, `P` opposite sign in the two orbits, `R_B` opposite sign in the two orbits, `R_G` same sign in both orbits.

## The post state against the held-out year (Q7)

`median(post) - median(outside)`, with the same seeded block bootstrap; a descriptive interval on a difference between two held-out classes, not a test, and not a verdict input.

| quantity | median(post) - median(held-out) | interval | excludes 0 |
| --- | --- | --- | --- |
| `gamma2` | +0.069 | [-0.358, +0.419] | no |
| `P` | +0.079 | [-0.203, +0.616] | no |
| `R_B` | +0.159 | [-0.164, +0.282] | no |
| `A` | +0.371 | [-0.133, +0.579] | no |
| `D` | -0.146 | [-0.551, +0.148] | no |
| `F` | +0.365 | [-0.365, +0.674] | no |
| `S` | +0.557 | [+0.207, +1.026] | yes |
| `R_G` | +0.165 | [-0.025, +0.290] | no |

## The environmental context (Q8)

The co-variates are reported beside the OCV quantities, not inside them: no weather correction is applied and no causal attribution is made. A pronounced month is one whose `R_B` or `R_G` exceeds its class's reference q95.

| bin | class | R_B | R_G | z(γ²) | z(P) | `precipitation_mm` | `temperature_c` | above q95 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 2023-07 | outside | +0.675 | +1.056 | -0.66 | +0.13 | +2.28 | +1.23 | none |
| 2023-08 | outside | +0.975 | +0.823 | -0.30 | +0.69 | +2.36 | +1.04 | none |
| 2023-09 | outside | +0.838 | +0.599 | -0.61 | +0.11 | -0.67 | +0.99 | none |
| 2023-10 | outside | +0.860 | +0.869 | -0.10 | -0.73 | +0.00 | +0.13 | none |
| 2023-11 | outside | +1.029 | +1.168 | +0.46 | -0.43 | +0.59 | -0.47 | R_G |
| 2023-12 | outside | +1.060 | +0.838 | -0.21 | -0.93 | +1.60 | -0.75 | none |
| 2024-01 | outside | +1.025 | +0.997 | +0.44 | -0.89 | +1.52 | -1.17 | none |
| 2024-02 | outside | +1.390 | +1.212 | -0.21 | -0.39 | +1.77 | -0.38 | R_B, R_G |
| 2024-03 | outside | +1.094 | +0.604 | -0.69 | -0.87 | -0.51 | -0.42 | none |
| 2024-04 | outside | +0.775 | +0.820 | -0.16 | -0.30 | +0.00 | +0.13 | none |
| 2024-05 | outside | +0.678 | +0.865 | +0.16 | -0.21 | -0.00 | +0.69 | none |
| 2024-06 | outside | +0.592 | +0.800 | +0.18 | -0.00 | +0.17 | +1.08 | none |
| 2024-07 | during | +0.552 | +1.151 | -0.17 | +0.13 | -0.51 | +1.28 | R_G |
| 2024-08 | during | +1.087 | +1.505 | +0.18 | -0.19 | -0.67 | +1.45 | R_G |
| 2024-09 | during | +1.010 | +0.941 | -0.15 | +0.29 | -0.42 | +0.69 | none |
| 2024-10 | during | +0.730 | +1.131 | -0.06 | -0.18 | -0.59 | +0.12 | R_G |
| 2024-11 | during | +0.749 | +0.739 | +0.29 | -0.01 | -0.51 | -0.62 | none |
| 2024-12 | during | +1.315 | +1.306 | -0.52 | -0.32 | +0.76 | -0.92 | R_B, R_G |
| 2025-01 | during | +0.735 | +0.879 | +0.37 | +0.42 | +1.69 | -1.08 | none |
| 2025-02 | during | +0.945 | +1.218 | +0.35 | -0.10 | -0.67 | -1.39 | R_G |
| 2025-03 | during | +1.255 | +1.238 | +0.24 | -0.74 | -0.67 | -0.62 | R_B, R_G |
| 2025-04 | during | +0.829 | +1.463 | -0.32 | -0.57 | -0.67 | +0.19 | R_G |
| 2025-05 | during | +1.081 | +0.911 | +0.16 | -0.12 | -0.34 | +0.29 | none |
| 2025-06 | during | +0.898 | +1.096 | -0.26 | -0.58 | +0.51 | +1.04 | R_G |
| 2025-07 | during | +1.030 | +1.148 | -0.21 | -0.66 | +1.01 | +1.05 | R_G |
| 2025-08 | during | +0.898 | +1.151 | -0.41 | -0.85 | -0.67 | +0.95 | R_G |
| 2025-09 (pre) | during | +0.420 | +0.889 | -0.19 | +0.19 | -0.51 | +1.30 | none |
| 2025-09 (post) | post | +0.717 | +1.116 | +0.13 | -0.23 | -0.67 | +0.57 | R_G |
| 2025-10 | post | +1.121 | +0.880 | +0.03 | +0.20 | +0.34 | -0.03 | none |
| 2025-11 | post | +0.864 | +0.991 | -0.09 | -0.15 | -0.59 | -0.60 | none |
| 2025-12 | post | +1.737 | +1.267 | +0.19 | +0.70 | -0.67 | -0.99 | R_B, R_G |

Pronounced months: 2023-11, 2024-02, 2024-07, 2024-08, 2024-10, 2024-12, 2025-02, 2025-03, 2025-04, 2025-06, 2025-07, 2025-08, 2025-09 (post), 2025-12. Of these, 2024-02, 2024-07, 2024-08, 2025-02, 2025-06, 2025-07 also carry a co-variate departure larger than one robust scale. This is a temporal coincidence, not an explanation.

## The panels

| panel | block | what it draws |
| --- | --- | --- |
| **a1** | behaviour | monthly z(gamma2) |
| **a2** | behaviour | monthly z(P) |
| **a3** | behaviour | behaviour plane, z(gamma2) x z(P) |
| **a4** | behaviour | monthly R_B |
| **b1** | geometry | monthly z(A) |
| **b2** | geometry | monthly z(F) |
| **b3** | geometry | geometry plane, z(A) x z(F) |
| **b4** | geometry | monthly R_G |
| **c1** | context | monthly co-variates, in their own reference frame |
| **c2** | context | cadence, chips per bin by orbit |
| **c3** | context | echo-bearing share per bin |
| **c4** | context | trajectory matrix, z per dimension and bin |

The plate is set at the column's own width -- `figsize = (7.205 in, 6.628 in)`, 630 dpi, so 183.0 mm x 168.4 mm -- and saved with `bbox_inches="tight"`, which crops the figure to its ink rather than to that rectangle, so the delivered file is a few millimetres wider than the setting (the colour bar's ticks) and shorter than it (the empty margin the grid leaves below the key). The header is packed to the panel block's own width, measured from column one's y label to column four's right spine, and starts at that same left edge; panel c4's title is broken in two. Either, left as one line, would be wider than the block, and the tight box would then widen the *file* to hold the text instead of the text to fit the block.

## What the plate does not claim

* **No damage, failure or safety verdict.** Every ordinate is a deviation from the crossing's *own* reference, in that reference's robust scales; two plates are not comparable in magnitude, only in shape and in where the movement sits.
* **No attribution to the works.** The interval is the record's; a departure that coincides with it in time is not evidence that the works caused it, and no weather correction is applied.
* **No significance on the curve.** The intervals are differences of medians, not bands on the curve: the monthly medians are drawn with their `n` and without one.
* **No directional score.** `R_B` and `R_G` keep their frozen definitions, and the verdict never reads the orbit split, the co-variates or any interval.
* **The reference maximum is not a boundary.** The envelope is the q05/q50/q95 of the monthly reference states; the maximum is printed only so a reader can see how far q95 sits from it.
* **No interpolation.** A boundary, a hole and a class change break the line: nothing is drawn across them, and no bin mixes two classes.
* **A month is a bin, not a date.** Its median can mix both pass directions; panel c2 says how many chips of each the bin holds.

## Pins and reproduction

| pin | value |
| --- | --- |
| `classes` | `{'reference': 119, 'outside': 122, 'during': 189, 'post': 61}` |
| `cut_rows` | `491` |
| `event` | `('2024-07-01', '2025-09-05')` |
| `monthly_bins` | `43` |
| `reference_months` | `12` |
| `reference_window` | `('2022-07', '2023-06')` |
| `state_windows` | `{'W2_bestand': ('2022-07-02', '2024-06-29'), 'W3_sanierung': ('2024-07-03', '2025-09-04'), 'W4_post': ('2025-09-05', '2025-12-04')}` |
| `table` | `data/werder_havel/support/werder_havel_ocv_channels.csv` |
| `table_rows` | `991` |
| `unassigned` | `0` |

`--verify` re-derives every number of this report from `data/werder_havel/support/werder_havel_ocv_channels.csv` and fails if the JSON, the report and the plate disagree with the table; `--selftest` checks the rules and the pins above without writing anything. The bootstrap is seeded (7, block 3, 2000 draws), no wall-clock value enters any output, and the plate is drawn twice and compared byte for byte, so the JSON, the report and the figure reproduce deterministically.

```bash
bash figures/werder_havel/fig_werder_havel_ocv_trajectory.sh            # JSON, report, plate
bash figures/werder_havel/fig_werder_havel_ocv_trajectory.sh --no-figures # numbers only
bash figures/werder_havel/fig_werder_havel_ocv_trajectory.sh --verify      # re-derive and compare
bash figures/werder_havel/fig_werder_havel_ocv_trajectory.sh --selftest    # rules and pins, nothing written
```

The generator of the cut is `code/werder_havel/werder_havel_ocv_channels_csv.py`; the template this plate follows is `code/cases/ocv_trajectory.py` (sha256 `2b87ffb0092fb597`), imported rather than copied, so the seed, the transforms, the block definitions and the monthly-bin machinery cannot drift between the two plates.
