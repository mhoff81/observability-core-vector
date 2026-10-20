# SDS of the OCV dimensions on the support-sized chips

Bautzen and LUMO, recomputed on the **support-sized chip** (odd, anchor-centred, guard 30 m, reference ring 30–200 m) and beside it on the committed table. The metric is the project's own: `SDS = 10 · |Cliff's delta|`, computed with `code/lumo/lumo_ocv_stats.cliffs_delta`; the interval is a percentile bootstrap (10,000 resamples, seed 7) of **one comparison**.

## 0. Reading rules

1. **The primary analysis is per channel and per state pair.** The observation channel is the orbit; ASC and DESC see different geometries and are never pooled into the primary number.
2. **A maximum always travels with its winning pair** — the signed delta, that pair's interval and its two sample sizes.
3. **The interval belongs to one comparison.** For a maximum over several comparisons the report prints the *winning pair's* interval and says so; it is not an interval for the selection.
4. **`pooled_all_damage` is a compatibility pin and context**, not the primary evidence: it is the basis the committed LUMO reference used.

## 1. Primary analysis — per channel, per state pair

Support chip. ★ marks the winning comparison of its dimension (the one the compact table in §2 carries).

### LUMO lattice tower (Hannover)

**ASCENDING** — rows per class: healthy 40, DAM 3 8, DAM 4 5, DAM 6 4

| dimension | state pair | SDS | band | delta | 95% CI | n1 / n2 | flags |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `gamma2` | healthy – DAM 3 | 0.06 | very weak | +0.006 | [-0.412, +0.419] | 40 / 8 | CI includes 0 |
| `gamma2` | healthy – DAM 4 | 1.90 | very weak | -0.190 | [-0.760, +0.400] | 40 / 5 | CI includes 0 |
| `gamma2` | healthy – DAM 6 ★ | 4.88 | moderate | +0.487 | [-0.013, +0.900] | 40 / 4 | CI includes 0, n<5 |
| `P` | healthy – DAM 3 ★ | 4.69 | moderate | -0.469 | [-0.744, -0.150] | 40 / 8 | – |
| `P` | healthy – DAM 4 | 0.60 | very weak | +0.060 | [-0.450, +0.590] | 40 / 5 | CI includes 0 |
| `P` | healthy – DAM 6 | 3.12 | weak | +0.312 | [-0.163, +0.775] | 40 / 4 | CI includes 0, n<5 |
| `D` | healthy – DAM 3 | 1.03 | very weak | +0.103 | [-0.297, +0.487] | 40 / 8 | CI includes 0 |
| `D` | healthy – DAM 4 ★ | 3.35 | weak | +0.335 | [-0.135, +0.760] | 40 / 5 | CI includes 0 |
| `D` | healthy – DAM 6 | 1.94 | very weak | +0.194 | [-0.600, +0.975] | 40 / 4 | CI includes 0, n<5 |
| `F` | healthy – DAM 3 | 0.72 | very weak | +0.072 | [-0.344, +0.475] | 40 / 8 | CI includes 0 |
| `F` | healthy – DAM 4 | 1.25 | very weak | -0.125 | [-0.490, +0.225] | 40 / 5 | CI includes 0 |
| `F` | healthy – DAM 6 ★ | 7.12 | strong | -0.713 | [-0.950, -0.406] | 40 / 4 | n<5 |
| `A` | healthy – DAM 3 | 1.16 | very weak | +0.116 | [-0.272, +0.491] | 40 / 8 | CI includes 0 |
| `A` | healthy – DAM 4 | 1.60 | very weak | +0.160 | [-0.415, +0.670] | 40 / 5 | CI includes 0 |
| `A` | healthy – DAM 6 ★ | 5.94 | moderate | -0.594 | [-0.963, -0.100] | 40 / 4 | n<5 |
| `S` | healthy – DAM 3 ★ | 2.06 | weak | +0.206 | [-0.219, +0.606] | 40 / 8 | CI includes 0 |
| `S` | healthy – DAM 4 | 0.10 | very weak | -0.010 | [-0.440, +0.420] | 40 / 5 | CI includes 0 |
| `S` | healthy – DAM 6 | 1.50 | very weak | -0.150 | [-0.969, +0.700] | 40 / 4 | CI includes 0, n<5 |


**DESCENDING** — rows per class: healthy 81, DAM 3 17, DAM 4 11, DAM 6 12

| dimension | state pair | SDS | band | delta | 95% CI | n1 / n2 | flags |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `gamma2` | healthy – DAM 3 | 1.42 | very weak | -0.142 | [-0.454, +0.182] | 81 / 17 | CI includes 0 |
| `gamma2` | healthy – DAM 4 ★ | 2.88 | weak | +0.288 | [-0.071, +0.623] | 81 / 11 | CI includes 0 |
| `gamma2` | healthy – DAM 6 | 1.48 | very weak | -0.148 | [-0.486, +0.212] | 81 / 12 | CI includes 0 |
| `P` | healthy – DAM 3 ★ | 2.50 | weak | -0.250 | [-0.521, +0.040] | 81 / 17 | CI includes 0 |
| `P` | healthy – DAM 4 | 1.34 | very weak | +0.134 | [-0.240, +0.493] | 81 / 11 | CI includes 0 |
| `P` | healthy – DAM 6 | 0.29 | very weak | +0.029 | [-0.342, +0.399] | 81 / 12 | CI includes 0 |
| `D` | healthy – DAM 3 | 1.36 | very weak | -0.136 | [-0.432, +0.163] | 81 / 17 | CI includes 0 |
| `D` | healthy – DAM 4 ★ | 3.66 | weak | +0.366 | [+0.104, +0.602] | 81 / 11 | – |
| `D` | healthy – DAM 6 | 2.49 | weak | -0.249 | [-0.541, +0.064] | 81 / 12 | CI includes 0 |
| `F` | healthy – DAM 3 | 0.54 | very weak | +0.054 | [-0.217, +0.326] | 81 / 17 | CI includes 0 |
| `F` | healthy – DAM 4 ★ | 3.08 | weak | -0.308 | [-0.543, -0.043] | 81 / 11 | – |
| `F` | healthy – DAM 6 | 2.63 | weak | +0.263 | [-0.009, +0.514] | 81 / 12 | CI includes 0 |
| `A` | healthy – DAM 3 | 0.47 | very weak | +0.047 | [-0.248, +0.344] | 81 / 17 | CI includes 0 |
| `A` | healthy – DAM 4 | 0.79 | very weak | -0.079 | [-0.420, +0.259] | 81 / 11 | CI includes 0 |
| `A` | healthy – DAM 6 ★ | 1.56 | very weak | +0.156 | [-0.168, +0.461] | 81 / 12 | CI includes 0 |
| `S` | healthy – DAM 3 | 0.12 | very weak | +0.012 | [-0.316, +0.339] | 81 / 17 | CI includes 0 |
| `S` | healthy – DAM 4 ★ | 2.27 | weak | -0.227 | [-0.586, +0.159] | 81 / 11 | CI includes 0 |
| `S` | healthy – DAM 6 | 2.25 | weak | +0.225 | [-0.072, +0.498] | 81 / 12 | CI includes 0 |


### Bautzen OpenLabs bridge

**ASC** — rows per class: RS 14, DS1 46, DS2 12

| dimension | state pair | SDS | band | delta | 95% CI | n1 / n2 | flags |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `gamma2` | RS – DS1 ★ | 0.89 | very weak | -0.089 | [-0.464, +0.292] | 12 / 32 | CI includes 0 |
| `gamma2` | RS – DS2 | 0.45 | very weak | -0.045 | [-0.530, +0.455] | 12 / 11 | CI includes 0 |
| `P` | RS – DS1 ★ | 1.52 | very weak | +0.152 | [-0.194, +0.488] | 13 / 40 | CI includes 0 |
| `P` | RS – DS2 | 0.64 | very weak | -0.064 | [-0.513, +0.397] | 13 / 12 | CI includes 0 |
| `D` | RS – DS1 | 0.06 | very weak | +0.006 | [-0.340, +0.362] | 13 / 40 | CI includes 0 |
| `D` | RS – DS2 ★ | 1.22 | very weak | -0.122 | [-0.564, +0.353] | 13 / 12 | CI includes 0 |
| `F` | RS – DS1 ★ | 2.02 | weak | -0.202 | [-0.497, +0.102] | 14 / 46 | CI includes 0 |
| `F` | RS – DS2 | 0.71 | very weak | +0.071 | [-0.310, +0.470] | 14 / 12 | CI includes 0 |
| `A` | RS – DS1 ★ | 1.13 | very weak | -0.113 | [-0.446, +0.225] | 14 / 46 | CI includes 0 |
| `A` | RS – DS2 | 0.71 | very weak | +0.071 | [-0.393, +0.518] | 14 / 12 | CI includes 0 |
| `S` | RS – DS1 | 0.58 | very weak | -0.058 | [-0.387, +0.283] | 13 / 40 | CI includes 0 |
| `S` | RS – DS2 ★ | 1.22 | very weak | +0.122 | [-0.353, +0.577] | 13 / 12 | CI includes 0 |


**DESC** — rows per class: RS 14, DS1 46, DS2 10

| dimension | state pair | SDS | band | delta | 95% CI | n1 / n2 | flags |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `gamma2` | RS – DS1 ★ | 0.31 | very weak | -0.031 | [-0.401, +0.342] | 14 / 46 | CI includes 0 |
| `gamma2` | RS – DS2 | 0.14 | very weak | -0.014 | [-0.471, +0.457] | 14 / 10 | CI includes 0 |
| `P` | RS – DS1 | 1.49 | very weak | -0.149 | [-0.509, +0.225] | 14 / 46 | CI includes 0 |
| `P` | RS – DS2 ★ | 2.29 | weak | -0.229 | [-0.657, +0.243] | 14 / 10 | CI includes 0 |
| `D` | RS – DS1 | 1.23 | very weak | +0.123 | [-0.248, +0.472] | 14 / 46 | CI includes 0 |
| `D` | RS – DS2 ★ | 2.14 | weak | +0.214 | [-0.257, +0.671] | 14 / 10 | CI includes 0 |
| `F` | RS – DS1 | 0.43 | very weak | +0.043 | [-0.280, +0.363] | 14 / 46 | CI includes 0 |
| `F` | RS – DS2 ★ | 2.14 | weak | +0.214 | [-0.229, +0.614] | 14 / 10 | CI includes 0 |
| `A` | RS – DS1 | 1.52 | very weak | +0.152 | [-0.207, +0.495] | 14 / 46 | CI includes 0 |
| `A` | RS – DS2 ★ | 3.57 | weak | +0.357 | [-0.100, +0.757] | 14 / 10 | CI includes 0 |
| `S` | RS – DS1 | 1.68 | very weak | +0.168 | [-0.171, +0.509] | 14 / 46 | CI includes 0 |
| `S` | RS – DS2 ★ | 2.00 | weak | +0.200 | [-0.272, +0.629] | 14 / 10 | CI includes 0 |


## 2. Compact summary — max pairwise

The paper table. Every row names the state pair that produced the maximum, its signed delta, **that pair's** interval, and the two sample sizes — a maximum without them would read as a general discrimination the number does not support.

| asset | channel | dimension | max SDS | winning pair | delta | 95% CI (winning pair) | n1 / n2 | flags |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| lumo | ASCENDING | `gamma2` | 4.88 | healthy – DAM 6 | +0.487 | [-0.013, +0.900] | 40 / 4 | CI includes 0, n<5 |
| lumo | ASCENDING | `P` | 4.69 | healthy – DAM 3 | -0.469 | [-0.744, -0.150] | 40 / 8 | – |
| lumo | ASCENDING | `D` | 3.35 | healthy – DAM 4 | +0.335 | [-0.135, +0.760] | 40 / 5 | CI includes 0 |
| lumo | ASCENDING | `F` | 7.12 | healthy – DAM 6 | -0.713 | [-0.950, -0.406] | 40 / 4 | n<5 |
| lumo | ASCENDING | `A` | 5.94 | healthy – DAM 6 | -0.594 | [-0.963, -0.100] | 40 / 4 | n<5 |
| lumo | ASCENDING | `S` | 2.06 | healthy – DAM 3 | +0.206 | [-0.219, +0.606] | 40 / 8 | CI includes 0 |
| lumo | DESCENDING | `gamma2` | 2.88 | healthy – DAM 4 | +0.288 | [-0.071, +0.623] | 81 / 11 | CI includes 0 |
| lumo | DESCENDING | `P` | 2.50 | healthy – DAM 3 | -0.250 | [-0.521, +0.040] | 81 / 17 | CI includes 0 |
| lumo | DESCENDING | `D` | 3.66 | healthy – DAM 4 | +0.366 | [+0.104, +0.602] | 81 / 11 | – |
| lumo | DESCENDING | `F` | 3.08 | healthy – DAM 4 | -0.308 | [-0.543, -0.043] | 81 / 11 | – |
| lumo | DESCENDING | `A` | 1.56 | healthy – DAM 6 | +0.156 | [-0.168, +0.461] | 81 / 12 | CI includes 0 |
| lumo | DESCENDING | `S` | 2.27 | healthy – DAM 4 | -0.227 | [-0.586, +0.159] | 81 / 11 | CI includes 0 |
| bautzen | ASC | `gamma2` | 0.89 | RS – DS1 | -0.089 | [-0.464, +0.292] | 12 / 32 | CI includes 0 |
| bautzen | ASC | `P` | 1.52 | RS – DS1 | +0.152 | [-0.194, +0.488] | 13 / 40 | CI includes 0 |
| bautzen | ASC | `D` | 1.22 | RS – DS2 | -0.122 | [-0.564, +0.353] | 13 / 12 | CI includes 0 |
| bautzen | ASC | `F` | 2.02 | RS – DS1 | -0.202 | [-0.497, +0.102] | 14 / 46 | CI includes 0 |
| bautzen | ASC | `A` | 1.13 | RS – DS1 | -0.113 | [-0.446, +0.225] | 14 / 46 | CI includes 0 |
| bautzen | ASC | `S` | 1.22 | RS – DS2 | +0.122 | [-0.353, +0.577] | 13 / 12 | CI includes 0 |
| bautzen | DESC | `gamma2` | 0.31 | RS – DS1 | -0.031 | [-0.401, +0.342] | 14 / 46 | CI includes 0 |
| bautzen | DESC | `P` | 2.29 | RS – DS2 | -0.229 | [-0.657, +0.243] | 14 / 10 | CI includes 0 |
| bautzen | DESC | `D` | 2.14 | RS – DS2 | +0.214 | [-0.257, +0.671] | 14 / 10 | CI includes 0 |
| bautzen | DESC | `F` | 2.14 | RS – DS2 | +0.214 | [-0.229, +0.614] | 14 / 10 | CI includes 0 |
| bautzen | DESC | `A` | 3.57 | RS – DS2 | +0.357 | [-0.100, +0.757] | 14 / 10 | CI includes 0 |
| bautzen | DESC | `S` | 2.00 | RS – DS2 | +0.200 | [-0.272, +0.629] | 14 / 10 | CI includes 0 |

## 3. Compatibility pin and context — pooled over both channels

**Not the primary evidence.** This is the aggregation the committed LUMO reference used (`sds_pooled_all_damaged`), so it is reproduced here as a pin and reported for context. Pooling ASC and DESC mixes two observation geometries.

| asset | table | dimension | SDS | delta | 95% CI | n1 / n2 | flags |
| --- | --- | --- | --- | --- | --- | --- | --- |
| lumo | committed | `gamma2` | 0.10 | +0.010 | [-0.165, +0.185] | 121 / 57 | CI includes 0 |
| lumo | committed | `P` | 1.05 | -0.105 | [-0.292, +0.089] | 120 / 57 | CI includes 0 |
| lumo | committed | `D` | 0.86 | -0.086 | [-0.272, +0.098] | 121 / 57 | CI includes 0 |
| lumo | committed | `F` | 0.67 | -0.067 | [-0.244, +0.107] | 121 / 57 | CI includes 0 |
| lumo | committed | `A` | 0.81 | -0.081 | [-0.255, +0.090] | 121 / 57 | CI includes 0 |
| lumo | committed | `S` | 0.10 | -0.010 | [-0.196, +0.172] | 121 / 57 | CI includes 0 |
| lumo | support | `gamma2` | 0.26 | +0.026 | [-0.157, +0.212] | 121 / 57 | CI includes 0 |
| lumo | support | `P` | 0.54 | -0.054 | [-0.237, +0.128] | 121 / 57 | CI includes 0 |
| lumo | support | `D` | 0.20 | +0.020 | [-0.159, +0.198] | 121 / 57 | CI includes 0 |
| lumo | support | `F` | 0.51 | -0.051 | [-0.219, +0.119] | 121 / 57 | CI includes 0 |
| lumo | support | `A` | 0.40 | -0.040 | [-0.215, +0.138] | 121 / 57 | CI includes 0 |
| lumo | support | `S` | 0.22 | -0.022 | [-0.205, +0.162] | 121 / 57 | CI includes 0 |
| bautzen | committed | `gamma2` | 0.17 | -0.017 | [-0.277, +0.238] | 27 / 101 | CI includes 0 |
| bautzen | committed | `P` | 0.32 | -0.032 | [-0.281, +0.216] | 27 / 107 | CI includes 0 |
| bautzen | committed | `D` | 1.40 | -0.140 | [-0.370, +0.099] | 27 / 107 | CI includes 0 |
| bautzen | committed | `F` | 0.39 | +0.039 | [-0.184, +0.252] | 28 / 114 | CI includes 0 |
| bautzen | committed | `A` | 0.47 | -0.047 | [-0.277, +0.185] | 28 / 114 | CI includes 0 |
| bautzen | committed | `S` | 0.84 | +0.084 | [-0.147, +0.315] | 27 / 107 | CI includes 0 |
| bautzen | support | `gamma2` | 0.68 | -0.068 | [-0.319, +0.192] | 26 / 99 | CI includes 0 |
| bautzen | support | `P` | 0.11 | -0.011 | [-0.242, +0.222] | 27 / 108 | CI includes 0 |
| bautzen | support | `D` | 0.43 | +0.043 | [-0.210, +0.298] | 27 / 108 | CI includes 0 |
| bautzen | support | `F` | 0.55 | -0.055 | [-0.280, +0.167] | 28 / 114 | CI includes 0 |
| bautzen | support | `A` | 0.25 | +0.025 | [-0.211, +0.262] | 28 / 114 | CI includes 0 |
| bautzen | support | `S` | 0.83 | +0.083 | [-0.154, +0.318] | 27 / 108 | CI includes 0 |

- **lumo**: the committed EMV reference (`data/lumo/reference/lumo_echo_mask_vector.json`) is reproduced exactly for all five dimensions.
- **bautzen**: nothing to pin — the registry names no reference for this site.

## 4. What the chip size changed

Max pairwise, committed → support, with the winning pair on each side. A changed pair is the point: the maximum moved to a different comparison, not just to a different number.

| asset | channel | dimension | committed | support | Δ SDS | winning pair committed → support |
| --- | --- | --- | --- | --- | --- | --- |
| lumo | ASCENDING | `gamma2` | 1.62 | 4.88 | +3.25 | healthy – DAM 3 → healthy – DAM 6 |
| lumo | ASCENDING | `P` | 7.12 | 4.69 | -2.43 | healthy – DAM 6 → healthy – DAM 3 |
| lumo | ASCENDING | `D` | 3.40 | 3.35 | -0.05 | healthy – DAM 4 → healthy – DAM 4 |
| lumo | ASCENDING | `F` | 4.62 | 7.12 | +2.50 | healthy – DAM 3 → healthy – DAM 6 |
| lumo | ASCENDING | `A` | 2.50 | 5.94 | +3.44 | healthy – DAM 3 → healthy – DAM 6 |
| lumo | ASCENDING | `S` | 1.40 | 2.06 | +0.66 | healthy – DAM 4 → healthy – DAM 3 |
| lumo | DESCENDING | `gamma2` | 2.77 | 2.88 | +0.12 | healthy – DAM 3 → healthy – DAM 4 |
| lumo | DESCENDING | `P` | 3.15 | 2.50 | -0.65 | healthy – DAM 6 → healthy – DAM 3 |
| lumo | DESCENDING | `D` | 1.76 | 3.66 | +1.90 | healthy – DAM 4 → healthy – DAM 4 |
| lumo | DESCENDING | `F` | 2.44 | 3.08 | +0.64 | healthy – DAM 3 → healthy – DAM 4 |
| lumo | DESCENDING | `A` | 2.37 | 1.56 | -0.81 | healthy – DAM 3 → healthy – DAM 6 |
| lumo | DESCENDING | `S` | 0.66 | 2.27 | +1.61 | healthy – DAM 3 → healthy – DAM 4 |
| bautzen | ASC | `gamma2` | 1.28 | 0.89 | -0.40 | RS – DS1 → RS – DS1 |
| bautzen | ASC | `P` | 1.97 | 1.52 | -0.45 | RS – DS1 → RS – DS1 |
| bautzen | ASC | `D` | 3.01 | 1.22 | -1.79 | RS – DS2 → RS – DS2 |
| bautzen | ASC | `F` | 2.50 | 2.02 | -0.48 | RS – DS2 → RS – DS1 |
| bautzen | ASC | `A` | 0.96 | 1.13 | +0.17 | RS – DS1 → RS – DS1 |
| bautzen | ASC | `S` | 3.08 | 1.22 | -1.86 | RS – DS2 → RS – DS2 |
| bautzen | DESC | `gamma2` | 0.57 | 0.31 | -0.26 | RS – DS2 → RS – DS1 |
| bautzen | DESC | `P` | 3.03 | 2.29 | -0.74 | RS – DS1 → RS – DS2 |
| bautzen | DESC | `D` | 3.79 | 2.14 | -1.64 | RS – DS2 → RS – DS2 |
| bautzen | DESC | `F` | 3.50 | 2.14 | -1.36 | RS – DS2 → RS – DS2 |
| bautzen | DESC | `A` | 0.86 | 3.57 | +2.71 | RS – DS2 → RS – DS2 |
| bautzen | DESC | `S` | 2.64 | 2.00 | -0.64 | RS – DS2 → RS – DS2 |

## 5. What this does not claim

- SDS measures **separability on the given contrast** — two labelled groups on one record — never damage.
- A channel with few damaged overpasses gets a wide interval and the `n<5` flag; the score is printed, the flag is not optional.
- The two tables are two measurements, not a better and a worse instrument: the chip changes which cells enter the mask.
- State, season and orbit are confounded in these records; a per-channel SDS does not remove that.
