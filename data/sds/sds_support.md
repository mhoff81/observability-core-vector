# SDS of the OCV dimensions on the support-sized chips

8 sites (bautzen, carola, cts, haltern, kdlo, lumo, morandi, ywf), recomputed on the **support-sized chip** (odd, anchor-centred, guard 30 m, reference ring 30–200 m) and beside it on the committed table. The metric is the project's own: `SDS = 10 · |Cliff's delta|`, computed with `code/lumo/lumo_ocv_stats.cliffs_delta`; the interval is a percentile bootstrap (10,000 resamples, seed 7) of **one comparison**.

## 0. Reading rules

1. **The primary analysis is per channel and per state pair.** The observation channel is the orbit; ASC and DESC see different geometries and are never pooled into the primary number.
2. **A maximum always travels with its winning pair** — the signed delta, that pair's interval and its two sample sizes.
3. **The interval belongs to one comparison.** For a maximum over several comparisons the report prints the *winning pair's* interval and says so; it is not an interval for the selection.
4. **`pooled_all_damage` is a compatibility pin and context**, not the primary evidence: it is the basis the committed LUMO reference used.

## 1. Primary analysis — per channel, per state pair

Support chip. ★ marks the winning comparison of its dimension (the one the compact table in §2 carries).

### Bautzen OpenLabs bridge

**ASC** — rows per class: DS1 46, DS2 12, RS 14

| dimension | state pair | SDS | band | delta | 95% CI | n1 / n2 | flags |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `A` | RS – DS1 ★ | 1.13 | very weak | -0.113 | [-0.446, +0.225] | 14 / 46 | CI includes 0 |
| `A` | RS – DS2 | 0.71 | very weak | +0.071 | [-0.393, +0.518] | 14 / 12 | CI includes 0 |
| `D` | RS – DS1 | 0.06 | very weak | +0.006 | [-0.340, +0.362] | 13 / 40 | CI includes 0 |
| `D` | RS – DS2 ★ | 1.22 | very weak | -0.122 | [-0.564, +0.353] | 13 / 12 | CI includes 0 |
| `F` | RS – DS1 ★ | 2.02 | weak | -0.202 | [-0.497, +0.102] | 14 / 46 | CI includes 0 |
| `F` | RS – DS2 | 0.71 | very weak | +0.071 | [-0.310, +0.470] | 14 / 12 | CI includes 0 |
| `P` | RS – DS1 ★ | 1.52 | very weak | +0.152 | [-0.194, +0.488] | 13 / 40 | CI includes 0 |
| `P` | RS – DS2 | 0.64 | very weak | -0.064 | [-0.513, +0.397] | 13 / 12 | CI includes 0 |
| `S` | RS – DS1 | 0.58 | very weak | -0.058 | [-0.387, +0.283] | 13 / 40 | CI includes 0 |
| `S` | RS – DS2 ★ | 1.22 | very weak | +0.122 | [-0.353, +0.577] | 13 / 12 | CI includes 0 |
| `gamma2` | RS – DS1 ★ | 0.89 | very weak | -0.089 | [-0.464, +0.292] | 12 / 32 | CI includes 0 |
| `gamma2` | RS – DS2 | 0.45 | very weak | -0.045 | [-0.530, +0.455] | 12 / 11 | CI includes 0 |


**DESC** — rows per class: DS1 46, DS2 10, RS 14

| dimension | state pair | SDS | band | delta | 95% CI | n1 / n2 | flags |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `A` | RS – DS1 | 1.52 | very weak | +0.152 | [-0.207, +0.495] | 14 / 46 | CI includes 0 |
| `A` | RS – DS2 ★ | 3.57 | weak | +0.357 | [-0.100, +0.757] | 14 / 10 | CI includes 0 |
| `D` | RS – DS1 | 1.23 | very weak | +0.123 | [-0.248, +0.472] | 14 / 46 | CI includes 0 |
| `D` | RS – DS2 ★ | 2.14 | weak | +0.214 | [-0.257, +0.671] | 14 / 10 | CI includes 0 |
| `F` | RS – DS1 | 0.43 | very weak | +0.043 | [-0.280, +0.363] | 14 / 46 | CI includes 0 |
| `F` | RS – DS2 ★ | 2.14 | weak | +0.214 | [-0.229, +0.614] | 14 / 10 | CI includes 0 |
| `P` | RS – DS1 | 1.49 | very weak | -0.149 | [-0.509, +0.225] | 14 / 46 | CI includes 0 |
| `P` | RS – DS2 ★ | 2.29 | weak | -0.229 | [-0.657, +0.243] | 14 / 10 | CI includes 0 |
| `S` | RS – DS1 | 1.68 | very weak | +0.168 | [-0.171, +0.509] | 14 / 46 | CI includes 0 |
| `S` | RS – DS2 ★ | 2.00 | weak | +0.200 | [-0.272, +0.629] | 14 / 10 | CI includes 0 |
| `gamma2` | RS – DS1 ★ | 0.31 | very weak | -0.031 | [-0.401, +0.342] | 14 / 46 | CI includes 0 |
| `gamma2` | RS – DS2 | 0.14 | very weak | -0.014 | [-0.471, +0.457] | 14 / 10 | CI includes 0 |


### Carolabrücke deck

**ASCENDING** — rows per class: post-collapse 166, pre-collapse (healthy) 124

| dimension | state pair | SDS | band | delta | 95% CI | n1 / n2 | flags |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `A` | pre-collapse (healthy) – post-collapse ★ | 0.92 | very weak | +0.092 | [-0.047, +0.226] | 124 / 166 | CI includes 0 |
| `D` | pre-collapse (healthy) – post-collapse ★ | 0.94 | very weak | +0.094 | [-0.047, +0.225] | 124 / 166 | CI includes 0 |
| `F` | pre-collapse (healthy) – post-collapse ★ | 0.76 | very weak | +0.076 | [-0.062, +0.209] | 124 / 166 | CI includes 0 |
| `P` | pre-collapse (healthy) – post-collapse ★ | 1.20 | very weak | -0.120 | [-0.259, +0.021] | 124 / 166 | CI includes 0 |
| `S` | pre-collapse (healthy) – post-collapse ★ | 0.96 | very weak | -0.096 | [-0.235, +0.043] | 124 / 166 | CI includes 0 |
| `gamma2` | pre-collapse (healthy) – post-collapse ★ | 0.25 | very weak | -0.025 | [-0.157, +0.110] | 124 / 166 | CI includes 0 |


**DESCENDING** — rows per class: post-collapse 86, pre-collapse (healthy) 60

| dimension | state pair | SDS | band | delta | 95% CI | n1 / n2 | flags |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `A` | pre-collapse (healthy) – post-collapse ★ | 6.00 | strong | -0.600 | [-0.733, -0.455] | 60 / 86 | – |
| `D` | pre-collapse (healthy) – post-collapse ★ | 5.74 | moderate | -0.574 | [-0.714, -0.422] | 60 / 86 | – |
| `F` | pre-collapse (healthy) – post-collapse ★ | 5.81 | moderate | -0.581 | [-0.717, -0.432] | 60 / 86 | – |
| `P` | pre-collapse (healthy) – post-collapse ★ | 2.90 | weak | -0.290 | [-0.465, -0.112] | 60 / 86 | – |
| `S` | pre-collapse (healthy) – post-collapse ★ | 4.21 | moderate | +0.421 | [+0.250, +0.581] | 60 / 86 | – |
| `gamma2` | pre-collapse (healthy) – post-collapse ★ | 5.50 | moderate | +0.550 | [+0.396, +0.694] | 60 / 86 | – |


### Champlain Towers South (target footprint)

**ASCENDING** — rows per class: post-collapse 42, pre-collapse (healthy) 136

| dimension | state pair | SDS | band | delta | 95% CI | n1 / n2 | flags |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `A` | pre-collapse (healthy) – post-collapse ★ | 1.18 | very weak | +0.118 | [-0.096, +0.332] | 136 / 42 | CI includes 0 |
| `D` | pre-collapse (healthy) – post-collapse ★ | 1.41 | very weak | +0.141 | [-0.075, +0.355] | 136 / 42 | CI includes 0 |
| `F` | pre-collapse (healthy) – post-collapse ★ | 1.18 | very weak | +0.118 | [-0.090, +0.328] | 136 / 42 | CI includes 0 |
| `P` | pre-collapse (healthy) – post-collapse ★ | 2.28 | weak | -0.228 | [-0.419, -0.029] | 136 / 42 | – |
| `S` | pre-collapse (healthy) – post-collapse ★ | 1.75 | very weak | -0.175 | [-0.370, +0.031] | 136 / 42 | CI includes 0 |
| `gamma2` | pre-collapse (healthy) – post-collapse ★ | 0.16 | very weak | +0.016 | [-0.178, +0.206] | 136 / 42 | CI includes 0 |


### Haltern AV9 wind turbine

**ASCENDING** — rows per class: post-collapse 5, pre-collapse (healthy) 44

| dimension | state pair | SDS | band | delta | 95% CI | n1 / n2 | flags |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `A` | pre-collapse (healthy) – post-collapse ★ | 1.32 | very weak | -0.132 | [-0.609, +0.400] | 44 / 5 | CI includes 0 |
| `D` | pre-collapse (healthy) – post-collapse ★ | 1.09 | very weak | -0.109 | [-0.591, +0.418] | 44 / 5 | CI includes 0 |
| `F` | pre-collapse (healthy) – post-collapse ★ | 2.73 | weak | +0.273 | [-0.345, +0.773] | 44 / 5 | CI includes 0 |
| `P` | pre-collapse (healthy) – post-collapse ★ | 1.23 | very weak | +0.123 | [-0.282, +0.500] | 44 / 5 | CI includes 0 |
| `S` | pre-collapse (healthy) – post-collapse ★ | 1.09 | very weak | +0.109 | [-0.282, +0.500] | 44 / 5 | CI includes 0 |
| `gamma2` | pre-collapse (healthy) – post-collapse ★ | 2.18 | weak | -0.218 | [-0.645, +0.309] | 44 / 5 | CI includes 0 |


**DESCENDING** — rows per class: post-collapse 6, pre-collapse (healthy) 52

| dimension | state pair | SDS | band | delta | 95% CI | n1 / n2 | flags |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `A` | pre-collapse (healthy) – post-collapse ★ | 4.04 | moderate | +0.404 | [-0.128, +0.833] | 52 / 6 | CI includes 0 |
| `D` | pre-collapse (healthy) – post-collapse ★ | 2.76 | weak | +0.276 | [-0.269, +0.744] | 52 / 6 | CI includes 0 |
| `F` | pre-collapse (healthy) – post-collapse ★ | 4.36 | moderate | +0.436 | [+0.038, +0.776] | 52 / 6 | – |
| `P` | pre-collapse (healthy) – post-collapse ★ | 1.41 | very weak | +0.141 | [-0.224, +0.500] | 52 / 6 | CI includes 0 |
| `S` | pre-collapse (healthy) – post-collapse ★ | 0.90 | very weak | +0.090 | [-0.449, +0.647] | 52 / 6 | CI includes 0 |
| `gamma2` | pre-collapse (healthy) – post-collapse ★ | 2.18 | weak | +0.218 | [-0.250, +0.654] | 52 / 6 | CI includes 0 |


### KDLO media tower

**DESCENDING** — rows per class: pre 16, rebuild 43

| dimension | state pair | SDS | band | delta | 95% CI | n1 / n2 | flags |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `A` | pre – rebuild ★ | 4.99 | moderate | +0.499 | [+0.173, +0.779] | 16 / 43 | – |
| `D` | pre – rebuild ★ | 2.01 | weak | -0.201 | [-0.532, +0.145] | 16 / 43 | CI includes 0 |
| `F` | pre – rebuild ★ | 1.21 | very weak | +0.121 | [-0.186, +0.423] | 16 / 43 | CI includes 0 |
| `P` | pre – rebuild ★ | 6.28 | strong | -0.628 | [-0.860, -0.337] | 16 / 43 | – |
| `S` | pre – rebuild ★ | 6.40 | strong | -0.640 | [-0.872, -0.363] | 16 / 43 | – |
| `gamma2` | pre – rebuild ★ | 0.47 | very weak | +0.047 | [-0.288, +0.384] | 16 / 43 | CI includes 0 |


### LUMO lattice tower (Hannover)

**ASCENDING** — rows per class: DAM 3 8, DAM 4 5, DAM 6 4, healthy 40

| dimension | state pair | SDS | band | delta | 95% CI | n1 / n2 | flags |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `A` | healthy – DAM 3 | 1.16 | very weak | +0.116 | [-0.272, +0.491] | 40 / 8 | CI includes 0 |
| `A` | healthy – DAM 4 | 1.60 | very weak | +0.160 | [-0.415, +0.670] | 40 / 5 | CI includes 0 |
| `A` | healthy – DAM 6 ★ | 5.94 | moderate | -0.594 | [-0.963, -0.100] | 40 / 4 | n<5 |
| `D` | healthy – DAM 3 | 1.03 | very weak | +0.103 | [-0.297, +0.487] | 40 / 8 | CI includes 0 |
| `D` | healthy – DAM 4 ★ | 3.35 | weak | +0.335 | [-0.135, +0.760] | 40 / 5 | CI includes 0 |
| `D` | healthy – DAM 6 | 1.94 | very weak | +0.194 | [-0.600, +0.975] | 40 / 4 | CI includes 0, n<5 |
| `F` | healthy – DAM 3 | 0.72 | very weak | +0.072 | [-0.344, +0.475] | 40 / 8 | CI includes 0 |
| `F` | healthy – DAM 4 | 1.25 | very weak | -0.125 | [-0.490, +0.225] | 40 / 5 | CI includes 0 |
| `F` | healthy – DAM 6 ★ | 7.12 | strong | -0.713 | [-0.950, -0.406] | 40 / 4 | n<5 |
| `P` | healthy – DAM 3 ★ | 4.69 | moderate | -0.469 | [-0.744, -0.150] | 40 / 8 | – |
| `P` | healthy – DAM 4 | 0.60 | very weak | +0.060 | [-0.450, +0.590] | 40 / 5 | CI includes 0 |
| `P` | healthy – DAM 6 | 3.12 | weak | +0.312 | [-0.163, +0.775] | 40 / 4 | CI includes 0, n<5 |
| `S` | healthy – DAM 3 ★ | 2.06 | weak | +0.206 | [-0.219, +0.606] | 40 / 8 | CI includes 0 |
| `S` | healthy – DAM 4 | 0.10 | very weak | -0.010 | [-0.440, +0.420] | 40 / 5 | CI includes 0 |
| `S` | healthy – DAM 6 | 1.50 | very weak | -0.150 | [-0.969, +0.700] | 40 / 4 | CI includes 0, n<5 |
| `gamma2` | healthy – DAM 3 | 0.06 | very weak | +0.006 | [-0.412, +0.419] | 40 / 8 | CI includes 0 |
| `gamma2` | healthy – DAM 4 | 1.90 | very weak | -0.190 | [-0.760, +0.400] | 40 / 5 | CI includes 0 |
| `gamma2` | healthy – DAM 6 ★ | 4.88 | moderate | +0.487 | [-0.013, +0.900] | 40 / 4 | CI includes 0, n<5 |


**DESCENDING** — rows per class: DAM 3 17, DAM 4 11, DAM 6 12, healthy 81

| dimension | state pair | SDS | band | delta | 95% CI | n1 / n2 | flags |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `A` | healthy – DAM 3 | 0.47 | very weak | +0.047 | [-0.248, +0.344] | 81 / 17 | CI includes 0 |
| `A` | healthy – DAM 4 | 0.79 | very weak | -0.079 | [-0.420, +0.259] | 81 / 11 | CI includes 0 |
| `A` | healthy – DAM 6 ★ | 1.56 | very weak | +0.156 | [-0.168, +0.461] | 81 / 12 | CI includes 0 |
| `D` | healthy – DAM 3 | 1.36 | very weak | -0.136 | [-0.432, +0.163] | 81 / 17 | CI includes 0 |
| `D` | healthy – DAM 4 ★ | 3.66 | weak | +0.366 | [+0.104, +0.602] | 81 / 11 | – |
| `D` | healthy – DAM 6 | 2.49 | weak | -0.249 | [-0.541, +0.064] | 81 / 12 | CI includes 0 |
| `F` | healthy – DAM 3 | 0.54 | very weak | +0.054 | [-0.217, +0.326] | 81 / 17 | CI includes 0 |
| `F` | healthy – DAM 4 ★ | 3.08 | weak | -0.308 | [-0.543, -0.043] | 81 / 11 | – |
| `F` | healthy – DAM 6 | 2.63 | weak | +0.263 | [-0.009, +0.514] | 81 / 12 | CI includes 0 |
| `P` | healthy – DAM 3 ★ | 2.50 | weak | -0.250 | [-0.521, +0.040] | 81 / 17 | CI includes 0 |
| `P` | healthy – DAM 4 | 1.34 | very weak | +0.134 | [-0.240, +0.493] | 81 / 11 | CI includes 0 |
| `P` | healthy – DAM 6 | 0.29 | very weak | +0.029 | [-0.342, +0.399] | 81 / 12 | CI includes 0 |
| `S` | healthy – DAM 3 | 0.12 | very weak | +0.012 | [-0.316, +0.339] | 81 / 17 | CI includes 0 |
| `S` | healthy – DAM 4 ★ | 2.27 | weak | -0.227 | [-0.586, +0.159] | 81 / 11 | CI includes 0 |
| `S` | healthy – DAM 6 | 2.25 | weak | +0.225 | [-0.072, +0.498] | 81 / 12 | CI includes 0 |
| `gamma2` | healthy – DAM 3 | 1.42 | very weak | -0.142 | [-0.454, +0.182] | 81 / 17 | CI includes 0 |
| `gamma2` | healthy – DAM 4 ★ | 2.88 | weak | +0.288 | [-0.071, +0.623] | 81 / 11 | CI includes 0 |
| `gamma2` | healthy – DAM 6 | 1.48 | very weak | -0.148 | [-0.486, +0.212] | 81 / 12 | CI includes 0 |


### Ponte Morandi deck

**ASCENDING** — rows per class: post-collapse 3, pre-collapse (healthy) 91

| dimension | state pair | SDS | band | delta | 95% CI | n1 / n2 | flags |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `A` | pre-collapse (healthy) – post-collapse ★ | 1.54 | very weak | +0.154 | [-0.121, +0.421] | 91 / 3 | CI includes 0, n<5 |
| `D` | pre-collapse (healthy) – post-collapse ★ | 7.29 | strong | -0.729 | [-1.000, -0.231] | 91 / 3 | n<5 |
| `F` | pre-collapse (healthy) – post-collapse ★ | 6.45 | strong | +0.645 | [+0.440, +0.813] | 91 / 3 | n<5 |
| `P` | pre-collapse (healthy) – post-collapse ★ | 1.06 | very weak | +0.106 | [-0.099, +0.304] | 91 / 3 | CI includes 0, n<5 |
| `S` | pre-collapse (healthy) – post-collapse ★ | 1.72 | very weak | -0.172 | [-0.436, +0.092] | 91 / 3 | CI includes 0, n<5 |
| `gamma2` | pre-collapse (healthy) – post-collapse ★ | 3.70 | weak | -0.370 | [-0.766, +0.143] | 91 / 3 | CI includes 0, n<5 |


**DESCENDING** — rows per class: post-collapse 6, pre-collapse (healthy) 107

| dimension | state pair | SDS | band | delta | 95% CI | n1 / n2 | flags |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `A` | pre-collapse (healthy) – post-collapse ★ | 2.04 | weak | -0.204 | [-0.565, +0.156] | 107 / 6 | CI includes 0 |
| `D` | pre-collapse (healthy) – post-collapse ★ | 2.24 | weak | -0.224 | [-0.592, +0.137] | 107 / 6 | CI includes 0 |
| `F` | pre-collapse (healthy) – post-collapse ★ | 2.24 | weak | -0.224 | [-0.576, +0.129] | 107 / 6 | CI includes 0 |
| `P` | pre-collapse (healthy) – post-collapse ★ | 3.80 | weak | +0.380 | [-0.103, +0.810] | 107 / 6 | CI includes 0 |
| `S` | pre-collapse (healthy) – post-collapse ★ | 3.02 | weak | +0.302 | [-0.153, +0.717] | 107 / 6 | CI includes 0 |
| `gamma2` | pre-collapse (healthy) – post-collapse ★ | 2.12 | weak | +0.212 | [-0.246, +0.660] | 107 / 6 | CI includes 0 |


### Yeongdeok Unit 21

**ASCENDING** — rows per class: post-collapse 0, pre-collapse (healthy) 1

| dimension | state pair | SDS | band | delta | 95% CI | n1 / n2 | flags |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `A` | pre-collapse (healthy) – post-collapse | – | – | – | – | 1 / 0 | n<5 |
| `D` | pre-collapse (healthy) – post-collapse | – | – | – | – | 1 / 0 | n<5 |
| `F` | pre-collapse (healthy) – post-collapse | – | – | – | – | 1 / 0 | n<5 |
| `P` | pre-collapse (healthy) – post-collapse | – | – | – | – | 1 / 0 | n<5 |
| `S` | pre-collapse (healthy) – post-collapse | – | – | – | – | 1 / 0 | n<5 |
| `gamma2` | pre-collapse (healthy) – post-collapse | – | – | – | – | 1 / 0 | n<5 |


**DESCENDING** — rows per class: post-collapse 2, pre-collapse (healthy) 11

| dimension | state pair | SDS | band | delta | 95% CI | n1 / n2 | flags |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `A` | pre-collapse (healthy) – post-collapse ★ | 4.55 | moderate | +0.455 | [-0.273, +1.000] | 11 / 2 | CI includes 0, n<5 |
| `D` | pre-collapse (healthy) – post-collapse ★ | 4.55 | moderate | +0.455 | [-0.273, +1.000] | 11 / 2 | CI includes 0, n<5 |
| `F` | pre-collapse (healthy) – post-collapse ★ | 0.00 | very weak | +0.000 | [-0.636, +0.636] | 11 / 2 | CI includes 0, n<5 |
| `P` | pre-collapse (healthy) – post-collapse ★ | 2.73 | weak | +0.273 | [-0.273, +0.818] | 11 / 2 | CI includes 0, n<5 |
| `S` | pre-collapse (healthy) – post-collapse ★ | 1.82 | very weak | +0.182 | [-0.636, +1.000] | 11 / 2 | CI includes 0, n<5 |
| `gamma2` | pre-collapse (healthy) – post-collapse ★ | 7.27 | strong | +0.727 | [+0.273, +1.000] | 11 / 2 | n<5 |


## 2. Compact summary — max pairwise

The paper table. Every row names the state pair that produced the maximum, its signed delta, **that pair's** interval, and the two sample sizes — a maximum without them would read as a general discrimination the number does not support.

| asset | channel | dimension | max SDS | winning pair | delta | 95% CI (winning pair) | n1 / n2 | flags |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| bautzen | ASC | `A` | 1.13 | RS – DS1 | -0.113 | [-0.446, +0.225] | 14 / 46 | CI includes 0 |
| bautzen | ASC | `D` | 1.22 | RS – DS2 | -0.122 | [-0.564, +0.353] | 13 / 12 | CI includes 0 |
| bautzen | ASC | `F` | 2.02 | RS – DS1 | -0.202 | [-0.497, +0.102] | 14 / 46 | CI includes 0 |
| bautzen | ASC | `P` | 1.52 | RS – DS1 | +0.152 | [-0.194, +0.488] | 13 / 40 | CI includes 0 |
| bautzen | ASC | `S` | 1.22 | RS – DS2 | +0.122 | [-0.353, +0.577] | 13 / 12 | CI includes 0 |
| bautzen | ASC | `gamma2` | 0.89 | RS – DS1 | -0.089 | [-0.464, +0.292] | 12 / 32 | CI includes 0 |
| bautzen | DESC | `A` | 3.57 | RS – DS2 | +0.357 | [-0.100, +0.757] | 14 / 10 | CI includes 0 |
| bautzen | DESC | `D` | 2.14 | RS – DS2 | +0.214 | [-0.257, +0.671] | 14 / 10 | CI includes 0 |
| bautzen | DESC | `F` | 2.14 | RS – DS2 | +0.214 | [-0.229, +0.614] | 14 / 10 | CI includes 0 |
| bautzen | DESC | `P` | 2.29 | RS – DS2 | -0.229 | [-0.657, +0.243] | 14 / 10 | CI includes 0 |
| bautzen | DESC | `S` | 2.00 | RS – DS2 | +0.200 | [-0.272, +0.629] | 14 / 10 | CI includes 0 |
| bautzen | DESC | `gamma2` | 0.31 | RS – DS1 | -0.031 | [-0.401, +0.342] | 14 / 46 | CI includes 0 |
| carola | ASCENDING | `A` | 0.92 | pre-collapse (healthy) – post-collapse | +0.092 | [-0.047, +0.226] | 124 / 166 | CI includes 0 |
| carola | ASCENDING | `D` | 0.94 | pre-collapse (healthy) – post-collapse | +0.094 | [-0.047, +0.225] | 124 / 166 | CI includes 0 |
| carola | ASCENDING | `F` | 0.76 | pre-collapse (healthy) – post-collapse | +0.076 | [-0.062, +0.209] | 124 / 166 | CI includes 0 |
| carola | ASCENDING | `P` | 1.20 | pre-collapse (healthy) – post-collapse | -0.120 | [-0.259, +0.021] | 124 / 166 | CI includes 0 |
| carola | ASCENDING | `S` | 0.96 | pre-collapse (healthy) – post-collapse | -0.096 | [-0.235, +0.043] | 124 / 166 | CI includes 0 |
| carola | ASCENDING | `gamma2` | 0.25 | pre-collapse (healthy) – post-collapse | -0.025 | [-0.157, +0.110] | 124 / 166 | CI includes 0 |
| carola | DESCENDING | `A` | 6.00 | pre-collapse (healthy) – post-collapse | -0.600 | [-0.733, -0.455] | 60 / 86 | – |
| carola | DESCENDING | `D` | 5.74 | pre-collapse (healthy) – post-collapse | -0.574 | [-0.714, -0.422] | 60 / 86 | – |
| carola | DESCENDING | `F` | 5.81 | pre-collapse (healthy) – post-collapse | -0.581 | [-0.717, -0.432] | 60 / 86 | – |
| carola | DESCENDING | `P` | 2.90 | pre-collapse (healthy) – post-collapse | -0.290 | [-0.465, -0.112] | 60 / 86 | – |
| carola | DESCENDING | `S` | 4.21 | pre-collapse (healthy) – post-collapse | +0.421 | [+0.250, +0.581] | 60 / 86 | – |
| carola | DESCENDING | `gamma2` | 5.50 | pre-collapse (healthy) – post-collapse | +0.550 | [+0.396, +0.694] | 60 / 86 | – |
| cts | ASCENDING | `A` | 1.18 | pre-collapse (healthy) – post-collapse | +0.118 | [-0.096, +0.332] | 136 / 42 | CI includes 0 |
| cts | ASCENDING | `D` | 1.41 | pre-collapse (healthy) – post-collapse | +0.141 | [-0.075, +0.355] | 136 / 42 | CI includes 0 |
| cts | ASCENDING | `F` | 1.18 | pre-collapse (healthy) – post-collapse | +0.118 | [-0.090, +0.328] | 136 / 42 | CI includes 0 |
| cts | ASCENDING | `P` | 2.28 | pre-collapse (healthy) – post-collapse | -0.228 | [-0.419, -0.029] | 136 / 42 | – |
| cts | ASCENDING | `S` | 1.75 | pre-collapse (healthy) – post-collapse | -0.175 | [-0.370, +0.031] | 136 / 42 | CI includes 0 |
| cts | ASCENDING | `gamma2` | 0.16 | pre-collapse (healthy) – post-collapse | +0.016 | [-0.178, +0.206] | 136 / 42 | CI includes 0 |
| haltern | ASCENDING | `A` | 1.32 | pre-collapse (healthy) – post-collapse | -0.132 | [-0.609, +0.400] | 44 / 5 | CI includes 0 |
| haltern | ASCENDING | `D` | 1.09 | pre-collapse (healthy) – post-collapse | -0.109 | [-0.591, +0.418] | 44 / 5 | CI includes 0 |
| haltern | ASCENDING | `F` | 2.73 | pre-collapse (healthy) – post-collapse | +0.273 | [-0.345, +0.773] | 44 / 5 | CI includes 0 |
| haltern | ASCENDING | `P` | 1.23 | pre-collapse (healthy) – post-collapse | +0.123 | [-0.282, +0.500] | 44 / 5 | CI includes 0 |
| haltern | ASCENDING | `S` | 1.09 | pre-collapse (healthy) – post-collapse | +0.109 | [-0.282, +0.500] | 44 / 5 | CI includes 0 |
| haltern | ASCENDING | `gamma2` | 2.18 | pre-collapse (healthy) – post-collapse | -0.218 | [-0.645, +0.309] | 44 / 5 | CI includes 0 |
| haltern | DESCENDING | `A` | 4.04 | pre-collapse (healthy) – post-collapse | +0.404 | [-0.128, +0.833] | 52 / 6 | CI includes 0 |
| haltern | DESCENDING | `D` | 2.76 | pre-collapse (healthy) – post-collapse | +0.276 | [-0.269, +0.744] | 52 / 6 | CI includes 0 |
| haltern | DESCENDING | `F` | 4.36 | pre-collapse (healthy) – post-collapse | +0.436 | [+0.038, +0.776] | 52 / 6 | – |
| haltern | DESCENDING | `P` | 1.41 | pre-collapse (healthy) – post-collapse | +0.141 | [-0.224, +0.500] | 52 / 6 | CI includes 0 |
| haltern | DESCENDING | `S` | 0.90 | pre-collapse (healthy) – post-collapse | +0.090 | [-0.449, +0.647] | 52 / 6 | CI includes 0 |
| haltern | DESCENDING | `gamma2` | 2.18 | pre-collapse (healthy) – post-collapse | +0.218 | [-0.250, +0.654] | 52 / 6 | CI includes 0 |
| kdlo | DESCENDING | `A` | 4.99 | pre – rebuild | +0.499 | [+0.173, +0.779] | 16 / 43 | – |
| kdlo | DESCENDING | `D` | 2.01 | pre – rebuild | -0.201 | [-0.532, +0.145] | 16 / 43 | CI includes 0 |
| kdlo | DESCENDING | `F` | 1.21 | pre – rebuild | +0.121 | [-0.186, +0.423] | 16 / 43 | CI includes 0 |
| kdlo | DESCENDING | `P` | 6.28 | pre – rebuild | -0.628 | [-0.860, -0.337] | 16 / 43 | – |
| kdlo | DESCENDING | `S` | 6.40 | pre – rebuild | -0.640 | [-0.872, -0.363] | 16 / 43 | – |
| kdlo | DESCENDING | `gamma2` | 0.47 | pre – rebuild | +0.047 | [-0.288, +0.384] | 16 / 43 | CI includes 0 |
| lumo | ASCENDING | `A` | 5.94 | healthy – DAM 6 | -0.594 | [-0.963, -0.100] | 40 / 4 | n<5 |
| lumo | ASCENDING | `D` | 3.35 | healthy – DAM 4 | +0.335 | [-0.135, +0.760] | 40 / 5 | CI includes 0 |
| lumo | ASCENDING | `F` | 7.12 | healthy – DAM 6 | -0.713 | [-0.950, -0.406] | 40 / 4 | n<5 |
| lumo | ASCENDING | `P` | 4.69 | healthy – DAM 3 | -0.469 | [-0.744, -0.150] | 40 / 8 | – |
| lumo | ASCENDING | `S` | 2.06 | healthy – DAM 3 | +0.206 | [-0.219, +0.606] | 40 / 8 | CI includes 0 |
| lumo | ASCENDING | `gamma2` | 4.88 | healthy – DAM 6 | +0.487 | [-0.013, +0.900] | 40 / 4 | CI includes 0, n<5 |
| lumo | DESCENDING | `A` | 1.56 | healthy – DAM 6 | +0.156 | [-0.168, +0.461] | 81 / 12 | CI includes 0 |
| lumo | DESCENDING | `D` | 3.66 | healthy – DAM 4 | +0.366 | [+0.104, +0.602] | 81 / 11 | – |
| lumo | DESCENDING | `F` | 3.08 | healthy – DAM 4 | -0.308 | [-0.543, -0.043] | 81 / 11 | – |
| lumo | DESCENDING | `P` | 2.50 | healthy – DAM 3 | -0.250 | [-0.521, +0.040] | 81 / 17 | CI includes 0 |
| lumo | DESCENDING | `S` | 2.27 | healthy – DAM 4 | -0.227 | [-0.586, +0.159] | 81 / 11 | CI includes 0 |
| lumo | DESCENDING | `gamma2` | 2.88 | healthy – DAM 4 | +0.288 | [-0.071, +0.623] | 81 / 11 | CI includes 0 |
| morandi | ASCENDING | `A` | 1.54 | pre-collapse (healthy) – post-collapse | +0.154 | [-0.121, +0.421] | 91 / 3 | CI includes 0, n<5 |
| morandi | ASCENDING | `D` | 7.29 | pre-collapse (healthy) – post-collapse | -0.729 | [-1.000, -0.231] | 91 / 3 | n<5 |
| morandi | ASCENDING | `F` | 6.45 | pre-collapse (healthy) – post-collapse | +0.645 | [+0.440, +0.813] | 91 / 3 | n<5 |
| morandi | ASCENDING | `P` | 1.06 | pre-collapse (healthy) – post-collapse | +0.106 | [-0.099, +0.304] | 91 / 3 | CI includes 0, n<5 |
| morandi | ASCENDING | `S` | 1.72 | pre-collapse (healthy) – post-collapse | -0.172 | [-0.436, +0.092] | 91 / 3 | CI includes 0, n<5 |
| morandi | ASCENDING | `gamma2` | 3.70 | pre-collapse (healthy) – post-collapse | -0.370 | [-0.766, +0.143] | 91 / 3 | CI includes 0, n<5 |
| morandi | DESCENDING | `A` | 2.04 | pre-collapse (healthy) – post-collapse | -0.204 | [-0.565, +0.156] | 107 / 6 | CI includes 0 |
| morandi | DESCENDING | `D` | 2.24 | pre-collapse (healthy) – post-collapse | -0.224 | [-0.592, +0.137] | 107 / 6 | CI includes 0 |
| morandi | DESCENDING | `F` | 2.24 | pre-collapse (healthy) – post-collapse | -0.224 | [-0.576, +0.129] | 107 / 6 | CI includes 0 |
| morandi | DESCENDING | `P` | 3.80 | pre-collapse (healthy) – post-collapse | +0.380 | [-0.103, +0.810] | 107 / 6 | CI includes 0 |
| morandi | DESCENDING | `S` | 3.02 | pre-collapse (healthy) – post-collapse | +0.302 | [-0.153, +0.717] | 107 / 6 | CI includes 0 |
| morandi | DESCENDING | `gamma2` | 2.12 | pre-collapse (healthy) – post-collapse | +0.212 | [-0.246, +0.660] | 107 / 6 | CI includes 0 |
| ywf | ASCENDING | `A` | – | – | – | – | – | – |
| ywf | ASCENDING | `D` | – | – | – | – | – | – |
| ywf | ASCENDING | `F` | – | – | – | – | – | – |
| ywf | ASCENDING | `P` | – | – | – | – | – | – |
| ywf | ASCENDING | `S` | – | – | – | – | – | – |
| ywf | ASCENDING | `gamma2` | – | – | – | – | – | – |
| ywf | DESCENDING | `A` | 4.55 | pre-collapse (healthy) – post-collapse | +0.455 | [-0.273, +1.000] | 11 / 2 | CI includes 0, n<5 |
| ywf | DESCENDING | `D` | 4.55 | pre-collapse (healthy) – post-collapse | +0.455 | [-0.273, +1.000] | 11 / 2 | CI includes 0, n<5 |
| ywf | DESCENDING | `F` | 0.00 | pre-collapse (healthy) – post-collapse | +0.000 | [-0.636, +0.636] | 11 / 2 | CI includes 0, n<5 |
| ywf | DESCENDING | `P` | 2.73 | pre-collapse (healthy) – post-collapse | +0.273 | [-0.273, +0.818] | 11 / 2 | CI includes 0, n<5 |
| ywf | DESCENDING | `S` | 1.82 | pre-collapse (healthy) – post-collapse | +0.182 | [-0.636, +1.000] | 11 / 2 | CI includes 0, n<5 |
| ywf | DESCENDING | `gamma2` | 7.27 | pre-collapse (healthy) – post-collapse | +0.727 | [+0.273, +1.000] | 11 / 2 | n<5 |

## 3. Compatibility pin and context — pooled over both channels

**Not the primary evidence.** This is the aggregation the committed LUMO reference used (`sds_pooled_all_damaged`), so it is reproduced here as a pin and reported for context. Pooling ASC and DESC mixes two observation geometries.

| asset | table | dimension | SDS | delta | 95% CI | n1 / n2 | flags |
| --- | --- | --- | --- | --- | --- | --- | --- |
| bautzen | committed | `A` | 0.47 | -0.047 | [-0.277, +0.185] | 28 / 114 | CI includes 0 |
| bautzen | committed | `D` | 1.40 | -0.140 | [-0.370, +0.099] | 27 / 107 | CI includes 0 |
| bautzen | committed | `F` | 0.39 | +0.039 | [-0.184, +0.252] | 28 / 114 | CI includes 0 |
| bautzen | committed | `P` | 0.32 | -0.032 | [-0.281, +0.216] | 27 / 107 | CI includes 0 |
| bautzen | committed | `S` | 0.84 | +0.084 | [-0.147, +0.315] | 27 / 107 | CI includes 0 |
| bautzen | committed | `gamma2` | 0.17 | -0.017 | [-0.277, +0.238] | 27 / 101 | CI includes 0 |
| bautzen | support | `A` | 0.25 | +0.025 | [-0.211, +0.262] | 28 / 114 | CI includes 0 |
| bautzen | support | `D` | 0.43 | +0.043 | [-0.210, +0.298] | 27 / 108 | CI includes 0 |
| bautzen | support | `F` | 0.55 | -0.055 | [-0.280, +0.167] | 28 / 114 | CI includes 0 |
| bautzen | support | `P` | 0.11 | -0.011 | [-0.242, +0.222] | 27 / 108 | CI includes 0 |
| bautzen | support | `S` | 0.83 | +0.083 | [-0.154, +0.318] | 27 / 108 | CI includes 0 |
| bautzen | support | `gamma2` | 0.68 | -0.068 | [-0.319, +0.192] | 26 / 99 | CI includes 0 |
| carola | committed | `A` | 1.00 | -0.100 | [-0.155, -0.047] | 896 / 821 | – |
| carola | committed | `D` | 0.63 | +0.063 | [+0.010, +0.119] | 896 / 821 | – |
| carola | committed | `F` | 0.74 | -0.074 | [-0.129, -0.021] | 896 / 821 | – |
| carola | committed | `P` | 6.61 | -0.661 | [-0.700, -0.618] | 896 / 821 | – |
| carola | committed | `S` | 1.14 | -0.114 | [-0.169, -0.059] | 896 / 821 | – |
| carola | committed | `gamma2` | 1.71 | +0.171 | [+0.117, +0.226] | 896 / 821 | – |
| carola | support | `A` | 1.63 | -0.163 | [-0.270, -0.054] | 184 / 252 | – |
| carola | support | `D` | 1.50 | -0.150 | [-0.257, -0.041] | 184 / 252 | – |
| carola | support | `F` | 1.27 | -0.127 | [-0.235, -0.019] | 184 / 252 | – |
| carola | support | `P` | 2.15 | -0.215 | [-0.321, -0.110] | 184 / 252 | – |
| carola | support | `S` | 1.23 | +0.123 | [+0.013, +0.231] | 184 / 252 | – |
| carola | support | `gamma2` | 1.70 | +0.170 | [+0.059, +0.277] | 184 / 252 | – |
| cts | committed | `A` | 1.67 | +0.167 | [+0.057, +0.275] | 543 / 168 | – |
| cts | committed | `D` | 2.76 | +0.276 | [+0.181, +0.367] | 543 / 168 | – |
| cts | committed | `F` | 1.53 | +0.153 | [+0.044, +0.262] | 543 / 168 | – |
| cts | committed | `P` | 1.68 | -0.168 | [-0.271, -0.063] | 543 / 168 | – |
| cts | committed | `S` | 1.18 | -0.118 | [-0.217, -0.018] | 543 / 168 | – |
| cts | committed | `gamma2` | 0.93 | -0.093 | [-0.190, +0.008] | 543 / 168 | CI includes 0 |
| cts | support | `A` | 1.18 | +0.118 | [-0.096, +0.332] | 136 / 42 | CI includes 0 |
| cts | support | `D` | 1.41 | +0.141 | [-0.075, +0.355] | 136 / 42 | CI includes 0 |
| cts | support | `F` | 1.18 | +0.118 | [-0.090, +0.328] | 136 / 42 | CI includes 0 |
| cts | support | `P` | 2.28 | -0.228 | [-0.419, -0.029] | 136 / 42 | – |
| cts | support | `S` | 1.75 | -0.175 | [-0.370, +0.031] | 136 / 42 | CI includes 0 |
| cts | support | `gamma2` | 0.16 | +0.016 | [-0.178, +0.206] | 136 / 42 | CI includes 0 |
| haltern | committed | `A` | 0.89 | +0.089 | [-0.661, +0.839] | 62 / 4 | CI includes 0, n<5 |
| haltern | committed | `D` | 0.97 | -0.097 | [-0.653, +0.492] | 62 / 4 | CI includes 0, n<5 |
| haltern | committed | `F` | 3.55 | +0.355 | [-0.266, +0.871] | 62 / 4 | CI includes 0, n<5 |
| haltern | committed | `P` | 3.67 | +0.367 | [-0.117, +0.758] | 62 / 4 | CI includes 0, n<5 |
| haltern | committed | `S` | 3.63 | +0.363 | [-0.286, +0.823] | 62 / 4 | CI includes 0, n<5 |
| haltern | committed | `gamma2` | 3.15 | -0.315 | [-0.677, +0.073] | 62 / 4 | CI includes 0, n<5 |
| haltern | support | `A` | 1.38 | +0.138 | [-0.276, +0.552] | 96 / 11 | CI includes 0 |
| haltern | support | `D` | 0.83 | +0.083 | [-0.326, +0.494] | 96 / 11 | CI includes 0 |
| haltern | support | `F` | 3.67 | +0.367 | [+0.032, +0.670] | 96 / 11 | – |
| haltern | support | `P` | 1.24 | +0.124 | [-0.188, +0.443] | 96 / 11 | CI includes 0 |
| haltern | support | `S` | 1.23 | +0.123 | [-0.239, +0.472] | 96 / 11 | CI includes 0 |
| haltern | support | `gamma2` | 0.34 | +0.034 | [-0.316, +0.388] | 96 / 11 | CI includes 0 |
| kdlo | committed | `A` | 2.81 | -0.281 | [-0.570, +0.029] | 16 / 43 | CI includes 0 |
| kdlo | committed | `D` | 0.70 | -0.070 | [-0.360, +0.233] | 16 / 43 | CI includes 0 |
| kdlo | committed | `F` | 2.67 | -0.267 | [-0.555, +0.038] | 16 / 43 | CI includes 0 |
| kdlo | committed | `P` | 2.35 | +0.235 | [-0.087, +0.544] | 16 / 43 | CI includes 0 |
| kdlo | committed | `S` | 1.45 | -0.145 | [-0.471, +0.180] | 16 / 43 | CI includes 0 |
| kdlo | committed | `gamma2` | 5.29 | +0.529 | [+0.262, +0.762] | 16 / 43 | – |
| kdlo | support | `A` | 4.99 | +0.499 | [+0.173, +0.779] | 16 / 43 | – |
| kdlo | support | `D` | 2.01 | -0.201 | [-0.532, +0.145] | 16 / 43 | CI includes 0 |
| kdlo | support | `F` | 1.21 | +0.121 | [-0.186, +0.423] | 16 / 43 | CI includes 0 |
| kdlo | support | `P` | 6.28 | -0.628 | [-0.860, -0.337] | 16 / 43 | – |
| kdlo | support | `S` | 6.40 | -0.640 | [-0.872, -0.363] | 16 / 43 | – |
| kdlo | support | `gamma2` | 0.47 | +0.047 | [-0.288, +0.384] | 16 / 43 | CI includes 0 |
| lumo | committed | `A` | 0.81 | -0.081 | [-0.255, +0.090] | 121 / 57 | CI includes 0 |
| lumo | committed | `D` | 0.86 | -0.086 | [-0.272, +0.098] | 121 / 57 | CI includes 0 |
| lumo | committed | `F` | 0.67 | -0.067 | [-0.244, +0.107] | 121 / 57 | CI includes 0 |
| lumo | committed | `P` | 1.05 | -0.105 | [-0.292, +0.089] | 120 / 57 | CI includes 0 |
| lumo | committed | `S` | 0.10 | -0.010 | [-0.196, +0.172] | 121 / 57 | CI includes 0 |
| lumo | committed | `gamma2` | 0.10 | +0.010 | [-0.165, +0.185] | 121 / 57 | CI includes 0 |
| lumo | support | `A` | 0.40 | -0.040 | [-0.215, +0.138] | 121 / 57 | CI includes 0 |
| lumo | support | `D` | 0.20 | +0.020 | [-0.159, +0.198] | 121 / 57 | CI includes 0 |
| lumo | support | `F` | 0.51 | -0.051 | [-0.219, +0.119] | 121 / 57 | CI includes 0 |
| lumo | support | `P` | 0.54 | -0.054 | [-0.237, +0.128] | 121 / 57 | CI includes 0 |
| lumo | support | `S` | 0.22 | -0.022 | [-0.205, +0.162] | 121 / 57 | CI includes 0 |
| lumo | support | `gamma2` | 0.26 | +0.026 | [-0.157, +0.212] | 121 / 57 | CI includes 0 |
| morandi | committed | `A` | 0.42 | +0.042 | [-0.336, +0.398] | 221 / 11 | CI includes 0 |
| morandi | committed | `D` | 0.14 | +0.014 | [-0.334, +0.393] | 221 / 11 | CI includes 0 |
| morandi | committed | `F` | 0.87 | +0.087 | [-0.342, +0.486] | 221 / 11 | CI includes 0 |
| morandi | committed | `P` | 1.30 | -0.130 | [-0.483, +0.221] | 221 / 11 | CI includes 0 |
| morandi | committed | `S` | 0.44 | +0.044 | [-0.358, +0.408] | 221 / 11 | CI includes 0 |
| morandi | committed | `gamma2` | 2.18 | +0.218 | [-0.164, +0.582] | 221 / 11 | CI includes 0 |
| morandi | support | `A` | 0.71 | +0.071 | [-0.230, +0.369] | 198 / 9 | CI includes 0 |
| morandi | support | `D` | 4.26 | -0.426 | [-0.722, -0.122] | 198 / 9 | – |
| morandi | support | `F` | 1.50 | +0.150 | [-0.089, +0.390] | 198 / 9 | CI includes 0 |
| morandi | support | `P` | 0.77 | +0.077 | [-0.257, +0.396] | 198 / 9 | CI includes 0 |
| morandi | support | `S` | 1.96 | +0.196 | [-0.218, +0.598] | 198 / 9 | CI includes 0 |
| morandi | support | `gamma2` | 0.84 | -0.084 | [-0.383, +0.208] | 198 / 9 | CI includes 0 |
| ywf | committed | `A` | 8.12 | +0.812 | – | 16 / 1 | n<5 |
| ywf | committed | `D` | 7.50 | -0.750 | – | 16 / 1 | n<5 |
| ywf | committed | `F` | 9.38 | +0.938 | – | 16 / 1 | n<5 |
| ywf | committed | `P` | 2.50 | +0.250 | – | 16 / 1 | n<5 |
| ywf | committed | `S` | 1.25 | -0.125 | – | 16 / 1 | n<5 |
| ywf | committed | `gamma2` | 3.75 | -0.375 | – | 16 / 1 | n<5 |
| ywf | support | `A` | 5.00 | +0.500 | [-0.167, +1.000] | 12 / 2 | CI includes 0, n<5 |
| ywf | support | `D` | 5.00 | +0.500 | [-0.167, +1.000] | 12 / 2 | CI includes 0, n<5 |
| ywf | support | `F` | 0.83 | +0.083 | [-0.500, +0.667] | 12 / 2 | CI includes 0, n<5 |
| ywf | support | `P` | 3.33 | +0.333 | [-0.167, +0.833] | 12 / 2 | CI includes 0, n<5 |
| ywf | support | `S` | 2.50 | +0.250 | [-0.500, +1.000] | 12 / 2 | CI includes 0, n<5 |
| ywf | support | `gamma2` | 7.50 | +0.750 | [+0.333, +1.000] | 12 / 2 | n<5 |

- **bautzen**: nothing to pin — the registry names no reference for this site.
- **carola**: nothing to pin — the registry names no reference for this site.
- **cts**: nothing to pin — the registry names no reference for this site.
- **haltern**: nothing to pin — the registry names no reference for this site.
- **kdlo**: nothing to pin — the registry names no reference for this site.
- **lumo**: the committed EMV reference (`data/lumo/reference/lumo_echo_mask_vector.json`) is reproduced exactly for all five dimensions.
- **morandi**: nothing to pin — the registry names no reference for this site.
- **ywf**: nothing to pin — the registry names no reference for this site.

## 4. What the chip size changed

Max pairwise, committed → support, with the winning pair on each side. A changed pair is the point: the maximum moved to a different comparison, not just to a different number.

| asset | channel | dimension | committed | support | Δ SDS | winning pair committed → support |
| --- | --- | --- | --- | --- | --- | --- |
| bautzen | ASC | `A` | 0.96 | 1.13 | +0.17 | RS – DS1 → RS – DS1 |
| bautzen | ASC | `D` | 3.01 | 1.22 | -1.79 | RS – DS2 → RS – DS2 |
| bautzen | ASC | `F` | 2.50 | 2.02 | -0.48 | RS – DS2 → RS – DS1 |
| bautzen | ASC | `P` | 1.97 | 1.52 | -0.45 | RS – DS1 → RS – DS1 |
| bautzen | ASC | `S` | 3.08 | 1.22 | -1.86 | RS – DS2 → RS – DS2 |
| bautzen | ASC | `gamma2` | 1.28 | 0.89 | -0.40 | RS – DS1 → RS – DS1 |
| bautzen | DESC | `A` | 0.86 | 3.57 | +2.71 | RS – DS2 → RS – DS2 |
| bautzen | DESC | `D` | 3.79 | 2.14 | -1.64 | RS – DS2 → RS – DS2 |
| bautzen | DESC | `F` | 3.50 | 2.14 | -1.36 | RS – DS2 → RS – DS2 |
| bautzen | DESC | `P` | 3.03 | 2.29 | -0.74 | RS – DS1 → RS – DS2 |
| bautzen | DESC | `S` | 2.64 | 2.00 | -0.64 | RS – DS2 → RS – DS2 |
| bautzen | DESC | `gamma2` | 0.57 | 0.31 | -0.26 | RS – DS2 → RS – DS1 |
| carola | ASCENDING | `A` | 1.38 | 0.92 | -0.46 | pre-collapse (healthy) – post-collapse → pre-collapse (healthy) – post-collapse |
| carola | ASCENDING | `D` | 0.91 | 0.94 | +0.03 | pre-collapse (healthy) – post-collapse → pre-collapse (healthy) – post-collapse |
| carola | ASCENDING | `F` | 1.08 | 0.76 | -0.33 | pre-collapse (healthy) – post-collapse → pre-collapse (healthy) – post-collapse |
| carola | ASCENDING | `P` | 6.10 | 1.20 | -4.90 | pre-collapse (healthy) – post-collapse → pre-collapse (healthy) – post-collapse |
| carola | ASCENDING | `S` | 0.79 | 0.96 | +0.17 | pre-collapse (healthy) – post-collapse → pre-collapse (healthy) – post-collapse |
| carola | ASCENDING | `gamma2` | 1.84 | 0.25 | -1.58 | pre-collapse (healthy) – post-collapse → pre-collapse (healthy) – post-collapse |
| carola | DESCENDING | `A` | 0.10 | 6.00 | +5.90 | pre-collapse (healthy) – post-collapse → pre-collapse (healthy) – post-collapse |
| carola | DESCENDING | `D` | 1.11 | 5.74 | +4.64 | pre-collapse (healthy) – post-collapse → pre-collapse (healthy) – post-collapse |
| carola | DESCENDING | `F` | 0.64 | 5.81 | +5.16 | pre-collapse (healthy) – post-collapse → pre-collapse (healthy) – post-collapse |
| carola | DESCENDING | `P` | 7.70 | 2.90 | -4.80 | pre-collapse (healthy) – post-collapse → pre-collapse (healthy) – post-collapse |
| carola | DESCENDING | `S` | 2.54 | 4.21 | +1.67 | pre-collapse (healthy) – post-collapse → pre-collapse (healthy) – post-collapse |
| carola | DESCENDING | `gamma2` | 1.72 | 5.50 | +3.79 | pre-collapse (healthy) – post-collapse → pre-collapse (healthy) – post-collapse |
| cts | ASCENDING | `A` | 1.67 | 1.18 | -0.49 | pre-collapse (healthy) – post-collapse → pre-collapse (healthy) – post-collapse |
| cts | ASCENDING | `D` | 2.76 | 1.41 | -1.35 | pre-collapse (healthy) – post-collapse → pre-collapse (healthy) – post-collapse |
| cts | ASCENDING | `F` | 1.53 | 1.18 | -0.35 | pre-collapse (healthy) – post-collapse → pre-collapse (healthy) – post-collapse |
| cts | ASCENDING | `P` | 1.68 | 2.28 | +0.60 | pre-collapse (healthy) – post-collapse → pre-collapse (healthy) – post-collapse |
| cts | ASCENDING | `S` | 1.18 | 1.75 | +0.56 | pre-collapse (healthy) – post-collapse → pre-collapse (healthy) – post-collapse |
| cts | ASCENDING | `gamma2` | 0.93 | 0.16 | -0.77 | pre-collapse (healthy) – post-collapse → pre-collapse (healthy) – post-collapse |
| haltern | ASCENDING | `A` | 1.77 | 1.32 | -0.46 | pre-collapse (healthy) – post-collapse → pre-collapse (healthy) – post-collapse |
| haltern | ASCENDING | `D` | 3.06 | 1.09 | -1.97 | pre-collapse (healthy) – post-collapse → pre-collapse (healthy) – post-collapse |
| haltern | ASCENDING | `F` | 6.45 | 2.73 | -3.72 | pre-collapse (healthy) – post-collapse → pre-collapse (healthy) – post-collapse |
| haltern | ASCENDING | `P` | 1.61 | 1.23 | -0.39 | pre-collapse (healthy) – post-collapse → pre-collapse (healthy) – post-collapse |
| haltern | ASCENDING | `S` | 7.74 | 1.09 | -6.65 | pre-collapse (healthy) – post-collapse → pre-collapse (healthy) – post-collapse |
| haltern | ASCENDING | `gamma2` | 1.94 | 2.18 | +0.25 | pre-collapse (healthy) – post-collapse → pre-collapse (healthy) – post-collapse |
| haltern | DESCENDING | `A` | 0.00 | 4.04 | +4.04 | pre-collapse (healthy) – post-collapse → pre-collapse (healthy) – post-collapse |
| haltern | DESCENDING | `D` | 1.29 | 2.76 | +1.47 | pre-collapse (healthy) – post-collapse → pre-collapse (healthy) – post-collapse |
| haltern | DESCENDING | `F` | 0.97 | 4.36 | +3.39 | pre-collapse (healthy) – post-collapse → pre-collapse (healthy) – post-collapse |
| haltern | DESCENDING | `P` | 6.61 | 1.41 | -5.20 | pre-collapse (healthy) – post-collapse → pre-collapse (healthy) – post-collapse |
| haltern | DESCENDING | `S` | 0.65 | 0.90 | +0.25 | pre-collapse (healthy) – post-collapse → pre-collapse (healthy) – post-collapse |
| haltern | DESCENDING | `gamma2` | 3.87 | 2.18 | -1.69 | pre-collapse (healthy) – post-collapse → pre-collapse (healthy) – post-collapse |
| kdlo | DESCENDING | `A` | 2.81 | 4.99 | +2.18 | pre – rebuild → pre – rebuild |
| kdlo | DESCENDING | `D` | 0.70 | 2.01 | +1.31 | pre – rebuild → pre – rebuild |
| kdlo | DESCENDING | `F` | 2.67 | 1.21 | -1.47 | pre – rebuild → pre – rebuild |
| kdlo | DESCENDING | `P` | 2.35 | 6.28 | +3.92 | pre – rebuild → pre – rebuild |
| kdlo | DESCENDING | `S` | 1.45 | 6.40 | +4.94 | pre – rebuild → pre – rebuild |
| kdlo | DESCENDING | `gamma2` | 5.29 | 0.47 | -4.83 | pre – rebuild → pre – rebuild |
| lumo | ASCENDING | `A` | 2.50 | 5.94 | +3.44 | healthy – DAM 3 → healthy – DAM 6 |
| lumo | ASCENDING | `D` | 3.40 | 3.35 | -0.05 | healthy – DAM 4 → healthy – DAM 4 |
| lumo | ASCENDING | `F` | 4.62 | 7.12 | +2.50 | healthy – DAM 3 → healthy – DAM 6 |
| lumo | ASCENDING | `P` | 7.12 | 4.69 | -2.43 | healthy – DAM 6 → healthy – DAM 3 |
| lumo | ASCENDING | `S` | 1.40 | 2.06 | +0.66 | healthy – DAM 4 → healthy – DAM 3 |
| lumo | ASCENDING | `gamma2` | 1.62 | 4.88 | +3.25 | healthy – DAM 3 → healthy – DAM 6 |
| lumo | DESCENDING | `A` | 2.37 | 1.56 | -0.81 | healthy – DAM 3 → healthy – DAM 6 |
| lumo | DESCENDING | `D` | 1.76 | 3.66 | +1.90 | healthy – DAM 4 → healthy – DAM 4 |
| lumo | DESCENDING | `F` | 2.44 | 3.08 | +0.64 | healthy – DAM 3 → healthy – DAM 4 |
| lumo | DESCENDING | `P` | 3.15 | 2.50 | -0.65 | healthy – DAM 6 → healthy – DAM 3 |
| lumo | DESCENDING | `S` | 0.66 | 2.27 | +1.61 | healthy – DAM 3 → healthy – DAM 4 |
| lumo | DESCENDING | `gamma2` | 2.77 | 2.88 | +0.12 | healthy – DAM 3 → healthy – DAM 4 |
| morandi | ASCENDING | `A` | 3.81 | 1.54 | -2.27 | pre-collapse (healthy) – post-collapse → pre-collapse (healthy) – post-collapse |
| morandi | ASCENDING | `D` | 5.75 | 7.29 | +1.53 | pre-collapse (healthy) – post-collapse → pre-collapse (healthy) – post-collapse |
| morandi | ASCENDING | `F` | 6.09 | 6.45 | +0.36 | pre-collapse (healthy) – post-collapse → pre-collapse (healthy) – post-collapse |
| morandi | ASCENDING | `P` | 2.42 | 1.06 | -1.36 | pre-collapse (healthy) – post-collapse → pre-collapse (healthy) – post-collapse |
| morandi | ASCENDING | `S` | 7.12 | 1.72 | -5.40 | pre-collapse (healthy) – post-collapse → pre-collapse (healthy) – post-collapse |
| morandi | ASCENDING | `gamma2` | 1.02 | 3.70 | +2.68 | pre-collapse (healthy) – post-collapse → pre-collapse (healthy) – post-collapse |
| morandi | DESCENDING | `A` | 1.36 | 2.04 | +0.69 | pre-collapse (healthy) – post-collapse → pre-collapse (healthy) – post-collapse |
| morandi | DESCENDING | `D` | 3.15 | 2.24 | -0.90 | pre-collapse (healthy) – post-collapse → pre-collapse (healthy) – post-collapse |
| morandi | DESCENDING | `F` | 2.06 | 2.24 | +0.19 | pre-collapse (healthy) – post-collapse → pre-collapse (healthy) – post-collapse |
| morandi | DESCENDING | `P` | 0.47 | 3.80 | +3.33 | pre-collapse (healthy) – post-collapse → pre-collapse (healthy) – post-collapse |
| morandi | DESCENDING | `S` | 2.91 | 3.02 | +0.11 | pre-collapse (healthy) – post-collapse → pre-collapse (healthy) – post-collapse |
| morandi | DESCENDING | `gamma2` | 4.08 | 2.12 | -1.96 | pre-collapse (healthy) – post-collapse → pre-collapse (healthy) – post-collapse |
| ywf | ASCENDING | `A` | 8.89 | – | – | pre-collapse (healthy) – post-collapse → – |
| ywf | ASCENDING | `D` | 10.00 | – | – | pre-collapse (healthy) – post-collapse → – |
| ywf | ASCENDING | `F` | 10.00 | – | – | pre-collapse (healthy) – post-collapse → – |
| ywf | ASCENDING | `P` | 7.78 | – | – | pre-collapse (healthy) – post-collapse → – |
| ywf | ASCENDING | `S` | 1.11 | – | – | pre-collapse (healthy) – post-collapse → – |
| ywf | ASCENDING | `gamma2` | 1.11 | – | – | pre-collapse (healthy) – post-collapse → – |
| ywf | DESCENDING | `A` | – | 4.55 | – | – → pre-collapse (healthy) – post-collapse |
| ywf | DESCENDING | `D` | – | 4.55 | – | – → pre-collapse (healthy) – post-collapse |
| ywf | DESCENDING | `F` | – | 0.00 | – | – → pre-collapse (healthy) – post-collapse |
| ywf | DESCENDING | `P` | – | 2.73 | – | – → pre-collapse (healthy) – post-collapse |
| ywf | DESCENDING | `S` | – | 1.82 | – | – → pre-collapse (healthy) – post-collapse |
| ywf | DESCENDING | `gamma2` | – | 7.27 | – | – → pre-collapse (healthy) – post-collapse |

## 5. What this does not claim

- SDS measures **separability on the given contrast** — two labelled groups on one record — never damage.
- A channel with few damaged overpasses gets a wide interval and the `n<5` flag; the score is printed, the flag is not optional.
- The two tables are two measurements, not a better and a worse instrument: the chip changes which cells enter the mask.
- State, season and orbit are confounded in these records; a per-channel SDS does not remove that.
