# Espoo mast (OSM 2233589804, event-free control) — two years, two orbits, no states

* table `data/espoo/espoo_ocv_channels.csv` (sha256 `0f167c108f7db8f5…`), 150 overpasses
* 2024-09-05 … 2026-09-07; orbits: ASCENDING 70, DESCENDING 80
* **no states**: the record carries no damage or structural label, so the figure has no state bands — the axis is time only
* **no chip**: Espoo has no mask layer; `A` is the pipeline's `coherence_masked_pixels`, and `D, F, S, P, gamma2_mask` do not exist here
* drawing rule: the Theil-Sen line is solid only where the ported stationarity flag says trend (block interval excludes zero *and* Mann-Kendall p < 0.05), dashed otherwise

## 1. Drift inside the two years

| dimension | orbit | n | slope / year | 95% CI | MK p | lag-1 | split-half delta | trend | serial flags | ties |
|---|---|---|---|---|---|---|---|---|---|---|
| `gamma2` | ASC afternoon | 70 | -0.0175 | [-0.0984, +0.0725] | 0.591 | +0.097 | -0.0057 | no | residual_autocorrelation | 0.01 at 0.484068996 |
| `gamma2` | DESC morning | 80 | +0.0071 | [-0.00127, +0.0172] | 0.088 | -0.066 | +0.14 | no | non_random_order | 0.01 at 0.006308362 |
| `A` | ASC afternoon | 70 | +0 | [-0.87, +1.08] | 0.850 | -0.12 | +0.011 | no | non_random_order | 0.17 at 4.0 |
| `A` | DESC morning | 80 | -1.43 | [-3.21, +1.52] | 0.310 | -0.2 | -0.08 | no | non_random_order | 0.09 at 3.0 |
| `displacement` | ASC afternoon | 70 | +0 | [+0, +0] | 0.135 | +0.011 | +0.074 | no | none | 0.64 at 0.0 **(degenerate)** |
| `displacement` | DESC morning | 80 | +0 | [-0.000651, +0.00115] | 0.933 | +0.063 | +0.11 | no | none | 0.24 at 0.0 |
| `brightness` | ASC afternoon | 70 | +0 | [+0, +0] | 0.138 | +0.045 | -0.074 | no | none | 0.64 at 1.0 **(degenerate)** |
| `brightness` | DESC morning | 80 | +0 | [-0.116, +0.0788] | 0.907 | +0.048 | -0.11 | no | none | 0.24 at 1.0 |
| `phase_coherence` | ASC afternoon | 39 | +0.000213 | [-0.0379, +0.0483] | 0.981 | -0.028 | +0.042 | no | none | 0.03 at 0.136927078 |
| `phase_coherence` | DESC morning | 41 | +0.0155 | [-0.0193, +0.0566] | 0.229 | -0.0071 | +0.22 | no | none | 0.02 at 0.129885374 |
| `phase_rms` | ASC afternoon | 39 | -0.0399 | [-1.61, +2.95] | 0.942 | +0.22 | -0.021 | no | none | 0.03 at 5.431609213 |
| `phase_rms` | DESC morning | 41 | +0.147 | [-1.84, +1.15] | 0.762 | +0.039 | -0.24 | no | none | 0.02 at 4.871973894 |

## 2. Year repeat — the same calendar month in both years

| dimension | orbit | months | median Δ (y2 − y1) | 95% CI | months + / − or the reason |
|---|---|---|---|---|---|
| `gamma2` | ASC afternoon | 12 | +0.0696 | [-0.0304, +0.143] | 8 / 4 |
| `gamma2` | DESC morning | 12 | +0.0156 | [-0.00145, +0.0165] | 8 / 4 |
| `A` | ASC afternoon | 12 | +0 | [-2, +1] | 4 / 8 |
| `A` | DESC morning | 12 | -3.5 | [-12, -0.5] | 2 / 10 |
| `displacement` | ASC afternoon | 12 | +0 | [+0, +0.000597] | 5 / 7 |
| `displacement` | DESC morning | 12 | +0.000719 | [-0.000987, +0.00189] | 7 / 5 |
| `brightness` | ASC afternoon | 12 | +0 | [-0.0751, +0] | 1 / 11 |
| `brightness` | DESC morning | 12 | -0.0985 | [-0.274, +0.0578] | 5 / 7 |
| `phase_coherence` | ASC afternoon | 2 | — | — | only 2 calendar month(s) carry a value in both campaign years (2024-09-05 … 2025-10-05) |
| `phase_coherence` | DESC morning | 2 | — | — | only 2 calendar month(s) carry a value in both campaign years (2024-09-05 … 2025-10-05) |
| `phase_rms` | ASC afternoon | 2 | — | — | only 2 calendar month(s) carry a value in both campaign years (2024-09-05 … 2025-10-05) |
| `phase_rms` | DESC morning | 2 | — | — | only 2 calendar month(s) carry a value in both campaign years (2024-09-05 … 2025-10-05) |

## 3. Orbit contrast, paired by date (no event in it)

| dimension | pairs | median ASC − DESC | 95% CI | share positive | drift of that difference |
|---|---|---|---|---|---|
| `gamma2` | 69 | +0.113 | [+0.0711, +0.224] | 0.93 | -0.0283 [-0.112, +0.0784] → stable |
| `A` | 69 | -4 | [-5, -1] | 0.28 | +1.35 [-4.06, +4.36] → stable |
| `displacement` | 69 | -0.00244 | [-0.00362, -0.000877] | 0.17 | +0.000977 [-0.00122, +0.0026] → stable |
| `brightness` | 69 | +0.261 | [+0.0854, +0.514] | 0.67 | -0.106 [-0.302, +0.135] → stable |
| `phase_coherence` | 39 | +0.00523 | [-0.00786, +0.0183] | 0.59 | -0.0173 [-0.0717, +0.0501] → stable |
| `phase_rms` | 39 | -0.565 | [-1.28, +0.253] | 0.38 | -0.269 [-2.78, +3.78] → stable |

## What this can and cannot say

* **Stability, not a signature.** There is no event here, so nothing in this file is an effect; it is the baseline the event sites lack.
* **No mask geometry.** `D, F, S, P` and `gamma2_mask` do not exist for Espoo (`data/espoo/asset_evidence.json`), so the stability statement is about coherence, phase, displacement and brightness. A chip-sized version needs a confirmed `base_width_m` first.
* **`phase_coherence` and `phase_rms_rad` cover half the rows** (80 of 150), so their `n` is smaller in every table.
* **`trend = no` is not `stable`**: the serial flags (autocorrelation, runs, Ljung-Box) are a different statement, and the split-half and the year-repeat columns are the second and third check. A dimension can pass one and fail another.
* **The solid line follows the ported flag rule**, which is stricter than the interval alone: it needs the block interval to exclude zero *and* Mann-Kendall p < 0.05 (`espoo_ocv_stats.stationarity_flags`).
* **A zero slope with a zero-width interval is not stability.** `brightness_ratio` is capped at 1.0 and `displacement_los_m` floored at 0.0; where that value fills at least half the rows the pairwise slopes are all zero and the bootstrap can only return zero. Those rows carry `(degenerate)` in the ties column: no power, not a result.

