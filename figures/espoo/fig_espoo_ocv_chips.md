# Espoo mast (OSM 2233589804) — support chips — the OCV vector of the support chips

* table `data/espoo/support/espoo_ocv_channels.csv` (sha256 `82a8cc04208e5105…`), 75 chips
* 2024-09-05 … 2026-09-07; orbits: ASCENDING 31, DESCENDING 44; **3 shared dates**
* **second instrument**: these are chip values, not the pipeline columns of `fig_espoo_ocv_orbits.py` (r(log A, log `coherence_masked_pixels`) = −0.24, r(`gamma2_mask`, `coherence_gamma2`) = 0.13)
* **two tracks**: ASC = relative orbit 160, DESC = 51; the other two tracks place the mast on the far range edge, where the 200 m ring runs past the swath
* **no pairing**: the tracks have three shared dates, so the orbit contrast below is a two-sample comparison
* **tautology**: `D = A / bbox_area` exactly, and `gamma2_mask` carries `n = A` in its denominator; `A` differs by a factor 7 between the orbits, so both dimensions are dominated by the mask size

## 1. Drift across the two years

| dimension | orbit | n | slope / year | 95% CI | MK p | lag-1 | split-half delta | trend | serial flags | ties |
|---|---|---|---|---|---|---|---|---|---|---|
| `A` | ASC R160 | 31 | +0 | [-25.8, +23.7] | 1.000 | +0.038 | +0.19 | no | none | 0.06 at 61.0 |
| `A` | DESC R51 | 44 | -24 | [-190, +97.7] | 0.739 | +0.19 | -0.012 | no | none | 0.05 at 392.0 |
| `D` | ASC R160 | 31 | +0.00016 | [-0.00599, +0.00581] | 0.865 | +0.021 | +0.22 | no | none | 0.03 at 0.018889858 |
| `D` | DESC R51 | 44 | -0.00531 | [-0.0421, +0.0216] | 0.739 | +0.19 | -0.012 | no | none | 0.05 at 0.086821705 |
| `F` | ASC R160 | 31 | -0.951 | [-3.38, +2.34] | 0.527 | +0.0075 | +0.0042 | no | none | 0.16 at 18.0 |
| `F` | DESC R51 | 44 | -8.69 | [-19.6, +3.38] | 0.056 | +0.16 | -0.21 | no | none | 0.09 at 133.0 |
| `S` | ASC R160 | 31 | -1.55 | [-7.01, +3.06] | 0.196 | +0.14 | -0.2 | no | none | 0.03 at 28.096039857 |
| `S` | DESC R51 | 44 | +0.78 | [-4.78, +2.95] | 0.754 | +0.11 | -0.079 | no | none | 0.02 at 6.920030672 |
| `P` | ASC R160 | 31 | -0.0216 | [-0.044, +0.000731] | 0.013 | +0.06 | -0.37 | no | none | 0.03 at 0.271179487 |
| `P` | DESC R51 | 44 | +0.00308 | [-0.00392, +0.0128] | 0.337 | -0.068 | +0.033 | no | none | 0.02 at 0.162346939 |
| `gamma2_mask` | ASC R160 | 31 | -0.00993 | [-0.0301, +0.00534] | 0.118 | +0.037 | -0.36 | no | none | 0.03 at 0.028037644 |
| `gamma2_mask` | DESC R51 | 44 | -0.000176 | [-0.001, +0.000275] | 0.298 | +0.056 | -0.21 | no | none | 0.02 at 0.001871815 |

## 2. Year repeat — the same calendar month in both campaign years

| dimension | orbit | months | median Δ (y2 − y1) | 95% CI | months + / − |
|---|---|---|---|---|---|
| `A` | ASC R160 | 11 | +6 | [-27.5, +13] | 7 / 4 |
| `A` | DESC R51 | 12 | +6 | [-55.5, +162] | 7 / 5 |
| `D` | ASC R160 | 11 | +0.00173 | [-0.00738, +0.0045] | 6 / 5 |
| `D` | DESC R51 | 12 | +0.00133 | [-0.0123, +0.0359] | 7 / 5 |
| `F` | ASC R160 | 11 | +0.5 | [-2.5, +3] | 6 / 5 |
| `F` | DESC R51 | 12 | -5.75 | [-12.8, +19] | 6 / 6 |
| `S` | ASC R160 | 11 | -1.76 | [-3.77, +0.929] | 4 / 7 |
| `S` | DESC R51 | 12 | -0.165 | [-5.58, +8.76] | 6 / 6 |
| `P` | ASC R160 | 11 | -0.0184 | [-0.0317, -0.00238] | 2 / 9 |
| `P` | DESC R51 | 12 | -0.00275 | [-0.00858, +0.00796] | 5 / 7 |
| `gamma2_mask` | ASC R160 | 11 | -0.0193 | [-0.0199, +0.00868] | 5 / 6 |
| `gamma2_mask` | DESC R51 | 12 | -9.11e-05 | [-0.0013, +0.000327] | 6 / 6 |

## 3. Orbit contrast, unpaired (ASC R160 vs DESC R51)

| dimension | n ASC/DESC | median ASC | median DESC | Cliff's δ | 95% CI | Mann-Whitney p |
|---|---|---|---|---|---|---|
| `A` | 31/44 | 76 | 528 | -1 | [-1, -0.99] | <0.001 |
| `D` | 31/44 | 0.021 | 0.117 | -1 | [-1, -0.99] | <0.001 |
| `F` | 31/44 | 19 | 133 | -1 | [-1, -1] | <0.001 |
| `S` | 31/44 | 25.7 | 14.4 | +0.84 | [+0.7, +0.95] | <0.001 |
| `P` | 31/44 | 0.234 | 0.147 | +0.94 | [+0.85, +1] | <0.001 |
| `gamma2_mask` | 31/44 | 0.0164 | 0.000704 | +0.89 | [+0.74, +0.99] | <0.001 |

## What this can and cannot say

* **It is a different instrument from the pipeline figure.** The two must not be mixed: they correlate at −0.24 (`A`) and 0.13 (`gamma2`).
* **The series is halved and selected by track.** 75 of 150 overpasses; the excluded tracks (ASC R58, DESC R153) are named in the site README.
* **No paired orbit contrast.** Three shared dates; the contrast is a two-sample comparison and is reported as such.
* **`D` and `gamma2_mask` are functions of `A`**, which itself differs by a factor 7 between the orbits.

