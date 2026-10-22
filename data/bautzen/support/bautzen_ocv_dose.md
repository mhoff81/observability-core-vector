# Bautzen bridge (RS -> DS1 -> DS2, with a control track) — the damage ladder as a paired contrast

* table `data/bautzen/support/bautzen_ocv_channels.csv` (sha256 `5e7fd543a5b2d7dc…`), 142 rows
* 71 complete `(date, orbit)` pairs, 0 without a partner; 2025-04-02 … 2025-10-29
* pairs per state: DS1 46, DS2 11, RS 14

## Why the difference, and not the state itself

The states are also the calendar: `RS` exists only in April/May, `DS2` only in September/October. A state effect and a season effect are therefore the same axis in this record, and the table carries no temperature. What it does carry is a **control track 500 m east on the same dates**, so every `(date, orbit)` exists twice and the analysis works on the difference bridge minus control. That removes whatever the two tracks share — season, acquisition day, processing — and leaves the structural part.

`gamma2` is the floor-corrected band coherence and can be **negative**, so it is analysed on the identity scale; `P, D, F` are logged. `D` is `A / bbox_area` exactly and `A == n_band` in every row, so neither can be read as evidence against the area.

## The area itself

| contrast | did log A | ratio | n (RS / damaged) |
|---|---|---|---|
| RS_vs_DS1 | +0.409 | 1.51 | 13 / 40 |
| RS_vs_DS2 | +0.675 | 1.96 | 13 / 11 |
| RS_vs_damaged | +0.466 | 1.59 | 13 / 51 |

## The contrasts

### RS_vs_DS1 — RS vs DS1

| target | Did (bridge−control) | + dlog A | within bridge | within control | Holm |
|---|---|---|---|---|---|
| `gamma2` | -0.093 [-0.283, +0.026] | -0.079 [-0.229, +0.022] | -0.037 [-0.11, +0.02] | +0.028 [-0.08, +0.18] | 0.350 |
| `P` | -0.003 [-0.347, +0.171] | -0.255 [-0.685, +0.033] | +0.009 [-0.05, +0.07] | -0.023 [-0.18, +0.25] | 0.350 |
| `D` | +0.194 [-0.120, +0.478] | +0.183 [-0.127, +0.492] | +0.127 [-0.02, +0.29] | -0.065 [-0.31, +0.14] | 0.542 |
| `F` | -0.284 [-0.637, +0.188] | -0.273 [-0.621, +0.204] | -0.201 [-0.48, +0.08] | +0.098 [-0.14, +0.27] | 0.542 |
| `S` | -0.420 [-0.898, +0.239] | -0.667 [-1.256, +0.025] | -0.079 [-0.66, +0.48] | +0.263 [-0.25, +0.66] | 0.350 |

Holm family: 5 targets within this contrast.

### RS_vs_DS2 — RS vs DS2

| target | Did (bridge−control) | + dlog A | within bridge | within control | Holm |
|---|---|---|---|---|---|
| `gamma2` | +0.179 [-0.055, +0.416] | +0.203 [+0.006, +0.428] | +0.036 [-0.05, +0.11] | -0.156 [-0.41, +0.09] | 0.126 |
| `P` | +0.228 [-0.261, +0.479] | -0.179 [-0.757, +0.163] | +0.024 [-0.12, +0.10] | -0.180 [-0.38, +0.17] | 0.205 |
| `D` | +0.495 [+0.217, +0.988] | +0.375 [+0.165, +0.827] | +0.246 [+0.10, +0.40] | -0.207 [-0.66, -0.01] | 0.008 |
| `F` | -0.647 [-1.061, -0.215] | -0.586 [-0.942, -0.217] | -0.251 [-0.57, +0.10] | +0.350 [+0.15, +0.66] | 0.005 |
| `S` | -0.413 [-1.402, +0.312] | -0.748 [-1.696, -0.033] | +0.127 [-0.45, +0.66] | +0.521 [-0.11, +1.45] | 0.126 |

Holm family: 5 targets within this contrast.

### RS_vs_damaged — RS vs DS1/DS2

| target | Did (bridge−control) | + dlog A | within bridge | within control | Holm |
|---|---|---|---|---|---|
| `gamma2` | -0.028 [-0.225, +0.111] | -0.013 [-0.178, +0.112] | -0.023 [-0.10, +0.03] | -0.016 [-0.15, +0.15] | 0.557 |
| `P` | +0.047 [-0.316, +0.195] | -0.247 [-0.703, +0.018] | +0.012 [-0.04, +0.06] | -0.057 [-0.19, +0.23] | 0.244 |
| `D` | +0.259 [-0.060, +0.551] | +0.230 [-0.066, +0.552] | +0.150 [+0.01, +0.31] | -0.096 [-0.35, +0.10] | 0.426 |
| `F` | -0.362 [-0.719, +0.122] | -0.340 [-0.709, +0.135] | -0.210 [-0.48, +0.06] | +0.152 [-0.07, +0.34] | 0.456 |
| `S` | -0.418 [-0.866, +0.231] | -0.706 [-1.288, -0.032] | -0.039 [-0.61, +0.53] | +0.318 [-0.19, +0.69] | 0.220 |

Holm family: 5 targets within this contrast.

## What this can and cannot say

* **No placebo.** The state is a campaign phase, so a pseudo event would be a placebo in time only. These cells are descriptive/supporting by construction and never reach the `generalisation-grade` label of the four-site comparison.
* **No temperature, no season control.** Both are absorbed by the pairing, not estimated.
* **`D` is the area.** `D = A / bbox_area` in every row, so a `D` effect after `log A` adjustment is a change of the bounding box, not an independent dimension.
* **`gamma2` is signed** (floor-corrected) and unlike the other sites not log-transformable.

