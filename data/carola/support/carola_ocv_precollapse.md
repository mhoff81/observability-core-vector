# Carolabrücke deck — pre-collapse trend and late level

* table: `data/carola/support/carola_ocv_channels.csv` (sha256 `3137d0b31fef80b5…`)
* event (documented collapse): **2024-09-11T00:00:00**
* window: last **3 months** before the event — baseline `t < 2024-06-11T00:00:00`, late-pre `2024-06-11T00:00:00 <= t < 2024-09-11T00:00:00`
* dimensions: `D`, `F`, `A` — the unit is the one the table carries
* the two channels are analysed **separately and never pooled**

**How to read this.** The numbers to trust are the intervals: the moving-block bootstrap CI of the slope, and the CI of Cliff's delta. The p-values are printed beside them because they were asked for, but Mann-Kendall assumes independent observations and consecutive overpasses of one deck are not independent, so its p is anti-conservative. `lag-1` and `runs p` say how far that assumption is from holding.

**Sign of the effect.** Cliff's delta is positive when the late-pre values are *larger* than the baseline, so a sinking dimension gives a **negative** delta. Nothing here is a new state: this script writes no state label, and the states stay the two the committed table carries.

## ASCENDING

*n* = 124 pre-collapse overpasses, 2022-09-03 … 2024-09-09, median revisit 5.0 d. Last pre-collapse overpass **2 d before the event**.

### Trend inside the pre-collapse series

| dim | u | n | Theil-Sen slope | CI95 (block 5) | rel. %/a | Kendall τ | MK p | Holm p | lag-1 | runs p |
|---|---|---|---|---|---|---|---|---|---|---|
| `D` | m | 124 | -0.0028 | [-0.0089, +0.0010] | -11.8 % | -0.071 | 0.246 | 0.738 | -0.54 | <0.001 |
| `F` | Hz | 124 | -8.8060 | [-25.9602, +3.2038] | -12.6 % | -0.067 | 0.271 | 0.738 | -0.65 | <0.001 |
| `A` | count | 124 | -27.0529 | [-82.6204, +8.7945] | -13.6 % | -0.064 | 0.295 | 0.738 | -0.56 | <0.001 |

### Late level — last 3 months versus the older baseline

| dim | n base / late | baseline median | late median | Δ median [CI95] | Cliff's δ [CI95] | MW p | Holm p | ties |
|---|---|---|---|---|---|---|---|---|
| `D` | 108 / 16 | 0.0241 | 0.0238 | -0.0004 [-0.0091, +0.0122] | -0.050 [-0.250, +0.152] | 0.751 | 1.000 | 0 |
| `F` | 108 / 16 | 70.0000 | 72.0000 | +2.0000 [-34.0000, +41.0000] | -0.039 [-0.231, +0.161] | 0.803 | 1.000 | 38 |
| `A` | 108 / 16 | 199.0000 | 199.5000 | +0.5000 [-84.5000, +108.0000] | -0.038 [-0.231, +0.160] | 0.811 | 1.000 | 15 |

Holm family: 3 tests (level), 3 tests (trend), alpha = 0.05. † = fewer than 8 values in the late window.

## DESCENDING

*n* = 60 pre-collapse overpasses, 2022-09-12 … 2024-08-20, median revisit 12.0 d. Last pre-collapse overpass **22 d before the event**.

### Trend inside the pre-collapse series

| dim | u | n | Theil-Sen slope | CI95 (block 5) | rel. %/a | Kendall τ | MK p | Holm p | lag-1 | runs p |
|---|---|---|---|---|---|---|---|---|---|---|
| `D` | m | 60 | +0.0051 | [-0.0069, +0.0240] | +14.1 % | +0.103 | 0.246 | 0.677 | -0.01 | 0.435 |
| `F` | Hz | 60 | +8.7382 | [-15.2188, +37.6870] | +11.6 % | +0.097 | 0.281 | 0.677 | +0.06 | 0.193 |
| `A` | count | 60 | +44.5364 | [-62.6388, +210.2955] | +14.4 % | +0.108 | 0.226 | 0.677 | -0.01 | 0.435 |

### Late level — last 3 months versus the older baseline

| dim | n base / late | baseline median | late median | Δ median [CI95] | Cliff's δ [CI95] | MW p | Holm p | ties |
|---|---|---|---|---|---|---|---|---|
| `D` † | 54 / 6 | 0.0370 | 0.0271 | -0.0100 [-0.0178, +0.0118] | -0.080 [-0.383, +0.284] | 0.758 | 1.000 | 1 |
| `F` † | 54 / 6 | 76.5000 | 63.0000 | -13.5000 [-25.5000, +16.0000] | -0.093 [-0.346, +0.241] | 0.721 | 1.000 | 9 |
| `A` † | 54 / 6 | 317.0000 | 233.5000 | -83.5000 [-153.0000, +102.0000] | -0.077 [-0.370, +0.284] | 0.767 | 1.000 | 3 |

Holm family: 3 tests (level), 3 tests (trend), alpha = 0.05. † = fewer than 8 values in the late window.

## Sensitivity — does the verdict depend on the boundary?

Delta per dimension and the Holm decision per window; the event date is also shifted by -12 / +0 / +12 days. No intervals here on purpose: the ladder answers *does the window drive the answer*, not how large the effect is.

| channel | window | cut | n late (D/F/A) | δ D | δ F | δ A | Holm-significant |
|---|---|---|---|---|---|---|---|
| ASCE | 2m-12d | 2024-06-30 | 12/12/12 | +0.02 | +0.02 | +0.03 | — |
| ASCE | 2m+0d | 2024-07-11 | 11/11/11 | -0.03 | -0.05 | -0.02 | — |
| ASCE | 2m+12d | 2024-07-23 | 9/9/9 | -0.10 | -0.11 | -0.09 | — |
| ASCE | 3m-12d | 2024-05-30 | 18/18/18 | -0.09 | -0.08 | -0.08 | — |
| ASCE | 3m+0d | 2024-06-11 | 16/16/16 | -0.05 | -0.04 | -0.04 | — |
| ASCE | 3m+12d | 2024-06-23 | 14/14/14 | -0.01 | -0.01 | +0.01 | — |
| ASCE | 4m-12d | 2024-04-30 | 23/23/23 | -0.22 | -0.18 | -0.19 | — |
| ASCE | 4m+0d | 2024-05-11 | 21/21/21 | -0.19 | -0.15 | -0.16 | — |
| ASCE | 4m+12d | 2024-05-23 | 19/19/19 | -0.11 | -0.10 | -0.10 | — |
| ASCE | 6m-12d | 2024-02-29 | 33/33/33 | -0.18 | -0.15 | -0.15 | — |
| ASCE | 6m+0d | 2024-03-11 | 31/31/31 | -0.16 | -0.13 | -0.14 | — |
| ASCE | 6m+12d | 2024-03-23 | 29/29/29 | -0.17 | -0.14 | -0.15 | — |
| DESC | 2m-12d | 2024-06-30 | 5/5/5 | +0.01 | -0.02 | +0.01 | — |
| DESC | 2m+0d | 2024-07-11 | 4/4/4 | +0.05 | -0.01 | +0.05 | — |
| DESC | 2m+12d | 2024-07-23 | 3/3/3 | +0.18 | +0.09 | +0.18 | — |
| DESC | 3m-12d | 2024-05-30 | 7/7/7 | -0.15 | -0.16 | -0.15 | — |
| DESC | 3m+0d | 2024-06-11 | 6/6/6 | -0.08 | -0.09 | -0.08 | — |
| DESC | 3m+12d | 2024-06-23 | 5/5/5 | +0.01 | -0.02 | +0.01 | — |
| DESC | 4m-12d | 2024-04-30 | 10/10/10 | -0.16 | -0.16 | -0.15 | — |
| DESC | 4m+0d | 2024-05-11 | 9/9/9 | -0.05 | -0.06 | -0.05 | — |
| DESC | 4m+12d | 2024-05-23 | 8/8/8 | -0.04 | -0.05 | -0.04 | — |
| DESC | 6m-12d | 2024-02-29 | 15/15/15 | -0.16 | -0.19 | -0.16 | — |
| DESC | 6m+0d | 2024-03-11 | 14/14/14 | -0.08 | -0.11 | -0.08 | — |
| DESC | 6m+12d | 2024-03-23 | 13/13/13 | -0.17 | -0.20 | -0.17 | — |

## What these numbers do not say

* a negative delta is a statement about the last three months before a date, with the `n` printed beside it — not a prediction and not a cause;
* **no state is created here.** A visible move and these tests can coexist with only two states; this script writes no state label at all, and the guard in `--selftest` checks that none reaches the report;
* a result that exists at one window width only is a finding about the window — the ladder is the check, not a search for the significant cell;
* `DESCENDING` has a gap directly before the collapse (see `gap_before_event_days` above), so its late window ends before the event does.

