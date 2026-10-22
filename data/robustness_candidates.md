# Robustness of the named candidates -- block lengths and size-matched placebos

Pre-registered: the model of each candidate is the one its site report already used, the block lengths are 2, 3 and 5, and the placebo subsamples each pseudo-post to the **real** post-group size. Nothing here searches for a specification that works.

**Gate order** (a cell must pass in this order; a failure is final):

* 1 ladder (A-free and A-adjusted both exclude no effect)
* 2 blocks (2, 3 and 5 all exclude no effect)
* 3 placebo (size-matched rank 1 of 8; rank <= 2 counts as supporting)
* 4 Holm (p < 0.05 within the site's family)

*generalisation-grade* = all four gates. *supporting* = ladder and block gates plus a placebo rank of at most 2. Blocks: 2, 3, 5; 1000 bootstrap draws; 20 subsamples per pseudo event.

| site | stratum | target | model | n post | block 2 | block 3 | block 5 | placebo rank | Holm | gate | verdict |
|---|---|---|---|---|---|---|---|---|---|---|---|
| morandi | ALL | `gamma2` | M1 | 9 | +0.06 [0.55, 1.46] | +0.06 [0.61, 1.41] | +0.06 [0.71, 1.38] | 9/8 | 0.756 | ladder, blocks, placebo_rank_1, placebo_rank_le_2, holm | **not rescued (ladder gate)** |
| morandi | ALL | `D` | M1 | 9 | -1.19 [0.17, 0.55] | -1.19 [0.23, 0.52] | -1.19 [0.24, 0.52] | 1/8 | <0.001 | all pass | **generalisation-grade** |
| morandi | ALL | `F` | M1 | 9 | +0.24 [1.13, 1.44] | +0.24 [1.14, 1.38] | +0.24 [1.15, 1.39] | 1/8 | <0.001 | all pass | **not rescued (pre-registered)** |
| haltern | ALL | `gamma2` | M1 | 11 | +0.33 [0.89, 2.86] | +0.33 [1.10, 2.53] | +0.33 [0.98, 2.66] | 2/8 | 0.090 | ladder, blocks, placebo_rank_1, holm | **not rescued (ladder gate)** |
| carola | ASC 36.9 | `gamma2` | M3 | 85 | +0.27 [0.86, 1.81] | +0.27 [0.83, 1.88] | +0.27 [0.82, 1.85] | 5/8 | 0.972 | ladder, blocks, placebo_rank_1, placebo_rank_le_2, holm | **not rescued (ladder gate)** |
| carola | DESC | `P` | M3 | 86 | -0.23 [0.68, 0.95] | -0.23 [0.68, 0.95] | -0.23 [0.68, 0.96] | 1/8 | 0.100 | holm | **candidate, not Holm-firm** |
| cts | ASC | `P` | M3 | 42 | -0.07 [0.88, 1.00] | -0.07 [0.87, 1.00] | -0.07 [0.88, 1.00] | 1/8 | 0.210 | holm | **candidate, not Holm-firm** |

## What each candidate may now be called

* **morandi / ALL / `gamma2`** — not rescued (ladder gate) (correction: gamma2 was missing only through a column mapping error).
* **morandi / ALL / `D`** — generalisation-grade (strongest independent candidate).
* **morandi / ALL / `F`** — not rescued (pre-registered) (ladder already disagreed; not rescuable).
* **haltern / ALL / `gamma2`** — not rescued (ladder gate) (thin, no geometry adjustment).
* **carola / ASC 36.9 / `gamma2`** — not rescued (ladder gate) (robust, but the placebos say not event-specific).
* **carola / DESC / `P`** — candidate, not Holm-firm (direction-consistent, never Holm-firm).
* **cts / ASC / `P`** — candidate, not Holm-firm (boundary case).

A placebo rank of 1 of 8 is an ordering, not a p-value: with eight pseudo events the floor is about 0.11. It says no pseudo event produced a stronger signature, nothing more.

## Freeze

After this table the analysis is frozen. **Not** searched, and not to be searched later: further covariates, change points, alternative windows, alternative block lengths, further sites, or a different specification for a candidate that failed a gate. A failed cell stays in this table and in the site report with its label; it is not removed.

