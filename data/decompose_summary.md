# Post-collapse decomposition across the sites

The frozen per-candidate verdicts, after the block-length and size-matched-placebo gate, are in `data/robustness_candidates.md` and `data/analysis_freeze.md`.

One row per site x stratum x target. `A-free` is the richest model the site could fit without the area, `A-adjusted` the same model with `log A`. Both come from the site's tier, and the tier is stated: a `thin` label marks a post group below twenty overpasses.

| site | tier | stratum | target | n post | raw | A-free | A-adjusted | 95% CI | Holm | category |
|---|---|---|---|---|---|---|---|---|---|---|
| carola | full | ASC 36.9 | `gamma2` | 85 | 1.71 | 1.72 | 1.32 | [0.82, 1.85] | 0.972 | co-moves with A |
| carola | full | ASC 36.9 | `P` | 85 | 1.03 | 1.06 | 1.05 | [1.00, 1.14] | 0.350 | co-moves with A |
| carola | full | ASC 36.9 | `D` | 85 | 0.83 | 0.83 | 1.02 | [0.99, 1.05] | 0.888 | no effect |
| carola | full | ASC 36.9 | `F` | 85 | 0.85 | 0.85 | 1.01 | [0.98, 1.04] | 0.972 | no effect |
| carola | full | ASC 36.9 | `S` | 85 | +0.67 | +0.92 | -0.03 | [-1.82, +1.15] | 0.972 | no effect |
| carola | full | ASC 44.8 | `gamma2` | 81 | 0.53 | 0.52 | 0.95 | [0.65, 1.38] | 1.000 | co-moves with A |
| carola | full | ASC 44.8 | `P` | 81 | 0.80 | 0.80 | 0.97 | [0.91, 1.04] | 1.000 | co-moves with A |
| carola | full | ASC 44.8 | `D` | 81 | 1.50 | 1.56 | 1.01 | [0.97, 1.04] | 1.000 | co-moves with A |
| carola | full | ASC 44.8 | `F` | 81 | 1.49 | 1.57 | 0.96 | [0.93, 0.99] | 0.110 | sign flip |
| carola | full | ASC 44.8 | `S` | 81 | -3.75 | -4.00 | -0.72 | [-2.69, +1.28] | 1.000 | co-moves with A |
| carola | full | DESC | `gamma2` | 86 | 4.82 | 4.98 | 1.32 | [0.89, 1.77] | 0.432 | co-moves with A |
| carola | full | DESC | `P` | 86 | 0.81 | 0.64 | 0.80 | [0.68, 0.96] | 0.100 | survives |
| carola | full | DESC | `D` | 86 | 0.32 | 0.44 | 0.98 | [0.88, 1.07] | 0.644 | co-moves with A |
| carola | full | DESC | `F` | 86 | 0.31 | 0.37 | 1.05 | [0.99, 1.10] | 0.432 | co-moves with A |
| carola | full | DESC | `S` | 86 | +9.29 | +9.01 | -1.19 | [-3.07, +0.99] | 0.644 | co-moves with A |
| cts | full | ASC | `gamma2` | 42 | 0.86 | 1.15 | 1.26 | [0.74, 1.90] | 1.000 | no effect |
| cts | full | ASC | `P` | 42 | 0.88 | 0.90 | 0.93 | [0.88, 1.00] | 0.210 | survives |
| cts | full | ASC | `D` | 42 | 1.19 | 1.09 | 1.01 | [0.98, 1.04] | 1.000 | no effect |
| cts | full | ASC | `F` | 42 | 1.22 | 1.07 | 1.01 | [0.98, 1.04] | 1.000 | no effect |
| cts | full | ASC | `S` | 42 | -3.25 | -2.31 | -1.91 | [-5.38, +1.71] | 1.000 | no effect |
| haltern | reduced | ALL | `gamma2` | 11 (thin) | 0.74 | 1.21 | 1.40 | [1.09, 2.52] | 0.090 | unstable |
| haltern | reduced | ALL | `P` | 11 (thin) | 1.02 | 1.00 | 1.00 | [0.97, 1.03] | 0.974 | no effect |
| haltern | reduced | ALL | `D` | 11 (thin) | 1.14 | 1.14 | 0.96 | [0.93, 1.00] | 0.200 | unstable |
| haltern | reduced | ALL | `F` | 11 (thin) | 1.10 | 1.17 | 1.07 | [1.00, 1.10] | 0.210 | co-moves with A |
| haltern | reduced | ALL | `S` | 11 (thin) | -2.75 | +3.80 | +4.22 | [-5.06, +13.45] | 0.592 | no effect |
| morandi | minimal | ALL | `gamma2` | 9 (thin) | 1.16 | 0.96 | 1.06 | [0.63, 1.45] | 0.756 | no effect |
| morandi | minimal | ALL | `P` | 9 (thin) | 1.06 | 1.07 | 1.09 | [0.98, 1.17] | 0.288 | no effect |
| morandi | minimal | ALL | `D` | 9 (thin) | 0.39 | 0.30 | 0.31 | [0.22, 0.52] | <0.001 | survives |
| morandi | minimal | ALL | `F` | 9 (thin) | 1.85 | 1.37 | 1.27 | [1.14, 1.38] | <0.001 | survives |
| morandi | minimal | ALL | `S` | 9 (thin) | +1.75 | +4.76 | +4.59 | [-0.46, +8.74] | 0.288 | co-moves with A |

## The tally: does the Carola pattern generalise?

| target | survives | co-moves with A | sign flip | no effect | n/a |
|---|---|---|---|---|---|
| `gamma2` | — | carola/ASC 36.9 (n=85), carola/ASC 44.8 (n=81), carola/DESC (n=86) | — | cts/ASC (n=42), morandi/ALL (n=9) | — |
| `P` | carola/DESC (n=86), cts/ASC (n=42) | carola/ASC 36.9 (n=85), carola/ASC 44.8 (n=81) | — | haltern/ALL (n=11), morandi/ALL (n=9) | — |
| `D` | morandi/ALL (n=9) | carola/ASC 44.8 (n=81), carola/DESC (n=86) | — | carola/ASC 36.9 (n=85), cts/ASC (n=42) | — |
| `F` | morandi/ALL (n=9) | carola/DESC (n=86), haltern/ALL (n=11) | carola/ASC 44.8 (n=81) | carola/ASC 36.9 (n=85), cts/ASC (n=42) | — |
| `S` | — | carola/ASC 44.8 (n=81), carola/DESC (n=86), morandi/ALL (n=9) | — | carola/ASC 36.9 (n=85), cts/ASC (n=42), haltern/ALL (n=11) | — |

`survives` = the A-adjusted interval excludes no effect and keeps the sign of the raw move; `co-moves with A` = the raw move is largely gone once `log A` is in the model; `no effect` = neither ladder separates from no effect. The categories are about the covariates, not about multiplicity -- the Holm column is the multiplicity answer -- and a cell marked `thin` rests on fewer than twenty post overpasses.

