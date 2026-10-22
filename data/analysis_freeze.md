# Analysis freeze v2 — the decomposition after the robustness gate

**v1 of this file (same commit series) is void.** It computed `gamma2` from the
pipeline's `coherence_gamma2` column, which is *not* the OCV dimension: the orbit
figures read `gamma2_mask` ("echo coherence") for every site, and the two columns
are essentially uncorrelated (r = 0.09 on Carola, 0.12 on CTS, −0.12 on Haltern;
medians 0.136 vs 0.003). Morandi has no `coherence_gamma2` at all, which is why
v1 wrongly declared its `gamma2` untestable. Everything below is computed on
`gamma2_mask`.

Scope: `carola`, `cts`, `haltern`, `morandi`; targets `gamma2, P, D, F, S`. Site
reports: `data/<site>/support/<site>_ocv_decompose.{json,md}`; cross-site table:
`data/decompose_summary.md`; candidate check: `data/robustness_candidates.md`.

### Why the two columns are not the same measurement (the chip)

The support tables measure a chip whose shape is set by the support geometry, not
the committed 80 x 80 export:

| site | chip rows x cols | support cells | ring cells |
|---|---|---|---|
| Carola | 69 x 125 | 585 | 7 140 |
| CTS | 39 x 145 | 225 | 5 010 |
| Haltern | 35 x 229 | 545 | 6 618 |
| Morandi | 51 x 137 | 861 | 5 448 |

(guard 30 m, reference ring 30-200 m). `tools/ocv_table.py` recomputes `A, D, F,
S, P` and `gamma2` on that chip -- `MASK_COLUMNS`, and `gamma2_column =
"gamma2_mask"` for every site -- and copies everything else (weather, state,
`coherence_gamma2`) verbatim from the committed table, keeping the committed mask
values beside the new ones as `*_committed`. Its own docstring states the
consequence: these numbers are "not comparable with the committed table row by
row". The evidence is in the columns: `coherence_gamma2 == gamma2_export` in
432 of 432 Carola rows and 107 of 107 Haltern rows, i.e. `coherence_gamma2` is
the untouched committed value. The r = 0.09 and the fifty-fold median gap are
therefore **structural** (different chip, different peak reference), not an error
in either column.

Consequences: a support-chip model must use the support columns, and any
comparison against the pre-migration analysis has to go through the
`*_committed` columns. `P` cannot cross chips at all -- it is a pixel-frequency
map computed per chip shape.

## The pre-registered rule (unchanged)

1. **ladder** — A-free and A-adjusted must both exclude no effect, same sign;
2. **blocks** — block lengths 2, 3 and 5 must all exclude no effect;
3. **placebo** — size-matched rank 1 of 8 (rank 2 counts as *supporting*);
4. **Holm** — p < 0.05 within the site's family.

*generalisation-grade* = all four, first failure decides, a failure is final.

## The frozen verdicts

| site | stratum | target | verdict |
|---|---|---|---|
| **morandi** | ALL | **`D`** | **generalisation-grade** |
| morandi | ALL | `F` | not rescued (pre-registered) |
| morandi | ALL | `gamma2` | not rescued (ladder gate) |
| haltern | ALL | `gamma2` | not rescued (ladder gate) |
| carola | ASC 36.9 | `gamma2` | not rescued (ladder gate) |
| carola | DESC | `P` | candidate, not Holm-firm |
| cts | ASC | `P` | candidate, not Holm-firm |

## The numbers to quote

| site | target | effect (log) | ratio | 95% CI | Holm | placebo | n post |
|---|---|---|---|---|---|---|---|
| morandi | `D` | −1.187 | 0.31 | [0.22, 0.52] | <0.001 | 1/8 | 9 |
| morandi | `F` | +0.236 | 1.27 | [1.14, 1.38] | <0.001 | 1/8 | 9 |
| morandi | `gamma2` | +0.060 | 1.06 | [0.63, 1.45] | 0.756 | 9/8 | 9 |
| haltern | `gamma2` | +0.334 | 1.40 | [1.09, 2.52] | 0.090 | 2/8 | 11 |
| carola ASC 36.9 | `gamma2` | +0.274 | 1.32 | [0.82, 1.85] | 0.972 | 5/8 | 85 |
| carola DESC | `P` | −0.227 | 0.80 | [0.68, 0.96] | 0.100 | 1/8 | 86 |
| cts | `P` | −0.070 | 0.93 | [0.87, 1.00] | 0.210 | 1/8 | 42 |

## What the column correction changed

* **Carola `gamma2` is no longer a result.** With the mask-derived coherence the
  A-free effect is large (+0.540 in ASC 36.9, CI [0.04, 1.02]; +1.605 in DESC)
  and it **vanishes once `log A` is in the model** (+0.274, CI [−0.20, 0.61]).
  That is the Carola pattern itself — co-movement with the observed area — not an
  exception to it.
* **Morandi `gamma2` exists and is null**: +0.060 [−0.47, 0.37], placebo rank 9/8.
* **Haltern `gamma2` shrinks** from 12.6× (pipeline column) to 1.40× with a ladder
  that disagrees (A-adj [1.09, 2.52], A-free [−0.09, 0.83]).

The generalisation answer therefore **flips**: with the correct `gamma2` the area
story generalises to three of the four sites (Carola all three strata, Haltern,
Morandi `gamma2`/`P`/`S`, CTS), and **Morandi `D` is the single cell that stays
independent of the observed area** — which is the opposite conclusion to v1.

## The two defects fixed before this freeze

1. **Placebo matching.** The pseudo-post groups were 3.5–17× larger than the real
   post group, so their ranks measured power. Both sides are now cut to
   `min(pseudo_post, real_post)`.
2. **Adaptive block length.** `ladder()` passed `block=5` explicitly, so Haltern
   and Morandi ran with block 5 on 11 and 9 post points. The rule
   `min(5, floor(n/3))` now applies: block 3 for both.

## What was not searched, and will not be

No further covariates, no change-point search, no alternative windows, no block
lengths beyond 2/3/5, no further sites, no different specification for a failed
candidate, and no re-specification after seeing a rank. A failed candidate stays
in `data/robustness_candidates.md` with its label.

## Bautzen — a separate design, outside the four-site tally

Bautzen cannot enter the pre/post decomposition, and the record says why: the
states are a dose ladder (`RS 28 / DS1 92 / DS2 22`); the states **are** the
calendar (RS only in April/May, DS2 only in September/October); there is no
temperature column; `gamma2` is the floor-corrected band coherence and **can be
negative** (DS1 median −0.021); and `D = A / bbox_area` exactly while
`A == n_band` in every row.

What Bautzen has that the others do not is a **control track 500 m east on the
same dates**: every `(date, orbit)` exists twice (71 pairs, 0 incomplete, 0
orphans). The analysis is therefore the paired difference bridge − control,
which removes season, acquisition day and processing — everything the missing
temperature control would otherwise have to absorb.

Estimates are the **A-adjusted** paired difference, Holm within the contrast
(`data/bautzen/support/bautzen_ocv_dose.md`):

| contrast | Δlog A (bridge−control) | ratio | `D` | `F` | `S` | `P` | `gamma2` |
|---|---|---|---|---|---|---|---|
| RS → DS1 | +0.409 | 1.51 | +0.183 | −0.273 | −0.667 | −0.255 | −0.079 |
| **RS → DS2** | **+0.675** | **1.96** | **+0.375** (Holm 0.008) | **−0.586** (0.005) | −0.748 (0.126) | −0.179 (0.205) | +0.203 (0.126) |
| RS → DS1/DS2 | +0.466 | 1.59 | +0.230 | −0.340 | −0.706 (0.220) | −0.247 | −0.013 |

The reading:

* the **area** signature is unambiguous and structural: the bridge's mask grows
  51 % (DS1) and 96 % (DS2) relative to the control, whose own area is flat to
  falling — the within-control contrast for `A` is negative;
* `D` and `F` survive Holm, but `D` **is** `A / bbox_area` and `F` tracks `A` at
  r ≈ 0.98, so both are the area in disguise. Their only non-tautological content
  is that the bounding box grew less than the cell count (Δlog bbox ≈ +0.30
  against Δlog A +0.675): the mask became **denser**;
* of the dimensions that are not area functions, `S` (−0.748) and `gamma2`
  (+0.203) have intervals excluding zero but **Holm 0.126 each**, and `P` shows
  nothing. Under the pre-registered gates none of them is a finding;
* **no cell can be generalisation-grade here**: there are no pseudo events,
  because the state *is* a campaign phase. Bautzen is descriptive/supporting by
  construction and stays outside the four-site tally.

Command and checks: `python3 code/orbits/ocv_dose.py` (~30 s), `--selftest`
(5 checks, including a planted bridge-only step that the pairing finds and an
identical step in both tracks that the pairing removes) and `--verify`
(30 fields recomputed from the table).
