# observability-core-vector

Eight standalone, offline figure/report packages for eight Sentinel-1 observability
sites, side by side in one repository — KDLO reports the same kind of contrast
for a *rebuild epoch* rather than for damage labels, and Carola for the two sides
of a **bridge collapse** (the same structure before and after 2024-09-11). Each
answers the same question about its own **observability core vector**

```
OCV = [gamma2, P, D] = [interferometric coherence, mask persistence, mask density]
```

for its own structural states — and none claims an answer: all eight **recompute**
the committed numbers of the site's analysis project from a small input table,
render their figures and report, and **pin every number against the committed
reference JSONs**: a single deviation makes the run exit with a non-zero status.
The claim of the repository is therefore not "we assert", but "we reproduce and
read off".

| site | folders | states | mask layer | results |
| --- | --- | --- | --- | --- |
| LUMO lattice tower (Hannover) | `code/lumo/`, `data/lumo/`, `figures/lumo/` | healthy, DAM3, DAM4, DAM6 | echo mask (`0.30 x peak`, `5 x median`) | [`figures/lumo/fig_lumo_ocv_paper.md`](figures/lumo/fig_lumo_ocv_paper.md) |
| Bautzen OpenLabs bridge | `code/bautzen/`, `data/bautzen/`, `figures/bautzen/` | RS, DS1, DS2 | deck-edge mask (geometrically anchored deck band) | [`figures/bautzen/fig_bautzen_ocv_paper.md`](figures/bautzen/fig_bautzen_ocv_paper.md), [`…fig_bautzen_ocv_amp_phase.md`](figures/bautzen/fig_bautzen_ocv_amp_phase.md) |
| KDLO media tower (Garden City SD) | `code/kdlo/`, `data/kdlo/`, `figures/kdlo/` | pre-collapse, during-rebuild (**epochs**) | echo mask recomputed from a committed 64-strip full-dwell cache (`0.30 x peak`, `5 x upper median`) | [`figures/kdlo/fig_kdlo_ocv_paper.md`](figures/kdlo/fig_kdlo_ocv_paper.md) |
| Carolabrücke (Dresden) | `code/carola/`, `data/carola/`, `figures/carola/` | pre-collapse (healthy), post-collapse | **two layers**: the echo mask of every 80 x 80 chip and the per-girder deck mask, both recomputed from a committed window-payload cache (`0.30 x peak`, `5 x np.median`) | [`figures/carola/fig_carola_ocv_paper.md`](figures/carola/fig_carola_ocv_paper.md), [`…fig_carola_ocv_months.md`](figures/carola/fig_carola_ocv_months.md) |
| Ponte Morandi / Polcevera (Genova) | `code/morandi/`, `data/morandi/`, `figures/morandi/` | pre-collapse (healthy), post-collapse | **two layers**: the echo mask of every 80 x 80 chip and the per-track deck mask, both recomputed from a committed window-payload cache (`0.30 x peak`, `5 x np.median`) | [`figures/morandi/fig_morandi_ocv_paper.md`](figures/morandi/fig_morandi_ocv_paper.md), [`…fig_morandi_ocv_months.md`](figures/morandi/fig_morandi_ocv_months.md) |
| Champlain Towers South (Surfside) | `code/cts/`, `data/cts/`, `figures/cts/` | pre-collapse (healthy), post-collapse | **two layers**: the echo mask of every 80 x 80 chip and the per-track deck mask, both recomputed from a committed window-cache (`0.30 x peak`, `5 x np.median`) — plus a first-class **difference-in-differences** layer (target vs. two control towers) | [`figures/cts/fig_cts_ocv_paper.md`](figures/cts/fig_cts_ocv_paper.md), [`…fig_cts_ocv_months.md`](figures/cts/fig_cts_ocv_months.md); DiD: [`data/cts/reference/cts_reference_did_ctn.md`](data/cts/reference/cts_reference_did_ctn.md), [`…did_placebo_ctn_vs_cte.md`](data/cts/reference/cts_reference_did_placebo_ctn_vs_cte.md) |
| Yeongdeok Wind Farm, Unit 21 (Samgye-ri) | `code/ywf/`, `data/ywf/`, `figures/ywf/` | pre-collapse (healthy), post-collapse | **one layer**: the echo mask of the committed 7 x 7 window-payload (`0.30 x peak`, `5 x np.median`) — the record is not segment-resolved, so a per-mast-section layer would be dishonest — plus a pinned **per-orbit** (ASCENDING / DESCENDING) echo-coverage split, reported as description (both orbits straddle the event), and two references resolved from public catalogues (the CDSE burst **sub-swath** of every chip, the OSM **road anchor**, whose `road_vs_window` block commits why the carriageway lies *outside* the 7 x 7 frame) | [`figures/ywf/fig_ywf_ocv_paper.md`](figures/ywf/fig_ywf_ocv_paper.md) |
| Espoo Kurttila mast (Espoo, Finland) | `code/espoo/`, `data/espoo/`, `figures/espoo/` | **none** — the record is event-free (ASCENDING / DESCENDING are the *orbit geometry*, not a damage state) | **no mask layer**: the record delivers only the mask *size* `A` (`coherence_masked_pixels`) and no raster at all, so the OCV is `[gamma2, A]` and `A` is **copied and proved row by row**, never recomputed | [`figures/espoo/fig_espoo_ocv_paper.md`](figures/espoo/fig_espoo_ocv_paper.md), [`…fig_espoo_ocv_months.md`](figures/espoo/fig_espoo_ocv_months.md) |

The Bautzen package is a **strict mirror** of the LUMO one — same module names,
same file layout, same JSON/Markdown shapes, same determinism — and only the
mask layer (and with it the site constants) differ, so the two packages can be
read as a diff of one another. The KDLO package keeps the same shape (stats
layer, data layer, generator, figure script, pin block) and the same 6D
echo-mask vector `[gamma2, P, D, A, F, S]` with the LUMO mask rule, but its mask
layer does not come from the record: the API stored only the mask *size*, so the
masks are recomputed from a committed 64-strip full-dwell SLC cache (see
`code/kdlo/kdlo_ocv_masks.py`), where `A` and the masked coherence are pinned to
the pipeline's own `coherence_masked_pixels` / `coherence_gamma2`. Its two states
are disjoint time windows, so every effect size is reported on all rows **and**
on a month-matched subset, and the classifier is run against a weather-only
competitor. The Carola package is the same shape again, on the per-girder chip
export of the collapsed bridge: the mask layer is a committed **window-payload
cache** (`code/carola/carola_ocv_masks.py`) that carries the echo mask of every
chip, and the package analyses *two* mask layers over one chip — the echo mask of
the whole 80 x 80 window and the same rule on the girder's own window (its
`A` likewise differs from the pipeline's `coherence_masked_pixels`, on 1,653 of
1,713 rows). The Morandi package is the same shape once more, on the per-track
chips (ASC/DESC) of the collapsed Genova viaduct: the same committed
window-payload cache, with the two tracks playing the role of the girders (the
constant `GIRDERS` and the helper `by_girder` are kept as aliases, so the figure
scripts read as one diff across the whole family); here the site's
`coherence_masked_pixels` is a *fixed 640-px quantile mask* and differs from the
recomputed `A` on all 232 rows, and only 11 chips survive the event, so its
figures are an **honest null** rather than a contrast. The CTS package (Champlain
Towers South, Surfside) is the one **difference-in-differences** site: its
per-track chip export carries four footprints of the *same* orbit (ASC rel-48
IW3) — the collapsed tower (`target`) and two standing control towers
(`control_ctn`, 165 m N; `control_cte`, ~90 m N) plus a beach moisture reference
— which share identical strata, dates and counts, so the package can contrast a
target against a control whose common atmosphere cancels (`code/cts/cts_ocv_did.py`,
pinned against the three committed `data/cts/reference/cts_reference_did_*.json`
to 1e-9, series and markdown). Its mask layer is the same committed window-cache
as Morandi, and the site's `coherence_masked_pixels` is again a *fixed 640-px
quantile mask* (its DiD interaction is pure noise and is pinned structurally
only). The seventh package, YWF (Yeongdeok Wind Farm, Unit 21 — a steel monopole
tower that collapsed onto a public road on 2026-02-02), is the family's **small
record**: 57 committed 7 x 7 windows on 32 dates, 33 unique chips after 24
byte-identical duplicates, all of them the *same* mast section — so the package
carries **one** mask layer and its headline is an **echo-coverage** statement
(16 of 27 pre-event chips carry an echo, 1 of 6 post-event chips does, Fisher
p = 0.085) with an explicit power caveat and **no classifier at all**; the
reference layer is not pinned against an upstream project (there is none) but
*defined* by `fig_ywf_ocv_paper.py --write-reference` and re-derived field by
field on every later run.

The Espoo package (the Kurttila communications mast near Espoo, Finland) is the
family's **event-free control**, and it is the only package with no damage axis at
all: its 150 acquisitions (2024-09-05 … 2026-09-07, 25 months at exactly six each)
carry an empty `damage_label`, `structural_state` and `condition_label` on every
row, and the committed upstream verdict says so itself ("This site has no
ground-truth state label"). There is no event to anchor a window on and no
severity to regress against, so the only contrast the record supports is the one
between its two **orbit geometries** — `ASCENDING` (ascending, afternoon pass,
n = 70) and `DESCENDING` (descending, morning pass, n = 80). The record also
delivers **no mask raster**: `peak_intensity`, `scatterers`, `amplitude`,
`phase_rad`, `mast_peak_row/col` and every other geometry column are empty on all
150 rows (56 of the export's 78), so `DIMS = ["A"]`, the OCV is `[gamma2, A]`, and
`D`, `F`, `S` and `P` are deliberately absent — a generator that produced them
would be inventing an echo-mask layer the record does not have. `A`
(`coherence_masked_pixels`, the delivered mask *size*) is therefore **copied and
proved**, never recomputed: `espoo_ocv_channels_csv.py` checks the copy against the
committed column *and* against the committed reference statistics row by row. The
headline is a geometry effect read as a baseline — median `gamma^2` 0.1746 on ASC
against 0.0258 on DESC, Cliff's delta −0.752, Welch t = +7.31, p = 2.7e-10 — and
the report's whole job is to show that this is **not** a damage signal: nothing
here is a damage signal, which is precisely what makes the site the reference the
other packages' effects have to beat. The phase ladder
(`phase_coherence` / `phase_rms_rad`) exists for only 80 of the 150 acquisitions
(53.3 %, nothing after 2025-10) and `phase_snr_db` for **0**, so it is reported as
*availability* in panel E4 and §A4, never as a third observability channel.


## What the LUMO figures show

| figure | question | content |
| --- | --- | --- |
| A | what do the channels look like? | raw distributions of gamma2, A, D, F, S, P over the 178 coherence overpasses, per state |
| B | which component carries a difference? | in-sample effect sizes: Cliff's delta, bootstrap CI, SDS = 10&#124;delta&#124;, for all six components |
| C | is the state recoverable at all? | 4-class LDA, leave-one-out, shrinkage grid, confusion matrix, label-permutation null |
| D | is the separation state-specific? | out-of-fold scores of the six binary state pairs for four feature sets (+ paired bootstrap vs. gamma2 alone) |
| E | is it the state driving it? | severity trend (Jonckheere-Terpstra, Spearman), season/orbit strata, mask brightness and wind covariates |

`P` belongs to the core vector; `A`, `F` and `S` are controls for the
*morphology* of the same echo mask (they are all derived from the identical mask
as `P` and `D`).

## What the Bautzen figure shows

The same five panels for the bridge: the six channels per series (ASC/DESC
bridge and control 500 m east), effect sizes, state recovery, pairwise
separation and the controls — but on the **deck-edge mask** instead of the echo
mask, and with the site's own two documented caveats:

* only the **ASC bridge** series passes the deck-edge detection gate (a global
  echo rule of the LUMO kind is not usable at this site: the scene clutter pulls
  the reference median up and the mask collapses to ~0 pixels);
* **state and season are perfectly confounded** at this site (RS is spring, DS1
  summer, DS2 autumn), so the environment covariates are a counter-check, not a
  control.

Read [`figures/bautzen/fig_bautzen_ocv_paper.md`](figures/bautzen/fig_bautzen_ocv_paper.md)
for the verdict; it is written by the script and contains no claim that the data
does not support.

## What the KDLO figure shows

The same five panels for the KDLO media tower, on the **60-acquisition committed
record** (16 `pre` before the 2022-12-11 collapse, 43 `rebuild`, one acquisition
in between that belongs to neither) and on the **6D echo-mask vector**
`[gamma2, P, D, A, F, S]` — the OCV proper is `[gamma2, P, D]`, recomputed from
the committed strip cache. The headline is a coherence contrast — median
`gamma^2` 0.003698 → 0.017673, Cliff's delta +0.529 on 688 pairs — and the report
spends its length on why that is *not* automatically damage:

* the second state is an **epoch**, not a label, and it is *during* the rebuild:
  every effect size is therefore also computed on the **month-matched** subset
  (delta +0.429) and the classifier is compared against a **weather-only**
  competitor (which never beats chance);
* the compact-echo count grows from 1 to 12, so the contrast is partly a change of
  **scattering geometry** — the record cannot separate the two readings;
* two of the four reference JSONs (`kdlo_s1_availability.json`,
  `kdlo_fem_frequencies.json`) are catalogue and solver output and are **pinned
  structurally**, not recomputed;
* with the echo-mask dimensions in the model the per-acquisition classifier does
  clear the label permutation (p = 0.023), but the 1D `gamma2` baseline alone
  scores *higher* (0.686 vs. 0.643), so the mask dimensions move the permutation
  null, not the decision boundary — stated in the report rather than hidden.

Read [`figures/kdlo/fig_kdlo_ocv_paper.md`](figures/kdlo/fig_kdlo_ocv_paper.md)
for the full verdict and [`code/kdlo/README.md`](code/kdlo/README.md) for the
module list.

## What the Carola figure shows

The same five panels for the collapsed **Carolabrücke** (Dresden), on the 1,717
committed 80 x 80 girder chips of the per-girder chip export (Girder A 1,189,
B and C 264 each; 896 `pre-collapse (healthy)` before 2024-09-11, 821
`post-collapse`) and on **two mask layers over the same chip** — the echo mask of
the whole window and the same rule on the girder's own window. The headline is a
coherence rise — median `gamma^2` 0.1228 → 0.2000, Cliff's delta +0.172, SDS 1.715
— on top of a mask-persistence collapse (`P`: delta −0.661, SDS 6.606), while the
mask *size* moves much less (`A`: delta −0.100, SDS 0.999). The report is explicit
about what that does **not** mean:

* both states are the **same bridge**, so season, weather, traffic and pipeline
  generation move with the collapse date as well; figures C–E report that
  co-variation rather than a mechanism (the stratified SDS stays between 1.27 and
  3.31 in every girder, orbit and season stratum);
* the nested OLS of the reference is recomputed and printed: the state explains
  `r2 = 0.028` of `gamma2` alone, the mask size `r2 = 0.361`, the inverse-N law
  `r2 = 0.595`, and adding `A` to the state adds `delta r2 = +0.349` — the
  reference's own scenario reading is an **indirect** channel
  (`state -> n_masked -> gamma2`);
* the pipeline's own `coherence_masked_pixels` is *not* the mask of these chips
  (it was computed on the pre-fix asset-point chip): it agrees with the
  recomputed `A` on 60 of 1,713 rows only, so the CSV keeps it as provenance and
  the package derives its own mask from the committed window cache;
* unlike the KDLO record, the multi-dimensional models here *do* separate better
  than the 1D baseline out of fold: SDS 1.638 for `gamma2` alone vs. 6.299 for the
  6D vector (paired delta +4.661, CI [4.079, 5.233]), and the 2-class LDA on the
  6D vector reaches balanced accuracy 0.741 (chance 0.5), which clears the
  label-permutation null (1000 draws, p = 0.001 for accuracy and balanced
  accuracy alike).

Read [`figures/carola/fig_carola_ocv_paper.md`](figures/carola/fig_carola_ocv_paper.md)
for the full verdict and [`code/carola/README.md`](code/carola/README.md) for the
module list.

A companion figure,
[`figures/carola/fig_carola_ocv_months.md`](figures/carola/fig_carola_ocv_months.md),
keeps the **time axis** instead of pooling the two states: it splits the chips
into the four rolling one-month windows around the collapse (M-3, M-2, M-1 before
the event and M+1 after it) and shows the same six echo-mask dimensions per
window (jittered strips for `gamma2`/`P`/`S`/`D`, the `(P, D)` plane and the 6D
fingerprint against the pre-baseline). Its reading is deliberately narrower than
the pooled headline: the three pre windows are **not** stationary, and against
that month-to-month noise only `D` clears the bootstrap CI in the adjacent-month
contrast — `P`, the sharpest pooled difference, moves by only −0.19 with a CI
spanning zero — so the month axis narrows the verdict rather than reproducing it.

## What the Morandi figure shows

The same five panels for the collapsed **Ponte Morandi / Polcevera** (Genova),
on the 232 committed 80 x 80 chips of the per-track chip export (`A_asc` 119,
`A_des` 113; 221 `pre-collapse (healthy)` before 2018-08-14, 11 `post-collapse`)
and on **two mask layers over the same chip** — the echo mask of the whole
window and the same rule on the track's own window. Unlike Carola, the headline
is an **honest null**: of the six echo-mask channels **none** clears the
bootstrap CI (gamma^2 +0.218, CI [-0.164, 0.582]; P -0.130, CI [-0.392, 0.221]);
the only channel that does is `temperature_c` (+0.635, CI [0.454, 0.806]) — a
co-variate, which is the point of panels C–E. The 6D LDA reaches 0.746 raw
accuracy but only 0.521 balanced accuracy (chance 0.5) on the 221 / 11
imbalance, and its label-permutation p is 0.077, so it does not clear the null;
the weather-only competitor reaches 0.685 / 0.662, and the out-of-fold OCV
scores carry SDS 0.06 against 1.64 for the 6D vector, so the states are not
separable in the analysed window. The per-track deck-mask table does separate
within one track (`A_asc`: D -0.58, S +0.71, CI≠0; `A_des`: nothing), but with
only 5–6 post chips per track the report carries that split as one more
co-variate to read, not a verdict.

The report is explicit about why no contrast is a damage verdict:

* only **11** chips survive the collapse (5 `A_asc` + 6 `A_des`) against 221
  before it, so every post-state interval is wide;
* the analysed 80 x 80 window is **not on the collapsed span**, and the archive
  holds no persistent scatterer at this footprint;
* both sides are the **same bridge**, so season, weather and pipeline generation
  move with the collapse date as well;
* the site's `coherence_masked_pixels` is a *fixed 640-px (10 %) quantile mask*,
  not the echo mask this package recomputes — the two agree on 0 of 232 rows
  (panel E), so the export is kept as provenance only.

Read [`figures/morandi/fig_morandi_ocv_paper.md`](figures/morandi/fig_morandi_ocv_paper.md)
for the full verdict and [`code/morandi/README.md`](code/morandi/README.md) for
the module list.

A companion figure,
[`figures/morandi/fig_morandi_ocv_months.md`](figures/morandi/fig_morandi_ocv_months.md),
keeps the **time axis** instead of pooling the two states: the four rolling
one-month windows around the collapse (M-3, M-2, M-1 before the event and M+1
after it), each carrying only 8–11 chips. Against that noise **no** dimension
clears the bootstrap CI in the adjacent-month contrast M+1 vs. M-1, so the month
axis shows the co-variate scatter the pooled figure already flags rather than a
resolvable step.

A companion figure,
[`figures/bautzen/fig_bautzen_ocv_amp_phase.md`](figures/bautzen/fig_bautzen_ocv_amp_phase.md),
reads the **amplitude** `A(t)` and the **phase** `phi(t)` of the deck edge **over
the deck-edge pixels only** (not the whole bridge) across the same four series and
the three states, and adds their **paired differential phase**
`dphi(t) = wrap(phi_bridge - phi_control)`; its `A(t)^2` reproduces the committed
`gamma2_band_raw` channel date by date, which is the only verification it needs.
The same report closes with the **state-to-state** shift
`wrap(circmean(state) - circmean(RS))` (RS → DS1, RS → DS2) for both the paired
`dphi` and the absolute `phi` per series, each with a seeded bootstrap interval —
and states plainly that none of those twelve deltas is resolvable at these sample
sizes.

## What the CTS figure shows

The same six panels for **Champlain Towers South** (Surfside), the one site whose
export is a **difference-in-differences** design: its four footprints are a single
orbit (ASC rel-48 IW3), so the collapsed tower (`target`) shares its orbit, dates
and atmosphere with two standing control towers (`control_ctn`, 165 m N;
`control_cte`, ~90 m N) and a beach moisture reference. Panels A–E mirror the
other sites on the 711 committed chips (543 `pre-collapse (healthy)`, 168
`post-collapse`); panel **F** is the DiD layer no single-structure site can have,
and it re-derives and pins the three committed
`data/cts/reference/cts_reference_did_*.json` (series and markdown) before the run
may pass.

The headline is an **honest null**: the post state is *bare sand* whose intensity
is moisture-driven, so a change in the footprint is not evidence about the
collapse, and **no** channel's interaction clears the 5 % level in either the CTN
or the CTE contrast (smallest p ≈ 0.06, `intensity` vs. CTN — a channel whose
reference *also* steps, p = 0.004, and whose pre-period is already incomparable,
p = 0.000). The pooled effect sizes do separate the states on five of the six
channels (P, D, A, F, S, CI ≠ 0), but that is season / moisture / bare-sand
movement on the *same* footprint; panel F removes exactly those common confounders
and returns a null. The site's `coherence_masked_pixels` is a *fixed 640-px
quantile mask*, constant per chip, so its interaction is pure floating-point noise
and is pinned structurally only.

Read [`figures/cts/fig_cts_ocv_paper.md`](figures/cts/fig_cts_ocv_paper.md) for
the full verdict and [`code/cts/README.md`](code/cts/README.md) for the module
list.

A companion figure,
[`figures/cts/fig_cts_ocv_months.md`](figures/cts/fig_cts_ocv_months.md), keeps
the **time axis** instead of pooling the two states and, uniquely among the sites,
reaches back to the **first images of the footprint**: `F0` is the record start
(`2015-09-21`, the ASC rel-48 IW3 stratum's first acquisition — **no** Sentinel-1
data for this footprint exist in 2014), then the three rolling one-month windows
before the collapse (`M-3`, `M-2`, `M-1`) and the first month after it (`M+1`).
It shows the same six echo-mask dimensions per window (jittered strips for
`gamma2`/`P`/`S`/`D`, the `(P, D)` plane and the 6D fingerprint against the `F0`
baseline). The reading is deliberately narrow: against the **first-images** epoch
`F0` the long baseline clears the bootstrap CI in **0** of the six dimensions, and
against the adjacent `M-1` it is again **0** of six — the `F0` baseline carries
only 8 chips and seven years of drift, and the post month is bare sand, so the
`F0 -> M+1` step conflates long-term drift, season and the collapse and is carried
as co-variation around the collapse cut, **never** as a damage verdict. The
companion adds no number to the A–F pin contract.

## What the YWF figure shows

The YWF site is the family's smallest record: **57** committed 7 x 7 windows on
**32** dates, **33** unique chips after **24** byte-identical duplicates — and all
of them the *same* mast section, because the two requested sections decoded the
same chip at the asset point (`segment_resolved = false`). The package therefore
carries **one** mask layer (the echo mask; a per-mast-section layer would be an
artefact of that duplication), and its six channels exist only where the mask rule
holds: **17** of the 33 chips.

The figure's headline is therefore the **echo coverage**, not a contrast: **16 of
27** pre-event chips carry an echo, **1 of 6** post-event chips does. The Fisher
exact test gives **p = 0.085** — not significant, and it cannot be, because that
single surviving echo sits in the post side; had all six been echo-free, four
chips would already have sufficed at alpha = 0.05, so what limits this site is the
surviving echo, not only the record length. That one post-event echo chip also
lies inside the pre-event cloud of the OCV plane, and **no classifier is fitted
at all** — one chip is not a class, so a leave-one-out accuracy would measure the
sample size, not the tower.

The record's only second cross-cut is the pass geometry, and §A2 of the report
states it: **10 of 20** ascending chips carry an echo against **7 of 13**
descending ones (Fisher exact **p = 1.0000** — the orbits are indistinguishable
here), and that single post-event echo is an *ascending* pass, so the descending
side holds none. `pass_label` (morning / afternoon) is the same cut. The block is
pinned field by field like the state block, but it is read as a description of
where the 17 echoes sit — both orbits straddle the collapse date, so it is not an
orbit contrast, and with 6 post-event chips it cannot repair the power limit above.

Panels A–E: **A** the six echo-mask channels per state over the 17 echo-bearing
chips (a single marker where a state holds one); **B** the coverage bars with the
Fisher test and the coverage timeline across the collapse; **C** the effect sizes
of the mask channels beside the co-variate controls (the season moves with the
state — pre is Aug–Jan, post Feb–Apr — and so do weather and acquisition
generation); **D** the OCV plane; **E** the reasons not to read the contrast as
damage: the pipeline's own `coherence_masked_pixels` (75, 133 … px, computed on the
different-generation asset-point chip) equals the recomputed echo `A` on **0 of
17** rows, plus the echo-mode mix, the season composition and the 24 duplicated
windows.

YWF is also the one site **without an upstream analysis project** (the record *is*
the site's own export), so its reference layer `data/ywf/reference/ywf_ocv_findings.json`
is not pinned against a foreign JSON but **defined** by
`code/ywf/fig_ywf_ocv_paper.py --write-reference` once and then re-derived field by
field (1179 checks) on every later run; a `nan` — a control test that cannot be
computed because the two states are nearly identical — is stored as JSON `null`.

Two facts the export does not carry are resolved from public catalogues and
committed beside that reference: the **pass geometry** of all 33 chips
(`data/ywf/reference/ywf_bursts.json`, resolved against the public CDSE burst
catalogue — ascending is IW3 on relative orbit 54 at 09:2x UT, descending IW2 on
relative orbit 61 at 21:2x UT, **VV and VH interleaved** in both) and the **road**
the tower fell on (`ywf_osm_road.json` + `ywf_geometry.json`, OSM way 577440814),
with `code/ywf/ywf_bursts_resolve.py` and `code/ywf/ywf_osm_anchor.py` as their own
offline `--check` steps and their digests in the figure's pin block. No edge or
band layer is claimed, for two independent reasons the anchor file commits
(`road_vs_window`): the road lies **outside** the 7 x 7 frame — its near edge is
8.10 px (ascending) / 9.33 px (descending) across range from the window centre,
against the 6 px the frame reaches, so containment would need 19-21 px per side —
and ~7 m of carriageway would still be 2 px across range but 0.5 px across
azimuth. The mask rule is itself only ever studied **offline**
(`code/ywf/ywf_mask_sensitivity.py`), because it is shared by the whole family.

Read [`figures/ywf/fig_ywf_ocv_paper.md`](figures/ywf/fig_ywf_ocv_paper.md) for the
full verdict and [`code/ywf/README.md`](code/ywf/README.md) for the module list.

## What the Espoo figure shows

Espoo is the **event-free** site, so its five panels are read as the *baseline*
rather than as a contrast: the same machinery, on a record whose 150 acquisitions
carry no damage label at all, splits almost perfectly along the **orbit geometry**
alone. The figure's job is to show that every number which would look like damage
under a classifier-first reading is accounted for by controls that have nothing to
do with the structure.

Panels A–E: **A** the raw distributions of the two measured channels per orbit
(median `gamma^2` 0.1746 ASC vs. 0.0258 DESC, mask size `A` 6 px vs. 11 px) plus
the phase ladder; **B** the in-sample effect sizes — Cliff's delta with a bootstrap
CI and SDS for `gamma2`, `A`, the two phase columns and the two weather
co-variates: `gamma2` gives delta −0.7521 [−0.8493, −0.6389], SDS 7.52, Welch
t = 7.31 (p = 2.7e-10), the largest effect in the record; **C** the 2-class LDA over
all 150 acquisitions (leave-one-out, shrinkage grid, label-permutation null) —
balanced accuracy 0.7259 for `gamma2` at a chance level of **0.5**, the two-state
counterpart of LUMO's four-state 0.25; **D** the same headline pair as out-of-fold
scores across the four feature sets, with the **out-of-time** split (75 / 75 by
acquisition order, cut at 2025-09-06) showing that `A` and the weather co-variates
add nothing that survives the split; **E** the four no-event controls.

The controls are the point. **E1** re-derives the 14 committed phase-coherence
series field by field — 13 of 14 are stationary (Theil-Sen slope with a
moving-block bootstrap CI, Mann-Kendall, half-split, annual harmonic, lag-1, runs
test, Ljung-Box), and the one flagged series is
`8-18 Hz · full_burst · ASCENDING`. **E2** recomputes the contrast per season,
since there is no event to cut on: every season keeps the same sign and roughly the
same size (`shoulder_Aug-Sep` −0.6356, `summer_Mar-Jul` −0.8653,
`winter_Oct-Feb` −0.6786). **E3** does not regress the weather out — the weather
columns are constant over long stretches — it recomputes the *same* effect size
**inside** each distinct wind value (wind = 5 m/s, 132 acquisitions, delta
−0.7479) and reports the pooled Spearman rho of every channel against wind and
temperature (all |rho| < 0.05); none of them inverts the sign. **E4** reports the
phase ladder as **availability** — `phase_coherence` and `phase_rms_rad` exist for
80 of 150 acquisitions (53.3 %, nothing after 2025-10) and `phase_snr_db` for 0 —
because a channel that is missing on half the record cannot carry a claim about the
record.

The report's reading is explicit: an effect size of t = 7.31 and a balanced
accuracy of 0.7259 would look like damage under a classifier-first reading, the
record's own no-event controls account for all of it, and therefore at every other
site of this repository the orbit-geometry contrast has to be reported *beside* the
damage contrast before the latter can be read — which is exactly what the LUMO,
Carola, Morandi and YWF packages do. Espoo is also the family's only package whose
`A` is not a mask layer at all: it is the delivered mask *size*, copied and proved
row by row (the pipeline's `coherence_masked_pixels`), because the record carries no
raster. Its reference layer is pinned field by field against the three committed
JSONs the upstream project left (`espoo_mast_observability.json`,
`espoo_phase_stationarity.json`, `fig_espoo_channels.json`) plus the generator's own
`output_sha256`, and — like YWF — the package's **own** lock
`data/espoo/reference/espoo_ocv_findings.json` is *defined* by
`code/espoo/fig_espoo_ocv_paper.py --write-reference` once and re-derived on every
later run. No wall-clock value appears in any output, so two runs produce identical
bytes.

Read [`figures/espoo/fig_espoo_ocv_paper.md`](figures/espoo/fig_espoo_ocv_paper.md)
for the full verdict and [`code/espoo/README.md`](code/espoo/README.md) for the
module list.

## Folder structure

```
observability-core-vector/
├── code/
│   ├── lumo/      LUMO analysis + plotting scripts — see code/lumo/README.md
│   ├── bautzen/   the Bautzen mirror (deck-edge mask) — see code/bautzen/README.md
│   ├── kdlo/      the KDLO record package (epoch contrast) — see code/kdlo/README.md
│   ├── carola/    the Carola chip export (pre/post collapse, two mask layers)
│   │              — see code/carola/README.md
│   ├── morandi/   the Morandi per-track chip export (pre/post collapse, two mask
│   │              layers) — see code/morandi/README.md
│   ├── cts/       the CTS chip export (four footprints of one orbit) — the
│   │              difference-in-differences site — see code/cts/README.md
│   ├── ywf/       the Yeongdeok wind-turbine tower record (one 7 x 7 window
│   │              payload, 33 unique chips, one mask layer) and the two resolved
│   │              references (CDSE burst pass geometry, OSM road anchor) — see
│   │              code/ywf/README.md
│   ├── espoo/     the Kurttila mast record — the **event-free control** (150
│   │              acquisitions, no damage label at all; OCV = [gamma2, A], and `A`
│   │              is the delivered mask *size*, copied and proved, never
│   │              recomputed) — see code/espoo/README.md
│   └── README.md
├── data/
│   ├── lumo/      the copied LUMO overpass table, the extended table generated
│   │              from it, and the five committed reference JSONs
│   ├── bautzen/   the four rect files, the generated channel table, the two
│   │              committed reference JSONs and the hourly weather series
│   ├── kdlo/      the two committed raw API dumps, the vendored nine-column DB
│   │              extraction, the generated channel table, the four committed
│   │              reference JSONs and the committed strip cache (64 full-dwell
│   │              SLC strips + manifest)
│   ├── carola/    the committed mask cache of the chip export (1717 echo masks in
│   │              mask coordinates), the committed measurement extract and
│   │              segments/index files, the generated channel table and the
│   │              three committed reference JSONs
│   ├── morandi/   the committed mask cache of the per-track chip export (232 echo
│   │              masks over the 80 x 80 crop + manifest), the committed weather
│   │              table and tracks file, the generated channel table and the
│   │              three committed reference JSONs
│   ├── cts/       the committed mask cache of the four-footprint chip export (711
│   │              echo masks + manifest), the committed measurement extract,
│   │              tracks and census files, the generated channel table and the
│   │              three committed difference-in-differences reference JSONs
│   ├── ywf/       the committed 7 x 7 window payload (57 windows, 38 kB) with its
│   │              index, segments and measurement extract, the committed mask
│   │              cache of its 17 echo windows (+ manifest), the generated channel
│   │              table, the analysis-layer reference the figure script defines and
│   │              the two references resolved from public catalogues (the CDSE
│   │              burst pass geometry, the OSM road extract + derived anchor)
│   ├── espoo/     the copied Kurttila mast export, the generated channel table
│   │              (`A` copied from the delivered `coherence_masked_pixels` and
│   │              proved row by row) and the four committed reference JSONs (the
│   │              asset observability brief, the phase-stationarity table, the
│   │              channel statistics and this package's own findings lock)
│   └── README.md
├── figures/
│   ├── lumo/      the .sh wrapper, the .md report and the .png figures
│   ├── bautzen/   the same for the bridge site
│   ├── kdlo/      the same for the KDLO tower
│   ├── carola/    the same for the collapsed Carola bridge
│   ├── morandi/   the same for the collapsed Morandi viaduct
│   ├── cts/       the same for the collapsed CTS condominium (the DiD site)
│   ├── ywf/       the same for the collapsed Yeongdeok tower (echo coverage)
│   ├── espoo/     the same for the event-free control (orbits, not damage labels)
│   └── README.md
├── .gitignore
└── README.md
```

Per site: `code/<site>/` reads `data/<site>/…`, writes `data/<site>/*.json` and
`figures/<site>/{*.md,*.png}`. See the eight READMEs for the file-by-file lists.

## Quick start

```bash
cd research/observability-core-vector

# LUMO: figures + report + JSON, verified against the five committed references
bash figures/lumo/fig_lumo_ocv_paper.sh            # ~5 min (two permutation tests)
bash figures/lumo/fig_lumo_ocv_paper.sh --quick    # ~2 min (permutation tests copied)

# Bautzen: the same for the bridge (deck-edge mask instead of the echo mask)
bash figures/bautzen/fig_bautzen_ocv_paper.sh            # ~6 min
bash figures/bautzen/fig_bautzen_ocv_paper.sh --quick    # ~1 min (placebos copied)
bash figures/bautzen/fig_bautzen_ocv_amp_phase.sh        # ~2 s (amplitude/phase/differential)

# KDLO: the same for the tower record (epoch contrast, mask recomputed from a
# committed strip cache)
bash figures/kdlo/fig_kdlo_ocv_paper.sh                  # ~1 min
bash figures/kdlo/fig_kdlo_ocv_paper.sh --quick          # ~35 s (no permutation null)

# Carola: the same for the collapsed bridge (per-girder chip export, echo mask and
# deck mask recomputed from a committed window-payload cache)
bash figures/carola/fig_carola_ocv_paper.sh              # ~21 min
bash figures/carola/fig_carola_ocv_paper.sh --quick      # ~6 min (no permutation null)
bash figures/carola/fig_carola_ocv_months.sh             # ~14 s (month-by-month companion)

# Morandi: the same for the collapsed viaduct (per-track chip export, echo mask and
# deck mask recomputed from a committed window-payload cache) — an honest null
bash figures/morandi/fig_morandi_ocv_paper.sh            # full run (permutation null)
bash figures/morandi/fig_morandi_ocv_paper.sh --quick    # ~23 s (no permutation null)
bash figures/morandi/fig_morandi_ocv_months.sh           # ~7 s (month-by-month companion)

# CTS: the difference-in-differences site (four footprints of one orbit, echo mask
# and per-track deck mask recomputed from a committed window cache)
python3 code/cts/cts_ocv_masks.py --check                # mask rule on every cache row
python3 code/cts/cts_ocv_channels_csv.py                 # build the channel table
python3 code/cts/cts_ocv_did.py --verify                 # pin the three DiD references
bash figures/cts/fig_cts_ocv_months.sh                   # ~10 s (F0 + M-3..M+1 companion)

# YWF: the Yeongdeok wind-turbine tower (one 7 x 7 window payload, 33 unique chips,
# one mask layer, echo coverage instead of a contrast — split by state and,
# descriptively, by orbit direction; plus the two references resolved from public
# catalogues: the CDSE pass geometry per chip and the OSM road anchor)
python3 code/ywf/ywf_ocv_masks.py --check                # mask rule on every cache row
python3 code/ywf/ywf_ocv_channels_csv.py                 # build + verify the channel table
python3 code/ywf/ywf_ocv_core.py                         # read-only brief of the committed layer
python3 code/ywf/ywf_bursts_resolve.py                   # --check the pass geometry (default)
python3 code/ywf/ywf_osm_anchor.py                       # --check the road anchor + grid (default)
python3 code/ywf/ywf_mask_sensitivity.py                 # the offline mask-rule sensitivity (writes nothing)
bash figures/ywf/fig_ywf_ocv_paper.sh                    # ~8 s (figures + pins)
bash figures/ywf/fig_ywf_ocv_paper.sh --write-reference  # re-define the reference file

# Espoo: the event-free control (no damage label on any row; the only contrast the
# record supports is the orbit geometry, and `A` is the delivered mask *size*,
# copied and proved row by row — the record carries no mask raster)
python3 code/espoo/espoo_ocv_channels_csv.py             # build + verify the channel table
python3 code/espoo/espoo_ocv_core.py                     # read-only brief of the committed layer
bash figures/espoo/fig_espoo_ocv_paper.sh                # ~4.5 min (figures A-E + pins)
bash figures/espoo/fig_espoo_ocv_paper.sh --quick        # ~1.5 min (no permutation nulls)
bash figures/espoo/fig_espoo_ocv_paper.sh --write-reference  # re-define the package lock
bash figures/espoo/fig_espoo_ocv_months.sh               # ~4 s (month-by-month companion)

# if `python3` of the PATH is a virtualenv without numpy/matplotlib:
PYTHON=/usr/bin/python3 bash figures/lumo/fig_lumo_ocv_paper.sh

# rebuild the channel tables (LUMO needs the local burst cache, Bautzen only
# the committed rect files, KDLO the two committed raw files plus the committed
# strip cache, Carola the committed mask cache of the chip export)
python3 code/lumo/lumo_ocv_channels_csv.py
python3 code/bautzen/bautzen_ocv_channels_csv.py
python3 code/kdlo/kdlo_ocv_channels_csv.py
python3 code/carola/carola_ocv_channels_csv.py
python3 code/morandi/morandi_ocv_channels_csv.py
python3 code/cts/cts_ocv_channels_csv.py
python3 code/ywf/ywf_ocv_channels_csv.py
python3 code/espoo/espoo_ocv_channels_csv.py
```

Everything is offline and deterministic: fixed seeds (RNG seed 11 and 2000
bootstrap draws for LUMO, 7 / 42 with 10000 bootstrap draws for Bautzen, 7 with
10000 bootstrap draws and 1000 permutations for KDLO, Carola and Morandi, 7 with
10000 bootstrap draws for YWF, and 7 with 10000 bootstrap draws — 2000 in the
strata — plus 1000 permutations for Espoo), no timestamps in any output, so
repeated runs produce identical JSON, Markdown and PNG/JPEG files.
Requirements: Python 3 with numpy, matplotlib and scipy (scipy only for the
p-values of the control tests). Only the KDLO generator can talk to a database,
and only when it is explicitly asked to re-extract nine columns
(`--fetch-db-extra`); the figure scripts never do, and neither the Carola package
nor its figures ever touch the 182 MB payload file (only the one-time
`carola_ocv_masks.py --extract` does); the Morandi package likewise reads only
its committed window cache, which the one-time `morandi_ocv_masks.py --extract`
writes from the site's 400 x 400 rect files. The Espoo package has neither a cache
nor a payload at all: its generator and both of its figure scripts read only
`data/espoo/` — the committed export, the generated channel table and the four
committed reference JSONs — so the whole package runs from a fresh clone.

