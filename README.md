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
| LUMO lattice tower (Hannover) | `code/lumo/`, `data/lumo/`, `figures/lumo/` | healthy, DAM3, DAM4, DAM6 | echo mask (`0.30 x peak`, `5 x median`) | [`figures/lumo/fig_lumo_ocv_paper.md`](figures/lumo/fig_lumo_ocv_paper.md), [`…fig_lumo_ocv_orbits.md`](figures/lumo/fig_lumo_ocv_orbits.md) |
| Bautzen OpenLabs bridge | `code/bautzen/`, `data/bautzen/`, `figures/bautzen/` | RS, DS1, DS2 | deck-edge mask (geometrically anchored deck band) | [`figures/bautzen/fig_bautzen_ocv_paper.md`](figures/bautzen/fig_bautzen_ocv_paper.md), [`…fig_bautzen_ocv_amp_phase.md`](figures/bautzen/fig_bautzen_ocv_amp_phase.md), [`…fig_bautzen_ocv_orbits.md`](figures/bautzen/fig_bautzen_ocv_orbits.md) |
| KDLO media tower (Garden City SD) | `code/kdlo/`, `data/kdlo/`, `figures/kdlo/` | pre-collapse, during-rebuild (**epochs**) | echo mask recomputed from a committed 64-strip full-dwell cache (`0.30 x peak`, `5 x upper median`); the 2015–2021 era has its own 162-strip cache, read only by the series companion | [`figures/kdlo/fig_kdlo_ocv_paper.md`](figures/kdlo/fig_kdlo_ocv_paper.md), [`…fig_kdlo_ocv_months.md`](figures/kdlo/fig_kdlo_ocv_months.md), [`…fig_kdlo_ocv_series.md`](figures/kdlo/fig_kdlo_ocv_series.md) |
| Carolabrücke (Dresden) | `code/carola/`, `data/carola/`, `figures/carola/` | pre-collapse (healthy), post-collapse | **two layers** over the same chip: the echo mask of the 80 x 80 window + the per-girder deck mask, both from a committed window-payload cache (`0.30 x peak`, `5 x np.median`) | [`figures/carola/fig_carola_ocv_paper.md`](figures/carola/fig_carola_ocv_paper.md), [`…fig_carola_ocv_months.md`](figures/carola/fig_carola_ocv_months.md) |
| Ponte Morandi / Polcevera (Genova) | `code/morandi/`, `data/morandi/`, `figures/morandi/` | pre-collapse (healthy), post-collapse | **two layers** as Carola, the two tracks in place of the girders (`0.30 x peak`, `5 x np.median`) | [`figures/morandi/fig_morandi_ocv_paper.md`](figures/morandi/fig_morandi_ocv_paper.md), [`…fig_morandi_ocv_months.md`](figures/morandi/fig_morandi_ocv_months.md) |
| Champlain Towers South (Surfside) | `code/cts/`, `data/cts/`, `figures/cts/` | pre-collapse (healthy), post-collapse | **two layers** as Morandi, plus a **difference-in-differences** layer against two control towers of the same orbit | [`figures/cts/fig_cts_ocv_paper.md`](figures/cts/fig_cts_ocv_paper.md), [`…fig_cts_ocv_months.md`](figures/cts/fig_cts_ocv_months.md); DiD: [`data/cts/reference/cts_reference_did_ctn.md`](data/cts/reference/cts_reference_did_ctn.md), [`…did_placebo_ctn_vs_cte.md`](data/cts/reference/cts_reference_did_placebo_ctn_vs_cte.md) |
| Haltern AV9 wind turbine (Kreis Recklinghausen) | `code/haltern/`, `data/haltern/`, `figures/haltern/` | pre-collapse (healthy), post-collapse | **one layer**: the echo mask of the committed 7 x 7 window payload (`0.30 x peak`, `5 x np.median`, ≥ 2 px) — the record is not segment-resolved and the tower base alone exceeds a 7 x 7 frame, so no second layer is claimed (reasons in the anchor file, `structure_vs_window`) | [`figures/haltern/fig_haltern_ocv_paper.md`](figures/haltern/fig_haltern_ocv_paper.md), [`…fig_haltern_ocv_orbits.md`](figures/haltern/fig_haltern_ocv_orbits.md) |
| Espoo Kurttila mast (Espoo, Finland) | `code/espoo/`, `data/espoo/`, `figures/espoo/` | **none** — event-free (ASCENDING / DESCENDING are orbit geometry, not a damage state) | **no mask layer**: the record delivers only the mask *size* `A` (`coherence_masked_pixels`), no raster, so the OCV is `[gamma2, A]` and `A` is **copied and proved**, never recomputed | [`figures/espoo/fig_espoo_ocv_paper.md`](figures/espoo/fig_espoo_ocv_paper.md), [`…fig_espoo_ocv_months.md`](figures/espoo/fig_espoo_ocv_months.md) |

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
competitor. The Carola package repeats the shape on the per-girder chip export:
**two mask layers** over one chip — the echo mask of the whole 80 x 80 window and
the same rule on the girder's own window — both from a committed window-payload
cache (`code/carola/carola_ocv_masks.py`). The Morandi package is the same shape
once more, on the per-track chips (ASC/DESC) of the collapsed Genova viaduct, the
two tracks playing the role of the girders (the constant `GIRDERS` and the helper
`by_girder` are kept as aliases, so the figure scripts read as one diff across the
whole family). The CTS package (Champlain Towers South, Surfside) is the one
**difference-in-differences** site: its per-track chip export carries four
footprints of the *same* orbit (ASC rel-48 IW3) — the collapsed tower (`target`),
two standing control towers (`control_ctn`, 165 m N; `control_cte`, ~90 m N) and a
beach moisture reference — with identical strata, dates and counts, so the common
atmosphere cancels (`code/cts/cts_ocv_did.py`, pinned to 1e-9 against the three
committed `data/cts/reference/cts_reference_did_*.json`).

The seventh package, Haltern AV9 (a Nordex N149/4.5 turbine of the Haltern AV9
wind park, collapsed on 2021-09-29), is the family's **one-layer record**: 168
committed 7 x 7 windows on 120 dates, all of the *same* mast section once the 48
byte-identical duplicates are collapsed, so its headline is an **echo-coverage**
statement with a power caveat and **no classifier at all**, and its anchor is a
**point** — the turbine's own OSM node, which makes the event date externally
verifiable rather than quoted.

The Espoo package (the Kurttila communications mast near Espoo, Finland) is the
family's **event-free control**, the only package with no damage axis at all: its
150 acquisitions (2024-09-05 … 2026-09-07, 25 months at exactly six each) carry an
empty `damage_label`, `structural_state` and `condition_label` on every row, and
the committed upstream verdict says so itself ("This site has no ground-truth
state label"). With no event to anchor a window on and no severity to regress
against, the only contrast the record supports is between its two **orbit
geometries** — `ASCENDING` (afternoon pass, n = 70) and `DESCENDING` (morning
pass, n = 80) — and the report's whole job is to show that the effect that
produces is **not** a damage signal. Read as a baseline rather than as a contrast,
it is the reference the other packages' effects have to beat. The record carries
no mask raster either (56 of its 78 columns are empty on every row), so
`DIMS = ["A"]` and `D`, `F`, `S` and `P` are deliberately absent.

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

Its amplitude/phase companion,
[`figures/bautzen/fig_bautzen_ocv_amp_phase.md`](figures/bautzen/fig_bautzen_ocv_amp_phase.md),
reads `A(t)` and `phi(t)` **over the deck-edge pixels only** across the same four
series and three states, adds the paired differential phase
`dphi(t) = wrap(phi_bridge - phi_control)`, and reproduces the committed
`gamma2_band_raw` date by date from `A(t)^2` — the only verification it needs. It
closes with the state-to-state shift `wrap(circmean(state) - circmean(RS))` (RS ->
DS1, RS -> DS2) for both `dphi` and the absolute `phi` per series, each with a
seeded bootstrap interval, and states plainly that none of those twelve deltas is
resolvable at these sample sizes.

The orbit companion,
[`figures/bautzen/fig_bautzen_ocv_orbits.md`](figures/bautzen/fig_bautzen_ocv_orbits.md),
keeps one panel per **OCV dimension and orbit** — a 6 x 2 grid in the formal order
`gamma2`, `P`, `D`, `F`, `A`, `S` for ASC and DESC — over the **actual local
overpass instant**, with the RS / DS1 / DS2 windows as coloured bands. It is
bridge-only (the control's ASC `gamma2` range would flatten the bridge signal),
the bridge series are complete on all 36 ASC and 35 DESC overpasses, and it pins
142 rows against the committed metadata (41 checks, 0 failures).

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
  structurally**, not recomputed — the first of them also carries the era block
  (`figure_0`) and the geometry break that separates it from the record;
* the **2015-2021 ascending IW2 era** (162 acquisition days, rel. orbit 136,
  ~00:2x UTC) is a *different* geometry and is **not** in this record: with no
  full-dwell strips for it, `P`/`D`/`F`/`S` cannot be recomputed, so it is kept as
  catalogue facts plus a committed, sha-pinned export of its decoded rows
  (`data/kdlo/fig_kdlo_0_raw.json`, evidence — never pooled into the 60 rows and
  read only by the era-to-collapse series companion, which re-verifies its pinned
  sha256), and the 69-day blind gap into the descending track (2021-12-19 →
  2022-02-26) is drawn, never interpolated;
* with the echo-mask dimensions in the model the per-acquisition classifier does
  clear the label permutation (p = 0.023), but the 1D `gamma2` baseline alone
  scores *higher* (0.686 vs. 0.643), so the mask dimensions move the permutation
  null, not the decision boundary — stated in the report rather than hidden.

The month companion,
[`figures/kdlo/fig_kdlo_ocv_months.md`](figures/kdlo/fig_kdlo_ocv_months.md),
keeps the **time axis** instead of pooling the two epochs — the role the month
companions of Carola, Morandi and CTS play by rolling 24 one-month windows before
the collapse. KDLO cannot: it has no dense pre-collapse stretch, being 60
acquisitions in two disjoint epochs (16 `pre` 2022-06-02…2022-12-11, 43 `rebuild`
2023-05-04…2024-09-19, one acquisition 2022-12-23 in neither), so the axis is the
whole record by **calendar month** — 28 buckets, four of them empty — taken from
`acquisition_ts`. The collapse month 2022-12 is **mixed** (it holds both the last
`pre` acquisition and the singleton), so the collapse cut is drawn inside that
bucket, and the four-month acquisition gap keeps its four empty buckets and a
dotted gap line rather than being filled.
The reading is descriptive: the six pure `pre` months carry 16 acquisitions, the
mixed collapse bucket 2, and **1** of six dimensions clears the bootstrap CI
against the last month before the collapse — the record is read as its own
co-variation, never as damage. The companion adds no number to the A–E pin
contract.

## What the KDLO era-to-collapse figure shows

The era companion,
[`figures/kdlo/fig_kdlo_ocv_series.md`](figures/kdlo/fig_kdlo_ocv_series.md),
puts the record back into the seven years that precede it: one dated axis,
`2015-03-20 … 2022-12-11`, read in **four stretches** (all four taken from the
committed availability reference, not typed). The ascending IW2 era
(`2015-03-20 … 2021-12-19`, 162 acquisition days, relative orbit 136, catalogue
hour 00) is the only stretch whose geometry differs from the record's
(`2022-06-02 … 2022-12-11`, descending IW1, relative orbit 12); between
them sit the 69-day blind gap and a stretch that exists in the CDSE catalogue but
has no committed row here. The figure is 3 x 2 and draws **all six** dimensions on
**both** dated stretches: `gamma2` and the echo-mask size `A` (the two backend
columns whose values both sides carry), and `P`, `D`, `F` and `S` recomputed with
the pipeline's own mask rule from a committed full-dwell strip cache — the
channel table's 64 strips for the record, a **second** cache for the era
(`data/kdlo/strips_era/`, 162 strips, one per era day, held to the era export's
`A` and `gamma2` before a value is drawn). Every panel draws the era's dots and
per-year medians beside the record's dots and line, each **divided by its own
stretch median** (log y). Because the two stretches are imaged on different tracks,
swaths and times of day, their levels are not comparable by construction — what each
stretch *does* is — and the report states that no statistic of this figure is
computed across the break. The era export is read **only** here (sha-verified against
the availability reference before a single value is used) and enters no row of the
channel table; the singleton `2022-12-23` is drawn hollow at the right edge, outside
the axis; and 125 self-checks recompute every drawn number from the four committed
inputs, aborting before any write.

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

The month companion,
[`figures/carola/fig_carola_ocv_months.md`](figures/carola/fig_carola_ocv_months.md),
splits the chips into the **24 rolling one-month windows before** the collapse
(M-24, 2022-09, through M-1, 2024-08) plus the month the event opens (M+1) and
shows the same six echo-mask dimensions per month (one dot per chip for
`gamma2`/`P`/`D`/`A` with the month median and IQR, the `(P, D)` plane and the 6D
fingerprint against the record-start `M-24` baseline). Its reading is deliberately
narrower than the pooled headline: the 24 pre months carry 892 chips (median 40
per month) and are **not** stationary (`gamma2` ranges 0.013 … 0.297 month to
month), and against
that noise only `D` clears the bootstrap CI in the adjacent-month contrast —
`P`, the sharpest pooled difference, moves by only −0.192 with a CI spanning zero —
so the month axis narrows the verdict rather than reproducing it.

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

The month companion,
[`figures/morandi/fig_morandi_ocv_months.md`](figures/morandi/fig_morandi_ocv_months.md),
covers the **24 rolling one-month windows before** the collapse (M-24, 2016-08,
through M-1, 2018-07) plus the month the event opens (M+1). The 24 pre months carry
221 chips (median 10 per month), the first post month 11, and `gamma2`'s monthly
median swings from 0.032 to 0.790 across the pre months alone. Against that noise **no** dimension clears the
bootstrap CI in the adjacent-month contrast M+1 vs. M-1, so the month axis shows
the co-variate scatter the pooled figure already flags rather than a resolvable
step.

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

The month companion,
[`figures/cts/fig_cts_ocv_months.md`](figures/cts/fig_cts_ocv_months.md), uniquely
among the sites reaches back to the **first images of the footprint**: `F0` is the
record start (`2015-09-21`, the ASC rel-48 IW3 stratum's first acquisition — **no**
Sentinel-1 data for this footprint exist in 2014), then the **24 rolling one-month
windows before the collapse** (`M-24`, 2019-06, through `M-1`, 2021-05) and the
first month after it (`M+1`). It shows the same per-month panels against the `F0`
baseline and draws the **acquisition gap** between `F0` and `M-24` as a labelled
dotted line (`gap 3.7 yr`, measured from the window days) instead of
interpolating it. The reading is deliberately narrow: the 24 pre months carry 236
chips (median 10 per month) and the first post month 8, and against the
**first-images** epoch `F0` the long
baseline clears the bootstrap CI in **0** of the six dimensions — and against the
adjacent `M-1` again **0** of six — the post month is bare sand, so the
`F0 -> M+1` step spans six years of drift and conflates that drift, season and the
collapse and is carried as co-variation around the collapse cut, **never** as a
damage verdict. The companion adds no number to the A–F pin contract.

## What the Haltern figure shows

The Haltern AV9 site is the family's **one-layer record with an externally anchored
event**: **168** committed 7 x 7 windows on **120** dates, **120** unique chips after
**48** byte-identical duplicates across 48 payload groups — and all of them the *same*
mast section, because the two requested sections decoded the same chip at the asset
point (`segment_resolved = false`). The package therefore carries **one** mask layer
(the echo mask; a per-mast-section layer would be an artefact of that duplication), and
its six channels exist only where the mask rule holds: **66** of the 120 chips.

The figure's headline is therefore the **echo coverage**, not a contrast: **62 of 108**
pre-event chips carry an echo, **4 of 12** post-event chips do. The Fisher exact test
gives **p = 0.1343** — not significant, and it cannot be, because four of the twelve
post-event chips still carry an echo; had all twelve been echo-free, four chips would
already have sufficed at alpha = 0.05, so what limits this site is the surviving echo,
not only the length of the post side. And **no classifier is fitted at all**: twelve
post-event chips, four of them echo-bearing, are not a class — a leave-one-out accuracy
would measure the sample size, not the turbine.

The record's second cross-cut is the pass geometry, and §A2 of the report states it
per orbit — both orbits straddle the collapse date, so the block is read as a
description of where the echoes sit, never as an orbit contrast. `pass_label`
(morning / afternoon) is the same cut, and the sub-swath is **not** fixed per
direction: the committed pass geometry alternates between IW1/IW3 on relative orbits
15/88 (afternoon) and 37/139 (morning), with **VV and VH interleaved** in both. The
mask rule is a per-window relative threshold and is therefore invariant to that gain,
while the export's own absolute columns (`intensity`, `brightness_ratio`) are not.

Panels A–E: **A** the six echo-mask channels per state over the 66 echo-bearing
chips; **B** the coverage bars with the Fisher test and the coverage timeline across
the collapse; **C** the effect sizes of the mask channels beside the co-variate
controls; **D** the OCV plane; **E** the reasons not to read the contrast as damage:
the pipeline's own `coherence_masked_pixels` (computed on the asset-point chip of the
whole-scene window) equals the recomputed echo `A` on **6 of 66** rows, plus the
echo-mode mix, the pass composition and the 48 duplicated windows. The 54 echo-free
chips are echo-free in **two** different ways, and neither is an empty payload: some
have exactly one pixel passing both thresholds (the rule requires `A >= 2`, so the
mask is refused), the others have none at all, because the peak does not stand 5x
above the median.

The **event date is externally verifiable**: the turbine's own OSM node gained
`damage:date=2021-09-29` and `destroyed=yes` within a day of the collapse, and the
package commits that node version, its changeset and its timestamp. The analysis-layer
reference `data/haltern/reference/haltern_ocv_findings.json` is **defined** by
`code/haltern/fig_haltern_ocv_paper.py --write-reference` once and then re-derived
field by field on every later run.

Two facts the export does not carry are resolved from public catalogues and committed
beside that reference: the **pass geometry** of the chips
(`data/haltern/reference/haltern_bursts.json`) and the **anchor** — the turbine's own
OSM node `8268709968` (`generator:source=wind`, `model=N149`, `operator=RAG Montan
Immobilien, Stadtwerke Haltern`) with the geometry derived from it
(`haltern_osm_node.json` + `haltern_geometry.json`), with
`code/haltern/haltern_bursts_resolve.py` and `code/haltern/haltern_osm_anchor.py` as
their own offline `--check` steps and their digests in the figure's pin block. That
anchor is a **point**, the strongest of this family, and it is checked against the
record's own asset point (`distance_from_asset_point_m`). No edge or band layer is
claimed, for two independent reasons the anchor file commits
(`structure_vs_window`, `structure_vs_pixel_grid`): the tower base alone is 2.4-2.6 px
across range and 0.6 px across azimuth against the 6 px a 7 x 7 frame reaches, so not
one committed pixel spans it; and the N149 rotor (149 m) *would* span a wider frame,
but a rotor is not a band/edge object — its ground projection is a line whose length
and direction are set by the yaw, which neither the export nor the catalogue carries.
The mask rule is itself only ever studied **offline**
(`code/haltern/haltern_mask_sensitivity.py`), because it is shared by the whole family.

The support-sized companion
[`figures/haltern/fig_haltern_ocv_orbits.md`](figures/haltern/fig_haltern_ocv_orbits.md)
draws the same six dimensions on the 107-overpass support cut (44 pre / 5 post ASC,
52 pre / 6 post DESC).

Two companions put the record back on a time axis without adding a number to the
A–E pin contract.
[`figures/haltern/fig_haltern_ocv_weeks.md`](figures/haltern/fig_haltern_ocv_weeks.md)
opens the **pre-collapse window** itself (108 chips on 108 acquisition days,
2021-01-01 -> 2021-09-17, 62 of them echo-bearing, 29 ISO weeks of which 26 hold an
echo) — there is no post-event chip in that figure, the window is the subject. And
[`figures/haltern/fig_haltern_ocv_series.md`](figures/haltern/fig_haltern_ocv_series.md)
draws the **whole record** — 120 chips on 120 acquisition days, 2021-01-01 ->
2021-10-17 — on one dated axis, each dimension divided by its own pre-window median
(taken over the 62 pre-event echo chips, so 1.0 is that median). It is the site's only
figure that shows *where* the collapse falls: no chip sits on 2021-09-29, the event
lies inside the 14-day step `2021-09-17 -> 2021-10-01`, the 9 steps longer than the
cadence are drawn as bands and never interpolated over, and the post side is drawn as
what it is — four echoes among twelve chips, no line, no median and no trend.

## What the Espoo figure shows

Espoo is the **event-free** site, so its five panels are read as the *baseline*
rather than as a contrast: the same machinery splits the 150 acquisitions almost
perfectly along the **orbit geometry** alone, and the figure has to show that
nothing in it needs a structure to be explained.

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

The report's reading is explicit, and it is a statement about the repository, not
just about this site: an effect size of t = 7.31 and a balanced accuracy of 0.7259
look like damage under a classifier-first reading, the record's own no-event
controls account for all of it, so at every other site the orbit-geometry contrast
has to be reported *beside* the damage contrast before the latter can be read —
which is exactly what the LUMO, Carola, Morandi and Haltern packages do. The
reference layer is pinned field by field against the three committed JSONs the
upstream project left (`espoo_mast_observability.json`,
`espoo_phase_stationarity.json`, `fig_espoo_channels.json`) plus the generator's
own `output_sha256`, and — like Haltern — the package's **own** lock
`data/espoo/reference/espoo_ocv_findings.json` is *defined* by
`code/espoo/fig_espoo_ocv_paper.py --write-reference` once and re-derived on every
later run. No wall-clock value appears in any output, so two runs produce identical
bytes.

## Folder structure

```
observability-core-vector/
├── code/          one package per site, plus the cross-site scripts
│   ├── <site>/    the eight records above, plus werder_havel (the family's
│   │              positive control) — see code/<site>/README.md
│   ├── orbits/    one OCV time-series figure per support-sized table, all sites
│   ├── cases/     the six-case overview page, the trajectory plate, the
│   │              correspondence and the prospective model
│   ├── sds/       the cross-site SDS report and the reconfiguration pairs
│   ├── audit_ocv_units.py   the unit audit over every channel table
│   └── README.md
├── data/
│   ├── <site>/    the committed inputs of that site (payload or mask cache,
│   │              reference JSONs) and the generated channel table
│   ├── ocv_channel_units.json, ocv_metadata_quarantine.json, …   the cross-site
│   │              audit files
│   └── README.md
├── figures/
│   ├── <site>/    that site's .sh wrapper, .md report and plate(s)
│   └── README.md
├── .gitignore
└── README.md
```

Per site: `code/<site>/` reads `data/<site>/…`, writes `data/<site>/*.json` and
`figures/<site>/{*.md,*.png}`. `code/README.md`, `data/README.md`,
`figures/README.md` and the site READMEs carry the file-by-file lists.

## Quick start

```bash
cd research/observability-core-vector

# LUMO
bash figures/lumo/fig_lumo_ocv_paper.sh            # ~5 min (two permutation tests)
bash figures/lumo/fig_lumo_ocv_paper.sh --quick    # ~2 min (permutation tests copied)
bash figures/lumo/fig_lumo_ocv_orbits.sh           # ~3 s (6 x 2 OCV time series + pins)

# Bautzen
bash figures/bautzen/fig_bautzen_ocv_paper.sh            # ~6 min
bash figures/bautzen/fig_bautzen_ocv_paper.sh --quick    # ~1 min (placebos copied)
bash figures/bautzen/fig_bautzen_ocv_amp_phase.sh        # ~2 s (amplitude/phase/differential)
bash figures/bautzen/fig_bautzen_ocv_orbits.sh           # ~3 s (OCV time series, 6 x 2)

# KDLO
bash figures/kdlo/fig_kdlo_ocv_paper.sh                  # ~1 min
bash figures/kdlo/fig_kdlo_ocv_paper.sh --quick          # ~35 s (no permutation null)
bash figures/kdlo/fig_kdlo_ocv_months.sh                 # ~11 s (calendar-month companion)
bash figures/kdlo/fig_kdlo_ocv_series.sh                 # <1 s (2015-03-20..2022-12-11 series)

# Carola
bash figures/carola/fig_carola_ocv_paper.sh              # ~21 min
bash figures/carola/fig_carola_ocv_paper.sh --quick      # ~6 min (no permutation null)
bash figures/carola/fig_carola_ocv_months.sh             # ~36 s (24-month companion)

# Morandi
bash figures/morandi/fig_morandi_ocv_paper.sh            # full run (permutation null)
bash figures/morandi/fig_morandi_ocv_paper.sh --quick    # ~23 s (no permutation null)
bash figures/morandi/fig_morandi_ocv_months.sh           # ~15 s (24-month companion)

# CTS (masks, channel table, DiD pins, month companion)
python3 code/cts/cts_ocv_masks.py --check                # mask rule on every cache row
python3 code/cts/cts_ocv_channels_csv.py                 # build the channel table
python3 code/cts/cts_ocv_did.py --verify                 # pin the three DiD references
bash figures/cts/fig_cts_ocv_months.sh                   # ~19 s (F0 + M-24..M+1 companion)

# Haltern (masks, table, the two resolved references, figures)
python3 code/haltern/haltern_ocv_masks.py --check            # mask rule on every cache row
python3 code/haltern/haltern_ocv_channels_csv.py             # build + verify the channel table
python3 code/haltern/haltern_ocv_core.py                     # read-only brief of the committed layer
python3 code/haltern/haltern_bursts_resolve.py               # --check the pass geometry (default)
python3 code/haltern/haltern_osm_anchor.py                   # --check the OSM anchor + grid (default)
python3 code/haltern/haltern_mask_sensitivity.py             # the offline mask-rule sensitivity (writes nothing)
bash figures/haltern/fig_haltern_ocv_paper.sh                # figures A-E + report + pins
bash figures/haltern/fig_haltern_ocv_series.sh               # whole-record series + pins
bash figures/haltern/fig_haltern_ocv_weeks.sh                # pre-collapse weeks companion + pins
bash figures/haltern/fig_haltern_ocv_orbits.sh               # support-chip trajectories + pins
bash figures/haltern/fig_haltern_ocv_paper.sh --write-reference  # re-define the reference file

# Espoo (event-free control)
python3 code/espoo/espoo_ocv_channels_csv.py             # build + verify the channel table
python3 code/espoo/espoo_ocv_core.py                     # read-only brief of the committed layer
bash figures/espoo/fig_espoo_ocv_paper.sh                # ~4.5 min (figures A-E + pins)
bash figures/espoo/fig_espoo_ocv_paper.sh --quick        # ~1.5 min (no permutation nulls)
bash figures/espoo/fig_espoo_ocv_paper.sh --write-reference  # re-define the package lock
bash figures/espoo/fig_espoo_ocv_months.sh               # ~4 s (month-by-month companion)

# the case pages (support-sized chip, nothing pooled across sites)
bash figures/cases/fig_ocv_trajectory.sh                 # ~17 s (the six trajectory plates)
bash figures/cases/fig_ocv_trajectory.sh --verify        # re-derive and compare
python3 code/cases/ocv_trajectory.py --selftest          # 201 checks, writes nothing
python3 code/cases/fig_cases_overview.py                 # the one-page case set

# if `python3` of the PATH is a virtualenv without numpy/matplotlib:
PYTHON=/usr/bin/python3 bash figures/lumo/fig_lumo_ocv_paper.sh
PYTHON=/usr/bin/python3 bash figures/lumo/fig_lumo_ocv_orbits.sh

# rebuild the channel tables (each site README lists what its table reads)
python3 code/lumo/lumo_ocv_channels_csv.py
python3 code/bautzen/bautzen_ocv_channels_csv.py
python3 code/kdlo/kdlo_ocv_channels_csv.py
python3 code/carola/carola_ocv_channels_csv.py
python3 code/morandi/morandi_ocv_channels_csv.py
python3 code/cts/cts_ocv_channels_csv.py
python3 code/haltern/haltern_ocv_channels_csv.py
python3 code/espoo/espoo_ocv_channels_csv.py
```

Everything is offline and deterministic: fixed seeds (RNG seed 11 and 2000
bootstrap draws for LUMO; 7 / 42 with 10000 bootstrap draws for Bautzen; 7 with
10000 bootstrap draws and 1000 permutations for KDLO, Carola, Morandi and Haltern;
7 with 10000 bootstrap draws — 2000 in the strata — plus 1000 permutations for
Espoo), no timestamps in any output, so repeated runs produce identical JSON,
Markdown and PNG/JPEG files.
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
