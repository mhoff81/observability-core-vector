# observability-core-vector

Four standalone, offline figure/report packages for four Sentinel-1 observability
sites, side by side in one repository — KDLO reports the same kind of contrast
for a *rebuild epoch* rather than for damage labels, and Carola for the two sides
of a **bridge collapse** (the same structure before and after 2024-09-11). Each
answers the same question about its own **observability core vector**

```
OCV = [gamma2, P, D] = [interferometric coherence, mask persistence, mask density]
```

for its own structural states — and none claims an answer: all four **recompute**
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
| Carolabrücke (Dresden) | `code/carola/`, `data/carola/`, `figures/carola/` | pre-collapse (healthy), post-collapse | **two layers**: the echo mask of every 80 x 80 chip and the per-girder deck mask, both recomputed from a committed window-payload cache (`0.30 x peak`, `5 x np.median`) | [`figures/carola/fig_carola_ocv_paper.md`](figures/carola/fig_carola_ocv_paper.md) |

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
1,713 rows).


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

## Folder structure

```
observability-core-vector/
├── code/
│   ├── lumo/      LUMO analysis + plotting scripts — see code/lumo/README.md
│   ├── bautzen/   the Bautzen mirror (deck-edge mask) — see code/bautzen/README.md
│   ├── kdlo/      the KDLO record package (epoch contrast) — see code/kdlo/README.md
│   ├── carola/    the Carola chip export (pre/post collapse, two mask layers)
│   │              — see code/carola/README.md
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
│   └── README.md
├── figures/
│   ├── lumo/      the .sh wrapper, the .md report and the .png figures
│   ├── bautzen/   the same for the bridge site
│   ├── kdlo/      the same for the KDLO tower
│   ├── carola/    the same for the collapsed Carola bridge
│   └── README.md
├── .gitignore
└── README.md
```

Per site: `code/<site>/` reads `data/<site>/…`, writes `data/<site>/*.json` and
`figures/<site>/{*.md,*.png}`. See the four READMEs for the file-by-file lists.

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

# if `python3` of the PATH is a virtualenv without numpy/matplotlib:
PYTHON=/usr/bin/python3 bash figures/lumo/fig_lumo_ocv_paper.sh

# rebuild the channel tables (LUMO needs the local burst cache, Bautzen only
# the committed rect files, KDLO the two committed raw files plus the committed
# strip cache, Carola the committed mask cache of the chip export)
python3 code/lumo/lumo_ocv_channels_csv.py
python3 code/bautzen/bautzen_ocv_channels_csv.py
python3 code/kdlo/kdlo_ocv_channels_csv.py
python3 code/carola/carola_ocv_channels_csv.py
```

Everything is offline and deterministic: fixed seeds (RNG seed 11 and 2000
bootstrap draws for LUMO, 7 / 42 with 10000 bootstrap draws for Bautzen, 7 with
10000 bootstrap draws and 1000 permutations for KDLO and Carola), no timestamps in
any output, so repeated runs produce identical JSON, Markdown and PNG/JPEG files.
Requirements: Python 3 with numpy, matplotlib and scipy (scipy only for the
p-values of the control tests). Only the KDLO generator can talk to a database,
and only when it is explicitly asked to re-extract nine columns
(`--fetch-db-extra`); the figure scripts never do, and neither the Carola package
nor its figures ever touch the 182 MB payload file (only the one-time
`carola_ocv_masks.py --extract` does).
