# `figures/` — one subfolder per site: the wrapper plus the report and the figures

```
figures/
├── lumo/      LUMO lattice tower  — fig_lumo_ocv_paper.sh|.md|[A-E].png
├── bautzen/   Bautzen OpenLabs bridge — fig_bautzen_ocv_paper.sh|.md|[A-E].png
│              and the companion fig_bautzen_ocv_amp_phase.sh|.md|.png
├── kdlo/      KDLO media tower — fig_kdlo_ocv_paper.sh|.md|[A-E].png
├── carola/    Carolabrücke, Dresden — fig_carola_ocv_paper.sh|.md|[A-E].png
└── README.md
```

The plotting scripts themselves live in `../code/<site>/` (they share that
folder's Python conventions and import each other); the wrapper in each site
folder runs them from the repository root, so that `data/<site>/…` and
`figures/<site>/…` resolve. Each script reads only its own `data/<site>/`
(channel table + reference JSONs) and writes its JSON artifact back to
`data/<site>/`, the report and the images here. No network, no database, no
burst cache, no 182 MB payload file.

```bash
bash figures/lumo/fig_lumo_ocv_paper.sh            # full run, ~5 min
bash figures/lumo/fig_lumo_ocv_paper.sh --quick    # ~2 min, permutation tests copied
bash figures/bautzen/fig_bautzen_ocv_paper.sh      # ~6 min, deck-edge mask
bash figures/kdlo/fig_kdlo_ocv_paper.sh            # ~1 min, epoch contrast
bash figures/carola/fig_carola_ocv_paper.sh        # ~21 min, bridge collapse
bash figures/carola/fig_carola_ocv_paper.sh --quick # ~6 min, no permutation null

# if the python3 of PATH lacks numpy/matplotlib, name the interpreter:
PYTHON=/usr/bin/python3 bash figures/lumo/fig_lumo_ocv_paper.sh
PYTHON=/usr/bin/python3 bash figures/carola/fig_carola_ocv_paper.sh
```

Images are git-ignored (`figures/**/*.png|jpg|jpeg`); the wrappers regenerate
them from the committed `data/<site>/`. The report (`*.md`) is committed, next
to the wrapper that writes it.

## `lumo/fig_lumo_ocv_paper.sh` → `../code/lumo/fig_lumo_ocv_paper.py`

Five figures on the question of whether the observability core vector
`OCV = [gamma2, P, D]` discriminates the LUMO structural states (healthy, DAM3,
DAM4, DAM6) over the 178 coherence overpasses, plus the English report
[`lumo/fig_lumo_ocv_paper.md`](lumo/fig_lumo_ocv_paper.md) and
[`../data/lumo/fig_lumo_ocv_paper.json`](../data/lumo/fig_lumo_ocv_paper.json):

| file | panel content |
| --- | --- |
| `lumo/fig_lumo_ocv_paper_A.png` | the six channels one panel each: raw per-overpass dots and the median bar for gamma^2, mask area A, density D, fragmentation F, centroid<->peak shift S and persistence P |
| `lumo/fig_lumo_ocv_paper_B.png` | in-sample effect sizes: SDS = 10&#124;Cliff's delta&#124; per component, pooled and per damage state, plus the ranking of the components |
| `lumo/fig_lumo_ocv_paper_C.png` | 4-class LDA: balanced accuracy over the shrinkage grid (6D vs. gamma^2 vs. single components), the leave-one-out confusion matrix and the label-permutation null against the observed accuracy |
| `lumo/fig_lumo_ocv_paper_D.png` | the six binary state pairs: out-of-fold Fisher-LDA scores for the four feature sets (gamma^2, [P,D,S], [gamma^2,P,D], 6D), with the SDS of the headline model and the permutation p per pair |
| `lumo/fig_lumo_ocv_paper_E.png` | controls: severity trend of gamma^2 (Spearman, Jonckheere-Terpstra), season/orbit strata, and the mask brightness / wind covariates per state |

The figures are 183 mm wide, 600 dpi PNGs with the same style as the Espoo
figures (DejaVu Sans + STIX). They are deterministic: fixed seeds, no timestamps
in any annotation, so two runs produce byte-identical files.

`lumo/fig_lumo_ocv_paper.md` is the same content in prose: question and design,
data and mask definition, the results per figure as tables, the pin verdict
against the committed references, the exact commands and the caveats.

## `bautzen/fig_bautzen_ocv_paper.sh` → `../code/bautzen/fig_bautzen_ocv_paper.py`

The bridge counterpart of the LUMO package, on the **deck-edge mask** of the
Bautzen bridge, one column per series (ASC/DESC bridge and the control rectangle
500 m east), plus the English report
[`bautzen/fig_bautzen_ocv_paper.md`](bautzen/fig_bautzen_ocv_paper.md) and
[`../data/bautzen/fig_bautzen_ocv_paper.json`](../data/bautzen/fig_bautzen_ocv_paper.json):

| file | panel content |
| --- | --- |
| `bautzen/fig_bautzen_ocv_paper_A.png` | the six OCV components one panel each: per-overpass dots (state-coloured) and the series/state median bar, per series |
| `bautzen/fig_bautzen_ocv_paper_B.png` | the deck-edge mask layer: band contrast per series and state against the detect gate, mask visibility (Wilson 95% CI), deck-edge pixel count, and the in-series gamma^2 effect size per state |
| `bautzen/fig_bautzen_ocv_paper_C.png` | the two cut dates against cuts shifted by +/-2..12 weeks (placebo sweeps), for band contrast and for gamma^2 |
| `bautzen/fig_bautzen_ocv_paper_D.png` | bridge minus control, date by date: the paired per-date difference with the RS / DS1 / DS2 bands and the step/placebo reading of that difference |
| `bautzen/fig_bautzen_ocv_paper_E.png` | seasonality and environment controls: monthly medians per series, and gamma^2 against air temperature / relative humidity at the overpass instant |

It is deliberately *not* a claim of the same result as LUMO: only the ASC bridge
series passes the deck-edge detection gate, and state and season are perfectly
confounded (RS / DS1 / DS2 are time intervals of one deployment) — which is why
the four series and the shifted-cut placebos are carried through every panel.
Both deviations are documented in `../code/bautzen/README.md` and in the report.

## `bautzen/fig_bautzen_ocv_amp_phase.sh` → `../code/bautzen/fig_bautzen_ocv_amp_phase.py`

The standalone companion of the Bautzen paper figure, on the **amplitude and the
phase of the deck edge over the deck-edge pixels only** (not the whole bridge)
and on their **paired differential phase** between the bridge and its control
rectangle, plus the English report
[`bautzen/fig_bautzen_ocv_amp_phase.md`](bautzen/fig_bautzen_ocv_amp_phase.md)
and [`../data/bautzen/fig_bautzen_ocv_amp_phase.json`](../data/bautzen/fig_bautzen_ocv_amp_phase.json):

| file | panel content |
| --- | --- |
| `bautzen/fig_bautzen_ocv_amp_phase.png` | three rows × one column per series: the coherence amplitude `A(t) = |sum z|/sqrt(N sum |z|^2)` in [0,1], the deck-edge phase `phi(t) = arg(sum z)` in (-pi, pi] and the paired differential phase `dphi(t) = wrap(phi_bridge - phi_control)` in (-pi, pi] (drawn under the bridge column of each geometry), per-overpass dots coloured by state (RS/DS1/DS2) with the state windows as background bands and the two cut dates as dashed lines |

`A(t)^2` is exactly the site's coherence channel, so the figure is verified
internally against the committed `gamma2_band_raw` column date by date (pixel
count, band contrast and the identity `A(t)^2 == gamma2_band_raw`; **515 checks**
counting the per-pair wrap identity `wrap_pi(dphi - (phi_bridge - phi_control))
== 0` and the state-to-state delta checks, exit code != 0 on any failure). The
`dphi` pairs are formed date by date without interpolation — ASC 28 pairs, DESC
35 pairs. The report also gives the **state-to-state** shift
`delta = wrap(circmean(state) - circmean(RS))` (RS → DS1, RS → DS2) for the paired
`dphi` and for each series' absolute `phi`, with a seeded bootstrap interval; all
twelve deltas turn out **not resolvable**. The figure adds one image to the site
and does **not** touch the committed reference JSONs or the figures A–E pin block.

```bash
bash figures/bautzen/fig_bautzen_ocv_amp_phase.sh            # ~2 s
PYTHON=/usr/bin/python3 bash figures/bautzen/fig_bautzen_ocv_amp_phase.sh
```

## `kdlo/fig_kdlo_ocv_paper.sh` → `../code/kdlo/fig_kdlo_ocv_paper.py`

Five figures for the KDLO media tower, on the committed 60-acquisition record and
its two analysis channels (`gamma2`, `A`), plus the English report
[`kdlo/fig_kdlo_ocv_paper.md`](kdlo/fig_kdlo_ocv_paper.md) and
[`../data/kdlo/fig_kdlo_ocv_paper.json`](../data/kdlo/fig_kdlo_ocv_paper.json):

| file | panel content |
| --- | --- |
| `kdlo/fig_kdlo_ocv_paper_A.png` | the ranked channels one panel each: per-acquisition dots (pre vs. rebuild) and the median bar for the 6D echo-mask vector gamma^2, P, D, A, F, S and the controls (the three phase channels, intensity, wind and temperature) |
| `kdlo/fig_kdlo_ocv_paper_B.png` | in-sample effect sizes: SDS = 10&#124;Cliff's delta&#124; per channel with the signed delta, plus the delta and its bootstrap CI for every channel on all rows **and** on the month-matched subset |
| `kdlo/fig_kdlo_ocv_paper_C.png` | 2-class LDA (pre vs. rebuild): balanced accuracy over the shrinkage grid (gamma^2, [P,D,S], the [gamma^2,P,D] headline, the 6D vector and a weather-only competitor), the leave-one-out confusion matrix and the label-permutation verdict |
| `kdlo/fig_kdlo_ocv_paper_D.png` | out-of-fold Fisher-LDA scores per acquisition and the SDS of the four models, all rows vs. month-matched, with the paired bootstrap against gamma^2 |
| `kdlo/fig_kdlo_ocv_paper_E.png` | the alternative explanations: the echo-mode split (compact/intermediate/distributed), the per-month delta table, the mask-size and weather medians, and the singleton acquisition |

Its two states are **epochs** (2022-06-02..2022-12-11 vs. 2023-05-04..2024-09-19,
with 2022-12-23 belonging to neither), so the report and the figures carry the
month-matched control and the weather competitor everywhere, state the
non-significance of the per-acquisition classifier, and separate the four
reference files into two that are recomputed exactly and two that are only pinned
structurally. 2884 checks in a full run, 2875 in `--quick`; exit code != 0 on any
failure.

## `carola/fig_carola_ocv_paper.sh` → `../code/carola/fig_carola_ocv_paper.py`

Five figures on the **collapsed Carolabrücke** (Dresden), on the 1,717 committed
80 x 80 girder chips and the two mask layers of the same chip — the echo mask of
the whole window and the per-girder deck mask — plus the English report
[`carola/fig_carola_ocv_paper.md`](carola/fig_carola_ocv_paper.md) and
[`../data/carola/fig_carola_ocv_paper.json`](../data/carola/fig_carola_ocv_paper.json):

| file | panel content |
| --- | --- |
| `carola/fig_carola_ocv_paper_A.png` | the six echo-mask channels one panel each (gamma^2, P, D, A, F, S): per-chip dots coloured by state (pre vs. post) and the median bar |
| `carola/fig_carola_ocv_paper_B.png` | in-sample effect sizes: the signed Cliff's delta with its bootstrap CI and SDS = 10&#124;delta&#124; for the six channels on all chips, beside the same delta per girder (the deck-mask layer, Girder A/B/C) |
| `carola/fig_carola_ocv_paper_C.png` | 2-class LDA (pre vs. post): balanced accuracy over the shrinkage grid (gamma^2, [P,D], [gamma^2,P,D], the 6D vector and a weather-only competitor), the leave-one-out confusion matrix of the headline model and the label-permutation null |
| `carola/fig_carola_ocv_paper_D.png` | out-of-fold Fisher-LDA scores per chip for the five models, plus the SDS of the OOF score with the paired bootstrap CI against gamma^2 alone |
| `carola/fig_carola_ocv_paper_E.png` | what else moves with the collapse date: the girder / orbit / season strata (SDS within each stratum), the mask-size co-variate (A vs. gamma^2 with the pooled and per-orbit rho) and the weather co-variate (wind vs. gamma^2) |

The figures are 183 mm wide, 600 dpi PNGs in the style of the other three sites
(DejaVu Sans, no matplotlib title box). They are deterministic: fixed seeds, no
timestamps in any annotation, so two runs produce byte-identical files.

The report is the same content in prose: the two mask layers, the state
definition, the results per figure as tables, the pin verdict against the three
committed references (**474 checks**, 414 exact, 60 within `1e-12`, exit code != 0
on any deviation), the reading of the co-variates and the exact commands.

Two caveats are carried through every panel: both states belong to the *same*
structure around the collapse date, so season, weather, traffic and pipeline
generation move with the state (each contrast is shown beside its co-variate
control), and the pipeline's own `coherence_masked_pixels` is not the mask of
these chips (it agrees with the recomputed `A` on 60 of 1,713 rows only).


