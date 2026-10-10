# `figures/` — one subfolder per site: the wrapper plus the report and the figures

```
figures/
├── lumo/      LUMO lattice tower  — fig_lumo_ocv_paper.sh|.md|[A-E].png
├── bautzen/   Bautzen OpenLabs bridge — fig_bautzen_ocv_paper.sh|.md|[A-E].png
│              and the companion fig_bautzen_ocv_amp_phase.sh|.md|.png
├── kdlo/      KDLO media tower — fig_kdlo_ocv_paper.sh|.md|[A-E].png
├── carola/    Carolabrücke, Dresden — fig_carola_ocv_paper.sh|.md|[A-E].png
│              and the companion fig_carola_ocv_months.sh|.md|.png
├── morandi/   Ponte Morandi, Genova — fig_morandi_ocv_paper.sh|.md|[A-E].png
│              and the companion fig_morandi_ocv_months.sh|.md|.png
├── cts/       Champlain Towers South, Surfside — fig_cts_ocv_paper.sh|.md|[A-F].png
│              the difference-in-differences site (panel F); its three contrast
│              results are the committed references in ../data/cts/reference/
│              and the companion fig_cts_ocv_months.sh|.md|.png
├── ywf/       Yeongdeok Wind Farm, Unit 21 — fig_ywf_ocv_paper.sh|.md|[A-E].png
│              the smallest record of the family (one mask layer, echo coverage
│              instead of a contrast) and the one site whose reference JSON is
│              defined by its own figure script (--write-reference)
│              plus the companion fig_ywf_ocv_weeks.md|.png (no wrapper)
├── espoo/     Espoo Kurttila mast — fig_espoo_ocv_paper.sh|.md|[A-E].png
│              the **event-free control**: no damage label on any of the 150
│              acquisitions, no mask raster either, so the only contrast is the
│              orbit geometry (ASC / DESC) and `A` is the delivered mask size,
│              copied and proved; plus the companion fig_espoo_ocv_months.sh|.md|.png
└── README.md
```

The plotting scripts themselves live in `../code/<site>/` (they share that
folder's Python conventions and import each other); the wrapper in each site
folder runs them from the repository root, so that `data/<site>/…` and
`figures/<site>/…` resolve. Each script reads only its own `data/<site>/`
(channel table + reference JSONs) and writes its JSON artifact back to
`data/<site>/`, the report and the images here. No network, no database, no
burst cache, no 182 MB payload file.

The CTS site is the one **difference-in-differences** site, so its figure script
carries a sixth panel: panels A–E mirror the shapes of the other sites and panel
**F** is the DiD layer, which re-derives and pins the three committed
`../data/cts/reference/cts_reference_did_*.json` (also reproducible on their own
with `../code/cts/cts_ocv_did.py --verify`).

The Espoo site is the family's **event-free control**, and it is the only site
whose figures have no damage axis at all: its 150 acquisitions carry an empty
`damage_label`, `structural_state` and `condition_label` on every row, the only
contrast the record supports is the one between its two **orbit geometries**
(ASCENDING n = 70 / DESCENDING n = 80), and the record delivers **no mask raster**
(56 of its 78 columns are empty), so `A` is the delivered
`coherence_masked_pixels` — **copied and proved**, never recomputed — and `D`, `F`,
`S` and `P` are deliberately absent. Its reference layer pins three committed JSONs
of the upstream project field by field *and* defines its own
`espoo_ocv_findings.json` lock with `--write-reference`, like YWF does.

```bash
bash figures/lumo/fig_lumo_ocv_paper.sh            # full run, ~5 min
bash figures/lumo/fig_lumo_ocv_paper.sh --quick    # ~2 min, permutation tests copied
bash figures/bautzen/fig_bautzen_ocv_paper.sh      # ~6 min, deck-edge mask
bash figures/kdlo/fig_kdlo_ocv_paper.sh            # ~1 min, epoch contrast
bash figures/carola/fig_carola_ocv_paper.sh        # ~21 min, bridge collapse
bash figures/carola/fig_carola_ocv_paper.sh --quick # ~6 min, no permutation null
bash figures/carola/fig_carola_ocv_months.sh       # ~14 s, month-by-month companion
bash figures/morandi/fig_morandi_ocv_paper.sh      # full run (permutation null)
bash figures/morandi/fig_morandi_ocv_paper.sh --quick # ~23 s, no permutation null
bash figures/morandi/fig_morandi_ocv_months.sh     # ~7 s, month-by-month companion
bash figures/cts/fig_cts_ocv_paper.sh              # ~15 min, permutation null
bash figures/cts/fig_cts_ocv_paper.sh --quick      # ~50 s, no permutation null
bash figures/cts/fig_cts_ocv_months.sh             # ~10 s, months companion (F0 + M-3..M+1)
bash figures/ywf/fig_ywf_ocv_paper.sh              # ~8 s, echo coverage + pins
python3 code/ywf/fig_ywf_ocv_weeks.py              # ~10 s, pre-window companion (no wrapper)
bash figures/espoo/fig_espoo_ocv_paper.sh          # ~4.5 min, event-free control + pins
bash figures/espoo/fig_espoo_ocv_paper.sh --quick  # ~1.5 min, no permutation nulls
bash figures/espoo/fig_espoo_ocv_months.sh         # ~4 s, month-by-month companion

# if the python3 of PATH lacks numpy/matplotlib, name the interpreter:
PYTHON=/usr/bin/python3 bash figures/lumo/fig_lumo_ocv_paper.sh
PYTHON=/usr/bin/python3 bash figures/carola/fig_carola_ocv_paper.sh
PYTHON=/usr/bin/python3 bash figures/morandi/fig_morandi_ocv_paper.sh
PYTHON=/usr/bin/python3 bash figures/cts/fig_cts_ocv_paper.sh
PYTHON=/usr/bin/python3 bash figures/ywf/fig_ywf_ocv_paper.sh
PYTHON=/usr/bin/python3 bash figures/espoo/fig_espoo_ocv_paper.sh
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

The figures are 183 mm wide, 600 dpi PNGs in the style of the other four sites
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

## `carola/fig_carola_ocv_months.sh` → `../code/carola/fig_carola_ocv_months.py`

The **standalone companion** of the Carola paper figure, on the *time axis*
around the collapse instead of the two pooled states: the same six echo-mask
dimensions `x = [gamma2, P, D, A, F, S]` but split into the four rolling
one-month windows that straddle the event (the collapse date `2024-09-11` opens
M+1, so M-1 is pure pre and M+1 pure post), plus the English report
[`carola/fig_carola_ocv_months.md`](carola/fig_carola_ocv_months.md) and
[`../data/carola/fig_carola_ocv_months.json`](../data/carola/fig_carola_ocv_months.json):

| window | range | state | chips |
| --- | --- | --- | --- |
| M-3 | [2024-06-11, 2024-07-11) | pre-collapse (healthy) | 42 |
| M-2 | [2024-07-11, 2024-08-11) | pre-collapse (healthy) | 52 |
| M-1 | [2024-08-11, 2024-09-11) | pre-collapse (healthy) | 56 |
| M+1 | [2024-09-11, 2024-10-11) | post-collapse | 60 |

| file | panel content |
| --- | --- |
| `carola/fig_carola_ocv_months.png` | 2 x 3: (a)-(d) one jittered strip per window for `gamma2`, `P`, `S`, `D` with the window median, the IQR and the collapse cut as a dashed line; (e) the `(P, D)` plane, one dot per chip coloured by window (symlog y, so the near-zero `D` cluster stays visible); (f) the 6D fingerprint — each window's median divided by the pre-baseline median on a log axis |

The figure adds **one** image to the site and does **not** touch the committed
reference JSONs or the figures A–E pin block. Its reading is deliberately
narrower than the pooled headline: the three pre windows are **not** stationary
(`gamma2`'s median roughly halves twice over M-3..M-1, 0.158 → 0.072 → 0.037),
and against that month-to-month noise only `D` clears the bootstrap CI in the
single adjacent-month contrast M+1 vs. M-1 (median 0.0073 → 0.0058, Cliff's
delta −0.36); `P` — the sharpest pooled difference in the A–E figures — moves by
only −0.19 here, with a CI that spans zero. The month axis therefore narrows the
pooled verdict rather than reproducing it, and the report says so.

```bash
bash figures/carola/fig_carola_ocv_months.sh            # ~14 s
PYTHON=/usr/bin/python3 bash figures/carola/fig_carola_ocv_months.sh
```



## `morandi/fig_morandi_ocv_paper.sh` → `../code/morandi/fig_morandi_ocv_paper.py`

Five figures on the **collapsed Ponte Morandi / Polcevera** (Genova), on the 232
committed 80 x 80 chips of the per-track chip export (`A_asc` 119, `A_des` 113)
and the two mask layers of the same chip — the echo mask of the whole window and
the per-track deck mask — plus the English report
[`morandi/fig_morandi_ocv_paper.md`](morandi/fig_morandi_ocv_paper.md) and
[`../data/morandi/fig_morandi_ocv_paper.json`](../data/morandi/fig_morandi_ocv_paper.json):

| file | panel content |
| --- | --- |
| `morandi/fig_morandi_ocv_paper_A.png` | the six echo-mask channels one panel each (gamma^2, P, D, A, F, S): per-chip dots coloured by state (pre vs. post) and the median bar, with log axes on A / D / S |
| `morandi/fig_morandi_ocv_paper_B.png` | in-sample effect sizes: the signed Cliff's delta with its bootstrap CI for the six channels on all chips, beside the same delta per track (the deck-mask layer, `A_asc` / `A_des`) |
| `morandi/fig_morandi_ocv_paper_C.png` | 2-class LDA (pre vs. post): leave-one-out accuracy for gamma^2, [P,D], [gamma^2,P,D], the 6D vector and a weather-only competitor, with the label-permutation p per model |
| `morandi/fig_morandi_ocv_paper_D.png` | the OCV plane (gamma^2 vs. P, gamma^2 vs. D, P vs. D), one dot per chip coloured by state |
| `morandi/fig_morandi_ocv_paper_E.png` | what else moves with the collapse date: the weather co-variates' effect sizes, the track composition per state, the echo-mode mix, and the export-vs-recomputed mask scatter |

The figures are 183 mm wide, 600 dpi PNGs in the style of the other sites
(DejaVu Sans, no matplotlib title box). They are deterministic: fixed seeds, no
timestamps in any annotation, so two runs produce byte-identical files.

The report is the same content in prose: the two mask layers, the state
definition, the results per figure as tables, the pin verdict against the three
committed references (**47 checks**, exit code != 0 on any deviation) and the
exact commands. Its headline is an **honest null**: of the six echo-mask
channels none clears the bootstrap CI (gamma^2 +0.218, P -0.130), the only
channel that does is `temperature_c` (+0.635), and the 6D LDA reaches 0.746 raw
but only 0.521 balanced accuracy on the 221 / 11 imbalance, with a
label-permutation p of 0.077 that does not clear the null. The per-track
deck-mask table does separate within one track (`A_asc`: D -0.58, S +0.71,
CI≠0; `A_des`: nothing), but with only 5–6 post chips per track the report
carries that split as one more co-variate, not a verdict.

Three caveats are carried through every panel: only 11 chips survive the event;
the analysed window is not on the collapsed span; and the site's
`coherence_masked_pixels` is a fixed 640-px quantile mask that agrees with the
recomputed `A` on 0 of 232 rows.


## `morandi/fig_morandi_ocv_months.sh` → `../code/morandi/fig_morandi_ocv_months.py`

The **standalone companion** of the Morandi paper figure, on the *time axis*
around the collapse instead of the two pooled states: the same six echo-mask
dimensions `x = [gamma2, P, D, A, F, S]` but split into the four rolling
one-month windows that straddle the event (the collapse date `2018-08-14` opens
M+1, so M-1 is pure pre and M+1 pure post), plus the English report
[`morandi/fig_morandi_ocv_months.md`](morandi/fig_morandi_ocv_months.md) and
[`../data/morandi/fig_morandi_ocv_months.json`](../data/morandi/fig_morandi_ocv_months.json):

| window | range | state | chips |
| --- | --- | --- | --- |
| M-3 | [2018-05-14, 2018-06-14) | pre-collapse (healthy) | 10 |
| M-2 | [2018-06-14, 2018-07-14) | pre-collapse (healthy) | 8 |
| M-1 | [2018-07-14, 2018-08-14) | pre-collapse (healthy) | 10 |
| M+1 | [2018-08-14, 2018-09-14) | post-collapse | 11 |

| file | panel content |
| --- | --- |
| `morandi/fig_morandi_ocv_months.png` | 2 x 3: (a)-(d) one jittered strip per window for `gamma2`, `P`, `S`, `D` with the window median, the IQR and the collapse cut as a dashed line; (e) the `(P, D)` plane, one dot per chip coloured by window (symlog y, so the near-zero `D` cluster stays visible); (f) the 6D fingerprint — each window's median divided by the pre-baseline median on a log axis |

The figure adds **one** image to the site and does **not** touch the committed
reference JSONs or the figures A–E pin block. Its reading is deliberately
narrower than the pooled one: with only 8–11 chips per window **no** dimension
clears the bootstrap CI in the adjacent-month contrast M+1 vs. M-1 or against
the pooled pre side, so the month axis shows the co-variate scatter the pooled
figure already flags rather than a resolvable step.

```bash
bash figures/morandi/fig_morandi_ocv_months.sh            # ~7 s
PYTHON=/usr/bin/python3 bash figures/morandi/fig_morandi_ocv_months.sh
```

## `cts/fig_cts_ocv_paper.sh` → `../code/cts/fig_cts_ocv_paper.py`

Six figures on the question of whether the observability core vector
`OCV = [gamma2, P, D]` discriminates the pre- and post-collapse states of the
**same** structure, plus the English report
[`cts/fig_cts_ocv_paper.md`](cts/fig_cts_ocv_paper.md) and
[`../data/cts/fig_cts_ocv_paper.json`](../data/cts/fig_cts_ocv_paper.json):

| file | panel content |
| --- | --- |
| `cts/fig_cts_ocv_paper_A.png` | the six echo-mask channels (`gamma2, P, D, A, F, S`) per state — one jittered box per state with `n` (`A`, `D`, `S` on a log axis) |
| `cts/fig_cts_ocv_paper_B.png` | in-sample effect sizes: (a) the pooled Cliff's delta per channel with its bootstrap 95 % CI and SDS; (b) the same per footprint (the four tracks) |
| `cts/fig_cts_ocv_paper_C.png` | is the state learnable from the chips — the 2-class leave-one-out LDA accuracy for `gamma2`, `P+D`, `gamma2+P+D`, the 6 dims and the weather block, against the 0.5 chance line |
| `cts/fig_cts_ocv_paper_D.png` | the OCV plane (`gamma2` vs `P`, `gamma2` vs `D`, `P` vs `D`), one dot per chip coloured by state |
| `cts/fig_cts_ocv_paper_E.png` | what else moves with the state: the weather co-variates' effect sizes, the track composition per state, the echo-mode mix, and the site's fixed 640-px `A` export against the recomputed echo `A` |
| `cts/fig_cts_ocv_paper_F.png` | the difference-in-differences layer: (a) the interaction t-statistic per channel for both controls and the CTN-vs-CTE placebo (±1.96 dashed); (b) the differencing assumption (reference-step p vs. pre-period comparability p); (c) parallel pre-trends (target rho vs. reference rho on the diagonal) |

CTS is the **difference-in-differences** site: its four footprints are one orbit
(the collapsed tower as `target`, the standing towers as `control_ctn` /
`control_cte`, the beach as `reference_beach`, all sharing dates and counts), so
panel **F** does what no single-structure site can — it differences the target
against controls that see the *same* atmosphere and orbit, and pins the three
committed `data/cts/reference/cts_reference_did_{ctn,cte,placebo_ctn_vs_cte}.json`
series-by-series. The script's `pin_block` fails the whole run (non-zero exit) on
any deviation, and it pins the channel table, the mask-cache manifest and the
per-track counts beside the references. Panels A–E mirror the other sites; only
`coherence_masked_pixels` is special-cased as the site's *fixed 640-px quantile
mask* (pinned structurally, never on its `t`/`p`).

The figures are 183 mm wide, 600 dpi PNGs in the style of the other sites
(DejaVu Sans, no matplotlib title box). They are deterministic: fixed seeds, no
timestamps in any annotation, so two runs produce byte-identical files.

The report is the same content in prose as tables. Its headline is the **honest
null of the DiD layer**: the post state is *bare sand* whose intensity is
moisture-driven, so a change in the footprint is not evidence about the collapse,
and **no** channel's interaction clears the 5 % level in either the CTN or the
CTE contrast (smallest p ≈ 0.06, `intensity` vs. CTN — a channel whose reference
*also* steps, p = 0.004, and whose pre-period is already incomparable, p = 0.000).
The pooled effect sizes do separate the states on five of the six channels (P, D,
A, F, S with CI ≠ 0), but that is season / moisture / bare-sand movement on the
*same* footprint, not a damage verdict — panel F removes exactly those common
confounders and returns a null.

```bash
bash figures/cts/fig_cts_ocv_paper.sh            # ~15 min, permutation null
bash figures/cts/fig_cts_ocv_paper.sh --quick    # ~50 s, no permutation null
PYTHON=/usr/bin/python3 bash figures/cts/fig_cts_ocv_paper.sh
```


## `cts/fig_cts_ocv_months.sh` → `../code/cts/fig_cts_ocv_months.py`

The time-axis **companion** of the CTS paper script (like
[`carola/fig_carola_ocv_months.sh`](carola/fig_carola_ocv_months.sh) and
[`morandi/fig_morandi_ocv_months.sh`](morandi/fig_morandi_ocv_months.sh) are for
their sites): instead of the two pooled states it keeps the time axis and shows
the **first Sentinel-1 images of the footprint** (`F0`, the 2015-09 record start),
the three one-month windows before the collapse and the **first month after** it
(`2021-06-24` opens `M+1`). It produces the English report
[`cts/fig_cts_ocv_months.md`](cts/fig_cts_ocv_months.md) and
[`../data/cts/fig_cts_ocv_months.json`](../data/cts/fig_cts_ocv_months.json):

| file | panel content |
| --- | --- |
| `cts/fig_cts_ocv_months.png` | 2 x 3: (a)-(d) one jittered strip per window for `gamma2`, `P`, `S`, `D` with the window median, the IQR, the collapse cut as a dashed line and the record gap (`F0` -> `M-3`, seven years) as a dotted line; (e) the `(P, D)` plane, one dot per chip coloured by window (symlog y, so the near-zero `D` cluster stays visible); (f) the 6D fingerprint — each window's median divided by the `F0` (first-images) baseline median on a log axis |

The five windows are `F0` `[2015-09-21, 2015-11-01)`, `M-3` `[2021-03-24,
2021-04-24)`, `M-2` `[2021-04-24, 2021-05-24)`, `M-1` `[2021-05-24, 2021-06-24)`
and `M+1` `[2021-06-24, 2021-07-24)`, so **48 chips** (8 / 12 / 8 / 12 / 8). The
`F0` baseline is 2015 rather than 2014 because **no** Sentinel-1 data for this
footprint exist before 2015-09-21 (the ASC rel-48 IW3 stratum's first
acquisition; `../data/cts/cts_acquisition_census.md`) — the report states that
rather than inventing a 2014 window. The month axis is descriptive and carries the
site's standing caveat: the post state is *bare sand*, and the long `F0 -> M+1`
baseline conflates seven years of drift, season and the collapse, so the panels are
read as co-variation, **never** as damage. This companion adds **no** number to
the paper script's A–F pin contract.

```bash
bash figures/cts/fig_cts_ocv_months.sh            # ~10 s
PYTHON=/usr/bin/python3 bash figures/cts/fig_cts_ocv_months.sh
```

## `ywf/fig_ywf_ocv_paper.sh` → `../code/ywf/fig_ywf_ocv_paper.py`

Five figures on the question of what can be *observed* on the collapsed Unit 21
tower of the Yeongdeok Wind Farm (Changpo Wind Power Complex, Samgye-ri) — the
steel monopole that fell onto a public road on **2026-02-02** — in the two sides
of that date (`pre-collapse (healthy)` Aug 2025 - Jan 2026, `post-collapse`
Feb - Apr 2026) of the *same* structure. The site's record is the smallest of the
family — **57** committed 7 x 7 windows on **32** dates, **33** unique chips
after **24** byte-identical duplicates, all of them the *same* mast section — so
the package carries **one** mask layer and its claim is deliberately narrow.

The script reads `../data/ywf/ywf_ocv_channels.csv` (33 rows, one per unique
chip, built by `../code/ywf/ywf_ocv_channels_csv.py`), recomputes every aggregate
and writes `../data/ywf/fig_ywf_ocv_paper.json`, the report
[`ywf/fig_ywf_ocv_paper.md`](ywf/fig_ywf_ocv_paper.md) and the five PNGs:

| figure | content |
| --- | --- |
| `ywf/fig_ywf_ocv_paper_A.png` | the six echo-mask channels (`gamma2, P, D, A, F, S`) per state over the **17 echo-bearing chips** — a box where a state holds ≥ 3 chips, a single marker for the one post-event chip (`A`, `D`, `S` on a log axis) |
| `ywf/fig_ywf_ocv_paper_B.png` | the **headline**: (a) the echo coverage per state (with/without echo, the rate and the Fisher exact p); (b) the coverage timeline across the collapse, marker area = unique chips of the day |
| `ywf/fig_ywf_ocv_paper_C.png` | in-sample effect sizes: (a) Cliff's delta + bootstrap 95 % CI and SDS for the six mask channels (no CI where the post side holds one value); (b) the same for the co-variate controls, which move with season, weather and generation |
| `ywf/fig_ywf_ocv_paper_D.png` | the OCV plane (`gamma2` vs `P`, `gamma2` vs `D`, `P` vs `D`) — 16 pre-event chips and 1 post-event chip, **no classifier**: one chip is not a class |
| `ywf/fig_ywf_ocv_paper_E.png` | what else moves with the collapse date: (a) the pipeline's `coherence_masked_pixels` (75, 133 … px, a different-generation asset-point chip) against the recomputed echo `A` — 0 of 17 rows coincide; (b) the echo-mode mix; (c) the season composition; (d) the 24 duplicated windows that make the record not segment-resolved |

The reading is a **coverage** statement, not a contrast: **16 of 27** pre-event
chips carry an echo and **1 of 6** post-event chips does, so the Fisher exact test
gives **p = 0.085** — not significant, and it cannot be, because the single
surviving echo is in the post side. Had all six post-event chips been echo-free,
four chips would already have sufficed at alpha = 0.05, so what limits the site is
the surviving echo and not only the record length. The report states that instead
of a damage verdict.

This is also the one site **without an upstream analysis project**, so the
reference layer is defined here: `bash figures/ywf/fig_ywf_ocv_paper.sh
--write-reference` writes `../data/ywf/reference/ywf_ocv_findings.json` (and proves
the round trip), and every later run re-derives it field by field — 1179 checks
against the reference, the channel metadata, the committed mask cache, the sha256
digests of the committed inputs and the digests of the two references resolved
from public catalogues (`ywf_bursts.json`, `ywf_osm_road.json`, `ywf_geometry.json`)
— exiting non-zero on the first deviation.

The report's four blockers are backed by the committed anchor file: its
`road_vs_window` block proves the carriageway lies **outside** the 7 x 7 frame
(its near edge is 8.10-9.33 px across range from the window centre, against the
6 px the frame reaches), and `../code/ywf/ywf_mask_sensitivity.py` sweeps the
mask rule **offline** on the committed payload — the rule is family-wide, so no
variant is switched on and no figure depends on either.

```bash
bash figures/ywf/fig_ywf_ocv_paper.sh                    # ~8 s, figures + report + pins
bash figures/ywf/fig_ywf_ocv_paper.sh --write-reference  # <1 s, re-define the reference
PYTHON=/usr/bin/python3 bash figures/ywf/fig_ywf_ocv_paper.sh
```

## `ywf/fig_ywf_ocv_weeks.md` → `../code/ywf/fig_ywf_ocv_weeks.py` (no wrapper)

The companion that opens the pre-collapse window instead of pooling it: the OCV
time series of the **27** pre-event chips on **26** acquisition days (23 ISO
weeks, **16** of them echo-bearing) from 2025-08-06 to the collapse line at
`W-0`, with **no** post-event chip drawn. It adds no number to the A-E pin
contract of the paper figure; it re-derives the pre window from the committed
channel table and mask cache and re-verifies the three reference digests of the
paper figure's `sources` block (422 checks in all).

| panel | content |
| --- | --- |
| a | the echo of every acquisition (`A` of 49) as a stem, grey where the chip has no echo, coloured by echo mode, with the chip `intensity` on a twin log axis |
| b | `gamma2` over the 16 echo-bearing chips plus its gap-aware 4-week rolling median |
| c | `P` (the same) |
| d | `D` (the same) |
| e | the three OCV channels normalised to their pre-window median, one log axis |
| f | what else moves over these weeks: temperature, wind speed and gust |

The time axis is weeks before the collapse (`x = -(2026-02-02 - day)/7`, 0 at the
event); the four **12-day repeat gaps** (12 Aug → 24 Aug, 24 Aug → 5 Sep,
3 Jan → 15 Jan, 15 Jan → 27 Jan, when the ascending and the descending pass
alternate) are drawn as bands and never interpolated over. The figure's point is
the *reading*: the coordinate is not stationary — temperature falls from 33.2 C to
0.6 C across the window — and the coverage moves with it, while the channels
themselves show no monotone drift (Spearman against acquisition time over the 16
echo chips: `gamma2` +0.135, p 0.62; `P` −0.17, p 0.52; `D` −0.25, p 0.35). A
"before vs. after" statement computed from the pooled window is therefore read
against the co-variates of panel f, never as a change at the event. `P` is a
**window-wide reference** (the mean pixel frequency inside the own mask over the
17 echo-bearing chips of the whole record), so these chips are not independent
draws of it; the report says so.

```bash
python3 code/ywf/fig_ywf_ocv_weeks.py              # ~10 s, figure + JSON + report
python3 code/ywf/fig_ywf_ocv_weeks.py --no-figures # cross-checks only
```

## `espoo/fig_espoo_ocv_paper.sh` → `../code/espoo/fig_espoo_ocv_paper.py`

Five figures on the **event-free control**: the Espoo Kurttila mast record has **no
damage label on any of its 150 acquisitions** (2024-09-05 … 2026-09-07, 25 months at
exactly 6 each) and the committed upstream verdict says so itself ("This site has no
ground-truth state label"), so there is no severity axis to regress against. The only
contrast the record supports is the one between its two **orbit geometries**
(ASCENDING / afternoon pass, n = 70; DESCENDING / morning pass, n = 80), and the
figure set is read as the **baseline** every damage claim elsewhere has to beat. The
site's OCV is `[gamma2, A]` where `A` = `coherence_masked_pixels`, the size of the
mask the pipeline found: this export carries **no mask raster** (`peak_intensity`,
`scatterers`, `amplitude`, `phase_rad`, `mast_peak_row/col` and every other geometry
column are empty on all 150 rows — 56 of 78), so `D`, `F`, `S` and `P` are absent and
`A` is copied from the committed column and proved row by row.

The script recomputes and then pins the effect sizes (median `gamma^2` 0.1746 ASC vs.
0.0258 DESC; Cliff's delta −0.7521 with CI [−0.8493, −0.6389], SDS 7.52; Welch
t = 7.31, p = 2.7e-10), the 2-class LDA over the whole record (leave-one-out,
balanced accuracy 0.7259 at a chance level of 0.5, label-permutation null) and the
four no-event controls. It writes `../data/espoo/fig_espoo_ocv_paper.json`, the report
[`espoo/fig_espoo_ocv_paper.md`](espoo/fig_espoo_ocv_paper.md) and the five PNGs, and
it reads **only** `../data/espoo/` (the channel table, the three committed references
of the upstream project and its own findings lock) — no network, no database, no
cache, no payload file.

| file | panel content |
| --- | --- |
| `espoo/fig_espoo_ocv_paper_A.png` | (a)-(b) the raw distributions of the two measured channels `gamma2` and `A` per orbit, per-acquisition dots under the box; (c) the **phase ladder** for the 80 of 150 acquisitions that carry it (and the note that `phase_snr_db` exists for 0) |
| `espoo/fig_espoo_ocv_paper_B.png` | (a) the pooled effect sizes ASC vs. DESC — Cliff's delta (DESC − ASC) with its bootstrap 95 % CI and SDS for `gamma2`, `A`, the two phase columns and the two weather co-variates; (b) the same contrast recomputed **per season**, the strata of the same sign |
| `espoo/fig_espoo_ocv_paper_C.png` | (a) the 2-class LDA over the whole record: leave-one-out balanced accuracy against the label-permutation null (chance **0.5**, the two-state counterpart of LUMO's four-state 0.25); (b) the LOO counts of the headline model — where the orbit is read back wrong |
| `espoo/fig_espoo_ocv_paper_D.png` | (a) the out-of-fold SDS of the four feature sets with the paired CI of the difference against `gamma2` alone; (b) the **out-of-time** split (75 / 75 by acquisition order, cut at 2025-09-06) — `A` and the weather co-variates add nothing that survives it |
| `espoo/fig_espoo_ocv_paper_E.png` | the four no-event controls: **E1** the 14 committed phase-coherence series (13 of 14 stationary, the flagged one marked), **E2** the month/season composition (nothing is cut on — there is no event), **E3** the same effect size recomputed *inside* each distinct wind stratum, with the pooled Spearman rho on the title line, **E4** the phase-ladder availability per month |

```bash
bash figures/espoo/fig_espoo_ocv_paper.sh                   # ~4.5 min, figures + report + pins
bash figures/espoo/fig_espoo_ocv_paper.sh --quick           # ~1.5 min, no permutation nulls
bash figures/espoo/fig_espoo_ocv_paper.sh --write-reference # (once) define the package lock
PYTHON=/usr/bin/python3 bash figures/espoo/fig_espoo_ocv_paper.sh
```

## `espoo/fig_espoo_ocv_months.sh` → `../code/espoo/fig_espoo_ocv_months.py`

The **standalone companion** of the Espoo paper figure, on the *time axis* instead of
the two pooled orbits — and, because the record is event-free, on the **whole
record**: nothing is cut, anchored or windowed, so all **150** acquisitions of the
**25** months (exactly 6 each, 2024-09 … 2026-09) appear. It reads only the committed
channel table and `../data/espoo/reference/espoo_phase_stationarity.json`, adds **no**
number to the A–E pin contract, and writes
[`espoo/fig_espoo_ocv_months.md`](espoo/fig_espoo_ocv_months.md),
[`../data/espoo/fig_espoo_ocv_months.json`](../data/espoo/fig_espoo_ocv_months.json)
and one PNG:

| file | panel content |
| --- | --- |
| `espoo/fig_espoo_ocv_months.png` | 2 x 3: (a) `gamma^2` by month, per orbit, per-acquisition dots with the monthly median; (b) the same for `A`; (c) the per-month Cliff's delta with its 2000-draw bootstrap CI, quoted only for the months in which both orbits carry ≥ 3 acquisitions (20 of 25 — the other 5 are reported as not evaluable rather than quoted on two rows); (d) the two primary stationarity series `3-8 Hz · strip400_pipeline` ASC / DESC; (e) coverage and the phase ladder by month (the two availability controls); (f) the weather co-variates by month (wind, temperature) |

Its own pin block re-derives the two primary series (`ASC` 0.2474 / `DESC` 0.2165 as
literal claims to four decimals), the 25 committed monthly means of both, the trend
signs, the per-orbit pooled medians, the 25 x 6 cadence and the phase ladder, and it
exits non-zero on any deviation.

```bash
bash figures/espoo/fig_espoo_ocv_months.sh            # ~4 s
python3 code/espoo/fig_espoo_ocv_months.py --figdir /tmp/x --json /tmp/y.json --md /tmp/z.md
```

