# `figures/` — one subfolder per site: the wrapper plus the report and the figures

```
figures/
├── lumo/      LUMO lattice tower  — fig_lumo_ocv_paper.sh|.md|[A-E].png
├── bautzen/   Bautzen OpenLabs bridge — fig_bautzen_ocv_paper.sh|.md|[A-E].png
└── README.md
```

The plotting scripts themselves live in `../code/<site>/` (they share that
folder's Python conventions and import each other); the wrapper in each site
folder runs them from the repository root, so that `data/<site>/…` and
`figures/<site>/…` resolve. Each script reads only its own `data/<site>/`
(channel table + reference JSONs) and writes its JSON artifact back to
`data/<site>/`, the report and the images here. No network, no database, no
burst cache.

```bash
bash figures/lumo/fig_lumo_ocv_paper.sh            # full run, ~5 min
bash figures/lumo/fig_lumo_ocv_paper.sh --quick    # ~2 min, permutation tests copied
bash figures/bautzen/fig_bautzen_ocv_paper.sh      # ~6 min, deck-edge mask

# if the python3 of PATH lacks numpy/matplotlib, name the interpreter:
PYTHON=/usr/bin/python3 bash figures/lumo/fig_lumo_ocv_paper.sh
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

