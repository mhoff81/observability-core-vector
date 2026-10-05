# `figures/` — the `.sh` wrapper plus the `.md` report and the `.png` figures

The plotting script itself lives in `../code/` (it shares that directory's
Python conventions and imports); the wrapper here runs it with the right working
directory. The script reads only `data/lumo_ocv_channels.csv` and the five
`data/reference/*.json`, writes its JSON artifact to `data/`, and the report and
the five PNGs here. No network, no database, no burst cache.

```bash
./figures/fig_lumo_ocv_paper.sh            # full run, ~5 min
./figures/fig_lumo_ocv_paper.sh --quick    # ~2 min, permutation tests copied

# if the python3 of PATH lacks numpy/matplotlib, name the interpreter:
PYTHON=/usr/bin/python3 ./figures/fig_lumo_ocv_paper.sh
```

### `fig_lumo_ocv_paper.sh` → `../code/fig_lumo_ocv_paper.py`

Five figures on the question of whether the observability core vector
`OCV = [gamma2, P, D]` discriminates the LUMO structural states (healthy, DAM3,
DAM4, DAM6) over the 178 coherence overpasses, plus the English report
[`fig_lumo_ocv_paper.md`](fig_lumo_ocv_paper.md) and
`../data/fig_lumo_ocv_paper.json`:

| file | panel content |
| --- | --- |
| `fig_lumo_ocv_paper_A.png` | the six channels one panel each: raw per-overpass dots and the median bar for gamma^2, mask area A, density D, fragmentation F, centroid<->peak shift S and persistence P |
| `fig_lumo_ocv_paper_B.png` | in-sample effect sizes: SDS = 10&#124;Cliff's delta&#124; per component, pooled and per damage state, plus the ranking of the components |
| `fig_lumo_ocv_paper_C.png` | 4-class LDA: balanced accuracy over the shrinkage grid (6D vs. gamma^2 vs. single components), the leave-one-out confusion matrix and the label-permutation null against the observed accuracy |
| `fig_lumo_ocv_paper_D.png` | the six binary state pairs: out-of-fold Fisher-LDA scores for the four feature sets (gamma^2, [P,D,S], [gamma^2,P,D], 6D), with the SDS of the headline model and the permutation p per pair |
| `fig_lumo_ocv_paper_E.png` | controls: severity trend of gamma^2 (Spearman, Jonckheere-Terpstra), season/orbit strata, and the mask brightness / wind covariates per state |

The figures are 183 mm wide, 600 dpi PNGs with the same style as the Espoo
figures (DejaVu Sans + STIX). They are deterministic: fixed seeds, no timestamps
in any annotation, so two runs produce byte-identical files. The `.png` files are
git-ignored — the wrapper regenerates them from the committed data.

`fig_lumo_ocv_paper.md` is the same content in prose: question and design, data
and mask definition, the results per figure as tables, the pin verdict against
the committed references, the exact commands and the caveats.
