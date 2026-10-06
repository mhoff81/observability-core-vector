# Carola Bruecke, Dresden — observability core vector, pre/post collapse

**Site** Carola Bruecke, Dresden · **event** 2024-09-11 · **rows** 1717 80x80 girder chips · **states** `pre-collapse (healthy)` / `post-collapse`

```
OCV  = [gamma2, P, D]              observability core vector
x    = [gamma2, P, D, A, F, S]      echo-mask / deck-mask vector
mask = intensity >= 0.30*peak AND >= 5.0*np.median(intensity), A >= 2
layers: echo mask (whole chip)  |  deck mask (per girder A/B/C)
```

## What is recomputed, and what is pinned

| input | role |
| --- | --- |
| `data/carola/carola_windows_mask_cache.txt` | the echo mask of every committed chip (the mask layer's input) |
| `data/carola/carola_measurements_full.txt` | the committed measurement extract (channels, weather, geometry) |
| `data/carola/carola_segments.txt` | the girder labels + FEM baselines |
| `data/carola/carola_ocv_channels.csv` | generated: one row per chip, 57 columns |

The mask cache holds **1717** 80x80 chips (265 7x7 tower chips are excluded), after dropping **472** byte-identical duplicates (the per-girder re-extraction resolves the segments: `segment_resolved = True`, 264 acquisition groups with differing payloads).

### The two mask layers

* **echo mask** — the whole 80x80 chip; the layer of the 6D vector.
* **deck mask** — the same rule on the girder's own chip (`segment_index` of the export index, labelled by `carola_segments.txt`).

Neither layer's `A` is the pipeline's own `coherence_masked_pixels`: that column was computed on the pre-fix *asset-point* chip while the committed export holds the per-girder re-extraction, so the two agree only by coincidence (60 of 1713 rows). The site's analysis project therefore recomputed the mask from the payloads, and this package re-derives it from the committed cache.

## A — the channels

Raw distributions of the six mask channels per state (fig. A).

| channel | pre median | post median | delta | SDS | CI95 excludes 0 |
| --- | --- | --- | --- | --- | --- |
| gamma2 | 0.1228 | 0.2000 | 0.171 | 1.715 | yes |
| P | 0.0351 | 0.0122 | -0.661 | 6.606 | yes |
| D | 0.0066 | 0.0075 | 0.063 | 0.630 | yes |
| A | 13.0000 | 12.0000 | -0.100 | 0.999 | yes |
| F | 5.0000 | 4.0000 | -0.074 | 0.745 | yes |
| S | 19.3763 | 14.6161 | -0.114 | 1.139 | yes |

## B — the two layers side by side

Cliff's delta (pre -> post) of the deck-mask channels inside each girder:

| girder | n | pre | post | gamma2 | P | D | A | F | S |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Girder A | 1189 | 878 | 311 | 0.250 | -0.749 | 0.108 | -0.160 | -0.138 | -0.209 |
| Girder B | 264 | 9 | 255 | 0.331 | -0.335 | 0.114 | -0.081 | -0.141 | -0.049 |
| Girder C | 264 | 9 | 255 | 0.271 | -0.226 | -0.156 | -0.312 | -0.280 | -0.038 |

The reference's own per-girder table (recomputed, fig. B/C):

| girder | n | healthy | damaged | gamma2 median | A median | rho(gamma2, A) | repair target |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Girder A | 1189 | 878 | 311 | 0.1391 | 13.0 | -0.753 | yes |
| Girder B | 264 | 9 | 255 | 0.1752 | 13.0 | -0.806 | no |
| Girder C | 264 | 9 | 255 | 0.1471 | 13.5 | -0.813 | no |

## C — is the state recoverable?

2-class LDA on the 6D vector, leave-one-out at `lambda = 0.5`: accuracy **0.740**, balanced accuracy **0.741** (chance 0.5).

Label permutation (1000 draws): p(accuracy) = 0.0010, p(balanced) = 0.0010; null mean accuracy 0.500.

## D — does the extra structure help?

| model | features | SDS (OOF) | delta vs gamma2 | CI95 |
| --- | --- | --- | --- | --- |
| gamma2 | gamma2 | 1.638 | - | [-, -] |
| P+D | P+D | 6.483 | 4.843 | [4.158, 5.515] |
| gamma2+P+D | gamma2+P+D | 6.381 | 4.746 | [4.132, 5.340] |
| 6 dims | gamma2+P+D+A+F+S | 6.299 | 4.661 | [4.079, 5.233] |

## E — what else moves with the collapse date?

Stratified SDS (= 10|Cliff's delta|) of `gamma2` within each stratum:

| stratifier | stratum | n pre | n post | SDS gamma2 |
| --- | --- | --- | --- | --- |
| girder | Girder A | 878 | 311 | 2.501 |
| girder | Girder B | 9 | 255 | 3.307 |
| girder | Girder C | 9 | 255 | 2.706 |
| orbit | ASCENDING | 607 | 541 | 1.838 |
| orbit | DESCENDING | 289 | 280 | 1.715 |
| season | shoulder_Aug-Sep | 161 | 147 | 2.995 |
| season | summer_Mar-Jul | 424 | 348 | 1.720 |
| season | winter_Oct-Feb | 311 | 326 | 1.269 |

The mask-size co-variate (fig. E panel b) and the nested OLS (the reference's `model_block`):

* rho(gamma2, A) = -0.771 pooled, -0.844 ascending, -0.361 descending.
* `gamma2 ~ damaged`: r2 = 0.028, beta(damaged) = 0.0826 (HC3 p = 0.000).
* `gamma2 ~ A`: r2 = 0.361, beta(A) = -0.00887.
* `gamma2 ~ damaged + A`: r2 = 0.377; delta r2 = +0.349 from A and +0.016 from the state.
* inverse-N law: r2 = 0.595.
* SDS(gamma2) = 1.715, SDS(A) = 0.999, scenario **B** — SDS(n_masked) < SDS(gamma2), so a state effect on n_masked could reach gamma2 through the mask (an indirect channel).

## Pins

474 checks against the committed references (414 exact, 60 within 1e-12); **0 deviations**.

| reference | what is pinned |
| --- | --- |
| `carola_coherence_states.json` | per-window `gamma2` and `n_masked` (925 of its 1454 stored windows are in the current export; 925 match exactly) |
| `carola_echo_mask_gamma2.json` | the state table, the per-girder and per-request tables, the SDS block and the nested OLS — number by number |
| `carola_bridge_osm.json` | the map extract structurally (Carolabrücke B 170, Q1044279, two 12-node bridge ways, two razed DVB ways) |

### Reading (not a claim of damage)

The echo-mask coherence rises across the collapse (SDS 1.715), while the mask size `A` moves much less (SDS 0.999). Both states sit on the same bridge, so season, weather, traffic and pipeline generation move with the state too — figures C–E report that co-variation rather than a mechanism. What the figures support is the *observability* statement: the vector separates the two epochs of the record, and the 6D deck-mask vector does so per girder.

```bash
python3 code/carola/carola_ocv_masks.py --check      # the mask rule
python3 code/carola/carola_ocv_channels_csv.py       # build + verify
python3 code/carola/fig_carola_ocv_paper.py          # figures + pins
```
