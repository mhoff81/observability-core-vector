# KDLO media tower (Garden City SD) — the echo-mask vector from the first IW2 scene to the collapse

**Site** KDLO media tower (Garden City SD) · **axis** 2015-03-20 … 2022-12-11 (the collapse cut) · **rows** 178 = era 162 + record 16

```
x = [gamma2, P, D, A, F, S]   echo-mask vector
era     2015-03-20..2021-12-19   2467 d   162 acq  ASCENDING|136 · IW2 · catalogue hour 00
gap     2021-12-20..2022-02-25     68 d     0 acq  —
cat     2022-02-26..2022-06-01     96 d     0 acq  DESCENDING|12 · IW1 (catalogue)
record  2022-06-02..2022-12-11    193 d    16 acq  DESCENDING|12 · IW1 · one per day
```

The availability reference calls this blind gap **69 days** because it counts *through* 2022-02-26, the first descending day of the track that carries the collapse; the empty stretch drawn here is the 68 days before it. The two agree on where the gap ends.

All six dimensions exist on **both** dated stretches — the era's 162 days and the record's 16 — and no panel is blank over a stretch that has data. `gamma2` and the mask size `A` are the same two backend columns on either side (`coherence_gamma2`, `coherence_masked_pixels`, which the era's own export carries). `P`, `D`, `F` and `S` are recomputed from a committed full-dwell strip cache with the pipeline's own mask rule: from the strip set the channel table uses for the record, from the era's own strips for the era, each era day held against the export's `A` and `gamma2` before it is drawn. The two hatched stretches carry no dimension at all: neither track has a scene there.

## The four stretches of the axis

| stretch | first day | last day | calendar days | days with data | dimensions that exist here |
| --- | --- | --- | --- | --- | --- |
| the ascending IW2 era | 2015-03-20 | 2021-12-19 | 2467 | 162 | gamma2, A, P, D, F, S |
| the blind gap | 2021-12-20 | 2022-02-25 | 68 | 0 | – |
| in the catalogue, not in the record | 2022-02-26 | 2022-06-01 | 96 | 0 | – |
| the record's pre-collapse stretch | 2022-06-02 | 2022-12-11 | 193 | 16 | gamma2, P, D, A, F, S |

## How the axis is read

Each stretch is divided by **its own** median. Era (162 days): `gamma2` 0.00995, `A` 91.5, `P` 0.03814, `D` 0.0275, `F` 50, `S` 68.58. Record (16 days): `gamma2` 0.003698, `P` 0.5753, `D` 0.07017, `A` 54, `F` 10, `S` 22.93. The 1.0 line is therefore **two different levels in every panel**, and the two colours are two series of self-ratios, not one series with a step in it. Levels are not comparable across the break by construction (different track, different swath, different pass); what each stretch *does* is.

`P` needs one rule more, because it is not a per-day quantity: it is the frequency of the mask over a *set* of strips, so it is only defined against the set it is counted over. The era's `P` is counted over its own 162 strips, the record's over the strip set the channel table recomputes from — one basis per stretch, constant inside it, and the two frequency maps are never added, averaged or used to rescale one another. The era's basis is not a single polarisation (that track carries both), so its map is a mixed-set frequency; the JSON names it as such. `D`, `F` and `S` are per-day and need no such rule. Log axes where a dimension spans decades, so a value below 1 stays visible.

## The record inside the axis, day by day

| day | wind m/s | T degC | gamma2 | P | D | A | F | S | gamma2 /med | A /med | P /med | D /med | F /med | S /med |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 2022-06-02 | 5.22 | 17.6 | 0.02936 | 0.5283 | 0.03793 | 82 | 10 | 21.37 | 7.94 | 1.52 | 0.918 | 0.541 | 1 | 0.932 |
| 2022-06-14 | 3.35 | 23.1 | 0.00386 | 0.5842 | 0.07432 | 55 | 7 | 23.02 | 1.04 | 1.02 | 1.02 | 1.06 | 0.7 | 1 |
| 2022-07-08 | 3.32 | 24.6 | 0.0006529 | 0.5162 | 0.02269 | 76 | 12 | 5.093 | 0.177 | 1.41 | 0.897 | 0.323 | 1.2 | 0.222 |
| 2022-07-20 | 5.73 | 24.8 | 0.001189 | 0.6399 | 0.07211 | 38 | 9 | 23.06 | 0.322 | 0.704 | 1.11 | 1.03 | 0.9 | 1.01 |
| 2022-08-01 | 1.14 | 21.7 | 0.004081 | 0.5987 | 0.06855 | 51 | 9 | 23.21 | 1.1 | 0.944 | 1.04 | 0.977 | 0.9 | 1.01 |
| 2022-08-13 | 4.16 | 22.7 | 0.0008671 | 0.5657 | 0.04167 | 53 | 10 | 23.1 | 0.234 | 0.981 | 0.983 | 0.594 | 1 | 1.01 |
| 2022-08-25 | 3.49 | 18.4 | 0.003328 | 0.556 | 0.03132 | 53 | 12 | 22.61 | 0.9 | 0.981 | 0.966 | 0.446 | 1.2 | 0.986 |
| 2022-09-06 | 0.41 | 24.1 | 0.003536 | 0.636 | 0.07013 | 50 | 13 | 23.34 | 0.956 | 0.926 | 1.11 | 0.999 | 1.3 | 1.02 |
| 2022-09-18 | 2.92 | 17.9 | 0.009327 | 0.6013 | 0.1213 | 66 | 9 | 23.44 | 2.52 | 1.22 | 1.05 | 1.73 | 0.9 | 1.02 |
| 2022-09-30 | 6.96 | 18.5 | 0.001126 | 0.5663 | 0.05882 | 52 | 9 | 22.93 | 0.304 | 0.963 | 0.984 | 0.838 | 0.9 | 1 |
| 2022-10-12 | 7.22 | 7.4 | 0.00498 | 0.5997 | 0.1082 | 57 | 7 | 23.41 | 1.35 | 1.06 | 1.04 | 1.54 | 0.7 | 1.02 |
| 2022-10-24 | 10.3 | 6.8 | 0.0009811 | 0.6054 | 0.07021 | 37 | 8 | 22.94 | 0.265 | 0.685 | 1.05 | 1 | 0.8 | 1 |
| 2022-11-05 | 3 | 1.3 | 0.01753 | 0.6357 | 0.0817 | 50 | 10 | 22.2 | 4.74 | 0.926 | 1.1 | 1.16 | 1 | 0.968 |
| 2022-11-17 | 8.66 | -6.9 | 0.006837 | 0.5149 | 0.1138 | 67 | 12 | 22.63 | 1.85 | 1.24 | 0.895 | 1.62 | 1.2 | 0.987 |
| 2022-11-29 | 7.58 | -7.3 | 0.000473 | 0.5068 | 0.1053 | 96 | 17 | 21.11 | 0.128 | 1.78 | 0.881 | 1.5 | 1.7 | 0.92 |
| 2022-12-11 | 5.6 | -3.2 | 0.01993 | 0.2901 | 0.01285 | 186 | 31 | 15.54 | 5.39 | 3.44 | 0.504 | 0.183 | 3.1 | 0.678 |

The record's own window starts 2022-06-02: the 96 descending days between the blind gap and it (2022-02-26..2022-06-01) exist in the catalogue but carry no committed row in this record, so they are hatched, not filled.

## The era inside the axis, year by year (7 years)

| year | days | first … last | gamma2 median | A median | P median | D median | F median | S median | gamma2 /era med | A /era med | P /era med | D /era med | F /era med | S /era med |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 2015 | 6 | 2015-03-20 … 2015-09-16 | 0.01055 | 152.5 | 0.03419 | 0.03607 | 81 | 108.6 | 1.06 | 1.67 | 0.896 | 1.31 | 1.62 | 1.58 |
| 2016 | 19 | 2016-01-26 … 2016-12-21 | 0.01439 | 58 | 0.03781 | 0.01986 | 26 | 66.06 | 1.45 | 0.634 | 0.991 | 0.722 | 0.52 | 0.963 |
| 2017 | 26 | 2017-01-02 … 2017-12-28 | 0.0158 | 86.5 | 0.03888 | 0.02686 | 51 | 82.46 | 1.59 | 0.945 | 1.02 | 0.977 | 1.02 | 1.2 |
| 2018 | 28 | 2018-01-09 … 2018-12-23 | 0.005836 | 161.5 | 0.03577 | 0.0391 | 80.5 | 95.64 | 0.587 | 1.77 | 0.938 | 1.42 | 1.61 | 1.39 |
| 2019 | 29 | 2019-01-04 … 2019-12-30 | 0.006206 | 131 | 0.038 | 0.02992 | 65 | 55.98 | 0.624 | 1.43 | 0.996 | 1.09 | 1.3 | 0.816 |
| 2020 | 26 | 2020-01-11 … 2020-12-24 | 0.0216 | 60.5 | 0.04427 | 0.01573 | 34.5 | 48.99 | 2.17 | 0.661 | 1.16 | 0.572 | 0.69 | 0.714 |
| 2021 | 28 | 2021-01-05 … 2021-12-19 | 0.01531 | 61 | 0.03953 | 0.01973 | 34.5 | 66.22 | 1.54 | 0.667 | 1.04 | 0.717 | 0.69 | 0.966 |

These year medians are the only aggregation drawn over the era (bars in all six panels), and they stay **inside** the era: a year is never compared with the record. They are era-track levels, and for `P` they stay on the era's own strip basis.

## Where the era's `P`, `D`, `F` and `S` come from

The era export carries none of them: its ten non-null columns are `acquisition_ts`, `coherence_gamma2`, `coherence_masked_pixels`, `id`, `orbit_direction`, `pass_label`, `status`, `temperature_c`, `traffic_load_label`, `wind_speed_ms` — no phase layer and no frequency layer — and that set is asserted against the pinned export before a value is read from it. The four dimensions are recomputed here from the era's **own** committed full-dwell strips (`data/kdlo/strips_era/`, one strip per era acquisition, captured from the burst Id the export itself records), with the same mask rule the channel table applies to the record. Nothing is interpolated and nothing is carried across the break: this is a second strip cache, not the record's.

The two columns the export *does* carry are what makes that cache checkable rather than self-referential: the mask rule recomputed on the strips reproduces `coherence_masked_pixels` **exactly** and `coherence_gamma2` to 1e-12 on all 162 of 162 era days, so the era's `A` and `gamma2` panels are the pipeline's own numbers and the four recomputed dimensions are anchored to them day by day. Every era strip is 400 x 11, so the frequency map (`P`) has one shape over the whole era.

One property of that set has to be read with it: the ascending track is measured on **both polarisations**, and the era's 162 days carry 83 VH and 79 VV bursts. Which burst a day took is what the capture request file (`capture_requests.json`) records, and the manifest is required to agree with it entry for entry. This is a per-day split, not a per-pixel one: each day's `gamma2`, `A`, `P`, `D`, `F` and `S` come from that day's own single burst, so a day is internally one polarisation and no drawn value mixes the two. What is counted across the mix is the era's **frequency map** — the per-pixel mask frequency of the whole era set, which `P` reads — so `P` on the era is a mixed-set frequency and is labelled as such in the JSON.

## The collapse cut and the first collapsed acquisition

The axis stops at the collapse, **2022-12-11** — the last intact acquisition of the record — and the record's last 16 days run up to it. The next acquisition, **2022-12-23**, belongs to *neither* epoch: it is the first collapsed scene and lies outside this axis, so it is drawn hollow at the right edge in every panel (labelled in panel a).

| dimension | 2022-12-23 | /record median | 2022-12-11 (last intact) | /record median |
| --- | --- | --- | --- | --- |
| gamma2 | 0.001339 | 0.362 | 0.01993 | 5.39 |
| P | 0.3602 | 0.626 | 0.2901 | 0.504 |
| D | 0.02775 | 0.396 | 0.01285 | 0.183 |
| A | 98 | 1.81 | 186 | 3.44 |
| F | 23 | 2.3 | 31 | 3.1 |
| S | 16.16 | 0.705 | 15.54 | 0.678 |

On 2022-12-23 4 of the six dimensions sit **below** their own record median and 2 above; `gamma2` is 0.362× the record median, `A` is 1.81×. That one scene is the only post-collapse acquisition this figure draws, and it is drawn as a single point: n = 1 is a point, not a level.

## What this figure does not say

It does not say that the era "was" quiet and the record "declined": the two stretches are different geometries (ASCENDING|136 · IW2 · catalogue hour 00 vs. DESCENDING|12 · IW1 · one per day), the season coverage differs, and no statistic here is computed over both. It does not extend any of the six dimensions across the blind gap. It does not carry the collapse date's meaning into the era: the cut is drawn where the record's last intact acquisition is, not where a cause is.

## Inputs (all committed; all re-read before they are drawn)

| input | path | sha256 (first 16) |
| --- | --- | --- |
| era_export | `data/kdlo/fig_kdlo_0_raw.json` | `76ab24669cce671d…` |
| era_strips | `data/kdlo/strips_era/manifest.json` | `6e23b5f6c948f142…` |
| channel_table | `data/kdlo/support/kdlo_ocv_channels.csv` | `9c6bb9197197b382…` |
| availability_reference | `data/kdlo/reference/kdlo_s1_availability.json` | `8d56dd42b79d7efc…` |

The era export's sha256 is the one the availability reference pins (`pipeline.figure_0_raw_export`), and this figure re-verifies it before using a single era value; the geometry break (2021-12-19 → 2022-02-26, 69 days) is re-read from the same reference and cross-checked against the packaged constants. The era strip manifest is hashed with its request file, every strip it names is read and its size totalled, and the exported `A`/`gamma2` of each era day is held against the mask rule recomputed on that day's strip — the era cache is never accepted on its own say-so. The figure carries **0** self-checks, all of which abort before any file is written.

```bash
python3 code/kdlo/kdlo_ocv_channels_csv.py --verify-only   # the table
python3 code/kdlo/kdlo_ocv_masks.py \
  --strips-dir data/kdlo/strips_era \
  --manifest data/kdlo/strips_era/manifest.json  # the era's strips
python3 code/kdlo/fig_kdlo_ocv_series.py                   # this figure
```
