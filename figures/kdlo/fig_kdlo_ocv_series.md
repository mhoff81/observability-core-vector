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

The six dimensions exist together on only **one** of those stretches. `gamma2` and the mask size `A` exist on both — the era's backend rows carry `coherence_gamma2` and `coherence_masked_pixels`, the same columns the record carries. `P`, `D`, `F` and `S` are recomputed from the committed full-dwell strips, which exist only for the record: the era has no mask geometry, and its rows have no phase or frequency layer either.

## The four stretches of the axis

| stretch | first day | last day | calendar days | days with data | dimensions that exist here |
| --- | --- | --- | --- | --- | --- |
| the ascending IW2 era | 2015-03-20 | 2021-12-19 | 2467 | 162 | gamma2, A |
| the blind gap | 2021-12-20 | 2022-02-25 | 68 | 0 | – |
| in the catalogue, not in the record | 2022-02-26 | 2022-06-01 | 96 | 0 | – |
| the record's pre-collapse stretch | 2022-06-02 | 2022-12-11 | 193 | 16 | gamma2, P, D, A, F, S |

## How the axis is read

Each stretch is divided by **its own** median: era `gamma2` 0.00995, era `A` 91.5 over the 162 era days; record `gamma2` 0.003698, record `A` 120.5 over the 16 record days of the axis. The 1.0 line is therefore **two different levels in one panel**, and the two colours of panels (a) and (b) are two series of self-ratios, not one series with a step in it. Levels are not comparable across the break by construction (different track, different swath, different pass); what each stretch *does* is.

The record-only ratio for the four mask-geometry dimensions uses the same 16-day median of the record (2022-06-02..2022-12-11). Log axes where a dimension spans decades, so a value below 1 stays visible.

## The record inside the axis, day by day

| day | wind m/s | T degC | gamma2 | P | D | A | F | S | gamma2 /med | A /med |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 2022-06-02 | 5.22 | 17.6 | 0.02936 | 0.2313 | 0.04545 | 5 | 2 | 4.952 | 7.94 | 0.0415 |
| 2022-06-14 | 3.35 | 23.1 | 0.00386 | 0.05692 | 0.0143 | 56 | 27 | 66.77 | 1.04 | 0.465 |
| 2022-07-08 | 3.32 | 24.6 | 0.0006529 | 0.04918 | 0.0402 | 176 | 108 | 105.5 | 0.177 | 1.46 |
| 2022-07-20 | 5.73 | 24.8 | 0.001189 | 0.04974 | 0.04466 | 196 | 94 | 17.45 | 0.322 | 1.63 |
| 2022-08-01 | 1.14 | 21.7 | 0.004081 | 0.04019 | 0.04153 | 180 | 114 | 168.1 | 1.1 | 1.49 |
| 2022-08-13 | 4.16 | 22.7 | 0.0008671 | 0.04131 | 0.0553 | 236 | 112 | 179.6 | 0.234 | 1.96 |
| 2022-08-25 | 3.49 | 18.4 | 0.003328 | 0.03877 | 0.04818 | 212 | 105 | 3.045 | 0.9 | 1.76 |
| 2022-09-06 | 0.41 | 24.1 | 0.003536 | 0.03906 | 0.0295 | 122 | 74 | 115.5 | 0.956 | 1.01 |
| 2022-09-18 | 2.92 | 17.9 | 0.009327 | 0.04015 | 0.05496 | 237 | 115 | 146.3 | 2.52 | 1.97 |
| 2022-09-30 | 6.96 | 18.5 | 0.001126 | 0.05188 | 0.03557 | 153 | 72 | 137.6 | 0.304 | 1.27 |
| 2022-10-12 | 7.22 | 7.4 | 0.00498 | 0.05606 | 0.0229 | 97 | 52 | 143.1 | 1.35 | 0.805 |
| 2022-10-24 | 10.3 | 6.8 | 0.0009811 | 0.06757 | 0.01847 | 77 | 45 | 25.17 | 0.265 | 0.639 |
| 2022-11-05 | 3 | 1.3 | 0.01753 | 0.07528 | 0.0133 | 55 | 39 | 53.53 | 4.74 | 0.456 |
| 2022-11-17 | 8.66 | -6.9 | 0.006837 | 0.1008 | 0.004735 | 20 | 14 | 101.7 | 1.85 | 0.166 |
| 2022-11-29 | 7.58 | -7.3 | 0.000473 | 0.05515 | 0.02803 | 119 | 77 | 84.45 | 0.128 | 0.988 |
| 2022-12-11 | 5.6 | -3.2 | 0.01993 | 0.05547 | 0.01377 | 60 | 42 | 46.8 | 5.39 | 0.498 |

The record's own window starts 2022-06-02: the 96 descending days between the blind gap and it (2022-02-26..2022-06-01) exist in the catalogue but carry no committed row in this record, so they are hatched, not filled.

## The era inside the axis, year by year (7 years)

| year | days | first … last | gamma2 median | A median | gamma2 /era med | A /era med |
| --- | --- | --- | --- | --- | --- | --- |
| 2015 | 6 | 2015-03-20 … 2015-09-16 | 0.01055 | 152.5 | 1.06 | 1.67 |
| 2016 | 19 | 2016-01-26 … 2016-12-21 | 0.01439 | 58 | 1.45 | 0.634 |
| 2017 | 26 | 2017-01-02 … 2017-12-28 | 0.0158 | 86.5 | 1.59 | 0.945 |
| 2018 | 28 | 2018-01-09 … 2018-12-23 | 0.005836 | 161.5 | 0.587 | 1.77 |
| 2019 | 29 | 2019-01-04 … 2019-12-30 | 0.006206 | 131 | 0.624 | 1.43 |
| 2020 | 26 | 2020-01-11 … 2020-12-24 | 0.0216 | 60.5 | 2.17 | 0.661 |
| 2021 | 28 | 2021-01-05 … 2021-12-19 | 0.01531 | 61 | 1.54 | 0.667 |

These year medians are the only aggregation drawn over the era (bars in panels a and b), and they stay **inside** the era: a year is never compared with the record. They are era-track levels.

## Why panels (c)-(f) are empty on the era side

A line, a dashed interpolation or a greyed continuation over 2015-03-20..2021-12-19 in panels (c)-(f) would be an invention. `P`, `D`, `F` and `S` are recomputed here from the committed full-dwell strips, and `data/kdlo/strips/` holds the 60 acquisitions of the **descending** record only; the era's backend rows carry no phase or frequency layer at all (their non-null columns are exactly the ten of the export, checked against the pinned reference before anything is drawn). The era is therefore drawn in those four panels as a hatched region read **not measurable here** — with the reason written inside it — and the figure carries no value of those four dimensions before 2022-06-02.

## The collapse cut and the first collapsed acquisition

The axis stops at the collapse, **2022-12-11** — the last intact acquisition of the record — and the record's last 16 days run up to it. The next acquisition, **2022-12-23**, belongs to *neither* epoch: it is the first collapsed scene and lies outside this axis, so it is drawn hollow at the right edge in every panel (labelled in panel a).

| dimension | 2022-12-23 | /record median | 2022-12-11 (last intact) | /record median |
| --- | --- | --- | --- | --- |
| gamma2 | 0.001339 | 0.362 | 0.01993 | 5.39 |
| P | 0.05556 | 1.04 | 0.05547 | 1.04 |
| D | 0.01812 | 0.557 | 0.01377 | 0.423 |
| A | 63 | 0.523 | 60 | 0.498 |
| F | 23 | 0.315 | 42 | 0.575 |
| S | 36.69 | 0.394 | 46.8 | 0.503 |

On 2022-12-23 5 of the six dimensions sit **below** their own record median and 1 above; `gamma2` is 0.362× the record median, `A` is 0.523×. That one scene is the only post-collapse acquisition this figure draws, and it is drawn as a single point: n = 1 is a point, not a level.

## What this figure does not say

It does not say that the era "was" quiet and the record "declined": the two stretches are different geometries (ASCENDING|136 · IW2 · catalogue hour 00 vs. DESCENDING|12 · IW1 · one per day), the season coverage differs, and no statistic here is computed over both. It does not extend any of the six dimensions across the blind gap. It does not carry the collapse date's meaning into the era: the cut is drawn where the record's last intact acquisition is, not where a cause is.

## Inputs (all committed; all re-read before they are drawn)

| input | path | sha256 (first 16) |
| --- | --- | --- |
| era_export | `data/kdlo/fig_kdlo_0_raw.json` | `76ab24669cce671d…` |
| channel_table | `data/kdlo/kdlo_ocv_channels.csv` | `21dfe5b649c535ad…` |
| availability_reference | `data/kdlo/reference/kdlo_s1_availability.json` | `8d56dd42b79d7efc…` |

The era export's sha256 is the one the availability reference pins (`pipeline.figure_0_raw_export`), and this figure re-verifies it before using a single era value; the geometry break (2021-12-19 → 2022-02-26, 69 days) is re-read from the same reference and cross-checked against the packaged constants. The figure carries **84** self-checks, all of which abort before any file is written.

```bash
python3 code/kdlo/kdlo_ocv_channels_csv.py --verify-only   # the table
python3 code/kdlo/fig_kdlo_ocv_series.py                   # this figure
```
