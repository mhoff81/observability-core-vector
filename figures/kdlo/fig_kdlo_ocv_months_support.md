# KDLO media tower (Garden City SD) — the echo-mask dimensions month by month over the record

**Site** KDLO media tower (Garden City SD) · **record** 2022-06 … 2024-09 · **acquisitions** 60

```
x     = [gamma2, P, D, A, F, S]      echo-mask vector
mask  = intensity >= 0.30*peak AND >= 5.0*np.median(intensity), A >= 2
pre     2022-06-02..2022-12-11  pre-rebuild (standing mast) 16 acquisitions
single  2022-12-23              in neither epoch    1 acquisition
rebuild 2023-05-04..2024-09-19  post-rebuild (reconstruction) 43 acquisitions
```

## Why this is not a 24-month series

The other site companions roll 24 one-month buckets before their event because those records carry a dense pre-collapse stretch. The KDLO record is 60 acquisitions in two **disjoint epochs** (``pre`` 2022-06-02..2022-12-11 and ``rebuild`` 2023-05-04..2024-09-19), so the axis here is the whole record by **calendar month** (28 buckets, 4 of them with no acquisition). Two properties of the record are kept visible rather than smoothed away: the collapse month **2022-12 is mixed** (it holds both the last ``pre`` acquisitions, to 2022-12-11, and the ``single`` 2022-12-23 that belongs to neither epoch), and the record has a four-month acquisition gap between the two epochs.

## The record, month by month (28 calendar buckets)

| month | epoch(s) | acquisitions | first day | last day | days |
| --- | --- | --- | --- | --- | --- |
| 2022-06 | pre | 2 | 2022-06-02 | 2022-06-14 | 2 |
| 2022-07 | pre | 2 | 2022-07-08 | 2022-07-20 | 2 |
| 2022-08 | pre | 3 | 2022-08-01 | 2022-08-25 | 3 |
| 2022-09 | pre | 3 | 2022-09-06 | 2022-09-30 | 3 |
| 2022-10 | pre | 2 | 2022-10-12 | 2022-10-24 | 2 |
| 2022-11 | pre | 3 | 2022-11-05 | 2022-11-29 | 3 |
| 2022-12 | pre+single | 2 | 2022-12-11 | 2022-12-23 | 2 |
| 2023-01 | (no acquisition) | 0 | – | – | 0 |
| 2023-02 | (no acquisition) | 0 | – | – | 0 |
| 2023-03 | (no acquisition) | 0 | – | – | 0 |
| 2023-04 | (no acquisition) | 0 | – | – | 0 |
| 2023-05 | rebuild | 3 | 2023-05-04 | 2023-05-28 | 3 |
| 2023-06 | rebuild | 2 | 2023-06-09 | 2023-06-21 | 2 |
| 2023-07 | rebuild | 3 | 2023-07-03 | 2023-07-27 | 3 |
| 2023-08 | rebuild | 2 | 2023-08-08 | 2023-08-20 | 2 |
| 2023-09 | rebuild | 3 | 2023-09-01 | 2023-09-25 | 3 |
| 2023-10 | rebuild | 3 | 2023-10-07 | 2023-10-31 | 3 |
| 2023-11 | rebuild | 2 | 2023-11-12 | 2023-11-24 | 2 |
| 2023-12 | rebuild | 3 | 2023-12-06 | 2023-12-30 | 3 |
| 2024-01 | rebuild | 2 | 2024-01-11 | 2024-01-23 | 2 |
| 2024-02 | rebuild | 3 | 2024-02-04 | 2024-02-28 | 3 |
| 2024-03 | rebuild | 2 | 2024-03-11 | 2024-03-23 | 2 |
| 2024-04 | rebuild | 3 | 2024-04-04 | 2024-04-28 | 3 |
| 2024-05 | rebuild | 2 | 2024-05-10 | 2024-05-22 | 2 |
| 2024-06 | rebuild | 3 | 2024-06-03 | 2024-06-27 | 3 |
| 2024-07 | rebuild | 2 | 2024-07-09 | 2024-07-21 | 2 |
| 2024-08 | rebuild | 3 | 2024-08-02 | 2024-08-26 | 3 |
| 2024-09 | rebuild | 2 | 2024-09-07 | 2024-09-19 | 2 |

`pre` months: 2022-06, 2022-07, 2022-08, 2022-09, 2022-10, 2022-11; collapse month: 2022-12; empty: 2023-01, 2023-02, 2023-03, 2023-04; `rebuild` months: 2023-05 … 2024-09.

## Median of each dimension per month

| month | n | epoch | gamma2 | P | D | A | F | S |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 2022-06 | 2 | pre | 0.01661 | 0.5562 | 0.05613 | 68.5 | 8.5 | 22.19 |
| 2022-07 | 2 | pre | 0.0009209 | 0.5781 | 0.0474 | 57 | 10.5 | 14.08 |
| 2022-08 | 3 | pre | 0.003328 | 0.5657 | 0.04167 | 53 | 10 | 23.1 |
| 2022-09 | 3 | pre | 0.003536 | 0.6013 | 0.07013 | 52 | 9 | 23.34 |
| 2022-10 | 2 | pre | 0.00298 | 0.6026 | 0.08918 | 47 | 7.5 | 23.17 |
| 2022-11 | 3 | pre | 0.006837 | 0.5149 | 0.1053 | 67 | 12 | 22.2 |
| 2022-12 | 2 | pre+single | 0.01063 | 0.3251 | 0.0203 | 142 | 27 | 15.85 |
| 2023-01 | 0 | – | – | – | – | – | – | – |
| 2023-02 | 0 | – | – | – | – | – | – | – |
| 2023-03 | 0 | – | – | – | – | – | – | – |
| 2023-04 | 0 | – | – | – | – | – | – | – |
| 2023-05 | 3 | rebuild | 0.02099 | 0.4533 | 0.06593 | 100 | 18 | 20.31 |
| 2023-06 | 2 | rebuild | 0.01215 | 0.5336 | 0.08488 | 73 | 11 | 13.3 |
| 2023-07 | 3 | rebuild | 0.006371 | 0.5933 | 0.05497 | 57 | 6 | 5.862 |
| 2023-08 | 2 | rebuild | 0.002206 | 0.4164 | 0.07368 | 134.5 | 33 | 33.43 |
| 2023-09 | 3 | rebuild | 0.01386 | 0.5074 | 0.07124 | 86 | 8 | 20.23 |
| 2023-10 | 3 | rebuild | 0.1912 | 0.5428 | 0.07212 | 76 | 10 | 21.94 |
| 2023-11 | 2 | rebuild | 0.1547 | 0.4996 | 0.04366 | 86.5 | 12 | 12.6 |
| 2023-12 | 3 | rebuild | 0.003263 | 0.4498 | 0.04905 | 83 | 14 | 18.4 |
| 2024-01 | 2 | rebuild | 0.02302 | 0.4584 | 0.04147 | 87 | 11.5 | 10.26 |
| 2024-02 | 3 | rebuild | 0.05936 | 0.476 | 0.05869 | 102 | 12 | 9.887 |
| 2024-03 | 2 | rebuild | 0.03592 | 0.4396 | 0.05433 | 97.5 | 14 | 6.813 |
| 2024-04 | 3 | rebuild | 0.04671 | 0.4975 | 0.07828 | 85 | 10 | 5.369 |
| 2024-05 | 2 | rebuild | 0.0108 | 0.4115 | 0.07947 | 72.5 | 12 | 5.649 |
| 2024-06 | 3 | rebuild | 0.04027 | 0.4575 | 0.09662 | 62 | 9 | 6.056 |
| 2024-07 | 2 | rebuild | 0.1705 | 0.4529 | 0.04559 | 63.5 | 15 | 6.06 |
| 2024-08 | 3 | rebuild | 0.05768 | 0.4842 | 0.05185 | 75 | 9 | 5.59 |
| 2024-09 | 2 | rebuild | 0.007545 | 0.4561 | 0.01451 | 112 | 18.5 | 28.84 |

`baseline` = the `2022-06` (first record month) medians — the 1.0 line of panel f: `gamma2` 0.01661, `P` 0.5562, `D` 0.05613, `A` 68.5, `F` 8.5, `S` 22.19.

## The 6D fingerprint at the cut and at the epoch edge (month / baseline)

| month | gamma2 | P | D | A | F | S |
| --- | --- | --- | --- | --- | --- | --- |
| 2022-06 | 1 | 1 | 1 | 1 | 1 | 1 |
| 2022-11 | 0.4116 | 0.9257 | 1.875 | 0.9781 | 1.412 | 1 |
| 2022-12 | 0.6401 | 0.5845 | 0.3617 | 2.073 | 3.176 | 0.7141 |
| 2023-05 | 1.263 | 0.815 | 1.175 | 1.46 | 2.118 | 0.9152 |
| 2024-09 | 0.4542 | 0.82 | 0.2585 | 1.635 | 2.176 | 1.299 |

## The collapse month against its baselines

| dimension | δ 2022-12 vs prev | CI≠0 | δ 2022-12 vs pooled pre | CI≠0 | δ rebuild vs pre | CI≠0 |
| --- | --- | --- | --- | --- | --- | --- |
| gamma2 | 0.3333 | no | 0.2812 | no | 0.5291 | yes |
| P | -1 | yes | -0.9062 | yes | -0.6279 | yes |
| D | -1 | yes | -0.8438 | yes | -0.2006 | no |
| A | 1 | yes | 0.9062 | yes | 0.4985 | yes |
| F | 1 | yes | 0.9062 | yes | 0.1206 | no |
| S | -1 | yes | -0.7812 | yes | -0.6395 | yes |

`prev` = the last month before the collapse (`2022-11`); `pooled pre` = the 6 pure `pre` months; `rebuild vs pre` is the paper's epoch contrast. Cliff's delta > 0 means the right-hand set is larger. Full CIs, medians and SDS are in the JSON (`deltas[dim].*`).

## Reading (descriptive — not a claim of damage)

The 6 pure `pre` months carry **16** acquisitions (two to three per month) and the collapse month **2**, so every monthly interval is very wide: this is a descriptive companion, not a test.

Relative to the first record month `2022-06` (= 1.0) the last pure `pre` month `2022-11` sits at `gamma2` 0.4116, `P` 0.9257, `D` 1.875, `A` 0.9781, `F` 1.412, `S` 1 and the first `rebuild` month `2023-05` at `gamma2` 1.263, `P` 0.815, `D` 1.175, `A` 1.46, `F` 2.118, `S` 0.9152. Against the last month before the collapse the mixed bucket clears the bootstrap CI in **5** of the six dimensions — and that bucket is not a clean post month. The pooled epoch contrast clears it in **4** of six. The axis is read as the record's own co-variation over two disjoint windows, **never** as damage: the reconstruction, the season, the acquisition geometry and the pre/post mix move with the collapse date as well.

```bash
python3 code/kdlo/kdlo_ocv_channels_csv.py       # build + verify the table
python3 code/kdlo/fig_kdlo_ocv_months.py          # this figure
```
