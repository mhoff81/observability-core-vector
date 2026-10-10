# KDLO media tower (Garden City SD) — the echo-mask dimensions month by month over the record

**Site** KDLO media tower (Garden City SD) · **record** 2022-06 … 2024-09 · **acquisitions** 60

```
x     = [gamma2, P, D, A, F, S]      echo-mask vector
mask  = intensity >= 0.30*peak AND >= 5.0*np.median(intensity), A >= 2
pre     2022-06-02..2022-12-11  standing mast       16 acquisitions
single  2022-12-23              in neither epoch    1 acquisition
rebuild 2023-05-04..2024-09-19  during the rebuild  43 acquisitions
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
| 2022-06 | 2 | pre | 0.01661 | 0.1441 | 0.02988 | 30.5 | 14.5 | 35.86 |
| 2022-07 | 2 | pre | 0.0009209 | 0.04946 | 0.04243 | 186 | 101 | 61.49 |
| 2022-08 | 3 | pre | 0.003328 | 0.04019 | 0.04818 | 212 | 112 | 168.1 |
| 2022-09 | 3 | pre | 0.003536 | 0.04015 | 0.03557 | 153 | 74 | 137.6 |
| 2022-10 | 2 | pre | 0.00298 | 0.06181 | 0.02069 | 87 | 48.5 | 84.13 |
| 2022-11 | 3 | pre | 0.006837 | 0.07528 | 0.0133 | 55 | 39 | 84.45 |
| 2022-12 | 2 | pre+single | 0.01063 | 0.05551 | 0.01595 | 61.5 | 32.5 | 41.75 |
| 2023-01 | 0 | – | – | – | – | – | – | – |
| 2023-02 | 0 | – | – | – | – | – | – | – |
| 2023-03 | 0 | – | – | – | – | – | – | – |
| 2023-04 | 0 | – | – | – | – | – | – | – |
| 2023-05 | 3 | rebuild | 0.02099 | 0.07422 | 0.009193 | 36 | 17 | 34.73 |
| 2023-06 | 2 | rebuild | 0.01215 | 0.04391 | 0.03312 | 144 | 78.5 | 60.69 |
| 2023-07 | 3 | rebuild | 0.006371 | 0.04335 | 0.01818 | 79 | 53 | 172.3 |
| 2023-08 | 2 | rebuild | 0.002206 | 0.04118 | 0.04408 | 192 | 102 | 63.82 |
| 2023-09 | 3 | rebuild | 0.01386 | 0.04172 | 0.0422 | 162 | 100 | 53.67 |
| 2023-10 | 3 | rebuild | 0.1912 | 0.2773 | 1 | 4 | 1 | 0.7071 |
| 2023-11 | 2 | rebuild | 0.1547 | 0.1363 | 0.2189 | 48.5 | 29.5 | 20.95 |
| 2023-12 | 3 | rebuild | 0.003263 | 0.05434 | 0.02689 | 113 | 65 | 128.6 |
| 2024-01 | 2 | rebuild | 0.02302 | 0.0605 | 0.0218 | 93.5 | 50.5 | 124.3 |
| 2024-02 | 3 | rebuild | 0.05936 | 0.1328 | 0.04166 | 6 | 4 | 59.71 |
| 2024-03 | 2 | rebuild | 0.03592 | 0.08454 | 0.00932 | 31.5 | 20.5 | 14.07 |
| 2024-04 | 3 | rebuild | 0.04671 | 0.06368 | 0.01285 | 53 | 33 | 32.93 |
| 2024-05 | 2 | rebuild | 0.0108 | 0.05866 | 0.02681 | 113 | 54 | 138.7 |
| 2024-06 | 3 | rebuild | 0.04027 | 0.08136 | 0.01336 | 29 | 16 | 68.48 |
| 2024-07 | 2 | rebuild | 0.1705 | 0.1124 | 0.1474 | 100 | 51.5 | 47.19 |
| 2024-08 | 3 | rebuild | 0.05768 | 0.05208 | 0.04708 | 30 | 18 | 39.08 |
| 2024-09 | 2 | rebuild | 0.007545 | 0.04072 | 0.04002 | 163.5 | 93 | 125.4 |

`baseline` = the `2022-06` (first record month) medians — the 1.0 line of panel f: `gamma2` 0.01661, `P` 0.1441, `D` 0.02988, `A` 30.5, `F` 14.5, `S` 35.86.

## The 6D fingerprint at the cut and at the epoch edge (month / baseline)

| month | gamma2 | P | D | A | F | S |
| --- | --- | --- | --- | --- | --- | --- |
| 2022-06 | 1 | 1 | 1 | 1 | 1 | 1 |
| 2022-11 | 0.4116 | 0.5225 | 0.4451 | 1.803 | 2.69 | 2.355 |
| 2022-12 | 0.6401 | 0.3853 | 0.5338 | 2.016 | 2.241 | 1.164 |
| 2023-05 | 1.263 | 0.5151 | 0.3077 | 1.18 | 1.172 | 0.9685 |
| 2024-09 | 0.4542 | 0.2826 | 1.34 | 5.361 | 6.414 | 3.496 |

## The collapse month against its baselines

| dimension | δ 2022-12 vs prev | CI≠0 | δ 2022-12 vs pooled pre | CI≠0 | δ rebuild vs pre | CI≠0 |
| --- | --- | --- | --- | --- | --- | --- |
| gamma2 | 0.3333 | no | 0.2812 | no | 0.5291 | yes |
| P | -0.3333 | no | 0.2188 | no | 0.2355 | no |
| D | 0.3333 | no | -0.5938 | yes | -0.06977 | no |
| A | 0.3333 | no | -0.4062 | no | -0.2805 | no |
| F | 0 | no | -0.5938 | yes | -0.2674 | no |
| S | -1 | yes | -0.4688 | no | -0.1453 | no |

`prev` = the last month before the collapse (`2022-11`); `pooled pre` = the 6 pure `pre` months; `rebuild vs pre` is the paper's epoch contrast. Cliff's delta > 0 means the right-hand set is larger. Full CIs, medians and SDS are in the JSON (`deltas[dim].*`).

## Reading (descriptive — not a claim of damage)

The 6 pure `pre` months carry **16** acquisitions (two to three per month) and the collapse month **2**, so every monthly interval is very wide: this is a descriptive companion, not a test.

Relative to the first record month `2022-06` (= 1.0) the last pure `pre` month `2022-11` sits at `gamma2` 0.4116, `P` 0.5225, `D` 0.4451, `A` 1.803, `F` 2.69, `S` 2.355 and the first `rebuild` month `2023-05` at `gamma2` 1.263, `P` 0.5151, `D` 0.3077, `A` 1.18, `F` 1.172, `S` 0.9685. Against the last month before the collapse the mixed bucket clears the bootstrap CI in **1** of the six dimensions — and that bucket is not a clean post month. The pooled epoch contrast clears it in **1** of six. The axis is read as the record's own co-variation over two disjoint windows, **never** as damage: the reconstruction, the season, the acquisition geometry and the pre/post mix move with the collapse date as well.

```bash
python3 code/kdlo/kdlo_ocv_channels_csv.py       # build + verify the table
python3 code/kdlo/fig_kdlo_ocv_months.py          # this figure
```
