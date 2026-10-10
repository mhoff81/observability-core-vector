# Carola Bruecke, Dresden — the echo-mask dimensions over the 24 months before the collapse

**Site** Carola Bruecke, Dresden · **event** 2024-09-11 · **chips** 952 (windows M-24 … M-1 plus M+1)

```
x     = [gamma2, P, D, A, F, S]      echo-mask vector
mask  = intensity >= 0.30*peak AND >= 5.0*np.median(intensity), A >= 2
M-24  [2022-09-11, 2022-10-11)       pre-collapse (healthy)  first record month
M-1   [2024-08-11, 2024-09-11)       pre-collapse (healthy)  last
M+1   [2024-09-11, 2024-10-11)       post-collapse
```

## The 24 pre-collapse months and the month after

Rolling one-month buckets anchored at the collapse date, from the record start (M-24) up to the event (M-1); the event itself opens M+1, so the pre side is the 24 M-* windows and the post side is M+1 alone. This is a companion of `fig_carola_ocv_paper.py` (figures A–E): it adds one figure and does not touch the A–E pin contract.

| month | range | state | chips | first day | last day | days |
| --- | --- | --- | --- | --- | --- | --- |
| M-24 | [2022-09-11, 2022-10-11) | pre | 16 | 2022-09-12 | 2022-10-09 | 8 |
| M-23 | [2022-10-11, 2022-11-11) | pre | 14 | 2022-10-14 | 2022-11-07 | 7 |
| M-22 | [2022-11-11, 2022-12-11) | pre | 16 | 2022-11-11 | 2022-12-08 | 8 |
| M-21 | [2022-12-11, 2023-01-11) | pre | 25 | 2022-12-13 | 2023-01-10 | 8 |
| M-20 | [2023-01-11, 2023-02-11) | pre | 35 | 2023-01-13 | 2023-02-06 | 7 |
| M-19 | [2023-02-11, 2023-03-11) | pre | 35 | 2023-02-11 | 2023-03-07 | 7 |
| M-18 | [2023-03-11, 2023-04-11) | pre | 40 | 2023-03-11 | 2023-04-07 | 8 |
| M-17 | [2023-04-11, 2023-05-11) | pre | 40 | 2023-04-12 | 2023-05-10 | 8 |
| M-16 | [2023-05-11, 2023-06-11) | pre | 35 | 2023-05-13 | 2023-06-06 | 7 |
| M-15 | [2023-06-11, 2023-07-11) | pre | 40 | 2023-06-11 | 2023-07-09 | 8 |
| M-14 | [2023-07-11, 2023-08-11) | pre | 40 | 2023-07-12 | 2023-08-10 | 8 |
| M-13 | [2023-08-11, 2023-09-11) | pre | 40 | 2023-08-14 | 2023-09-10 | 8 |
| M-12 | [2023-09-11, 2023-10-11) | pre | 35 | 2023-09-15 | 2023-10-09 | 7 |
| M-11 | [2023-10-11, 2023-11-11) | pre | 40 | 2023-10-13 | 2023-11-09 | 8 |
| M-10 | [2023-11-11, 2023-12-11) | pre | 35 | 2023-11-14 | 2023-12-08 | 7 |
| M-9 | [2023-12-11, 2024-01-11) | pre | 40 | 2023-12-12 | 2024-01-08 | 8 |
| M-8 | [2024-01-11, 2024-02-11) | pre | 40 | 2024-01-13 | 2024-02-10 | 8 |
| M-7 | [2024-02-11, 2024-03-11) | pre | 38 | 2024-02-13 | 2024-03-08 | 7 |
| M-6 | [2024-03-11, 2024-04-11) | pre | 48 | 2024-03-13 | 2024-04-10 | 8 |
| M-5 | [2024-04-11, 2024-05-11) | pre | 42 | 2024-04-13 | 2024-05-07 | 7 |
| M-4 | [2024-05-11, 2024-06-11) | pre | 48 | 2024-05-12 | 2024-06-09 | 8 |
| M-3 | [2024-06-11, 2024-07-11) | pre | 42 | 2024-06-12 | 2024-07-06 | 7 |
| M-2 | [2024-07-11, 2024-08-11) | pre | 52 | 2024-07-11 | 2024-08-08 | 8 |
| M-1 | [2024-08-11, 2024-09-11) | pre | 56 | 2024-08-11 | 2024-09-09 | 7 |
| M+1 | [2024-09-11, 2024-10-11) | post | 60 | 2024-09-13 | 2024-10-10 | 8 |

## Median of each dimension per month

| month | n | gamma2 | P | D | A | F | S |
| --- | --- | --- | --- | --- | --- | --- | --- |
| M-24 | 16 | 0.06874 | 0.03428 | 0.005454 | 18 | 9 | 19.3 |
| M-23 | 14 | 0.1908 | 0.0381 | 0.004589 | 9 | 5 | 21.97 |
| M-22 | 16 | 0.1073 | 0.03809 | 0.005667 | 17 | 8 | 18.27 |
| M-21 | 25 | 0.03728 | 0.05057 | 0.006944 | 12 | 4 | 18.72 |
| M-20 | 35 | 0.1611 | 0.03626 | 0.01292 | 7 | 3 | 15.6 |
| M-19 | 35 | 0.05066 | 0.03243 | 0.004872 | 24 | 11 | 12.34 |
| M-18 | 40 | 0.05284 | 0.04633 | 0.007777 | 9.5 | 4 | 12.26 |
| M-17 | 40 | 0.08806 | 0.03592 | 0.006252 | 12 | 5 | 17.03 |
| M-16 | 35 | 0.09045 | 0.03998 | 0.006259 | 15 | 6 | 19.75 |
| M-15 | 40 | 0.01327 | 0.03209 | 0.00603 | 18.5 | 7.5 | 17.84 |
| M-14 | 40 | 0.02251 | 0.04059 | 0.01167 | 20.5 | 8 | 9.751 |
| M-13 | 40 | 0.06156 | 0.02766 | 0.005985 | 13.5 | 6.5 | 26.3 |
| M-12 | 35 | 0.04962 | 0.03439 | 0.00565 | 22 | 10 | 16.59 |
| M-11 | 40 | 0.1121 | 0.02425 | 0.005396 | 18 | 8.5 | 23.3 |
| M-10 | 35 | 0.1806 | 0.03663 | 0.007822 | 22 | 9 | 18.88 |
| M-9 | 40 | 0.2975 | 0.04191 | 0.008024 | 8 | 2 | 13.05 |
| M-8 | 40 | 0.2506 | 0.03304 | 0.009723 | 12 | 2.5 | 17.32 |
| M-7 | 38 | 0.1546 | 0.03955 | 0.007532 | 13 | 6 | 23.16 |
| M-6 | 48 | 0.1519 | 0.03417 | 0.006222 | 12.5 | 4 | 23.42 |
| M-5 | 42 | 0.239 | 0.04565 | 0.005634 | 11 | 2 | 18.89 |
| M-4 | 48 | 0.197 | 0.03986 | 0.00642 | 11 | 3 | 18.08 |
| M-3 | 42 | 0.1578 | 0.02199 | 0.004717 | 13 | 7 | 25 |
| M-2 | 52 | 0.07215 | 0.02563 | 0.005169 | 13 | 7 | 15.11 |
| M-1 | 56 | 0.0365 | 0.0327 | 0.007286 | 20 | 8 | 19.09 |
| M+1 | 60 | 0.06237 | 0.03182 | 0.005789 | 17 | 9 | 10.65 |

`baseline` = the M-24 monthly medians (the 1.0 line of panel f): `gamma2` 0.06874, `P` 0.03428, `D` 0.005454, `A` 18, `F` 9, `S` 19.3.

## The 6D fingerprint at the cut (month / baseline)

| month | gamma2 | P | D | A | F | S |
| --- | --- | --- | --- | --- | --- | --- |
| M-24 | 1 | 1 | 1 | 1 | 1 | 1 |
| M-12 | 0.7219 | 1.003 | 1.036 | 1.222 | 1.111 | 0.8598 |
| M-1 | 0.531 | 0.9537 | 1.336 | 1.111 | 0.8889 | 0.9895 |
| M+1 | 0.9073 | 0.9281 | 1.061 | 0.9444 | 1 | 0.5521 |

## What moved across the cut (M+1 vs M-1)

| dim | median M-1 | median M+1 | Cliff's delta | SDS | CI95 excludes 0 | CI95 |
| --- | --- | --- | --- | --- | --- | --- |
| gamma2 | 0.0365 | 0.06237 | 0.07798 | 0.7798 | no | [-0.15, 0.3036] |
| P | 0.0327 | 0.03182 | -0.1917 | 1.917 | no | [-0.4024, 0.02292] |
| D | 0.007286 | 0.005789 | -0.3554 | 3.554 | yes | [-0.5503, -0.1571] |
| A | 20 | 17 | 0.1821 | 1.821 | no | [-0.02827, 0.3839] |
| F | 8 | 9 | 0.1664 | 1.664 | no | [-0.04673, 0.3723] |
| S | 19.09 | 10.65 | 0.01786 | 0.1786 | no | [-0.2021, 0.2327] |

Cliff's delta > 0 means M+1 larger than M-1. The same comparison against all 24 pre months pooled is in the JSON (`deltas[dim].vs_pre_pool`).

## Reading (descriptive — not a claim of damage)

The 24 pre-collapse months carry **892** chips (median 40 per month) and the first post month **60**, so the monthly intervals are wide and this is a descriptive companion, not a test.

Relative to the record-start baseline `M-24` (= 1.0) the **last** pre month M-1 sits at `gamma2` 0.531, `P` 0.9537, `D` 1.336, `A` 1.111, `F` 0.8889, `S` 0.9895; the first post month M+1 sits at `gamma2` 0.9073, `P` 0.9281, `D` 1.061, `A` 0.9444, `F` 1, `S` 0.5521. The pre-collapse record is therefore **not** a flat line — season, weather, traffic and pipeline generation move the six dimensions month to month well before the event. Against that spread the single adjacent contrast (M+1 vs M-1) clears the bootstrap CI in **1** of the six dimensions; the full per-month trajectory is in the JSON (`blocks`, `fingerprint`, `deltas`).

The month axis therefore does **not** reproduce the pooled A–E headline by itself; it shows how much of the pre-collapse record already moves on its own and where the collapse sits inside it. Both sides live on the *same* bridge, so season, weather, traffic and pipeline generation move with the collapse date as well — the panels are read as that co-variation, **never** as damage.

```bash
python3 code/carola/carola_ocv_channels_csv.py       # build + verify the table
python3 code/carola/fig_carola_ocv_months.py          # this figure
```
