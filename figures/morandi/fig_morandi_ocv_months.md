# Ponte Morandi (Polcevera), Genova — the echo-mask dimensions over the 24 months before the collapse

**Site** Ponte Morandi (Polcevera), Genova · **event** 2018-08-14 · **chips** 232 in the one-month windows M-24 … M-1 plus M+1

```
x     = [gamma2, P, D, A, F, S]      echo-mask vector
mask  = intensity >= 0.30*peak AND >= 5.0*np.median(intensity), A >= 2
M-24  [2016-08-14, 2016-09-14)   pre-collapse   first record month
M-1   [2018-07-14, 2018-08-14)   pre-collapse   last
M+1   [2018-08-14, 2018-09-14)   post-collapse
```

## Windows

| month | interval | state | n | days |
| --- | --- | --- | --- | --- |
| M-24 | [2016-08-14, 2016-09-14) | pre | 4 | 2016-08-24 … 2016-09-13 |
| M-23 | [2016-09-14, 2016-10-14) | pre | 7 | 2016-09-17 … 2016-10-13 |
| M-22 | [2016-10-14, 2016-11-14) | pre | 9 | 2016-10-17 … 2016-11-12 |
| M-21 | [2016-11-14, 2016-12-14) | pre | 10 | 2016-11-16 … 2016-12-12 |
| M-20 | [2016-12-14, 2017-01-14) | pre | 10 | 2016-12-16 … 2017-01-11 |
| M-19 | [2017-01-14, 2017-02-14) | pre | 8 | 2017-01-15 … 2017-02-10 |
| M-18 | [2017-02-14, 2017-03-14) | pre | 8 | 2017-02-16 … 2017-03-12 |
| M-17 | [2017-03-14, 2017-04-14) | pre | 9 | 2017-03-16 … 2017-04-11 |
| M-16 | [2017-04-14, 2017-05-14) | pre | 10 | 2017-04-15 … 2017-05-11 |
| M-15 | [2017-05-14, 2017-06-14) | pre | 9 | 2017-05-15 … 2017-06-10 |
| M-14 | [2017-06-14, 2017-07-14) | pre | 9 | 2017-06-14 … 2017-07-10 |
| M-13 | [2017-07-14, 2017-08-14) | pre | 11 | 2017-07-14 … 2017-08-13 |
| M-12 | [2017-08-14, 2017-09-14) | pre | 10 | 2017-08-15 … 2017-09-12 |
| M-11 | [2017-09-14, 2017-10-14) | pre | 10 | 2017-09-14 … 2017-10-12 |
| M-10 | [2017-10-14, 2017-11-14) | pre | 11 | 2017-10-14 … 2017-11-13 |
| M-9 | [2017-11-14, 2017-12-14) | pre | 10 | 2017-11-17 … 2017-12-13 |
| M-8 | [2017-12-14, 2018-01-14) | pre | 10 | 2017-12-17 … 2018-01-12 |
| M-7 | [2018-01-14, 2018-02-14) | pre | 10 | 2018-01-16 … 2018-02-11 |
| M-6 | [2018-02-14, 2018-03-14) | pre | 10 | 2018-02-15 … 2018-03-13 |
| M-5 | [2018-03-14, 2018-04-14) | pre | 9 | 2018-03-17 … 2018-04-12 |
| M-4 | [2018-04-14, 2018-05-14) | pre | 9 | 2018-04-16 … 2018-05-12 |
| M-3 | [2018-05-14, 2018-06-14) | pre | 10 | 2018-05-16 … 2018-06-11 |
| M-2 | [2018-06-14, 2018-07-14) | pre | 8 | 2018-06-17 … 2018-07-11 |
| M-1 | [2018-07-14, 2018-08-14) | pre | 10 | 2018-07-15 … 2018-08-10 |
| M+1 | [2018-08-14, 2018-09-14) | post | 11 | 2018-08-14 … 2018-09-13 |

Rolling one-month buckets anchored at the collapse date: the pre side is the 24 M-* months and the post side is M+1 alone. This is a companion of `fig_morandi_ocv_paper.py` (figures A–E); it adds no number to the A–E pin contract.

## Median of each dimension per month

| dimension | M-24 | M-23 | M-22 | M-21 | M-20 | M-19 | M-18 | M-17 | M-16 | M-15 | M-14 | M-13 | M-12 | M-11 | M-10 | M-9 | M-8 | M-7 | M-6 | M-5 | M-4 | M-3 | M-2 | M-1 | M+1 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| gamma2 | 0.07128 | 0.09801 | 0.03207 | 0.04993 | 0.06899 | 0.08584 | 0.06269 | 0.0554 | 0.09024 | 0.05625 | 0.08629 | 0.03201 | 0.05023 | 0.0924 | 0.1398 | 0.1772 | 0.1574 | 0.1448 | 0.1175 | 0.09103 | 0.7895 | 0.08577 | 0.1082 | 0.258 | 0.1016 |
| P | 0.07518 | 0.06342 | 0.03843 | 0.06144 | 0.07685 | 0.08621 | 0.08985 | 0.1275 | 0.08804 | 0.09914 | 0.08621 | 0.07266 | 0.08409 | 0.1107 | 0.1169 | 0.06304 | 0.0625 | 0.07694 | 0.07042 | 0.05603 | 0.08966 | 0.05981 | 0.07085 | 0.04918 | 0.07256 |
| D | 0.004773 | 0.00644 | 0.02579 | 0.02532 | 0.02247 | 0.02282 | 0.01259 | 0.005335 | 0.005073 | 0.004755 | 0.0053 | 0.005757 | 0.006104 | 0.005642 | 0.005725 | 0.007352 | 0.006929 | 0.01003 | 0.01688 | 0.003791 | 0.5833 | 0.004879 | 0.003505 | 0.008443 | 0.005159 |
| A | 18.5 | 9 | 12 | 10.5 | 11 | 11.5 | 12.5 | 14 | 15.5 | 14 | 19 | 19 | 22.5 | 11.5 | 7 | 6.5 | 9 | 7.5 | 8 | 10 | 5 | 11.5 | 9.5 | 8.5 | 13 |
| F | 6.5 | 3 | 4 | 2.5 | 4 | 3.5 | 4 | 5 | 7 | 7 | 8 | 9 | 8 | 5.5 | 4 | 3 | 3 | 3 | 3.5 | 4 | 1 | 5 | 5 | 4 | 7 |
| S | 13.71 | 6.184 | 5.077 | 4.376 | 5.098 | 5.973 | 6.169 | 8.157 | 13.89 | 25.05 | 11.07 | 20.11 | 9.106 | 7.61 | 4.693 | 14.57 | 21.46 | 12.7 | 12.18 | 11.47 | 1.088 | 7.767 | 5.466 | 14.47 | 13.78 |

`baseline` = the M-24 monthly medians (the 1.0 line of panel f): `gamma2` 0.07128, `P` 0.07518, `D` 0.004773, `A` 18.5, `F` 6.5, `S` 13.71.

## The 6D fingerprint at the cut (month / baseline)

| month | gamma2 | P | D | A | F | S |
| --- | --- | --- | --- | --- | --- | --- |
| M-24 | 1 | 1 | 1 | 1 | 1 | 1 |
| M-12 | 0.7046 | 1.119 | 1.279 | 1.216 | 1.231 | 0.664 |
| M-1 | 3.619 | 0.6542 | 1.769 | 0.4595 | 0.6154 | 1.055 |
| M+1 | 1.425 | 0.9651 | 1.081 | 0.7027 | 1.077 | 1.005 |

## M+1 against the pre side

| dimension | δ vs M-1 | CI≠0 vs M-1 | δ vs pooled pre | CI≠0 vs pooled pre |
| --- | --- | --- | --- | --- |
| gamma2 | -0.1091 | no | 0.2184 | no |
| P | -0.09091 | no | -0.1304 | no |
| D | -0.1182 | no | 0.01399 | no |
| A | 0.1455 | no | 0.04237 | no |
| F | 0.07273 | no | 0.08721 | no |
| S | -0.07273 | no | 0.0436 | no |

Cliff's delta > 0 means M+1 larger than M-1.

## What the figure shows (descriptive — not a claim of damage)

The 24 pre-collapse months carry **221** chips (median 10 per month) and the first post month **11**, so its intervals are wide and this is a descriptive companion, not a test.

Relative to the record-start baseline `M-24` (= 1.0) the **last** pre month M-1 sits at `gamma2` 3.619, `P` 0.6542, `D` 1.769, `A` 0.4595, `F` 0.6154, `S` 1.055; the first post month M+1 sits at `gamma2` 1.425, `P` 0.9651, `D` 1.081, `A` 0.7027, `F` 1.077, `S` 1.005. The pre-collapse record is **not** a flat line: season, weather, traffic and pipeline generation move the six dimensions month to month well before the event. Against that spread **0** of the six dimensions show a Cliff's-delta 95 % CI that excludes zero against M-1. The (P, D) trajectory (panel e) and the 6D fingerprint (panel f) show *when* the mask moved over the whole 24-month axis; the full per-month values are in the JSON (`blocks`, `fingerprint`, `deltas`).

Both sides of the event live on the *same* bridge, so season, weather, traffic and pipeline generation move with the collapse date too; the panels are read as that co-variation, **never** as damage.

```bash
python3 code/morandi/morandi_ocv_channels_csv.py    # build + verify the table
python3 code/morandi/fig_morandi_ocv_months.py     # this figure + report
```
