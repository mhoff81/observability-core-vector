# Ponte Morandi (Polcevera), Genova — the observability core vector around the collapse

**Site** Ponte Morandi (Polcevera), Genova · **event** 2018-08-14 · **chips** 232 (221 pre-collapse (healthy) / 11 post-collapse)

```
state = pre-collapse (healthy)  day <  2018-08-14
        post-collapse            day >= 2018-08-14
x     = [gamma2, P, D, A, F, S]   echo-mask vector
OCV   = [gamma2, P, D]            observability core vector
mask  = intensity >= 0.30*peak AND >= 5.0*np.median(intensity), A >= 2
```

## Honest null — read this first

Only **11** chips survive the collapse (5 `A_asc` + 6 `A_des`) against **221** before it, so every post-state number carries a wide interval. The analysed 80x80 window is **not on the collapsed span** and the archive holds no persistent scatterer at this footprint. This figure therefore asks whether the states are *separable* here — it must never be read as a damage verdict.

## Effect sizes of the six channels (post − pre)

| channel | Cliff's delta | 95 % CI | CI≠0 | SDS |
| --- | --- | --- | --- | --- |
| gamma2 | 0.2184 | [-0.1641, 0.5821] | no | 2.184 |
| P | -0.1304 | [-0.4825, 0.2209] | no | 1.304 |
| D | 0.01399 | [-0.3344, 0.3933] | no | 0.1399 |
| A | 0.04237 | [-0.3365, 0.3982] | no | 0.4237 |
| F | 0.08721 | [-0.3422, 0.4862] | no | 0.8721 |
| S | 0.0436 | [-0.3579, 0.4077] | no | 0.436 |

## Per-track deck-mask effect sizes

**A_asc**

| channel | Cliff's delta | 95 % CI | CI≠0 |
| --- | --- | --- | --- |
| gamma2 | -0.1018 | [-0.5317, 0.2632] | no |
| P | -0.2421 | [-0.7088, 0.214] | no |
| D | -0.5754 | [-0.786, -0.328] | yes |
| A | 0.3807 | [0.2026, 0.5606] | yes |
| F | 0.6088 | [0.4421, 0.7632] | yes |
| S | 0.7123 | [0.456, 0.8982] | yes |

**A_des**

| channel | Cliff's delta | 95 % CI | CI≠0 |
| --- | --- | --- | --- |
| gamma2 | 0.4081 | [-0.142, 0.849] | no |
| P | 0.04673 | [-0.4191, 0.5327] | no |
| D | 0.3146 | [-0.222, 0.8092] | no |
| A | -0.1355 | [-0.6194, 0.3958] | no |
| F | -0.2056 | [-0.7976, 0.3825] | no |
| S | -0.2913 | [-0.8039, 0.2274] | no |

## Is the state learnable? (2-class LDA, leave-one-out)

| model | dims | LOO acc. | balanced acc. | p (perm.) |
| --- | --- | --- | --- | --- |
| gamma2 | 1 | 0.793 | 0.589 | 0.063 |
| PD | 2 | 0.772 | 0.535 | 0.108 |
| gamma2PD | 3 | 0.772 | 0.578 | 0.096 |
| x_6d | 6 | 0.746 | 0.521 | 0.077 |
| weather | 4 | 0.685 | 0.662 | 0.134 |

## What this does and does not say

Both sides of the event live on the *same* bridge, so season, weather and pipeline generation move with the state — every post-state interval is wide and no contrast here is a damage verdict. The point of the figure is the **honest null**: with post = 11 chips the states are not cleanly separable in the analysed window, and the `A` column of the site's export is a *fixed 640-px quantile mask*, not the echo mask this package recomputes — the two never coincide (panel E).

```bash
python3 code/morandi/morandi_ocv_channels_csv.py   # build + verify the table
python3 code/morandi/fig_morandi_ocv_paper.py      # figures A–E + this report
```
