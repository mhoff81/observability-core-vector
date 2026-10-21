# Haltern AV9 wind turbine — OCV trajectories before and after the collapse

Trajectories of the six OCV dimensions across the collapse, one panel per **dimension x observation channel**, drawn from the support-sized chip: `gamma2, P, D, F, A, S` down the rows (the order the LUMO and Bautzen figures use), ASC / DESC across the columns. `P, D, F` and `A` are logarithmic, `gamma2` and `S` linear.

The thin dark line in each panel is a **rolling median over 5 consecutive overpasses** (centred, at least 3 valid values). It summarises local evolution — a dip in it is a summary, not a state, and whether it is a drift or a level shift is a question this figure does not answer.

## States

Two states, as the record defines them: **pre-collapse** and **post-collapse** (table labels `pre-collapse (healthy)` / `post-collapse`).

## Rows and refusals

| channel | pre | post | runs |
| --- | --- | --- | --- |
| ASC | 44 | 5 | 2 |
| DESC | 52 | 6 | 2 |

Table `data/haltern/support/haltern_ocv_channels.csv` (107 overpasses, sha256 `cb39c7473b29`), cut manifest `data/haltern/support/haltern_rects_support.txt.manifest.json`: chip [35, 229], band [5, 109], support 545 cells, ring 6618 cells, guard 30.0 m, ring 200.0 m, **13 of 120 dates refused**.

## Per-state medians

**ASC**

| dimension | pre-collapse | post-collapse |
| --- | --- | --- |
| `gamma2` | 0.004418 (n=44) | 0.00244 (n=5) |
| `P` | 0.06676 (n=44) | 0.0671 (n=5) |
| `D` | 0.04519 (n=44) | 0.0438 (n=5) |
| `F` | 186 (n=44) | 204 (n=5) |
| `A` | 359 (n=44) | 348 (n=5) |
| `S` | 21.51 (n=44) | 27.82 (n=5) |

**DESC**

| dimension | pre-collapse | post-collapse |
| --- | --- | --- |
| `gamma2` | 0.003275 (n=52) | 0.005925 (n=6) |
| `P` | 0.07472 (n=52) | 0.07407 (n=6) |
| `D` | 0.06139 (n=52) | 0.09126 (n=6) |
| `F` | 198.5 (n=52) | 223.5 (n=6) |
| `A` | 455.5 (n=52) | 686 (n=6) |
| `S` | 17.69 (n=52) | 15.08 (n=6) |

## What this does not claim

- The panels show the **support-sized chip**, not the committed table: the chip is larger and the reference is a disjoint ring, so these values are not the committed ones.
- A band with few overpasses is a handful of points, not a distribution; the median carries its n.
- The rolling median is a summary over five overpasses, not a state and not a trend test; it is drawn so that the local evolution is visible beside the global pre/post difference.
- The state boundary is drawn between runs, not at an event date: the record labels the runs, and the band follows the labels.
