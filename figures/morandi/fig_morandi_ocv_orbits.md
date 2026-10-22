# Ponte Morandi deck — OCV trajectories before and after the collapse

Trajectories of the six OCV dimensions before and after the collapse, one panel per **dimension x observation channel**, drawn from the support-sized chip: `gamma2, P, D, F, A, S` down the rows (the order the LUMO and Bautzen figures use), `ASC` / `DESC` across the columns. `P, D, F` and `A` are logarithmic, `gamma2` and `S` linear.

## States

Two states, as the record defines them: **pre-collapse** and **post-collapse** (table labels `pre-collapse (healthy)` / `post-collapse`).

## Rows and refusals

| channel | pre-collapse | post-collapse | runs |
| --- | --- | --- | --- |
| ASC | 91 | 3 | 2 |
| DESC | 107 | 6 | 2 |

Table `data/morandi/support/morandi_ocv_channels.csv` (207 overpasses, sha256 `249c9ee6630c`), cut manifest `data/morandi/support/morandi_rects_support.txt.manifest.json`: chip [51, 137], band [21, 41], support 861 cells, ring 5448 cells, guard 30.0 m, ring 200.0 m, **25 of 232 dates refused**.

## Per-state medians

**ASC**

| dimension | pre-collapse | post-collapse |
| --- | --- | --- |
| `gamma2` | 0.09039 (n=91) | 0.05069 (n=3) |
| `P` | 0.1682 (n=91) | 0.1752 (n=3) |
| `D` | 0.1143 (n=91) | 0.004079 (n=3) |
| `F` | 4 (n=91) | 6 (n=3) |
| `A` | 11 (n=91) | 12 (n=3) |
| `S` | 2.637 (n=91) | 2.091 (n=3) |

**DESC**

| dimension | pre-collapse | post-collapse |
| --- | --- | --- |
| `gamma2` | 0.002662 (n=107) | 0.0112 (n=6) |
| `P` | 0.06533 (n=107) | 0.0825 (n=6) |
| `D` | 0.02169 (n=107) | 0.01286 (n=6) |
| `F` | 54 (n=107) | 36.5 (n=6) |
| `A` | 142 (n=107) | 83.5 (n=6) |
| `S` | 10.37 (n=107) | 19.72 (n=6) |

## What this does not claim

- The panels show the **support-sized chip**, not the committed table: the chip is larger and the reference is a disjoint ring, so these values are not the committed ones.
- A band with few overpasses is a handful of points, not a distribution; the median carries its n.
- The state boundary is drawn between runs, not at an event date: the record labels the runs, and the band follows the labels.
