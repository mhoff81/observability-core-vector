# Yeongdeok Unit 21 — OCV trajectories before and after the collapse

Trajectories of the six OCV dimensions before and after the collapse, one panel per **dimension x observation channel**, drawn from the support-sized chip: `gamma2, P, D, F, A, S` down the rows (the order the LUMO and Bautzen figures use), `ASC` / `DESC` across the columns. `P, D, F` and `A` are logarithmic, `gamma2` and `S` linear.

## States

Two states, as the record defines them: **pre-collapse** and **post-collapse** (table labels `pre-collapse (healthy)` / `post-collapse`).

## Rows and refusals

| channel | pre-collapse | post-collapse | runs |
| --- | --- | --- | --- |
| ASC | 1 | 0 | 1 |
| DESC | 11 | 2 | 2 |

Table `data/ywf/support/ywf_ocv_channels.csv` (14 overpasses, sha256 `ef6690e06267`), cut manifest `data/ywf/support/ywf_rects_support.txt.manifest.json`: chip [35, 189], band [5, 69], support 345 cells, ring 5658 cells, guard 30.0 m, ring 200.0 m, **19 of 33 dates refused**.

## Per-state medians

**ASC**

| dimension | pre-collapse | post-collapse |
| --- | --- | --- |
| `gamma2` | 0.0007654 (n=1) | – (n=0) |
| `P` | 0.1746 (n=1) | – (n=0) |
| `D` | 0.08889 (n=1) | – (n=0) |
| `F` | 127 (n=1) | – (n=0) |
| `A` | 588 (n=1) | – (n=0) |
| `S` | 7.79 (n=1) | – (n=0) |

**DESC**

| dimension | pre-collapse | post-collapse |
| --- | --- | --- |
| `gamma2` | 0.003112 (n=11) | 0.01003 (n=2) |
| `P` | 0.3496 (n=11) | 0.3535 (n=2) |
| `D` | 0.09675 (n=11) | 0.1014 (n=2) |
| `F` | 146 (n=11) | 146.5 (n=2) |
| `A` | 640 (n=11) | 670.5 (n=2) |
| `S` | 18.67 (n=11) | 25.75 (n=2) |

## What this does not claim

- The panels show the **support-sized chip**, not the committed table: the chip is larger and the reference is a disjoint ring, so these values are not the committed ones.
- A band with few overpasses is a handful of points, not a distribution; the median carries its n.
- The state boundary is drawn between runs, not at an event date: the record labels the runs, and the band follows the labels.
