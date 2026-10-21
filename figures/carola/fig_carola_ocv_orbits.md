# Carolabrücke deck — OCV trajectories before and after the collapse

Trajectories of the six OCV dimensions across the collapse, one panel per **dimension x observation channel**, drawn from the support-sized chip: `gamma2, P, D, F, A, S` down the rows (the order the LUMO and Bautzen figures use), ASC / DESC across the columns. `P, D, F` and `A` are logarithmic, `gamma2` and `S` linear.

## States

Two states, as the record defines them: **pre-collapse** and **post-collapse** (table labels `pre-collapse (healthy)` / `post-collapse`).

## Rows and refusals

| channel | pre | post | runs |
| --- | --- | --- | --- |
| ASC | 124 | 166 | 2 |
| DESC | 60 | 86 | 2 |

Table `data/carola/support/carola_ocv_channels.csv` (436 overpasses, sha256 `3137d0b31fef`), cut manifest `data/carola/support/carola_rects_support.txt.manifest.json`: chip [69, 125], band [39, 15], support 585 cells, ring 7140 cells, guard 30.0 m, ring 200.0 m, **0 of 436 dates refused**.

## Per-state medians

**ASC**

| dimension | pre-collapse | post-collapse |
| --- | --- | --- |
| `gamma2` | 0.003712 (n=124) | 0.00361 (n=166) |
| `P` | 0.1031 (n=124) | 0.09496 (n=166) |
| `D` | 0.02413 (n=124) | 0.0261 (n=166) |
| `F` | 70 (n=124) | 77 (n=166) |
| `A` | 199 (n=124) | 219 (n=166) |
| `S` | 18.59 (n=124) | 17.27 (n=166) |

**DESC**

| dimension | pre-collapse | post-collapse |
| --- | --- | --- |
| `gamma2` | 0.0008537 (n=60) | 0.004132 (n=86) |
| `P` | 0.1138 (n=60) | 0.0919 (n=86) |
| `D` | 0.03606 (n=60) | 0.01141 (n=86) |
| `F` | 75.5 (n=60) | 23.5 (n=86) |
| `A` | 308.5 (n=60) | 70 (n=86) |
| `S` | 10.93 (n=60) | 20.22 (n=86) |

## What this does not claim

- The panels show the **support-sized chip**, not the committed table: the chip is larger and the reference is a disjoint ring, so these values are not the committed ones.
- A band with few overpasses is a handful of points, not a distribution; the median carries its n.
- The state boundary is drawn between runs, not at an event date: the record labels the runs, and the band follows the labels.
