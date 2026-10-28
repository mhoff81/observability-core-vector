# Prospective read — the case and its control, side by side

| | morandi (case) | carola (control) |
|---|---|---|
| event | 2018-08-14 | 2024-09-11 |
| window source | the stated rule selected it: the earliest run of twelve consecutive calendar months in which every month holds at least four chips | the registry declares it, and its descriptive plate is published against it: the registry's window for this asset |
| trajectory rule | **CHANGE** — sustained multidimensional observability departure | **WATCH** — isolated unusual departure |
| pinned policy | **NOMINAL** | **NOMINAL** |
| policy, every firing rule | **INSUFFICIENT** | **INSUFFICIENT** |
| first month outside | 2017-09 (2017-09-15, 333 days) | 2024-01 (2024-01-15, 240 days) |
| first persistence | 2017-12 (2017-12-15, 242 days) | 2024-02 (2024-02-15, 209 days) |
| first CHANGE | 2017-12 (2017-12-15, 242 days) | did not happen |
| signed lead (days) | 333 | 240 |
| sustained lead (days) | 242 | — |

Both readers, both assets, the same layer. The readers disagree on both, and the disagreement is printed as two words rather than resolved into one — no word on this page is a reconciliation of the two readers, here or anywhere.

| stream | case run | case months outside | control run | control months outside | case median | control median |
|---|---|---|---|---|---|---|
| `R_B` | 6 | 7 / 12 | 1 | 2 / 9 | 1.11 | 0.81 |
| `R_G` | 2 | 4 / 12 | 1 | 2 / 9 | 1.66 | 0.63 |
| `gamma2` | 10 | 11 / 12 | 2 | 4 / 9 | -1.34 | 0.26 |
| `P` | 1 | 1 / 12 | 3 | 5 / 9 | 0.18 | -0.37 |
| `A` | 10 | 10 / 12 | 0 | 0 / 9 | 1.45 | 0.05 |
| `D` | 0 | 0 / 12 | 0 | 0 / 9 | 1.31 | 0.05 |
| `F` | 9 | 9 / 12 | 1 | 1 / 9 | 1.16 | 0.05 |
| `S` | 1 | 1 / 12 | 2 | 3 / 9 | 0.06 | 0.06 |

Streams whose longest run differs: `R_B`, `R_G`, `gamma2`, `P`, `A`, `F`, `S`.

Read in this order: the case, then the control, then the case again. The third reading is compared with the first and must be identical; the control is never an input to the case. The two assets are not pooled and no number here is an average of them.
