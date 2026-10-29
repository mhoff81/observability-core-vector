# `concept/fig_ocv_schematic` — the observability core vector, drawn

![the observability core vector](fig_ocv_schematic.jpg)

One plate, three panels, no data. Every other figure in `figures/` is computed from a record and pinned to it; this one is computed from the **definition** the records are measured with, and its job is to say what the core vector is before any site's numbers are read. It therefore has no pins to a table — it has structural checks instead (88 of them, listed below and in `data/concept/fig_ocv_schematic.json`).

```
physical asset -> SAR observation geometry -> asset-specific measurement
support -> radar echo representation -> OCV
```

## The three panels

**(a) the asset and the acquisition.** A generic lattice asset `G_a^{3D}` in its own frame (the grey triad), observed from a satellite along the orange line of sight `LOS` at the look angle `θ` measured from nadir; the operator `Π_c` projects it into the radar frame `(r, x)`. The asset is drawn as four legs, five belts and alternating braces — something with a real cross-section, because that is what makes the projection lossy.

**(b) the measurement support.** The projected asset `Π_c(G_a^{3D})` sits inside the asset- and channel-specific support `W_{a,c}` (blue), which sits inside the read chip `C_{a,c}` (grey, outermost). Between them is the **background ring** — dashed, lightly tinted — from which the no-echo reference is read. The panel states the set algebra `W_{a,c} ⊂ C_{a,c}` once at its foot and nothing else: no box, no legend, no count.

**(c) the observable representation and the vector.** Inside the same `W_{a,c}` the echo mask is a raster of green cells: 4 connected regions of 8, 6, 4, 3 cells in a 16 × 10 grid (density 0.131), one closed hole, and a dashed centroid-to-peak shift (2.63 cells) between `p̄` and `p⋆`. That is the extent, the density, the fragmentation and the spatial shift the vector has to summarise, made visible before it is made numeric. Below it:

```
OCV(t) = [ γ², P, D, F, A, S ]
           ⏟           ⏟
        cross-acquisition  spatial geometry
```

The two groups are carried by two restrained colours and two thin under-braces — not by boxes and not by a second legend. Each brace spans the measured extent of exactly the tokens it names: the row is laid out from the renderer's own text metrics at the output dpi, and the inter-group gap is *solved* so that the two names clear each other.

| group | components | reads as | colour |
| --- | --- | --- | --- |
| cross-acquisition | `γ²`, `P` | what changes between acquisitions of the same place | blue `#1f4e79` |
| spatial geometry | `D`, `F`, `A`, `S` | how the echo is laid out inside `W_{a,c}` | green `#3f7d3f` |

## What the plate does not show

It is a schematic, so it carries no site, no date, no state and no number that could be read as a result. Deliberately absent:

- structural health, damage probability, failure detection or any other diagnosis — the figure defines a measurement, it does not deliver a verdict;
- displacement, deformation or InSAR geometry — the OCV is a coherence and echo-geometry vector, and nothing in the plate moves;
- `R_B` and `R_G` — those belong to the temporal reconfiguration analysis (`figures/sds/`), not to the definition of the vector;
- a score, a threshold, a weighting or a count.

The prohibition is mechanical, not editorial: a forbidden-token scan runs over every drawn string, and the whole text inventory is compared against a written-out whitelist, so the plate cannot grow a `damage` label without failing the run.

## The plates, and how to reproduce them

| file | what it is |
| --- | --- |
| `concept/fig_ocv_schematic.jpg` | 183 mm × 61 mm, 630 dpi JPEG (quality 95, 4:4:4 subsampling) |
| `concept/fig_ocv_schematic.pdf` | the same plate as a true vector page (FontType 42; no raster object anywhere in the file) |
| `concept/fig_ocv_schematic.md` | this report |
| `../data/concept/fig_ocv_schematic.json` | the drawn geometry, the text inventory, every check and the sha256 of both images |

```bash
bash figures/concept/fig_ocv_schematic.sh                  # both images
python3 code/concept/fig_ocv_schematic.py --no-figures     # checks only
PYTHON=/usr/bin/python3 bash figures/concept/fig_ocv_schematic.sh
```

## The guard block

**88 of 88 checks pass.** A failure exits with status 2 *after* the artefacts are written, so a broken plate is reported rather than silently delivered. (`--no-figures` builds, measures and checks, writes the JSON and this report, and skips both images and the eleven file checks that go with them.)

| family | what it pins |
| --- | --- |
| canvas | 183 mm × 61 mm at 630 dpi, ratio 3.000, `savefig.bbox` `standard` (a tight bbox would move the 183 mm edge); every panel isotropic, so one unit is one unit in x and in y and the raster cells are square |
| text inventory | exactly the declared 28 strings and no others; no title; no legend; all axes off |
| legibility | every text ≥ 7.0 pt at 183 mm and ≥ 6.5 pt at 170 mm; every text inside its own panel; no two texts overlap |
| vocabulary | no forbidden token in any drawn string (health, damage, fail, displace, deform, InSAR, score, threshold, px, `R_B`, `R_G`, …) |
| set algebra | `W_{a,c}` inside the ring inside `C_{a,c}` with a visible margin; the band between `W` and the ring is ≥ 0.03 units; the projected asset lies inside `W` |
| one W, two panels | the support is **one** polygon definition mapped into (b) and into (c); normalised to its own box it is identical to 1 × 10⁻¹² in both |
| raster integrity | cells square on the drawn canvas and ≥ 1.5 mm; every echo cell inside `W` and inside the raster; exactly 4 connected regions of 8, 6, 4 and 3 cells; exactly one enclosed hole; density in (5 %, 20 %); centroid and peak ≥ 1.5 cells apart |
| the vector | six components in the order `γ², P, D, F, A, S`; the partition into (γ², P) and (D, F, A, S) proved from the token table; each brace spans its own group's measured extent; the two names clear each other; the row sits inside the panel |
| determinism | two independent builds render identical pixels and lay the row out byte-identically; the PDF is written without a creation date, so both files are byte-identical across separate runs |
| the files | JPEG 4539 × 1513 px ± 2 at 630 dpi; the PDF has no `/Subtype /Image`, embeds `/FontFile2`, carries no Type-3 font, and its page box is 183 mm × 61 mm |

### Every check, in the order it runs

| # | check | result | detail |
| --- | --- | --- | --- |
| 1 | canvas width 183 mm | ok | 183.0000 mm |
| 2 | canvas height 61 mm | ok | 61.0000 mm |
| 3 | canvas ratio 3.000:1 | ok | 3.00000 |
| 4 | figure dpi is the output dpi | ok | 630 |
| 5 | panel a is isotropic (unit = unit) | ok | dx=1379.8488px dy=1379.8488px |
| 6 | panel a box matches the layout | ok | 0.304000x0.960000 |
| 7 | panel b is isotropic (unit = unit) | ok | dx=1379.8488px dy=1379.8488px |
| 8 | panel b box matches the layout | ok | 0.304000x0.960000 |
| 9 | panel c is isotropic (unit = unit) | ok | dx=1379.8488px dy=1379.8488px |
| 10 | panel c box matches the layout | ok | 0.304000x0.960000 |
| 11 | panel a/b boxes do not overlap | ok | gap 0.0320 |
| 12 | panel b/c boxes do not overlap | ok | gap 0.0320 |
| 13 | panel unit is 55.632 mm | ok | 55.632 mm |
| 14 | text inventory is exactly the declared set | ok | 28 objects |
| 15 | text count <= 30 | ok | 28 |
| 16 | no duplicate panel label | ok |  |
| 17 | no title anywhere | ok |  |
| 18 | no legend anywhere | ok |  |
| 19 | axes are off (no spines, no ticks) | ok |  |
| 20 | every text >= 7.0 pt | ok | min 7.00 pt |
| 21 | every text >= 7.0 pt at 170 mm (scale 0.929) = 6.5 pt | ok | 6.50 pt |
| 22 | no forbidden vocabulary in any drawn text | ok |  |
| 23 | R_B absent | ok |  |
| 24 | R_G absent | ok |  |
| 25 | the coherence channel is mathtext | ok |  |
| 26 | panel a: every text inside the panel | ok |  |
| 27 | panel a: no two texts overlap | ok |  |
| 28 | panel b: every text inside the panel | ok |  |
| 29 | panel b: no two texts overlap | ok |  |
| 30 | panel c: every text inside the panel | ok |  |
| 31 | panel c: no two texts overlap | ok |  |
| 32 | W has eight vertices in both panels | ok |  |
| 33 | W in (b) sits inside C with a >= 20 um margin | ok |  |
| 34 | W in (b) sits strictly inside the background ring | ok |  |
| 35 | the background ring is a band >= 0.03 wide | ok | 0.060 units |
| 36 | the ring sits inside C, so it cannot be the chip edge | ok |  |
| 37 | C sits inside the panel | ok |  |
| 38 | the projected asset lies inside W | ok |  |
| 39 | the projected asset is inside C, not straddling it | ok |  |
| 40 | the ring leader ends between W and the ring | ok | y=0.762 |
| 41 | W in (b) and W in (c) are the same polygon | ok |  |
| 42 | W in (b) and W in (c) are the same shape, not the same place | ok |  |
| 43 | raster cells are square on the drawn canvas | ok | 0.03625 units |
| 44 | raster cell >= 1.5 mm on the printed page | ok | 2.02 mm |
| 45 | the raster sits inside W in (c) | ok |  |
| 46 | the raster sits inside W with a visible margin | ok |  |
| 47 | every echo cell lies inside W in (c) | ok |  |
| 48 | every echo cell lies inside the raster | ok |  |
| 49 | the echo mask has exactly four connected regions | ok | 4 |
| 50 | the four regions are 8, 6, 4 and 3 cells | ok | [8, 6, 4, 3] |
| 51 | every region is 4-connected (a fragment, not a dust field) | ok |  |
| 52 | the mask encloses exactly one hole | ok | 1 |
| 53 | the hole is the centre of the 3x3 ring | ok | [(4, 4)] |
| 54 | echo density is sparse and legible (5 % < d < 20 %) | ok | 0.1313 |
| 55 | the peak lies on an echo cell | ok |  |
| 56 | the centroid and the peak are >= 1.5 cells apart | ok | 2.626 cells |
| 57 | the centroid is not inside an echo cell (the mask is fragmented) | ok |  |
| 58 | the vector has six components | ok |  |
| 59 | the order is gamma2, P, D, F, A, S | ok | gamma2,P,D,F,A,S |
| 60 | the partition covers every component exactly once | ok |  |
| 61 | group one is exactly gamma2 and P | ok | $\gamma^2,\ P$ |
| 62 | group two is exactly D, F, A and S | ok | $D,\ F,\ A,\ S$ |
| 63 | the row starts inside the panel | ok | 0.0339 |
| 64 | the row ends inside the panel | ok | 0.9661 |
| 65 | the row is centred in its panel | ok | 0.5000 |
| 66 | the cross brace spans its own group only | ok | text 0.312..0.441 brace 0.312..0.441 |
| 67 | the cross name sits under its brace | ok | top 0.279 < brace 0.285 |
| 68 | the spatial brace spans its own group only | ok | text 0.648..0.922 brace 0.648..0.922 |
| 69 | the spatial name sits under its brace | ok | top 0.279 < brace 0.285 |
| 70 | the separating comma sits between the two groups | ok | 0.453..0.471 |
| 71 | the separating comma does not touch either group | ok | 0.012 / 0.177 |
| 72 | the two group names do not overlap | ok | 0.562 < 0.597 |
| 73 | the two group names do not run into each other's brace | ok | gap 0.035 |
| 74 | the equation sits above the braces | ok |  |
| 75 | the equation does not touch the echo raster | ok |  |
| 76 | JPEG width is 183 mm at 630 dpi (4539 px +/- 2) | ok | 4538 px |
| 77 | JPEG height is 61 mm at 630 dpi (1513 px +/- 2) | ok | 1512 px |
| 78 | JPEG pixel ratio is 3.000 | ok | 3.00132 |
| 79 | JPEG carries 630 dpi | ok | (630, 630) |
| 80 | JPEG is a full-colour raster | ok | RGB |
| 81 | PDF page contains no raster image | ok |  |
| 82 | PDF embeds its fonts as TrueType (FontType 42) | ok |  |
| 83 | PDF has no Type-3 fonts | ok |  |
| 84 | PDF carries no timestamp (it must be reproducible) | ok |  |
| 85 | PDF page box is 183 x 61 mm | ok | [0.0, 0.0, 518.7401574803, 172.9133858268] |
| 86 | re-saving the same figure gives the same JPEG bytes | ok |  |
| 87 | two builds render identical pixels | ok | (1512, 4538, 4) |
| 88 | two builds lay the OCV row out identically | ok |  |

## Files

- `fig_ocv_schematic.jpg` — sha256 `8a06dfbd0e497040c4431c003f3ba9e0a67c73b3cdffb09c366e5fadf53301fc`
- `fig_ocv_schematic.pdf` — sha256 `e5053ea80833ae4a3338986537f5dcbcea5b995806afd6fe35654b2e52bd3162`

Environment: Python 3.13.5, matplotlib 3.11.1, numpy 2.4.6. The plate is deterministic: no timestamp, no seed-dependent value and no machine-dependent path reaches any output, so two runs produce byte-identical JSON, Markdown, JPEG and PDF.
