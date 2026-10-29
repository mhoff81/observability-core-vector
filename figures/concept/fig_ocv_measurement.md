# `concept/fig_ocv_measurement` — one chip, three stages, one number set

![one chip, three stages, one number set](fig_ocv_measurement.jpg)

The companion plate to `concept/fig_ocv_schematic`. That one defines the vector from the asset outwards; this one answers the question a reader asks next — **what exactly is measured, and where?** It takes one radar read chip and shows it three times, each time admitting less:

```
C_{a,c}  ->  W_{a,c}  ->  M_{a,c}(t)  ->  OCV(t)
native       measurement     echo mask      the six
read chip    support                        components
```

One idea, three panels, one frame: **restriction**. Panel (a) holds everything the sensor returned; (b) keeps the part that belongs to this asset in this channel; (c) keeps the part of that which carries echo. Nothing in (b) is absent from (a), and nothing in (c) is absent from (b) — the panels are one observation at three levels of admission, not three pictures.

Like the first schematic it carries no data, so it has no pins to a table; it has 160 structural checks instead, listed at the foot of this page and in `../data/concept/fig_ocv_measurement.json`.

## The three panels

**(a) the native read chip.** `C_{a,c}` as it comes off the acquisition: the intensity field over the whole chip, 52 × 39 = 2028 native samples at 0.81 mm each, drawn as 2028 filled vector quads on a 16-step grey ramp — not as an embedded image, so the PDF stays a true vector page. The asset footprint `Π_c(G_a^{3D})` is a thin outline in the middle of it (the symbol itself is set in (b), on the footprint), the chip's own boundary carries `C_{a,c}` on its lower edge, and the radar frame `az` / `r` sits at the chip's lower-left corner. Nothing is subdued: this is the read, with more cells in it than any measurement will use.

**(b) the measurement support.** The same field, restricted. Outside the background ring it is pulled 65 % of the way to white; the ring band itself is a flat tint (it is the no-echo reference, not structure); the guard between the ring and the support is pure white; and inside `W_{a,c}` the field returns at full contrast under a blue veil (18 % blue). The asset is dark again, as on the first schematic, but solid: what is drawn is the *projection* `Π_c(G_a^{3D})`, a silhouette, not the tower with its storeys, so it carries no floor lines. Its name is set in paper white on the silhouette itself: the first schematic sets the same symbol under its tower, and here there is no band between the footprint and the support's own edge deep enough to hold a 9 pt label clear of both. The restriction is therefore *visible*: the reader sees what is kept and what is dropped, in the same picture, at the same place, rather than being told about it in a caption.

**(c) the echo mask.** Inside `W_{a,c}` the mask: every support cell is drawn, the 130 without echo in very light grey and the 20 with echo in green. The footprint is repeated as a thin outline, so the echo can be seen to sit on the asset rather than beside it, and no echo cell crosses that outline. Outside the support there is no chip plate at all — bare paper, because the vector row runs through the band where the chip's lower edge would fall, and nothing outside the support is evaluated. Below the support the mask becomes the vector, `OCV(t) = [γ², P, D, F, A, S]`, in its two conceptual groups.

## One frame, drawn three times

All three panels read one frame, and the plate says so in a way a reader can check by eye: the chip rectangle, the asset footprint, the radar frame and the scale are identical in (a), (b) and (c) — the same 0.760 × 0.570 panel-unit rectangle (42.3 × 31.7 mm), drawn at the same device coordinates in all three panels, and the checks compare the drawn rectangles in pixels rather than trusting the constants.

The nesting is `W_{a,c} ⊂ guard ⊂ background ring ⊂ C_{a,c}`:

| element | what it is | drawn as |
| --- | --- | --- |
| `C_{a,c}` | the read chip | grey rectangle, 42.3 × 31.7 mm, 2028 native samples |
| background ring | the no-echo reference | flat grey band 2.2 mm wide, dashed boundary, `background ring` on a leader |
| guard | keeps the ring out of the support | pure white band, 1.1 mm |
| `W_{a,c}` | the asset- and channel-specific support | the same chamfered octagon as the first schematic, blue, 27.8 × 19.5 mm |

`W_{a,c}` is **one** polygon definition mapped once into the two panels that draw it, so the panels cannot drift apart; a check normalises both to their own box and compares them to the definition at 1 × 10⁻¹².

## The mask is derived, not decorated

The field is written down in closed form: a smooth base, two texture terms and seven declared responses (one annulus — the tower seen as a ring — three compact blobs on the asset, one diffuse glow that never reaches the reference level, and two clutter returns *outside* the support). There is no random number anywhere in the file, so the field is a constant and so is the mask.

The mask is the set of support cells whose **peak native sample exceeds one declared reference level** (0.490 on a field that runs 0.186 … 0.782). That definition produces, and the checks insist on:

| property | value |
| --- | --- |
| support lattice | 15 × 10 = 150 cells, one cell = exactly 2 × 2 native samples (1.6 mm) |
| echo cells | 20 (13.3 % of the support) |
| connected regions | 4, of 8, 6, 4, 2 cells |
| enclosed holes | 1 — the dark centre of the annulus, and the only place the mask is 4-connected around nothing |
| margin | weakest echo cell 0.551, strongest inactive cell 0.429 — the mask survives ±0.02 on the reference |

The lattice is not decoration: one mask cell is exactly two native samples, on the native sample grid, so *the mask is a coarsened read of the same chip* is something the plate proves rather than something the caption claims. And the mask is not a threshold of the picture either: the chip outside `W_{a,c}` does carry returns above the reference level (peak 0.685) and they produce no cell at all, because the mask is evaluated on the support and nowhere else.

The centroid and the peak are measured and reported — they are 1.95 cells apart, and that distance is what `S` measures — and neither is drawn. No dot, no arrow, no symbol: the reader is not asked to define a glyph the plate cannot justify.

## The vector, in its two blocks

| block | components | colour | what it is about |
| --- | --- | --- | --- |
| cross-acquisition | `γ²`, `P` | blue `#1f4e79` | what changes between passes of the same chip |
| spatial geometry | `D`, `F`, `A`, `S` | green `#3f7d3f` | what the echo looks like inside the support |

The two blocks are separated three ways — colour, a hairline rule in the gap, and an under-brace with a name — and each brace spans the measured extent of its own group, so a brace can never be wider than what it names. The row is the **same row** the first schematic draws, to 10⁻⁶ of a panel unit: a check reads `fig_ocv_schematic.json` and compares.

## What the plate does not show

It is a schematic, so it carries no site, no date, no state and no number that could be read as a result. Deliberately absent:

- structural health, damage probability, failure detection or any other diagnosis — the figure defines a measurement, it does not deliver a verdict;
- displacement, deformation or InSAR geometry — nothing in the plate moves;
- `R_B` and `R_G` — those belong to the temporal reconfiguration analysis (`figures/sds/`), not to the definition of the vector;
- any shift marker (`p̄`, `p⋆`) — the vector has no shift parameter on this plate;
- a centroid or a peak symbol — `S` is a distance the mask yields, not a glyph the plate has to explain;
- the reference level itself, any count, any weight and any score: the number is in this report and in the JSON, not on the plate.

The prohibition is mechanical, not editorial: a forbidden-token scan runs over every drawn string, the whole text inventory is compared against a written-out whitelist, and a separate scan rejects the shift markers, so the plate cannot grow a `damage` label or a `p̄` without failing the run.

## The plates, and how to reproduce them

| file | what it is |
| --- | --- |
| `concept/fig_ocv_measurement.jpg` | 183 mm × 61 mm, 630 dpi JPEG (quality 95, 4:4:4 subsampling) |
| `concept/fig_ocv_measurement.pdf` | the same plate as a true vector page (FontType 42; no raster object anywhere in the file — the intensity field is 2,028 filled quads) |
| `concept/fig_ocv_measurement.md` | this report |
| `../data/concept/fig_ocv_measurement.json` | the drawn geometry, the field, the mask, the text inventory, every check and the sha256 of both images |

```bash
bash figures/concept/fig_ocv_measurement.sh                  # both images
python3 code/concept/fig_ocv_measurement.py --no-figures     # checks only
PYTHON=/usr/bin/python3 bash figures/concept/fig_ocv_measurement.sh
```

For the manuscript:

```latex
\includegraphics[width=170mm]{figures/concept/fig_ocv_measurement.pdf}
```

## The checks

**160 of 160 checks pass.** A failure exits with status 2 *after* the artefacts are written, so a broken plate is reported rather than silently delivered. (`--no-figures` builds, measures and checks, writes the JSON and this report, and skips both images and the file checks that go with them.)

| family | what it pins |
| --- | --- |
| canvas | 183 mm × 61 mm at 630 dpi, ratio 3.000, `savefig.bbox` `standard`; every panel isotropic, so one unit is one unit in x and in y and the native cells stay square |
| one frame | the chip, the asset footprint, the radar frame and the scale are identical in all three panels — compared on the *drawn* geometry in device pixels, and on the chip's offset inside each panel |
| text inventory | exactly the declared strings and no others; one 7 pt title per panel; no plate title; no legend; all axes off |
| legibility | every text ≥ 7.0 pt at 183 mm and ≥ 6.5 pt at 170 mm; every text inside its own panel; no two texts overlap; every title fits its run |
| vocabulary | no forbidden token in any drawn string, and no shift marker (`p̄`, `p⋆`); the radar frame is `az` / `r`; the physical frame is absent |
| set algebra | `W_{a,c}` inside the guard inside the ring inside `C_{a,c}`, each with a visible margin; the guard ≥ 1.1 mm and the band outside it ≥ 1.7 mm; the asset footprint lies inside `W` |
| the lattice | one mask cell is exactly 2 × 2 native samples, on the native grid; the lattice is inside `W`; every mask cell is inside `W`; `C` holds more than twice the cells `W` does |
| raster integrity | 52 × 39 = 2,028 quads that tile the chip exactly once and touch every edge; cells square and ≥ 0.8 mm; the ramp is 16 monotone steps, neither ink nor paper; the field is bit-identical on a rebuild |
| the mask | 4 connected regions of 8, 6, 4 and 2 cells, strictly decreasing; 20 cells, 18.5 % density; exactly one enclosed hole, off the border; every echo cell above the reference level and every inactive cell below it, with 0.02 to spare; the mask is the threshold set of the field cell for cell |
| the mask stops at W | the chip outside the support carries returns above the reference level and produces no cell at all — the mask is a measurement, not a threshold of the picture |
| the four zones | read off the composed canvas, not the geometry: the outer zone is the panel-(a) colour pulled 65 % to white, the band and the guard are exact composites, the inside of `W` is the field under 18 % blue, and every one of them to within 2/255 |
| nothing annotated | no centroid or peak symbol, no marker artist anywhere in panel (c); the centroid and the peak are measured, reported and ≥ 1.5 cells apart — the quantity `S` measures |
| the vector | six components in the order `γ², P, D, F, A, S`; the partition into (γ², P) and (D, F, A, S) proved from the token table; each brace spans its own group's measured extent; the row is centred, clear of the mask, and identical to the first schematic's row |
| determinism, files | two independent builds render identical pixels and lay the row out byte-identically; JPEG 4 539 × 1 513 px ± 2 at 630 dpi; the PDF has no `/Subtype /Image`, embeds `/FontFile2`, carries no Type-3 font, no timestamp, and a 183 mm × 61 mm page box |

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
| 14 | the chip rectangle maps to the same size in all three panels | ok | [1048.685  786.514] px |
| 15 | the chip sits at the same place inside all three panels | ok | corner offset [193.179 400.156] px from each panel's origin |
| 16 | one chip rectangle, shared by all three panels | ok |  |
| 17 | all three panels draw the same asset footprint | ok |  |
| 18 | panels (a) and (b) draw the same radar frame in the same place | ok |  |
| 19 | panel (c) names no axes, because it labels nothing on the mask | ok |  |
| 20 | panel (c) draws no chip plate: the vector row runs through that band | ok |  |
| 21 | the frame encloses the chip on all four sides | ok | x 0.100..0.920 y 0.170..0.900 |
| 22 | the frame stays inside the panel | ok | {'x0': 0.1, 'y0': 0.17, 'x1': 0.92, 'y1': 0.9} |
| 23 | the frame's arm ends are clear of the chip's own edge | ok | r arm 0.120 below the chip |
| 24 | panel (c)'s ground is the support alone, not a second frame | ok |  |
| 25 | panels (b) and (c) draw the same W_{a,c} | ok | 8 vertices |
| 26 | panel (a) draws no support, so the restriction is visible | ok |  |
| 27 | the chip covers the same fraction of each panel | ok | 0.760000 of the panel width |
| 28 | the native sample grid is the same grid in (a) and (b) | ok | 52 x 39 cells |
| 29 | text inventory is exactly the declared set | ok | 24 objects |
| 30 | text count <= 30 | ok | 24 |
| 31 | no duplicate panel label | ok |  |
| 32 | every panel carries exactly one title | ok |  |
| 33 | no title anywhere | ok |  |
| 34 | no legend anywhere | ok |  |
| 35 | axes are off (no spines, no ticks) | ok |  |
| 36 | every text >= 7.0 pt | ok | min 7.00 pt |
| 37 | every text >= 7.0 pt at 170 mm (scale 0.929) = 6.5 pt | ok | 6.50 pt |
| 38 | panel a: every text inside the panel | ok |  |
| 39 | panel a: no two texts overlap | ok |  |
| 40 | panel a: the title fits its own run | ok | 0.115..0.577 |
| 41 | panel b: every text inside the panel | ok |  |
| 42 | panel b: no two texts overlap | ok |  |
| 43 | panel b: the title fits its own run | ok | 0.115..0.928 |
| 44 | panel c: every text inside the panel | ok |  |
| 45 | panel c: no two texts overlap | ok |  |
| 46 | panel c: the title fits its own run | ok | 0.115..0.681 |
| 47 | (b) names the asset on the asset: the symbol sits inside the footprint | ok | label 0.437..0.603 x 0.577..0.635 |
| 48 | (b) sets the asset symbol in paper white, on the dark footprint | ok |  |
| 49 | (b)'s W_{a,c} label and its knockout stay clear of the footprint | ok | knockout to 0.405, footprint from 0.396 at its base |
| 50 | no forbidden vocabulary in any drawn text | ok |  |
| 51 | R_B absent | ok |  |
| 52 | R_G absent | ok |  |
| 53 | no shift marker (p-bar, p-star) anywhere on the plate | ok |  |
| 54 | the coherence channel is mathtext | ok |  |
| 55 | the radar frame is az / r and nothing else | ok |  |
| 56 | the physical frame (X, Y, Z) is not on this plate at all | ok |  |
| 57 | every channel-dependent symbol carries c | ok | C_{a,c}, W_{a,c}, W_{a,c} subset C_{a,c}, Pi_c(G_a^{3D}) |
| 58 | no channel-dependent symbol is spelled out instead | ok |  |
| 59 | W has eight vertices (the chamfered rectangle of the first plate) | ok |  |
| 60 | W sits inside C with a visible margin | ok | margin 0.080 |
| 61 | W sits inside the guard | ok |  |
| 62 | the guard sits inside the ring | ok |  |
| 63 | the guard is a white band >= 0.015 wide | ok | 0.020 units = 1.11 mm |
| 64 | the ring band outside the guard is >= 0.03 wide | ok | 0.040 units = 2.23 mm |
| 65 | the ring sits inside C, so it cannot be the chip edge | ok | margin 0.020 |
| 66 | C sits inside the panel | ok |  |
| 67 | the asset footprint lies inside W | ok |  |
| 68 | the asset footprint is inside C, not straddling it | ok |  |
| 69 | the ring leader ends on the ring band, clear of the guard | ok | 0.800 < 0.812 < 0.840 |
| 70 | W in (b) and W in (c) are the same shape, in the same box | ok |  |
| 71 | one mask cell is exactly 2 x 2 native samples | ok | 0.029 / 0.015 = 2.0 |
| 72 | the mask lattice sits on the native sample grid | ok | origin at sample (11, 12) |
| 73 | the mask lattice sits inside W | ok | margin 0.022 |
| 74 | every mask cell is inside W | ok |  |
| 75 | C holds far more cells than W (the read is not the measurement) | ok | 800 of 2028 native cells in W |
| 76 | the mask has fewer cells than W has native samples | ok | 20 mask cells, 800 native cells |
| 77 | native cells are square on the drawn canvas | ok | 0.0146 units |
| 78 | native cell >= 0.8 mm on the printed page | ok | 0.81 mm |
| 79 | the chip holds 52 x 39 native samples | ok | 2028 |
| 80 | the raster is drawn as quads, not as an image | ok |  |
| 81 | the quads cover the chip exactly once | ok | 2028 quads |
| 82 | the quad union spans the chip from edge to edge | ok | [0.1396 0.2896 0.9004 0.8604] |
| 83 | no native cell is counted twice | ok |  |
| 84 | the intensity ramp has 16 steps | ok |  |
| 85 | the ramp is monotone and has no repeated step | ok | 87..240 |
| 86 | the darkest ramp step is grey, not ink | ok | 87/255 |
| 87 | the lightest ramp step is not paper white | ok | 240/255 |
| 88 | the field is deterministic (no random number anywhere) | ok |  |
| 89 | the field uses at least half the ramp | ok | 14 of 16 steps used |
| 90 | the brightest sample in the chip lies inside W | ok | (0.483, 0.648) = 0.782 |
| 91 | the chip outside W does carry bright returns | ok | peak outside W 0.685 > reference 0.490 |
| 92 | the strongest return in the chip is on the support | ok | 0.782 inside vs 0.685 outside |
| 93 | the echo mask has exactly four connected regions | ok | 4 |
| 94 | the four regions are 8, 6, 4 and 2 cells | ok | [8, 6, 4, 2] |
| 95 | the regions are strictly decreasing in size | ok |  |
| 96 | every region is 4-connected (a fragment, not a dust field) | ok |  |
| 97 | the mask has 20 cells | ok | 20 |
| 98 | the mask encloses exactly one hole | ok | 1 |
| 99 | the hole is enclosed, not on the lattice border | ok | [(5, 5)] |
| 100 | no echo cell touches the lattice border | ok |  |
| 101 | echo density is sparse and legible (5 % < d < 22 %) | ok | 0.1333 |
| 102 | every echo cell is above the reference level | ok | min echo 0.551 > 0.490 |
| 103 | no inactive cell reaches the reference level | ok | max inactive 0.429 < 0.490 |
| 104 | the mask is not knife-edge: it clears the level by 0.02 either way | ok | margins 0.061 / 0.061 |
| 105 | the mask is exactly the threshold set of the field, cell for cell | ok |  |
| 106 | the asset footprint carries the bulk of the echo | ok | 20 of 20 echo cells on the asset |
| 107 | no echo cell crosses the asset outline | ok | [] |
| 108 | the strongest echo cell lies inside W | ok | (6, 6) = 0.782 |
| 109 | echo stops at the support: bright returns outside W make no cell | ok | 20 cells, all inside W |
| 110 | the centroid and the peak are >= 1.5 cells apart | ok | 1.946 cells, the quantity S measures |
| 111 | the centroid of the mask lies inside the support | ok | (0.516, 0.604) in panel units |
| 112 | nothing is annotated on the mask: no marker, no dot, no label | ok | 7 line objects, all markers empty |
| 113 | (a) paints the field cell at its own ramp step | ok | rgb 0.502 0.502 0.502 |
| 114 | (b) shows the same cell outside the ring, pulled 65 % to white | ok | rgb 0.827 0.827 0.827 |
| 115 | (b) keeps the field inside W at full contrast under the blue veil | ok | rgb 0.396 0.427 0.459 |
| 116 | the veil is visible: (b) inside W is not (a) | ok |  |
| 117 | (b)'s guard is pure white, knocked out of the ring tint | ok | rgb 1.000 1.000 1.000 |
| 118 | (b)'s ring band is the flat tint, not the field | ok | rgb 0.914 0.914 0.914 |
| 119 | the ring tint never reaches inside W (annulus, not a fill) | ok |  |
| 120 | (b)'s footprint is a dark plate, so the white symbol reads on it | ok | rgb 0.141 0.145 0.149 |
| 121 | (c) draws no chip plate: outside the support is bare paper | ok | rgb 1.000 1.000 1.000 |
| 122 | (c)'s support ground is the pale blue over paper | ok | rgb 0.910 0.925 0.945 |
| 123 | (c)'s inactive support cells are very light grey | ok | rgb 0.914 0.914 0.914 |
| 124 | (c)'s echo cells are green | ok | rgb 0.247 0.490 0.247 |
| 125 | green appears in (c), and only there | ok | 37255 exact-green pixels on the whole plate |
| 126 | the green is only the mask: (a) and (b) carry no green at all | ok | left two thirds of the canvas |
| 127 | the vector has six components | ok |  |
| 128 | the order is gamma2, P, D, F, A, S | ok | gamma2,P,D,F,A,S |
| 129 | the partition covers every component exactly once | ok |  |
| 130 | group one is exactly gamma2 and P | ok | $\gamma^2,\ P$ |
| 131 | group two is exactly D, F, A and S | ok | $D,\ F,\ A,\ S$ |
| 132 | the row starts inside the panel | ok | 0.0339 |
| 133 | the row ends inside the panel | ok | 0.9661 |
| 134 | the row is centred in its panel | ok | 0.5000 |
| 135 | the cross brace spans its own group only | ok | text 0.312..0.441 brace 0.312..0.441 |
| 136 | the cross name sits under its brace | ok | top 0.279 < brace 0.285 |
| 137 | the spatial brace spans its own group only | ok | text 0.648..0.922 brace 0.648..0.922 |
| 138 | the spatial name sits under its brace | ok | top 0.279 < brace 0.285 |
| 139 | the separating comma sits between the two groups | ok | 0.453..0.471 |
| 140 | the separating comma does not touch either group | ok | 0.012 / 0.177 |
| 141 | the two group names do not overlap | ok | 0.562 < 0.597 |
| 142 | the equation sits above the braces | ok |  |
| 143 | the equation does not touch the mask lattice | ok | row top 0.378 < lattice 0.465 |
| 144 | the arrow leaves the support's lower edge and points at the row | ok | 0.382 -> 0.418 under support 0.430 |
| 145 | the whole vector row sits below the support | ok | row top 0.378 < support 0.430 |
| 146 | JPEG width is 183 mm at 630 dpi (4539 px +/- 2) | ok | 4538 px |
| 147 | JPEG height is 61 mm at 630 dpi (1513 px +/- 2) | ok | 1512 px |
| 148 | JPEG pixel ratio is 3.000 | ok | 3.00132 |
| 149 | JPEG carries 630 dpi | ok | (630, 630) |
| 150 | JPEG is a full-colour raster | ok | RGB |
| 151 | PDF page contains no raster image | ok | the raster is 2028 vector quads, not an image object |
| 152 | PDF embeds its fonts as TrueType (FontType 42) | ok |  |
| 153 | PDF has no Type-3 fonts | ok |  |
| 154 | PDF carries no timestamp (it must be reproducible) | ok |  |
| 155 | PDF is a manageable vector page (< 1.5 MB) | ok | 0.10 MB |
| 156 | PDF page box is 183 x 61 mm | ok | [0.0, 0.0, 518.7401574803, 172.9133858268] |
| 157 | re-saving the same figure gives the same JPEG bytes | ok |  |
| 158 | two builds render identical pixels | ok | (1512, 4538, 4) |
| 159 | two builds lay the OCV row out identically | ok |  |
| 160 | the vector row agrees with the first schematic | ok | left 0.033885 vs 0.033885, width 0.932230 vs 0.932230 |

## Files

- `fig_ocv_measurement.jpg` — sha256 `2a4974802f9f2affef734094d4f72911f3089212b1374a382a11a527ad0d8630`
- `fig_ocv_measurement.pdf` — sha256 `f66976d1823b7964f9b6a38b05bcd8ca8237950017ae27a902eaad40c741cb7c`

Environment: Python 3.13.5, matplotlib 3.11.1, numpy 2.4.6. The plate is deterministic: no timestamp, no seed-dependent value and no machine-dependent path reaches any output, so two runs produce byte-identical JSON, Markdown, JPEG and PDF.
