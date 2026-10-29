# `concept/fig_ocv_schematic` — the observability core vector, drawn

![the observability core vector](fig_ocv_schematic.jpg)

One plate, three panels, no data. Every other figure in `figures/` is computed from a record and pinned to it; this one is computed from the **definition** the records are measured with, and its job is to say what the core vector is before any site's numbers are read. It therefore has no pins to a table — it has structural checks instead (97 of them, listed below and in `data/concept/fig_ocv_schematic.json`).

```
physical asset -> SAR observation geometry -> asset-specific measurement
support -> radar echo representation -> OCV
```

## The three panels

**(a) the asset and the acquisition.** A generic lattice asset `G_a^{3D}` in its own frame (the grey triad `X`, `Y`, `Z`), observed from a satellite along the orange line of sight `LOS_c` at the look angle `θ_c` measured from nadir; the operator `Π_c` maps that frame into the radar frame `(az, r)`. The asset is drawn as four legs, five belts and alternating braces — something with a real cross-section, because that is what makes the projection lossy.

The subscript on `LOS_c` and `θ_c` is the acquisition channel: the same `c` that `Π_c`, `W_{a,c}` and `C_{a,c}` carry, because the line of sight and the incidence angle belong to the pass rather than to the asset. Nothing is spelled out about ascending or descending — the subscript says it, and a check fails the build if a channel-dependent symbol ever loses it or a frame symbol picks it up.

The mapping is written out exactly once, at the foot of the panel, left-aligned to the radar frame's own left edge: `(X, Y, Z) ⟶^{Π_c} (az, r)`. The operator appears nowhere else on the plate, and the two frames share no symbol — the physical frame is `X`, `Y`, `Z`, the radar frame is `az`, `r`, and a check fails the build if the two symbol sets ever intersect. The radar frame is drawn 0.040 panel units above the foot so that the statement has clear air beneath it, and a check pins both the gap and the alignment.

**(b) the measurement support.** The projected asset `Π_c(G_a^{3D})` sits inside the asset- and channel-specific support `W_{a,c}` (blue), which sits inside the read chip `C_{a,c}` (grey, outermost). Outward from the support, in this order, come the **guard** — a white band 0.020 units wide (1.1 mm) — and then the **background ring**, dashed and lightly tinted, from which the no-echo reference is read.

The ring is drawn as a true annulus: its tint is knocked out over the guard, so no tint reaches inside `W_{a,c}` and no background cell can be read as support. That order — support, guard, ring, chip — is the order the audit walks (`code/lumo/lumo_ocv_masks.py`), and two pixel checks prove it on the composed canvas rather than on the geometry: W's interior must equal the blue fill at 0.13 over white and the outer band must equal the ring tint, each to within 2/255.

`C_{a,c}` is set *on* the chip's own lower edge with a white knockout, so the edge breaks behind the glyphs and the label cannot be read as naming the ring it used to float between; a check pins the label's own box across that edge. The panel states the set algebra `W_{a,c} ⊂ C_{a,c}` once at its foot and nothing else: no box, no legend, no count.

**(c) the observable representation and the vector.** Inside the same `W_{a,c}` the echo mask is a raster of green cells: 4 connected regions of 8, 6, 4, 3 cells in a 16 × 10 grid (density 0.131) and one closed hole. That is the extent, the density and the fragmentation the vector has to summarise, made visible before it is made numeric.

Nothing is annotated on the mask itself — no centroid, no peak, no arrow between them, no symbol naming either. The centroid-to-peak distance that `S` measures is a number the raster yields, not an element the plate has to define: every mark in this panel can be read without a caption, and the only words on it are the vector's and its two group names. Below the mask:

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
- a centroid or a peak symbol — `S` is a distance the mask yields, not a glyph the plate has to explain;
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

## The checks

**97 of 97 checks pass.** A failure exits with status 2 *after* the artefacts are written, so a broken plate is reported rather than silently delivered. (`--no-figures` builds, measures and checks, writes the JSON and this report, and skips both images and the eleven file checks that go with them.)

| family | what it pins |
| --- | --- |
| canvas | 183 mm × 61 mm at 630 dpi, ratio 3.000, `savefig.bbox` `standard` (a tight bbox would move the 183 mm edge); every panel isotropic, so one unit is one unit in x and in y and the raster cells are square |
| text inventory | exactly the declared 26 strings and no others; no title; no legend; all axes off |
| legibility | every text ≥ 7.0 pt at 183 mm and ≥ 6.5 pt at 170 mm; every text inside its own panel; no two texts overlap |
| vocabulary | no forbidden token in any drawn string (health, damage, fail, displace, deform, InSAR, score, threshold, px, `R_B`, `R_G`, …) |
| set algebra | `W_{a,c}` inside the guard inside the ring inside `C_{a,c}`, each with a visible margin; the guard is ≥ 0.015 units and the tinted band outside it ≥ 0.03; the projected asset lies inside `W` |
| the annulus | the composed canvas is sampled, not the geometry: W's interior must be the blue fill over white and the outer band the ring tint, each within 2/255, so the tint provably never runs behind the support |
| the two frames | the physical frame (`X`, `Y`, `Z`) and the radar frame (`az`, `r`) share no symbol; the mapping statement carries the operator and both frames; `C_{a,c}` sits on the chip's own edge |
| the channel | every channel-dependent symbol carries the `c` subscript — `LOS_c`, `θ_c`, `Π_c`, `W_{a,c}`, `C_{a,c}` — and no frame symbol does; ascending and descending are never spelled out |
| one W, two panels | the support is **one** polygon definition mapped into (b) and into (c); normalised to its own box it is identical to 1 × 10⁻¹² in both |
| raster integrity | cells square on the drawn canvas and ≥ 1.5 mm; every echo cell inside `W` and inside the raster; exactly 4 connected regions of 8, 6, 4 and 3 cells; exactly one enclosed hole; density in (5 %, 20 %); the mask's centroid and its peak ≥ 1.5 cells apart — the quantity `S` measures, measured but not drawn |
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
| 14 | text inventory is exactly the declared set | ok | 26 objects |
| 15 | text count <= 30 | ok | 26 |
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
| 26 | the physical frame and the radar frame share no symbol | ok | [] |
| 27 | the mapping statement names the operator and both frames | ok | no drift |
| 28 | every channel-dependent symbol carries c, and no frame symbol does | ok | 7 channel symbols |
| 29 | panel a: every text inside the panel | ok |  |
| 30 | panel a: no two texts overlap | ok |  |
| 31 | panel b: every text inside the panel | ok |  |
| 32 | panel b: no two texts overlap | ok |  |
| 33 | panel c: every text inside the panel | ok |  |
| 34 | panel c: no two texts overlap | ok |  |
| 35 | C_{a,c} sits on the chip's own lower edge, not between chip and ring | ok | y 0.216..0.264 across the edge 0.240 |
| 36 | the mapping statement is the radar frame's own caption, below it | ok | top 0.088 < frame 0.140, left 0.560 vs 0.560 |
| 37 | W has eight vertices in both panels | ok |  |
| 38 | W in (b) sits inside C with a >= 20 um margin | ok |  |
| 39 | W in (b) sits strictly inside the background ring | ok |  |
| 40 | W in (b) sits inside the guard, and the guard inside the ring | ok |  |
| 41 | the guard is a white band >= 0.015 wide between W and the ring | ok | 0.020 units |
| 42 | the ring band outside the guard is >= 0.03 wide | ok | 0.040 units |
| 43 | the ring sits inside C, so it cannot be the chip edge | ok |  |
| 44 | C sits inside the panel | ok |  |
| 45 | the projected asset lies inside W | ok |  |
| 46 | the projected asset is inside C, not straddling it | ok |  |
| 47 | the ring leader ends on the ring band, clear of the guard | ok | y=0.770 |
| 48 | the ring's tint does not reach inside W (annulus, not a fill) | ok | W rgb 0.886 0.910 0.929 |
| 49 | the band between the guard and the ring's boundary is tinted | ok | band rgb 0.914 0.914 0.914 |
| 50 | W in (b) and W in (c) are the same polygon | ok |  |
| 51 | W in (b) and W in (c) are the same shape, not the same place | ok |  |
| 52 | raster cells are square on the drawn canvas | ok | 0.03625 units |
| 53 | raster cell >= 1.5 mm on the printed page | ok | 2.02 mm |
| 54 | the raster sits inside W in (c) | ok |  |
| 55 | the raster sits inside W with a visible margin | ok |  |
| 56 | every echo cell lies inside W in (c) | ok |  |
| 57 | every echo cell lies inside the raster | ok |  |
| 58 | the echo mask has exactly four connected regions | ok | 4 |
| 59 | the four regions are 8, 6, 4 and 3 cells | ok | [8, 6, 4, 3] |
| 60 | every region is 4-connected (a fragment, not a dust field) | ok |  |
| 61 | the mask encloses exactly one hole | ok | 1 |
| 62 | the hole is the centre of the 3x3 ring | ok | [(4, 4)] |
| 63 | echo density is sparse and legible (5 % < d < 20 %) | ok | 0.1313 |
| 64 | the peak lies on an echo cell | ok |  |
| 65 | the centroid and the peak are >= 1.5 cells apart | ok | 2.626 cells, the quantity S measures |
| 66 | the centroid is not inside an echo cell (the mask is fragmented) | ok |  |
| 67 | the vector has six components | ok |  |
| 68 | the order is gamma2, P, D, F, A, S | ok | gamma2,P,D,F,A,S |
| 69 | the partition covers every component exactly once | ok |  |
| 70 | group one is exactly gamma2 and P | ok | $\gamma^2,\ P$ |
| 71 | group two is exactly D, F, A and S | ok | $D,\ F,\ A,\ S$ |
| 72 | the row starts inside the panel | ok | 0.0339 |
| 73 | the row ends inside the panel | ok | 0.9661 |
| 74 | the row is centred in its panel | ok | 0.5000 |
| 75 | the cross brace spans its own group only | ok | text 0.312..0.441 brace 0.312..0.441 |
| 76 | the cross name sits under its brace | ok | top 0.279 < brace 0.285 |
| 77 | the spatial brace spans its own group only | ok | text 0.648..0.922 brace 0.648..0.922 |
| 78 | the spatial name sits under its brace | ok | top 0.279 < brace 0.285 |
| 79 | the separating comma sits between the two groups | ok | 0.453..0.471 |
| 80 | the separating comma does not touch either group | ok | 0.012 / 0.177 |
| 81 | the two group names do not overlap | ok | 0.562 < 0.597 |
| 82 | the two group names do not run into each other's brace | ok | gap 0.035 |
| 83 | the equation sits above the braces | ok |  |
| 84 | the equation does not touch the echo raster | ok |  |
| 85 | JPEG width is 183 mm at 630 dpi (4539 px +/- 2) | ok | 4538 px |
| 86 | JPEG height is 61 mm at 630 dpi (1513 px +/- 2) | ok | 1512 px |
| 87 | JPEG pixel ratio is 3.000 | ok | 3.00132 |
| 88 | JPEG carries 630 dpi | ok | (630, 630) |
| 89 | JPEG is a full-colour raster | ok | RGB |
| 90 | PDF page contains no raster image | ok |  |
| 91 | PDF embeds its fonts as TrueType (FontType 42) | ok |  |
| 92 | PDF has no Type-3 fonts | ok |  |
| 93 | PDF carries no timestamp (it must be reproducible) | ok |  |
| 94 | PDF page box is 183 x 61 mm | ok | [0.0, 0.0, 518.7401574803, 172.9133858268] |
| 95 | re-saving the same figure gives the same JPEG bytes | ok |  |
| 96 | two builds render identical pixels | ok | (1512, 4538, 4) |
| 97 | two builds lay the OCV row out identically | ok |  |

## Files

- `fig_ocv_schematic.jpg` — sha256 `94c94bddce41ba3537adc34594d69e5b4a83005de2e11443f3348a275a242b42`
- `fig_ocv_schematic.pdf` — sha256 `eafc0e62bfa0843fe5b33dcfad68ae24a292cb662466719c8f8edfe66ad1a175`

Environment: Python 3.13.5, matplotlib 3.11.1, numpy 2.4.6. The plate is deterministic: no timestamp, no seed-dependent value and no machine-dependent path reaches any output, so two runs produce byte-identical JSON, Markdown, JPEG and PDF.
