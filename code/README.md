# `code/` — the eight site packages

One subfolder per site, all standalone, offline and built the same way (stats
layer, data layer, generator, figure script). They differ in what the observability
core vector contains: LUMO masks the **echo** of the tower, Bautzen masks the
**deck edge** of the bridge, KDLO consumes the mask size `A` that the
committed record delivers (it has no mask layer of its own), and Carola, Morandi
and CTS each mask both the **echo of every chip** and the **per-girder /
per-track deck** of a collapsed structure. CTS is the one
**difference-in-differences** site: four footprints of the *same* orbit (target
plus two control towers) share strata, dates and counts, so its package carries a
first-class DiD layer (`cts_ocv_did.py`) next to the mask layer. YWF is the one
site whose record is too small and too duplicated for a contrast at all: one
mask layer, echo **coverage** instead of effect sizes, and a reference file that
the package itself defines and then pins. Espoo is the one site with **no event at
all**: 150 acquisitions with an empty `damage_label`, `structural_state` and
`condition_label` on every row and **no mask raster** either (56 of its 78 columns
are empty), so its OCV is `[gamma2, A]` with `A` = the delivered mask *size*,
**copied and proved** rather than recomputed, and the only contrast left is the
orbit geometry (ASC / DESC) — the family's event-free control.

| folder | site | states | mask layer | package README |
| --- | --- | --- | --- | --- |
| [`lumo/`](lumo/README.md) | LUMO lattice tower, Hannover | healthy, DAM3, DAM4, DAM6 | echo mask (`0.30 x peak`, `5 x median`, ≥ 2 px) | `lumo/README.md` |
| [`bautzen/`](bautzen/README.md) | Bautzen OpenLabs bridge | RS, DS1, DS2 | deck-edge mask (geodetically anchored deck band) | `bautzen/README.md` |
| [`kdlo/`](kdlo/README.md) | KDLO media tower, Garden City SD | pre-collapse, during-rebuild (epochs) | echo mask recomputed from a committed 64-strip full-dwell cache (LUMO rule, upper median) | `kdlo/README.md` |
| [`carola/`](carola/README.md) | Carolabrücke, Dresden (collapsed 2024-09-11) | pre-collapse (healthy), post-collapse | two layers: the echo mask of every chip **and** the per-girder deck mask, both from a committed window-payload cache (`0.30 x peak`, `5 x np.median`, ≥ 2 px) | `carola/README.md` |
| [`morandi/`](morandi/README.md) | Ponte Morandi / Polcevera, Genova (collapsed 2018-08-14) | pre-collapse (healthy), post-collapse | two layers: the echo mask of every chip **and** the per-track deck mask, both from a committed window-payload cache (`0.30 x peak`, `5 x np.median`, ≥ 2 px) | `morandi/README.md` |
| [`cts/`](cts/README.md) | Champlain Towers South, Surfside (collapsed 2021-06-24) | pre-collapse (healthy), post-collapse | two layers over four footprints of one orbit (target + two controls + a beach reference), the echo mask **and** the per-track deck mask from a committed window cache (`0.30 x peak`, `5 x np.median`, ≥ 2 px), **plus a difference-in-differences layer** | `cts/README.md` |
| [`ywf/`](ywf/README.md) | Yeongdeok Wind Farm, Unit 21, Samgye-ri (collapsed 2026-02-02) | pre-collapse (healthy), post-collapse | **one layer**: the echo mask of the committed 7 x 7 window payload (`0.30 x peak`, `5 x np.median`, ≥ 2 px) — the record is not segment-resolved, so a per-mast-section layer would be dishonest; its second cross-cut (the pass geometry, ASCENDING / DESCENDING) is reported as a pinned, explicitly descriptive echo-coverage split; two further references are *resolved* from public catalogues (the CDSE burst pass of every chip, the OSM road anchor) and committed beside it | `ywf/README.md` |
| [`espoo/`](espoo/README.md) | Espoo Kurttila mast, Espoo, Finland (**event-free** 2024-09-05 … 2026-09-07) | none — the record has no damage label at all; the contrast is the orbit geometry (ASCENDING n = 70 / DESCENDING n = 80) | **no mask layer**: the record delivers only the mask *size* `A` (`coherence_masked_pixels`) and no raster (56 of 78 columns empty), so the OCV is `[gamma2, A]` and `A` is **copied and proved**, never recomputed | `espoo/README.md` |

The scripts are started from the **repository root** (each wrapper in
`../figures/<site>/` does that), so that their defaults resolve to
`data/<site>/…` in and `data/<site>/…`, `figures/<site>/…` out. No package:
the script's own folder is on `sys.path` (a script run puts its directory first),
which is why the modules import each other by their `lumo_…` / `bautzen_…` /
`kdlo_…` / `carola_…` / `morandi_…` / `cts_…` / `ywf_…` / `espoo_…` names.

Adding another site means copying one package and swapping the mask layer (if the
record has one) plus the `DATA`/`FIGDIR` lines of `…_ocv_core.py` — nothing
outside the new `code/<site>/`, `data/<site>/` and `figures/<site>/` changes. The
KDLO package is the example of a site *without* a mask layer: the two OCV
components are `gamma2` and the delivered mask size. The Espoo package is the
example of a site *without an event*: the same two components (`gamma2` and the
delivered `coherence_masked_pixels`), no damage label on any row, no mask raster
behind the size, and the orbit geometry (ASC / DESC) as the only contrast the
record supports — so its `espoo_ocv_channels_csv.py` **proves** the copied `A`
against the committed column and the committed reference statistics instead of
recomputing a mask, and its reports read every effect size against the record's own
no-event controls (stationarity, season, weather, coverage). The Carola package is the
example of a site *with two* mask layers over one chip: the echo mask of the
whole 80 x 80 window and the same rule on the girder's own window. The Morandi
package is the same shape one more time, with the two Sentinel-1 tracks
(`A_asc` / `A_des`) in place of the girders; because only 11 of its 232 chips
survive the collapse, its figures are an honest null rather than a contrast. The
CTS package keeps the same shape but adds what no single-structure site can: the
four footprints of its one orbit (the collapsed tower as `target`, two standing
towers as `control_ctn` / `control_cte`, the beach as `reference_beach`) let it
difference the target against a control that shares the atmosphere, so
`code/cts/cts_ocv_did.py` pins the three committed
`data/cts/reference/cts_reference_did_*.json` files (series and markdown, 1e-9).
Its figure script `code/cts/fig_cts_ocv_paper.py` (run by
`figures/cts/fig_cts_ocv_paper.sh`) mirrors the other sites' panels A–E and adds
panel **F**, the DiD forest — the interaction t per channel for both controls and
the placebo, the differencing assumption and the pre-trend parallelism — so the
whole contrast, and its honest null, is readable in one place. The YWF package is
the family's smallest and the only one with a single mask layer: its 7 x 7 window
record holds 33 unique chips (24 windows are byte-identical duplicates across the
requested mast sections, `segment_resolved = false`), only 17 of them carry an
echo at all, and 1 of the 6 post-event chips does — so `fig_ywf_ocv_paper.py` fits
no classifier, leads with the echo-coverage table and its Fisher test, adds §A2
with the same coverage split by pass geometry (10 of 20 ascending vs. 7 of 13
descending chips, p = 1.0000 — pinned, but read as a description of where the
echoes sit), and
*defines* `data/ywf/reference/ywf_ocv_findings.json` with `--write-reference`
before pinning it on every later run.

Two further references are not *defined* but *resolved* against public catalogues
and committed beside it — `reference/ywf_bursts.json` (the CDSE pass of all 33
chips: sub-swath, polarisation, relative orbit, azimuth time, with ASCENDING =
IW3 / rel-orbit 54 at 09:2x UT and DESCENDING = IW2 / rel-orbit 61 at 21:2x UT, VV
and VH interleaved in both) and `reference/ywf_osm_road.json` +
`reference/ywf_geometry.json` (the OSM way the tower fell on, its centreline and
the estimated burst grid) — each with its own offline `--check`
(`ywf_bursts_resolve.py`, `ywf_osm_anchor.py`) and its sha256 digest in the
figures' pin block. The package also states what the record **cannot** do: a
7 x 7 px window per acquisition (`TOWER_WINDOW_SIZE`), no rect spec of its own, no
per-section coordinates and no deck edge — the road lies *outside* the frame (its
near edge is 8.10-9.33 px across range from the window centre against a 6 px
reach, so containment needs 19-21 px per side) and a ~7 m carriageway would still
be 2.0 px across range but 0.5 px across azimuth — which is why the echo mask
stays the only layer it claims, and why the mask rule itself is only ever
*studied*: `ywf_mask_sensitivity.py` sweeps it offline on the committed payload
(the rule is family-wide — identical in the backend's `sarsolve` and in Carola's
`gamma2`).

