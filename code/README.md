# `code/` — the four site packages

One subfolder per site, all standalone, offline and built the same way (stats
layer, data layer, generator, figure script). They differ in what the observability
core vector contains: LUMO masks the **echo** of the tower, Bautzen masks the
**deck edge** of the bridge, KDLO consumes the mask size `A` that the
committed record delivers (it has no mask layer of its own), and Carola masks
both the **echo of every chip** and the **per-girder deck** of the collapsed
bridge.

| folder | site | states | mask layer | package README |
| --- | --- | --- | --- | --- |
| [`lumo/`](lumo/README.md) | LUMO lattice tower, Hannover | healthy, DAM3, DAM4, DAM6 | echo mask (`0.30 x peak`, `5 x median`, ≥ 2 px) | `lumo/README.md` |
| [`bautzen/`](bautzen/README.md) | Bautzen OpenLabs bridge | RS, DS1, DS2 | deck-edge mask (geodetically anchored deck band) | `bautzen/README.md` |
| [`kdlo/`](kdlo/README.md) | KDLO media tower, Garden City SD | pre-collapse, during-rebuild (epochs) | echo mask recomputed from a committed 64-strip full-dwell cache (LUMO rule, upper median) | `kdlo/README.md` |
| [`carola/`](carola/README.md) | Carolabrücke, Dresden (collapsed 2024-09-11) | pre-collapse (healthy), post-collapse | two layers: the echo mask of every chip **and** the per-girder deck mask, both from a committed window-payload cache (`0.30 x peak`, `5 x np.median`, ≥ 2 px) | `carola/README.md` |

The scripts are started from the **repository root** (each wrapper in
`../figures/<site>/` does that), so that their defaults resolve to
`data/<site>/…` in and `data/<site>/…`, `figures/<site>/…` out. No package:
the script's own folder is on `sys.path` (a script run puts its directory first),
which is why the modules import each other by their `lumo_…` / `bautzen_…` /
`kdlo_…` / `carola_…` names.

Adding another site means copying one package and swapping the mask layer (if the
record has one) plus the `DATA`/`FIGDIR` lines of `…_ocv_core.py` — nothing
outside the new `code/<site>/`, `data/<site>/` and `figures/<site>/` changes. The
KDLO package is the example of a site *without* a mask layer: the two OCV
components are `gamma2` and the delivered mask size. The Carola package is the
example of a site *with two* mask layers over one chip: the echo mask of the
whole 80 x 80 window and the same rule on the girder's own window.

