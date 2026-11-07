# Individual shrub outputs (`shrub_outputs.py`)

Every shrub in a plot gets its own set of outputs, mirroring the per-tree ones,
under `individual_shrubs/<plot_uid>_shrub_<uid>/`. These are written by
`understory.py` right after the understory cloud (shrubs are generated there),
so they need no separate stage. Disable with `understory.py --no-shrub-outputs`.

## Layout

```
individual_shrubs/<plot>_shrub_<uid>/
  las/
    <tag>_raw.laz           the shrub's own points only (semantic 8)
    <tag>_understory.laz    the shrub + surrounding context (bbox crop)
    <tag>_normalized.laz     written if a normalize run is performed
  deep_learning/
    volume/labels/<tag>_volume.json
    canopy_map/features/<tag>_canopy_bounding_box_raw.tif   (crown CHM)
    canopy_map/labels/<tag>_canopy_mask_raw.tif
    canopy_map/labels/<tag>_canopy_bounding_box.json
    canopy_map/labels/<tag>_canopy_convex_hull.json
    stem_map/features/<tag>_stem_bounding_box_raw.tif       (DBH-band count)
    stem_map/labels/<tag>_stem_mask_raw.tif
```

A plot-level `plot_placement/textual/<plot>_shrub_index.json` lists every shrub
with its id, point count, and the three volumes.

## Clouds

- **raw** — only this shrub's points (`semantic_code == 4` and its `shrub_id`).
- **understory** — every point (any class) inside the shrub's footprint bbox
  padded by one shrub-width, so the shrub sits in its neighbourhood (grass,
  ground, and any overlapping tree points). Mirrors the tree understory crop.
- **normalized** — written when `normalize.py` runs, using the same steps and
  grid (`--tree-l`/`--tree-z`) as trees. Shrubs have no `trimmed`/`terrain`
  source, so an `understory` source maps to the shrub's understory cloud and any
  other source maps to its raw cloud.

## Volume JSON

Three volume estimates per shrub, plus the generation parameters:

| field | meaning |
|---|---|
| `volume_m3.primitive` | nominal ellipsoid `(4/3)pi (w/2)^2 (h/2)` from the shrub's height/width, before jitter |
| `volume_m3.mesh` | sampled-shape volume: the actual (jittered) primitive mesh that produced the shrub's points (trimesh `.volume`, or the convex-hull volume of the mesh vertices as a fallback) — the shrub volume per the notes.md ruling |
| `height_m`, `width_m`, `center_xy` | the generated shrub geometry |

## Canopy and DBH maps

- **Canopy map** — the shrub's crown, built exactly like the tree canopy map: a
  max-height raster (CHM) over the shrub points, a filled crown mask, and a
  bbox + convex-hull JSON (hull volume/area when the shrub has >= 4 points). Like
  trees, the CHM stores absolute z.
- **DBH (stem) map** — the cross-section of shrub points in the breast-height
  band, using the **same DBH band parameters as trees** (`dbh_height`,
  `threshold`), measured relative to the shrub's own base. A shrub has no trunk,
  so this is usually empty or a small patch -- that's expected and kept only for
  parity with the tree stem map, so downstream code can treat shrubs and trees
  uniformly.

## Parameters

Set on `understory.py` (and surfaced through `run_pipeline.py`):

- `--shrub-canopy-cell-size` (default 0.25), `--shrub-stem-cell-size` (0.25) —
  raster resolution for the shrub canopy/DBH maps.
- `--no-shrub-outputs` — skip per-shrub outputs entirely.

Shrub volume is calculated from the sampled shape, never from voxels
(notes.md ruling, 2026-11-07); the earlier occupied-voxel estimate and its
`--occupied-voxel-*` flags were removed.

## Notes / caveats

- `mesh` volume falls back to the convex hull of the mesh vertices when
  trimesh's own `.volume` isn't available (e.g. a non-watertight mesh); the hull
  slightly overestimates relative to the true jittered surface.
- The stub trimesh used in the offline dev environment has no `.volume`, so the
  hull fallback is what runs there; with real trimesh installed the watertight
  mesh volume is used directly.
