# Output restructuring (canonical layout)

Stage 6: reorganise the flat stage outputs into the canonical directory
tree, and generate the derived artifacts that tree implies.

```
# after stages 2-5:
python restructure.py plot_11519.csv
```

Two modules:

- **`layout.py`** — the single source of truth for every output path. Stage
  scripts and the restructure step both address files through it, so the layout
  is defined in exactly one place.
- **`restructure.py`** — moves the flat files into the layout and produces the
  new derived artifacts.

## The layout

```
PLOT_UID/
  individual_trees/
    PLOT_UID_TREE_UID/
      textual/        generation.json, raw_skeleton.json, trimmed_skeleton.json
      las/            raw.laz, trimmed.laz
      deep_learning/
        volume/labels/            volume.json
        stem_map/features/        stem_bounding_box_raw.tif
        stem_map/labels/          stem_mask_raw.tif
        canopy_map/features/      canopy_bounding_box_raw.tif
        canopy_map/labels/        canopy_mask_raw.tif
                                  canopy_bounding_box.json, canopy_convex_hull.json
  plot_placement/
    textual/          manifest.csv, PLOT_UID_scene.json
    las/              raw.laz, trimmed.laz
    deep_learning/
      volume/labels/            PLOT_UID_volume.json
      stem_map/features/        PLOT_UID_stem_map_raw.tif
      stem_map/labels/          PLOT_UID_stem_mask_raw.tif
      canopy_map/features/      PLOT_UID_canopy_map_raw.tif
      canopy_map/labels/        PLOT_UID_canopy_mask_raw.tif
```

`TREE_UID` files are prefixed with `PLOT_UID_TREE_UID` (e.g.
`plot_11519_003_raw.laz`); plot files with `PLOT_UID`.

## What moves vs. what's derived

**Moved** (renamed into place, not copied): the per-tree generation/skeleton
JSONs, raw/trimmed clouds, and the stem/canopy crop + mask TIFFs; the plot-level
clouds, manifest, scene JSON, and the stem-map and canopy-map feature rasters.

**Derived here** (new artifacts the layout calls for that the stages didn't
emit):

- **Per-tree `volume.json`** — that tree's convex-hull canopy volume, `raw` and
  (if the trimmed cloud exists) `trimmed`, plus its bbox volume. Recomputed from
  the branch points.
- **Per-tree `canopy_bounding_box.json` / `canopy_convex_hull.json`** — split
  out of the plot-level combined files the canopy stage produced (which have no
  home in the new layout and are removed after splitting).
- **Plot-level `volume.json`** — total canopy volume summed over trees, `raw`
  and `trimmed` (equal to the sum of the per-tree files).
- **Plot-level `stem_mask_raw.tif` / `canopy_mask_raw.tif`** — the label rasters
  that pair with the existing feature rasters (see below).

The plot-level stem bounding-box list (no dedicated file in the new layout) is
folded into the scene JSON under `stem_map.bounding_boxes_inline`.

## Plot-level masks (features ↔ labels)

Each `deep_learning/*/features/*_raw.tif` (a continuous raster) now has a
matching `labels/*_mask_raw.tif`. The labels are **instance masks keyed by
`tree_id`** (uint32, 0 = background), written to the exact grid of their feature
raster (origin, cell size and shape read back from the feature's GeoTIFF tags),
so they align pixel-for-pixel. `(mask > 0)` gives a binary semantic mask.

- **Stem label** — `tree_id` of the stem points in each breast-height-band cell.
  Collisions (rare, since stems are separated) resolve to the lowest `tree_id`.
- **Canopy label** — `tree_id` of the **highest** point in each cell, matching
  the CHM feature, which stores that same point's height: labelled cells are
  exactly the feature's occupied cells.

Note the plot-level canopy label differs from the **per-tree** canopy masks: the
per-tree masks mark *all* of a tree's branch cells (including understory ones
hidden under a taller crown, as specified earlier), whereas the plot-level label
can only name one owner per cell and names the visible top. Both are intentional
for their respective purposes.

A subtlety worth recording: the CHM feature is built with a zero-initialised
max, so a cell whose only points sit at z ≤ 0 (possible from the near-ground
foliage z-jitter) reads as empty. The label mirrors that — cells the feature
treats as empty are cleared — so feature and label never disagree. Without this
guard the canopy label had 694 phantom cells over ground-grazing points.

## Scene JSON as the index

The rewritten `PLOT_UID_scene.json` keeps all the assembly/stem/canopy metadata
and adds:

- `plot_files` — relative links to every plot-level artifact.
- `trees` — per tree, relative links to each artifact **that exists** for it.
  Trees with no stem in the DBH band (5 of 49 here) simply omit their stem
  links rather than pointing at missing files; all 439 links in the example
  resolve.

Links are forward-slash relative paths from the collection root (the parent of
the plot directory), e.g. `plot_11519/plot_placement/las/plot_11519_raw.laz`.

## Re-running

`restructure.py` moves files, so it is a one-shot per plot. Re-running detects
that the flat placement cloud is already gone and skips cleanly (prints
"already restructured") rather than erroring, so it is safe in a re-run of a
batch pipeline.

## Design note

The stage scripts still write the flat layout and `restructure.py` packages it,
rather than each stage writing into the tree directly. This kept five tested
scripts untouched. Because all paths are centralised in `layout.py`, migrating
the stages to emit the canonical layout directly (removing the move step) is a
later change confined to swapping their output-path calls — the layout itself
would not move.
