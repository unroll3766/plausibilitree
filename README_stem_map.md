# Stem map raster

Stage 4: rasterise a breast-height cross-section of the placed plot into a stem
map for deep learning, with per-tree bounding boxes and semantic masks.

```
python stem_map.py plot_11519.csv --stem-resolution 0.01
```

```
plot_11519/
  plot_placement/
    plot_11519_raw_placement_stem_map.tif                  count raster (GeoTIFF)
    plot_11519_raw_placement_stem_map_bounding_boxes.json  per-tree stem bboxes
    plot_11519_scene.json                                  gains a "stem_map" link
  individual_trees/
    plot_11519_tree_<UID>_stem_map_bounding_box.tif        raster cropped to bbox
    plot_11519_tree_<UID>_stem_map_semantic_mask.tif       stem mask, same size
```

## What the raster is

A top-down (x-y) image of the point cloud in a thin horizontal band at DBH
(breast height, 1.37 m) ± a threshold. Band 1 holds the **count of all
points** falling in each cell. Band 2 is the stem **mask**: 1 only where a
**branch-order-0** tree point (`branch_code == 1`, i.e. a trunk-stem point on
an order-0 element per spec Section 7.4) lies in the slab — never ground,
shrub/grass, foliage or higher-order branch-stem returns. Built from the
placed scene `plot_placement/<PLOT_UID>_raw.laz` (use `--source trimmed` for
the trimmed cloud instead; the raw plot cloud carries `branch_code` from
assembly, so the order-0 rule applies on either source).

Two mask rules (spec Section 5, item W28):

- **Order 0 only.** Only branch-order-0 points enter the mask. The
  skeleton-derived masks (per-tree semantic masks, plot instance label) use
  only order-0 trunk polylines — a fork daughter with order ≥ 1 never
  contributes, even where it crosses the band.
- **Below-band forks excluded.** A tree with any fork below the DBH range
  (split height < DBH_HT − STEM_SECTION, here < 1.12 m) is excluded from the
  stem map mask entirely: no band-2 cells, no per-tree crop/mask. It is
  listed under `excluded_trees` in the bounding-box JSON (and the scene
  block) with its lowest below-band fork height, so the exclusion is
  auditable rather than silent.

Because every tree is generated with its base at z = 0 and placement only
shifts x-y, breast height is a single global z-band across the whole plot — no
terrain to account for.

Defaults: cell size 0.05 m (`--stem-resolution`), band half-thickness 0.25 m
so the slice is 0.5 m thick (`--stem-section`), DBH 1.37 m (`--dbh-height`).
The raster
value type is uint16 (counts stay small — max 5 per cell in the example).

## Bounding boxes

`_bounding_boxes.json` lists, per tree that has an order-0 stem in the band
(plus `n_trees_excluded` / `excluded_trees` for below-band-fork trees):

- `tree_id`, `n_stem_points`, `padding_px`
- `bbox_pixel` — the **padded** box (inclusive col/row min/max plus
  width/height) in the stem map's pixel grid. This is the box used for the crop
  and mask.
- `bbox_pixel_tight` — the tight box around the stem points, before padding.
- `bbox_world` — the padded box in metres (cell corners).

Boxes come from the **order-0 stem cross-section** (skeleton discs), while the
raster's band 1 counts **all** band points, so a box is the stem's footprint
inside the fuller image.

### Padding

Each tree box is grown by `--bbox-padding` pixels on every side (default 1), so
a tight 5×10 stem box becomes a 7×12 crop. Padding is **clipped to the raster**:
a tree whose stem touches a raster edge gets no padding on that side (in the
example, three trees sit on the grid boundary and pad only on their interior
sides). `--bbox-padding 0` reproduces the tight box exactly. The tight box is
always kept in `bbox_pixel_tight`, so padding is reversible.

The file also records `raster_origin_world` (rounded, for reading) and
`raster_origin_world_exact` (full precision). Use the exact origin to
reconstruct pixels from world coordinates: `col = floor((x - x_min)/cell)`,
`row = floor((y_max - y)/cell)`, then clip to `[0, width/height - 1]`. The clip
matters — this is a `PixelIsArea` grid, so a point exactly on the top or right
edge (which can be the point that *defined* the extent) lands one index past the
grid before clipping.

## Per-tree crops and masks

For each non-excluded tree with an order-0 stem in the band:

- `_stem_map_bounding_box.tif` — the stem map cropped to that tree's **padded**
  pixel bbox.
- `_stem_map_semantic_mask.tif` — a binary (0/1) raster the **same size** as the
  crop: the filled order-0 stem cross-section (skeleton discs), 1 inside the
  stem disc. With padding, the stem cells sit inside an empty border the width
  of the padding.

The two align pixel-for-pixel, so they're a ready image/label pair: every mask
cell has a positive count in the matching crop.

## Scene link

`_scene.json` gains a `stem_map` block with relative links in the requested
form (`<plot>/plot_placement/<plot>_raw_placement_stem_map.tif` etc.), the band
z-range, cell size, source cloud and raster size.

## GeoTIFF writing

Written with `rasterio` when installed, otherwise with `tifffile`, embedding the
GeoTIFF tags (ModelPixelScale, ModelTiepoint, a minimal GeoKeyDirectory) so the
transform is preserved either way. Rasters use DEFLATE compression (a standard
GeoTIFF codec, no extra package needed), which collapses the mostly-empty
full-plot raster from ~55 MB to ~65 KB; if the codec is somehow unavailable the
writer falls back to uncompressed rather than emitting an unreadable file.
There is no CRS — coordinates are plot-local metres — so the projected-CS key
is left user-defined; readers get the pixel size and origin but no datum, which
is correct for a synthetic local frame.

## Two things to be aware of

**Thin bands lose trees.** At the 0.10 m default the example keeps 44 of 49
trees; at `--stem-section 0.03` only 22 survive, because a 6 cm slice of a sparse
synthetic cloud simply contains no stem points for many trees. Thinner slices
give rounder stem discs but fewer usable stems — pick per how dense your clouds
are.

**Below-band forks are excluded, not boxed.** A tree that has forked below
the DBH range (fork height < DBH_HT − STEM_SECTION) gets no box, crop or
mask — it is recorded in `excluded_trees` instead, because a cross-section
above a low fork is not a measurable DBH stem. A fork *within* the band keeps
its box, but only the order-0 daughter's cross-section enters the mask, so
the box is the dominant stem's, not the gap between both apical stems.
