# Canopy height map + 3D canopy geometry

Stage 5: a canopy height model (CHM) raster of the placed plot, plus per-tree 3D
canopy bounding boxes, convex hulls, crops, masks, and a total canopy volume.

```
python canopy_map.py plot_11519.csv --chm-resolution 0.01
```

```
plot_11519/
  plot_placement/
    plot_11519_raw_placement_canopy_map.tif                 CHM (GeoTIFF, float32)
    plot_11519_raw_placement_canopy_bounding_box.json       per-tree 3D bbox
    plot_11519_raw_placement_canopy_convex_hull.json        per-tree 3D hull
    plot_11519_scene.json                                   gains a "canopy_map" block
  individual_trees/
    plot_11519_tree_<UID>_canopy_bounding_box.tif           CHM cropped to bbox+buffer
    plot_11519_tree_<UID>_canopy_semantic_mask.tif          canopy mask, same size
```

## The CHM raster

A 2-band top-down (x-y) raster. Band 1's value in each cell is the **height of
the highest point** falling in it, and **0** where no point falls (no tree).
Band 2 is the canopy mask (see below). Built from all
points of the placed scene `plot_placement/<PLOT_UID>_raw.laz`, so it is the top
surface of the whole plot. Since every tree sits with its base at z = 0 and
placement only shifts x-y, the value is height above ground directly. Written as
a 2-band float32 GeoTIFF (DEFLATE-compressed, geotransform embedded, same
writer as the stem map). Default cell size 0.05 m (`--chm-resolution`). Each band-1
cell holds the max z of the points in it, and the raster max equals the cloud's
max z.

## 3D canopy bounding box

`_canopy_bounding_box.json` gives, per tree with branch points, the axis-aligned
rectangular prism enclosing that tree's **branch points** (`semantic == 0`):
`bbox_world` with x/y/z min/max and dx/dy/dz, `bbox_volume_m3`, and the buffered
x-y pixel box used for the crop. World min/max are rounded **outward** (floor
the mins, ceil the maxes) so the stored box always encloses every branch point
despite 4-dp rounding.

## 3D convex hull

`_canopy_convex_hull.json` gives, per tree, the convex hull of its branch
points: `vertices` (hull points), `faces` (triangles indexing into vertices),
`volume_m3`, `surface_area_m2`. Faces are **wound outward** (normals agree with
the hull's outward equations), so a naive signed-volume sum over the stored mesh
reproduces the reported volume. A tree with fewer than 4 branch points,
or a coplanar set, is marked `degenerate: true` with volume 0 and empty
vertices/faces; none occur in the example.

## Per-tree crops and masks

For each tree with a canopy:

- `_canopy_bounding_box.tif` — the CHM cropped to that tree's canopy **x-y bbox
  plus `--buffer` pixels** (default 2), clipped to the raster.
- `_canopy_semantic_mask.tif` — a binary (0/1) raster the **same size** as the
  crop: the tree's **filled canopy footprint** — the 2D convex hull of its
  foliage points with the interior filled (item W29), unioned with the point
  cells so a hull-boundary point is never lost — so cells inside the crown
  are 1 even where no foliage point falls. A foliage set too small or
  collinear to hull contributes its point cells only, so a thin canopy never
  vanishes.

The plot-level CHM band 2 is the same filled footprint unioned over all trees,
and the plot-level instance label (`canopy_mask_raw.tif`) fills each tree's
footprint with its `tree_id`, giving the **taller crown** the cell where
footprints overlap. Like the stem map's masks, these are deliberately not
clipped to the feature raster's point support — the canopy interior is part of
the mask, so band 2 may mark cells where band 1 is 0.

## Total canopy volume

`_scene.json`'s `canopy_map.total_canopy_volume` sums the per-tree convex hull
volumes over the successfully placed trees — a single post-trim value
(`volume_m3`); there is no raw/trimmed pair (spec 3.2: no record of removed
points is kept).

Convex hull volume (not bounding-box volume) is used as the canopy volume, since
the hull is the tighter solid the task defines around the branch points.

## Two things to be aware of

**Volume is resolution-independent.** It comes from the 3D points and hulls, not
the raster, so it does not change with `--chm-resolution` (confirmed 4786.0 m³ at 1 cm
and 5 cm). Only the raster and pixel bboxes scale with cell size.

**Convex hulls overestimate real canopy volume.** A hull wraps the outermost
foliage points into one solid, filling gaps and concavities a real crown has. It
is the right answer to "the convex hull that contains the branch points," but do
not read it as occupied leaf volume; `build_info.py`'s occupied-voxel volume is
the tighter estimate.
