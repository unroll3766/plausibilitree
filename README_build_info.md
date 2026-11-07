# Info JSONs (`build_info.py`)

`build_info.py` emits the plot- and tree-level `*_info.json` files, replacing the
old `*_volume.json` files. It runs last in the pipeline (after every cloud,
skeleton, and canopy map exists) and only reads, so its order is flexible as long
as those inputs are present. Shrub `*_volume.json` files are **not** touched.

```
plot_placement/textual/<plot>_info.json
individual_trees/<plot>_<uid>/textual/<plot>_<uid>_info.json   (one per tree)
```

## Plot-level schema

```
meta: plot_uid, coordinate_system, units, plot_shape, plot_center, bounding_box
vegetation_characteristics:
  trees:
    n_trees, basal_area_m2
    total_volume_m3                        = trunk + branch + foliage
    trunk_volume_m3
    branch_volume_m3
    foliage_volume_m3
    canopy_volume_m3
  shrubs:
    n_shrubs, total_volume_m3, primitive_reference_m3
  grass:
    percent_cover, total_volume_m3
```

`meta.plot_center` is `{lat, lon, elev_m}`, read from the scene JSON's
`plot_lat`/`plot_lon`/`plot_elev` (see README_assembly.md's "Scene metadata"),
which are themselves a straight pass-through of the FIA input's LAT/LON/ELEV
by `fia_helper.py` (README_fia_helper.md). All three are `null` if the FIA
input didn't carry them. This is plot-center metadata, not per-tree — tree-
and shrub-level info JSONs don't repeat it.

## Tree-level schema

```
meta: tree_uid, plot_uid, coordinate_system, units, bounding_box
tree_characteristics:
  center [x, y], dbh_m, height_m
  canopy:
    convex_hull_volume_m3
    major_axis {length_m, orientation_deg}
    minor_axis {length_m, orientation_deg}
  volume_m3:
    trunk
    branch
    foliage
neighborhood:
  n_other_trees_in_bbox, n_shrubs_in_bbox
  target_voxel_volume_m3
  proportion_target_vs_other_trees
  proportion_target_vs_shrubs
  proportion_target_vs_grass
  proportion_target_vs_all_other_veg
```

## How each value is derived

**Trunk volume** — frustum (truncated-cone) sum along the skeleton trunk
centerline. For each pair of consecutive polyline nodes the radius comes from the
taper model `r(h) = base_radius * (1 - h / height)`, and the frustum volume is
`(pi/3) L (r0^2 + r0 r1 + r1^2)` with `L` the 3D segment length. The (placed)
skeleton carries only the surviving geometry (spec 3.2).

**Branch volume** — each branch is a **cylinder**: `pi * r^2 * length`, where `r`
is the taper radius at the branch's base height; trimmed branches are physically
absent from the skeleton.

**Foliage volume** — occupied-voxel volume of the foliage points
(`semantic_code == 3`): voxelize at `occupied_voxel_size`, count voxels holding
`>= occupied_voxel_threshold` points, times the voxel cell volume, over the
surviving (trimmed) cloud.

**Canopy volume** — convex-hull volume of the foliage points, read from the
per-tree `canopy_convex_hull.json` (single post-trim value).

**Total volume (plot trees)** — `trunk + branch + foliage` (canopy is a separate
envelope measure and is reported alongside, not summed into the total).

**Basal area** — `sum over trees of pi * (dbh/2)^2`, in m^2 per plot. DBH is the
taper diameter at breast height (1.37 m).

**Crown major / minor axis** — computed on the crown's convex hull in x-y (foliage
points, trimmed preferred). The **major axis** is the longest chord between hull
vertices; the **minor axis** is the widest extent of the hull measured
perpendicular to the major direction. Lengths in metres; orientation in degrees
CCW from +x in the local frame (minor = major +/- 90 deg by construction).

**Shrub total volume** — sums each of the three per-shrub methods (primitive,
mesh, voxel) from the shrub index. NaN mesh values (non-watertight meshes) are
skipped.

**Grass** — `total_volume_m3` is the occupied-voxel volume of grass points.
`percent_cover` is the **realized** cover: the fraction of plot footprint cells
(0.5 m) that contain any grass point. This is measured from the cloud, not the
input `grass_cover` parameter, because the generator uses that parameter only as
a thinning weight (modulated by the vigor field and canopy attenuation), so the
realized cover differs from the input.

**Neighborhood** — within the target tree's bounding box: counts of other trees
and shrubs, and volume proportions `target / (target + other)`. All volumes here
(target and every neighbour class) are occupied-voxel volumes at the same voxel
params, so the proportions are on one consistent basis. `target_voxel_volume_m3`
is the target tree's own occupied-voxel volume inside its bbox, exposed so the
proportions are interpretable. A proportion is `null` when both volumes are zero.

## Post-trim values

Every volume is a single post-trim value (spec 3.2): trimmed elements are
physically absent from the (placed) skeleton, and the clouds are the
surviving points — so there are no raw/trimmed pairs and no record of
removed material anywhere.

## Parameters

Set on `build_info.py` (and surfaced through `run_pipeline.py`):

- `--occupied-voxel-size` (default 0.1) — voxel edge (m) for foliage/grass/
  neighborhood voxel volumes.
- `--occupied-voxel-threshold` (default 1) — min points for an occupied voxel.
- `--coord-system` (default `local`), `--units` (default `meters`) — recorded in
  `meta`.

## Notes / caveats

- The offline dev environment stubs trimesh, so shrub `mesh` volumes in that
  environment are not the true watertight volumes; with real trimesh installed
  they are correct. Trunk/branch/foliage/canopy/axis math uses real numpy/scipy
  and is exact.
- Coordinate system is `local` (plot metres, z = 0 at ground) unless a real EPSG
  is supplied; there is no north in the local frame, so axis orientation is
  measured from +x.
