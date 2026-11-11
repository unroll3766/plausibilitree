# Info JSONs (`build_info.py`)

`build_info.py` emits the plot- and tree-level `*_info.json` files, replacing the
old tree `*_volume.json` files. It runs last in the pipeline (after every
cloud, skeleton, and canopy map exists) and only reads, so its order is
flexible as long as those inputs are present. Per-shrub values come from the
plot shrub index (`understory.py` writes it); there are no per-shrub
`*_volume.json` files anymore (removed in the 2026-11 review).

```
plot_placement/textual/<plot>_info.json
individual_trees/<plot>_<uid>/textual/<plot>_<uid>_info.json   (one per tree)
```

## Plot-level schema

notes.md field names at the top level (spec 8.4), extension blocks below them:

```
NUM_TREES                      int
BASAL_AREA                     float, m^2 (sum of pi * (TREE_DBH/2)^2)
NUM_SHRUBS                     int
VOL_TREES                      float, sum of every tree's VOL_TREE
VOL_TRUNK  VOL_BRANCH  VOL_FOLIAGE   float, sums over tree info JSONs
VOL_SHRUB                      float, sum of per-shrub mesh volumes
VOL_GRASS                      float, occupied-voxel at GRASS_VOXEL
PLOT_SHAPE                     "circle" | "rectangle"
PLOT_AREA                      float, m^2 (pi*r^2 or length*width)
BOUNDING_BOX                   [x_min, y_min, z_min, x_max, y_max, z_max]
CHM_PATH                       relative path, null if the stage did not run
STEM_PATH                      relative path, null if the stage did not run
PLAUSIBILITREE_CONFIGURATION   relative path to the scene JSON (settings record)
PLOT_CENTER                    optional [lon, lat, elev]; present only when the
                               FIA input carried the pass-through columns
DTM_PATH                       optional relative path; present only when a plot
                               DTM was generated
meta: plot_uid, coordinate_system, units
shrubs: primitive_reference_m3        (analytic reference; VOL_SHRUB is the mesh)
grass: percent_cover                  (realized; measured from the cloud)
```

`PLOT_CENTER` is read from the scene JSON's `plot_lat`/`plot_lon`/`plot_elev`
(see README_assembly.md), a straight pass-through of the FIA input's
LAT/LON/ELEV by `fia_helper.py`. notes.md's spelling of the configuration
field is "PLAUSIBILITREE_CONFIGURAITION"; the corrected spelling is used
(owner-approved).

## Tree-level schema

```
TREE_UID    string              TREE_HT   float, m
TREE_DBH    float, m (taper)    TREE_X    float, m (base, matches tree list)
TREE_Z      float, m            TREE_Y    float, m (base, matches tree list)
VOL_TREE    = VOL_TRUNK + VOL_BRANCH + VOL_FOLIAGE
VOL_TRUNK   VOL_BRANCH   VOL_FOLIAGE    float, m^3
BOUNDING_BOX                  [x_min, y_min, z_min, x_max, y_max, z_max]
meta: plot_uid, coordinate_system, units
canopy: major_axis, minor_axis        (convex-hull volume removed in the
                                      2026-11 review: it was only populated
                                      when the optional canopy-map stage ran)
neighborhood: (unchanged; see below)
```

`TREE_Z` is the centroid height above the tree's base: the mean
height-above-ground of the tree's points in the understory cloud, falling
back to the centroid height above the lowest point of the tree cloud when
no understory cloud exists. When the plot is anchored (scene `geo` block,
spec Section 10): BOUNDING_BOX becomes `[lon_min, lat_min, elev_min,
lon_max, lat_max, elev_max]` (elev = plot-center height + local z, EPSG
4326/5773) in both info levels, and each tree info gains the optional
`TREE_CENTER_REAL [lon, lat, elev]` at its base (the bottom of the tree
is the elevation, per notes.md). Unanchored plots keep local metres.

`canopy.major_axis` / `canopy.minor_axis` are `{length_m, orientation_deg}`
(crown hull in x-y; major = longest chord, minor = widest span
perpendicular to it).

## Shrub-level schema

`meta` (shrub_uid / plot_uid / coordinate_system / units / bounding_box) +
`shrub_characteristics` (center, max_width_m, height_of_max_width_m,
height_m, volume_m3.primitive, volume_m3.mesh) + `neighborhood`. The
2026-11 review removed the ellipsoid `dbh_m` field and the always-null
`volume_m3.voxel` key (shrub volume is never voxel-based per the notes.md
ruling).

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
(`semantic_code == 3`): voxelize at `--foliage-voxel` (FOLIAGE_VOXEL), count voxels holding
`>= --voxel-threshold` points, times the voxel cell volume, over the
surviving (trimmed) cloud.

The canopy convex-hull volume extension was removed in the 2026-11 review
(it was only populated when the optional canopy-map stage ran); the crown
axes remain.

**Total volume (VOL_TREE / VOL_TREES)** — `trunk + branch + foliage`.

**Basal area** — `sum over trees of pi * (dbh/2)^2`, in m^2 per plot. DBH is the
taper diameter at breast height (1.37 m).

**Crown major / minor axis** — computed on the crown's convex hull in x-y (foliage
points, trimmed preferred). The **major axis** is the longest chord between hull
vertices; the **minor axis** is the widest extent of the hull measured
perpendicular to the major direction. Lengths in metres; orientation in degrees
CCW from +x in the local frame (minor = major +/- 90 deg by construction).

**Shrub total volume (VOL_SHRUB)** — sums the `volume_m3.mesh` values from
the per-shrub info JSONs (sampled shape per the notes.md ruling; never
voxels). The nominal primitive (ellipsoid) volumes are the
`shrubs.primitive_reference_m3` extension.

**Grass (VOL_GRASS)** — occupied-voxel volume of grass points at
`--grass-voxel` (GRASS_VOXEL); `grass.percent_cover` is the **realized** cover: the fraction of plot footprint cells
(0.5 m) that contain any grass point. This is measured from the cloud, not the
input `grass_cover` parameter, because the generator uses that parameter only as
a thinning weight (modulated by the vigor field and canopy attenuation), so the
realized cover differs from the input.

**Neighborhood** — within the target tree's bounding box: counts of other trees
and shrubs, and volume proportions `target / (target + other)`. All volumes here
(target and every neighbour class) are occupied-voxel volumes at
FOLIAGE_VOXEL / the shared threshold, so the proportions are on one
consistent basis. `target_voxel_volume_m3`
is the target tree's own occupied-voxel volume inside its bbox, exposed so the
proportions are interpretable. A proportion is `null` when both volumes are zero.

## Post-trim values

Every volume is a single post-trim value (spec 3.2): trimmed elements are
physically absent from the (placed) skeleton, and the clouds are the
surviving points — so there are no raw/trimmed pairs and no record of
removed material anywhere.

## Parameters

Set on `build_info.py` (and surfaced through `run_pipeline.py`):

- `--foliage-voxel` (default 0.1) — FOLIAGE_VOXEL: voxel edge (m) for
  foliage and neighborhood voxel volumes.
- `--grass-voxel` (default 0.1) — GRASS_VOXEL: voxel edge (m) for the
  grass voxel volume.
- `--voxel-threshold` (default 1) — min points for an occupied voxel.
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
