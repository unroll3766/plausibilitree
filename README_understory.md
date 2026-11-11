# Understory (shrubs + grass) — `understory.py`

Optional stage that adds procedural **shrubs** and **grass** to a plot's terrain
cloud. Runs **after `terrain.py`**, as the final cloud-writing stage (it needs the draped
terrain cloud and the DTM).

```
python understory.py plots/plot_11519.csv --seed 42 \
    --shrub-cover 0.15 --grass-cover 0.6 --grass-candidate-density 200
```

Writes:

```
plot_placement/las/<PLOT_UID>_understory.laz
individual_trees/<PLOT_UID_TREE_UID>/las/<..>_understory.laz
```

Implemented from `vegetation_design_sketch.py`; the two generators live in
`vegetation.py` (`gen_shrub_layer`, `gen_grass_layer`, plus the primitives
`gen_shrub_ellipsoid_mesh` and `gen_vigor_field`).

## Shrubs

Discrete jittered-ellipsoid shrubs, one mesh each. Grid the plot at
`--shrub-cell-size`; exclude cells within `canopy_r + --shrub-margin` of any
tree (cKDTree); keep the rest with probability `--shrub-cover`; per kept cell
jitter a centre, sample height in `[--shrub-height-min, --shrub-height-max]`,
width in `[--shrub-width-min, --shrub-width-max]`, build a warped icosphere,
sample `--points-per-shrub` points, drape onto the DTM. Because cells near
canopies are excluded up front, shrubs sit in the gaps between trees.

## Grass

Continuous point field via Poisson thinning of a diamond-square **vigor field**
`D(x,y) ∈ [0,1]` (reusing the DTM generator, so no second noise implementation).
Scatter `--grass-candidate-density` candidates/m², accept each with probability
`cover · D · attenuation`, where attenuation `= --canopy-attenuation ** (#
overlapping tree/shrub footprints)` — under both a tree and a shrub it's applied
twice (multiplicative, per the sketch's first option). Height interpolates in
`[--grass-height-min, --grass-height-max]` by `D`; `z = ground + height · u²`
(u²-biased toward the blade base). No per-plant geometry.

## Schema

The plot understory cloud is the terrain cloud (trees + ground, unchanged) with
shrub and grass points appended, in the canonical schema (README_schema):
`tree_id`, `shrub_id`, `semantic_code`, `classification`, `branch_code`,
`height_above_ground`.

- **`semantic_code`** (object identity): shrubs → **4**, grass → **5**, on top
  of the existing tree (1 trunk, 2 branch, 3 foliage) and ground (6) codes.
- **LAS `classification`** stays orthogonal, set by height above ground:
  ≤0.5 m → 3 (low veg), 0.5–2 m → 4 (med veg), >2 m → 5; ground stays 2. So a
  shrub is classed by its height like a real return, independent of its identity.
- **`shrub_id` / `tree_id`**: trees keep `tree_id`, `shrub_id=0`; shrubs get
  `tree_id=0`, `shrub_id≥1`; grass and ground get both `=0`.
- **`branch_code`** is 0 on shrub/grass/ground (they're not tree parts).

The ID rules hold exactly; classification matches the HAG bands (a shrub spans
classes 3–5 by height); every terrain point carries through unchanged.

Per-tree `_understory.laz` crops contain all points inside that tree's canopy
x-y bbox — so a crop includes nearby shrubs/grass and the ground patch, not just
the tree.

## In the pipeline

`run_pipeline.py --understory` enables it (requires terrain via `--dtm` or
`--generate-dtm`; otherwise it's skipped with a note). Runs per plot after
terrain, as the final cloud-writing stage. The normalization stage is retired
(2026-11 review): the DL side owns it. Knobs: `--shrub-cover`,
`--shrub-cell-size`, `--grass-cover`, `--grass-candidate-density`,
`--canopy-attenuation`.
