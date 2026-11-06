# Individual tree generation

Stage 2: turn a PlausibiliTree plot table into one labelled point cloud per
tree, with full provenance.

```
python generate_trees.py plot_11519.json --seed 42   # canonical tree list (default)
python generate_trees.py plot_11519.csv  --seed 42   # legacy SST parameter table
```

```
plot_11519/
  manifest.csv
  individual_trees/
    001_raw.laz
    001_generation.json
    001_raw_skeleton.json
    ...
```

The output directory is derived from the table name (`plot_001.json` or
`plot_001.csv` →
`plot_001/`) unless `-o` is given. Rows with `successful_placement=False` are
dropped before generation.

## Point budget

`Np = points_per_metre * ceil(HT_HI)`, default 1024, so points per vertical
metre is constant across the plot rather than points per tree. A 4 m sapling
gets 5,120 points and a 40 m dominant gets 40,960. Override with
`--points-per-metre`.

Note this normalises by *height only* — a wide-crowned oak and a narrow fir of
the same height get the same budget, so the oak's surface density is lower.

## Reproducibility

Every tree gets a deterministic seed derived from the run seed and its
`TREE_UID` (`crc32`, not `hash()`, which is salted per process). Re-running
with the same `--seed` reproduces the plot exactly, structures and surface
point placements alike; adding or removing trees does not disturb the
others. trimesh >= 5 samples mesh surfaces with its own unseeded generator
that ignores `np.random.seed`, so `treesim.py` installs a wrapper that
derives the sampler seed from the global RNG state without consuming a draw
-- without it the same seed would build the same skeleton but scatter
different surface points on every run. The seed
actually used is recorded in each `_generation.json`, so any single tree can
be regenerated in isolation.

`_generation.json` carries the full `model_params` dict, the point budget, the
realised geometry, the source plot-table row, and the seed — enough to rebuild
the tree without the plot table.

## Skeleton format

`_raw_skeleton.json` records the underlying structure, before point sampling:

```json
{
  "tree_uid": "015", "plot_uid": "plot_11519", "tree_id": 15,
  "coordinate_frame": "plot metres, z = 0 at tree base",
  "base_xy": [-3.8434, 12.2217],
  "height_m": 16.72, "base_radius_m": 0.324,
  "fork": 2, "fork_height_m": 3.92,
  "min_canopy_height_m": 3.2, "max_canopy_width_height_m": 7.96,
  "trunk": {
    "total_length_m": 21.4, "n_segments": 3,
    "segments": [{"z_range": [0.0, 3.92], "length_m": 3.93}],
    "polylines": [[[x, y, z], "... 20 vertices per segment"]]
  },
  "branches": [{"base": [x, y, z], "tip": [x, y, z],
                "height_m": 4.88, "azimuth_rad": 0.76, "length_m": 3.26}],
  "n_branches": 197
}
```

Changes from the prior skeleton snippet: kept the trunk polyline/arc-length
logic essentially as-is (it reads straight off `stem_data`, which is exactly
the right source), dropped the underscore-prefixed names since these are now
module-level rather than scratch variables, and added the branch axes. Branch
segments are cheap to record here and expensive to recover from the point
cloud later, which is the whole argument for writing a skeleton at all —
it gives you noise-free ground truth for branch detection, DBH-from-cloud, and
QSM comparisons.

A forked tree has 2–3 trunk polylines; `fork_height_m` is `null` when
`fork == 1`. Note that `n_branches` can exceed `NUM_BRANCH_HI` for forked
trees — `treesim` scales the upper-canopy branch count by the fork factor.
That is upstream behaviour, left unchanged.

## Changes to treesim.py

Patched rather than forked, so there is one implementation:

1. `gen_single_simtree(..., return_skeleton=False)` — when true, returns
   `(points, skeleton)` instead of `points`.
2. Branch axes are recorded in the branch loop (base = first vertex, tip =
   last vertex of each 6-vertex branch strip).
3. `gen_simtree(..., tx=None, ty=None, return_skeleton=False)` — passes the
   placement offset through to `gen_single_simtree`, which already accepted it
   but had no way to be told.

Default behaviour is unchanged; `gen_simtree(Np=4096)` still returns a plain
`(N, 4)` array.

## LAS output

Each cloud is written via the project's `las.py` as `(x, y, z, semantic,
tree_id)`, semantic 1 = stem, 0 = foliage. `tree_id` is the integer form of
`TREE_UID`, so per-tree files concatenate into a plot-level cloud without
further bookkeeping. Coordinates are already in plot metres — `X`/`Y` from the
table are applied during generation, not afterwards.

## Failure handling

| status | meaning |
|---|---|
| `ok` | written |
| `invalid_params` | non-finite or non-positive height/diameter/crown width, or inverted LO/HI — caught before simulating |
| `failed` | simulator raised (or returned non-finite points) on every attempt |
| `skipped_exists` | `--no-overwrite` and output already present |

Degenerate randomisations do occur — if every branch happens to land in the
azimuth gap the foliage mesh is empty and `trimesh` raises — so each tree gets
up to `--max-attempts` reseeded tries before being recorded as failed. The
run continues either way, and `manifest.csv` is the record of what happened.
