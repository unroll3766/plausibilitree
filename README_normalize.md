# Normalization (optional DL preprocessing)

Stage 7: optional normalizations of the placed clouds, writing the layout's
`*_normalized.laz` slots. Users pick any combination of three steps, always
applied in this fixed order.

```
# after restructure.py:
python normalize.py plot_11519.csv --max-density 5000 --spatial \
    --density cube --plot-l 64 --tree-l 32
```

> **Run `restructure.py` first.** `normalize.py` reads its source clouds from the
> canonical layout (`plot_placement/las/*_trimmed.laz` and the per-tree `las/`
> dirs), which only exist after restructuring. Run it on the flat outputs and it
> raises a clear error naming the cause. Correct order: generate → assemble →
> stem → canopy → **restructure** → normalize.

## The three steps (in order)

1. **`--max-density N`** — hard cap of N points per 1 m³. A fixed 1 m grid is
   keyed off `floor(coord)`; any voxel over N has the excess randomly dropped
   (survivors keep every column).

2. **`--spatial`** — per-axis rescale of x, y, z each into `[-1, 1]`. Aspect
   ratio is **not** preserved (a tall thin tree fills the cube). Applied after
   the cap, before density.

3. **`--density square|cube`** — replaces the cloud with a fixed lattice of
   points, one per cell/voxel centre, **including empty cells**. Terminal step.

Square and cube are mutually exclusive and are the shared base shape across the
plot and tree levels; only `L` (and `Z`) may differ per level (`--plot-l/-z`,
`--tree-l/-z`). If a step is omitted the others still run; with no density step
the (capped/rescaled) point cloud is written as-is with columns
`x,y,z,semantic,tree_id`.

## Lattice columns (density modes)

Every lattice point (occupied or empty) carries:

| column | meaning |
|---|---|
| `x,y,z` | cell/voxel centre |
| `density` | the pseudocode value (see below) |
| `num_trees` | distinct `tree_id` count among that cell/voxel's points |
| `<cat>_k` | one column **per category value in the whole run**, ratio of that cell/voxel's points in category k |

The categorical `<cat>` adapts to the source cloud (see "Source" below):
`branch_code_k` for raw/trimmed clouds, `semantic_k` for terrain/understory
clouds. `num_trees` and the ratio columns are always computed **per cell/voxel**
(its own points). The ratio columns sum to 1 where the cell has points, and are
all 0 where empty.

### density value

- **cube (`--plot-l L`)**: `L³` voxels. `density = S_voxel / S_total`
  (whole-cloud total). All voxels sum to **1**.
- **square (`--plot-l L --plot-z Z`)**: `L×L` columns × `Z` height bins =
  `L·L·Z` points. `density = S_bin / S_column` — the **per-column** fraction,
  so each non-empty column's Z values sum to **1** (a vertical profile-shape
  descriptor, useful for mimicking aerial-LiDAR returns). Each non-empty
  column's values sum to 1, empty columns to 0.

Note the deliberate **two-denominator asymmetry** in square: `density` uses the
column total `S_column`, while the ratio columns use the bin's own total. Both
describe the same bin, normalized differently.

### grid extent

The `L` grid spans the cloud's own bounding box, or the unit cube `[-1,1]³` when
`--spatial` was applied. Plot and tree levels can use different `L`/`Z` but must
share square-vs-cube.

## Uniform schema across the run

The ratio-column set is the **union of the categorical's values over every
cloud in the run** (scanned up front), so all normalized files — plot and every
tree — share one schema. A single-stem tree that lacks a category still gets
that column, all zeros, so the schema is identical across plot and trees.

## Source cloud and the categorical

`--source` defaults to the **furthest-processed cloud present**, resolved
understory → terrain → trimmed → raw. Every canonical cloud carries both
`branch_code` and `semantic_code`, so normalize chooses the meaningful one:

- **raw / trimmed** source → categorical is **`branch_code`** (branch structure:
  0, 1, 2, …), columns `branch_code_k`.
- **terrain / understory** source → categorical is **`semantic_code`** (object
  identity incl. ground/shrub/grass, which branch_code doesn't distinguish),
  columns `semantic_k`.

Trimming removes points upstream (in `assemble_plot.py`), so there is no
`--drop-trimmed` option; the trimmed cloud is already the culled geometry. If a run's plots resolve to
mixed source kinds (some understory, some only trimmed), normalize errors rather
than mixing `semantic_*` and `branch_code_*` schemas.

## Applied to both levels

Runs on the plot cloud and every per-tree cloud, writing each level's
`las/*_normalized.laz`. The scene JSON gains a `normalization` block recording
the steps, source, run branch orders, per-level L/Z, and the column schema.

## LAS writer

`las.py`'s `export_points_las` infers float32 for floating extra dims (the
`density` and ratio columns) in addition to the uint8/int32 paths; `num_trees`
is uint16.
