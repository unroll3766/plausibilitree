# Point-cloud schema (`schema.py`)

Every cloud in the pipeline carries **one canonical set of extra dimensions**,
defined once in `schema.py` and written/read through `las.export_canonical_las`
/ `las.read_canonical_las`. Plot-level and tree-level clouds share the schema
from generation through the derived stages, so consumers never have to
branch on which stage produced a file. A stage that hasn't computed a dim yet writes its
zero/default (e.g. a raw cloud has `shrub_id = 0` everywhere); the column still
exists, so the schema never changes shape across the chain.

Point rows are always `x, y, z` plus these extra dims.

## Plot-level clouds

| dim | dtype | meaning |
|---|---|---|
| `tree_id` | int32 | 0 if not a tree; otherwise the FIA-helper tree id |
| `shrub_id` | int32 | 0 if not a shrub; assigned during shrub generation |
| `semantic_code` | int8 | object identity (below) |
| `classification` | uint8 | LAS classification (below) |
| `branch_code` | int8 | 0 if not a tree; 1 trunk, 2 connected-to-trunk, 3 connected-to-2, … |
| `branch_id` | int32 | 0 if not a branch; otherwise a compact plot-unique id from the tree skeleton |
| `height_above_ground` | float32 | point z minus DTM(x, y) |

`tree_id`/`shrub_id`/`branch_id` are all int32 rather than int16 because each
is a per-plot sequential count, and the id space has to outgrow the actual
count. int16's 32767 ceiling is easy to exceed once plot area grows: FIA's
placeholder TPA expansion scales the expected tree count with plot area, so a
100 m plot can need tree ids well past 32767 where the default ~7.3 m FIA
subplot never would. `branch_id` is remapped to a compact plot-unique
sequential id at assembly time (see README_assembly.md) rather than using the
full `tree_id * 100000 + element` value directly, but that remap counter is
cumulative *across every tree in the plot* — a plot with many large,
many-branched trees can assign more than 32767 distinct branch ids in one
plot well before running out of trees, so it needs the same int32 headroom.

## Tree-level clouds

Identical to the plot-level schema (`branch_id` is on both, so it survives the
terrain/understory chain, whose per-tree crops are built from the plot cloud).

## Semantic code (int8) — object identity

| code | meaning |
|---|---|
| 0 | ground — also "no class needed" (ground/noise are the LAS classification's job) |
| 1 | tree (trunk stem) |
| 2 | branch (branch/fork stem) |
| 3 | foliage |
| 8 | shrub |
| 9 | grass |

There are no separate unassigned or noise semantic codes: a point that
needs no object identity is 0, and noise is the traditional LAS
classification's concern (classification 7).

## LAS classification (uint8)

| code | meaning |
|---|---|
| 2 | ground (bare earth) |
| 3 | low vegetation (≤ 0.5 m above ground) |
| 4 | medium vegetation (0.5–2 m) |
| 5 | high vegetation (> 2 m) |
| 7 | noise (reserved) |

## Two orthogonal axes

`semantic_code` and `classification` are **independent**. `semantic_code` is the
object's identity (a foliage point is always semantic 3, wherever it sits);
`classification` is the height-above-ground LAS band (that same foliage point is
class 3/4/5 depending on how high it is). The terrain step sets `classification`
by height and leaves `semantic_code` untouched, so both are available.

`branch_code` mirrors `semantic_code` for tree points: trunk stems (semantic 1)
have `branch_code` 1; fork/branch stems (semantic 2) have `branch_code` ≥ 2;
foliage (semantic 3) carries the branch element's code. Non-tree points
(shrub/grass/ground) have `branch_code` 0.

## How stages populate it

- **generate_trees** — writes tree clouds with provisional `semantic_code`
  (stem → 1, foliage → 3); `branch_code`/`branch_id` stay 0.
- **assemble_plot** — resolves trunk (1) vs branch/fork (2) from the skeleton,
  sets `branch_code = order + 1`, assigns `branch_id`. **Trimming removes the
  collision points** (no flag): the trimmed cloud has fewer points than raw.
- **terrain** — adds ground points (semantic 0, class 2), sets `classification`
  by height above ground, fills `height_above_ground`. Tree `semantic_code`
  passes through unchanged.
- **understory** — appends shrubs (semantic 8) and grass (semantic 9) with
  `shrub_id` on shrubs; classifies them by height.

## Writer / reader

```python
from las import export_canonical_las, read_canonical_las

export_canonical_las(path, xyz, {"tree_id": ..., "semantic_code": ...},
                     level="plot")   # or level="tree"
xyz, cols = read_canonical_las(path, level="plot")
```

`export_canonical_las` fills any dim you don't supply with its schema default,
casts each to its canonical dtype, and sets the LAS `classification` field from
the `classification` dim. `read_canonical_las` returns `(xyz, {name: array})`
for the level's fields, defaulting any missing dim to zeros.
