# Plot-level assembly and trimming

Stage 3: rotate and place the per-tree clouds into one scene, then use the
skeletons to find inter-tree intersections and trim the losing parts.

```
python assemble_plot.py plot_11519.csv --seed 42
```

```
plot_11519/
  plot_placement/
    plot_11519_raw.laz        merged scene, rotated + placed (all points)
    plot_11519_trimmed.laz    same scene with collision points removed
    plot_11519_scene.json     bounds, per-tree rotations
  individual_trees/
    <UID>_trimmed.laz             collision-culled points, canonical tree schema
    <UID>_trimmed_skeleton.json   placed skeleton of the surviving geometry
```

Both per-tree outputs are written for every tree, including untouched ones, so
the file set is always complete. All clouds use the canonical schema
(README_schema): `tree_id`, `shrub_id`, `semantic_code`, `classification`,
`branch_code`, `height_above_ground`, and — on tree-level clouds — `branch_id`.
Stages that haven't run yet leave their dims at 0 (shrub_id, classification,
height_above_ground are 0 out of assembly).

## Trimming removes points

Collision losers are **removed**, not flagged. The `_raw.laz` cloud is the full
placed geometry; `_trimmed.laz` is the same scene with the trimmed points
deleted, so it has fewer points (e.g. a plot might drop ~0.4% to collisions).
There is no per-point `trimmed` column and no record of removed elements
anywhere (spec 3.2) — downstream stages that want the culled
geometry simply read the trimmed cloud or the trimmed skeleton, which
contain only the surviving structure.

## Semantic and branch code

Assembly resolves the provisional generation labels into the canonical
`semantic_code` and computes `branch_code`:

- stem points on the **main trunk** → `semantic_code` 1 (trunk), `branch_code` 1
- stem points on a **fork** (order ≥ 1) → `semantic_code` 2 (branch stem),
  `branch_code` ≥ 2
- **foliage** points → `semantic_code` 3, carrying their branch element's
  `branch_code`

`branch_code` is the old branch-order idea shifted so 0 means "not a tree": 1
trunk, 2 connected-to-trunk, 3 connected-to-2, and so on. It equals the
skeleton element's order + 1 on every assigned tree point.

## Branch id

Every point gets a branch id identifying the skeleton element it was assigned
to — a specific branch, or a specific trunk sub-segment. Points on the same
branch share an id; the trunk is distinguishable from every branch.

Two encodings, because one integer can't be both plot-global and fit the LAS
int32 field (real FIA `PLT_CN` values are ~14 digits):

| form | where | value | scope |
|---|---|---|---|
| `branch_id` (int32) | LAS extra dim | compact plot-unique sequential id | unique within the plot |
| `branch_uid` (string) | skeleton JSON | `<plot_uid>_t<tree_id>_e<element>` | globally unique |

The point-cloud `branch_id` is a compact sequential id (1, 2, 3, …) assigned per
plot rather than the full `tree_id * 100000 + element` value — nothing
downstream depends on its exact value, only that it's unique per plot. Even
compacted, the id space needs int32, not int16: the compaction counter is
cumulative across every tree in the plot, so a plot with many large,
many-branched trees can assign more than int16's 32767 distinct branch ids
well before running out of trees. The full `branch_id_int` values still appear in the skeleton JSON. Unassignable points
get `branch_id = 0`. Note `branch_id` is a **tree-level** dim; the plot-level
cloud omits it (per the canonical schema).

In the trimmed skeleton, each branch element carries its own `branch_id` and
`branch_uid`, and each trunk segment carries `subsegment_branch_ids` /
`subsegment_branch_uids` lists (one entry per vertex pair). Every point-cloud
`branch_id` resolves back to an element in its skeleton via the remap.

## Branch code

Every point also gets a `branch_code` (int8) = the skeleton element's order + 1
on tree points, 0 on non-tree points:

- `0` — not a tree point (foliage still gets its branch element's code; this 0
  is for shrub/grass/ground added by later stages).
- `1` — the **main stem** (base trunk, up to the first fork).
- `2` — a **fork** or a branch attached to the main stem (1st-order structure).
- `3` — a branch attached to a fork.
- `n` — a structure attached to an (n−1)th-order one.

This is the old branch-order value shifted up by one so that 0 is a clean
"not-a-tree" sentinel. Order is computed by walking the skeleton's parent graph,
so it is correct to arbitrary depth. The current simulator forks the main stem
at most once, so in practice you see:

- single-stem trees: trunk stem `branch_code` 1, foliage 2.
- forked trees: main stem 1, forks and main-stem branches 2, fork branches 3.

`semantic_code` tracks the same split: trunk stems are semantic 1, fork/branch
stems semantic 2, foliage semantic 3.
Each branch element in the trimmed skeleton carries its own `branch_order`, and
each trunk segment carries the order of its polyline (`0` for the base stem, `1`
for each fork). The skeleton keeps the raw `branch_order`; the LAS `branch_code`
column is that value **+ 1** (so LAS trunk = 1 where the skeleton says order 0).
Unassigned points (which
class-aware assignment makes effectively impossible) fall back to 0.

If you later want a fork to count as order 0 again (part of the trunk rather
than a 1st-order structure), the rule lives entirely in `TreeSkeleton.order`.

## Rotation, reflection and placement

Each tree gets a rotation drawn in steps of `--rotation-step` degrees
(default 30, so 12 possible orientations) from [0, 360) about **its own base**,
not about the plot origin — the stems were already generated at their `X`/`Y`,
so rotating about the origin would move them off their placement. The rotation
is applied to both the cloud and the skeleton, then the base is re-anchored on
the table `X`/`Y` (a no-op unless the table changed since generation).
On top of the rotation, each tree is mirrored with probability 0.5 across the
vertical plane through its base (spec Section 3.1's random reflection) — the
mirror changes the tree's handedness while preserving all distances, radii and
skeleton topology, so collision trimming is unaffected beyond the new geometry.
Rotations in degrees and the reflection flag are recorded per tree in
`_scene.json` and in each trimmed skeleton.

## Intersection test

Skeleton elements are trunk sub-segments (one per polyline vertex pair) and
branch axes. Two elements intersect when their minimum segment-to-segment
distance is less than the sum of their radii plus `--tolerance`.

Trunk radii come from the simulator's taper, `r_base * (1 - z/h)`. **Branch
axes have no radius in the skeleton**, so one is synthesised as
`--branch-radius-frac` × branch length (default 0.10, floor 0.02 m). The
simulator builds each branch as a ribbon of half-width 0.2 × length, so 0.2 is
the geometrically literal value; 0.1 is a more forgiving reading of
"intersecting". This is by far the most influential parameter here:

| `--branch-radius-frac` | collisions | points trimmed |
|---|---|---|
| 0.05 | 36 | 0.29% |
| 0.10 (default) | 63 | 0.60% |
| 0.20 | 139 | 1.31% |

Tree pairs are AABB-prefiltered, then all element pairs are tested at once with
a vectorised clamped-parametric distance, so the cost is manageable despite
~250 elements per tree.

## Trim rules

Trees are visited in ascending `TREE_UID`. Within a contact, the loser is
chosen by:

| contact | trimmed |
|---|---|
| branch × branch | shorter path length back to the trunk; ties broken randomly (seeded) |
| branch × trunk | the branch |
| trunk × trunk | shorter stem; ties broken by higher `TREE_UID` |

Contacts are resolved deepest-interpenetration first, so the decisive contact
is settled before incidental ones nearby, and already-trimmed elements are
skipped.

Trimming an element takes its whole sub-tree with it — the "nearest joint
before the collision" is the element's start node, and everything distal to
that node goes: the rest of that stem, any forks branching above it, and every
branch attached above it. Path length back to the trunk is computed generically
by walking parents, so a higher-order branching model would not need this code
changed, even though every branch currently hangs directly off a stem.

**Interpretation worth checking:** I read "nearest joint" as the nearest
*polyline vertex* below the contact, not the nearest structural joint (base or
fork). The structural reading would mean any contact on a single-stem tree
decapitates it from ground level. The vertex reading gives 20 possible cut
points per stem and much more proportionate damage. Easy to switch if you
meant the stricter version.

## Point-to-element association

Points carry no element membership, so each is assigned to its nearest
skeleton element and inherits that element's trimmed flag. The assignment is
class-aware: stem points (semantic 1) only match trunk elements and foliage
points only match branches, because that is how they were sampled. Letting
foliage snap to a trunk element would delete canopy every time a stem was cut.

## Stem geometry notes

`generate_trees.py` rejects and reseeds any tree whose stem strays further
than `TOP_DIST + 3.5 m + slack` from its base (the fork-offset allowance
scales with the row's `SPLIT_NUM_HI`, since every cascading split can move a
stem top sideways), so a divergent stem fails loudly instead of quietly
producing bad training data.

`treesim`'s stem spline (`get_spline_points`) differs from the standard cubic
Hermite basis in the sign of one term:

```python
b01x = -kx[1]*(z1-z0) - (x1-x0)
```

where the standard Hermite form uses `+ (x1-x0)`. Endpoints are unaffected
(the term is multiplied by `t(1-t)`, zero at both knots), but interior points
are displaced: for a natural spline through x = 0, 1, 0 at z = 0, 1, 2, the
code returns 0.4375 at z = 0.5 where the Hermite form gives 0.6875. The
published equations (14)–(17) in the Simple Synth Tree paper carry the same
sign, so this is a characteristic of the published model rather than a
transcription slip — its effect is a modestly different stem curvature, not
divergence.

## Scene metadata

`_scene.json` records bounds, plot radius, per-tree rotations and reflections, the parameters
used and point counts — enough to reproduce the assembly. It also carries `plot_lat`,
`plot_lon`, `plot_elev` — the plot-center values passed straight through from
the FIA input by `fia_helper.py` (see README_fia_helper.md), read back out of
the plot table's `plot_lat`/`plot_lon`/`plot_elev` columns here and copied
into the scene JSON unchanged. All three are `null` if the FIA input didn't
have LAT/LON/ELEV to pass through. `build_info.py` reads them back out of the
scene JSON into `plot_info.json`'s `meta.plot_center` (see
README_build_info.md).

Note the scene extends past the placement radius: stems are inside the
boundary but crowns overhang, giving a 32.4 m radial extent for a 30 m plot in
the example. Nothing is clipped, since a real scan of a plot also sees the
overhanging parts of edge trees.
