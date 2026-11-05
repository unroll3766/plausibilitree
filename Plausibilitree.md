# PlausibilTree — Design Specification

PlausibiliTree is a modular assortment of software that generates realistic and
labelled point clouds of forested scenes.

This document is the authoritative design. Code and this document must agree;
where they disagree, one of them is wrong and must be fixed. Development
proceeds by changing this document first, then changing the code to match.

Each section is tagged with its implementation status:

- **[implemented]** — the codebase behaves as described
- **[partial]** — partially implemented; deviations are called out inline
- **[planned]** — designed but not yet implemented

Section 11 lists every known code-to-design gap as actionable work items.

## 1. Design principles

1. **Generator-agnostic core.** The core pipeline expresses every operation
   against the generator module contract (Section 7.1). Specific simulators
   (currently Simple Synth Tree) plug in behind that contract and may carry
   their own parameter formats and skeleton formats. No core stage may depend
   on a named simulator.
2. **Minimal, module-independent tree list.** Every plot is described by a
   tree list holding only inventory facts: position, DBH, height, tree code.
   Simulator parameters are never part of the tree list; they live in a
   separate tree-code key parameter dictionary supplied per generator module.
3. **Placement is a frontend concern.** The core never invents tree positions.
   Position sampling (expansion-factor population, spacing rules, edge
   effects, retry logic) belongs to the FIA helper, which converts inventory
   data into a placed tree list. A user who supplies their own field inventory
   with positions skips that machinery entirely.
4. **Derived outputs are separate stages.** The core run produces the labelled
   plot cloud, terrain, understory and info JSONs. Everything downstream
   (individual-tree crops, normalization, stem/canopy rasters, HELIOS++ scan
   simulation, polygon meshes) is an opt-in derived stage, never bundled into
   the core run.
5. **Metres internally, EPSG:4326 on export.** The user thinks and specifies
   in metres in a local plot frame; every georeferenced output file is written
   in EPSG:4326 (with heights in EPSG:5773) whenever the plot center is known
   (Section 10).
6. **Reproducibility.** Every stochastic draw (tree shape, rotation,
   placement, terrain, understory) is derived from explicit seeds recorded in
   the scene JSON, so any plot can be regenerated exactly.

## 2. Inputs

**[partial]** — all input kinds exist; the tree list contract below is the
target design and the param-lookup at generation time is partially in place
(see Section 7 and Section 11, items W1-W3).

The core consumes:

- A **tree list** table representing one plot (contract in Section 8.1).
- A **tree parameter dictionary**: tree code -> generator parameters
  (Section 8.6). Module specific; currently the functional-group table
  consumed by SST (Section 7.3).
- A **generator module** selection. Currently SST only; the module must
  satisfy the contract in Section 7.1.
- A **topography dict** (Section 8.7).
- A **vegetation dict** (Section 8.8).
- **Point cloud settings**: points-per-metre budget for generation, and
  derived-stage settings such as max density for normalization (Section 9).

A plot with its own measured positions does not run the FIA helper at all; a
hand-authored tree list in the Section 8.1 format is a first-class input.

## 3. Core process

### 3.1 Tree generation and placement — **[partial]**

For each tree in the tree list, in order:

- The generator module builds a tree from the tree's (X, Y, DBH, HT,
  TREE_CODE) plus that code's parameter dictionary entry. Shape parameters are
  drawn randomly (seeded; seed is a hash of plot seed + TREE_UID) from the
  parameter dictionary's ranges. Realized DBH and HT may therefore differ
  slightly from the input values (differences recorded in tree info JSON).
- The tree is placed at (X, Y) with a random rotation about its base, drawn
  from [0, 360) degrees, and a random reflection.

Implementation status: rotation is continuous and seeded [implemented];
reflection is not implemented; SST consumes frozen per-tree parameter ranges
in the plot CSV rather than looking up DBH/HT/TREE_CODE at generation time
(item W1).

### 3.2 Collision trimming — **[implemented]**

Uses the per-tree skeletons to find intersections between different trees.
For each tree, in ascending order of TREE_UID:

| Collision | Resolution |
|---|---|
| branch x branch | the one whose path length back to trunk is shorter is trimmed; tie: random (seeded) |
| branch x trunk | the branch is trimmed |
| trunk x trunk | the shorter stem is trimmed; tie: the higher TREE_UID is trimmed |

When an element loses, the whole sub-tree distal to its nearest joint before
the collision is trimmed with it (child forks, attached branches, and all
their points). Collision geometry: branch axes carry no radius in the
skeleton, so branches use `BRANCH_RADIUS_FRAC * branch_length`
(`MIN_BRANCH_RADIUS` floor) as their collision radius; trunks use their taper
radii; a `COLLISION_TOLERANCE` may relax both.

Trimmed points are removed from the cloud at trim time (final clouds are
surviving points only). Trimming is additionally recorded structurally: the
trimmed tree skeletons carry `node_trimmed` flags, and the scene JSON records
per-tree removed-point counts. There is no per-point `trimmed` dimension.

### 3.3 Topography — **[partial]**

Topography is resolved from the topography dict (Section 8.7):

- **DTM filepath**: crop to the plot boundary, sample ground points from it,
  drape the plot onto it.
- **Random terrain parameters**: generate a DTM from them, then proceed as
  above.
- **None**: sample ground points from a uniform flat plane at z = 0.

Resolution order when multiple are present: None wins if present, else
filepath wins, else random parameters; a warning describes the choice. A
malformed or missing dict is treated as None. If no sampling density is
provided with None, a default density is used and a warning is returned.

Every point receives `height_above_ground` = height above the DTM surface
(0 for ground points), and an LAS classification from it: ground = 2,
low veg (<= 0.5 m) = 3, medium veg (0.5-2 m) = 4, high veg (> 2 m) = 5.
Each tree is dropped vertically so its lowest point sits on the DTM at its
(X, Y); the DTM must cover the plot footprint (no tiling or stretching).
The DTM is an ESRI ASCII grid (.asc) in metres.

Implementation status: filepath and random generation, draping, ground
sampling, height_above_ground and classification are implemented. The None ->
flat-plane path is not: when no topography source is given, the terrain stage
is skipped entirely and no ground points or HAG are produced (item W4).

### 3.4 Vegetation (understory) — **[implemented]**

Generated from the vegetation dict (Section 8.8), draped onto the same DTM
the terrain step used (regenerated from the scene JSON record so the surface
is identical).

- **Shrubs first.** Each shrub is procedurally generated and placed at random
  within cover/pattern constraints; tree proximity is enforced via canopy
  attenuation. Any shrub point intersecting a tree point is removed.
  Shrubs are labelled with `shrub_id` and semantic code 4.
- **Grass next.** Any grass point intersecting a tree point or a shrub point
  is removed. Grass carries no per-plant identity (semantics only).

### 3.5 Info JSONs — **[implemented]**

After all clouds exist, plot-level and tree-level info JSONs are written
(Section 8.4), including volumes with raw and trimmed variants:
trunk (frustum sum along the skeleton with taper radii), branch (cylinder
per branch element), foliage (occupied-voxel volume), canopy (convex-hull
volume of foliage).

## 4. Core outputs

**[partial]** — core writes plot-level cloud, per-tree intermediate clouds,
scene JSON and info JSONs; the derived stages are currently always run
(Section 11, item W6).

All point clouds are LAZ (compressed LAS, lazrs backend). All artifacts are
written into the canonical per-plot layout (Section 8.2), and one
`<PLOT_UID>_scene.json` links every artifact with a relative path, plus
bounds, rotations, and the collision log.

- **Plot cloud** (raw, trimmed, terrain, understory variants) with the
  canonical extra dimensions (Section 8.3).
- **Per-tree intermediate clouds** (required by trimming and derived stages;
  cached in `individual_trees/`).
- **Plot info JSON** and **per-tree info JSONs** (Section 8.4).
- **Skeleton JSON** per tree, in the generator module's format
  (Section 7.4).
- **DTM** (GeoTIFF, EPSG:4326 per Section 10) when one was generated.

## 5. Derived output stages (opt-in)

Each is a separate stage over the core outputs, never part of the core run.

- **Individual-tree crops** — final per-tree clouds with surrounding-plot
  points. **[partial]**: written by every core run today.
- **Stem map** **[implemented as an always-on stage]** — a point-count raster
  of the breast-height slab (DBH_HT +/- STEM_THRESHOLD, cell
  STEM_CELL_SIZE), plus per-tree stem bounding boxes (JSON) and per-tree
  cropped maps with semantic masks for DL labels.
- **Canopy map (CHM)** **[implemented as an always-on stage]** — a canopy
  raster at CANOPY_CELL_SIZE (max height above ground per cell; foliage
  presence), plus per-tree crown bounding boxes and convex hulls.
- **Normalization** **[partial]** — subsample / rescale / voxelize clouds for
  DL consumption.
- **HELIOS++ scan simulation** — **[planned]**. Re-simulate acquisition with
  HELIOS++: input a voxelized point cloud or the polygon meshes of
  Section 5 item below, plus scanner parameters; output a .laz with
  realistic scan effects.
- **Polygon mesh export** — **[planned]**. Plot-level mesh with per-element
  attributes: semantic code (Section 8.3 values), tree UID, tree code,
  shrub UID. Required as input for the HELIOS++ polygon pathway and for
  downstream mesh consumers.

## 6. Helper tools

### 6.1 FIA inventory -> tree list — **[partial]**

Converts widely available FIA inventory data into placed tree lists. This is
where placement lives; the core never places trees.

Inputs:

- FIA TREE table (plus optional PLOT join for lat/lon/elev) for a state.
- Plot shape (circle or rectangle) and dimensions.
- Expansion factor modifier: multiplied into FIA-recorded TPA. 0 = only
  documented trees are generated; between 0 and 1 (exclusive) = global
  reduction; > 1 = more trees than recorded TPA.
- Minimum distance between trees: negative = trunks may overlap; 0 = trunks
  cannot touch; positive = required buffer.
- Maximum number of tries per tree.
- Edge effects, 0.0-1.0: the allowed proportion of a trunk's cross-sectional
  area outside the plot. The helper never places a tree 100% outside the
  plot; the full 0-1 range lets users ignore edge-effect geometry.
- Basis: whole plots or subplots; seed; optional output limit.

Process:

1. Extract the geographic coordinate of the plot center (passes through
   to outputs as PLOT_LAT/PLOT_LON/PLOT_ELEV).
2. Remove dead trees (DECAYCD >= 1) and rows missing DIA/HT/TPA/condition.
3. Resolve each tree's code to a functional group (species code ->
   functional group table, coarse group fallback, DEFAULT fallback).
4. Populate the expected tree count from the expansion factor and the
   user-supplied modifier for the user-supplied plot dimensions.
5. Placement: sort descending by height with DBH as tiebreaker (ties broken
   randomly, seeded); place each tree at a random position inside the plot;
   the trunk disc (diameter DBH + minimum distance) must not overlap
   previously placed trunks, and the proportion of the disc outside the plot
   must be <= edge effects; failed placements retry up to the maximum number
   of tries; placed trees get ascending UIDs.
6. Write the placed tree list (Section 8.1) and a JSON report of parameters
   used, trees placed, trees not placed, and the proportion.

Implementation status: unit conversion, dead-tree removal, functional-group
resolution, `--allowed-overlap` spacing and the CSV/JSON outputs are
implemented. TPA expansion is a placeholder (round-to-nearest); the
expansion factor modifier, minimum distance, max tries and edge-effects
controls are not implemented (items W2a-W2c). A user-supplied inventory with
positions replaces this helper entirely.

### 6.2 Tree parameter editor — **[partial]**

A tool for establishing generator parameters that produce realistic trees for
a given tree code: load a tree-code parameter dictionary (or a tree list to
derive one), edit per code, sample trees from a row at chosen HT/DBH, inspect
in a 3D viewer, save the dictionary. **[partial]**: today this is the SST
functional-group editor (`param_editor.py`), which edits a per-code CSV and
samples trees per row; it satisfies the design against the SST dict format.

### 6.3 Environment / vegetation parameter editor — **[implemented]**

An interactive editor over an already-generated plot: 3D viewer for the
trimmed plot; edit the topography and vegetation dicts (defaults loaded if
absent); regenerate terrain/understory from the edits and toggle between
original and edited clouds; save the dicts. (`scene_editor.py`.)

### 6.4 Environment replicates — **[implemented]**

Replicates a tree list N times, each with a randomly drawn terrain/shrub/
grass sidecar JSON that the pipeline applies per plot (`replicate.py`).

## 7. Generator module contract (with Simple Synth Tree)

### 7.1 Module contract — **[planned as a formal interface]**

A generator module must provide:

- `build(rng, DBH, HT, tree_code, params) -> TreeModel` where `TreeModel`
  holds:
  - points (N x 3, metres, origin at tree base) with a trunk / branch / foliage
    label per point;
  - a skeleton with trunk and branch centerlines carrying radii (or a
    documented proxy) and joint structure sufficient for collision trimming;
  - realized dimensions (realized DBH, realized height);
- a parameter dictionary format (Section 8.6) mapping tree code ->
  parameters.

The core calls the module through this contract only. **[partial]**: today
the core (generate/assemble/build_info) calls SST conventions directly;
formalizing the interface is item W1.

### 7.2 Conceptual model: from tree list entry to synthetic tree — **[implemented]**

A synthetic tree is produced from three separable sources of information:

| Source | Supplies | Example |
|---|---|---|
| Tree list (Section 8.1) | dimensions and position: X, Y, DBH, HT, TREE_CODE | "a 0.35 m DBH, 28 m tall tree of code DF at (12.4, 3.1)" |
| Parameter dictionary (per tree code, Section 7.3) | form: allometry exponents, forking behaviour, crown shape, point budget split | "Douglas-fir-like: rarely forks, high crown base, narrow crown" |
| Seeded per-tree draws | individual variation within the code's parameter ranges | "this particular tree forked twice and leans 1.2 m" |

Nothing about the tree's size lives in the dictionary and nothing about the
code's form lives in the tree list; the two meet only inside the generator
module call (Section 7.1). Realized DBH and HT may deviate slightly from the
input values; the deviations are recorded in the tree info JSON.

Synthesis chain for one tree (H = HT, DBH in metres, DBH_HT = 1.37 m):

1. **Trunk diameter profile.** DBH is measured at 1.37 m; the trunk swells
   below that point. Base diameter, with b = TAPER_B:

   d(0) = DBH * (H / (H - 1.37))^(1/b)

   Small b gives strong butt swell (a large d(0) relative to DBH). From the
   base radius r0 = d(0)/2 the trunk tapers linearly to zero at the tip:

   r(h) = r0 * (1 - h/H)

   This taper also defines the collision radii used by trimming (Section
   3.2) and the frustum volumes in the info JSONs.
2. **Crown width (allometry).** With DBH converted to inches (DBH_in):

   CW_ft = exp(CW_A0 + CW_A1 * ln(DBH_in)),  MAX_CAN_W (m) = 0.3048 * CW_ft

   A0 sets crown width at any given size; A1 sets how strongly width grows
   with DBH.
3. **Stem path.** The stem runs base -> mid control point -> tip. TOP_DIST
   (m) is the maximum horizontal wander of the tip (the whole tree's lean);
   MID_DIST (m) is the wander of the mid-stem control point (stem bow).
   These change stem geometry only, not the crown.
4. **Forking.** With probability SPLIT_PROB the stem forks at all. If it
   forks, a count drawn from SPLIT_NUM_LO/HI cascades splits: each split
   targets the tallest remaining stem, so forks climb the dominant leader
   and new crowns form in the upper canopy, interlocked with their
   neighbours. Each split is placed at a height drawn as a fraction
   (SPLIT_HT_LO/HI) of the stem segment being split.
5. **Crown envelope.** Branches exist between the crown base
   (MIN_CAN_HT fraction of H) and the tip; the crown is widest
   (MAX_CAN_W) at the height given by MAX_CAN_W_HT (fraction of H).
   Branch lengths ramp from zero at the crown base to full width at the
   widest height, then taper toward the tip. A count drawn from
   NUM_BRANCH_LO/HI is split between a lower crown band and an upper
   crown/fork band.
6. **Point sampling.** The total budget scales with the tree size
   (POINTS_PER_METRE). A fraction drawn from SPR_LO/HI is spent on the
   trunk and branch stem meshes, the rest on foliage meshes. Foliage
   points are allocated across 1 m height slices in proportion to the
   number of branches rooted in each slice (within a slice, by mesh
   area), so better-branched crown sections are sampled more heavily.
   FP_RAND_HT (m) adds upward-only vertical jitter to foliage points,
   softening the sharp geometric canopy surface.
7. **Individual variation.** Every *_LO/*_HI pair above is drawn once per
   tree from a seeded RNG, so a dictionary entry describes a distribution
   of forms, not a single tree. (The editor reseeds at seed + i * 7919
   for sample i; the design moves these draws to pipeline generation time
   keyed by plot seed + TREE_UID, item W1.)

Because raw trunk volume is a frustum sum along this profile,
V = pi/3 * L * (r0^2 + r0*r1 + r1^2) for each 3D segment between skeleton
nodes at radii r0, r1, the inventory inputs (DBH, HT) alone largely
determine trunk volume; the dictionary parameters change how that volume is
distributed and how much crown accompanies it.

### 7.3 SST parameter dictionary

The SST parameter dictionary is the functional-group table (one row per tree
code). Each entry is read at generation time, shape parameters are drawn from
the LO/HI ranges, and crown width / base diameter are derived from DBH and
HT by the entry's allometry (Section 7.2).

| Column(s) | Meaning | Raising it ... |
|---|---|---|
| FUNC_GROUP | tree code this entry is keyed by | — |
| CW_A0, CW_A1 | crown-width allometry: CW (ft) = exp(A0 + A1 ln(DBH_in)), then to metres | widens the crown at a given DBH (A0) / makes big trees disproportionately wide (A1); wider crowns mean more collision trimming in dense plots |
| TAPER_B | butt-swell exponent, Section 7.2 step 1 | reduces butt swell (base diameter closer to DBH); lowering it thickens the trunk base |
| SPLIT_PROB | probability the tree forks at all (0 = single straight stem) | more forked trees |
| SPLIT_NUM_LO/HI | number of cascading trunk splits, drawn per tree (default 1:1 = single fork) | more levels of forking and bushier multi-leader crowns; each leader's crown is airier at the same point count, so raise the point budget for dense multi-leader crowns |
| SPLIT_HT_LO/HI | where each split happens, as a fraction of the stem being split | pushes splits higher up the stems |
| NUM_BRANCH_LO/HI | nominal branch count, drawn per tree and split between the lower crown band and the upper crown/fork band | denser canopy and more foliage points |
| MIN_CAN_HT_LO/HI | crown base: height of the lowest branch, as a fraction of H | raises the crown base (pruned-plantation look); more bare, visible stem |
| MAX_CAN_W_HT_LO/HI | height (fraction) where the crown is widest | moves the widest part of the crown upward |
| TOP_DIST | maximum horizontal wander of the stem tip, m (tree lean) | more leaning trees |
| MID_DIST | horizontal wander of the mid-stem control point, m (stem bow) | more bowed, curved stems |
| FP_RAND_HT | upward vertical jitter added to foliage points, m | fuzzier, less clearly bounded canopy |
| SPR_LO/HI | stem-point fraction of the point budget, drawn per tree | more points on trunk and branches, at the cost of canopy detail |

Useful silhouettes: a **conifer** look is low SPLIT_PROB, small CW_A0/A1,
high crown base variation and high NUM_BRANCH; a **hardwood** look is high
SPLIT_PROB, high CW_A0/A1, a low crown base, and SPLIT_NUM_LO/HI above
1:1 for cascading multi-leader crowns. TOP_DIST and MID_DIST only change
stem geometry, so they add variety without changing the crown.

**[partial]**: keyed-by-code consumption already exists in the editor and in
the helper's derivations; the pipeline path still carries frozen per-tree
ranges in the plot CSV, and the tree's height/diameter tolerance ranges are
computed by the FIA helper around the inventory values rather than drawn at
generation time (item W1).

### 7.4 SST skeleton JSON — **[implemented]**

Per tree, a JSON describing the generated structure:

- `tree_uid`, `plot_uid`, `tree_id`, `coordinate_frame`, `base_xy`,
  `height_m`, `base_radius_m`, `fork`, `fork_height_m`, `n_splits`,
  `splits`, `min_canopy_height_m`, `max_canopy_width_height_m`;
- `trunk`: centerline polylines per segment with per-segment z-range,
  length, and (trimmed variants) `node_trimmed` flags;
- `branches`: same organization for branch elements, each with its length
  and the branch element id that its points carry.

Branch connectivity (joint structure) is expressed by the polylines and
splits; there is no separate SOURCE/DESTINATION table.

### 7.5 SST codebase modifications

- `treesim.py`: optional skeleton export (Section 7.4); a translation factor
  so tree tops can grow towards any point.
- `las.py`: canonical LAZ writer writing x/y/z plus the extra dims of
  Section 8.3 at their schema dtypes (LAS 1.2, point format 0, extra bytes
  records; LAZ compression via lazrs).

## 8. Data structures

### 8.1 Tree list (input; emitted by the FIA helper or authored by hand)

Required attributes:

| Attribute | Type | Description | Range | Notes |
|---|---|---|---|---|
| TREE_UID | string | unique tree id in the list | | ex: 0001, 0002 |
| X | float | local x of tree centroid, metres from lower-left | [0, PLOT_LENGTH] rect / [0, 2 PLOT_RADIUS] circle | |
| Y | float | local y of tree centroid, metres from lower-left | [0, PLOT_WIDTH] rect / [0, 2 PLOT_RADIUS] circle | |
| DBH | float | trunk diameter at breast height, metres | (MIN_DBH, inf) | |
| HT | float | tree height above ground, metres | (0, 120) | tallest known tree is ~116 m |
| TREE_CODE | string | vegetation type code (species, functional group) | | resolution order: species code -> functional group -> DEFAULT |

Optional attributes:

| Attribute | Type | Description | Range | Notes |
|---|---|---|---|---|
| Z | float | local z of tree base, metres | [0, inf) | typically 0 |
| LON | float | centroid longitude, EPSG:4326 | [-180.0, 180.0] | |
| LAT | float | centroid latitude, EPSG:4326 | [-90.0, 90.0] | |
| ELEV | float | base elevation, metres (EPSG:5773) | [-430, 5000] | Dead Sea shore to Andes treeline |

**[partial]**: the current CSV contract is the SST per-tree parameter table
(TREE_UID, FUNC_GROUP, X, Y, HT_LO/HI, BASE_DIA_LO/HI, SPLIT_*, ...,
PLOT_LAT/LON/ELEV); the minimal contract above is the target (item W1).

### 8.2 Canonical output layout

```
<PLOT_UID>/
  plot_placement/
    las/            raw / trimmed / terrain / understory / normalized clouds
    textual/        scene JSON, manifest, plot info JSON, DTM crop
    deep_learning/  stem_map / canopy_map feature+label rasters    (derived)
  individual_trees/<PLOT_UID>_<TREE_UID>/
    las/ textual/ deep_learning/          (intermediates + derived crops)
  individual_shrubs/<PLOT_UID>_shrub_<SHRUB_UID>/  (understory stage)
```

### 8.3 Point cloud dimensions — **[implemented]**

x, y, z plus:

| Dimension | Dtype | Meaning |
|---|---|---|
| tree_id | int32 | 0 = not a tree; per-plot sequential otherwise |
| shrub_id | int32 | 0 = not a shrub |
| semantic_code | int8 | object identity, below |
| classification | uint8 | LAS class by height above ground: 2 ground, 3 low (<=0.5 m), 4 medium (0.5-2 m), 5 high (>2 m), 7 noise |
| branch_code | int8 | 0 not-tree; 1 trunk; 2+ = branch element code |
| branch_id | int32 | branch element identity within the plot (0 = unassigned) |
| height_above_ground | float32 | z minus DTM at (x, y) |

Semantic codes: 0 unassigned, 1 trunk, 2 branch stem, 3 foliage, 4 shrub,
5 grass, 6 ground, 9 noise.

### 8.4 Info JSONs — **[implemented]**

Plot info (`<PLOT_UID>_info.json`): plot_uid, plot_shape (with radius or
length/width), n_trees, n_shrubs, bounding box, volumes
{trunk, branch, foliage, canopy} x {raw, trimmed}, heights (voxel counts /
percentiles), vegetation characteristics (percent cover), neighborhood
statistics, and `meta` (plot center lon/lat/elev when known, coordinate
system, units).

Tree info (per tree): tree_uid, dbh_m, height_m, volumes as above, bounding
box, crown major/minor axis, orientation, max width and its height,
and neighborhood counts. Shrub info analogous per shrub.

Volume methods: trunk = frustum sum along skeleton centerlines with taper
radii; branch = cylinder per branch (pi r^2 L); foliage = occupied-voxel
volume at OCCUPIED_VOXEL_SIZE; canopy = convex-hull volume of foliage.

### 8.5 Scene JSON — **[implemented]**

One per plot; relative links to every artifact; plot bounds; per-tree
rotations; collision log with removed-point counts; terrain record (DTM
source, cell size, generated-ness) sufficient to reproduce the DTM; topo
dict + vegetation settings.

### 8.6 Tree parameter dictionary

Module-specific mapping tree code -> parameters (SST format in Section 7.3);
edited with the Section 6.2 helper. Height and diameter come from the tree
list, never from this dictionary.

### 8.7 Topography dict

One of:

- DTM filepath (ESRI ASCII .asc, metres)
- parameters for random terrain generation (cell size, bumpiness, relief)
- None (flat plane)

plus ground sampling density (points per DTM cell) and ground noise (metres,
vertical jitter). Resolution and conflict rules in Section 3.3.

### 8.8 Vegetation dict

- Shrubs: cover fraction, horizontal cell size, min/max height (m), min/max
  width (m), points per shrub, placement margin (m), tree proximity
  (canopy attenuation 0-1).
- Grass: cover fraction, min/max height (m), patch scale (m), candidate
  point density (per m^2), tree proximity (canopy attenuation 0-1).
- Shared: occupied-voxel size and threshold for volume accounting.

## 9. Variables

| Variable | Type | Default | Description |
|---|---|---|---|
| PLOT_SHAPE | string | CIRCLE | plot geometry, CIRCLE or RECTANGLE |
| PLOT_LENGTH, PLOT_WIDTH | (float, float) | (30.0, 30.0) | rectangle sides, m; only if RECTANGLE |
| PLOT_RADIUS | float | 15.0 | circle radius, m; only if CIRCLE. FIA subplot ladder uses 7.315 m |
| MIN_DBH | float | 0.127 | minimum DBH for a tree to be placed, m |
| DBH_HT | float | 1.37 | breast height, m |
| POINTS_PER_METRE | int | 4096 | generation point budget per metre of tree height |
| MAX_ATTEMPTS | int | 10 | retries when resolving tree generation |
| BRANCH_RADIUS_FRAC | float | 0.10 | branch collision radius as a fraction of branch length |
| MIN_BRANCH_RADIUS | float | 0.02 | floor for the above, m |
| COLLISION_TOLERANCE | float | 0.0 | trim radius relaxation, m |
| STEM_CELL_SIZE | float | 0.1 | stem map cell size, m (was STEM_RESOLUTION) |
| STEM_THRESHOLD | float | 0.25 | half-height of the breast-height stem slab, m (was STEM_SECTION buffer) |
| CANOPY_CELL_SIZE | float | 0.1 | canopy map cell size, m (was CHM_RESOLUTION) |
| CANOPY_BUFFER | int | 2 | canopy raster growth padding, cells |
| OCCUPIED_VOXEL_SIZE | float | 0.1 | voxel size for foliage/grass/shrub volume, m (was FOLIAGE_VOXEL / GRASS_VOXEL) |
| OCCUPIED_VOXEL_THRESHOLD | int | 1 | points needed to count a voxel occupied |
| DTM_CELL_SIZE | float | 1.0 | generated DTM cell size, m (was DTM_RESOLUTION) |
| DTM_BUMPINESS | float | 0.5 | generated terrain micro-relief scale |
| DTM_RELIEF | float | 10.0 | generated terrain total relief, m |
| GROUND_DENSITY | float | 1.0 | ground points per DTM cell |
| GROUND_NOISE | float | 0.0 | ground point vertical jitter, m |
| SHRUB_COVER | float | 0.15 | fraction of plot area eligible for shrubs |
| SHRUB_CELL_SIZE | float | 2.0 | shrub placement cell size, m |
| SHRUB_HEIGHT_RANGE | (float, float) | (0.4, 1.8) | shrub height range, m |
| SHRUB_WIDTH_RANGE | (float, float) | (0.8, 1.6) | shrub width range, m |
| POINTS_PER_SHRUB | int | 300 | points sampled per shrub |
| GRASS_COVER | float | 0.6 | fraction of plot area in grass |
| GRASS_HEIGHT_RANGE | (float, float) | (0.03, 0.35) | grass height range, m |
| GRASS_PATCH_SCALE | float | 5.0 | grass vigor patch scale, m |
| GRASS_CANDIDATE_DENSITY | float | 200.0 | grass candidate density, per m^2 |
| CANOPY_ATTENUATION | float | 0.3 | vegetation reduction near tree canopies, 0-1 |
| MAX_DENSITY | int | 5000 | normalization: max points per m^3/cube (derived stage) |

## 10. Coordinate reference policy — **[planned]**

- All computation, parameters and user-facing lengths are in metres in a
  local plot frame (z = 0 at DTM base), independent of any CRS.
- Every exportable georeferenced artifact (DTM GeoTIFF, stem map, canopy
  map) is written in EPSG:4326 with heights in EPSG:5773, when the plot
  center is known from the tree list's LON/LAT/ELEV. Cell sizes are
  specified in metres and converted to degrees at the plot's latitude
  (1 deg latitude ~ 111320 m; longitude metres scale by 1/cos(lat)).
- When no geographic position is known, artifacts stay in the local
  metre frame and are marked as un-georeferenced.
- Limits of agreement: lat in [-90, 90], lon in [-180, 180], elevation in
  [-430, 5000] m (EPSG:5773; Dead Sea shore to Andes treeline).

**[partial]**: rasters currently carry a local-metre geotransform with no
CRS; the DTM is ASCII with optional EPSG sidecar. Reprojecting exports is
item W5.

## 11. Code-to-design gap register

Work items to bring code into agreement with this document, in rough
dependency order:

| # | Item | Gap |
|---|---|---|
| W1 | Minimal tree list contract | generate_trees consumes per-tree frozen SST ranges; switch to (X, Y, DBH, HT, TREE_CODE) + param-dict lookup at generation, drawing seeded shapes per TREE_UID. FIA helper emits the minimal list (DBH/HT per tree), keeping placement there |
| W2 | FIA helper placement controls | (a) implement expansion factor modifier; (b) minimum tree distance (generalize --allowed-overlap); (c) max tries and edge effects (spec Section 6.1). Replace the placeholder TPA rounding |
| W3 | Formal generator module interface | extract the Section 7.1 contract; move SST specifics behind it |
| W4 | Flat-plane topography | implement topography dict None -> uniform plane ground sampling + warnings |
| W5 | EPSG:4326 exports | write GeoTIFFs in EPSG:4326 (Section 10), metres-thinking preserved; keep local frame internal |
| W6 | Derived stages opt-in | make individual-tree crops, stem/canopy maps and normalization opt-in stages; core run stops after understory + info |
| W7 | Random reflection | add mirror reflection to placement draws |
| W8 | HELIOS++ stage | planned; needs scanner parameter design and voxel export |
| W9 | Polygon mesh export | planned; per-element semantics/ids; feeds HELIOS++ polygon pathway |

Registry of resolved decisions (for provenance):

- Tree list contract: minimal list + tree-code parameter dict (module
  parameters never in the tree list); placement belongs to the FIA helper.
- SST-specific formats are quarantined in Section 7.
- Semantic codes: 1 trunk, 2 branch, 3 foliage, 4 shrub, 5 grass, 6 ground,
  9 noise (the earlier 8/9/0 draft values are withdrawn).
- Trimmed points are removed at trim time; no per-point trimmed dimension
  (structural record in skeletons and scene JSON counts).
- HELIOS++ export and polygon meshes are planned (not current) work.
- Raster georeferencing: metres thinking, EPSG:4326 exports (Section 10).
