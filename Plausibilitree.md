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
   (individual-tree extraction, normalization, stem/canopy rasters, HELIOS++ scan
   simulation, polygon meshes) is an opt-in derived stage, never bundled into
   the core run.
5. **Metres in thought, WGS84 at the edges, UTM inside.** Users specify
   everything in metres of local plot space; ingestion and export are
   EPSG:4326 (heights in EPSG:5773), and processing internally converts to
   metric UTM (Section 10).
6. **Reproducibility.** Every stochastic draw (tree shape, rotation,
   placement, terrain, understory) is derived from explicit seeds recorded in
   the scene JSON, so any plot can be regenerated exactly.

## 2. Inputs

**[partial]** — the JSON tree list contract (Section 8.1) and generation-time
param-dict lookup are implemented; the parameter dictionary's JSON
serialization (W13) and the formal module interface (W3) are not.

The core consumes:

- A **tree list** table representing one plot (contract in Section 8.1).
- A **tree parameter dictionary**: tree code -> generator parameters
  (Section 8.6). Module specific; currently the functional-group table
  consumed by SST (Section 7.3).
- A **generator module** selection. Currently SST only; the module must
  satisfy the contract in Section 7.1.
- A **topography dict** (Section 8.7) and a **vegetation dict**
  (Section 8.8).
- **Point cloud settings**: points-per-metre budget for generation, and
  derived-stage settings such as max density for normalization (Section 9).

All interchange data between tools and stages — tree lists, parameter
dictionaries, scene/info records — is JSON; CSV appears only as a legacy
input form (item W13). A plot with its own measured positions does not run
the FIA helper at all; a hand-authored tree list in the Section 8.1 format
is a first-class input.

## 3. Core process

### 3.1 Tree generation and placement — **[implemented]**

For each tree in the tree list, in order:

- The generator module builds a tree from the tree's (X, Y, DBH, HT,
  TREE_CODE) plus that code's parameter dictionary entry. Shape parameters are
  drawn randomly (seeded; seed is a stable hash of plot seed + TREE_UID) from
  the parameter dictionary's ranges. Realized DBH and HT may therefore differ
  slightly from the input values (differences recorded in tree info JSON).
- The tree is placed at (X, Y) with a random rotation about its base, drawn
  from [0, 360) degrees, and a random reflection.

Implementation status: canonical JSON tree lists are generated per the
contract — the parameter dictionary lookup, allometric bounds and per-tree
seeded draws all happen at generation time (legacy SST parameter CSVs remain
accepted as input). Rotation is continuous and seeded. Reflection is not
implemented (item W7); the parameter dictionary is still the functional-group
CSV rather than JSON (item W13).

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

### 3.5 Downed woody debris (DWD) — **[planned]**

Final core step, after vegetation: populate downed woody debris (future work
item W10; details TBD when the stage is designed):

- downed logs of varying radius and length;
- tree stumps of varying radius and height.

Interaction rules (collision trimming against debris, shrub/grass removal
against debris points, semantic codes, volume accounting) will be specified
as part of W10.

### 3.6 Info JSONs — **[implemented]**

After all clouds exist, plot-level and tree-level info JSONs are written
(Section 8.4), including volumes with raw and trimmed variants:
trunk (frustum sum along the skeleton with taper radii), branch (cylinder
per branch element), foliage (occupied-voxel volume), canopy (convex-hull
volume of foliage).

## 4. Core outputs

**[partial]** — the design target is exactly one point cloud per plot, plus
per-tree and per-shrub JSONs carrying everything needed to extract
individuals (position, bounding box, tree_id/shrub_id, volumes). Per-tree
sub-clouds are never written as outputs. Today the pipeline writes per-tree
clouds at every stage instead (Section 11, item W6).

All point clouds are LAZ (compressed LAS, lazrs backend). All artifacts are
written into the canonical per-plot layout (Section 8.2), and one
`<PLOT_UID>_scene.json` links every artifact with a relative path, plus
bounds, rotations, and the collision log.

- **Plot cloud** (raw, trimmed, terrain, understory variants) with the
  canonical extra dimensions (Section 8.3). The `tree_id` and `shrub_id`
  dimensions plus the per-object JSONs make extracting individuals from this
  single cloud possible without any per-tree sub-clouds.
- **Plot info JSON** and, at the individual level, one **tree info JSON**
  and one **skeleton JSON** per tree, one **shrub info JSON** per shrub
  (Sections 8.4 and 7.4); positions, bounding boxes and realized
  dimensions so any consumer can crop individuals themselves.
- **DTM** (GeoTIFF, EPSG:4326 per Section 10) when one was generated.

(Per-tree clouds still exist as internal runtime artifacts — trimming
operates on them — but they are not deliverables.)

## 5. Derived output stages (opt-in)

Each is a separate stage over the core outputs, never part of the core run.

- **Individual-tree extraction** — optional convenience stage: crop per-tree
  clouds (optionally with surrounding-plot points) from the single plot
  cloud using `tree_id` and the per-tree JSONs. The core outputs alone are
  sufficient for anyone to do this themselves. **[partial]**: today per-tree
  sub-clouds are instead written directly by every core run (item W6).
- **Stem map** **[implemented as an always-on stage]** — a 2-band raster of
  the breast-height slab (DBH_HT +/- STEM_THRESHOLD, cell STEM_CELL_SIZE),
  both bands covering the whole plot: band 1 counts all points in the
  slab; band 2 is a binary mask, 1 where a tree point lies in the slab, 0
  elsewhere. Per-tree stem bounding boxes (JSON) accompany it.
  **[partial]**: the implementation writes a single-band count raster and
  expresses masks as per-tree cropped rasters (item W6).
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

Scope: FIA trees do not come with individual coordinates, and variable plot
size means trees are not recorded exhaustively over the plot. The helper
works at the scale of one unit of analysis (a whole plot, or a plot +
subplot combination).

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

1. Extract the geographic coordinate of the plot center by joining the
   FIA PLOT table (or equivalent). This join is optional: because FIA TREE
   data no longer carries tree coordinates, it is the only route to
   real-world positions. When supplied, the plot center combines with each
   tree's local placement to produce per-tree LON/LAT/ELEV (passed through
   as PLOT_LAT/PLOT_LON/PLOT_ELEV and written into the tree list's optional
   geo columns). When absent, everything stays in local coordinates and
   the geospatial ingest/export policy (Section 10) is moot — outputs are
   local and un-georeferenced.
2. Remove dead trees (DECAYCD >= 1; snags are out of scope) and rows
   missing any of: DIA, HT, TPA_UNADJ, CCLCD, or the selected
   functional-group column (SPCD, SPGRPCD, ...).
3. Convert units: DIA and HT to metres; TPA_UNADJ to trees per m^2.
4. Resolve each tree's code to a functional group (species code ->
   functional group table, coarse group fallback, DEFAULT fallback).
5. Expansion: for each tree, the expected count is
   N = ceil(TPA [trees/m^2] * plot area [m^2] * expansion factor modifier),
   with a floor of 1. A modifier of 0 therefore generates exactly the
   documented trees (no expansion), values in (0, 1) thin globally, and
   values > 1 multiply the recorded stocking.
6. Per tree row: N > 1 duplicates the row N times (each copy becomes an
   independent tree); N = 1 keeps it; N = 0 removes it.
7. Placement: sort descending by height with DBH as tiebreaker (ties broken
   randomly, seeded); place each tree at a random position inside the plot;
   the trunk disc (diameter DBH + minimum distance) must not overlap
   previously placed trunks, and the proportion of the disc outside the plot
   must be <= edge effects; failed placements retry up to the maximum number
   of tries; placed trees get ascending UIDs.
8. Write the placed tree list (Section 8.1) as JSON and a JSON report of
   parameters used, trees placed, trees not placed, and the proportion.

Implementation status: process steps 1-7 are implemented — the plot-center
join is a passthrough of pre-joined LAT/LON/ELEV columns, unit conversion,
dead-tree removal, functional-group resolution, height-descending placement
(DBH tiebreak), the expansion modifier with round-up, `min_distance`
spacing, `max_tries` and `edge_effects` controls, and the per-plot report
JSON. The placed tree list is written as canonical JSON (Section 8.1,
default) or the legacy SST parameter CSV (`--format sst-csv`);
per-tree geo derivation from plot center plus local placement is not yet
computed (item W5). A hand-authored tree list with positions replaces this
helper entirely.

### 6.2 Tree parameter editor — **[partial]**

A tool for establishing generator parameters that produce realistic trees for
a given tree code: load a tree-code parameter dictionary (or a tree list to
derive one), edit per code, sample trees from a row at chosen HT/DBH, inspect
in a 3D viewer, save and export the dictionary as JSON. **[partial]**: today
this is the SST functional-group editor (`param_editor.py`), which edits a
per-code CSV; JSON import/export is item W13.

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
formalizing the interface is item W3.

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
6. **Point sampling.** The point budget scales with both tree height and
   trunk thickness, so that vertical point density is constant across the
   plot and thick trunks carry proportionally more points than thin ones:

   Np = POINTS_PER_METRE * ceil(H) * base_dia(m)

   POINTS_PER_METRE is the number of points per 1 m height slice; the
   base-diameter factor keeps point density per unit stem surface
   comparable across DBH classes. A fraction drawn from SPR_LO/HI is spent
   on the
   trunk and branch stem meshes, the rest on foliage meshes. Foliage
   points are allocated across 1 m height slices in proportion to the
   number of branches rooted in each slice (within a slice, by mesh
   area), so better-branched crown sections are sampled more heavily.
   FP_RAND_HT (m) adds upward-only vertical jitter to foliage points,
   softening the sharp geometric canopy surface.
7. **Individual variation.** Every *_LO/*_HI pair above is drawn once per
   tree from a seeded RNG, so a dictionary entry describes a distribution
   of forms, not a single tree. The pipeline seeds each tree from a
   stable hash of the plot seed + TREE_UID (crc32-based, so regeneration
   is exact across runs and row orderings); the editor reseeds at
   seed + i * 7919 for sample i.

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
| DBH_TOL, HT_TOL | proportional size tolerances (default 0.01, i.e. +/- 1%). SST requires bounded height/base-diameter ranges, so the module bounds each tree's realized size around its tree list values: HT is drawn in [HT(1-tol), HT(1+tol)] and base diameter is derived from taper then bounded likewise | wider tolerances = more realized-size variation around the same DBH/HT |
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

**[partial]**: keyed-by-code consumption is implemented at generation time
(rows are expanded from the tree list's DBH/HT with tolerance bounds from
DBH_TOL/HT_TOL, allometry-derived crown width and base diameter); the
dictionary is a CSV table today rather than JSON (item W13).

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
  Section 8.3 at their schema dtypes via extra-bytes records (LAZ
  compression via lazrs) as **LAS 1.4**, point format 0 — **[implemented]**
  (formats 0-5 remain valid in 1.4, so only the version tag changed;
  merged-file point counts are no longer uint32-capped).

## 8. Data structures

### 8.1 Tree list (input; emitted by the FIA helper or authored by hand) — **[implemented]**

Serialized as **JSON**: one object per plot (PLOT_UID, plot shape and
dimensions, optional plot center geo) containing a `trees` array whose
elements carry the attributes below. Every tool in the pipeline exchanges
tree lists as JSON; CSV is the legacy form (item W13).
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
| LON | float | centroid longitude, EPSG:4326 | [-180.0, 180.0] | plot center + local placement; absent without a known plot center |
| LAT | float | centroid latitude, EPSG:4326 | [-90.0, 90.0] | plot center + local placement; absent without a known plot center |
| ELEV | float | base elevation, metres (EPSG:5773) | [-500, 5000] | plot center + local placement; absent without a known plot center |

**Status**: the canonical JSON contract is implemented — the FIA helper
emits `plot_<UID>.json` by default (`--format sst-csv` emits the legacy SST
parameter table) and every stage loads either contract through
`tree_list.py`; DBH/HT carry full inventory precision. Legacy SST CSVs are
accepted as input for generation. Parameter-dictionary serialization
(CSV -> JSON) remains item W13.

### 8.2 Canonical output layout

```
<PLOT_UID>/
  plot_placement/
    las/            raw / trimmed / terrain / understory / normalized clouds
    textual/        scene JSON, manifest, plot info JSON, DTM crop
    deep_learning/  stem_map / canopy_map feature+label rasters    (derived)
  individual_trees/<PLOT_UID>_<TREE_UID>/
    textual/        generation, skeleton and info JSONs  (no sub-clouds;
                    extract individuals from the plot cloud)
    deep_learning/  per-tree derived artifacts only            (derived)
  individual_shrubs/<PLOT_UID>_shrub_<SHRUB_UID>/
    textual/        shrub info JSONs         (understory stage)
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

Module-specific mapping tree code -> parameters, serialized as **JSON**
(the SST entry format is Section 7.3); edited with the Section 6.2 helper.
Height and diameter originate from the tree list, never from this
dictionary — the dictionary supplies only the proportional tolerances
(DBH_TOL, HT_TOL) that bound each tree's realized size around its
tree list values. **[partial]**: today the dictionary is the
functional-group CSV (item W13).

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
| PLOT_RADIUS | float | 7.315 | circle radius, m; only if CIRCLE (the FIA subplot radius) |
| MIN_DBH | float | 0.127 | minimum DBH for a tree to be placed, m |
| DBH_HT | float | 1.37 | breast height, m |
| POINTS_PER_METRE | int | 4096 | points per 1 m height slice, scaled by base diameter: Np = POINTS_PER_METRE x ceil(HT) x BASE_DIA |
| MAX_ATTEMPTS | int | 10 | retries when resolving tree generation |
| BRANCH_RADIUS_FRAC | float | 0.10 | branch collision radius as a fraction of branch length |
| MIN_BRANCH_RADIUS | float | 0.02 | floor for the above, m |
| COLLISION_TOLERANCE | float | 0.0 | trim radius relaxation, m |
| STEM_CELL_SIZE | float | 0.05 | stem map cell size, m (was STEM_RESOLUTION) |
| STEM_THRESHOLD | float | 0.25 | half-height of the breast-height stem slab, m (was STEM_SECTION buffer) |
| CANOPY_CELL_SIZE | float | 0.05 | canopy map cell size, m (was CHM_RESOLUTION) |
| CANOPY_BUFFER | int | 2 | canopy raster growth padding, cells |
| OCCUPIED_VOXEL_SIZE | float | 0.1 | voxel size for foliage/grass/shrub volume, m (was FOLIAGE_VOXEL / GRASS_VOXEL) |
| OCCUPIED_VOXEL_THRESHOLD | int | 1 | points needed to count a voxel occupied |
| DTM_CELL_SIZE | float | 0.25 | generated DTM cell size, m (was DTM_RESOLUTION) |
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

- **Ingest and export geospatial:** EPSG:4326 (lat/lon) with elevations in
  EPSG:5773, universally, for consistency.
- **Users think in metres.** Every user-facing length — buffers, cell
  sizes, ranges — is specified in metres of local plot space, and the tree
  list's X, Y are local metres relative to the plot origin (its lower-left
  or center).
- **Processing runs in UTM metres.** Internally, EPSG:4326 inputs are
  converted to a metric UTM projection (via the `utm` Python library, zone
  derived from the plot center) and the local plot frame is anchored there;
  every geometric stage (placement, trimming, DTM sampling, rasters) runs
  in metres of that frame.
- **Exports convert back.** Every georeferenced artifact (DTM GeoTIFF,
  stem map, canopy map) is written in EPSG:4326 with heights in EPSG:5773
  when the plot center is known. Cell sizes are specified in metres and
  converted to degrees at the plot's latitude (1 deg latitude ~ 111320 m;
  longitude metres scale by 1/cos(lat)).
- Geographic position is optional. It reaches the pipeline as the plot
  center — via the FIA helper's PLOT-table join (Section 6.1) or a
  hand-authored tree list's LAT/LON/ELEV — and per-tree real-world
  positions are the plot center plus the tree's local placement. Without
  it, the entire WGS84 ingest/export leg is moot: every artifact stays in
  the local metre frame and is marked as un-georeferenced.
- Limits of agreement: lat in [-90, 90], lon in [-180, 180], elevation in
  [-500, 5000] m (EPSG:5773).

**[partial]**: the codebase today computes in a CRS-less local metre frame
(rasters carry a local-metre geotransform with no CRS; the DTM is ASCII
with an optional EPSG sidecar; no UTM conversion is performed). Internal
UTM processing and EPSG:4326 exports are item W5.

## 11. Code-to-design gap register

Work items to bring code into agreement with this document, in rough
dependency order. Status: **open** (not yet started), **closed** (code
matches the design).

| # | Item | Gap | Status |
|---|---|---|---|
| W1 | Minimal tree list contract | generate_trees consumes per-tree frozen SST ranges; switch to (X, Y, DBH, HT, TREE_CODE) + param-dict lookup at generation, drawing seeded shapes per TREE_UID. FIA helper emits the minimal list (DBH/HT per tree), keeping placement there | closed |
| W2 | FIA helper placement controls | (a) implement expansion factor modifier; (b) minimum tree distance (generalize --allowed-overlap); (c) max tries and edge effects (spec Section 6.1). Replace the placeholder TPA rounding | closed |
| W3 | Formal generator module interface | extract the Section 7.1 contract; move SST specifics behind it | open |
| W4 | Flat-plane topography | implement topography dict None -> uniform plane ground sampling + warnings | open |
| W5 | UTM internals + EPSG:4326 exports | convert EPSG:4326 inputs to UTM metres internally (`utm` library, zone from plot center; Section 10); export GeoTIFFs back to EPSG:4326/5773; metres-thinking preserved | open |
| W6 | One plot cloud; derived stages opt-in | stop writing per-tree sub-clouds at every stage (plot cloud + per-tree/shrub JSONs are the core outputs); individual-tree extraction becomes an opt-in derived stage; stem/canopy maps and normalization become opt-in; stem map gains a second (mask) band (Sections 4-5) | open |
| W7 | Random reflection | add mirror reflection to placement draws | open |
| W8 | HELIOS++ stage | planned; needs scanner parameter design and voxel export | open |
| W9 | Polygon mesh export | planned; per-element semantics/ids; feeds HELIOS++ polygon pathway | open |
| W10 | Downed woody debris | design and implement the DWD stage (logs, stumps) per Section 3.5 | open |
| W11 | Design defaults | align code defaults with Section 9: PLOT_RADIUS 7.315, DTM_CELL_SIZE 0.25, CANOPY_CELL_SIZE 0.05, STEM_CELL_SIZE 0.05 | closed |
| W12 | LAS 1.4 writer | switch the canonical writer to LAS version 1.4 (point format 0 + extra bytes unchanged; Section 7.5) | closed |
| W13 | JSON interchange | all tool/stage interchange becomes JSON: tree list in/out of every stage, FIA helper outputs, generator parameter dictionary (CSV -> JSON), editor import/export | open |

Registry of resolved decisions (for provenance):

- Tree list contract: minimal list + tree-code parameter dict (module
  parameters never in the tree list); placement belongs to the FIA helper.
- SST-specific formats are quarantined in Section 7.
- Semantic codes: 1 trunk, 2 branch, 3 foliage, 4 shrub, 5 grass, 6 ground,
  9 noise (the earlier 8/9/0 draft values are withdrawn).
- Trimmed points are removed at trim time; no per-point trimmed dimension
  (structural record in skeletons and scene JSON counts).
- HELIOS++ export and polygon meshes are planned (not current) work.
- Coordinates: EPSG:4326/5773 ingest and export; internal processing in
  UTM metres via the `utm` library, local lengths in metres (Section 10).
- Core output is a single plot-level cloud; individual-tree information
  lives in per-tree JSONs (info + skeleton), so per-tree sub-clouds are
  never written (Section 4).
- FIA expansion: TPA converted to trees per m^2, scaled by plot area and
  the modifier, rounded up; modifier 0 keeps exactly the documented trees.
- Point budget: Np = POINTS_PER_METRE x ceil(HT) x base diameter —
  per-height-slice density normalized by trunk thickness (Section 7.2).
- Downed woody debris (logs, stumps) is planned future work (W10).
- Stem map: 2-band design — point count + tree-point mask in the breast
  height slab (Section 5).
- Defaults: CIRCLE plots, PLOT_RADIUS 7.315 m (FIA subplot); finer rasters
  (DTM 0.25, CHM 0.05, stem cell 0.05) with a robust +/-0.25 m stem slab
  (W11).
- LAS 1.4 is the target point cloud format (W12).
- DBH_TOL / HT_TOL in the parameter dictionary are proportional tolerances
  (default +/- 1%) bounding realized sizes around the tree list values;
  sizes themselves never come from the dictionary (Section 8.6).
- All interchange data (tree lists, parameter dictionaries, reports) is
  JSON; CSV is legacy (W13).
- First code-alignment round (2026-11): W2, W11 and W12 closed — FIA
  placement controls (expansion modifier with round-up, signed
  min-distance replacing the fractional allowed-overlap, max tries, edge
  effects, height-descending placement order, per-plot report JSON),
  design defaults in code (radius 7.315, DTM 0.25, CHM 0.05, stem cell
  0.05), and the LAS 1.4 canonical writer.
- Second code-alignment round (2026-11): W1 closed — canonical JSON tree
  list (Section 8.1) via tree_list.py; every stage accepts it or the
  legacy CSV; FIA helper emits JSON by default; generation looks up the
  parameter dictionary per TREE_CODE at generation time with allometric
  bounds around the tree list's DBH/HT and stable hashed per-tree seeds;
  contracts are pipeline-equivalent (identical plot clouds for the same
  trees and seeds).
