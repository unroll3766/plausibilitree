# PlausibilTree — Design Specification

PlausibiliTree is a modular assortment of software that generates realistic and
labelled point clouds of forested scenes.

`notes.md` is the upstream design document, written and commented by the
project owner; wherever this document disagrees with `notes.md`, `notes.md`
wins. This document is the specification of record: it tracks `notes.md`
section by section, tags implementation status, and turns every deviation
between `notes.md`, this spec, and the code into a register entry
(Section 11). Development proceeds by reconciling `notes.md` here first,
then changing the code to match.

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

**[implemented]** — the tree list contract (Section 8.1, both
serializations — CSV canonical per notes.md since W17, JSON alternate),
generation-time param-dict lookup, and the parameter dictionary's JSON
serialization (8.6).

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

Tree lists interchange as CSV and are also accepted as JSON (Section 8.1);
parameter dictionaries, scene and info records are JSON. A plot with its
own measured positions does not run the FIA
helper at all; a hand-authored tree list in the Section 8.1 format is a
first-class input.

## 3. Core process

### 3.1 Tree generation and placement — **[implemented]**

For each tree in the tree list, in order:

- The generator module builds a tree from the tree's (X, Y, DBH, HT,
  TREE_CODE) plus that code's parameter dictionary entry. Shape parameters are
  drawn randomly (seeded; seed is a stable hash of plot seed + TREE_UID) from
  the parameter dictionary's ranges. Realized DBH and HT may therefore differ
  slightly from the input values (differences recorded in tree info JSON).
- The tree is placed at (X, Y) with a random rotation about its base, drawn
  from [0, 360) in steps of N degrees (ROTATION_STEP, seeded), and a random
  reflection.

Implementation status: canonical JSON tree lists are generated per the
contract — the parameter dictionary lookup, allometric bounds and per-tree
seeded draws all happen at generation time (legacy SST parameter CSVs remain
accepted as input). Rotation is quantized: a multiple of ROTATION_STEP
degrees is drawn from [0, 360) (default 30; the scene JSON records
rotation_step_deg); reflection is a
seeded 50/50 mirror draw across the vertical plane through the tree's
base, applied with the rotation at plot assembly. The parameter dictionary
is JSON (`data/functional_groups.json`, Section 8.6).

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

Trimmed points are removed from the cloud at trim time. No record of a
trimmed point or branch is retained anywhere (owner ruling 2026-11-07):
there is no per-point `trimmed` dimension, no skeleton node or element
flag, no trim summary, no collision log, no removed-point count, and no
pre-trim ("raw") volume variant — the surviving geometry is the only
output.

### 3.3 Topography — **[implemented]**

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
The DTM is an ESRI ASCII grid (.asc) in metres (internal working
format); for plot-anchored runs the generated plot DTM is additionally
exported as a GeoTIFF with compound CRS EPSG:9707 and absolute heights
(Section 10).

Implementation status: all three topography resolutions are implemented —
DTM filepath (draped with a coverage check, no tiling/stretching), random
parameters (generated, saved as an .asc, cell DTM_RESOLUTION), and None (a
uniform flat plane at z = 0; the terrain stage always runs, sampling ground
points at the default density and warning when none was configured). Conflicts
resolve per the None > filepath > random order with a warning describing the
ignored choice; a malformed or missing topography behaves as None.

### 3.4 Vegetation (understory) — **[implemented]**

Generated from the vegetation dict (Section 8.8), draped onto the same DTM
the terrain step used (regenerated from the scene JSON record so the surface
is identical).

- **Shrubs first.** Each shrub is procedurally generated and placed at random
  within cover/pattern constraints; tree proximity is enforced via canopy
  attenuation. Any shrub point intersecting a tree point is removed.
  Shrubs are labelled with `shrub_id` and the shrub semantic code (8).
  Candidate shrub cells are constrained to the plot's own footprint
  (a circular plot seeds shrubs only inside the disc, not in the bounding
  square's corners), and a shrub that retains no points after clipping is
  dropped from the per-shrub record set with a warning — no shrub ever has
  a record without points (2026-11 review, item W23).
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
(Section 8.4) under the notes.md field names: trunk (frustum sum along
the skeleton with taper radii), branch (cylinder per branch element),
foliage (occupied-voxel volume at FOLIAGE_VOXEL), grass (occupied-voxel
volume at GRASS_VOXEL) — one post-trim value each, no pre-trim variants
(Section 3.2). The canopy convex-hull volume extension is withdrawn
(2026-11 review, item 10).

## 4. Core outputs

**[implemented]** — exactly one point cloud per plot, plus per-tree and
per-shrub JSONs carrying everything needed to extract individuals
(position, bounding box, tree_id/shrub_id, volumes). Per-tree sub-clouds
are never written as outputs (individual extraction is the opt-in stage
of Section 5).

All point clouds are LAZ (compressed LAS, lazrs backend). All artifacts are
written into the canonical per-plot layout (Section 8.2), and one
`<PLOT_UID>_scene.json` links every artifact with a relative path, plus
bounds, rotations and reflections.

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
operates on them — but they are not deliverables: generate/assemble stage
them in the flat working area and restructure deletes them, so nothing
per-tree ever lands in `individual_trees/*/las/` from a core run. The same
applies to per-shrub clouds; the extraction stage is the only stage that
crops per-object clouds out of the plot cloud.)

## 5. Derived output stages (opt-in)

Each is a separate stage over the core outputs, never part of the core run.

- **Individual-tree extraction** **[implemented (item W6)]** — optional
  convenience stage (`extract_trees.py`, orchestrated as `--extract-trees`):
  crop per-object clouds from the single plot cloud using `tree_id` /
  `shrub_id` — trees and shrubs alike, each written to the canonical
  `individual_trees/<TREE_UID>/las/` / `individual_shrubs/<SHRUB_UID>/las/`
  positions under the source name (`--source` resolves understory ->
  terrain -> trimmed -> raw; the extracted file is named for that source
  so the scene links, normalize and the editor read it unchanged).
  `--surroundings N` widens each crop to an N-metre window of the plot so
  the extracted cloud carries background points; without it a crop is the
  object's own points only. The core outputs alone are sufficient for
  anyone to do this themselves.
- **Stem map** **[implemented as an opt-in stage (`--stem-map`)]** — a 2-band raster of
  the breast-height slab (DBH_HT +/- STEM_SECTION, cell STEM_RESOLUTION),
  both bands covering the whole plot: band 1 (uint16) counts all points
  in the slab; band 2 (uint16, 0/1) is a binary mask, 1 where a *tree*
  point (`tree_id != 0`) lies in the slab, 0 elsewhere — so band 2 is
  tree points only, never ground or shrub/shrub-and-grass returns.
  Per-tree stem bounding boxes (JSON) accompany it. Per-tree instance
  masks (integer `tree_id` labels, skeleton-derived) remain separate
  per-tree rasters.
- **Canopy map (CHM)** **[implemented as an opt-in stage (`--canopy-map`)]** — a 2-band
  raster at CHM_RESOLUTION: band 1 (float32) is the max height above
  ground per cell over vegetation points (semantic codes 1, 2, 3, 8; grass
  code 9 is ground cover and contributes nothing to a canopy height map),
  0 where a cell holds no vegetation; band 2 (float32, 0/1) is 1 where a
  foliage point (code 3) lies in the cell. Plus per-tree crown bounding
  boxes and convex hulls.
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
8. Write the placed tree list (Section 8.1) as CSV and a JSON report of
   parameters used, trees placed, trees not placed, and the proportion.

Implementation status: process steps 1-7 are implemented — the plot-center
join is a passthrough of pre-joined LAT/LON/ELEV columns, unit conversion,
dead-tree removal, functional-group resolution, height-descending placement
(DBH tiebreak), the expansion modifier with round-up, `min_distance`
spacing, `max_tries` and `edge_effects` controls, and the per-plot report
JSON. The placed tree list is written as the minimal per-tree CSV (default;
tree rows only — plot geometry and plot center are recorded in the report
JSON, which `tree_list` reads back when loading CSV lists). The JSON
tree-list form remains the `--format tree-list-json` alternate and
`--format sst-csv` the legacy parameter table.
Per-tree geo derivation from plot center plus local placement is not yet
computed (item W5). A hand-authored tree list with positions replaces this
helper entirely.

### 6.2 Tree parameter editor — **[partial]** (item W20)

A tool for establishing generator parameters that produce realistic trees for
a given tree code: load a tree-code parameter dictionary (or a tree list to
derive one), edit per code, sample trees from a row at chosen HT/DBH, inspect
in a 3D viewer, save and export the dictionary as JSON. **[partial]**: JSON
import/export is implemented (`param_editor.py`, default
`data/functional_groups.json`; legacy CSVs open/save too); deriving a
dictionary from a tree list is not (item W20).

### 6.3 Environment / vegetation parameter editor — **[implemented]**

An interactive editor over an already-generated plot: 3D viewer for the
trimmed plot; edit the topography and vegetation dicts (defaults loaded if
absent); regenerate terrain/understory from the edits and toggle between
original and edited clouds; save the dicts. (`scene_editor.py`.)

### 6.4 Environment replicates — **[implemented]**

Replicates a tree list N times, each with a randomly drawn terrain/shrub/
grass sidecar JSON that the pipeline applies per plot (`replicate.py`).

### 6.5 Tree skeleton viewer — **[planned]** (item W26)

A 3D GUI helper that loads a tree skeleton JSON (raw or trimmed). Toggles:
branch radii (cylinders vs polylines), branch ID display, branch angle
display, branch-order colorization; zoom/rotate/pan; render the current
view to .png/.jpeg. Input: one skeleton JSON. Output: TBD (2026-11
review, item 11 — helper module noted; full implementation TBD).

## 7. Generator module contract (with Simple Synth Tree)

### 7.1 Module contract — **[implemented]**

A generator module must provide (implemented in `generator_module.py`; SST
in `sst_generator.py`, selectable by name from the registry —
`generate_trees.py -g` and `run_pipeline.py -g`):

- `build(rng, DBH, HT, tree_code, params) -> TreeModel` where `TreeModel`
  holds:
  - points (N x 4: x, y, z in metres plus one semantic label per point),
    **origin at the tree base** (z = 0 at the base);
  - a placement helper `place(model, x, y)` that positions the base-origin
    model in plot coordinates (points plus every positional field of the
    skeleton);
  - a skeleton (the module's documented format; SST: Section 7.4) with
    trunk and branch centerlines carrying radii (or a documented proxy)
    and joint structure sufficient for collision trimming;
  - realized dimensions and provenance (model params, stem fraction,
    point budget). The realized record carries the height, the base
    radius and the taper-implied DBH (diameter at DBH_HT) so
    `<TREE_UID>_generation.json` reports the realized DBH alongside
    `base_radius_m` (2026-11 review, item 9; item W25).
- `load_param_dict(path)`: the parameter dictionary format (Section 8.6)
  mapping tree code -> parameters.
- `resolve_tree_code(code)`: a tree code -> dictionary key.
- `entry_for(code, dbh, ht)`: the resolved dictionary entry for a tree,
  gated by the entry's applicability bounds when the tree's DBH/HT are
  given (Section 7.3), so the entry the core records for a tree is the
  one the build consumes.

The core calls the module through this contract only: it seeds each tree
from a stable hash of the plot seed + TREE_UID, hands the module a fresh
numpy Generator, and retries a degenerate build with a rehashed seed
(MAX_ATTEMPTS). The module must draw ALL randomness from that rng (SST's
treesim uses numpy's global generator, so its adapter bridges the rng
with one deterministic draw). The trimming and volume stages consume the
configured module's documented skeleton format, which is therefore part
of each module's documentation (SST: Section 7.4).

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
   below that point. Base diameter, with b = TAPER:

   d(0) = DBH * (H / (H - 1.37))^(1/b)

   Small b gives strong butt swell (a large d(0) relative to DBH). From the
   base radius r0 = d(0)/2 the trunk tapers linearly to zero at the tip:

   r(h) = r0 * (1 - h/H)

   This taper also defines the collision radii used by trimming (Section
   3.2) and the frustum volumes in the info JSONs.
2. **Crown width (allometry).** With DBH converted to inches (DBH_in):

   CW_ft = exp(A0 + A1 * ln(DBH_in)),  MAX_CAN_W (m) = 0.3048 * CW_ft

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
| A0, A1 | crown-width allometry (notes.md's a0/a1): CW (ft) = exp(A0 + A1 ln(DBH_in)), then to metres | widens the crown at a given DBH (A0) / makes big trees disproportionately wide (A1); wider crowns mean more collision trimming in dense plots |
| TAPER | butt-swell exponent of the taper function, Section 7.2 step 1 (notes.md: "taper function (base diameter)") | reduces butt swell (base diameter closer to DBH); lowering it thickens the trunk base |
| DBH_MIN, DBH_MAX, HT_MIN, HT_MAX | applicability bounds on the trees this entry covers (notes.md: bounds on the tree's DBH / HT). A tree's list DBH/HT must lie inside them; the realized SST ranges are derived around the tree list values | a tree outside the bounds falls through the code resolution order (species -> functional group -> DEFAULT); wider bounds admit more trees to this entry |
| SPLIT_PROB | probability the tree forks at all (0 = single straight stem) | more forked trees |
| SPLIT_NUM_LO/HI | number of cascading trunk splits, drawn per tree (default 1:1 = single fork); confirmed spec extension (owner ruling 2026-11-10), not in notes.md's parameter list | more levels of forking and bushier multi-leader crowns; each leader's crown is airier at the same point count, so raise the point budget for dense multi-leader crowns |
| SPLIT_HT_LO/HI | where each split happens, as a fraction of the stem being split | pushes splits higher up the stems |
| NUM_BRANCH_LO/HI | nominal branch count, drawn per tree and split between the lower crown band and the upper crown/fork band | denser canopy and more foliage points |
| MIN_CAN_HT_LO/HI | crown base: height of the lowest branch, as a fraction of H | raises the crown base (pruned-plantation look); more bare, visible stem |
| MAX_CAN_W_HT_LO/HI | height (fraction) where the crown is widest | moves the widest part of the crown upward |
| TOP_DIST | maximum horizontal wander of the stem tip, m (tree lean) | more leaning trees |
| MID_DIST | horizontal wander of the mid-stem control point, m (stem bow) | more bowed, curved stems |
| FP_RAND_HT | upward vertical jitter added to foliage points, m | fuzzier, less clearly bounded canopy |
| SPR_LO/HI | stem-point fraction of the point budget, drawn per tree | more points on trunk and branches, at the cost of canopy detail |

Useful silhouettes: a **conifer** look is low SPLIT_PROB, small A0/A1,
high crown base variation and high NUM_BRANCH; a **hardwood** look is high
SPLIT_PROB, high A0/A1, a low crown base, and SPLIT_NUM_LO/HI above
1:1 for cascading multi-leader crowns. TOP_DIST and MID_DIST only change
stem geometry, so they add variety without changing the crown.

**[implemented]**: keyed-by-code consumption at generation time — rows
are expanded from the tree list's DBH/HT, with allometry-derived crown
width (A0/A1) and taper-derived base diameter (TAPER). The dictionary is
serialized as JSON (`data/functional_groups.json`; legacy CSV tables
still load — pre-W18 spellings CW_A0/CW_A1/TAPER_B are renamed at load).
The per-code DBH_MIN/DBH_MAX/HT_MIN/HT_MAX applicability bounds are
enforced at entry resolution: a tree list DBH/HT outside an entry's
bounds falls through the resolution order (species -> functional group
-> DEFAULT), and DEFAULT always applies. Realized sizes are drawn within
the generation constant of +/-1% around the tree list values (the
DBH_TOL/HT_TOL dictionary keys and their per-code proportional
tolerances are gone from the schema).

### 7.4 SST skeleton JSON — **[implemented]**

Design target (per notes.md): a JSON where each element represents one
portion of the generated tree — an element graph:

| Attribute | Type | Description |
|---|---|---|
| BRANCH_UID | string | unique id of the element within the tree |
| BRANCH_ORDER | int | 0 = trunk, 1 = connected to the trunk, 2 = connected to a 1st order element, ... (5 is a soft upper limit). **-1 indicates the point does not belong to a branch** (non-branch points: ground, grass, shrubs; on cloud points this is `branch_code` 0 = order -1) |
| BRANCH_LENGTH | float | length of the element, m |
| BRANCH_RADIUS | float | radius of the element, m |
| BRANCH_VOLUME | float | cylinder volume pi * r^2 * L from the element's length and radius |
| SOURCE | string | UID of the element this element connects to (exactly one; null for the root trunk element) |
| DESTINATION | string | UID(s) of the elements connected to this one (multiple at a fork) |

Implementation (W15; split rule corrected in the 2026-11 review, item W22):
the graph is computed at generation time and
embedded in every skeleton JSON under the `elements` key (raw and trimmed
variants alike — the elements carry no coordinates, only UIDs/parents/
lengths/radii/volumes, so they are invariant under place() rotation and
reflection and the trimmed skeleton copies them verbatim). SST realization
decisions of the notes.md table:

- **Split rule (item 5 of the 2026-11 review):** when an element splits,
  one of the splits carries the parent order. The split daughter with the
  larger base radius keeps the parent's order; ties are broken by the
  daughter's own length, then deterministically (the first-listed
  daughter keeps the parent order). In SST the global taper gives both
  daughters of a split the same base radius (they start at the split
  height), so the larger-by-length daughter stays the trunk (order 0) and
  the shorter becomes order 1 — a fork of a fork can produce orders 2 and
  higher. A *branch* (a side element, not a split) is always
  parent order + 1, so branches hanging off an order-1 stem are order 2
  (SST branches never sub-branch).
- Partition: one element per trunk stem (the main stem is order 0; each
  fork daughter stem's order follows the split rule above) and one
  element per branch (source stem order + 1; SST branches do not
  sub-branch).
- UID scheme: `T00`, `T01`, ... for stems in stem order, `B000`, `B001`,
  ... for branches in skeleton order.
- SOURCE: the root trunk element's is null; each fork daughter stem's is
  the parent stem element (from the split log); each branch's is its
  recorded stem element (`stem_index`, recorded at generation).
- RADIUS: a stem element uses the taper radius of Section 7.2 step 1 at
  the element's base height; a branch element uses its recorded ribbon
  half-thickness. BRANCH_VOLUME is the notes.md cylinder; the frustum
  accounting of Section 8.4 is separate and unchanged.
- DESTINATION is the inverse of SOURCE (children in the same order).

The same order rule drives the point-cloud codes (assemble): a point on
an element carries `branch_code` = element order + 1, and a point on no
branch carries `branch_code` 0 = BRANCH_ORDER -1; stem points on an
order-0 element are semantic 1 (trunk), stem points on order >= 1
elements are semantic 2 (branch stem) — under the split rule the
dominant fork stem stays trunk.

The polyline-and-splits form below remains the trimming stage's native
format; consumers of the skeleton can use either.

- `tree_uid`, `plot_uid`, `tree_id`, `coordinate_frame`, `base_xy`,
  `height_m`, `base_radius_m`, `fork`, `fork_height_m`, `n_splits`,
  `splits`, `min_canopy_height_m`, `max_canopy_width_height_m`;
- `trunk`: centerline polylines per segment with per-segment z-range
  and length;
- `branches`: same organization for branch elements, each with its length
  and the branch element id that its points carry.

In the polyline form, branch connectivity (joint structure) is expressed
by the polylines and splits; there is no separate SOURCE/DESTINATION
table.

### 7.5 SST codebase modifications

- `treesim.py`: optional skeleton export (Section 7.4); a translation factor
  so tree tops can grow towards any point.
- `las.py`: canonical LAZ writer writing x/y/z plus the extra dims of
  Section 8.3 at their schema dtypes via extra-bytes records (LAZ
  compression via lazrs) as **LAS 1.4**, point format 0 — **[implemented]**
  (formats 0-5 remain valid in 1.4, so only the version tag changed;
  merged-file point counts are no longer uint32-capped). The added columns
  notes.md lists (`tree_ID`, `semantic`, `branch_ID`, `branch_order`) map
  to the Section 8.3 dimensions; `trimmed` is not implemented — trimmed
  points are removed and recorded structurally only (owner ruling
  2026-11-07, Section 3.2).

## 8. Data structures

### 8.1 Tree list (input; emitted by the FIA helper or authored by hand) — **[implemented]**

Serialized as **CSV**, one row per tree (per notes.md: the FIA helper's
placed tree list is a .csv and the parameter helper reads a tree-list
.csv). The JSON form — one object per plot (PLOT_UID, plot shape and
dimensions, optional plot center geo) with a `trees` array — remains an
accepted alternate serialization carrying the same attributes.
Required attributes:

| Attribute | Type | Description | Range | Notes |
|---|---|---|---|---|
| TREE_UID | string | unique tree id in the list | | ex: 0001, 0002 |
| X | float | local x of tree centroid, metres from lower-left | [0, PLOT_LENGTH] rect / [0, 2 PLOT_RADIUS] circle | |
| Y | float | local y of tree centroid, metres from lower-left | [0, PLOT_WIDTH] rect / [0, 2 PLOT_RADIUS] circle | |
| Z | float | local z of tree base, metres | [0, inf) | typically 0 |
| DBH | float | trunk diameter at breast height, metres | (MIN_DBH, inf) | |
| HT | float | tree height above ground, metres | (0, 120) | tallest known tree is ~116 m |
| TREE_CODE | string | vegetation type code (species, functional group) | | resolution order: species code -> functional group -> DEFAULT |

Optional attributes:

| Attribute | Type | Description | Range | Notes |
|---|---|---|---|---|
| LON | float | centroid longitude, EPSG:4326 | [-180.0, 180.0] | plot center + local placement; absent without a known plot center |
| LAT | float | centroid latitude, EPSG:4326 | [-90.0, 90.0] | plot center + local placement; absent without a known plot center |
| ELEV | float | base elevation, metres (EPSG:5773) | [-500, 5000] | plot center + local placement; absent without a known plot center |

**Status**: the attribute contract is implemented in both serializations
and every stage loads either through `tree_list.py`; DBH/HT carry full
inventory precision and CSV loads parse with
`float_precision="round_trip"` so CSV/JSON round-trips are bit-exact.
For CSV lists the plot-level fields come from the sibling
`plot_<UID>_report.json`; the JSON form carries them in the plot object.
Legacy SST parameter CSVs remain accepted as input for generation. The
parameter dictionary is JSON
(Section 8.6).

### 8.2 Canonical output layout

```
<PLOT_UID>/
  plot_placement/
    las/            raw / trimmed / terrain / understory clouds
    textual/        scene JSON, manifest, plot info JSON, DTM crop,
                    plot shrub index JSON (understory stage)
    rasters/        plot DTM GeoTIFF (EPSG:9707; anchored plots only)
    deep_learning/  stem_map / canopy_map feature+label rasters    (derived)
  individual_trees/<PLOT_UID>_<TREE_UID>/
    textual/        generation, skeleton and info JSONs  (no sub-clouds;
                    extract individuals from the plot cloud)
    deep_learning/  per-tree derived artifacts only            (derived)
  individual_shrubs/<PLOT_UID>_shrub_<SHRUB_UID>/
    textual/        shrub info JSONs         (understory stage)
```

Per-shrub volume labels under `deep_learning/volume/labels/` were removed
(2026-11 review, item W24; the per-shrub volume lives in the shrub info JSON
and the plot shrub index). Per-tree and plot-level derived artifacts under
`deep_learning/` (volume labels, stem/canopy map rasters) remain (owner
decision, 2026-11 review round).

### 8.3 Point cloud dimensions — **[implemented]**

x, y, z plus:

| Dimension | Dtype | Meaning |
|---|---|---|
| tree_id | int32 | 0 = not a tree; per-plot sequential otherwise |
| shrub_id | int32 | 0 = not a shrub |
| semantic_code | int8 | object identity, below |
| classification | uint8 | LAS class by height above ground: 2 ground, 3 low (<=0.5 m), 4 medium (0.5-2 m), 5 high (>2 m), 7 noise |
| branch_code | int8 | 0 not a tree/branch point (BRANCH_ORDER -1 per notes.md); 1 trunk (order 0); 2+ = element order + 1 (order 1 = connected to trunk, 3 = connected to a 2, ...) |
| branch_id | int32 | branch element identity within the plot (0 = unassigned) |
| height_above_ground | float32 | z minus DTM at (x, y) |

Semantic codes (per notes.md, owner ruling 2026-11-07): 0 ground, 1 trunk,
2 branch stem, 3 foliage, 8 shrub, 9 grass. Code 0 doubles as "no semantic
class needed" — ground and noise are the traditional LAS classification's
job (classification 2 ground, 7 noise); there are no separate unassigned
or noise semantic codes.

### 8.4 Info JSONs — **[implemented]**

Plot info (`<PLOT_UID>_info.json`), required attributes (per notes.md;
top-level keys, no nesting):

| Attribute | Type | Description |
|---|---|---|
| NUM_TREES | int | number of trees contained in the plot |
| NUM_SHRUBS | int | number of shrubs contained in the plot |
| BASAL_AREA | float | basal area of the plot's trees, m² (sum of pi * (DBH/2)²) |
| PLOT_SHAPE | string | "circle" or "rectangle" |
| PLOT_AREA | float | plot area, m² (pi * r² or length * width, from the scene JSON) |
| BOUNDING_BOX | 6 floats | [x_min, y_min, z_min, x_max, y_max, z_max] in local metres when the plot is unanchored; [lon_min, lat_min, elev_min, lon_max, lat_max, elev_max] in EPSG:4326/5773 when the plot is anchored (Section 10; elev = anchor elev + local z) |
| VOL_TREES | float | sum of every tree's volume |
| VOL_TRUNK | float | sum of trunk volumes (from the tree/skeleton JSONs) |
| VOL_BRANCH | float | sum of branch volumes (from the tree/skeleton JSONs) |
| VOL_FOLIAGE | float | sum of foliage volumes (voxelized at FOLIAGE_VOXEL) |
| VOL_SHRUB | float | sum of shrub volumes (from the per-shrub JSONs) |
| VOL_GRASS | float | volume of all grass points (voxelized at GRASS_VOXEL) |
| CHM_PATH | string | relative path to the generated canopy height map; null when the stage did not run |
| STEM_PATH | string | relative path to the generated stem map; null when the stage did not run |
| PLAUSIBILITREE_CONFIGURATION | string | relative path to the JSON recording the settings used (the scene JSON). notes.md's spelling is "PLAUSIBILITREE_CONFIGURAITION"; the corrected spelling is used (owner-approved 2026-11) |

Optional: PLOT_CENTER (lon, lat, elev in EPSG:4326/5773, present only
when the FIA input carried the pass-through columns) and DTM_PATH
(relative path to the plot DTM, present only when it was generated).

Tree info (per tree), required attributes (per notes.md): TREE_UID,
TREE_DBH (may differ from the tree list DBH), TREE_HT (may differ from the
tree list HT), TREE_X, TREE_Y (base position; match the tree list),
TREE_Z (centroid height above the tree's base: mean height-above-ground
of the tree's points from the understory cloud, falling back to the
centroid height above the lowest point of the tree cloud when no
understory cloud exists; varies from the tree list z), VOL_TREE (trunk +
branch + foliage), VOL_TRUNK, VOL_BRANCH, VOL_FOLIAGE, BOUNDING_BOX (6
floats, same form as the plot info). Optional: TREE_CENTER_REAL
[lon, lat, elev] at the tree's base in EPSG:4326/5773 — written when the
plot is anchored (Section 10; the bottom of the tree is the elevation
per notes.md), absent otherwise.

Fields beyond the notes.md set are retained as extensions: `meta`
(plot_uid / tree_uid, coordinate_system, units), `canopy` (per-tree crown
major/minor axes; the convex-hull volume extension was removed in the
2026-11 review, item 10 — it was only ever populated when the optional
canopy-map stage ran), `grass` (realized percent cover, plot info only),
`shrubs` (primitive-volume analytic reference, plot info only), and
per-tree / per-shrub `neighborhood` voxel-compaction records (extension
volumes measured at FOLIAGE_VOXEL so all class proportions share one
basis). The old tree_characteristics / vegetation_characteristics
nestings are gone.

Shrub info (notes.md defines no per-shrub info schema): meta
(shrub_uid / plot_uid / coordinate_system / units / bounding_box),
`shrub_characteristics` (center, max_width_m, height_of_max_width_m,
height_m, volume_m3 {primitive, mesh}) and `neighborhood`. The
2026-11 review (items 7/8) removed the ellipsoid `dbh_m` field and the
`volume_m3.voxel` key (shrub volume is never voxel-based per the
notes.md ruling); VOL_SHRUB is summed from the per-shrub info JSONs.

Volume methods: trunk = frustum sum along skeleton centerlines with taper
radii; branch = cylinder per branch (pi r^2 L); foliage = occupied-voxel
volume at FOLIAGE_VOXEL; grass at GRASS_VOXEL; shrub = sampled-shape
volume — the actual (jittered) primitive mesh that produced the shrub's
points (per the notes.md ruling, never occupied-voxel; the nominal
ellipsoid (4/3) pi (w/2)^2 (h/2) is the analytic reference). FOLIAGE_VOXEL and
GRASS_VOXEL are separate knobs (--foliage-voxel, --grass-voxel; both
default 0.1 per notes.md) with a shared occupancy threshold
(--voxel-threshold, default 1).

### 8.5 Scene JSON — **[implemented]**

One per plot; relative links to every artifact; plot bounds; per-tree
rotations and reflections; terrain record (DTM
source, cell size, generated-ness) sufficient to reproduce the DTM; topo
dict + vegetation settings.

### 8.6 Tree parameter dictionary

Module-specific mapping tree code -> parameters, serialized as **JSON**
— **[implemented]**: `data/functional_groups.json`, one object keyed by
tree code whose values are that code's parameters (the Section 7.3
columns: A0/A1 crown-width allometry, TAPER, the DBH_MIN/DBH_MAX /
HT_MIN/HT_MAX applicability bounds, and the shape parameters
SPLIT_HT_LO/HI, SPLIT_PROB, SPLIT_NUM_LO/HI (confirmed spec extension),
NUM_BRANCH_LO/HI, MIN_CAN_HT_LO/HI, MAX_CAN_W_HT_LO/HI, TOP_DIST,
MID_DIST, FP_RAND_HT, SPR_LO/HI). notes.md's required SST parameters
HT_LO/HT_HI, BASE_DIA_LO/BASE_DIA_HI and MAX_CAN_W are derived at
generation time from the tree list's DBH/HT (realized-size ranges within
the +/-1% generation constant; crown width from the A0/A1 allometry) —
they enter SST's model_params, never the per-code dictionary. The
Section 6.2 helper edits it; legacy CSV tables are still accepted
(pre-W18 spellings CW_A0/CW_A1/TAPER_B rename at load). Height and
diameter always originate in the tree list, never from this dictionary;
the dictionary supplies each code's applicability bounds, checked at
entry resolution (Section 7.3). The DBH_TOL/HT_TOL proportional
tolerance keys are withdrawn from the schema (W18 closed).

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
- Shared: the grass voxel (GRASS_VOXEL) is grass volume accounting;
  shrub volumes come from the sampled shape, not voxels.

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
| ROTATION_STEP | float | 30.0 | placement rotation quantization, degrees (notes.md's steps of N; 360 must be divisible by it) |
| BRANCH_RADIUS_FRAC | float | 0.10 | branch collision radius as a fraction of branch length |
| MIN_BRANCH_RADIUS | float | 0.02 | floor for the above, m |
| COLLISION_TOLERANCE | float | 0.0 | trim radius relaxation, m |
| STEM_RESOLUTION | float | 0.05 | stem map cell size, m |
| STEM_SECTION | float | 0.25 | buffer applied to DBH_HT for the breast-height stem map slab, m |
| CHM_RESOLUTION | float | 0.05 | canopy height map cell size, m |
| CANOPY_BUFFER | int | 2 | canopy raster growth padding, cells |
| FOLIAGE_VOXEL | float | 0.1 | voxel size for foliage volume, m |
| GRASS_VOXEL | float | 0.1 | voxel size for grass volume, m |
| OCCUPIED_VOXEL_THRESHOLD | int | 1 | points needed to count a voxel occupied |
| DTM_RESOLUTION | float | 0.25 | generated DTM cell size, m |
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

These surface as per-script CLI flags fanned out by run_pipeline; the
stem/canopy/DTM knobs carry the notes.md names (--stem-resolution,
--stem-section, --chm-resolution, --dtm-resolution; W19). build_info's
--occupied-voxel flags cover foliage and neighbor-class accounting until
W16 splits FOLIAGE_VOXEL / GRASS_VOXEL. Shrub volume
is never voxel-based (sampled-shape ruling, Section 8.4).

## 10. Coordinate reference policy — **[implemented]**

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

Implementation (item W5, `georef.py`):

- **Anchor at assembly.** The scene JSON carries a `geo` block (an
  explicit null when the plot carries no center — the un-georeferenced
  marker):
  epsg, utm_zone, utm_hemisphere, easting, northing, elev, local anchor
  point (always `[0, 0]` — placement centers every shape on the origin)
  and the source lat/lon. The UTM conversion is the `utm` library with
  the zone derived from the plot center (EPSG 326zz north / 327zz
  south); the limits of agreement are validated there and violations
  raise. The UTM system itself spans only 80 deg S .. 84 deg N — a plot
  center inside the ingest limits of agreement but outside that
  envelope (a polar plot) raises a plain `ValueError` at anchoring;
  polar plots cannot be georeferenced until the owner picks a polar
  CRS. `elev` is the EPSG:5773 height of the plot center.
- **Internals unchanged.** Local coordinates are metres offset from the
  anchor's UTM position, so placement, trimming, DTM sampling and every
  raster's pixel grid operate exactly as before — georeferencing only
  adds metadata at the edges (same-seed clouds are identical with and
  without a known plot center).
- **GeoTIFF exports** (stem map, canopy map, plot-level instance
  masks, per-tree / per-shrub crops): when anchored, written with an EPSG:4326 transform whose
  west/north corner comes from the anchor and whose pixel sizes are the
  metre cells converted to degrees at the plot latitude (x by
  1/cos(lat)); values (counts, masks, heights above ground) are
  unchanged. Unanchored plots keep local-metre transforms and no CRS.
- **DTM GeoTIFF**: a generated plot DTM is additionally exported as
  `plot_placement/rasters/<PLOT_UID>_dtm.tif` with the compound CRS
  EPSG:9707 (horizontal EPSG:4326, vertical NAVD88 height EPSG:5773;
  z = anchor elev + DTM height, so the exported band carries absolute
  heights). The ESRI ASCII `.asc` remains the internal working format;
  the GeoTIFF is written only for anchored plots (unanchored plots are
  un-georeferenced by design and keep the .asc only).
- **Info JSONs** (Section 8.4): when anchored, BOUNDING_BOX becomes
  [lon_min, lat_min, elev_min, lon_max, lat_max, elev_max] (elev = anchor
  + local z, EPSG:4326/5773), each tree gains TREE_CENTER_REAL
  [lon, lat, elev] at its base (plot center + (TREE_X, TREE_Y); bottom
  of tree for elev per notes.md), and PLOT_CENTER is written from the
  same anchor; unanchored output is local metres, exactly as before.

## 11. Code-to-design gap register

Work items to bring code into agreement with this document, in rough
dependency order. Status: **open** (not yet started), **closed** (code
matches the design).

| # | Item | Gap | Status |
|---|---|---|---|
| W1 | Minimal tree list contract | generate_trees consumes per-tree frozen SST ranges; switch to (X, Y, DBH, HT, TREE_CODE) + param-dict lookup at generation, drawing seeded shapes per TREE_UID. FIA helper emits the minimal list (DBH/HT per tree), keeping placement there | closed |
| W2 | FIA helper placement controls | (a) implement expansion factor modifier; (b) minimum tree distance (generalize --allowed-overlap); (c) max tries and edge effects (spec Section 6.1). Replace the placeholder TPA rounding | closed |
| W3 | Formal generator module interface | extract the Section 7.1 contract; move SST specifics behind it | closed |
| W4 | Flat-plane topography | implement topography dict None -> uniform plane ground sampling + warnings | closed |
| W5 | UTM internals + EPSG:4326 exports | convert EPSG:4326 inputs to UTM metres internally (`utm` library, zone from plot center; Section 10); export GeoTIFFs back to EPSG:4326/5773; metres-thinking preserved | closed |
| W6 | One plot cloud; derived stages opt-in | stop writing per-tree sub-clouds at every stage (plot cloud + per-tree/shrub JSONs are the core outputs); individual-tree extraction becomes an opt-in derived stage; stem/canopy maps and normalization become opt-in; stem map gains a second (mask) band (Sections 4-5) | closed (CHM band 1 now takes the max over vegetation codes 1/2/3/8 — grass 9 excluded as ground cover — and band 2 is the foliage mask; both maps' band contracts live in Section 5; per-tree DTM crops remain, map-gated) |
| W7 | Random reflection | add mirror reflection to placement draws | closed |
| W8 | HELIOS++ stage | planned; needs scanner parameter design and voxel export | open |
| W9 | Polygon mesh export | planned; per-element semantics/ids; feeds HELIOS++ polygon pathway | open |
| W10 | Downed woody debris | design and implement the DWD stage (logs, stumps) per Section 3.5 | open |
| W11 | Design defaults | align code defaults with Section 9: PLOT_RADIUS 7.315, DTM_CELL_SIZE 0.25, CANOPY_CELL_SIZE 0.05, STEM_CELL_SIZE 0.05 | closed |
| W12 | LAS 1.4 writer | switch the canonical writer to LAS version 1.4 (point format 0 + extra bytes unchanged; Section 7.5) | closed |
| W13 | JSON interchange | all tool/stage interchange becomes JSON: tree list in/out of every stage, FIA helper outputs, generator parameter dictionary (CSV -> JSON), editor import/export | closed (tree-list part reversed; see W17) |
| W14 | notes.md semantic codes | semantic codes per notes.md: 0 ground (doubling as "no class needed" — ground/noise are the traditional LAS classification's job, so no separate unassigned/noise codes), 8 shrub, 9 grass (Sections 3.4, 8.3). The per-point `trimmed` dimension was withdrawn from the design (owner ruling 2026-11-07): trimmed points are removed, structural record only | closed |
| W15 | Element-graph skeleton | export the notes.md per-element skeleton graph (BRANCH_UID/ORDER/LENGTH/RADIUS/VOLUME, SOURCE/DESTINATION) from the polyline form (Section 7.4) | closed (per-element trimmed flags on trimmed variants were to ride with W16, but the 2026-11-07 trim ruling removed trim flags entirely; trimmed skeletons subset the element graph with scrubbed DESTINATION and carry no trim records) |
| W16 | Info JSON schema per notes.md | plot info: NUM_TREES, BASAL_AREA, NUM_SHRUBS, VOL_TREES/TRUNK/BRANCH/FOLIAGE/SHRUB/GRASS, CHM_PATH/STEM_PATH, config-record path, optional PLOT_CENTER; tree info: TREE_UID/DBH/HT/X/Y/Z/VOL_* — adopt the notes.md field names over the current tree_characteristics/vegetation_characteristics nesting; split build_info's shared occupied-voxel accounting into FOLIAGE_VOXEL / GRASS_VOXEL (Section 8.4) | closed |
| W17 | CSV tree lists | tree lists are CSV-canonical per notes.md (FIA helper emits the placed tree list as .csv by default; JSON is the alternate); reverses W13's JSON-canonical tree lists (Sections 6.1, 8.1). The parameter dictionary stays JSON | closed |
| W18 | Parameter dictionary per notes.md | a0/a1 + taper naming; per-code DBH/HT applicability bounds (DBH_MIN/MAX, HT_MIN/MAX) replacing the proportional DBH_TOL/HT_TOL; confirm the SPLIT_NUM_LO/HI extension (Sections 7.3, 8.6) | closed |
| W20 | Parameter editor: dictionary from a tree list | derive a parameter dictionary from a tree list (propose per-code parameters from the list's species/DBH/HT) so the editor's load path covers tree lists as well (Section 6.2) | open |
| W19 | Rotation stepping + variable names | quantize placement rotation to steps of N (ROTATION_STEP, seeded); rename code variables to notes.md names: FOLIAGE_VOXEL, GRASS_VOXEL, DTM_RESOLUTION, CHM_RESOLUTION, STEM_RESOLUTION, STEM_SECTION (Sections 3.1, 9) | closed (FOLIAGE_VOXEL/GRASS_VOXEL split deferred to W16) |
| W21 | Normalization retired (2026-11 review item 2) | deep learning is out of scope: the normalization step happens on the DL side. Remove normalize.py, the --normalize pipeline stage and every normalized path key (Sections 5, 8.2) | closed |
| W22 | Branch-order split rule (2026-11 review items 4, 5) | -1 = a point on no branch (branch_code 0); at a split the daughter with the larger base radius keeps the parent order (ties: longer, then first-listed), the other is parent + 1; branches are parent + 1 — applied to both the element graph and the per-point codes (Section 7.4) | closed |
| W23 | Shrubs inside the plot footprint (2026-11 review items 6, 7) | circular plots seeded shrubs over the bounding square, and the disc clip removed whole corner shrubs while their records persisted (0-point shrubs, null bb/neighborhood). Shrub cells are constrained to the plot disc; a shrub with no surviving points is dropped with a warning (Section 3.4) | closed |
| W24 | Shrub volume label + shrub info cleanup (2026-11 review items 7, 8) | remove per-shrub deep_learning volume labels (volume lives in the shrub info JSON and the plot shrub index); drop dbh_m and the always-null volume_m3.voxel key; VOL_SHRUB is summed from the per-shrub info JSONs (Section 8.4) | closed |
| W25 | Tree info convex hull + realized DBH (2026-11 review items 9, 10) | generation.json's realized record carries the taper-implied DBH; the per-tree canopy convex-hull volume extension (null unless the optional canopy-map stage ran) is removed (Sections 7.1, 8.4) | closed |
| W26 | Tree skeleton viewer helper (2026-11 review item 11) | 3D GUI skeleton JSON viewer: radius/ID/angle toggles, branch-order colorization, zoom/rotate, render to png/jpeg (Section 6.5). Full implementation TBD | open |
| W27 | Review-round housekeeping (2026-11 review items 1, 3) | pin utm in requirements.txt (notes.md comments say the pipeline uses the utm library); drop the redundant --basis stem (it is the default) from quickstart instructions | closed |

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
  [Superseded by the W18 closing: the tolerance keys are withdrawn from
  the dictionary schema; the +/-1% is a generation-time constant, and the
  dictionary carries DBH/HT applicability bounds instead.]
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
- Third code-alignment round (2026-11): W3 closed — the Section 7.1 module
  contract is formalized in generator_module.py with an SST implementation
  (sst_generator.py) behind it; the core builds through build(rng, ...)
  with base-origin TreeModels and module place() positioning; the rng is
  the sole entropy source (treesim's global generator bridged through one
  rng draw); both input contracts drive the same build path.
- Fourth code-alignment round (2026-11): W4 closed — the None topography
  is a flat plane at z = 0 (terrain always runs; ground sampled at the
  default density with a spec warning when unconfigured); conflicts
  resolve None > filepath > random with a warning; the terrain stage's
  standalone gen-cell default aligned with DTM_CELL_SIZE 0.25.
- Fifth code-alignment round (2026-11): W13 closed — the parameter
  dictionary is JSON (data/functional_groups.json, one object keyed by
  tree code; legacy CSV tables still load) and the parameter editor
  imports/exports JSON; interchange is JSON end to end with CSV as a
  legacy input form only.
- Design-alignment ruling (2026-11-07): notes.md is the upstream design
  document and overrides this spec wherever they disagree. Reversed by
  that ruling: semantic codes return to 0 ground / 8 shrub / 9 grass (the
  earlier 4/5/6 decision is itself withdrawn; superseded by the W14 closing
  rulings — 0 doubles as "no class needed", no separate unassigned/noise
  codes); the per-point `trimmed` dimension is reinstated (the earlier
  structural-record-only decision is withdrawn — itself superseded by the
  W14 closing ruling: no trimmed dimension, structural record only, W14);
  placed tree lists are
  CSV again, with the parameter dictionary remaining JSON (W17); the
  parameter dictionary carries per-code DBH/HT applicability bounds and the
  a0/a1/taper form rather than proportional tolerances (W18); placement
  rotation is quantized in steps of N (W19); the skeleton JSON design
  target is the per-element graph with SOURCE/DESTINATION and cylinder
  volumes (W15). Section 8.4 now carries notes.md's plot/tree info schemas
  (W16).
- Sixth code-alignment round (2026-11-07): W17 closed — the FIA helper
  emits the minimal tree-list CSV by default (tree rows only; plot geometry
  and plot center recorded in the report sidecar, read back by tree_list),
  every stage loads both serializations through tree_list.py, and CSV
  loads parse with round_trip float precision so CSV/JSON round-trips are
  bit-exact; JSON is the alternate form, the SST table the legacy one.
- Owner rulings for W19 (2026-11-07): placement rotation is quantized in
  steps of N degrees with ROTATION_STEP defaulting to 30 (12
  orientations); shrub volume is calculated from the sampled shape (the
  actual jittered primitive mesh, nominal ellipsoid as the analytic
  reference) — shrub occupied-voxel accounting is withdrawn.
- Seventh code-alignment round (2026-11-07): W19 closed — rotations are
  drawn on the ROTATION_STEP grid (asserted to divide 360; scene JSON
  records rotation_step_deg; standing suite check guarantees the grid),
  the stem/canopy/DTM knobs carry the notes.md names end to end
  (--stem-resolution, --stem-section, --chm-resolution, --dtm-resolution,
  env-sidecar key dtm_resolution), and shrub volume JSONs carry
  primitive + mesh only. W16 inherits the FOLIAGE_VOXEL / GRASS_VOXEL
  split of build_info's shared occupied-voxel accounting.
- Owner rulings closing W14 (2026-11-07): semantic code 0 = ground and
  doubles as "no semantic class needed" — ground and noise are the
  traditional LAS classification's job, so there are no separate
  unassigned or noise semantic codes (the earlier interim proposal of
  255/7 is withdrawn; 255 was invalid for int8 anyway). The per-point
  `trimmed` dimension is out of the design for good: trimmed points are
  removed and recorded structurally (skeleton node_trimmed flags, scene
  JSON counts) only. W14 closed with codes 0/1/2/3/8/9 in code and the
  standing suites asserting the mapping.
- Eighth code-alignment round (2026-11-07): W15 closed — every skeleton
  JSON (raw and trimmed) carries the notes.md element graph under
  `elements`: one order-0 element per trunk stem (root SOURCE null;
  split daughters source their parent stem), one order-1 element per
  branch (sourced to its stem via `stem_index`, recorded at generation),
  cylinder volumes and UIDs per the Section 7.4 decisions (treesim also
  records each branch ribbon's half-thickness as `radius_m`; both
  additions are rng-neutral — same-seed clouds are bit-identical). The
  elements carry no coordinates, so the trimmed variant copies them
  verbatim. Standing suite section J re-derives the graph independently
  of the writer (root/reciprocity/partition/lengths/radii/volumes) and
  cascade-forced builds verified granddaughter chains.
- Owner ruling, trim records (2026-11-07): once a point or branch has
  been trimmed, no record of it is kept — the per-point dimension, the
  skeleton node/element trim flags, the trim summary, the scene-JSON
  collision log and removal counts, and all raw-vs-trimmed volume pairs
  were removed. The trimmed skeleton is now physically the surviving
  geometry (polylines cut at removed sub-segments, trimmed branches
  dropped, element graph subsetted), build_info computes single
  post-trim volumes, and the standing suites gained a forced-collision
  section proving record absence. The raw/trimmed cloud pair remains a
  pre-/post-trim stage variant until W6; whether it survives is W6's
  call. A latent bug also fixed here: the plot info's shrub totals read
  a shrub_index.json no stage writes (always zero); they now sum the
  per-shrub mesh volumes.
- Ninth code-alignment round (2026-11): W16 closed — the plot and tree
  info JSONs carry the notes.md field names as top-level keys (NUM_TREES,
  BASAL_AREA, NUM_SHRUBS, VOL_TREES/TRUNK/BRANCH/FOLIAGE/SHRUB/GRASS,
  PLOT_SHAPE, PLOT_AREA, BOUNDING_BOX as 6 floats, CHM_PATH/STEM_PATH,
  PLAUSIBILITREE_CONFIGURATION; per-tree TREE_UID/DBH/HT/X/Y/Z and
  VOL_TREE/TRUNK/BRANCH/FOLIAGE), the extension fields moved out of the
  old tree_characteristics/vegetation_characteristics nestings into
  meta/canopy/grass/shrubs/neighborhood blocks, and build_info's shared
  occupied-voxel knob split into FOLIAGE_VOXEL and GRASS_VOXEL
  (--foliage-voxel/--grass-voxel, shared --voxel-threshold; neighborhood
  extension proportions measure at FOLIAGE_VOXEL). TREE_Z = mean
  height-above-ground of the tree's points (understory cloud), with the
  tree-cloud fallback documented in Section 8.4. Owner-approved
  deviation: the config-record field is spelled
  PLAUSIBILITREE_CONFIGURATION (notes.md's "CONFIGURAITION" is a typo).
  Optional DTM_PATH rides along when a plot DTM exists.
- Tenth code-alignment round (2026-11-10): W18 closed — the parameter
  dictionary carries the notes.md names (A0/A1 crown-width allometry,
  TAPER; pre-W18 spellings CW_A0/CW_A1/TAPER_B rename at load, so legacy
  CSV/JSON files still open), every entry gained DBH_MIN/DBH_MAX and
  HT_MIN/HT_MAX applicability bounds enforced at entry resolution (a
  tree list DBH/HT outside an entry's bounds falls through species ->
  functional group -> DEFAULT, with DEFAULT the unconditional last
  resort), and the DBH_TOL/HT_TOL tolerance keys are withdrawn from the
  schema (realized sizes keep the +/-1% generation constant, now fixed
  in sst_generator rather than per-entry). Owner rulings: the taper
  parameter is named TAPER; SPLIT_NUM_LO/HI is a confirmed spec
  extension. The rename is numerically inert — a full demo regeneration
  produced byte-identical clouds, rasters and skeletons (only the
  generation JSON timestamps moved).
- Eleventh code-alignment round (2026-11-10): W5 closed — georef.py owns
  the plot anchor (utm library, zone from the plot center; EPSG 326zz
  north / 327zz south; limits of agreement validated), the scene JSON
  carries the geo block (an explicit null for unanchored plots), and
  internal geometry stays in local metres (offsets from the anchor's UTM
  position, so placement/trimming/DTM/raster grids are unchanged and
  same-seed clouds are identical either way). GeoTIFF exports (stem map,
  canopy map, per-tree and per-shrub crops) convert to EPSG:4326 through
  the anchor with degree pixel sizes at the plot latitude (pixel grids
  and values unchanged); the generated plot DTM gains a
  plot_placement/rasters GeoTIFF with compound CRS EPSG:9707 and
  absolute heights; info JSONs write geographic BOUNDING_BOX values and
  TREE_CENTER_REAL when anchored, and PLOT_CENTER from the same anchor.
  Standing suite section K runs an anchored mini plot through the whole
  pipeline (UTM anchor, raster CRS/degree cells, compound DTM, geo info
  JSONs, round-trip precision, limits of agreement) and asserts the
  unanchored fixture stays CRS-less. The committed demo regenerated
  byte-identically except its scene geo:null marker.
