# Check Log

Running log of unique checks, per the global agent rules. One to three
lines per entry: date, what was verified, command/method, actual result.
Re-entries only when behavior changes. Long-form narrative lives in
`VERIFICATION_LOG.md`; this file is the short pool.

- 2026-11-07 — FIA tree-list CSV/JSON equivalence (W17): `fia_helper.py`
  twice on CA_TREE_SMALL.csv (seed 42, r 20) in both formats, then
  `assert_frame_equal` on loaded trees + plot fields via report sidecar —
  exit 0, 111 trees, geometry+center identical; CSV column order differs
  (Z canonical position), content identical.
- 2026-11-07 — CSV float round-trip must be exact: pandas' default fast
  parser was 1 ULP off on 17-sig-digit values (DBH 46 in -> 1.1684 vs
  1.1683999999999999); loader fixed with `float_precision="round_trip"`,
  proven by exact `==` of all 111 `tree_requests` from CSV vs JSON — pass.
- 2026-11-07 — W17 regression run: `conda run -n ssf python
  tests/run_all.py` — exit 0, all suites passed.
- 2026-11-07 — W19 rotation grid: `assemble_plot` across seeds 2000-2007 —
  every recorded rotation a multiple of 30 deg in [0,360), all 12
  orientations occur; promoted into the standing suite (placement §7).
- 2026-11-07 — W19 regression: check.py exit 0; tests/run_all.py exit 0
  (both suites, renamed flags exercised end-to-end: stem/canopy maps at
  5 cm, terrain DTM 0.25, 45-shrub understory).
- 2026-11-07 — Shrub volume JSON schema: per-shrub *_volume.json after the
  invariant-suite run (--keep) carries volume_m3 {primitive, mesh} only —
  voxel key and voxel_volume_params absent (sampled-shape rule).
- 2026-11-07 — W19 review pass: rotation-grid assertion (25 deg rejected
  with the explicit message, 90 deg accepted); stem/canopy/terrain CLIs
  list the new flags; understory CLI has no occupied-voxel refs; replicate
  sidecar round-trip executed (dtm_resolution key) through a full
  run_pipeline on a legacy SST replicate (1 ok, exit 0); scene JSON
  records rotation_step_deg 30.0 with rotations on the grid.
- 2026-11-07 — W14 semantic codes (notes.md values): schema 0 ground
  (doubles as "no class needed"; no separate unassigned/noise codes —
  LAS classification handles them), 8 shrub, 9 grass; check.py exit 0;
  tests/run_all.py exit 0 with new standing invariants (understory
  domain subset {0,1,2,3,8,9}, shrubs=8, grass=9, ground=0); observed
  plot cloud reads semantic_code=[0,1,2,3,8,9].
- 2026-11-07 — W14 review pass: found+fixed two stale-code bugs the suites
  don't cover — colorize_plot (ground 6/grass 5/shrub 4 masks; now schema
  constants) and scene_editor SEM_COLORS/SEM_NAMES (4/5/6 keys; now
  schema constants). Colorize executed on a kept fixture: grass 61,263
  grey, ground 21,357 black, woody brown, shrub/foliage randomized
  (101,578 pts). re-run tests/run_all.py exit 0 (119 PASS).
- 2026-11-07 — W15 element graph: suite section J (60 checks) re-derives
  the graph from the polyline skeleton independent of the writer
  (root/reciprocity/partition/lengths/radii/volumes; fixture includes a
  fork tree) — exit 0; cascade-forced builds (SPLIT_PROB=1, num_splits 2
  and 3: 5- and 7-stem trees, granddaughter chains) hold reciprocity and
  acyclicity; run_all.py exit 0, 179 PASS total, same-seed fixture values
  unchanged (rng-neutral treesim additions).
- 2026-11-07 — Trim-record removal (owner ruling): placement suite section 8
  forces heavy collisions (8,430 -> 2,067 pts) — no collision log/removal
  counts in scene, no node_trimmed/trim_summary/"trimmed" in trimmed
  skeletons, surviving graph matches raw element-wise (DESTINATION subset),
  branch lists shrink — all PASS; cross suite J restated to subset
  semantics; colliding two-tree plot through the full pipeline: info
  volumes single-valued, scene clean (exit 0); shrub totals check added
  (n>0, mesh volume>0) pinning the shrub_index bug fix; run_all 210 PASS.
- 2026-11-07 — Trim-removal review pass: consumer grep clean (no code reads
  the removed fields; only the suite's ban-list names remain); cut-geometry
  probe on 4 forced-collision trees (seed 5, frac 45): every surviving
  polyline run is a contiguous span of the raw stem's edge lengths
  (tree 001 trunk 52.93 -> 50.60 m across 7 runs; graph 174 -> 7),
  subsets internally consistent; branch_order las_field doc drift fixed
  (real field is branch_code); re-run run_all.py exit 0, 210 PASS.
- 2026-11-10 — W16 info-JSON schema: build_info emits notes.md top-level
  fields (plot NUM_TREES/BASAL_AREA/NUM_SHRUBS/VOL_*/PLOT_SHAPE/PLOT_AREA/
  CHM_PATH/STEM_PATH/PLAUSIBILITREE_CONFIGURATION; tree TREE_UID/DBH/HT/
  X/Y/Z/VOL_*), BOUNDING_BOX as 6 floats, FOLIAGE_VOXEL/GRASS_VOXEL split
  (--foliage-voxel/--grass-voxel/--voxel-threshold, run_pipeline wired) —
  run_all.py exit 0 (221 PASS, 0 FAIL); invariant suite re-derives
  BASAL_AREA, VOL_TREES, VOL_SHRUB (per-shrub mesh sum), PLOT_AREA and
  resolves CHM/STEM/config paths.
- 2026-11-10 — W16 review pass: consumer grep clean (no code reads the old
  tree_characteristics/vegetation_characteristics nestings or dict-style
  bbox from info JSONs; terrain/understory reads are canopy_bbox_json's own
  schema); the committed examples/plot_demo still predated W14/W19/trim
  removal (old semantic codes 4/5/6, raw/trimmed volume pairs) —
  regenerated with the README command (seed 42, exit 0): normalize reports
  codes [0,1,2,3,8,9], no trim-record markers in scene/skeletons (only the
  raw/trimmed cloud-pair paths), VOL_SHRUB == per-shrub mesh sum, plot
  info sums == per-tree values, paths resolve; new flags on build_info
  --help.
- 2026-11-10 — W18 parameter dictionary: data/functional_groups.json uses
  A0/A1/TAPER + DBH/HT applicability bounds; resolution falls through on
  out-of-bounds sizes (DEFAULT last resort); tolerance keys withdrawn;
  legacy CSV/JSON rename at load. run_all.py exit 0 (233 PASS); invariant
  suite section Hb asserts the schema + fall-through. Demo regen with the
  same seed: byte-identical LAZ/TIF/skeleton/scene/info artifacts — rename
  numerically inert (only generation.json timestamps moved, reverted).
- 2026-11-10 — W18 review pass: resolution probes (202→DF in bounds; out of
  bounds → DEFAULT; key-only ignores bounds; unknown → DEFAULT); legacy CSV
  with CW_* / TAPER_B and no bounds renames to A0/A1/TAPER and applies
  unlimited; param_editor headless load+sample works (A0 allometry,
  TAPER extrapolation); fia_helper.crown_width vector fallback to DEFAULT.
  Found+fixed: entry_for was missing from the GeneratorModule contract
  (core calls it; only SST implemented it) — added to the abstract class
  and spec 7.1. check.py + run_all.py exit 0 (233 PASS).
- 2026-11-10 — Full-repo pass: git clean (main ahead 9, nothing pushed);
  check.py + run_all.py exit 0 (233 PASS). Spec section markers vs
  register cross-checked: 7.4's [partial] was stale (W15 closed it) —
  fixed; 6.2's unimplemented tree-list derivation was untracked — added
  register item W20. README-vs-CLI flag sweep across 12 scripts: one
  genuine drift (README_understory described the retired
  --shrub-width-cell-frac width model) — fixed; other hits were
  cross-script references or explicit "replaced/does not exist" notes.
- 2026-11-10 — W5 georeferencing: anchored mini plot through run_pipeline
  (plot_center in a JSON tree list): scene geo block (EPSG 32610,
  zone from center), stem/canopy/shrub GeoTIFFs in EPSG:4326 with degree
  pixel sizes, DTM GeoTIFF compound EPSG:9707 with absolute heights, info
  JSONs geo bbox + TREE_CENTER_REAL + PLOT_CENTER, unanchored fixture/demo
  rasters stay CRS-less and skip the DTM GeoTIFF; run_all.py exit 0
  (260 PASS, suite section K; two bugs found en route: dtm_tif layout key
  landed in _TREE (KeyError in terrain), and load_anchor missed the
  pre-restructure scene location — stem/canopy maps run before
  restructure files the scene).

2026-11-11 — W5 review: restructure's plot-level instance masks were built from the
  anchored feature tif's degree transform (metre geometry vs degree grid -> wrong labels,
  no/double CRS). Fix: georef.local_grid_from_geo inverts the geo transform; masks compute
  in local metres, write EPSG:4326 with anchor; re-run smoke: labels 1-3 present, mask
  transform within 8e-10 deg of the feature raster.
2026-11-11 — W5 review: polar plot centers (lat outside UTM's -80..84) now raise a plain
  ValueError at anchoring (utm lib would die with its own OutOfRangeError past the ingest
  limits of agreement); -80/84 edges inclusive. Probe: conda run -n ssf georef edge
  battery (sydney anchor EPSG 32756/S, no-elev anchor, limit edges, scene round-trip).
2026-11-11 — W5 review: CSV tree list + plot_<UID>_report.json sidecar plot_center anchors
  end-to-end (geo block, PLOT_CENTER, BOUNDING_BOX all set; /tmp smoke run).
2026-11-11 — W5 review: suite section K asserts both label masks are EPSG:4326, transform
  within 1e-6 deg of the feature raster, and label all three trees; test suite exit 0.
2026-11-11 — W5 review: run_all.py exit 0, 266 PASS 0 FAIL; demo regen (seed 42) differs
  from examples/plot_demo only in the four generated_utc timestamps (byte-identical
  clouds/rasters/JSONs).

2026-11-11 — W6: stem map and CHM are 2-band per notes.md (stem band 2 = tree-point
  mask in the slab, uint16; CHM band 1 = max vegetation height (codes 1/2/3/8,
  grass 9 excluded as ground cover), band 2 = foliage mask, float32); smoke-run +
  suite section K asserts count/dtypes/band-2-implies-band-1; run_all exit 0.
2026-11-11 — W6: derived stages opt-in (--stem-map/--canopy-map/--normalize,
  default off). Core-only smoke: no map stages, no deep_learning dirs, scene
  without map links/blocks, CHM_PATH/STEM_PATH null, no normalized cloud.
2026-11-11 — W6: core run writes zero per-object clouds (generate/assemble flat
  staging deleted by restructure; terrain/understory/shrub_outputs per-object
  crops removed); build_info per-tree geometry = tree_id slices of the plot
  trimmed/raw clouds; extracted crops equal the slices (allclose 1e-9; float
  noise ~4e-15 from LAS offset arithmetic), suite sections core + L.
2026-11-11 — W6: extract_trees crops by tree_id/shrub_id from the resolved plot
  cloud, links into the scene (trees + new shrubs section); --surroundings
  window verified (crop gains neighbouring tree ids); shrub normalized clouds
  flow through normalize after extraction (demo: 7 normalized clouds).
2026-11-11 — W6: demo regenerated with the full derived set (--stem-map
  --canopy-map --normalize --extract-trees); second regen byte-identical
  except the four generated_utc timestamps; run_all.py exit 0.

2026-11-11 — W6 review: extraction battery — core-only extract resolves terrain
  (chain order), re-extract idempotent (byte-identical), missing plot cloud raises
  FileNotFoundError, --surroundings crops carry background tree ids; committed demo
  scene link audit: every link resolves (7 plot / 28 tree / 2 shrub).
2026-11-11 — W6 review: build_info slice values recompute independently from the
  plot cloud (tree 001 bbox exact, VOL_FOLIAGE 2.044 == 2.044); placement clouds
  carry only codes 1/2/3 so the CHM vegetation filter is inert there (probed).
2026-11-11 — W6 review: fixed terrain's "N tree DTM crops" counter (counted trees
  with canopy bboxes even when no DTM was written); widened the core-only
  no-per-object-clouds check to every individual_*/**/*.la[sz]; spec's "extraction
  is the only writer" claim corrected (normalize writes per-object _normalized
  variants of extracted crops); run_all exit 0.

2026-11-11 — Review round (W21): normalize.py retired with --normalize and every
  normalized path key; check.py parses all 29 modules, run_all.py exit 0 with no
  normalized artifacts in fixtures.
2026-11-11 — Review round (W22): fork-order split rule (dominant daughter by
  taper radius, ties by stem length, keeps parent order; branches parent + 1)
  in both element graphs and per-point branch_code; demo hardwoods now show
  orders 0..3; run_all.py exit 0.
2026-11-11 — Review round (W23): shrub zero-point diagnosis confirmed by
  synthetic repro (37/225 bounding-square corner shrubs lost every point to the
  circular-plot disc clip while records persisted); after footprint-constrained
  cell placement the same repro yields 0 zero-point shrubs.
2026-11-11 — Review round (W24/25): demo contracts verified by script — realised
  dbh_m present in generation JSONs, no dbh_m / voxel key in shrub info, no
  per-shrub volume JSON, no convex_hull in tree info, VOL_SHRUB == sum of
  per-shrub info mesh volumes, no normalized artifacts (demo regen, examples,
  exit 0).
2026-11-11 — Review round: ground-on-DTM invariant bound now derived from the
  DTM's own max slope (scale*slope + scale = 4.6 mm for the fixture) instead of
  the fixture-specific 2.5 mm; the shrub-placement change shifted the cloud
  bbox (LAS offset grid) and tripped the old bound (full suite exit 0 after).
2026-11-11 — Review-round audit: the two claims not yet probed directly now hold —
  kept full-pipeline fixture contains 0 *_normalized* artifacts (find,
  suite --keep); demo regen in a fresh dir is byte-identical to the committed
  examples/plot_demo except the four generated_utc timestamps (diff -rq).
  The normalized-artifact check was reinstated as a standing invariant
  (stronger than the deleted core-only one: any *_normalized* anywhere in the
  core-only or full run fails); suite exit 0 with both new checks PASS.
2026-11-11 — Per-commit historical audit (doc-first round): diff review of all
  10 commits vs messages — faithful, no strays. Worktree state runs: check.py +
  placement suite pass at 2d7dfe8/f70da55/280fb6b; full run_all exit 0 at
  f9a45a3. Known gap: commits 2-4 cannot run the cross suite (fixture.py still
  passed normalize=True until f9a45a3 — TypeError). Demo commit deleted all
  7 *_normalized.laz + 2 shrub volume JSONs; 0 *_normalized* tracked at HEAD.
2026-11-11 — Full end-to-end repo pass (post review round): check.py (29
  modules) + run_all.py exit 0; git clean, no pyc/scratch/secrets tracked,
  largest tracked file 252 KB; retired formats (normalized, per-shrub volume
  labels, dbh_m/voxel keys, tree-info hull) have zero live readers; demo
  scene link audit 44/44 resolve (links are relative to the plot dir's
  parent); colorize_plot smoke on the demo understory cloud OK (26,988
  points); demo plot/tree/shrub info keys match spec 8.4 exactly (canopy
  axes only, no dbh_m/voxel, DTM_PATH optional present). Doc gaps found and
  fixed: --extract-source/--extract-surroundings were undocumented (W6-era);
  the W23 shrub footprint placement was undocumented in README_understory.
2026-11-11 — Promoted the review-round probes into a standing third suite
  (tests/test_review_round_regression.py, in run_all.py): shrub footprint/zero-
  point invariants (W23), synthetic fork-order battery incl. radius-first
  dominance, ties and subordinate cascades (W22), and the committed-demo
  contract (links, retired artifacts, item 9/10/7 schemas). 21 checks,
  ~10 s, exit 0; mutation-tested (inverted dominance -> 5 FAIL, disabled disc
  filter -> 3 FAIL, both reverted); run_all.py exit 0 with three suites.
