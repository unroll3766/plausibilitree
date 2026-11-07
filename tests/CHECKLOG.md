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
