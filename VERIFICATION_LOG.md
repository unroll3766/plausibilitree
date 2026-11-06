# Verification Log

Running record of the tests and checks used to verify code changes against
the design spec (`Plausibilitree.md`). Append a dated entry for every
verification round; keep entries copy-pasteable so any check can be re-run.

## Environment

- Interpreter: the `ssf` conda env (numpy, pandas, laspy 2.7 with lazrs,
  trimesh 5.1, rasterio, scipy).
  Run scripts with: `conda run -n ssf python <script>`
  (Directly invoking the env's python binary as a standalone command is
  blocked in this environment; always wrap with `conda run`.)
- Repo check: `conda run -n ssf python check.py` — parses every module,
  reports undefined names. Fast static sanity check; run after every edit.
- Sample data: `CA_TREE_SMALL.csv` (one plot, 111 documented trees,
  `TPA_UNADJ` = 1.0 everywhere, so expansion only engages at radii > ~32 m
  or with modifiers; see entries below).

## Checks

### Standing checks (run every round)

1. `check.py` — static, whole-repo: `conda run -n ssf python check.py`.
   (The runtime blocks `python3` as a standalone command; run every check
   through the conda wrapper.)
2. Stale identifiers: `grep -rn "allowed_overlap\|allowed-overlap\|7\.32" *.py *.md`
   (expect: no hits outside intentional provenance mentions); during the
   W3/W1 rounds also grep for API names removed by refactors
   (`load_rows`, `rows_from_tree_list`, ...).
3. Doc cross-references resolve:
   `conda run -n ssf python -c "import re; t=open('Plausibilitree.md').read(); refs=set(re.findall(r'Section (\d+\.\d+)', t)); heads=set(re.findall(r'^### (\d+\.\d+)', t, re.M)); print([r for r in refs if r not in heads])"`
   (expect: empty list)
4. Exit-code discipline: NEVER pipe a verification command through
   `tail`/`grep` without checking its exit code — `cmd | tail -2` masks
   a nonzero exit with tail's 0. Either run unfiltered, add
   `echo "exit: $?"`, or use `set -o pipefail`. (Added after the W3
   review caught a regression hidden this way; see the W3 review entry.)

### 2026-11-05 — W2 / W11 / W12 (commit e1868e6)

All snippets below assume `import sys; sys.path.insert(0, repo_dir)` if run
from outside the repo. Result is noted per check; all passed.

**FIA CLI smoke** (W2 end-to-end, plot CSV + report JSON):

```
conda run -n ssf python fia_helper.py CA_TREE_SMALL.csv -o /tmp/fia_test -r 20 --seed 42 --limit 3
```

Result: `plot_11839.csv` (111 trees placed, 0 not placed) +
`plot_11839_report.json` with parameters, trees_expanded/placed/not_placed,
proportion. Report keys verified by eye against spec Section 6.1 step 8.

**Expansion modifier semantics** (`expand_by_tpa`, spec 6.1 step 5), on
`CA_TREE_SMALL.csv` plot with `plot_radius_m=100` (expected/tree = 7.74:

```python
base = fia_to_sst_plot(sub, plot_radius_m=100, seed=42, return_intermediate=True)[1]
mod0 = fia_to_sst_plot(sub, plot_radius_m=100, expansion_modifier=0.0, ...)
mod05 = ... expansion_modifier=0.5 ... ; mod2 = ... expansion_modifier=2.0 ...
assert: mod0 = 111 (documented only), mod05 = 444, base = 888, mod2 = 1776
```

Result: 111 / 444 / 888 / 1776 — matches ceil(TPA x area x modifier) with
floor 1. (At r=20 all modifiers give 111: expected/tree = 0.31 < 1 floors.
This is correct, not a bug — verified before trusting later results.)

**Placement ordering invariant** (spec 6.1 step 7): on the placed frame,
`placed["HT"].diff() > 0` count must be 0 (height non-increasing).
Result: 0 violations.

**Min-distance spacing invariant**: `fia_to_sst_plot(..., radius 100,
min_distance=1.0, overlap_basis="stem")`; pairwise surface-to-surface gap
`hypot(dxy) - (DIA_i + DIA_j)/2` over all placed pairs.
Result: min gap 1.029 >= 1.0. (Negative min_distance sanity: crown basis,
min_distance -2 places more trees than +2.)

**Edge-effects quota**: disc-outside-proportion for every placed disc must
be <= quota + eps. Cases: circle r=40 edge=0.0 (worst = 0.0, strict
containment), edge=0.4 (worst <= 0.4); circle r=100 edge=0.25 shows 3 discs
legitimately straddling the boundary; rectangle 90x90 edge=0.3 (worst 0.000).
Result: invariant held in every case.

**Geo math cross-check**: analytic `_disc_outside_proportion` vs 4096x4096
midpoint grid, circle plot, 6 configurations (origin, edge, partial straddle,
fully outside, disc > plot, off-axis).
Result: agreement to 5 decimal places on all 6. Note: a disc centered exactly
on the boundary is ~0.51 outside, NOT 0.5 — a convex plot boundary curves
inward around the disc; the analytic and grid agree on 0.52124.
(`_disc_outside_proportion` rectangles use a 64x64 midpoint grid.)

**Same-seed reproducibility**: two `fia_to_sst_plot` runs, identical args
(r=40, seed=123, edge=0.2), `pd.testing.assert_frame_equal`.
Result: identical. Caution: the RNG stream changed vs the pre-W2 code
(ordering tiebreak draws precede placement draws), so plots generated with
old code will not regenerate identically from the same seed.

**LAS 1.4 roundtrip** (W12): `export_canonical_las` a 50-point cloud with
PLOT_FIELDS extras; `laspy.open` header reports version 1.4 / point
format 0; `read_canonical_las` returns xyz + all extra dims intact
(tree_id, height_above_ground values match what was written).
Result: pass.

**Defaults check** (W11): assert `run_pipeline.run_pipeline.__kwdefaults__`
radius 7.315 / stem 0.05 / canopy 0.05 / gen DTM 0.25; `stem_map`
DEFAULT_CELL_SIZE 0.05 / DEFAULT_THRESHOLD 0.25; `generate_dtm`
--cell-size default 0.25.
Result: pass.

**End-to-end pipeline smoke** with the new defaults (6-tree truncation of
the FIA output; 111 trees is ~120s+ of generation even at 64 points/m):

```
head -7 /tmp/fia_test/plot_11839.csv > /tmp/plots_small/plot_11839.csv
conda run -n ssf python run_pipeline.py /tmp/plots_small --radius 20 --seed 42 --points-per-metre 64 --limit 1
```

Result: 1 ok / 0 failed in 15 s; stem map "597x653 px @ 5 cm, band
1.12-1.62 m" (= DBH_HT +/- STEM_THRESHOLD with the W11 defaults), CHM
@ 5 cm, plot + 6 tree info JSONs, normalize fell back to 'trimmed'
(no terrain in this run — expected without --generate-dtm).

### Verification gaps / notes

- The report JSON and CLI smoke were verified by inspection, not asserts.
- No automated regression suite yet — the snippets above are the checks;
  consider promoting them into `tests/` as work proceeds.
- Pre-W2 seeds do not regenerate identical plots (see reproducibility entry).

### 2026-11-05 — W1: minimal JSON tree list contract (commit f26e2d6)

**Unit checks** (`conda run -n ssf python /tmp/test_w1a.py`, all pass):

1. FIA CLI default format emits the canonical `plot_11839.json` (111 trees)
   + report; `tree_list.load_plot` verifies plot_uid/shape/radius, null
   plot_center, required attrs, DBH/HT plausibility, and that every
   TREE_CODE resolves in the functional-group database.
2. `load_rows` + `rows_from_tree_list` produce generator rows numerically
   identical to the legacy SST table (`--format sst-csv`, same seed) on
   HT/BASE_DIA bounds, MAX_CAN_W, SPLIT_*, SPR, TOP_DIST, X, Y
   (atol 1e-4). Found+fixed during this check: round(4) on the tree
   list's DBH/HT shifted allometry by ~1 mm; DBH/HT now carry full
   precision.
3. `resolve_func_group` resolution order: species code (202 -> DF) ->
   functional group -> DEFAULT (unknown code).
4. `validate_plot` rejects non-physical lists (negative DBH).
5. `tree_list.trees_frame` parity: same TREE_UID/X/Y/geo columns from
   both contracts; identical X/Y values across them.

**End-to-end**:

- JSON contract: 6-tree subset of plot_11839 through the full
  `run_pipeline.py` -> 1 ok / 0 failed, plot cloud 8430 points (stem map
  @5 cm band 1.12-1.62 m, CHM @5 cm).
- Legacy CSV regression: the same 6 trees via the SST CSV -> 8430 points,
  identical — the two contracts are pipeline-equivalent for the same
  trees and seeds.
- Replicate on JSON: `plot_11839_r01/r02.json` + env sidecars written;
  `tree_list.find_plots` excludes `_report.json`/`_env.json` sidecars.

**Check**: `check.py` passes (26 files, tree_list.py included). Design doc
statuses updated (W1 closed; 3.1/8.1 implemented-except-W13; 6.1 output
line; Section 2 partial note).

### 2026-11-05 — W1 review pass (commit f4087c4; post-commit adversarial check)

Found and fixed two real bugs the first round's tests missed (both in
plot-table discovery, both caught by writing tests for what WASN'T
covered):

1. `tree_list.find_plots` had no `plot_` prefix check — any `*.csv` /
   `*.json` in a plot directory qualified as a table. `replicates_manifest.json`
   (written into replicate output dirs, which run_pipeline is documented
   to consume) and stray `manifest.csv` files would be picked up as plots.
   Fixed: only `plot_`-prefixed files qualify.
2. `is_tree_list_path` excluded `_env.json` (underscore) but replicate's
   sidecars are named `plot_<UID>_rNN.env.json` (DOT env) — they matched
   the tree-list pattern and would enter the pipeline as plots. Found by
   asserting the discovered-table list instead of eyeballing it (the
   earlier e2e replicate run listed the dir by hand and missed this).
   Fixed: exclude `.env.json`.

Review suite (`conda run -n ssf python /tmp/rev_check2.py`, passes):
(1) a dir containing manifest.csv, replicates_manifest.json, notes.txt and
one real tree list discovers exactly the tree list; (2) a replicate dir
discovers only `plot_*_rNN.json`; (3) both contracts load via load_rows
with X-parity on shared trees.

Also fixed this round: stale spec references (7.1's W1 -> W3; 7.2 step 7
and 7.3 status paragraphs now describe the implemented generation-time
path); README_fia_helper / README_generation / README_replicate document
the JSON default and the `--format sst-csv` legacy path.

Known gap (noted, unchanged): tree Z from the tree list is not yet
consumed by generation (Z is 0 for all current sources; will matter for
terrain work, track with W5).

### 2026-11-05 — W3: generator module interface (commit b0da89e)

**Unit checks** (`conda run -n ssf python /tmp/test_w3.py`, all pass):

1. Registry: `get_generator("sst")` instantiates; unknown names raise
   KeyError; `register` used by the SST module.
2. `resolve_tree_code` single-tree form: RW -> RW, SPCD 202 -> DF,
   unknown -> DEFAULT.
3. Canonical `_expand_entry` vs `legacy_row_params`: model_params
   equivalent (atol 5e-4; the residual is the legacy CSV's 4dp rounding).
4. build determinism: identical rng -> bit-identical TreeModel.points;
   base at origin (base_xy [0,0], crown within its reach of the base).
5. `place()` vs generating at (x, y): replaying the module's rng bridge
   and calling gen_simtree with tx,ty directly matches the
   place()-shifted model to 3.6e-15 — placement via shift is exactly
   equivalent to generating positioned.

   (An earlier draft of check 5 compared against the WRONG global seed
   and reported an 8.96 m "mismatch" — the module bridges the rng
   through one draw, so the comparison must replay that bridge. Test
   bug, not code bug; noted here so the 8.96 number is not misread
   later.)

**Found+fixed during first run**: SSTGenerator.load_param_dict/_spcd
passed `None` paths into load_func_group_db/load_spcd_map (crash);
`entry_for` called np.isnan on the DESCRIPTION text column. Both caught
by the unit suite before any e2e run.

**End-to-end** (`/tmp/test_w3e2e.py` + `/tmp/w3_repro.sh`):

- JSON contract, 6 trees, full pipeline: 1 ok / 0 failed, 8430-point
  trimmed cloud; legacy CSV contract, same 6 trees: 8430 points.
- Cross-contract parity: identical point count and ALL extra dims;
  sorted-coordinate max diff 1.996e-03 m (the legacy table rounds
  HT/BASE_DIA ranges to 4dp, which propagates through the stem spline).
- Reproducibility: same table + seed, two runs -> byte-identical
  trimmed cloud (md5 match).

**Behavior notes**:
- Seeding bridge: treesim draws from numpy's GLOBAL generator; the
  module now bridges the core rng with one rng.integers draw. Point
  counts are unchanged for these plots, but clouds generated pre-W3
  will NOT regenerate byte-identically from the same seed (the global
  stream changed). Within-version reproducibility is exact (md5 check).
- The rng bridges per attempt: the core rehashes the seed per retry
  (base + attempt*7919, spec MAX_ATTEMPTS), so retry paths stay
  deterministic.
- Test-suite supersession: the W1 unit suite (`test_w1a.py`) imported
  `load_rows` / `rows_from_tree_list`, which the W3 refactor replaced
  with `tree_requests` / `SSTGenerator._expand_entry`+`legacy_row_params`.
  Checks 3-5 of the W1 suite are covered by this round's checks;
  contract/discovery/parity checks (rev_check2.py) still apply
  unchanged.

### 2026-11-05 — W3 review pass (adversarial check on commit b0da89e)

Found one real regression the W3 verification missed, plus one important
coverage gap:

1. **Masked test failure** (real regression in the test suite, not the
   code): the W3 entry claimed "rev_check2.py still applies unchanged",
   but its check 3 imported the removed `generate_trees.load_rows` —
   it has been failing with AttributeError since the refactor. The
   failure was masked because the rerun was piped through `tail -2`,
   whose exit code is 0. Fixed: check 3 now uses `tree_requests` /
   `get_generator` (X/Y/UID parity across contracts); suite green with
   `exit: 0` verified explicitly. Standing check 4 (exit-code
   discipline) was added so this cannot recur silently.
2. **Collision trimming had never been exercised since W1**: every e2e
   run used the sparse 6-tree subset where raw == trimmed volume (zero
   collisions), so the W3 refactor's riskiest path — trimming against
   the module-placed skeletons — was untested. Now verified on a dense
   subset: 40 trees (r=20 plot), `-g sst`, 32 pts/m through generate +
   assemble -> 40/40 generated, **376 collisions, 202 points trimmed
   (2.6%)**, per-tree element trims recorded (e.g. tree 040: 5 elements,
   2 points removed). Trimming works correctly on the placed module
   skeletons.

Also verified this round: no stale symbol references remain
(param_editor's `stem_frac_from_row` is its own local copy, not an
import); the manifest columns dropped by the refactor
(fork/n_branches/n_stem_points) are consumed by nothing (restructure
only moves manifest.csv); run_pipeline `-g sst` pass-through works —
exercised through a complete 6-tree pipeline run after an earlier CLI
probe (which mutated the shared /tmp e2e fixture; the fixture was
regenerated and parity re-verified exit-0). Test-hygiene lesson: run
ad-hoc CLI probes against scratch directories, never against fixture
dirs other checks depend on.

**Verification gaps updated**: dense-plot timing — 40 trees at 32 pts/m
generate+assemble completed well within a 300 s budget; the earlier
"6-tree e2e only" gap is now "full-pipeline e2e still sparse-subset".
Trimming coverage exists; stem/canopy/normalize on a dense plot remain
untested (derived stages, unchanged by W3).

### 2026-11-05 — final pre-push sweep (W1+W3 stack)

The last adversarial pass before pushing, on the whole unpushed stack
(b0da89e W3, 27711d2 hash pin, plus this round's fixes):

**New coverage — geo passthrough** (`conda run -n ssf python
/tmp/test_geo.py`, passes with exit 0): every earlier round used
CA_TREE_SMALL.csv, which has NO PLOT join, so the LAT/LON/ELEV
passthrough had never been exercised. A synthetic 60-tree FIA table
with PLOT-join geo columns goes through `fia_to_tree_list` ->
`plot_center` = {lat 44.1234, lon -121.5678, elev 950.0}; inventory-only
tree arrays (no SST leakage); JSON round trip preserves the center;
`trees_frame` exposes plot_lat/lon/elev (what assemble writes into the
scene JSON); and the georeferenced list generates through the module
interface. Found no bugs — the path was correct but unverified until
now.

**Found+fixed: README.md staleness** — four places still described the
pre-W1/pre-W3 world: the quickstart held "plot_*.csv tables", the
pipeline diagram showed `fia_helper -> plots/plot_<UID>.csv` and
"one labelled cloud per tree (treesim.py)", module 1's description
said CSV output, and the output-layout listing had no tree-list entry.
All updated to the JSON default (with legacy mentions where relevant)
and the generator-module call path.

**Full-sweep results** (each with an explicit exit-code check,
per standing check 4):
check.py PASS, W3 unit suite PASS, W1 review suite PASS, cross-contract
e2e parity PASS, geo passthrough PASS. Stale-symbol grep: only the
intentional "legacy SST table" mention remains.

**Still-open verification debts (unchanged, tracked)**: tree Z
unconsumed (W5); scene_editor GUI untested; full-pipeline runs on dense
plots (stem/canopy/normalize after trimming) untested; /tmp test
scripts are ephemeral — promoting the standing suites into a repo
`tests/` directory is the natural next hygiene step.

### 2026-11-05 — coverage-debt round 2 (dense derived stages, terrain, replicates)

Round goal: close the three open coverage debts from the prior sweeps.
**No notable bugs found this round** — all three previously untested
workflows pass.

1. **Dense-plot derived stages** (the /tmp/w3_dense 40-tree plot, 376
   collisions / 202 points trimmed, left at generate+assemble): ran
   stem_map, canopy_map, restructure, build_info, normalize -- each
   with an explicit exit check. ALL PASS; normalized plot cloud 7427
   points (7629 raw - 202 trimmed). Post-trim per-tree artifacts
   (stem bboxes, crown hulls, info JSONs for trimmed trees) compute
   without errors.
2. **Terrain + understory on the JSON contract** (never run since W1's
   refactor; all prior e2e skipped the terrain stage):
   `run_pipeline.py <6-tree JSON list> --generate-dtm --understory` ->
   1 ok in 28 s; 29,943-point terrain cloud, 21,513 ground points;
   asserted the spec Section 3.3 invariant — ground-point HAG exactly
   0.0, no NaN anywhere; three per-shrub dirs written. Understory's
   DTM reproduction from the scene record works through the JSON path.
3. **Replicate + env-sidecar workflow E2E** (the documented consumer of
   replicate directories): replicate x2 of a 6-tree JSON list ->
   /tmp/w3_reps/{plot_11839_r01,r02}.json + .env.json sidecars ->
   `run_pipeline.py /tmp/w3_reps --generate-dtm --understory` -> 2 ok /
   0 failed in 48 s. Matched each plot's scene-JSON terrain record
   against its sidecar draw: r01 bumpiness 0.5619 / relief 12.3304 and
   r02 0.79 / 12.4153 — exact per-plot env application confirmed, and
   the two plots' generated DTMs differ accordingly.

Coverage debts remaining after this round: scene_editor GUI (no
display), tree Z consumption (W5), full 111-tree scale runs
(performance, not correctness), and the /tmp suite ephemerality
(promote to tests/). None of these blocked any workflow tested.

### 2026-11-05 — W4: flat-plane topography (commit c456162)

**Design** (spec 3.3): the None topography is a uniform flat plane at
z = 0 implemented as `DTM.flat(...)` — a zero grid with the full DTM
interface (sample -> 0 everywhere), so draping/HAG/classification/ground
sampling run unchanged; terrain is no longer skippable (the stage always
runs; absence of --dtm/--generate-dtm means flat); conflicts resolve
None > filepath > random with a warning; unconfigured density in flat
mode warns and defaults to 1.0 per cell.

**Unit/e2e checks** (all exit-0):

1. Flat full pipeline + understory on the 6-tree JSON list: PASS;
   scene record reads "flat plane (topography = None, z = 0)" with
   dtm_generated false; orchested run writes 6 per-tree terrain files.
2. Flat invariants (/tmp/w4_inv.py, on /tmp/w4_flat2): understory cloud
   102,806 pts with 21,513 ground at |z| == 0.0 EXACTLY, ground HAG ==
   0 exactly, no NaN; 10,704 shrub + 62,159 grass points draped on the
   plane (grass z within [0, 0.35]); vegetation |z - HAG| <= 5.05e-4 —
   i.e. HAG == z up to the LAZ 1 mm coordinate quantization, NOT a bug
   (my first assert used 1e-5, below the storage scale; loosened and
   documented).
3. Standalone flat terrain with no --density: exit 0 AND the spec warning
   present ("no ground sampling density provided ... default 1.0");
   flat terrain cloud 24,665 points (0 tree files — expected, that run
   had no canopy_map so no per-tree bboxes; pre-existing behavior, the
   orchestrated run writes all 6).
4. Conflict resolution: pipeline run with both --dtm <generated.asc>
   and --generate-dtm: exit 0, orchestrator warning printed, provided
   DTM used (filepath wins). Standalone terrain.py resolves and warns
   in build_terrain the same way.
5. Generated-DTM regression through the new wiring (--generate-dtm
   --understory): identical 29,943-point terrain cloud, 6 tree files,
   scene record "generated ... cell 0.25", ground HAG == 0 exactly,
   z range 0.73..53.71 m.

**Notable non-bug**: terrain.py standalone --gen-cell-size default
aligned 1.0 -> 0.25 (DTM_CELL_SIZE, W11 leftover the standalone CLI had
kept). My own edit duplicated the --gen-cell-size argparse line during
the CLI rewrite (argparse would use the last, silently restoring the
old default) — caught by rereading the diff before running anything;
guarded by the density-warning e2e which exercises the flat grid's cell
semantics. (Process lesson, related to the log-edit mistakes of earlier
rounds: insert new entries by targeting a UNIQUE anchor string, never
a heading used as a section boundary for content to be preserved.)

### Verification gaps / notes (W1 round)

- scene_editor's JSON acceptance is path-level (`find_plots`); its GUI was
  not exercised (no display).
- FIA geo passthrough (per-tree LON/LAT/ELEV from a PLOT join) is not yet
  derived from plot center + placement (item W5); this sample has no PLOT
  join (plot_center all-null).
