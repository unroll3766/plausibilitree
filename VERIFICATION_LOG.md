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
- Repo regression suites: `conda run -n ssf python tests/run_all.py` — the
  promoted W7/whole-repo invariant suites (placement & frame; cross-artifact
  invariants). Run after every code change; see `tests/README.md`.

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
   (2.6%)** [numbers computed against the place() frame bug later found
   by the W7 review; collision counts from the W3 era onward are
   phantom-trunk-tainted until that fix — see the W7 review entry],
   per-tree element trims recorded (e.g. tree 040: 5 elements,
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

### 2026-11-05 — W4 review pass (adversarial check on commit c456162+40cf20d)

**No notable bugs found this round.** Checked the surfaces the W4 round
did not cover:

1. **scene_editor compatibility**: its terrain invocation always passes
   --dtm or --generate-dtm (never relied on the old skip), so
   always-on terrain in the orchestrator does not affect it.
2. **Flat run WITHOUT --understory** (untested in the main round): full
   pipeline exit 0; normalize consumed the TERRAIN cloud directly
   (semantic codes [1,2,3,6] including ground; plot cloud 29,943 =
   terrain count, no fallback warning). Behavior note now explicit:
   a plain run (no flags) previously normalized the trimmed cloud
   (8,430 pts); with always-on terrain every default run's products
   include flat-plane ground points — spec-conformant (terrain always
   runs) but a visible default-output change.
3. **Mixed replicate dirs**: a sidecar-free plot in a replicate run now
   gets flat terrain instead of being skipped — intended Section 3.3
   behavior; the sidecar plots still use their drawn DTMs.
4. Absent-DTM handling: a --dtm path that doesn't exist still fails at
   read (same as before, now after the conflict warning — user asked
   for the file).

**Found+fixed: README/doc staleness**: the main README quickstart still
said terrain needs "--generate-dtm --understory" (changed to "terrain
always runs, flat by default"); the run_pipeline docstring still showed
the pre-W4 six-stage list without terrain (updated with the always-on
terrain stage and the None topography note) and the --dtm help still
said "enables the terrain step".

Standing checks after fixes: check.py PASS, doc cross-references clean.

### 2026-11-05 — W13: parameter dictionary JSON + editor import/export (commit 1939061)

**Design** (spec 8.6): the parameter dictionary is a JSON object keyed by
tree code (`data/functional_groups.json`), one entry per code with the
Section 7.3 columns. Converted from the shipped CSV with a one-off script
(preserving every field; NaN -> null); the CSV reader stays as the legacy
input path (`load_func_group_db` auto-detects by extension, returning the
same FUNC_GROUP-indexed frame either way — callers unchanged).
param_editor.py: `load_table`/new `save_table` handle JSON (default) and
legacy CSV; open/save pickers and CLI strings updated.
`data/functional_groups.csv` removed (CSV is a legacy *input* form per
spec, not a shipped artifact).

**Unit checks** (`conda run -n ssf python /tmp/w13_unit.py`, exit 0):

1. Default dictionary resolves to the JSON; 9 groups; DEFAULT present;
   index named FUNC_GROUP.
2. **Value parity**: every group x every column identical between the JSON
   and the pre-conversion CSV (recovered via
   `git show HEAD:data/functional_groups.csv`), NaN handling included.
3. `entry_for` through the JSON returns DESCRIPTION (string path) and
   numerics; unknown code resolves to DEFAULT's parameters.
4. Legacy CSV dictionary still loads through the same loader.
5. Editor round trip: load JSON -> save JSON -> load; frame-equal and
   JSON-equal (9 groups); legacy CSV save path still writes/loads CSV.

**E2E** (`/tmp/w13_e2e.sh`, exit 0): a fresh full 6-tree pipeline run
reading the JSON dictionary produces a md5-identical plot terrain cloud
(68c29645...) to the fixture run generated from the CSV dictionary — the
conversion is value-exact end to end.

**Found+fixed during this round (before any test run)**: my first
JSON-loader draft reset the FUNC_GROUP index into a column, breaking
`_fg_param`'s `lut.loc["DEFAULT"]` — caught while rereading the diff.
No runtime bugs found after that.

### 2026-11-05 — W13 review pass (adversarial check on commit 1939061+e72fafa)

**No notable bugs found this round.** Checked the surfaces the W13 round
did not cover (all exit-0, /tmp/w13_review.py + /tmp/w13_review2.sh):

1. **The FIA helper core path** (crown-width allometry, taper, group
   resolution — all `_fg_param` consumers) never ran through the JSON
   dictionary in the W13 round (only generation did). FIA helper via JSON
   dict: 111 trees placed, canonical JSON output valid.
2. **User-added groups**: a custom dictionary with a grafted "ZZ" group
   (CW_A0 modified) loads via the module, resolves, and BUILDS — the
   hot-swap path (load_param_dict(path)) works through JSON.
3. Unknown-code DEFAULT fallback through the JSON dictionary.
4. **Stale-reference grep**: no `functional_groups.csv` mentions remain
   in .py/.md outside legacy-captioned contexts and this log.
5. **Legacy FIA output + legacy pipeline contract through the JSON
   dictionary**: `fia_helper --format sst-csv` and a full legacy-CSV
   pipeline run both pass. The legacy-contract trimmed/terrain cloud
   differs from the JSON-contract one only by the documented cross-
   contract variance (frozen ranges rounded to 4dp in the legacy CSV;
   the W3 review measured 1.996e-03 m) — expected, not a regression.

### 2026-11-05 — W7: random mirror reflection in placement draws (commit b3023b9)

**Design** (spec 3.1): placement draws = uniform rotation on [0, 360)
about the tree's own base PLUS a 50/50 mirror reflection across the
vertical plane through the base. Implemented in `assemble_plot` (where
rotation already lived): a new `mirror_x(xy, cx)` in-place mirror,
`TreeSkeleton.reflect()` for element geometry (distances/radii/topology
are reflection-invariant, so collision trimming is unchanged beyond the
new geometry), the per-tree Bernoulli draw in the placement loop between
the rotation draw and the collision-tiebreak stream, and a `reflections`
dict in the scene JSON next to `rotations_deg` plus a `reflected` flag in
each trimmed skeleton. The serialized skeleton geometry (polyline nodes,
branch base/tip) mirrors exactly like the point cloud, so the skeleton
JSON keeps matching the LAZ beside it.

**Unit checks** (`conda run -n ssf python /tmp/w7_unit.py`, exit 0),
fixture = 6-tree subset of plot_11839 (r=20, seed 42, 64 pts/m):

1. Determinism: same seed twice -> bit-identical trimmed cloud +
   identical reflections/rotations records.
2. Scene JSON: reflections covers all 6 trees, values bool; both
   states occur within the plot.
3. **Cloud-vs-records parity**: for every tree, applying the recorded
   rotation+mirror manually to the per-tree raw cloud reproduces the
   merged `_raw.laz` points for that tree -- max |d| 5.0-5.1e-4 on every
   tree, exactly half the 1 mm LAS scale (the merged cloud re-quantizes
   at ITS offset); my first tolerances (1e-6 xy, 0 z) failed on this
   re-quantization -- fixed by documenting it, not loosening blindly.
   (The z diff, up to 4.5e-4, is the same effect: z is not transformed
   but is re-quantized when merged.)
4. **Skeleton-vs-records parity**: replaying the recorded 4dp
   rotation + mirror on each raw skeleton reproduces the trimmed
   skeleton's polylines and branch base/tip to 2e-5 m. The residual
   (~1.2e-5 m at 18 m reach) is the 4dp angle record itself (5e-5 deg =
   8.8e-7 rad); a real transform error would be metres. (First run had
   this check pass on 5/6 trees and fail on 003 by 1.18e-5 -- np.allclose's
   DEFAULT rtol=1e-5 was silently covering the others; now rtol=0 with an
   explicit atol. Test bug, not code bug.)
5. Draw distribution: seeds 1000..1007 -> both states occur, aggregate
   rate ~50% (37/48 expected range).

**Dense-plot collision coverage** (the riskiest untested path — trimming
against mirrored skeletons; 40-tree subset, 32 pts/m, seed 42):
originally recorded as 319 collisions / 250 points trimmed (21/40
reflected). [SUPERSEDED — those runs, and the "sparse 6-tree plot now
has 6 collisions / 125 points trimmed" observation below, were measured
against the frame bug this round's review uncovered in
`sst_generator.place`; with the bug fixed the 6-tree plot is back to
0 collisions / 0 trimmed and the dense plot trims 337 collisions /
314 points (4.1%). See the W7 review entry.]

**E2E** (fresh scratch dirs, never the unit fixture): full pipeline on
the 6-tree list -- 1 ok / 0 failed in 15 s; scene JSON carries the
reflections; two independent full runs -> **byte-identical trimmed
cloud** (recorded md5 7c6dffd85d145b0dab36db4a6051f463 on the
frame-bugged code; post-fix fixture md5
891ad001acbd45b38d5ea38c6ccf398e, still byte-identical across runs) and
identical reflections. Derived stages (stem map from mirrored skeleton
discs, CHM, restructure, flat terrain, build_info, normalize) all ran
green through the pipeline.

**RNG-stream caution (expected behavior change)**: the reflection
Bernoulli draw precedes every subsequent draw, so same-seed outputs
differ vs pre-W7 code (tree 001's rotation is unchanged — it is still
the first draw — but no other draw reproduces). Pre-W7 seeds do NOT
regenerate identical clouds.

**Process notes (own bugs, caught before/at first run)**: my first edit
line made `rotations` a list in the tuple-unpack (`[], [], [], {}`),
crashing at the first placement — caught by the determinism check before
any e2e; another edit dropped the `if not trees:` guard in
assemble_plot's placement loop while adding the reflections record,
restored immediately after re-reading the hunk. Both reinforce the
diff-reread-before-run discipline.

### 2026-11-05 — W7 review pass (adversarial check on commit b3023b9+91d959d; fix in 82aabe6)

**Found one REAL, significant bug — pre-existing since W3, not a W7
regression.** The W7 code is innocent: the review's same-seed
collision-change observation (sparse plot 0 -> 6 collisions) was the
thread that unraveled it, and the frame consistency invariant now
guards it.

1. **`SSTGenerator.place()` silently dropped the trunk-polyline shift**
   (`sst_generator.py`). `gen_simtree` emits polylines as plain lists
   (`.tolist()`), and place() did `arr = np.asarray(poly); arr[:, 0] += x`
   — np.asarray COPIES a list, so the shift landed on a discarded array.
   Points, `base_xy` and every branch base/tip used in-place list
   mutation and shifted correctly. Consequence: since commit b0da89e
   (W3), every module-generated raw skeleton declared
   `coordinate_frame: "plot metres"` (spec Section 7.4) but carried
   trunk polylines in the LOCAL base-origin frame while everything else
   was global. `assemble_plot.TreeSkeleton` built MIXED-frame collision
   geometry: trunk elements (local) were then rotated about the global
   base and drifted ~|base_xy| (4-14 m on these plots) from the true
   stems — phantom trunks. Impact: all collision counts and
   trunk-element trims since W3 were computed against phantom geometry
   (branch-vs-branch collisions were always real; trunk-involving ones
   were not); the sparse-plot "6 collisions / 125 points trimmed" that
   the W7 round attributed to mirror-induced crown overlap was ENTIRELY
   phantom. The W3-review equivalence check (place vs generating
   positioned, 3.6e-15) compared model.points only — the points shift
   was always correct, which is why dozens of cloud-level parity checks
   passed while the skeleton frame was wrong. **Fix**: write the shifted
   array back (`polys[i] = arr.tolist()`), with a comment explaining
   the list/ndarray trap.
2. **How it was caught**: a review check the W7 round didn't have — the
   per-tree trimmed skeleton's centerline vs the per-tree trimmed
   cloud's stem points. Reflected tree 001's stem points were 9.57 m
   off-axis from "its own" serialized centerline; the raw skeleton then
   showed polyline0 node0 = (0,0,0) with base_xy = (-3.10, 8.45) and
   branch bases already global. (Lesson, again: check the INVARIANT,
   not just stage exit codes — every stage "ran green" on phantom
   geometry because both sides were consistently wrong together.)

**Re-verification after the fix** (fixture regenerated; point counts,
rotations and reflection draws are identical pre/post fix — the bug
never touched rng or the point cloud):

- Frame consistency: every raw polyline node0 == base_xy; all trunk
  nodes within crown reach of the base (new permanent checks 0.x of
  /tmp/w7_unit.py, exit 0).
- W7 unit suite: ALL PASS (determinism, records replay, chirality,
  draw distribution).
- Sparse 6-tree plot: **0 collisions / 0 trimmed** — the phantom
  collisions gone; pre-W7 behavior for THIS fixture restored, so the
  reflection does not create real collisions at sparse spacing.
- Dense 40-tree @32 pts/m: **337 collisions / 314 points trimmed
  (4.1%)**, 21/40 mirrored — real trimming against real (mirrored)
  skeleton geometry, including trunk-involving collisions.
- Full pipeline e2e x2 (fresh scratch dirs): 1 ok / 0 failed; trimmed
  cloud **byte-identical across runs, md5 891ad001acbd45b...** (differs
  from the pre-fix fixture, as expected with real skeleton geometry);
  all derived stages green.
- Stem-vs-skeleton coax check (the one that caught the bug): max
  off-axis 0.18-0.58 m = trunk radius, on reflected and unreflected
  trees alike. The stem MASK stage now genuinely sees its tree.
- W7-specific claims re-verified post-fix: **RNG-stream comparison**
  vs the febb5ae worktree at the same seed: tree 001's rotation is
  identical (it is the first draw), all others differ, pre-W7 scene has
  no `reflections` key — the W7 entry's caution is exactly right.
  **Chirality**: signed-area over ordered branch bases flips sign
  exactly on reflected trees (001, 005: |S_out + S_raw| ~ 1e-15) and is
  preserved on the others — the mirror is a true handedness flip, not a
  substitute rotation.
- Invariants on a mirrored plot with generated DTM + understory:
  ground HAG == 0 exactly, no NaN (spec 3.3); scene_editor has no
  skeleton/rotation consumers (cloud display only).

**Also superseded**: the W7 entry's dense numbers (319/250) and sparse
"now trims" claim are annotated in place; W3-review collision counts
(376/202) likewise annotated as phantom-tainted. Historical request to
future rounds: any collision/trim counts recorded before this fix are
not reproducible on fixed code.

### 2026-11-05 — whole-repo adversarial review (fixes in 16e3998)

Triggered by finding a pre-existing bug hidden across many rounds; scope =
every module. Method: (a) static sweep for the same bug class
(np.asarray-then-mutate copies, swallowed exceptions), (b) frame/contract
audit of geometry-bearing writers/readers, (c) a NEW cross-artifact
invariant suite recomputing each artifact from a DIFFERENT artifact
(/tmp/wrev_invariants.py, 80 checks: skeleton-vs-cloud coax, info-JSON
bbox/dbh/frustum-volume recomputation, CHM, branch-id encoding, DTM
invariants, understory domains, generation records, scene records vs
table), (d) a scene-JSON link audit, (e) re-runs of every workflow e2e
never yet validated on post-frame-fix code.

**Three REAL bugs found and fixed** (all in the draping/scene-link
family — none exercised by the pre-W7 fixtures, which were flat or
zero-noise):

1. **`terrain.py` tree draping violated spec 3.3's anchor rule.** Spec:
   "Each tree is dropped vertically so its lowest point sits on the DTM
   at its (X, Y)." Code sampled the DTM at the LOWEST POINT'S OWN xy —
   on rugged terrain a downhill branch tip anchored the whole tree to
   a gully elevation, burying trunk bases below local ground (measured:
   trunk points with HAG to -0.86 m) and floating them on slopes. All
   prior terrain fixtures were flat (both samples identical), so the
   deviation was invisible. Fix: anchor shift at the tree's (X, Y) read
   from the plot table (tree_id == int(TREE_UID)), legacy fallback to
   the old behavior if the table is unreadable.
2. **`vegetation.py` shrubs were center-anchored, not draped.** The
   whole shrub mesh was placed at its CENTER's DTM elevation, burying
   upslope flanks (measured: 1019/10,800 shrub points with HAG to
   -1.22 m) — while the grass layer already draped per point
   correctly. Fix: per-point draping (sample the DTM at every point's
   own xy), matching grass; `shrub_geom.base_z_m` records the center
   ground elevation.
3. **`restructure.py` wrote dangling scene-JSON links.** (a)
   `plot_files.volume` pointed at the plot-level volume label JSON,
   which this stage computes-and-prints but RETIRED from writing
   ("superseded by build_info's _info.json") — a dead link in every
   scene JSON, violating the spec's every-artifact-resolves role; the
   link is now retired with the artifact. (Correction to this entry's
   first draft, which claimed "spec never listed plot-level volumes":
   the spec DOES list them — Section 8.4's Plot info JSON, which
   build_info duly writes under
   `vegetation_characteristics.trees.*_volume_m3` {raw, trimmed} — the
   retired artifact is the spec-unlisted deep-learning volume LABEL,
   a different thing entirely. Found while answering a user question
   about this entry.)
   (b) `stem_map.source_cloud`/`canopy_map.source_cloud` recorded the
   FLAT-layout cloud paths, which restructure itself then MOVED to
   `las/` — stale since the canonical-layout refactor; now rewritten
   to the canonical paths (same files). Module docstring updated.

**Verified green after the fixes** (each with explicit exit codes):

- Whole-repo invariant suite on a fresh full fixture (6-tree, r=20,
  seed 42, generated DTM + understory): 80/80 PASS. New spec-exact
  anchors asserted: every tree's lowest point == DTM at its (X, Y)
  (1.5 mm); shrub+grass HAG >= -1 mm; per-tree trunk z-extent ==
  height_m; trunk frustum volume == info JSON to 1e-3; scene records
  cover the table.
- Scene-link audit: every relative path in the scene JSON resolves
  (62/62 after the fix; 3 dangled before).
- **Flat-mode regression**: both draping fixes are provably no-ops on
  flat terrain — flat full pipeline (with understory) byte-identical
  (md5 terrain 26f7a95c..., understory cf62bf8e...) on both sides of
  the fix (verified via a 9361d61 worktree run of the same input).
- Rugged-regression behavior change: quantified on the fixture
  terrain — all 6 trees re-anchored (largest shift 1.83 m on tree
  002); shrub HAG min went from -0.866 m to 0.000 m.
- Replicate workflow x2 with drawn rugged envs (bumpiness 0.67, relief
  5.7, ground_noise 0.026): shrub/grass HAG >= 0, tree anchors exact,
  ground HAG == drawn noise +/- 1 mm (the noise is by design per the
  env sidecar; the "HAG == 0 exactly" invariant holds only for
  zero-noise runs — documented here so it is not misread as a bug).
- Legacy SST-CSV contract pipeline (6 trees, DTM + understory):
  1 ok / 0 failed.
- Dense 40-tree derived stages on the post-frame-fix assembly
  (replacing the W3-era tainted validation): stem/canopy maps
  (685x758 px — the phantom skeletons no longer inflate extents),
  restructure, build_info, normalize all exit 0; normalized plot
  cloud 7,315 pts (7,629 - 314 real trims).
- CLI smokes of the not-yet-run tools: colorize_plot.py (103,297 pts
  colored from the understory cloud), generate_example_tree.py (PLY +
  LAS written in a scratch dir), scene_editor remains GUI-untestable
  (no display) — static check only.
- Standing checks: check.py PASS (it caught a dangling-name of my own
  mid-round: my first vegetation.py edit left `base_z` referenced in
  shrub_geom — fixed before any run), doc cross-refs clean, stale
  identifier grep clean.

**Test-side findings (documented, not bugs)**: foliage discs extend up
to ~1.1 m above the trunk tip on DF forks (model behavior; the trunk
z-extent == height_m check is the exact one); CHM datum is
height-above-ground, not draped z; ground-on-DTM tolerance must be
~2.5 mm (1 mm xy-quantization x slope + 0.5 mm z scale); per-tree
`*_dtm.asc` files are CROPS (validating against them is only valid
within each tree's bbox); np.asarray-copy traps found nowhere else in
the repo (remaining `except` sites are legitimate best-effort
fallbacks).

**Late find while answering a user question about the volume link
(correction folded into item 3a above): Section 8.4's plot info field
list overclaims its [implemented] status.** Present in
`<PLOT_UID>_info.json` and spec-conformant: plot_uid, plot_shape,
n_trees, n_shrubs, bounding box, the {trunk, branch, foliage, canopy}
x {raw, trimmed} volume set, percent cover, and meta (plot center,
coordinate system, units). MISSING: "heights (voxel counts /
percentiles)" and "neighborhood statistics" — neither exists at plot
level (per-tree info carries height_m and neighbor counts, but no
plot aggregate). Decision needed: implement plot-level height/neighbor
aggregates in build_info, or amend Section 8.4 to say what build_info
actually writes. Until then this is a tracked doc-vs-code gap, not a
validated claim.

**Process lesson reinforced**: every bug found lived in a code path
exercised ONLY by configurations the fixtures never used (rugged DTM,
negative-noise, post-restructure link consumers). Flat + zero-noise
fixtures masked spec deviations for months of rounds; the W7-round
skeleton-frame bug was the same shape (no fixture ever checked the
invariant). Cross-artifact invariant checks (recompute from a
different artifact than the writer) are the pattern that finds them.

### 2026-11-05 — tests/ promotion: standing suites moved out of /tmp (commit ab0077c)

Closed the long-tracked /tmp-suite ephemerality debt: the two live
suites behind the W7 and whole-repo rounds are now repo files,
self-contained (they build their own fixtures via the programmatic
entry points `fia_to_tree_list`/`generate_plot`/`run_pipeline` — no CLI
subprocesses), deterministic (seed 42, radius 20, 6 trees at 64 pts/m,
temp dirs via mkdtemp, `--keep` for debugging), and exit-0 gated.

**Layout**: `tests/fixture.py` (shared fixture builder),
`test_placement_reflection.py` (frame-consistency guards for the
place() polyline bug, same-seed determinism, records replay onto cloud
and skeleton, ~50% mirror draw distribution, chirality via signed
area), `test_cross_artifact_invariants.py` (scene records vs table;
skeleton-vs-cloud coax; info-JSON bbox/DBH/frustum-volume recompute;
branch_id/branch_code domains; spec 3.3 draping anchors + shrub/grass
per-point drape on a RUGGED generated DTM; understory domains;
generation records; scene-JSON link audit), `run_all.py` (aggregates
exit codes), `README.md` (what each suite guards + conventions).

**Verification of the promotion itself** (each with explicit exit
codes, per standing check 4):

- `conda run -n ssf python tests/test_placement_reflection.py` ->
  exit 0, all checks PASS (frame guards, determinism bit-identical,
  cloud replay max |d| 5.0-5.1e-4 on every tree, skeleton replay
  within the 4dp-angle tolerance, draw rate 0.58 across seeds
  1000-1007, chirality flips exactly on the two reflected trees).
- `conda run -n ssf python tests/test_cross_artifact_invariants.py` ->
  exit 0: full-pipeline fixture built in-process (28 s, 1 ok / 0
  failed), then 88 checks all PASS (scene 30/30 links resolve;
  stem coax 0.18-0.58 m; frustum volumes equal to 1e-3; all six
  trees anchored on the DTM at their (X, Y); shrub+grass HAG min
  0.0000).
- `conda run -n ssf python tests/run_all.py` -> exit 0, 113 PASS
  total.
- Fixture parity: the in-suite fixture reproduces the logged /tmp
  fixture numbers exactly (8,430 pts, 0 collisions / 0 trimmed,
  understory 103,297, stem map 597x654 px) — the promotion is
  value-identical to the suites the W7/whole-repo entries ran.
- Standing checks: check.py PASS (28 file(s) — check.py globs only
  top-level *.py, so tests/ is exercised by the suites themselves,
  not by the static pass).

**Convention going forward**: new verification rounds run
`tests/run_all.py` as part of their standing checks and add any new
invariant into these suites (same check/FAIL style); the
verification-log entry then references the suite check that covers
the round's behavior, so nothing lives only in /tmp again.

### Verification gaps / notes (W1 round)

- scene_editor's JSON acceptance is path-level (`find_plots`); its GUI was
  not exercised (no display).
- FIA geo passthrough (per-tree LON/LAT/ELEV from a PLOT join) is not yet
  derived from plot center + placement (item W5); this sample has no PLOT
  join (plot_center all-null).

### 2026-11-07 — W17: CSV-canonical tree lists (commit ed8838e)

notes.md design applied (it wins): the FIA helper's placed tree list is a
minimal CSV — tree rows only (TREE_UID, X, Y, Z, DBH, HT, TREE_CODE;
optional LON/LAT/ELEV never written until per-tree geo exists, item W5) —
and plot geometry lives in the report JSON the helper writes with it.

**Code**: `tree_list.py` gained `save_plot_csv`, `is_tree_list_csv` (header
sniff: required attrs present; legacy SST tables carry FUNC_GROUP/HT_LO and
never match), `_load_plot_csv` (plot fields from the sibling
`plot_<UID>_report.json` — parameters + a new top-level `plot_center`;
hand-authored CSVs without the sidecar default to circle/unknown and
callers override), with `load_plot`/`trees_frame` dispatching both
serializations. `fia_helper.py --format tree-list-csv` is the new default
(`tree-list-json` alternate, `sst-csv` legacy). `generate_trees.tree_requests`
and `replicate` accept CSV lists (kind label now `tree-list`).

**Bug found and fixed (root cause)**: first round-trip failed on DBH —
pandas' default fast float parser returns a value 1 ULP off for
17-significant-digit decimals in CSV text (46 in tree -> 1.1684 parsed
from "1.1683999999999999", vs 1.1683999999999999 in the JSON path), which
would break bit-exact regeneration (design principle 6). Fix:
`float_precision="round_trip"` in the CSV loader.

**Checks** (each exit-0 gated):

- `conda run -n ssf python check.py` -> exit 0 (28 files, no undefined
  names).
- FIA run x2 (CA_TREE_SMALL.csv, seed 42, radius 20) in both formats ->
  111 placed each; `load_plot` on both: trees frames `assert_frame_equal`,
  plot fields identical via the report sidecar (circle, 20.0);
  `trees_frame` passthrough identical; format sniff true for the minimal
  CSV, false for JSON -> exit 0.
- `tree_requests` CSV vs JSON with an SSTGenerator: kind `tree-list` both,
  111 requests, exact `==` on every field including the resolved `params`
  dict (only correct with the round_trip parser) -> exit 0.
- `conda run -n ssf python tests/run_all.py` -> exit 0 (all suites pass,
  standing run of the W7/placement+frame and cross-artifact invariants on
  the JSON-form fixture).

Known notes: report sidecar for CSV lists now carries `plot_center`
(extension to the notes.md report: parameters + counts + proportion);
scene_editor's GUI was not exercised (no display). Column order in the CSV
form is canonical (Z before DBH), but every consumer reads by name.

### 2026-11-07 — W19: rotation stepping + variable renames (commit 2ad83a9)

Owner rulings applied: placement rotation quantized "in steps of N" with
ROTATION_STEP default 30 degrees (12 orientations); shrub volume from the
sampled shape, not voxels.

**Code**:
- `assemble_plot.py`: rotation drawn from the discrete grid
  rng.integers(0, 360/step) * step degrees (asserted to divide 360; default
  30); `--rotation-step` CLI, `rotation_step_deg` recorded in the scene
  JSON, fanned out by `run_pipeline.py --rotation-step`.
- Vocabulary renames (defaults unchanged): run_pipeline/stem_map
  `--cell-size`+`--threshold` -> `--stem-resolution`+`--stem-section`
  (constants STEM_RESOLUTION/STEM_SECTION), canopy `--cell-size` ->
  `--chm-resolution` (CHM_RESOLUTION, now defined in canopy_map instead of
  borrowed from stem_map), terrain `--gen-cell-size` -> `--dtm-resolution`
  (parameter dtm_resolution), env-sidecar key `gen_cell_size` ->
  `dtm_resolution` (replicate.py, run_pipeline.py, scene_editor.py).
  Importers (restructure, shrub_outputs) updated.
- Shrub volume: `shrub_outputs.py` writes volume_m3 {primitive, mesh} only
  (mesh = the actual jittered primitive, the shrub volume per the ruling;
  voxel entry and voxel_volume_params removed); understory's
  --occupied-voxel flags removed (they served only the shrub voxel
  estimate); run_pipeline's --occupied-voxel flags remain, fanning to
  build_info only (foliage + neighbor-class accounting; the FOLIAGE_VOXEL/
  GRASS_VOXEL split is W16).
- Spec corrected first (doc-first): the earlier Section 9 note claiming
  the code used STEM_CELL_SIZE-style constants was wrong — the knobs are
  per-script CLI flags.

**Checks** (exit-0 gated):
- `conda run -n ssf python check.py` -> exit 0 (28 files).
- `conda run -n ssf python tests/run_all.py` -> exit 0: placement suite
  with a new standing section 7 (rotation grid: every recorded rotation a
  multiple of 30 deg in [0, 360); all 12 orientations across seeds
  2000-2007), plus the pre-existing determinism/replay/chirality checks
  all passing on the quantized draw; cross-artifact suite runs the full
  pipeline through the renamed flags (stem/canopy 5 cm, DTM 0.25,
  45-shrub understory), all invariants PASS.
- Shrub volume schema: invariant-suite fixture (--keep) shows per-shrub
  *_volume.json = {primitive, mesh} only.

Known notes: same-seed fixture values shifted vs pre-W19 runs (the draw
mechanism changed from uniform float to integer grid draws — expected, and
per-seed determinism is preserved); scene_editor's GUI was not exercised
(no display); CLI flag renames are breaking — old spellings are not
aliased (checked no other caller in-repo).

### 2026-11-07 — W19 review pass (commit df5ee88)

No bugs found beyond what the round closed. Checks (exit-0/exit-1 as
expected, each observed):

- Rotation grid assertion: `--rotation-step 25` -> AssertionError
  ("does not divide 360 evenly") before any tree processing; 90 accepted.
- `--help` on stem_map / canopy_map / terrain / understory -> all parse;
  the notes.md flags are listed; understory has zero occupied-voxel
  references.
- Replicate env-sidecar round-trip actually executed for the first time:
  `replicate.py` wrote `dtm_resolution: 1.0` into the sidecar, and
  `run_pipeline.py --understory --limit 1` consumed it end-to-end on the
  legacy SST table replicate (generation -> assemble -> stem/canopy ->
  restructure -> generated-DTM terrain from the sidecar's drawn params ->
  understory -> build_info -> normalize; 1 ok / 0 failed, 22 s). This also
  exercised the legacy-table path through the W17 minimal-CSV loader.
- Scene JSON from that run: `rotation_step_deg: 30.0` present;
  `rotations_deg` values on the grid.
- Stale-name grep: no live references to gen_cell_size / --gen-cell-size /
  --stem-cell-size / --stem-threshold / --canopy-cell-size /
  DEFAULT_THRESHOLD / DEFAULT_CELL_SIZE outside VERIFICATION_LOG's
  append-only history. scene_editor's terrain/understory CLI builders
  carry only current flags (static review; GUI not exercised, no display).
- Confirmed in passing: tests/fixture.py runs `generate_dtm=True`, so the
  armed cross-artifact suite already executes the renamed
  --generate-dtm/--dtm-resolution terrain path on a rugged DTM.

### 2026-11-07 — W14: notes.md semantic codes; trimmed dimension withdrawn (commit 77b8ac7)

Owner rulings applied: semantic 0 = ground and doubles as "no semantic
class needed" — ground and noise are the traditional LAS classification's
job, so the separate unassigned and noise semantic codes are gone (the
interim 255/7 proposal is withdrawn; 255 was invalid for int8 anyway);
no per-point `trimmed` dimension in final clouds — trimmed points are
removed, structural record only (the earlier W14 reinstatement of the
trimmed dim is itself withdrawn).

**Code**: schema.py only — SEM_GROUND 0, SEM_SHRUB 8, SEM_GRASS 9;
SEM_UNASSIGNED/SEM_NOISE deleted (their only consumers were the constant
definitions and the default fill, now SEM_GROUND). All ~50 name-based
consumers inherit; writers/readers iterate schema.PLOT_FIELDS. The
hardcoded numeric sites found and updated: tests' understory masks
(grass 6->0/9, understory pair [4,5]->[8,9], ground 6->0). colorize/
normalize pass values through generically (palette by observed ids).

**Docs**: spec Sections 3.2/3.4/7.5/8.3 + W14 row closed, registry
rulings; README_schema code table + stage bullets; CHECKLIST Section 0;
AGENTS conventions + register snapshot.

**Checks** (exit-0 gated):
- `conda run -n ssf python check.py` -> exit 0 (28 files; no dangling
  SEM_UNASSIGNED/SEM_NOISE references).
- `conda run -n ssf python tests/run_all.py` -> exit 0; new standing
  invariants in cross-suite section F (understory semantic domain subset
  {0,1,2,3,8,9}; shrubs carry 8; grass carries 9; ground carries 0) pass
  on fixture data; the armed run's normalize output reads
  semantic_code=[0, 1, 2, 3, 8, 9] from the on-disk plot cloud.

**Regeneration caveat**: clouds generated before this change carry the
old codes (4/5/6) — no migration; regenerate.

**Not verified**: scene_editor's GUI colorization paths (no display);
they read semantic_code generically per code review.

### 2026-11-07 — W14 review pass: two stale-code bugs found and fixed (commit 0d0f8b0)

The standing suites could not catch these (colorize_plot and the
scene_editor GUI are outside both), found by grepping for numeric-code
keys/masks:

- `colorize_plot.py` still mapped the old codes — ground 6, grass 5,
  shrub 4 — so post-W14 grass would stay black instead of grey and shrubs
  would get no per-id palette; ground was black only by luck (zeroed
  array). Fixed at the root: the masks now use the schema.SEM_*
  constants, so a future code change re-points the colors automatically.
- `scene_editor.py` SEM_COLORS/SEM_NAMES keyed 4/5/6 for shrub/grass/
  ground: shrub/grass points would have had no display color or layer
  toggle after the change. Same fix: keyed by schema constants
  (GUI itself still not exercised — no display; corrected by static
  review + import).
- Also refreshed normalize's documented categorical column names in
  README_understory (dynamic in code: `semantic_{c}` from observed codes,
  so post-W14 the columns are semantic_0/1/2/3/8/9) and README_shrub_
  outputs' "(semantic 4)" -> 8. build_info's `isin(sem, (1,2,3))` masks
  survive unchanged (codes 1/2/3 are the same in both schemes).

Checks (exit-0 gated):
- `conda run -n ssf python check.py` -> exit 0 (28 files) after both fixes.
- Colorize functional: cross-artifact suite run with `--keep`, then
  `colorize_plot.py` on the on-disk understory cloud (101,578 pts) and an
  independent point-order-paired assertion of the output RGB: grass
  (61,263) fully grey, ground (21,357) fully black, woody (2,067) fully
  brown, shrub (10,528) + foliage (6,363) fully randomized — exit 0.
- `conda run -n ssf python tests/run_all.py` re-run after the fixes ->
  exit 0, 119 PASS.

### 2026-11-07 — W14 second review pass: registry supersession annotation (commit cc7ff3c)

Aimed at anything the first pass could have missed, in particular the
0-collision (semantic 0 used to mean "unassigned"):

- `sst_generator.py` `old_sem == 0`: false alarm — that is treesim's own
  module-local raw label column (1 stem / 0 foliage), never the schema
  semantic code; the inline comment disambiguates.
- No live references anywhere to the old codes remain (numeric masks,
  dict keys, prose); `branch_id`'s "0 = unassigned" is a different
  dimension and correct.
- Media/ examples hold images only.
- One documentation-internal contradiction fixed: the historical
  Design-alignment ruling bullet in the Section 11 registry still said
  the per-point trimmed dimension "is reinstated" and unassigned/noise
  "need new codes" — now annotated as superseded by the W14 closing
  rulings (pointer annotations, W13 precedent).

No code changed in this pass, so no suite re-run (the standing run_all
exit 0 / 119 PASS remains the current result for the code).

### 2026-11-07 — W15: notes.md element graph embedded in skeleton JSONs (commit 42618e1)

The notes.md element graph (spec Section 7.4) is computed at generation
by sst_generator.element_graph and embedded under `elements` in every
skeleton JSON; trimmed variants carry the identical graph (the trimmed
writer already deep-copies the raw skeleton, and the elements carry no
coordinates — UIDs, SOURCE/DESTINATION, lengths, radii, cylinder
volumes — so they are place()-invariant by construction).

**Code**:
- `treesim.py` records two facts the graph needs, both from values
  already drawn (rng-neutral): each branch's `stem_index` (the drawn
  stem it was placed on) and `radius_m` (the ribbon half-thickness).
- `sst_generator.py` `element_graph()`: one order-0 element per trunk
  stem (root T00 SOURCE null; split daughters 2k+1/2k+2 source their
  parent stem), one order-1 element per branch sourced to its stem;
  radii per Section 7.4 (taper at stem base; ribbon half-thickness);
  volumes = pi r^2 L; DESTINATION = inverse of SOURCE.
- Polyline-and-splits form unchanged (trimming keeps consuming it).

**Checks** (exit-0 gated):
- check.py -> exit 0 (28 files; one caught slip: the new function was
  initially inserted under SSTGenerator's @register decorator, caught
  immediately by the module failing to register).
- Direct build smoke: root/branch links, reciprocity, volumes.
- Standing suite section J (new, 60 checks on the fixture): graph
  presence, single null-SOURCE root, SOURCE resolution and exact
  DESTINATION reciprocity, stem/branch partition, element count ==
  trunk segments + branches, lengths recomputed from polylines and
  branch base-tip pairs, stem radii follow the taper profile, volumes
  = pi r^2 L, trimmed variant identical. Fixture includes a forked
  tree (005: 3 stems). -> exit 0.
- Cascade coverage: forced SPLIT_PROB=1 with num_splits 2 and 3 -> 5-
  and 7-stem trees with granddaughter chains; reciprocity and
  acyclicity hold. (First attempt forced the wrong key names; the
  entry's fork probability gates the cascade — forcing SPLIT_PROB=1
  engages it.)
- tests/run_all.py -> exit 0, 179 PASS; the placement suite's
  bit-identical same-seed check and unchanged fixture point counts
  (e.g. 001 RW 3,534 pts) confirm the treesim additions drew no new rng
  values.

Pending (rides with W16): per-element `trimmed` flags on trimmed
variants; the structural record remains node_trimmed + scene counts.

### 2026-11-07 — Trim records removed everywhere (owner ruling) (commit de88920)

Ruling: once a point or branch has been trimmed, no record of it is
kept. Consequences implemented:

-assemble_plot: resolve_collisions returns nothing (no collision log);
  the scene JSON drops n_collisions/collisions/points_flagged_trimmed;
  trimmed skeletons lose node_trimmed, per-segment/per-branch "trimmed"
  flags and trim_summary, and are now physically the surviving geometry —
  trunk polylines are cut into runs at removed sub-segments (rebuilt
  segments with recomputed z-range/length/subsegment ids), trimmed
  branches are dropped, and the element graph is subsetted (dropped
  elements and scrubbed DESTINATION links). Cabin summary print no
  longer reports trim counts.
- build_info: single post-trim volumes everywhere (trunk frustum, branch
  cylinders, foliage voxel over the surviving cloud, canopy hull);
  trunk_volume/branch_volume lose their trimmed_only machinery; plot and
  tree infos carry plain floats. Latent bug fixed while in here: the plot
  info's shrub totals read a shrub_index.json that no stage writes
  (n_shrubs was always 0); they now sum the per-shrub volume JSONs
  (mesh = the shrub volume per the W19 ruling; primitive kept as
  analytic reference).
- canopy_map/restructure: single post-trim canopy-volume totals
  (canopy_map's _trimmed_volume and restructure's raw/trimmed pairs
  removed); run_pipeline's --drop-trimmed removed (dead: normalize never
  defined it); normalize drops the phantom "trimmed" extra field.

Checks (exit-0 gated):
- check.py exit 0; tests/run_all.py exit 0, 210 PASS.
- Placement suite section 8 (new): forced heavy collisions
  (branch_radius_frac 45, tolerance 3; 8,430 -> 2,067 points) — scene
  carries no collision log or removal counts; trimmed skeletons contain
  no node_trimmed/trim_summary/"trimmed"; surviving graph elements match
  their raw counterparts (DESTINATION may only shrink); trimmed branch
  lists shrink. (First version of this check wrongly asserted full
  element equality and caught its own error — DESTINATION legitimately
  shrinks — before being fixed.)
- Colliding two-tree plot (trees 0.5 m apart) through the full pipeline:
  tree/plot info volumes all single-valued floats, scene free of trim
  records; exit 0.
- Cross suite: section J restated to the subset semantics; new standing
  check that the plot info counts its shrubs (n=45, mesh total > 0) —
  this pins the shrub_index fix.
- Zero-collision standing fixture unchanged (identical values), so the
  removals altered no surviving geometry.

Known scope note: the raw/trimmed CLOUD pair (pre-/post-trim stage
variants) remains — deleting it is W6's output-design call, not a trim
record per se.

### 2026-11-07 — Trim-removal review pass: cut geometry verified (commit c8dfe17)

- Consumer sweep: nothing in the codebase reads the removed fields
  (collision log, removal counts, node_trimmed, trim_summary, per-branch
  flags); the only remaining mentions are the suite's intentional
  ban-list. scene_editor reads neither skeletons nor the collision log.
- Cut-geometry probe (the standing fixture has zero collisions, so this
  path was otherwise untested on real mid-cuts): 4 trees assembled with
  branch_radius_frac 45 / tolerance 3 on the legacy demo table. Every
  trimmed polyline run proved to be a contiguous span of the raw stem's
  edge lengths (rotation-invariant comparison); tree 001 lost 2.33 m of
  trunk across a 7-run split, its 167 branches dropped, its element
  graph 174 -> 7 with consistent SOURCE/DESTINATION; others lost whole
  crowns to the branch-vs-trunk rule while keeping intact trunks (runs
  1/1). The aggressive trimming also spot-checked the runs edge cases
  (middle cuts, distal cuts).
- Doc drift found and fixed: the scene's branch_order_encoding block and
  the assemble docstring claimed a las_field "branch_order" that no LAS
  carries; the real order field is branch_code (order info also lives in
  the skeleton JSONs). Pre-existing drift, corrected while here.
- No code-behavior change in this pass beyond the docstring/scene-doc fix
  (scene "branch_order_encoding" now names las_field branch_code and
  points at the skeleton fields). check.py exit 0; run_all.py re-run exit 0
  (210 PASS, 0 FAIL).

### 2026-11-10 — W16 closed: info JSON schema per notes.md (commit a5a686c)

- Plot and tree info JSONs now carry the notes.md field names as
  top-level keys; the old tree_characteristics / vegetation_characteristics
  nestings are gone. Extension data (plot uid/coord-system/units, canopy
  hull total, shrub primitive reference, realized grass cover,
  neighborhoods) lives in separate meta/canopy/grass/shrubs/neighborhood
  blocks. BOUNDING_BOX is the notes.md 6-float form
  [x_min, y_min, z_min, x_max, y_max, z_max].
- Shrub/grass/object volumes: VOL_SHRUB sums the per-shrub mesh volumes
  (sampled shape, as ruled), VOL_GRASS voxelizes at GRASS_VOXEL,
  VOL_FOLIAGE at FOLIAGE_VOXEL -- the shared --occupied-voxel-size knob
  split into --foliage-voxel/--grass-voxel with a shared
  --voxel-threshold (build_info CLI + run_pipeline flags; the invariant
  suite asserts VOL_SHRUB against an independent mesh recompute).
- Path fields resolve against the plot dir and the suite checks they
  exist: CHM_PATH (canopy_map_raw.tif), STEM_PATH (stem_map_raw.tif),
  PLAUSIBILITREE_CONFIGURATION (scene JSON; notes.md's
  "CONFIGURAITION" spelling corrected, owner-approved). PLOT_CENTER is
  present only when the FIA input carried lat/lon/elev (absent in the
  local fixture); DTM_PATH rides along when a plot DTM exists.
- TREE_Z = mean height-above-ground of the tree's points from the
  understory cloud (fallback: centroid z above the tree cloud's lowest
  point); TREE_X/TREE_Y stay the tree-list base position, TREE_DBH the
  taper diameter at 1.37 m, TREE_HT the skeleton height; VOL_TREE =
  trunk + branch + foliage (suite asserts the identity).
- Verification: check.py exit 0; tests/run_all.py exit 0 (221 PASS,
  0 FAIL); a --keep run of the invariant suite was inspected by eye
  (plot 11839: NUM_TREES 6, NUM_SHRUBS 45, PLOT_AREA = pi r^2 exact,
  VOL_TREES 114.013 == sum of tree VOL_TREE 69.299 + ...), then removed.
  The first full-suite run caught a residual dict-style bbox index in
  _shrub_neighborhood (fixed; suite then green end to end).

### 2026-11-10 — W16 review pass: consumer sweep + stale demo fixed (commit 0c6aefe)

- Consumer grep: nothing in the codebase reads the old
  tree_characteristics / vegetation_characteristics nestings or the
  dict-style info bbox; terrain.py/understory.py's dict-style bbox reads
  are canopy_bbox_json's own schema (bbox_world), unrelated to info
  JSONs. The only remaining old-name mentions are the spec's W16 register
  row (the gap description) and provenance entries.
- Cross-writer check against the pre-W16 numbers (same inputs, old
  committed demo): NUM_TREES, BASAL_AREA, VOL_TREES/TRUNK/BRANCH/FOLIAGE,
  canopy hull and bbox were bit-identical under the rename; the DIFFs
  found were all pre-W16-rounds staleness in the demo itself, not schema
  bugs.
- Root cause found and fixed: examples/plot_demo predated the W14/W19/
  trim-removal rounds (semantic codes 4/5/6, raw-vs-trimmed volume pairs,
  unquantized rotations, pre-W15 skeletons), so the new info writer was
  reading artifacts it was never specified against (VOL_GRASS collapsed
  to 0 because the old understory cloud has no semantic-9 points).
  Regenerated with the documented examples/README command (seed 42,
  radius 15, 256 ppm, --generate-dtm --understory; run_pipeline exit 0;
  24 s). Post-regen checks: normalize reports codes [0,1,2,3,8,9];
  scene/skeletons carry no trim records (only the raw/trimmed cloud-pair
  paths); VOL_SHRUB == per-shrub mesh sum (0.4571); plot sums ==
  per-tree values; VOL_TREE identity holds; BASAL_AREA == sum of
  pi(TREE_DBH/2)^2; PLOT_AREA == pi r^2; CHM/STEM/config/DTM paths
  resolve; TREE_Z in (0, TREE_HT]; element graphs present; PLOT_CENTER
  correctly absent (no FIA lat/lon in the demo input); build_info --help
  lists --foliage-voxel/--grass-voxel/--voxel-threshold.
- No code changes in this pass; run_all.py's last green run (221 PASS)
  still covers the unchanged code.

### 2026-11-10 — W18 closed: parameter dictionary per notes.md (commit 3b078eb)

- Dictionary schema: data/functional_groups.json now carries A0/A1/TAPER
  (owner ruling: TAPER over TAPER_B) plus DBH_MIN/DBH_MAX, HT_MIN/HT_MAX
  per code (DF/PP/TF/CE/RW/OK/SW/HW/DEFAULT; FIA-typical metre ranges).
  DBH_TOL/HT_TOL are withdrawn; realized sizes keep the fixed +/-1%
  generation constant in sst_generator (identical behavior, since no
  shipped entry ever carried tolerance keys).
- Applicability bounds enforced at entry resolution: species ->
  functional group -> DEFAULT, first in-bounds candidate wins, DEFAULT
  unconditional. generate_trees passes each tree's DBH/HT. Legacy
  CSV/JSON files with pre-W18 spellings rename at load
  (fia_helper.normalize_param_columns; param_editor uses the same), so
  old tables still open.
- SPLIT_NUM_LO/HI confirmed as a spec extension (owner ruling) — the
  cascading-hardwood look and suite section J's forced-cascade builds
  keep working.
- Verification: check.py exit 0; run_all.py exit 0 (233 PASS, 0 FAIL)
  with new invariant-section Hb (dictionary schema, withdrawn-key ban,
  ordered bounds, per-tree resolved entries, in/out-of-bounds and
  unknown-code resolution). Numerical inertness: regenerating the
  committed demo (seed 42) produced byte-identical LAZ/TIF/skeleton/
  scene/info files — only the four generation.json timestamps moved
  (reverted; plot_demo.csv is a legacy frozen-range table unaffected by
  the rename). First full-suite run caught a misdirected check
  (source_row exists only for legacy requests); replaced with a
  direct-resolution check.
