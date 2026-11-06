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

1. `check.py` — static, whole-repo.
2. Lint-for-mistakes greps used during review rounds:
   - stale identifiers: `grep -rn "allowed_overlap\|allowed-overlap\|7\.32" *.py *.md`
     (expect: no hits outside intentional provenance mentions)
   - doc cross-references resolve:
     `python3 -c "import re; t=open('Plausibilitree.md').read(); refs=set(re.findall(r'Section (\d+\.\d+)', t)); heads=set(re.findall(r'^### (\d+\.\d+)', t, re.M)); print([r for r in refs if r not in heads])"`
     (expect: empty list)

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

### Verification gaps / notes (W1 round)

- The JSON-tree-list e2e used 6 trees at 64 points/m (a 111-tree run is
  minutes of generation even at that density).
- scene_editor's JSON acceptance is path-level (`find_plots`); its GUI was
  not exercised (no display).
- FIA geo passthrough (per-tree LON/LAT/ELEV from a PLOT join) is not yet
  derived from plot center + placement (item W5); this sample has no PLOT
  join (plot_center all-null).
