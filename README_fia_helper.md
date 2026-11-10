# FIA helper for PlausibiliTree

Converts raw FIA `TREE` records into one placed tree list per plot, ready to
drive a multi-tree version of SimpleSynthTree. By default it writes the
minimal tree-list CSV (`plot_<UID>.csv`, spec Section 8.1 — tree rows only:
TREE_UID, X, Y, Z, DBH, HT, TREE_CODE); plot geometry lives in the sibling
`plot_<UID>_report.json`, which also records the placement parameters and
counts. `--format tree-list-json` emits the alternate canonical JSON object;
`--format sst-csv` emits the legacy per-tree SST parameter table instead
(frozen HT/BASE_DIA ranges plus per-group shape parameters). For what each parameter controls (the `CW_*`
allometry, `TAPER`, crown shape ranges, etc.), see the parameter reference
in [README_param_editor.md](README_param_editor.md).

```
fia_helper.py
data/functional_groups.json   per-group allometry + SST parameters
                              (the parameter dictionary, spec Section 8.6)
data/spcd_func_group.csv      SPCD -> functional group lookup
tree_list.py                  the canonical tree-list contract (Section 8.1)
```

## Usage

```python
from fia_helper import fia_to_tree_list, split_plots
import pandas as pd

df = pd.read_csv("CA_TREE_SMALL.csv", low_memory=False)
for plot_id, sub in split_plots(df, "PLT_CN"):
    plot = fia_to_tree_list(sub, plot_radius_m=30.0, expansion_modifier=1.0,
                             min_distance=1.0, edge_effects=0.25,
                             seed=42, plot_uid=f"plot_{plot_id}")
    # plot: dict with plot shape/dims, plot_center, and a minimal trees array
```

or from the shell:

```
python fia_helper.py CA_TREE_SMALL.csv -o plots -r 20 --min-distance 1.0 --edge-effects 0.25
#    -> plots/plot_<UID>.csv (minimal tree list) + plot_<UID>_report.json
python fia_helper.py CA_TREE_SMALL.csv -o plots -r 20 --format tree-list-json
#    -> plots/plot_<UID>.json (canonical JSON tree list) + report
python fia_helper.py CA_TREE_SMALL.csv -o plots -r 20 --format sst-csv
#    -> plots/plot_<UID>.csv (legacy SST parameter table) + report
```


This is a standalone step: run it first to produce the tree lists (plus a
`plot_*_report.json` per plot), then point `run_pipeline.py` at that
directory to run the rest of the stages over those lists (see
[README_run_pipeline.md](README_run_pipeline.md)). `run_pipeline.py` no
longer runs the FIA helper itself. For per-row SST parameter access (the
legacy table) the `fia_to_sst_plot` API remains available and takes the
same placement controls.

### Placement controls

- `--expansion-modifier` (default 1.0): multiplied into FIA-recorded TPA
  before rounding **up** to the expected per-plot count. 0 = only the
  documented trees; values in (0, 1) = global reduction; > 1 = more trees
  than the recorded TPA. (Every tallied record always places at least one
  tree.)
- `--min-distance` (default 0.0): required surface-to-surface separation
  between two trees' discs, in metres. Negative = discs may interpenetrate
  by up to `|d|`; 0 = trees cannot touch; positive = a required buffer.
  The disc diameter is `DBH + min_distance` for the default `--basis stem`
  (the design spec's trunk disc). This generalizes and replaces the old
  fractional `--allowed-overlap`.
- `--basis` (default `stem`): `stem` packs trunk discs (spec default);
  `crown` packs crown-width discs instead.
- `--max-tries` (default 100): placement attempts per tree before it is
  dropped from the output.
- `--edge-effects` (default 0.0, range 0-1): allowed proportion of a tree's
  disc area outside the plot. 0.0 = the disc must fall entirely inside the
  plot; 1.0 = any sampled position (positions are always sampled inside the
  plot, so a disc is never fully outside).
- `--seed`: seeds the random placement draws and the ordering tiebreak.

Trees are placed tallest-first (DBH tiebreak, seeded random for exact
ties). Only successfully placed trees are written (with fresh sequential
`TREE_UID`s); each plot's `plot_*_report.json` records the parameters used
and the placed / not-placed / proportion counts. Unplaced trees in dense
plots are normal: raise `--min-distance`'s negative range, `--max-tries`,
`--edge-effects`, or the plot size if you need more tries to succeed.

### Plot shapes

`--plot-shape` is `circle` (default), `rectangle`, or `square`:

- `circle` uses `-r/--radius` (defaults to 7.315 m, the FIA subplot radius).
- `rectangle` uses `--length` (x) and `--width` (y); both are required.
- `square` uses `--side`, a single side length; required.

`--radius` only applies to `circle` and is rejected (a CLI error) alongside
`--plot-shape rectangle`/`square`, since the plot area/placement geometry for
those instead comes from `--length`/`--width` or `--side`. Likewise
`--length`/`--width`/`--side` are rejected with `--plot-shape circle`.
Programmatically, `plot_area_m2`/`plot_bbox`/`preprocess_fia`/`place_trees`/
`fia_to_sst_plot` all take a `plot_side_m` argument alongside the existing
`plot_radius_m`/`plot_length_m`/`plot_width_m` for the `square` case.

## Pipeline

Stages below run for all output formats; the postprocess stage diverges:
`--format tree-list-csv` (default) and `--format tree-list-json` both run
`build_tree_list` and emit the minimal tree list (Section 8.1 — one CSV row
per tree, or the JSON object); `--format sst-csv` runs
`build_sst_table` and emits the frozen-parameter columns below. In the
tree-list paths, the tolerance/allometry values are computed at *generation*
time by
`generate_trees.py` from the tree list's DBH/HT, so the tree list carries
only inventory facts.

| stage | function | notes |
|---|---|---|
| preprocess | `preprocess_fia` | drop `DECAYCD >= 1`, drop rows missing `DIA/HT/TPA_UNADJ/CCLCD/<group col>`, convert units, add `X`, `Y`, `successful_placement` |
| expansion | `expand_by_tpa` | `N = ceil(TPA * plot area * expansion modifier)`, floored at 1; duplicates rows to match |
| allometry | `crown_width` | `CW = exp(a0 + a1 ln DBH)` per functional group |
| placement | `place_trees` | height-descending placement (DBH tiebreak), disc spacing with `min_distance`, `edge_effects` quota, `max_tries` attempts |
| taper | `dbh_to_base_diameter` | DBH extrapolated to ground line |
| postprocess | `build_tree_list` / `build_sst_table` | minimal tree list (CSV default, JSON alternate) / the legacy columns below |

Legacy (sst-csv) output columns: `TREE_UID, FUNC_GROUP, X, Y, HT_LO, HT_HI,
BASE_DIA_LO, BASE_DIA_HI, SPLIT_HT_LO, SPLIT_HT_HI, SPLIT_PROB,
NUM_BRANCH_LO, NUM_BRANCH_HI, MIN_CAN_HT_LO, MIN_CAN_HT_HI, MAX_CAN_W,
MAX_CAN_W_HT_LO, MAX_CAN_W_HT_HI, TOP_DIST, MID_DIST, FP_RAND_HT, SPR_LO,
SPR_HI`, plus `plot_lat, plot_lon, plot_elev` when the input has them
(see below).

`X`/`Y` are metres relative to plot centre. `HT_*` and `BASE_DIA_*` are metres;
`SPLIT_HT_*`, `MIN_CAN_HT_*` and `MAX_CAN_W_HT_*` are fractions of tree height,
matching `treesim.gen_simtree`'s `model_params`.

### Plot-center pass-through (`plot_lat`/`plot_lon`/`plot_elev`)

If the input has FIA `LAT`/`LON`/`ELEV` columns (e.g. the TREE table was
joined with the PLOT table before being passed in), `preprocess_fia` copies
them straight through, unchanged, as `plot_lat`/`plot_lon`/`plot_elev`. No
unit conversion and no offset is applied — these are the plot-center
location, unlike `X`/`Y`, which are local placement coordinates relative to
that center. Columns not present in the input are simply omitted from the
output rather than added as empty/NaN columns.

`assemble_plot.py` reads these three columns back out of the plot table and
carries them into `_scene.json` (see README_assembly.md), and `build_info.py`
reads them from there into `plot_info.json`'s `meta.plot_center` (see
README_build_info.md) — so a plot's lat/lon/elev reach the pipeline's
textual/JSON output, not just the intermediate CSV.

## Functional groups

Resolution order is SPCD lookup → SPGRPCD coarse split (1–24 softwood,
25–48 hardwood) → `DEFAULT`. Any unmatched key silently falls back to
`DEFAULT`, so an incomplete table is never fatal. Add a row to
`data/functional_groups.json` and reference it from `data/spcd_func_group.csv`
to add a group; the dictionary is JSON and the lookup is CSV, both
hot-swappable via
`load_func_group_db(path)` (JSON or legacy CSV) / `load_spcd_map(path)`.

The shipped coefficients are plausible starting values, not fitted ones —
substitute published crown-width and taper coefficients for the species pool
you care about before using results quantitatively.

## Deviations from the written spec

1. **Height conversion.** The spec gives `HT * 3.6576`. FIA stores `HT` in
   feet, so the factor is `0.3048`; `3.6576` is 12× too large and produces
   150 m trees. `HT_TO_M` defaults to `0.3048`; pass `ht_factor=3.6576` to
   restore the specified behaviour.
2. **Crown width before placement.** The spec description computes crown
   width after placement, but the placement primitive is one of the two
   `overlap_basis` discs, so crown width is computed first. `overlap_basis
   ="stem"` (trunk discs, the design default) is now the default;
   `"crown"` packs crown-width discs instead.
3. **`MAX_CAN_W_HT_HT`** in the example output is read as a typo for
   `MAX_CAN_W_HT_HI`.
4. **TPA expansion** ignores the variable-radius/macroplot design; see
   Bechtold & Patterson (2005), GTR SRS-80.

## Known behaviour worth watching

Standard FIA subplot `TPA_UNADJ` is 23.6576 trees/acre, which over a 7.315 m
subplot is ~1.0 tree — records land at the floor of one stem each, and the
expansion modifier only starts adding stems at larger simulated plot radii
(or with per-acre TPA above the subplot constant). Use `--expansion-modifier
> 1` or a bigger plot to exercise duplication.

Strict disc non-overlap is unsatisfiable in closed-canopy stands: with
`min_distance >= 0` some trees will exhaust `max_tries` and be dropped
(they are omitted from the output and counted in the report JSON). Real
stems and canopies interpenetrate, so a small negative `--min-distance`
(or `--basis crown` with a negative `--min-distance`) is usually the more
realistic setting.
