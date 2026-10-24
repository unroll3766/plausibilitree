# FIA helper for PlausibiliTree

Converts raw FIA `TREE` records into one parameter table per plot, ready to
drive a multi-tree version of SimpleSynthTree. For what each parameter column
controls (the `CW_*` allometry, `TAPER_B`, crown shape ranges, etc.), see the
parameter reference in [README_param_editor.md](README_param_editor.md).

```
fia_helper.py
data/functional_groups.csv    per-group allometry + SST parameters
data/spcd_func_group.csv      SPCD -> functional group lookup
```

## Usage

```python
from fia_helper import fia_to_sst_plot, split_plots
import pandas as pd

df = pd.read_csv("CA_TREE_SMALL.csv", low_memory=False)
for plot_id, sub in split_plots(df, "PLT_CN"):
    fia_to_sst_plot(sub, plot_radius_m=30.0, allowed_overlap=0.3,
                    seed=42, out_csv=f"plot_{plot_id}.csv")
```

or from the shell:

```
python fia_helper.py CA_TREE_SMALL.csv -o plots -r 20 --allowed-overlap 0.3
```

This is a standalone step: run it first to produce `plots/plot_*.csv`, then
point `run_pipeline.py` at that directory to run the rest of the stages over
those tables (see README_run_pipeline.md). `run_pipeline.py` no longer runs
the FIA helper itself.

### Plot shapes

`--plot-shape` is `circle` (default), `rectangle`, or `square`:

- `circle` uses `-r/--radius` (defaults to 7.32 m, the FIA subplot radius).
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

| stage | function | notes |
|---|---|---|
| preprocess | `preprocess_fia` | drop `DECAYCD >= 1`, drop rows missing `DIA/HT/TPA_UNADJ/CCLCD/<group col>`, convert units, add `X`, `Y`, `successful_placement` |
| expansion | `expand_by_tpa` | **placeholder** — round to nearest int and duplicate |
| allometry | `crown_width` | `CW = exp(a0 + a1 ln DBH)` per functional group |
| placement | `place_trees` | random non-overlapping circles inside the plot |
| taper | `dbh_to_base_diameter` | DBH extrapolated to ground line |
| postprocess | `build_sst_table` | emits the columns below |

Output columns: `TREE_UID, FUNC_GROUP, X, Y, HT_LO, HT_HI, BASE_DIA_LO,
BASE_DIA_HI, SPLIT_HT_LO, SPLIT_HT_HI, SPLIT_PROB, NUM_BRANCH_LO,
NUM_BRANCH_HI, MIN_CAN_HT_LO, MIN_CAN_HT_HI, MAX_CAN_W, MAX_CAN_W_HT_LO,
MAX_CAN_W_HT_HI, TOP_DIST, MID_DIST, FP_RAND_HT, SPR_LO, SPR_HI`, plus
`plot_lat, plot_lon, plot_elev` when the input has them (see below).

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
`data/functional_groups.csv` and reference it from `data/spcd_func_group.csv`
to add a group; both files are plain CSV and hot-swappable via
`load_func_group_db(path)` / `load_spcd_map(path)`.

The shipped coefficients are plausible starting values, not fitted ones —
substitute published crown-width and taper coefficients for the species pool
you care about before using results quantitatively.

## Deviations from the written spec

1. **Height conversion.** The spec gives `HT * 3.6576`. FIA stores `HT` in
   feet, so the factor is `0.3048`; `3.6576` is 12× too large and produces
   150 m trees. `HT_TO_M` defaults to `0.3048`; pass `ht_factor=3.6576` to
   restore the specified behaviour.
2. **Crown width before placement.** The spec places trees first, then
   computes crown width, but the placement step is described as crown-cover
   based. Crown width is computed first so the packing primitive can be either
   the crown (`overlap_basis="crown"`, default) or the stem
   (`overlap_basis="stem"`).
3. **`MAX_CAN_W_HT_HT`** in the example output is read as a typo for
   `MAX_CAN_W_HT_HI`.
4. **TPA expansion** ignores the variable-radius/macroplot design; see
   Bechtold & Patterson (2005), GTR SRS-80.

## Known behaviour worth watching

Standard FIA subplot `TPA_UNADJ` is 23.6576 trees/acre, which over a 7.32 m
subplot is 0.983 trees — every such record rounds to exactly one stem, and
smaller expansion factors round to zero and disappear. The rounding rule is
therefore close to a no-op at subplot scale and only starts behaving like an
expansion at larger simulated plot radii. This is the main thing to replace
when the real estimation procedure goes in.

Strict crown non-overlap is unsatisfiable in closed-canopy stands: with
`allowed_overlap=0.0` some trees will exhaust `max_tries` and be dropped
(`successful_placement=False`). Real canopies interpenetrate, so
`allowed_overlap` around 0.3–0.5 is usually the more realistic setting.
