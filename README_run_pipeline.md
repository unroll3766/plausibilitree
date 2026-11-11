# Pipeline orchestrator (`run_pipeline.py`)

Runs the PlausibiliTree pipeline over a directory of already-generated
per-plot SST tables: generate_trees, assemble_plot, restructure,
(always) terrain, (optional) understory (grass + shrubs), then build_info —
for every plot, all with the same parameters. The derived stages — stem
map and canopy map — are OPT-IN (spec Section 5): they run only
when `--stem-map` and/or `--canopy-map` are given. The FIA
helper (`fia_helper.py`) is a separate, standalone step that must be run
first to produce those tables; see README_fia_helper.md.

```
python fia_helper.py CA_TREE_SMALL.csv -o plots -r 15 --seed 42

python run_pipeline.py plots \
    --radius 15 --seed 42 \
    --points-per-metre 4096 --max-attempts 10 \
    --stem-map --stem-resolution 0.1 --stem-section 0.25 \
    --canopy-map --chm-resolution 0.1 \
    --generate-dtm --understory
```

This reproduces the documented hand-run sequence:

```
fia_helper.py CA_TREE_SMALL.csv -o plots -r 15 --seed 42
  for each plots/plot_*.csv:
    generate_trees.py  --seed 42 --points-per-metre 4096 --max-attempts 10
    assemble_plot.py   --seed 42 --plot-radius 15
    stem_map.py        --stem-resolution 0.1 --stem-section 0.25
                       (opt-in: --stem-map)
    canopy_map.py      --chm-resolution 0.1 --source trimmed
                       (opt-in: --canopy-map)
    restructure.py
    terrain.py          --generate-dtm ...       (optional; --dtm or --generate-dtm)
    understory.py        ...                     (optional; requires terrain, --understory)
    build_info.py        ...                     (always; plot + per-tree info JSONs)
```

Every default matches the first example (terrain/understory and the derived
stages off by default), so the flags above are only needed to change them;
`python run_pipeline.py plots` alone runs the CORE pipeline (plot clouds +
JSONs, no maps) once `plots/` has been populated by
`fia_helper.py`.

## How it works

Each stage is invoked as its own subprocess — exactly like running it by hand —
so each stage's own arguments, one-line summary and error handling apply
unchanged. The orchestrator discovers the plot tables already sitting in the
given directory (`<plotdir>/plot_*.csv`, produced beforehand by `fia_helper.py`)
and feeds each through the stages in the required order (restructure before
terrain/understory/build_info, which the stages depend on).
It errors out with a reminder to run `fia_helper.py` first if no `plot_*.csv`
is found there.

Extraction (`--extract-trees`) also runs per plot, after build_info; its two
knobs are `--extract-source` (plot cloud the extraction slices; 'understory'
resolves to the furthest-processed present, default) and
`--extract-surroundings` (metres of context added around each object's own
points, default 0). Normalization is retired (2026-11 review): the DL side
owns it.

## Terrain and understory (grass + shrubs)

Both are optional per-plot stages that slot in after `restructure`, before
`build_info`. Full detail lives in README_terrain.md and
README_understory.md; this is the `run_pipeline.py` surface.

**Terrain** drapes the plot onto a DTM and adds ground points. Enable it with
either `--dtm <file.asc>` (use a provided DTM) or `--generate-dtm` (synthesize
one per plot, auto-sized to that plot's footprint) — they're mutually
exclusive. Omit both and the step is skipped entirely.

**Understory** (grass + shrubs) is enabled with `--understory`, but it reads
the draped terrain cloud, so it **requires terrain**: pass `--understory`
without `--dtm`/`--generate-dtm` and the orchestrator prints a note and skips
it rather than failing the run. Shrubs are discrete jittered-ellipsoid meshes
placed in the gaps between tree canopies (`--shrub-cover` of the eligible grid
cells, sized by `--shrub-cell-size`/`--shrub-width-min`/`--shrub-width-max`);
grass is a continuous point field thinned by `--grass-cover`,
`--grass-candidate-density`, and `--canopy-attenuation` (how much cover a
tree/shrub canopy suppresses beneath it).

`--foliage-voxel`/`--grass-voxel`/`--voxel-threshold` are used by `build_info.py`
(FOLIAGE_VOXEL and GRASS_VOXEL volume accounting), which always runs last
(after every cloud/skeleton/map, including terrain/understory when present)
to emit the plot- and tree-level `*_info.json` files — see
README_build_info.md.

`--radius` (and `--plot-shape`/`--length`/`--width`) should match whatever was
passed to `fia_helper.py` when the plot tables were produced, since it feeds
`assemble_plot --plot-radius` and the assembly geometry: the placement radius
used there needs to agree with the tree budget baked into the tables by the
FIA helper.

## Parameters

Grouped by the stage they feed. Anything omitted uses the default in
parentheses.

- shared / assemble_plot geometry: `--radius (15)`, `--plot-shape (circle)`,
  `--length`, `--width` (rectangle only), `--seed (42)`,
  `--limit` (first N plot tables found in plotdir)
- generate_trees: `--points-per-metre (4096)`, `--max-attempts (10)`
- assemble_plot: `--branch-radius-frac (0.10)`, `--tolerance (0.0)`,
  `--rotation-step (30 — rotations are quantized to multiples of it)`
- stem_map: `--stem-resolution (0.05)`, `--stem-section (0.25)`,
  `--bbox-padding (1)`
- canopy_map: `--chm-resolution (0.05)`, `--canopy-buffer (2)`,
  `--canopy-source (trimmed)`
- terrain (optional; off unless `--dtm` or `--generate-dtm` is given):
  `--dtm PATH` (provided DTM .asc), `--generate-dtm` (synthesize one per plot;
  mutually exclusive with `--dtm`), `--dtm-resolution (0.25)`,
  `--gen-bumpiness (0.5)`, `--gen-relief (10.0)`, `--ground-density (1.0)`,
  `--ground-noise (0.0)`, `--terrain-source (trimmed)`
- understory / grass + shrubs (optional; `--understory` flag, requires
  terrain): `--shrub-cover (0.15)`, `--shrub-cell-size (2.0)`,
  `--shrub-width-min (0.8)`, `--shrub-width-max (1.6)`, `--grass-cover (0.6)`,
  `--grass-candidate-density (200.0)`, `--canopy-attenuation (0.3)`
- build_info (always runs, last):
  `--foliage-voxel (0.1)`, `--grass-voxel (0.1)`, `--voxel-threshold (1)`

## Error handling and control

- **Default: skip and continue.** A stage failure logs the error, skips the rest
  of that plot, and moves to the next. The run ends with a summary
  (`N ok, M failed`) and exits non-zero (3) if any plot failed, so a batch of
  hundreds isn't lost to one bad plot.
- **`--stop-on-error`** aborts the whole batch on the first stage failure
  (exit 2).
- **`--dry-run`** prints every command without running the stages. Stale
  `plot_*/` stage directories from a previous orchestrator run are still
  cleared unless `--dry-run` is also given.
- **`--limit N`** processes only the first N plot tables found in `plotdir`.
- **`--understory` without `--dtm`/`--generate-dtm`** doesn't fail the run —
  it prints a note (`--understory requires terrain...; skipping understory.`)
  and proceeds without the understory stage.

## Notes

- `fia_helper.py` is a separate step and owns the `plot_*.csv` tables in
  `plotdir`; the orchestrator only reads them. It does clear stale `plot_*/`
  stage directories from a previous orchestrator run on each invocation, but
  never touches `plot_*.csv` — rerun `fia_helper.py` to regenerate those.
- The normalization stage is retired (2026-11 review): the DL side owns it.
