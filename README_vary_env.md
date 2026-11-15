# Per-plot environment draws (batch variation)

`vary_env.py` is the batch counterpart of `replicate.py`. Where
`replicate.py` makes N replicates of **one** plot (same trees, re-drawn
environment), `vary_env.py` varies the environment **across** a directory of
distinct plot tables: every plot gets its own terrain / shrub / grass draw.

```
python vary_env.py -i plots -o plots_var --seed 7
python run_pipeline.py plots_var --seed 42     # sidecars are applied per plot
```

Inputs go to `-i/--input` (repeatable, may be tables or directories;
positional paths work too and are combined with `-i`), output to
`-o/--out-dir` (default `plots_var`). A table that fails -- unreadable,
not a canonical plot table, missing plot geometry -- is skipped with a
warning and recorded in the manifest's `errors` list; one bad table does
not cancel the batch, but the run exits nonzero if any table failed.

## What it writes

```
plots_var/
  plot_<UID>.csv|json       byte copy of the source table
  plot_<UID>_report.json    byte copy of the FIA report sidecar (CSV tables
                            read their plot geometry from it)
  plot_<UID>.env.json       this plot's environment draw (picked up
                            automatically by run_pipeline.py's _load_env)
  vary_env_manifest.json    sampler config + every draw + provenance
```

The output directory must be empty unless `--force` is given. Tables are
copied byte-for-byte; the trees are not changed, and `run_pipeline.py`'s
global `--seed` still drives tree generation/assembly — variation is
environment-only (terrain, shrubs, grass), exactly as in replicate.py.

## The draws

- **Terrain / ground** (`gen_bumpiness`, `gen_relief`, `ground_density`,
  `ground_noise`, `dtm_resolution`, optional `--dtm-dir` DTM draw,
  `shrub_cell_size`, `shrub_width` pair, `grass_candidate_density`,
  `canopy_attenuation`, per-plot `seed`): uniform over CLI `--min:max`
  ranges, reusing replicate.py's `DEFAULT_RANGES` / `parse_range` /
  `draw_env`. `--shrub-cover` / `--grass-cover` ranges are deliberately not
  offered — those two keys come from the joint draw below.
- **Shrub + grass cover (joint, tree-conditioned).** With
  `BA = Σ π(DBH/2)²` from the plot table's DBH column and `A` the plot area
  (from the table / its report sidecar; CLI `--plot-shape`/`--radius`/
  `--plot-length`/`--plot-width` are the fallback for geometry-less
  tables), the total understory cover
  `T = shrub_cover + grass_cover ~ U[min(BA/A, 0.5), max(1 − BA/A, 0.5)]`
  (clamped at BA/A = 0.5, where the interval would be empty), split
  uniformly: `shrub_cover = T·b`, `grass_cover = T·(1−b)`, `b ~ U(0,1)`.
  Rounding keeps the two covers summing to T exactly at 4 dp. The bound is
  on the **nominal** parameters — realized understory area sits below T
  because grass is attenuated under canopies (`canopy_attenuation`) and
  shrubs keep a tree margin.

## Provenance and determinism

Each plot's seed derives from the run `--seed` plus the plot UID
(`per_plot_seed`, SHA-256 based), so draws are reproducible and independent
of table order (spec principle 6). The manifest records, per plot: source
table, copied table/report/env paths, BA, plot area, BA/A, the cover bounds
used, T, the split b, both covers, the full env draw and the geometry used.
Same seed → byte-identical sidecars; different seed → different draws.

## Tests

`tests/test_review_round_regression.py` section F (W30): bounds/sum/clamp
of the joint draw over several BA/A values, per-plot seed determinism, env
keys within run_pipeline's ENV_KEYS, and a two-plot end-to-end (copied
tables + report sidecars + sidecars + manifest, `_load_env` pickup,
manifest BA/A recomputed from the table, seed reproducibility, non-empty
out-dir refusal).
