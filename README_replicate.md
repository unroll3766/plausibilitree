# Environment replicates (`replicate.py`)

`replicate.py` makes scene-level augmentation cheap: every replicate keeps the
**same inventory rows** (same trees, same positions, same placement seed) and
re-draws the **environment around them** — terrain shape, ground sampling, and
the shrub / grass layers. A model trained on the set then sees one canopy
under many different grounds.

```
# 3 environment replicates of one plot (or pass a directory to replicate all)
python replicate.py plots/plot_11574.json -o plots_rep -n 3 --seed 7
# (legacy SST plot_*.csv tables work identically)

# draw real DTMs from a directory instead of synthesizing terrain
python replicate.py plots -o plots_rep -n 4 --dtm-dir dtms

# run the pipeline exactly as usual; sidecars are picked up per plot
python run_pipeline.py plots_rep --radius 15 --seed 42
```

## What it writes

Per source table `plot_<UID>.<json|csv>` (canonical tree list or legacy SST
table) and replicate `k`:

| File | Contents |
|---|---|
| `plot_<UID>_r<NN>.<json|csv>` | the source table, copied byte-for-byte (same format; trees unchanged) |
| `plot_<UID>_r<NN>.env.json` | this replicate's environment draw (see below) |
| `replicates_manifest.json` | the sampler config + every draw, for provenance |

The sidecar is plain JSON with a flat `env` object whose keys mirror
`run_pipeline.py`'s kwarg names, so it doubles as a hand-run record — the
`terrain.py` / `understory.py` command for any single replicate can be
reconstructed from it directly.

## What varies, what doesn't

**Varies per replicate** (uniform draws over `min:max` ranges; a single value
means fixed):

| Setting | Default range | Key(s) in the env JSON |
|---|---|---|
| DTM bumpiness | 0.2 : 0.8 | `gen_bumpiness` |
| DTM relief | 3 : 15 m | `gen_relief` |
| DTM cell size | fixed 1.0 | `dtm_resolution` |
| Ground density | 0.5 : 2.0 pts/cell | `ground_density` |
| Ground noise | 0 : 0.03 m | `ground_noise` |
| Shrub cover | 0.05 : 0.35 | `shrub_cover` |
| Shrub width | 0.8 : 1.6 m (sorted pair) | `shrub_width_min` / `shrub_width_max` |
| Grass cover | 0.3 : 0.7 | `grass_cover` |
| Grass density | 100 : 300 pts/m2 | `grass_candidate_density` |
| Canopy attenuation | 0.2 : 0.5 | `canopy_attenuation` |
| Terrain/understory seed | drawn once per replicate | `seed` |
| DTM file (with `--dtm-dir`) | uniform choice of the `.asc` files | `dtm` |

Every range is overridable on the CLI, e.g.
`--bumpiness 0.4:0.9 --shrub-cover 0.1:0.5 --relief 5:25`.

**Shared across replicates** (identical by design): the tree tables
themselves, `generate_trees` / `assemble_plot` / `stem_map` / `canopy_map` /
`restructure` (all driven by `run_pipeline`'s global `--seed`), and the
dataset-wide `normalize` call (one seed for the whole batch, so subsampling
does not differ between replicates).

## How the orchestrator applies it

`run_pipeline.py` looks for `plot_<UID>.env.json` next to each table and, if
present, uses that draw **for that plot's terrain and understory stages only**
— every other stage, and every sidecar-free plot, keeps the global settings.
A sidecar `dtm` entry replaces DTM synthesis for that plot. A mixed directory
(some replicated, some not) works unchanged; stages that aren't configured
for a plot are skipped with a note.

One caveat on mixed directories: `normalize.py` requires every plot in a
batch to resolve to the same source kind (understory vs. trimmed vs. raw), so
a batch mixing plots with and without an understory layer errors out with
instructions. Keep replicates in their own run directory, or pass
`--norm-source` explicitly to force one kind for the whole batch.

## Reproducibility

Given the same source tables and the same `--seed`, `replicate.py` produces
identical draws (the manifest records the seed and ranges). Each replicate's
own seed then drives its terrain and understory stages, so a re-run of
`run_pipeline.py` reproduces the same environments. To re-draw, change
`replicate.py --seed` and re-run both commands.
