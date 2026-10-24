# Scene editor (`scene_editor.py`)

An interactive tool for tuning the environment of already-generated plots —
topography, ground, shrubs and grass — with immediate visual feedback. The
companion of the parameter editor: that one tunes the trees before they are
generated; this one tunes the world around them after.

```
python scene_editor.py plots                  # a run directory of plots
python scene_editor.py plots/plot_11574.csv   # a single plot table
```

## What it needs

Plots must be **generated and restructured** (the canonical layout with
`plot_placement/las/*_trimmed.laz`). `terrain.py` reads from the canonical
layout, so unrestructured plots are listed but cannot be edited; the status
bar says so.

## Workflow

1. Pick a plot (combobox, or Left/Right keys to page through the directory).
2. Edit the settings on the left.
3. **Generate** runs `terrain.py` and then `understory.py` on the selected
   plot as subprocesses — the same commands `run_pipeline.py` issues — and
   shows the resulting cloud. Stage output and errors land in the log pane at
   the bottom.
4. The **cloud** selector switches the view between `understory`, `terrain`
   and the trees-only `trimmed` cloud for comparison.
5. **Export settings** writes a `plot_<UID>.env.json` sidecar next to the plot
   table, in the same schema as `replicate.py`'s sidecars, so the environment
   can be re-applied by future pipeline runs of that directory.

Generate overwrites the plot's terrain / understory outputs — they are
derived from the trees and regenerable by definition. Tree outputs are never
modified.

## The settings

**Topography** (`terrain.py`):

| Control | Flag | Meaning |
|---|---|---|
| generate DTM / provided .asc | `--generate-dtm` / `--dtm` | synthesize a DTM per plot, or drape onto a supplied ESRI ASCII grid |
| cell / bump / relief | `--gen-cell-size` / `--gen-bumpiness` / `--gen-relief` | DTM resolution, roughness and height range (generated mode) |
| ground pts/cell, noise | `--density` / `--noise` | bare-earth sampling density and vertical jitter |
| source | `--source` | which tree cloud the terrain drapes (`trimmed`/`raw`) |

**Shrubs + grass** (`understory.py`):

| Control | Flag | Meaning |
|---|---|---|
| shrub cover, cell | `--shrub-cover` / `--shrub-cell-size` | fraction of cells hosting a shrub, cell grid size |
| width / height min-max | `--shrub-width-*` / `--shrub-height-*` | shrub size draw ranges |
| pts/shrub, margin | `--points-per-shrub` / `--shrub-margin` | point budget per shrub; clearance kept from tree stems |
| grass cover, patch | `--grass-cover` / `--grass-patch-scale` | fraction of candidates kept; vigour-field patch size |
| grass heights | `--grass-height-min/max` | grass point height range |
| cand/m2 | `--grass-candidate-density` | candidate points per m2 before thinning |
| attenuation | `--canopy-attenuation` | how strongly grass thins under canopy |
| seed | `--seed` | one seed drives both stages |
| skip per-shrub outputs | `--no-shrub-outputs` | don't write `individual_shrubs/` |

## Viewer

3D scatter of the plot cloud colored by semantic code — trunk brown, branch
stem tan, foliage green, shrub olive, grass yellow, ground gray — with the
legend showing per-class counts over the full cloud. A **layers** row of
checkboxes toggles each class on and off: hidden classes leave both the
scatter and the display subsample, so the remaining classes stay densely
sampled. Large clouds are subsampled to 120k points for display (seeded, so
the same cloud always looks the same; the loaded cloud is cached, so toggling
redraws instantly). Drag to rotate, scroll to zoom, **R** resets the camera.

## Environment sidecar

`Export settings` writes the same flat `env` object `replicate.py` uses
(see [README_replicate.md](README_replicate.md)), with a `dtm` entry when a
provided DTM is selected. `run_pipeline.py` picks sidecars up per plot
automatically on future runs of that directory.
