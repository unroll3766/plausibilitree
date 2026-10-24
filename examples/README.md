# Examples

`plot_demo` is a small, fully-processed demo plot (committed to the repo) so
you can explore the pipeline's outputs — or the [scene editor](README_scene_editor.md)
— without generating anything first:

```
python scene_editor.py examples        # open the demo in the scene editor
python colorize_plot.py examples/plot_demo/plot_placement/las/plot_demo_understory.laz out.laz
```

It contains four trees — two cascading hardwoods (`SPLIT_NUM 1:5`) and two
conifers — placed at readable spacing, with the complete pipeline run over
them: per-tree clouds and skeletons, assembled and trimmed scene, stem and
canopy rasters, generated terrain, shrubs and grass, info JSONs, and
normalized clouds. The full canonical layout is present; see
[README_schema.md](../README_schema.md) for what each file contains and
[README_run_pipeline.md](../README_run_pipeline.md) to reproduce it:

```
python run_pipeline.py examples --radius 15 --seed 42 --points-per-metre 256 \
    --generate-dtm --understory
```

Pipeline outputs you generate elsewhere (e.g. a `plots/` run directory) are
git-ignored; `examples/` is the only committed plot.
