# Feature checklist

Everything a from-scratch reimplementation needs to cover, module by module.
Every item below exists in the code today; the checkbox is for tracking a
rebuild or port. Read [README.md](README.md) first for the overall workflow.

## 0. Core vocabulary (build first — everything depends on it)

- [ ] **Canonical per-point schema** (`schema.py`): x/y/z plus 7 extra dims —
      `tree_id` (i32), `shrub_id` (i32), `semantic_code` (i8: 0 ground / no
      class needed, 1 trunk, 2 branch stem, 3 foliage, 8 shrub, 9 grass;
      no separate unassigned/noise codes — LAS classification handles them),
      `classification` (u8: 2 ground, 3/4/5 low/med/high veg, 7 noise),
      `branch_code` (i8, branch order + 1), `branch_id` (i32),
      `height_above_ground` (f32). One schema for plot and tree clouds; dims
      not yet computed are written as their defaults so the schema never
      changes shape across the chain.
- [ ] **Canonical file layout** (`layout.py`): key -> (subdir, filename
      template) for plot / tree / shrub artifacts, plus path and relative-link
      helpers.
- [ ] **LAS I/O** (`las.py`): canonical export/read at the schema's dtypes;
      generic export with caller-specified extra dims; point format 0
      (no RGB), scale 0.001, offsets at the cloud minimum.
- [ ] **Scene JSON**: one file per plot recording geometry, seeds,
      id-encoding documentation, and relative links to every
      artifact.

## 1. Tree list preparation (optional entry point)

- [ ] FIA `TREE` table -> one `plot_<UID>.csv` per plot (`fia_helper.py`):
      preprocess -> TPA expansion (explicit placeholder) -> crown width from
      functional groups -> placement -> postprocess; each step usable alone.
- [ ] Species -> functional-group lookups from companion CSVs
      (`data/functional_groups.json`, `data/spcd_func_group.csv`).
- [ ] Plot-center LAT/LON/ELEV passthrough columns; `--radius` and `--basis`
      options; seeded placement.
- [ ] Any other source can write the same `plot_*.csv` column contract instead
      (FIA data is not required).
- [ ] **Environment replicates** (`replicate.py`): N copies of a plot table,
      each with a randomly drawn terrain/understory sidecar JSON
      (`plot_<UID>_r<NN>.env.json`) consumed automatically by the pipeline
      orchestrator; draws are uniform over CLI-given `min:max` ranges, with an
      optional `--dtm-dir` to draw real DTMs, and a manifest records every
      draw for provenance.

## 2. Point-cloud generation

- [ ] **Tree simulator** (`treesim.py`, vendored Simple Synth Tree): spline
      trunk + branch mesh, surface-sampled points, stem/foliage per-point
      labels, tunable model parameters. Forked trees support **cascading
      trunk splits** (the `SPLIT_NUM_LO/HI` functional-group columns, default
      1 = the original single fork): each split truncates a live stem and
      grows two daughters, with the skeleton recording the split cascade.
- [ ] **Per-tree generation** (`generate_trees.py`): table row -> simulator
      params, density scaled by tree size (points per metre), retry loop,
      writes per-tree `_raw.las` + raw skeleton JSON (trunk polylines and
      segments, branch base/tip, and the notes.md element graph under
      `elements`: one order-0 element per trunk stem, one order-1 element
      per branch, SOURCE/DESTINATION links, cylinder pi*r^2*L volumes;
      trimmed variants carry the same, place()-invariant graph)
      + generation JSON.
- [ ] **Interactive parameter editor** (`param_editor.py`): GUI to edit any
      cell of a parameter CSV (functional groups or similar), generate N
      sample trees from the selected row, page through them in a 3D viewer
      (keys/buttons), and save the edited table to a new CSV; heights and
      diameters are sampled from GUI ranges and completed via the same
      crown-width / taper allometry as fia_helper.
- [ ] **Plot assembly** (`assemble_plot.py`): random z-rotation per tree about
      its own base, quantized in `--rotation-step` degree steps (default 30),
      placement from table X/Y, merged raw scene LAS.
- [ ] **Collision trimming**: skeleton segment-segment distances with branch
      radii; trim rules (branch loses to trunk, shorter branch/stem loses,
      seeded and UID-based tiebreaks); the whole distal subtree is trimmed;
      trimmed points are removed, not flagged.
- [ ] **Per-point attribution**: points assigned to skeleton elements (stem
      points -> trunk elements, foliage -> branch elements); plot-unique
      compact `branch_id` + globally-unique string uid; branching order;
      semantic refinement (main trunk vs fork stems).
- [ ] **DTM synthesis** (`generate_dtm.py`, standalone CLI and library):
      diamond-square fractal terrain, ESRI ASCII grid output.
- [ ] **Terrain** (optional, `terrain.py`): drape onto a provided or generated
      DTM, per-tree rigid vertical drop, ground-point sampling with noise, LAS
      classification by height above ground, HAG column, per-tree terrain crops
      (bbox + neighboring points).
- [ ] **Understory** (optional, requires terrain; `understory.py` +
      `vegetation.py`): seeded shrub layer (jittered ellipsoid meshes, avoided
      around tree canopies) and grass layer (Poisson thinning of a
      diamond-square vigor field, attenuated under canopy), per-shrub ids,
      per-tree understory crops.
- [ ] **Interactive scene editor** (`scene_editor.py`): GUI to tune terrain,
      shrub and grass settings on already-generated (restructured) plots,
      regenerate the environment via the stage CLIs, view the resulting cloud
      colored by semantic code, and export the settings as an env sidecar in
      the replicate schema.
- [ ] **Per-shrub outputs** (`shrub_outputs.py`): shrub raw/understory clouds,
      volume JSONs (primitive / mesh, the sampled-shape volumes), per-shrub
      canopy and stem maps.

## 3. Labels and outputs for ML/DL

- [ ] **Stem map** (`stem_map.py`): breast-height (1.37 m +/- threshold)
      point-count raster (GeoTIFF), per-tree stem bounding boxes (pixel +
      world) JSON, per-tree crops and semantic masks.
- [ ] **Canopy map** (`canopy_map.py`): CHM raster (max height per cell),
      per-tree 3D canopy bounding boxes, convex hulls with volumes, per-tree
      crops and masks.
- [ ] **Restructure** (`restructure.py`): move flat stage outputs into the
      canonical layout; derive plot-level instance masks (stem: ground-truth
      skeleton discs with nearest-center arbitration; canopy: highest point's
      tree_id, zeroed where the CHM is empty); rewrite the scene JSON with
      relative links.
- [ ] **Info JSONs** (`build_info.py`, always runs): per-tree / per-shrub /
      plot metadata — occupied-voxel volumes, DBH and height, generation
      params, neighborhood statistics (other trees/shrubs/grass within the
      bbox, with voxel-volume proportions).
- [ ] **Normalization** (optional DL preprocessing, `normalize.py`): max
      points-per-m3 cap, rescale to the unit cube, optional density
      voxelization (cube L^3 or LxLxZ lattices with density / num_trees /
      per-order ratio columns); categorical is `branch_code` or
      `semantic_code` depending on the source cloud; run-wide category union so
      every `_normalized.las` shares one schema; runs at plot, tree and shrub
      levels.

## 4. Deep learning (separate companion project)

- [x] Lives in the companion project `../SimpleSynthForest-DL` (see its
      README): cloud cache building, autoencoder + encoder registry,
      training, Optuna tuning, memorisation probes, blobbiness diagnostics,
      and the trial tree-count task. It consumes this project's canonical
      layout and normalized clouds; nothing here imports it.

## 5. Orchestration

- [ ] **Pipeline runner** (`run_pipeline.py`): discovers `plot_*.csv` tables,
      runs all stages per plot as isolated subprocesses with shared
      parameters, gates terrain/understory on flags, defers normalize to one
      dataset-wide call (uniform schema), runs build_info last, cleans stale
      output dirs, supports `--limit`, `--stop-on-error`, `--dry-run`, failure
      logging with meaningful exit codes.

## 6. Auxiliaries

- [ ] `generate_example_tree.py` (upstream simulator demo).
- [ ] `colorize_plot.py` (PLY visualization colorization).
- [ ] `requirements.txt` (base deps; torch/optuna are training-only extras).
- [ ] `runs/` output convention (checkpoints, history JSON, reconstructions).

## Notes for a rebuild

- The pipeline is feature-complete per stage; the pluggable generator seam,
  config files, and a plain tree-list entry point are interface work rather
  than missing functionality.
- About a third of module 3 (five cloud variants per object,
  restructure-as-a-stage) is redundant storage: a rebuild could write only the
  canonical layout directly and drop the intermediate tiers without losing
  any information.
