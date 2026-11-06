# Terrain integration (`terrain.py`)

Stage that drapes a placed plot onto a DTM (or, when no DTM is configured,
onto a **flat plane at z = 0**), adds bare-earth ground points and a
height-above-ground column, and classifies points per the LAS spec. Runs
**after `restructure.py`, before `normalize.py`** — the stage always runs;
the three topography modes are chosen per the design spec (Section 3.3):
a provided DTM wins over `--generate-dtm` (a warning describes the ignored
choice), and when neither is given the topography is None: a uniform flat
plane at z = 0, with ground sampled at the default density (warning when
`--density` was not provided):

```
python terrain.py plots/plot_11519.json --dtm sample_dtm.asc \
    --density 2 --noise 0.05
python terrain.py plots/plot_11519.json                # flat plane (None)
```

Writes:

```
plot_placement/las/<PLOT_UID>_terrain.laz
individual_trees/<PLOT_UID_TREE_UID>/las/<..>_terrain.laz
```

## The DTM

Either **provided** or **generated per plot**.

**Provided** (`--dtm file.asc`): an ESRI ASCII grid in metres, used exactly as
given — no tiling, stretching, or resampling. Its horizontal centre is aligned
to the plot centre `(0,0)`, and it **must cover the whole plot footprint** or
`terrain.py` raises a clear error naming both extents:

```
plot_11519: DTM does not cover the plot footprint. DTM 5.0x5.0 m centred at
(0.0,0.0) spans x[-2.5,2.5] y[-2.5,2.5]; plot needs x[-15.1,13.7] y[-17.3,12.3].
Provide a larger DTM or use --generate-dtm (no tiling/stretching is applied).
```

**Generated** (`--generate-dtm`): synthesize a DTM per plot with the project's
`generate_dtm` (diamond-square), so provided and generated DTMs share one
terrain model. You supply `--gen-cell-size`, `--gen-bumpiness` and
`--gen-relief`; the **bounds are auto-sized to the plot's point-cloud bounding
box** (symmetric about the centre, plus one cell of margin), so coverage is
guaranteed and no coverage check is needed. Each plot gets its own DTM sized to
its own footprint. If both `--dtm` and `--generate-dtm` are given, the
provided DTM wins and a warning describes the ignored generation.

**Flat** (neither flag): the None topography — a uniform plane at z = 0
with the same DTM interface (sampling returns 0 everywhere), so draping is
a no-op (tree bases already sit at z = 0), height_above_ground equals z,
and ground points are sampled across the plot bbox at `--density` (default
1.0 per cell, applied with a warning when unconfigured).

Heights (either source) are read by bilinear interpolation at each point's x,y.

## Steps

1. **Align** the DTM centre to the plot centre `(0,0)`.
2. **Drape each tree**: shift the whole tree vertically so its lowest point sits
   on the ground at that tree's x,y (`z += DTM(x_low,y_low) − z_low`). The tree
   is dropped rigidly — geometry is unchanged, only its base lands on the
   terrain (every tree's lowest point ends at height-above-ground ≈ 0).
3. **Sample ground**: `--density` points per DTM cell across the plot's x-y
   bounding box, `z = DTM(x,y) + noise`, `semantic_code` 6 (ground),
   `classification` 2. `--noise` is the ± half-range of uniform vertical noise
   (m). Ground HAG is just that noise (±`noise`).
4. **Classify by height above ground** (LAS spec): ≤ 0.5 m → 3 (low veg),
   0.5–2 m → 4 (med veg), > 2 m → 5 (high veg). Ground stays 2. This is written
   to the LAS `classification` field: each class's HAG falls exactly in its
   band.
5. **Height-above-ground column**: every point gets `height_above_ground =
   z − DTM(x,y)`, a float32 extra dim. For a draped tree point this equals its
   local height within the tree; for ground it's the noise.

## Per-tree terrain files

Each tree's `_terrain.laz` contains the tree **and its surrounding points** —
every point (any class, any tree, plus ground) inside that tree's canopy x-y
bounding box (the same bbox used elsewhere). So an understorey or edge tree's
file includes overlapping neighbour canopy and the ground patch beneath it.

## Two orthogonal labels

`--source` defaults to the furthest-processed cloud present; `trimmed`/`raw`
also work. The terrain cloud uses the canonical schema (README_schema), and its
two class axes stay **orthogonal**:

- **`semantic_code`** — object identity, preserved unchanged from the source for
  tree points (1 trunk, 2 branch, 3 foliage) and set to **6 (ground)** for the
  new ground points.
- **`classification`** — the LAS height-above-ground band (2 ground, 3/4/5 low/
  med/high veg), independent of identity: a foliage point (semantic 3) is class
  3, 4, or 5 depending on how high it sits. Both the plot and per-tree terrain
  files carry the full canonical schema.

## In the pipeline

The terrain stage always runs (per plot, after restructure); `--dtm` /
`--generate-dtm` choose the DTM source, and neither means the flat plane:

```
--dtm PATH             provided DTM .asc to drape onto (wins over --generate-dtm)
--generate-dtm         synthesize a per-plot DTM instead (bounds auto-fit)
--gen-cell-size N      generated/flat DTM cell size, m (default 0.25)
--gen-bumpiness N      generated DTM roughness (default 0.5)
--gen-relief M         generated DTM vertical relief, m (default 10)
--ground-density N     ground points per DTM cell (default 1)
--ground-noise M       vertical noise half-range, m (default 0)
--terrain-source S     raw|trimmed (default trimmed)
```

Terrain slots between restructure and normalize: it is skipped when
`--dtm`/`--generate-dtm` is absent, and errors cleanly (per plot, so one
uncovered plot is skipped rather than sinking the batch) when the DTM is too
small.
