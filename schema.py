"""
schema.py: the canonical point-cloud column schema for the whole pipeline.

One definition of the extra dimensions every cloud carries, so plot-level and
tree-level clouds stay consistent from generation through normalization. A stage
that hasn't computed a given dim yet writes its zero/default fill (e.g. a raw
cloud has shrub_id = 0 everywhere); the column still exists so the schema never
changes shape across the chain.

Point rows are always x, y, z, plus these extra dims. Two schemas:
  PLOT  : tree_id, shrub_id, semantic_code, classification, branch_code, hag
  TREE  : the above + branch_id

Semantic code (int8):
  0 ground (also "no class needed" — ground/noise are the traditional LAS
  classification's job, so there are no separate unassigned/noise codes),
  1 tree(trunk stem), 2 branch(branch stem), 3 foliage, 8 shrub, 9 grass
LAS classification (uint8): 2 ground, 3 low veg, 4 med veg, 5 high veg, 7 noise
Branch code (int8): 0 not-tree, 1 trunk, 2 connected-to-trunk, 3 connected-to-2,
  ... (== the old branch_order + 1 on tree points)

tree_id / shrub_id / branch_id are all int32: each is a per-plot sequential
count, so the id space must outgrow whatever a plot can actually produce.
int16 (max 32767) was the original choice but is easy to exceed once plot
area grows -- FIA's placeholder TPA expansion scales the expected tree count
with plot area, so a 100 m plot can legitimately need tree ids past 32767
where the default ~7.3 m FIA subplot never would. branch_id is compacted to a
small plot-unique sequential id at assembly time (see assemble_plot.py's
`_branch_id_remap`), but that counter is cumulative *across every tree in the
plot* -- a plot with many large, many-branched trees can still assign more
than 32767 distinct branch ids in one plot, so it needs the same headroom as
tree_id/shrub_id, not less.
"""

from __future__ import annotations

import numpy as np

# ---- semantic codes -------------------------------------------------------
# 0 doubles as "no semantic class needed": ground and noise are the
# traditional LAS classification's job (Section 8.3), so there are no
# separate unassigned or noise semantic codes.
SEM_GROUND = 0
SEM_TREE = 1        # trunk stem
SEM_BRANCH = 2      # branch stem
SEM_FOLIAGE = 3
SEM_SHRUB = 8
SEM_GRASS = 9

# ---- LAS classification codes --------------------------------------------
LAS_GROUND = 2
LAS_LOW_VEG = 3
LAS_MED_VEG = 4
LAS_HIGH_VEG = 5
LAS_NOISE = 7

# ---- branch code ----------------------------------------------------------
BRANCH_NONE = 0     # not a tree point

# ---- the canonical extra-dim schema --------------------------------------
# name -> (numpy dtype, default fill). Order here is the on-disk column order.
_PLOT_DIMS = [
    ("tree_id",             np.int32,   0),
    ("shrub_id",            np.int32,   0),
    ("semantic_code",       np.int8,    SEM_GROUND),
    ("classification",      np.uint8,   0),
    ("branch_code",         np.int8,    BRANCH_NONE),
    ("branch_id",           np.int32,   0),
    ("height_above_ground", np.float32, 0.0),
]
_TREE_EXTRA = []

PLOT_DIMS = _PLOT_DIMS
TREE_DIMS = _PLOT_DIMS      # plot and tree schemas are now identical

# convenience name lists (on-disk order)
PLOT_FIELDS = [n for n, _, _ in PLOT_DIMS]
TREE_FIELDS = [n for n, _, _ in TREE_DIMS]

# dtype / default lookup (union of both schemas)
_DT = {n: dt for n, dt, _ in PLOT_DIMS + _TREE_EXTRA}
_DEF = {n: d for n, _, d in PLOT_DIMS + _TREE_EXTRA}


def dtype_of(name):
    return _DT[name]


def default_of(name):
    return _DEF[name]


def empty_extra(fields, n):
    """{name: zero/default array of length n} for the given field names."""
    return {name: np.full(n, _DEF[name], dtype=_DT[name]) for name in fields}


def build_extra(fields, n, **values):
    """Canonical extra-dim dict: start from defaults, overlay provided arrays.

    Each provided array is cast to the field's canonical dtype. Unspecified
    fields are filled with their default (so a stage only supplies what it
    knows; the rest are zeros of the right dtype).
    """
    out = empty_extra(fields, n)
    for name, arr in values.items():
        if name not in out:
            raise KeyError(f"{name!r} is not in this schema ({fields})")
        a = np.asarray(arr)
        if len(a) != n:
            raise ValueError(f"{name}: length {len(a)} != {n}")
        out[name] = a.astype(_DT[name])
    return out
