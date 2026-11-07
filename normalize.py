"""
normalize.py: optional deep-learning normalizations of the placed clouds.

Users pick any combination of three steps, always applied in this order:

  1. max_density   hard cap of N points per cubic metre (random subsample,
                   all semantics preserved on survivors)
  2. spatial       per-axis rescale to the unit cube [-1, 1]^3 (aspect ratio
                   NOT preserved)
  3. density       replace the cloud with a fixed lattice of points, one per
                   cell/voxel centre (empties included), carrying:
                     - density        the pseudocode value
                     - num_trees      distinct tree_id count in that cell/voxel
                     - branch_order_k one column per branch order in the RUN,
                                      ratio of the cell/voxel's points with
                                      that order (padded with 0 where absent)

density has two mutually exclusive modes, shared as the base shape across the
plot and tree levels; L (and Z) may differ per level:

  square (L, Z):  L x L columns on x-y, each split into Z height bins.
                  L*L*Z lattice points.
                  density = S_bin / S_column   (profile-shape descriptor;
                            each column's Z values sum to 1)
  cube (L):       L x L x L voxels. L^3 lattice points.
                  density = S_voxel / S_total  (all voxels sum to 1)

  num_trees and the branch_order_k ratios are always computed per cell/voxel
  (that cell's own points), independent of the density denominator.

Density-normalized clouds are written to the layout's *_normalized.las slots.
The branch_order column set is the union over every cloud written in the run,
so all normalized files share one schema (absent orders are all-zero columns).

Run after restructure.py:

    python normalize.py plot_11519.csv --max-density 5000 --spatial \\
        --density cube --plot-l 64 --tree-l 32

Requires: numpy, pandas, laspy (via las.py), layout.py
"""

from __future__ import annotations

import argparse
import json
import os

import numpy as np

import layout
from las import read_points_las, export_points_las, export_canonical_las, read_canonical_las
import schema

# source columns after read_points_las(..., extra_fields=SRC_EXTRA)
SRC_EXTRA = ["branch_id", "branch_order"]
C_X, C_Y, C_Z, C_SEM, C_TID, C_TRIM, C_BID, C_ORDER = range(8)

# The per-voxel categorical ratio columns adapt to the source cloud. Every
# canonical cloud carries both branch_code and semantic_code, so the choice is
# purely about which is meaningful for that stage:
#   raw / trimmed       -> branch_code   (branch structure: 0,1,2,3,...)
#   terrain / understory -> semantic_code (object identity incl. ground/shrub/
#                          grass, which branch_code doesn't distinguish)
# value: (canonical field name, column prefix)
CATEGORICAL = {
    "branch_code": ("branch_code", "branch_code"),
    "semantic_code": ("semantic_code", "semantic"),
}


def _categorical_kind(source):
    """Which categorical the given source cloud carries."""
    return "semantic_code" if source in ("terrain", "understory") else "branch_code"


# ---------------------------------------------------------------------------
# step 1: max density
# ---------------------------------------------------------------------------

def cap_density(pts, n_per_m3, rng):
    """Randomly drop points so no 1 m^3 voxel holds more than n_per_m3.

    Fixed 1 m grid keyed off floor(coord). Survivors keep every column.
    """
    if n_per_m3 is None or len(pts) == 0:
        return pts
    keys = np.floor(pts[:, :3]).astype(np.int64)
    # pack the 3 integer voxel coords into one key
    k = np.ascontiguousarray(keys).view(
        np.dtype((np.void, keys.dtype.itemsize * 3))).ravel()
    keep = np.ones(len(pts), dtype=bool)
    order = np.argsort(k, kind="stable")
    ks = k[order]
    # boundaries of runs of equal voxel key
    bounds = np.flatnonzero(np.r_[True, ks[1:] != ks[:-1], True])
    for a, b in zip(bounds[:-1], bounds[1:]):
        if b - a > n_per_m3:
            grp = order[a:b]
            drop = rng.choice(grp, size=(b - a) - n_per_m3, replace=False)
            keep[drop] = False
    return pts[keep]


# ---------------------------------------------------------------------------
# step 2: spatial normalization to the unit cube
# ---------------------------------------------------------------------------

def to_unit_cube(pts):
    """Per-axis rescale x,y,z into [-1, 1] (aspect ratio not preserved)."""
    if len(pts) == 0:
        return pts
    out = pts.copy()
    lo = pts[:, :3].min(axis=0)
    hi = pts[:, :3].max(axis=0)
    span = np.where(hi > lo, hi - lo, 1.0)
    out[:, :3] = 2.0 * (pts[:, :3] - lo) / span - 1.0
    return out


# ---------------------------------------------------------------------------
# step 3: density normalization
# ---------------------------------------------------------------------------

def _extent(pts, unit_cube):
    """(lo, hi) 3-vectors of the grid extent."""
    if unit_cube:
        return np.array([-1.0, -1.0, -1.0]), np.array([1.0, 1.0, 1.0])
    lo = pts[:, :3].min(axis=0)
    hi = pts[:, :3].max(axis=0)
    return lo, hi


def _bin_index(vals, lo, hi, n):
    """Bin index in [0, n-1] for vals over [lo, hi]."""
    span = hi - lo if hi > lo else 1.0
    idx = np.floor((vals - lo) / span * n).astype(int)
    return np.clip(idx, 0, n - 1)


def _order_ratio_columns(order_of_pts, cell_ids, n_cells, run_orders):
    """(n_cells, len(run_orders)) ratio of each cell's points per branch order."""
    ncol = len(run_orders)
    counts = np.zeros((n_cells, ncol), dtype=np.float64)
    order_col = {o: j for j, o in enumerate(run_orders)}
    if len(cell_ids):
        for j, o in enumerate(run_orders):
            m = order_of_pts == o
            if m.any():
                counts[:, j] = np.bincount(cell_ids[m], minlength=n_cells)
    totals = counts.sum(axis=1, keepdims=True)
    ratios = np.divide(counts, totals, out=np.zeros_like(counts),
                       where=totals > 0)
    return ratios


def _num_trees_per_cell(tid_of_pts, cell_ids, n_cells):
    """Distinct tree_id count per cell."""
    out = np.zeros(n_cells, dtype=np.int64)
    if not len(cell_ids):
        return out
    # unique (cell, tree) pairs, then count per cell
    pair = np.stack([cell_ids, tid_of_pts.astype(np.int64)], axis=1)
    pair = np.unique(pair, axis=0)
    np.add.at(out, pair[:, 0], 1)
    return out


def density_cube(pts, L, unit_cube, cat_values, categories, tid_col=4):
    """L^3 voxel lattice. density = S_voxel / S_total.

    cat_values: per-point categorical value (branch order or object type).
    categories: run-wide sorted category list -> one ratio column each.
    """
    lo, hi = _extent(pts, unit_cube)
    n_cells = L * L * L
    centers = _cube_centers(lo, hi, L)

    density = np.zeros(n_cells, dtype=np.float64)
    num_trees = np.zeros(n_cells, dtype=np.int64)
    ratios = np.zeros((n_cells, len(categories)), dtype=np.float64)

    if len(pts):
        ix = _bin_index(pts[:, C_X], lo[0], hi[0], L)
        iy = _bin_index(pts[:, C_Y], lo[1], hi[1], L)
        iz = _bin_index(pts[:, C_Z], lo[2], hi[2], L)
        cell = (ix * L + iy) * L + iz
        s_total = len(pts)
        density = np.bincount(cell, minlength=n_cells) / float(s_total)
        num_trees = _num_trees_per_cell(pts[:, tid_col], cell, n_cells)
        ratios = _order_ratio_columns(cat_values, cell, n_cells, categories)

    return centers, density, num_trees, ratios


def density_square(pts, L, Z, unit_cube, cat_values, categories, tid_col=4):
    """L x L x Z cell-bin lattice. density = S_bin / S_column."""
    lo, hi = _extent(pts, unit_cube)
    n_cells = L * L * Z
    centers = _square_centers(lo, hi, L, Z)

    density = np.zeros(n_cells, dtype=np.float64)
    num_trees = np.zeros(n_cells, dtype=np.int64)
    ratios = np.zeros((n_cells, len(categories)), dtype=np.float64)

    if len(pts):
        ix = _bin_index(pts[:, C_X], lo[0], hi[0], L)
        iy = _bin_index(pts[:, C_Y], lo[1], hi[1], L)
        iz = _bin_index(pts[:, C_Z], lo[2], hi[2], Z)
        column = ix * L + iy                       # 0 .. L*L-1
        cell = column * Z + iz                     # 0 .. L*L*Z-1
        s_column = np.bincount(column, minlength=L * L).astype(np.float64)
        s_bin = np.bincount(cell, minlength=n_cells).astype(np.float64)
        col_of_cell = np.arange(n_cells) // Z
        denom = s_column[col_of_cell]
        density = np.divide(s_bin, denom, out=np.zeros_like(s_bin),
                            where=denom > 0)
        num_trees = _num_trees_per_cell(pts[:, tid_col], cell, n_cells)
        ratios = _order_ratio_columns(cat_values, cell, n_cells, categories)

    return centers, density, num_trees, ratios


def _cube_centers(lo, hi, L):
    cx = _centers_1d(lo[0], hi[0], L)
    cy = _centers_1d(lo[1], hi[1], L)
    cz = _centers_1d(lo[2], hi[2], L)
    gx, gy, gz = np.meshgrid(cx, cy, cz, indexing="ij")
    return np.stack([gx.ravel(), gy.ravel(), gz.ravel()], axis=1)


def _square_centers(lo, hi, L, Z):
    cx = _centers_1d(lo[0], hi[0], L)
    cy = _centers_1d(lo[1], hi[1], L)
    cz = _centers_1d(lo[2], hi[2], Z)
    gx, gy, gz = np.meshgrid(cx, cy, cz, indexing="ij")
    return np.stack([gx.ravel(), gy.ravel(), gz.ravel()], axis=1)


def _centers_1d(lo, hi, n):
    span = hi - lo if hi > lo else 1.0
    edges = lo + (np.arange(n) + 0.5) * span / n
    return edges


# ---------------------------------------------------------------------------
# per-cloud driver
# ---------------------------------------------------------------------------

def normalize_cloud(src_las, dst_las, *, max_density=None, spatial=False,
                    density=None, L=None, Z=None,
                    run_categories, categorical="branch_code", seed=0):
    """Read a canonical source cloud, apply the selected steps, write dst_las.

    categorical: which per-point class drives the ratio columns -
      'branch_code' (raw/trimmed) or 'semantic_code' (terrain/understory).
      run_categories is the run-wide union for that categorical.

    Returns a small dict describing the written file (schema + point count).
    """
    field, prefix = CATEGORICAL[categorical]
    level = ("tree" if ("individual_trees" in src_las
                        or "individual_shrubs" in src_las) else "plot")
    xyz, cols = read_canonical_las(src_las, level=level)
    cat_full = cols[field].astype(int)
    rng = np.random.default_rng(seed)

    # carry EVERY canonical dim through, aligned to the points, so nothing is
    # dropped by normalization (shrub_id, classification, branch_code, branch_id,
    # height_above_ground were previously lost). Columns after xyz:
    #   3..(3+K-1) canonical dims (schema order), then the categorical at the end.
    fields = schema.PLOT_FIELDS
    extra_stack = [cols[f].astype(float) for f in fields]
    pts_aug = np.column_stack([xyz, *extra_stack, cat_full.astype(float)])
    pts_aug = cap_density(pts_aug, max_density, rng)
    if spatial:
        pts_aug = to_unit_cube(pts_aug)
    cat_col = 3 + len(fields)
    cat_values = pts_aug[:, cat_col].astype(int)
    pts = pts_aug

    cat_cols = [f"{prefix}_{c}" for c in run_categories]

    if density is None:
        os.makedirs(os.path.dirname(dst_las), exist_ok=True)
        canon = {f: pts[:, 3 + i].astype(schema.dtype_of(f))
                 for i, f in enumerate(fields)}
        export_canonical_las(dst_las, pts[:, :3], canon, level=level)
        return {"mode": "points", "n_points": int(len(pts)),
                "columns": ["x", "y", "z"] + (schema.TREE_FIELDS
                                              if level == "tree"
                                              else schema.PLOT_FIELDS)}

    if density == "cube":
        centers, dens, ntree, ratios = density_cube(
            pts, L, spatial, cat_values, run_categories,
            tid_col=3 + schema.PLOT_FIELDS.index("tree_id"))
    elif density == "square":
        centers, dens, ntree, ratios = density_square(
            pts, L, Z, spatial, cat_values, run_categories,
            tid_col=3 + schema.PLOT_FIELDS.index("tree_id"))
    else:
        raise ValueError(f"unknown density mode {density!r}")

    lattice = np.column_stack([centers, np.zeros((len(centers), 2))])
    extra = {"density": dens.astype(np.float32),
             "num_trees": ntree.astype(np.uint16)}
    for j, name in enumerate(cat_cols):
        extra[name] = ratios[:, j].astype(np.float32)

    os.makedirs(os.path.dirname(dst_las), exist_ok=True)
    export_points_las(dst_las, lattice, extra_fields=extra)
    return {"mode": density, "L": L, "Z": (Z if density == "square" else None),
            "n_points": int(len(centers)), "categorical": categorical,
            "columns": ["x", "y", "z", "density", "num_trees"] + cat_cols}


# ---------------------------------------------------------------------------
# run over a plot
# ---------------------------------------------------------------------------

def _scan_run_categories(plot_dirs, plot_uids, source):
    """Union of the source's categorical values across every plot in the run.

    branch_code for raw/trimmed; semantic_code for terrain/understory.
    """
    kind = _categorical_kind(source)
    field, _ = CATEGORICAL[kind]
    vals = set()
    for pd_, pu in zip(plot_dirs, plot_uids):
        src = layout.plot_path(pd_, pu, f"{source}_las", create=False)
        if os.path.exists(src):
            _, cols = read_canonical_las(src, level="plot")
            vals.update(np.unique(cols[field].astype(int)).tolist())
    return sorted(vals) if vals else [0]


# backward-compatible alias
def _scan_run_orders(plot_dirs, plot_uids, source):
    return _scan_run_categories(plot_dirs, plot_uids, source)


def _resolve_source(plot_dir, plot_uid, source):
    """Return the source that actually exists on disk.

    Precedence walks from the requested source toward the furthest-processed
    cloud available: understory -> terrain -> trimmed -> raw. So with the
    default the last cloud in the chain (understory if the step ran, else
    terrain, else trimmed) is picked up automatically. branch_order is only on
    trimmed/raw; terrain/understory carry the object-type categorical instead.
    Returns (effective_source, plot_src_path) or (None, None).
    """
    chain = ["understory", "terrain", "trimmed", "raw"]
    # start at the requested source, then fall through the rest of the chain
    if source in chain:
        order = chain[chain.index(source):] + [s for s in chain
                                               if chain.index(s) < chain.index(source)]
    else:
        order = [source] + chain
    seen = []
    for s in order:
        if s in seen:
            continue
        seen.append(s)
        p = layout.plot_path(plot_dir, plot_uid, f"{s}_las", create=False)
        if os.path.exists(p):
            return s, p
    return None, None


def normalize_plot(csv_path, plot_dir=None, *, max_density=None, spatial=False,
                   density=None, plot_l=None, plot_z=None, tree_l=None,
                   tree_z=None, source="understory",
                   run_orders=None, seed=0, verbose=True):
    plot_uid = os.path.splitext(os.path.basename(csv_path))[0]
    if plot_dir is None:
        plot_dir = os.path.join(os.path.dirname(os.path.abspath(csv_path)), plot_uid)

    effective_source, plot_src = _resolve_source(plot_dir, plot_uid, source)
    if effective_source is None:
        expected = layout.plot_path(plot_dir, plot_uid, "trimmed_las", create=False)
        flat = layout.cloud_file(os.path.join(plot_dir, layout.PLOT_PLACEMENT),
                                 f"{plot_uid}_trimmed")
        hint = ("The flat cloud exists but has not been moved into the canonical "
                "layout yet -- run restructure.py first." if os.path.exists(flat)
                else "Run the pipeline (generate -> assemble -> stem -> canopy -> "
                     "restructure) for this plot before normalize.")
        raise FileNotFoundError(
            f"{plot_uid}: no source cloud at {expected} (or the raw equivalent). "
            f"normalize.py reads from the canonical layout, so restructure.py "
            f"must run first. {hint}")
    if effective_source != source and verbose:
        print(f"{plot_uid}: '{source}' cloud not found, using "
              f"'{effective_source}' instead"
              + (" (branch_order will be all-zero on raw clouds)"
                 if effective_source == "raw" else ""))
    source = effective_source

    if run_orders is None:
        run_orders = _scan_run_categories([plot_dir], [plot_uid], source)

    categorical = _categorical_kind(source)
    prefix = CATEGORICAL[categorical][1]

    scene_path = layout.plot_path(plot_dir, plot_uid, "scene", create=False)
    scene = json.load(open(scene_path)) if os.path.exists(scene_path) else {}
    tree_uids = list(scene.get("trees", {}).keys()) or _tree_uids_from_disk(plot_dir)

    written = []

    # plot level
    src = plot_src
    if os.path.exists(src):
        info = normalize_cloud(
            src, layout.plot_path(plot_dir, plot_uid, "normalized_las"),
            max_density=max_density, spatial=spatial, density=density,
            L=plot_l, Z=plot_z,
            run_categories=run_orders, categorical=categorical, seed=seed)
        written.append(("plot", plot_uid, info))

    # tree level
    for uid in tree_uids:
        src = layout.tree_path(plot_dir, plot_uid, uid, f"{source}_las", create=False)
        if not os.path.exists(src):
            continue
        info = normalize_cloud(
            src, layout.tree_path(plot_dir, plot_uid, uid, "normalized_las"),
            max_density=max_density, spatial=spatial, density=density,
            L=tree_l, Z=tree_z,
            run_categories=run_orders, categorical=categorical, seed=seed)
        written.append(("tree", uid, info))

    # shrub level -- normalize each individual shrub the same way as trees.
    # Shrubs only have raw / understory clouds, so map the requested source to
    # the nearest shrub cloud: understory-source -> understory, else raw.
    shrub_source = "understory_las" if source == "understory" else "raw_las"
    for uid in _shrub_uids_from_disk(plot_dir):
        src = layout.shrub_path(plot_dir, plot_uid, uid, shrub_source, create=False)
        if not os.path.exists(src):
            continue
        info = normalize_cloud(
            src, layout.shrub_path(plot_dir, plot_uid, uid, "normalized_las"),
            max_density=max_density, spatial=spatial, density=density,
            L=tree_l, Z=tree_z,
            run_categories=run_orders, categorical=categorical, seed=seed)
        written.append(("shrub", uid, info))

    if not written:
        raise RuntimeError(
            f"{plot_uid}: found the source cloud but produced 0 outputs. "
            f"This usually means the tree list is empty; checked scene 'trees' "
            f"and {os.path.join(plot_dir, layout.INDIVIDUAL_TREES)}.")

    # record schema in the scene
    scene["normalization"] = {
        "steps": {"max_density": max_density, "spatial": bool(spatial),
                  "density": density},
        "source": source,
        "categorical": categorical,
        "run_categories": run_orders,
        "plot": {"L": plot_l, "Z": plot_z},
        "tree": {"L": tree_l, "Z": tree_z},
        "columns": (["x", "y", "z", "density", "num_trees"]
                    + [f"{prefix}_{o}" for o in run_orders])
                   if density else ["x", "y", "z", "semantic", "tree_id"],
    }
    if os.path.exists(scene_path):
        with open(scene_path, "w") as f:
            json.dump(scene, f, indent=2)

    if verbose:
        n = len(written)
        pl = next((w[2]["n_points"] for w in written if w[0] == "plot"), 0)
        print(f"{plot_uid}: normalized {n} clouds "
              f"(mode={density or 'points'}, {categorical}={run_orders}); "
              f"plot cloud {pl} points")
    return written


def _tree_uids_from_disk(plot_dir):
    d = os.path.join(plot_dir, layout.INDIVIDUAL_TREES)
    if not os.path.isdir(d):
        return []
    out = []
    for name in sorted(os.listdir(d)):
        if os.path.isdir(os.path.join(d, name)) and "_" in name:
            out.append(name.rsplit("_", 1)[-1])
    return out


def _shrub_uids_from_disk(plot_dir):
    """Shrub uids from individual_shrubs/<plot>_shrub_<uid> dir names."""
    d = os.path.join(plot_dir, layout.INDIVIDUAL_SHRUBS)
    if not os.path.isdir(d):
        return []
    out = []
    for name in sorted(os.listdir(d)):
        if os.path.isdir(os.path.join(d, name)) and "_shrub_" in name:
            out.append(name.rsplit("_", 1)[-1])
    return out


def main():
    ap = argparse.ArgumentParser(description="Normalize placed clouds for DL")
    ap.add_argument("csv", nargs="+", help="plot table csv(s)")
    ap.add_argument("-d", "--plot-dir", default=None)
    ap.add_argument("--source", default="understory",
                    choices=["raw", "trimmed", "terrain", "understory"],
                    help="which cloud to normalize; default resolves to the "
                         "furthest-processed cloud present (understory -> "
                         "terrain -> trimmed -> raw). terrain/understory use "
                         "the object-type categorical; trimmed/raw use branch_order")
    ap.add_argument("--max-density", type=int, default=None,
                    help="cap of N points per cubic metre")
    ap.add_argument("--spatial", action="store_true",
                    help="rescale to the unit cube [-1,1]^3 (per-axis)")
    ap.add_argument("--density", choices=["square", "cube"], default=None)
    ap.add_argument("--plot-l", type=int, default=None)
    ap.add_argument("--plot-z", type=int, default=None, help="square only")
    ap.add_argument("--tree-l", type=int, default=None)
    ap.add_argument("--tree-z", type=int, default=None, help="square only")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    if args.density == "square" and (args.plot_z is None or args.tree_z is None):
        ap.error("square density requires --plot-z and --tree-z")
    if args.density and (args.plot_l is None or args.tree_l is None):
        ap.error("density requires --plot-l and --tree-l")

    plot_dirs = [args.plot_dir or os.path.join(
        os.path.dirname(os.path.abspath(c)),
        os.path.splitext(os.path.basename(c))[0]) for c in args.csv]
    plot_uids = [os.path.splitext(os.path.basename(c))[0] for c in args.csv]

    # resolve each plot's effective source; they must share a categorical kind
    eff = [_resolve_source(pd_, pu, args.source)[0]
           for pd_, pu in zip(plot_dirs, plot_uids)]
    eff = [e for e in eff if e is not None]
    kinds = {_categorical_kind(e) for e in eff}
    if len(kinds) > 1:
        ap.error(f"plots resolve to mixed source kinds {sorted(set(eff))}; "
                 f"the categorical schema (branch_order vs object-type) would "
                 f"differ. Run the same stages for all plots, or pass --source "
                 f"explicitly to force one.")
    scan_source = eff[0] if eff else args.source

    orders = _scan_run_categories(plot_dirs, plot_uids, scan_source)

    for path in args.csv:
        normalize_plot(path, plot_dir=args.plot_dir,
                       max_density=args.max_density, spatial=args.spatial,
                       density=args.density, plot_l=args.plot_l,
                       plot_z=args.plot_z, tree_l=args.tree_l, tree_z=args.tree_z,
                       source=args.source,
                       run_orders=orders, seed=args.seed)


if __name__ == "__main__":
    main()
