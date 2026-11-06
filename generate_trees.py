"""
generate_trees.py: generate one SimpleSynthTree point cloud per row of a
PlausibiliTree plot table.

    python generate_trees.py plot_11519.json --seed 42   # canonical tree list
    python generate_trees.py plot_11519.csv  --seed 42   # legacy SST table

For each successfully placed tree this writes, under
<outdir>/individual_trees/:

    <TREE_UID>_raw.las            labelled point cloud, plot coordinates
    <TREE_UID>_generation.json    every parameter needed to reproduce it
    <TREE_UID>_raw_skeleton.json  trunk polylines + branch axes

plus a manifest.csv summarising the run.

Point budget is normalised by tree height AND DBH so that point density is
constant across the plot and thick trees carry proportionally more points than
thin ones: Np = points_per_metre * ceil(HT_HI) * BASE_DIA_HI.

Requires: numpy, pandas, trimesh, laspy (the latter two via treesim / las).
"""

from __future__ import annotations

import argparse
import datetime as _dt
import json
import math
import os
import zlib

import numpy as np
import pandas as pd

from treesim import gen_simtree
from las import export_canonical_las
import layout
import tree_list

POINTS_PER_METRE = 1024

# treesim offsets fork tops by up to 2.0*rand + 1.5 metres from the main stem
MAX_FORK_OFFSET = 3.5


def _max_stem_excursion(skeleton: dict) -> float:
    """Furthest any trunk vertex strays horizontally from the tree base."""
    bx, by = skeleton["base_xy"]
    worst = 0.0
    for poly in skeleton["trunk"]["polylines"]:
        arr = np.asarray(poly, dtype=float)
        worst = max(worst, float(np.hypot(arr[:, 0] - bx, arr[:, 1] - by).max()))
    return worst

# plot-table column -> treesim model_params key
_RANGE_PARAMS = {
    "height_range": ("HT_LO", "HT_HI"),
    "diam_range": ("BASE_DIA_LO", "BASE_DIA_HI"),
    "split_height_range": ("SPLIT_HT_LO", "SPLIT_HT_HI"),
    "num_branches": ("NUM_BRANCH_LO", "NUM_BRANCH_HI"),
    "min_can_height": ("MIN_CAN_HT_LO", "MIN_CAN_HT_HI"),
    "max_can_width_height": ("MAX_CAN_W_HT_LO", "MAX_CAN_W_HT_HI"),
}
_SCALAR_PARAMS = {
    "split_prob": "SPLIT_PROB",
    "max_can_width": "MAX_CAN_W",
    "tree_top_dist": "TOP_DIST",
    "tree_mid_dist": "MID_DIST",
    "foliage_noise": "FP_RAND_HT",
}


# ---------------------------------------------------------------------------
# Row -> simulator parameters
# ---------------------------------------------------------------------------

def params_from_row(row: pd.Series) -> dict:
    """Translate one plot-table row into a treesim model_params dict."""
    p = {}
    for key, (lo, hi) in _RANGE_PARAMS.items():
        p[key] = [float(row[lo]), float(row[hi])]
    for key, col in _SCALAR_PARAMS.items():
        p[key] = float(row[col])
    # treesim expects integer branch counts
    p["num_branches"] = [int(round(v)) for v in p["num_branches"]]
    # cascading trunk splits for forked trees (optional column pair; the
    # simulator's default [1, 1] reproduces the original single fork)
    if "SPLIT_NUM_LO" in row.index and "SPLIT_NUM_HI" in row.index:
        p["num_splits"] = [int(round(float(row["SPLIT_NUM_LO"]))),
                           int(round(float(row["SPLIT_NUM_HI"])))]
    return p


def stem_frac_from_row(row: pd.Series):
    """Stem point ratio (lambda); a list is sampled uniformly by gen_simtree."""
    lo, hi = float(row["SPR_LO"]), float(row["SPR_HI"])
    return [lo, hi] if hi > lo else lo


# ---------------------------------------------------------------------------
# Minimal JSON tree list (design spec Section 8.1) -> generation rows
# ---------------------------------------------------------------------------

DEFAULT_SIZE_TOL = 0.01        # +/- 1% around the tree list's DBH/HT


def _fg_param_series(fg_db, func_group: pd.Series, param: str, default=None):
    """One functional-group parameter as a Series, DEFAULT-fallback."""
    if param in fg_db.columns:
        lut = fg_db[param]
        out = func_group.map(lut)
        return out.fillna(lut.loc["DEFAULT"]).astype(float)
    if default is None:
        raise KeyError(f"functional-group table has no {param} column")
    return pd.Series(float(default), index=func_group.index)


def resolve_func_group(codes: pd.Series, fg_db, spcd_map: dict | None = None
                       ) -> pd.Series:
    """TREE_CODE -> functional group (spec Section 8.1 resolution order:
    species code -> functional group -> DEFAULT)."""
    fg = pd.Series("DEFAULT", index=codes.index, dtype=object)
    known = set(fg_db.index)
    codes_s = codes.astype(str)
    # direct functional-group label
    fg[codes_s.isin(known)] = codes_s[codes_s.isin(known)]
    # species code -> group lookup (codes may be numeric strings)
    if spcd_map:
        str_map = {str(k): v for k, v in spcd_map.items()}
        mapped = codes_s.map(str_map)
        hit = mapped.notna() & mapped.isin(known)
        fg[hit] = mapped[hit]
    # anything unknown falls back to DEFAULT
    fg = fg.where(fg.isin(known), "DEFAULT")
    return fg


def rows_from_tree_list(plot: dict, fg_db=None, spcd_map: dict | None = None,
                        dbh_tol: float | None = None,
                        ht_tol: float | None = None) -> pd.DataFrame:
    """Expand a minimal JSON tree list into generator rows.

    For every tree, the code's parameter dictionary entry supplies the form:
    crown width from the CW allometry, base diameter from the TAPER_B
    extrapolation, shape parameters copied from the functional group. Sizes
    (DBH, HT) come from the tree list; the tolerance bounds (DBH_TOL,
    HT_TOL columns if present, DEFAULT_SIZE_TOL otherwise) put the SST
    ranges around the list values [section 7.3]. The output rows are
    column-compatible with the legacy SST parameter CSV, so the generation
    loop below is shared by both contracts.
    """
    import fia_helper

    trees = plot["trees"].copy()
    if "TREE_UID" in trees.columns:
        trees["TREE_UID"] = trees["TREE_UID"].astype(str)
    if not len(trees):
        return pd.DataFrame()

    if fg_db is None:
        fg_db = fia_helper.load_func_group_db()
    if spcd_map is None:
        spcd_map = fia_helper.load_spcd_map()

    trees["FUNC_GROUP"] = resolve_func_group(trees["TREE_CODE"], fg_db,
                                             spcd_map)
    fg = trees["FUNC_GROUP"].astype(str)
    if ht_tol is None:
        ht_tol_series = _fg_param_series(fg_db, fg, "HT_TOL",
                                         default=DEFAULT_SIZE_TOL)
    else:
        ht_tol_series = pd.Series(float(ht_tol), index=trees.index)
    if dbh_tol is None:
        dbh_tol_series = _fg_param_series(fg_db, fg, "DBH_TOL",
                                         default=DEFAULT_SIZE_TOL)
    else:
        dbh_tol_series = pd.Series(float(dbh_tol), index=trees.index)

    ht = pd.to_numeric(trees["HT"], errors="coerce")
    # taper crown/base allometry reuse: fia_helper computes base diameter
    # and crown width from DBH/HT per functional group
    base_dia = fia_helper.dbh_to_base_diameter(trees["DBH"], ht, fg, fg_db)
    crown_w = fia_helper.crown_width(trees["DBH"], fg, fg_db)

    rows = pd.DataFrame(index=trees.index)
    rows["TREE_UID"] = trees["TREE_UID"]
    rows["FUNC_GROUP"] = fg
    rows["X"] = pd.to_numeric(trees["X"], errors="coerce").round(4)
    rows["Y"] = pd.to_numeric(trees["Y"], errors="coerce").round(4)
    rows["HT_LO"] = (ht * (1 - ht_tol_series)).round(4)
    rows["HT_HI"] = (ht * (1 + ht_tol_series)).round(4)
    rows["BASE_DIA_LO"] = (base_dia * (1 - dbh_tol_series)).round(4)
    rows["BASE_DIA_HI"] = (base_dia * (1 + dbh_tol_series)).round(4)
    rows["MAX_CAN_W"] = crown_w.round(4)

    for p in ["SPLIT_HT_LO", "SPLIT_HT_HI", "SPLIT_PROB",
              "NUM_BRANCH_LO", "NUM_BRANCH_HI",
              "MIN_CAN_HT_LO", "MIN_CAN_HT_HI",
              "MAX_CAN_W_HT_LO", "MAX_CAN_W_HT_HI",
              "TOP_DIST", "MID_DIST", "FP_RAND_HT", "SPR_LO", "SPR_HI"]:
        rows[p] = _fg_param_series(fg_db, fg, p)

    if "SPLIT_NUM_LO" in fg_db.columns and "SPLIT_NUM_HI" in fg_db.columns:
        rows["SPLIT_NUM_LO"] = _fg_param_series(
            fg_db, fg, "SPLIT_NUM_LO").astype(int)
        rows["SPLIT_NUM_HI"] = _fg_param_series(
            fg_db, fg, "SPLIT_NUM_HI").astype(int)
    else:
        rows["SPLIT_NUM_LO"] = 1
        rows["SPLIT_NUM_HI"] = 1

    rows["NUM_BRANCH_LO"] = rows["NUM_BRANCH_LO"].astype(int)
    rows["NUM_BRANCH_HI"] = rows["NUM_BRANCH_HI"].astype(int)
    rows["successful_placement"] = True

    # provenance from the tree list (kept on the source_row the generator
    # records): DBH/HT/TREE_CODE as given
    rows["DBH"] = pd.to_numeric(trees["DBH"], errors="coerce").round(4)
    rows["HT"] = ht.round(4)
    rows["TREE_CODE"] = trees["TREE_CODE"].astype(str)
    return rows.reset_index(drop=True)


def target_points(ht_hi: float, points_per_metre: int = POINTS_PER_METRE,
                  dbh: float = 1.0) -> int:
    """Point budget normalised by BOTH height and DBH.

    Np = points_per_metre * ceil(HT_HI) * dbh

    The height term keeps vertical point density constant; the dbh scalar keeps
    thin trees from being oversampled relative to thick ones (a thick trunk has
    more surface, so it should carry proportionally more points). dbh is the
    base diameter in metres.
    """
    n = points_per_metre * math.ceil(float(ht_hi)) * float(dbh)
    return max(1, int(round(n)))


_POSITIVE_COLS = ["HT_LO", "HT_HI", "BASE_DIA_LO", "BASE_DIA_HI", "MAX_CAN_W"]


def validate_row(row: pd.Series) -> str | None:
    """Return a reason string if the row can't be simulated, else None.

    Worth doing up front: trimesh will happily build a mesh from NaN
    geometry and hand back a cloud of NaN points instead of raising, which
    then fails silently all the way into the LAS file.
    """
    for col in _POSITIVE_COLS:
        v = pd.to_numeric(pd.Series([row.get(col)]), errors="coerce").iloc[0]
        if not np.isfinite(v) or v <= 0:
            return f"{col}={row.get(col)!r}"
    if float(row["HT_HI"]) < float(row["HT_LO"]):
        return "HT_HI < HT_LO"
    if float(row["BASE_DIA_HI"]) < float(row["BASE_DIA_LO"]):
        return "BASE_DIA_HI < BASE_DIA_LO"
    if not (np.isfinite(float(row["X"])) and np.isfinite(float(row["Y"]))):
        return "X/Y not finite"
    return None


def tree_seed(base_seed: int, tree_uid: str) -> int:
    """Deterministic per-tree seed, stable across runs and row ordering.

    Uses crc32 rather than hash() because hash() of a str is salted per
    process, which would silently break reproducibility.
    """
    return int((base_seed * 1000003 + zlib.crc32(str(tree_uid).encode())) % (2 ** 32))


def _tree_id_int(tree_uid: str, fallback: int) -> int:
    try:
        return int(str(tree_uid).lstrip("0") or 0)
    except ValueError:
        return int(fallback)


# ---------------------------------------------------------------------------
# Single tree
# ---------------------------------------------------------------------------

def generate_tree(row: pd.Series, plot_uid: str, out_dir: str,
                  base_seed: int = 0,
                  points_per_metre: int = POINTS_PER_METRE,
                  row_index: int = 0,
                  max_attempts: int = 5,
                  excursion_slack: float = 1.0,
                  overwrite: bool = True) -> dict:
    """Generate, label and write one tree. Returns a manifest record."""
    uid = str(row["TREE_UID"])
    las_path = layout.cloud_file(out_dir, f"{uid}_raw")
    gen_path = os.path.join(out_dir, f"{uid}_generation.json")
    skel_path = os.path.join(out_dir, f"{uid}_raw_skeleton.json")

    rec = {"tree_uid": uid, "plot_uid": str(plot_uid),
           "x": float(row["X"]), "y": float(row["Y"]),
           "func_group": str(row.get("FUNC_GROUP", "DEFAULT")),
           "n_points": 0, "status": "ok", "las": os.path.basename(las_path)}

    if not overwrite and os.path.exists(las_path):
        rec["status"] = "skipped_exists"
        return rec

    bad = validate_row(row)
    if bad is not None:
        rec["status"] = "invalid_params"
        rec["last_error"] = bad
        return rec

    model_params = params_from_row(row)
    stem_frac = stem_frac_from_row(row)
    np_target = target_points(row["HT_HI"], points_per_metre,
                              dbh=row["BASE_DIA_HI"])
    tid = _tree_id_int(uid, row_index + 1)

    # The simulator can fail on a degenerate randomisation (e.g. every branch
    # lands in the azimuth gap, leaving an empty foliage mesh). Reseed and
    # retry rather than losing the tree.
    points = skeleton = None
    seed = used_attempts = None
    for attempt in range(max_attempts):
        seed = tree_seed(base_seed + attempt * 7919, uid)
        np.random.seed(seed)
        try:
            points, skeleton = gen_simtree(
                Np=np_target, stem_frac_lambda=stem_frac,
                model_params=model_params, tx=float(row["X"]), ty=float(row["Y"]),
                return_skeleton=True)
            if not np.all(np.isfinite(points[:, :3])):
                raise ValueError("simulator produced non-finite point coordinates")
            # Guard against divergent stem splines (see the tx/ty fork fix in
            # treesim.py). A stem should never wander further from its base
            # than its own lean plus the fork offsets allows. Each cascading
            # split can offset a fork top by up to MAX_FORK_OFFSET, so the
            # limit scales with the row's configured split maximum.
            split_hi = 1
            if "SPLIT_NUM_HI" in row.index:
                split_hi = max(1, int(round(float(row["SPLIT_NUM_HI"]))))
            limit = (float(row["TOP_DIST"]) + split_hi * MAX_FORK_OFFSET
                     + excursion_slack)
            excursion = _max_stem_excursion(skeleton)
            if excursion > limit:
                raise ValueError(
                    f"stem excursion {excursion:.1f} m exceeds {limit:.1f} m")
            used_attempts = attempt + 1
            break
        except Exception as exc:                     # noqa: BLE001
            points = skeleton = None                 # don't keep a bad partial
            rec["last_error"] = f"{type(exc).__name__}: {exc}"

    if points is None:
        rec["status"] = "failed"
        return rec

    # provisional semantic: stem points -> tree (1), foliage -> foliage (3).
    # trunk(1) vs branch(2) is resolved later at assembly, once branch_code is
    # known; here every stem point is the coarse "tree" class. branch_code /
    # branch_id stay 0 and are assigned during assembly.
    import schema
    old_sem = points[:, 3].astype(int)               # 1 stem, 0 foliage
    sem = np.where(old_sem == 1, schema.SEM_TREE, schema.SEM_FOLIAGE).astype(np.int8)
    export_canonical_las(
        las_path, points[:, :3],
        {"tree_id": np.full(len(points), tid, dtype=np.int32),
         "semantic_code": sem},
        level="tree")

    generation = {
        "tree_uid": uid,
        "plot_uid": str(plot_uid),
        "tree_id": tid,
        "func_group": rec["func_group"],
        "x": float(row["X"]),
        "y": float(row["Y"]),
        "target_points": int(np_target),
        "points_per_metre": int(points_per_metre),
        "ht_hi_used_for_budget": float(row["HT_HI"]),
        "stem_frac_lambda": stem_frac,
        "model_params": model_params,
        "rng": {"library": "numpy.random (global)", "seed": int(seed),
                "base_seed": int(base_seed), "attempts": int(used_attempts)},
        "realised": {
            "n_points": int(len(points)),
            "n_stem_points": int((points[:, 3] == 1).sum()),
            "n_foliage_points": int((points[:, 3] == 0).sum()),
            "height_m": float(skeleton["height_m"]),
            "base_radius_m": float(skeleton["base_radius_m"]),
            "fork": int(skeleton["fork"]),
        },
        "outputs": {"las": os.path.basename(las_path),
                    "skeleton": os.path.basename(skel_path)},
        "source_row": {k: (None if pd.isna(v) else
                           (v.item() if hasattr(v, "item") else v))
                       for k, v in row.items()},
        "generated_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"),
        "generator": "SimpleSynthTree / treesim.gen_simtree (skeleton patch)",
    }
    with open(gen_path, "w") as f:
        json.dump(generation, f, indent=2)

    skeleton_out = {"tree_uid": uid, "plot_uid": str(plot_uid), "tree_id": tid,
                    "coordinate_frame": "plot metres, z = 0 at tree base",
                    **skeleton}
    with open(skel_path, "w") as f:
        json.dump(skeleton_out, f, indent=2)

    rec.update(n_points=int(len(points)),
               n_stem_points=int((points[:, 3] == 1).sum()),
               height_m=round(float(skeleton["height_m"]), 3),
               n_branches=int(skeleton["n_branches"]),
               fork=int(skeleton["fork"]),
               seed=int(seed))
    return rec


# ---------------------------------------------------------------------------
# Whole plot
# ---------------------------------------------------------------------------

def load_rows(table_path: str) -> tuple:
    """Rows for either contract: (kind, DataFrame).

    'tree-list-json' -> the minimal JSON contract (Section 8.1), expanded via
    the tree-code parameter dictionary at generation time (per-tree bounds
    from the tolerances; allometric crown width / base diameter from DBH/HT;
    seeded shape draws happen per TREE_UID in generate_tree below).
    'legacy-sst-csv' -> the frozen per-tree parameter CSV, consumed as-is.
    """
    if tree_list.is_tree_list_path(table_path):
        plot = tree_list.load_plot(table_path)
        return "tree-list-json", rows_from_tree_list(plot)
    return "legacy-sst-csv", pd.read_csv(table_path, dtype={"TREE_UID": str})


def generate_plot(table_path: str, out_dir: str | None = None,
                  seed: int = 0,
                  points_per_metre: int = POINTS_PER_METRE,
                  max_attempts: int = 5,
                  excursion_slack: float = 1.0,
                  overwrite: bool = True,
                  verbose: bool = True) -> pd.DataFrame:
    """Generate every tree in a plot table.

    Accepts a canonical JSON tree list (plot_*.json, Section 8.1) or a legacy
    SST parameter CSV. The output directory is derived from the table name
    (plot_001.json -> plot_001/) unless out_dir is given.
    """
    # TREE_UID must stay a string or "015" is read as 15 and the file names
    # stop matching the plot table.
    kind, df = load_rows(table_path)

    plot_uid = os.path.splitext(os.path.basename(table_path))[0]
    if out_dir is None:
        out_dir = os.path.join(os.path.dirname(os.path.abspath(table_path)), plot_uid)
    trees_dir = os.path.join(out_dir, "individual_trees")
    os.makedirs(trees_dir, exist_ok=True)

    n_in = len(df)
    if "successful_placement" in df.columns:
        keep = df["successful_placement"].astype(str).str.lower().isin(
            ["true", "1", "yes", "t"])
        df = df[keep]
    df = df.reset_index(drop=True)

    records = []
    for i, row in df.iterrows():
        rec = generate_tree(row, plot_uid, trees_dir, base_seed=seed,
                            points_per_metre=points_per_metre, row_index=i,
                            max_attempts=max_attempts,
                            excursion_slack=excursion_slack, overwrite=overwrite)
        records.append(rec)
        if verbose:
            note = (f"{rec['n_points']:>7d} pts  h={rec.get('height_m', 0):5.1f} m"
                    if rec["status"] == "ok" else f"  [{rec['status']}]")
            print(f"  {rec['tree_uid']}  {rec['func_group']:<8s}{note}")

    manifest = pd.DataFrame(records)
    manifest_path = os.path.join(out_dir, "manifest.csv")
    manifest.to_csv(manifest_path, index=False)

    if verbose:
        ok = int((manifest["status"] == "ok").sum()) if len(manifest) else 0
        print(f"{plot_uid}: {n_in} rows in, {len(df)} placed, {ok} generated "
              f"({int(manifest.get('n_points', pd.Series(dtype=int)).sum()):,} points) "
              f"-> {out_dir}")
    return manifest


def main():
    ap = argparse.ArgumentParser(description="Generate SimpleSynthTree clouds for a plot table")
    ap.add_argument("table", nargs="+",
                    help="plot tree-list JSON(s) (canonical contract) or "
                         "legacy SST parameter csv(s)")
    ap.add_argument("-o", "--outdir", default=None,
                    help="output root (default: alongside each table, named after it)")
    ap.add_argument("-s", "--seed", type=int, default=0)
    ap.add_argument("-p", "--points-per-metre", type=int, default=POINTS_PER_METRE)
    ap.add_argument("--max-attempts", type=int, default=5)
    ap.add_argument("--no-overwrite", action="store_true")
    args = ap.parse_args()

    for path in args.table:
        out = (os.path.join(args.outdir, os.path.splitext(os.path.basename(path))[0])
               if args.outdir else None)
        generate_plot(path, out_dir=out, seed=args.seed,
                      points_per_metre=args.points_per_metre,
                      max_attempts=args.max_attempts,
                      overwrite=not args.no_overwrite)


if __name__ == "__main__":
    main()