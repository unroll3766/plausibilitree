"""
generate_trees.py: generate one SimpleSynthTree point cloud per row of a
PlausibiliTree plot table.

    python generate_trees.py plot_11519.csv --seed 42

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

def generate_plot(csv_path: str, out_dir: str | None = None,
                  seed: int = 0,
                  points_per_metre: int = POINTS_PER_METRE,
                  max_attempts: int = 5,
                  excursion_slack: float = 1.0,
                  overwrite: bool = True,
                  verbose: bool = True) -> pd.DataFrame:
    """Generate every successfully placed tree in a plot table.

    The output directory is derived from the csv name (plot_001.csv ->
    plot_001/) unless out_dir is given.
    """
    # TREE_UID must stay a string or "015" is read as 15 and the file names
    # stop matching the plot table.
    df = pd.read_csv(csv_path, dtype={"TREE_UID": str})

    plot_uid = os.path.splitext(os.path.basename(csv_path))[0]
    if out_dir is None:
        out_dir = os.path.join(os.path.dirname(os.path.abspath(csv_path)), plot_uid)
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
    ap.add_argument("csv", nargs="+", help="plot table csv(s) from fia_helper")
    ap.add_argument("-o", "--outdir", default=None,
                    help="output root (default: alongside each csv, named after it)")
    ap.add_argument("-s", "--seed", type=int, default=0)
    ap.add_argument("-p", "--points-per-metre", type=int, default=POINTS_PER_METRE)
    ap.add_argument("--max-attempts", type=int, default=5)
    ap.add_argument("--no-overwrite", action="store_true")
    args = ap.parse_args()

    for path in args.csv:
        out = (os.path.join(args.outdir, os.path.splitext(os.path.basename(path))[0])
               if args.outdir else None)
        generate_plot(path, out_dir=out, seed=args.seed,
                      points_per_metre=args.points_per_metre,
                      max_attempts=args.max_attempts,
                      overwrite=not args.no_overwrite)


if __name__ == "__main__":
    main()