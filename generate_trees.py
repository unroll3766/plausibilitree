"""
generate_trees.py: generate one labelled point cloud per tree of a
PlausibiliTree tree list, through the generator module contract
(design spec Section 7.1).

    python generate_trees.py plot_11519.json --seed 42   # canonical tree list
    python generate_trees.py plot_11519.csv  --seed 42   # legacy SST table
    python generate_trees.py plot_11519.json -g sst     # module selection

The core selects a generator module (-g/--generator, registry in
generator_module.py; currently SST), loads its parameter dictionary,
resolves each tree's code, and builds the tree via module.build() with a
per-tree rng seeded from a stable hash of the plot seed + TREE_UID
(retried with rehashed seeds up to MAX_ATTEMPTS on degenerate builds).
The module returns a base-origin TreeModel; the core places it at the
tree's (X, Y) via module.place() and writes, per tree, under
<outdir>/individual_trees/:

    <TREE_UID>_raw.las            labelled point cloud, plot coordinates
    <TREE_UID>_generation.json    every parameter needed to reproduce it
    <TREE_UID>_raw_skeleton.json  the module's skeleton format (SST: 7.4)

plus a manifest.csv summarising the run.

Requires: numpy, pandas, trimesh, laspy (the latter two via treesim / las).
"""

from __future__ import annotations

import argparse
import datetime as _dt
import json
import os
import zlib

import numpy as np
import pandas as pd

from las import export_canonical_las
from generator_module import get_generator
import layout
import tree_list

POINTS_PER_METRE = 1024


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
# Tree requests: both contracts -> one per-tree request
# ---------------------------------------------------------------------------

def tree_requests(table_path: str, module, points_per_metre: int) -> tuple:
    """A plot table -> (kind, [request...]).

    Each request is a dict: TREE_UID, X, Y, DBH, HT, TREE_CODE and the
    resolved build `params` (a canonical entry for tree-list JSONs, an
    explicit frozen-range dict for legacy SST rows -- both feed the same
    module.build path).
    """
    if tree_list.is_tree_list_path(table_path) or tree_list.is_tree_list_csv(table_path):
        plot = tree_list.load_plot(table_path)
        trees = plot["trees"]
        reqs = []
        for _, t in trees.iterrows():
            code = str(t["TREE_CODE"])
            reqs.append({"TREE_UID": str(t["TREE_UID"]),
                         "X": float(t["X"]), "Y": float(t["Y"]),
                         "DBH": float(t["DBH"]), "HT": float(t["HT"]),
                         "TREE_CODE": code,
                         "params": dict(module.entry_for(code),
                                        points_per_metre=points_per_metre)})
        return "tree-list", reqs
    df = pd.read_csv(table_path, dtype={"TREE_UID": str})
    if "successful_placement" in df.columns:
        df = df[df["successful_placement"].astype(str).str.lower()
                .isin(["true", "1", "yes", "t"])]
    df = df.reset_index(drop=True)
    reqs = []
    for _, row in df.iterrows():
        reqs.append({"TREE_UID": str(row["TREE_UID"]),
                     "X": float(row["X"]), "Y": float(row["Y"]),
                     "DBH": float(row["BASE_DIA_HI"]),
                     "HT": float(row["HT_HI"]),
                     "TREE_CODE": str(row.get("FUNC_GROUP", "DEFAULT")),
                     "params": module.legacy_row_params(row,
                                                        points_per_metre),
                     "source_row": {k: (None if pd.isna(v) else
                                        (v.item() if hasattr(v, "item") else v))
                                    for k, v in row.items()}})
    return "legacy-sst-csv", reqs


# ---------------------------------------------------------------------------
# Single tree
# ---------------------------------------------------------------------------

def generate_tree(module, req: dict, plot_uid: str, out_dir: str,
                  base_seed: int = 0,
                  points_per_metre: int = POINTS_PER_METRE,
                  row_index: int = 0,
                  max_attempts: int = 10,
                  overwrite: bool = True) -> dict:
    """Build, place and write one tree. Returns a manifest record."""
    uid = str(req["TREE_UID"])
    las_path = layout.cloud_file(out_dir, f"{uid}_raw")
    gen_path = os.path.join(out_dir, f"{uid}_generation.json")
    skel_path = os.path.join(out_dir, f"{uid}_raw_skeleton.json")

    rec = {"tree_uid": uid, "plot_uid": str(plot_uid),
           "x": float(req["X"]), "y": float(req["Y"]),
           "func_group": str(req["TREE_CODE"]),
           "n_points": 0, "status": "ok", "las": os.path.basename(las_path)}

    if not overwrite and os.path.exists(las_path):
        rec["status"] = "skipped_exists"
        return rec

    tid = _tree_id_int(uid, row_index + 1)
    attempts = []
    model = None
    for attempt in range(max_attempts):
        seed = tree_seed(base_seed + attempt * 7919, uid)
        rng = np.random.default_rng(seed)
        attempts.append({"seed": int(seed)})
        try:
            model = module.build(rng, req["DBH"], req["HT"],
                                 req["TREE_CODE"], req["params"])
            break
        except Exception as exc:                     # noqa: BLE001
            model = None
            rec["last_error"] = f"{type(exc).__name__}: {exc}"

    if model is None:
        rec["status"] = "failed"
        return rec

    model = module.place(model, req["X"], req["Y"])
    points = model.points                       # (N, 4): x, y, z, semantic

    export_canonical_las(
        las_path, points[:, :3],
        {"tree_id": np.full(len(points), tid, dtype=np.int32),
         "semantic_code": model.labels},
        level="tree")

    generation = {
        "tree_uid": uid,
        "plot_uid": str(plot_uid),
        "tree_id": tid,
        "generator_module": getattr(module, "name", "unknown"),
        "func_group": str(req["TREE_CODE"]),
        "x": float(req["X"]),
        "y": float(req["Y"]),
        "tree_list": {"DBH": float(req["DBH"]), "HT": float(req["HT"]),
                      "TREE_CODE": str(req["TREE_CODE"])},
        "attempts": attempts,
        "realised": dict(model.realized),
        "provenance": dict(model.provenance),
        "outputs": {"las": os.path.basename(las_path),
                    "skeleton": os.path.basename(skel_path)},
        "source_row": req.get("source_row"),
        "generated_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"),
    }
    with open(gen_path, "w") as f:
        json.dump(generation, f, indent=2)

    skeleton_out = {"tree_uid": uid, "plot_uid": str(plot_uid), "tree_id": tid,
                    "coordinate_frame": "plot metres, z = 0 at tree base",
                    **model.skeleton}
    with open(skel_path, "w") as f:
        json.dump(skeleton_out, f, indent=2)

    rec.update(n_points=int(model.realized.get("n_points", len(points))),
               height_m=round(float(model.realized.get("height_m", 0.0)), 3),
               seed=int(attempts[-1]["seed"]))
    return rec


# ---------------------------------------------------------------------------
# Whole plot
# ---------------------------------------------------------------------------

def generate_plot(table_path: str, out_dir: str | None = None,
                  seed: int = 0,
                  points_per_metre: int = POINTS_PER_METRE,
                  max_attempts: int = 10,
                  generator: str = "sst",
                  overwrite: bool = True,
                  verbose: bool = True) -> pd.DataFrame:
    """Generate every tree in a plot table (canonical JSON tree list or
    legacy SST parameter CSV) through the selected generator module.

    The output directory is derived from the table name (plot_001.json ->
    plot_001/) unless out_dir is given.
    """
    module = get_generator(generator)
    kind, reqs = tree_requests(table_path, module, points_per_metre)

    plot_uid = os.path.splitext(os.path.basename(table_path))[0]
    if out_dir is None:
        out_dir = os.path.join(os.path.dirname(os.path.abspath(table_path)), plot_uid)
    trees_dir = os.path.join(out_dir, "individual_trees")
    os.makedirs(trees_dir, exist_ok=True)

    records = []
    for i, req in enumerate(reqs):
        rec = generate_tree(module, req, plot_uid, trees_dir,
                            base_seed=seed,
                            points_per_metre=points_per_metre,
                            row_index=i, max_attempts=max_attempts,
                            overwrite=overwrite)
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
        print(f"{plot_uid} [{kind}, generator={generator}]: "
              f"{len(reqs)} trees, {ok} generated "
              f"({int(manifest.get('n_points', pd.Series(dtype=int)).sum()):,} points) "
              f"-> {out_dir}")
    return manifest


def main():
    ap = argparse.ArgumentParser(
        description="Generate labelled tree clouds for a plot via a "
                    "generator module (spec Section 7.1)")
    ap.add_argument("table", nargs="+",
                    help="plot tree-list JSON(s) (canonical contract) or "
                         "legacy SST parameter csv(s)")
    ap.add_argument("-o", "--outdir", default=None,
                    help="output root (default: alongside each table, named after it)")
    ap.add_argument("-s", "--seed", type=int, default=0)
    ap.add_argument("-p", "--points-per-metre", type=int, default=POINTS_PER_METRE)
    ap.add_argument("-g", "--generator", default="sst",
                    help="generator module (registry in generator_module.py)")
    ap.add_argument("--max-attempts", type=int, default=10)
    ap.add_argument("--no-overwrite", action="store_true")
    args = ap.parse_args()

    for path in args.table:
        out = (os.path.join(args.outdir, os.path.splitext(os.path.basename(path))[0])
               if args.outdir else None)
        generate_plot(path, out_dir=out, seed=args.seed,
                      points_per_metre=args.points_per_metre,
                      max_attempts=args.max_attempts,
                      generator=args.generator,
                      overwrite=not args.no_overwrite)


if __name__ == "__main__":
    main()
