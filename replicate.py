"""
replicate.py: environment replicates of plot inventory tables.

A replicate keeps the inventory rows (the same trees, at the same positions)
and re-draws the environment around them: terrain shape (or a choice of DTM
from a directory), ground sampling, and the shrub / grass layers. Each
replicate is a new plot table plus a JSON "env sidecar" holding its draw;
run_pipeline.py reads the sidecar automatically and applies that plot's
settings to terrain.py and understory.py, leaving every other stage alone.

    # 3 environment replicates of one plot, drawn from the default ranges
    python replicate.py plots/plot_11574.csv -o plots_rep -n 3 --seed 7

    # replicate every plot table in a directory (replicate tables excluded)
    python replicate.py plots -o plots_rep -n 4 --seed 7

    # draw real DTMs from a directory instead of synthesizing terrain
    python replicate.py plots -o plots_rep -n 4 --dtm-dir dtms

    # then run the pipeline exactly as usual; sidecars are picked up per plot
    python run_pipeline.py plots_rep --radius 15 --seed 42

Outputs in out_dir:
    plot_<UID>_r<NN>.csv          the source table, rows unchanged (byte copy)
    plot_<UID>_r<NN>.env.json     the environment draw for this replicate
    replicates_manifest.json      sampler config + every draw (provenance)

Design notes:
  - Trees are identical across replicates: generate_trees / assemble_plot
    run with run_pipeline's global --seed and are unaffected by the sidecar.
    Variation comes from terrain.py / understory.py consuming the sidecar's
    seed and parameters, so replicates differ in ground shape, shrub cover
    and grass -- the canopy over them is shared.
  - Every draw is uniform over a [min:max] range given on the CLI; a single
    value means fixed. Shrub width is drawn as a sorted pair, so
    width_min <= width_max always holds.
  - The env JSON keys mirror run_pipeline.py's kwarg names, so a sidecar is
    also a convenient hand-run record: the terrain.py / understory.py command
    for one replicate can be reconstructed from it directly.

Requires: numpy
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import re
import shutil
import time

import numpy as np

# keys written into each env sidecar; run_pipeline.py reads exactly these
ENV_KEYS = [
    "seed",
    "gen_cell_size", "gen_bumpiness", "gen_relief",
    "ground_density", "ground_noise",
    "shrub_cover", "shrub_cell_size", "shrub_width_min", "shrub_width_max",
    "grass_cover", "grass_candidate_density", "canopy_attenuation",
]
# "dtm" is added only when --dtm-dir is given (a path, not a draw).

# default draw ranges: uniform over (min, max); equal values = fixed
DEFAULT_RANGES = {
    "gen_cell_size": (1.0, 1.0),
    "gen_bumpiness": (0.2, 0.8),
    "gen_relief": (3.0, 15.0),
    "ground_density": (0.5, 2.0),
    "ground_noise": (0.0, 0.03),
    "shrub_cover": (0.05, 0.35),
    "shrub_cell_size": (2.0, 2.0),
    "shrub_width": (0.8, 1.6),           # drawn as a sorted min/max pair
    "grass_cover": (0.3, 0.7),
    "grass_candidate_density": (100.0, 300.0),
    "canopy_attenuation": (0.2, 0.5),
}

# replicate tables written by an earlier run, excluded when a directory is
# the input (so replicating a directory never re-replicates old replicates)
_RE_REPLICATE = re.compile(r"plot_.+_r\d+\.csv$")


def parse_range(text):
    """'min:max' (or a single fixed value) -> (float, float)."""
    parts = str(text).split(":")
    if len(parts) == 1:
        lo = hi = float(parts[0])
    elif len(parts) == 2:
        lo, hi = float(parts[0]), float(parts[1])
    else:
        raise ValueError(f"range must be 'min:max' or 'value', got {text!r}")
    if hi < lo:
        lo, hi = hi, lo
    return lo, hi


def draw_env(rng, ranges, dtm_files=None):
    """One environment draw. Keys mirror run_pipeline.py's kwargs."""
    def u(key):
        lo, hi = ranges[key]
        return round(float(rng.uniform(lo, hi)), 4)

    env = {
        "seed": int(rng.integers(1, 2**31 - 1)),
        "gen_cell_size": u("gen_cell_size"),
        "gen_bumpiness": u("gen_bumpiness"),
        "gen_relief": u("gen_relief"),
        "ground_density": u("ground_density"),
        "ground_noise": u("ground_noise"),
        "shrub_cover": u("shrub_cover"),
        "shrub_cell_size": u("shrub_cell_size"),
        "grass_cover": u("grass_cover"),
        "grass_candidate_density": u("grass_candidate_density"),
        "canopy_attenuation": u("canopy_attenuation"),
    }
    w = sorted((u("shrub_width"), u("shrub_width")))
    env["shrub_width_min"] = w[0]
    env["shrub_width_max"] = w[1]
    if dtm_files:
        env["dtm"] = str(rng.choice(sorted(dtm_files)))
    return env


def discover_tables(inputs):
    """Resolve inputs (csvs and/or directories) into source table paths."""
    csvs = []
    for p in inputs:
        if os.path.isdir(p):
            found = sorted(glob.glob(os.path.join(p, "plot_*.csv")))
            fresh = [f for f in found if not _RE_REPLICATE.match(os.path.basename(f))]
            if not fresh and found:
                raise RuntimeError(
                    f"{p}: only replicate tables (plot_*_rNN.csv) found; pass the "
                    f"original tables or a directory holding them")
            csvs += fresh
        else:
            if not os.path.exists(p):
                raise FileNotFoundError(p)
            csvs.append(p)
    if not csvs:
        raise RuntimeError("no plot tables found in the given input(s)")
    return csvs


def make_replicates(inputs, out_dir, n, *, seed=0, ranges=None, dtm_dir=None,
                    verbose=True):
    """Write n replicate tables + env sidecars per source table. Returns the
    manifest dict (also written to <out_dir>/replicates_manifest.json)."""
    if n < 1:
        raise ValueError("n must be >= 1")
    ranges = {**DEFAULT_RANGES, **(ranges or {})}
    csvs = discover_tables(inputs)
    dtm_files = []
    if dtm_dir:
        dtm_files = sorted(glob.glob(os.path.join(dtm_dir, "*.asc")))
        if not dtm_files:
            raise FileNotFoundError(f"no .asc DTMs under {dtm_dir}")

    os.makedirs(out_dir, exist_ok=True)
    rng = np.random.default_rng(seed)
    pad = max(2, len(str(n)))

    manifest = {
        "tool": "replicate.py",
        "created": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "seed": int(seed),
        "n_replicates": int(n),
        "ranges": {k: [float(v[0]), float(v[1])] for k, v in ranges.items()},
        "dtm_dir": (os.path.abspath(dtm_dir) if dtm_dir else None),
        "replicates": [],
    }

    for src in csvs:
        src = os.path.abspath(src)
        uid = os.path.splitext(os.path.basename(src))[0]
        # tables are copied byte-for-byte: the rows (trees) are the point of a
        # replicate, and only the plot uid (the filename) changes
        with open(src) as f:
            header = f.readline()
        if "TREE_UID" not in header:
            raise RuntimeError(f"{src}: not a plot table (no TREE_UID column)")

        for k in range(1, n + 1):
            new_uid = f"{uid}_r{k:0{pad}d}"
            table = os.path.join(out_dir, f"{new_uid}.csv")
            env_path = os.path.join(out_dir, f"{new_uid}.env.json")
            env = draw_env(rng, ranges, dtm_files)
            shutil.copyfile(src, table)
            with open(env_path, "w") as f:
                json.dump({"plot_uid": new_uid, "source": src, "env": env},
                          f, indent=2)
            manifest["replicates"].append({
                "plot_uid": new_uid,
                "source": src,
                "table": os.path.abspath(table),
                "env_path": os.path.abspath(env_path),
                "env": env,
            })
            if verbose:
                dtm_note = f", dtm {os.path.basename(env['dtm'])}" if "dtm" in env \
                    else ""
                print(f"  {new_uid}: bump={env['gen_bumpiness']} "
                      f"relief={env['gen_relief']} "
                      f"shrub_cover={env['shrub_cover']} "
                      f"grass={env['grass_cover']}{dtm_note} "
                      f"seed={env['seed']}")

    manifest_path = os.path.join(out_dir, "replicates_manifest.json")
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)
    if verbose:
        print(f"\nwrote {len(manifest['replicates'])} replicate(s) from "
              f"{len(csvs)} table(s) -> {out_dir}")
        print(f"manifest: {manifest_path}")
        print("next: python run_pipeline.py "
              f"{out_dir} --radius 15 --seed 42  (sidecars are applied per plot)")
    return manifest


def main():
    ap = argparse.ArgumentParser(
        description="Make environment replicates of plot inventory tables "
                    "(same trees, randomly drawn terrain/shrub/grass settings)")
    ap.add_argument("inputs", nargs="+",
                    help="plot table csv(s) and/or directories of them")
    ap.add_argument("-o", "--out-dir", default="plots_rep",
                    help="directory for the replicate tables + sidecars "
                         "(default plots_rep)")
    ap.add_argument("-n", "--replicates", type=int, default=2,
                    help="replicates per source table (default 2)")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--dtm-dir", default=None,
                    help="directory of .asc DTMs; each replicate draws one "
                         "instead of synthesizing terrain")
    r = ap.add_argument_group("draw ranges", "uniform [min:max]; a single "
                                                "value is fixed")
    r.add_argument("--gen-cell-size", default=None, help="default 1.0:1.0")
    r.add_argument("--bumpiness", default=None,
                   help="DTM bumpiness, default 0.2:0.8")
    r.add_argument("--relief", default=None, help="DTM relief (m), default 3:15")
    r.add_argument("--ground-density", default=None,
                   help="ground points per DTM cell, default 0.5:2.0")
    r.add_argument("--ground-noise", default=None,
                   help="vertical ground noise half-range (m), default 0:0.03")
    r.add_argument("--shrub-cover", default=None,
                   help="fraction of cells hosting a shrub, default 0.05:0.35")
    r.add_argument("--shrub-cell-size", default=None, help="default 2.0:2.0")
    r.add_argument("--shrub-width", default=None,
                   help="range both width endpoints are drawn from; the pair "
                        "is sorted into min/max. default 0.8:1.6")
    r.add_argument("--grass-cover", default=None,
                   help="fraction of grass candidate cells kept, default 0.3:0.7")
    r.add_argument("--grass-density", default=None,
                   help="grass candidate points per m2, default 100:300")
    r.add_argument("--attenuation", default=None,
                   help="grass attenuation under canopy, default 0.2:0.5")
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()

    range_args = {
        "gen_cell_size": args.gen_cell_size,
        "gen_bumpiness": args.bumpiness,
        "gen_relief": args.relief,
        "ground_density": args.ground_density,
        "ground_noise": args.ground_noise,
        "shrub_cover": args.shrub_cover,
        "shrub_cell_size": args.shrub_cell_size,
        "shrub_width": args.shrub_width,
        "grass_cover": args.grass_cover,
        "grass_candidate_density": args.grass_density,
        "canopy_attenuation": args.attenuation,
    }
    ranges = {k: parse_range(v) for k, v in range_args.items() if v is not None}
    try:
        make_replicates(args.inputs, args.out_dir, args.replicates,
                        seed=args.seed, ranges=ranges, dtm_dir=args.dtm_dir,
                        verbose=not args.quiet)
    except (RuntimeError, FileNotFoundError, ValueError) as e:
        ap.exit(1, f"replicate.py: {e}\n")


if __name__ == "__main__":
    main()
