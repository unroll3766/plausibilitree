"""
vary_env.py: per-plot environment draws for a batch of plot tables.

The batch counterpart of replicate.py (spec Section 6.5): replicate.py
makes N replicates of ONE plot (same trees, re-drawn environment); this
tool varies the environment ACROSS plots. Given a directory of distinct
plot tables, every plot gets its own terrain / shrub / grass draw, written
as a plot_<UID>.env.json sidecar next to a byte-copy of its table (plus its
*_report.json sidecar, which CSV tables need for plot geometry) in the
output directory. run_pipeline.py picks the sidecars up automatically (one
env per plot), so a single pipeline run over the output directory generates
each plot with different topography and understory values while the trees
stay on the global --seed.

    python vary_env.py -i plots -o plots_var --seed 7

Terrain/ground keys are uniform draws over CLI ranges (reusing
replicate.py's machinery). The shrub/grass covers are a JOINT draw
conditioned on the plot's own trees: with BA = sum pi (DBH/2)^2 from the
table's DBH column and A the plot area, the total understory cover
T = shrub_cover + grass_cover is drawn uniformly from
[min(BA/A, 0.5), max(1 - BA/A, 0.5)] (clamped at BA/A = 0.5, where the
interval would be empty) and split uniformly, shrub = T*b, grass = T*(1-b)
with b ~ U(0,1), so the two always sum to T. The bound is on the nominal
parameters: realized understory area sits below T because grass is
attenuated under canopies (canopy_attenuation) and shrubs keep a tree
margin.

Inputs are positional paths and/or -i/--input (repeatable); output is
-o/--out-dir. A table that fails (unreadable, not a canonical table,
missing plot geometry) is skipped and recorded in the manifest's
`errors` list -- one bad table does not cancel the batch; the run
exits nonzero if any table failed.

Outputs in out_dir:
    plot_<UID>.csv|json        byte copy of the source table
    plot_<UID>_report.json     byte copy of the FIA report sidecar (CSVs)
    plot_<UID>.env.json        this plot's environment draw
    vary_env_manifest.json     sampler config + every draw + provenance
                               (+ "errors": failed input tables)

All draws derive from the run seed plus the plot UID, so they are
order-independent and reproducible (spec principle 6).

Requires: numpy
"""

from __future__ import annotations

import argparse
import glob
import hashlib
import json
import os
import shutil
import sys
import time

import numpy as np

import replicate
import tree_list
from fia_helper import plot_area_m2


def per_plot_seed(seed: int, plot_uid: str) -> int:
    """Deterministic per-plot seed from the run seed + plot UID, stable
    across runs and independent of table order (spec principle 6)."""
    digest = hashlib.sha256(f"{int(seed)}:{plot_uid}".encode()).digest()
    return int.from_bytes(digest[:4], "little")


def plot_basics(table_path, *, default_radius_m=7.315, default_shape="circle",
                default_length_m=None, default_width_m=None):
    """(basal_area_m2, plot_area_m2, geometry) for one plot table.

    BA comes from the table's DBH column (sum pi (DBH/2)^2, the same formula
    as build_info's BASAL_AREA but on list DBH). Plot geometry comes from
    the table itself (JSON tree list, or CSV + its *_report.json sidecar);
    CLI defaults fill in whatever the table does not carry.
    """
    plot = tree_list.load_plot(table_path)
    trees = plot["trees"]
    if "DBH" not in trees.columns:
        raise RuntimeError(f"{table_path}: no DBH column -- the joint cover "
                           f"bound needs per-tree DBH")
    dbh = trees["DBH"].astype(float)
    ba = float((np.pi * (dbh / 2.0) ** 2).sum())

    shape = plot.get("plot_shape") or default_shape
    radius = plot.get("plot_radius_m") or default_radius_m
    length = plot.get("plot_length_m") or default_length_m
    width = plot.get("plot_width_m") or default_width_m
    side = plot.get("plot_side_m")
    try:
        area = plot_area_m2(shape, radius_m=radius, length_m=length,
                            width_m=width, side_m=side)
    except TypeError:
        raise RuntimeError(
            f"{table_path}: incomplete plot geometry for shape {shape!r} "
            f"(radius={radius!r}, length={length!r}, width={width!r}, "
            f"side={side!r}); pass --plot-shape/--radius/--plot-length/"
            f"--plot-width or add a *_report.json sidecar") from None
    if area <= 0:
        raise RuntimeError(f"{table_path}: plot area {area} <= 0 (shape "
                           f"{shape!r}, radius {radius!r})")
    geom = {"plot_shape": shape, "plot_radius_m": radius,
            "plot_length_m": length, "plot_width_m": width,
            "plot_side_m": side}
    return ba, float(area), geom


def joint_cover_draw(rng, ba_frac):
    """One joint shrub/grass cover draw: returns (env_patch, record).

    T = shrub + grass ~ U[min(ba_frac, 0.5), max(1 - ba_frac, 0.5)], split
    uniformly (b ~ U(0,1)); the two covers are rounded so they sum to T
    exactly at 4 dp.
    """
    lo = min(ba_frac, 0.5)
    hi = max(1.0 - ba_frac, 0.5)
    total = round(float(rng.uniform(lo, hi)), 4)
    b = float(rng.uniform(0.0, 1.0))
    shrub = round(total * b, 4)
    grass = round(total - shrub, 4)
    patch = {"shrub_cover": shrub, "grass_cover": grass}
    rec = {"ba_frac": round(ba_frac, 6),
           "cover_bounds": [round(lo, 6), round(hi, 6)],
           "cover_total": total,
           "shrub_share": round(b, 4),
           "shrub_cover": shrub, "grass_cover": grass}
    return patch, rec


def draw_plot_env(rng, ranges, ba_frac, dtm_files=None):
    """One plot's full env draw: replicate.py's terrain/understory ranges,
    with the two cover keys replaced by the joint BA-conditioned draw."""
    env = replicate.draw_env(rng, ranges, dtm_files)
    patch, rec = joint_cover_draw(rng, ba_frac)
    env.update(patch)
    return env, rec


def vary(inputs, out_dir, *, seed=0, ranges=None, dtm_dir=None,
         default_radius_m=7.315, default_shape="circle",
         default_length_m=None, default_width_m=None,
         force=False, verbose=True):
    """Draw one environment per source plot table into out_dir. Returns the
    manifest dict (also written to <out_dir>/vary_env_manifest.json)."""
    ranges = {**replicate.DEFAULT_RANGES, **(ranges or {})}
    csvs = replicate.discover_tables(inputs)
    dtm_files = []
    if dtm_dir:
        dtm_files = sorted(glob.glob(os.path.join(dtm_dir, "*.asc")))
        if not dtm_files:
            raise FileNotFoundError(
                f"--dtm-dir {dtm_dir!r} contains no .asc DTMs "
                f"(looked for {os.path.join(dtm_dir, '*.asc')})")

    if os.path.exists(out_dir) and not os.path.isdir(out_dir):
        raise RuntimeError(f"output path {out_dir!r} exists and is not a "
                           f"directory; choose a different -o/--out-dir")
    if os.path.isdir(out_dir) and os.listdir(out_dir) and not force:
        raise RuntimeError(f"output directory {out_dir!r} is not empty "
                           f"({len(os.listdir(out_dir))} entries); pass "
                           f"--force to overwrite it")
    try:
        os.makedirs(out_dir, exist_ok=True)
    except OSError as e:
        raise RuntimeError(f"cannot create output directory {out_dir!r}: "
                           f"{e.strerror or e} (check that the parent path "
                           f"exists)") from None

    manifest = {
        "tool": "vary_env.py",
        "created": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "seed": int(seed),
        "ranges": {k: [float(v[0]), float(v[1])] for k, v in ranges.items()},
        "dtm_dir": (os.path.abspath(dtm_dir) if dtm_dir else None),
        "note": "shrub_cover + grass_cover are a joint draw: total T ~ "
                "U[min(BA/A,0.5), max(1-BA/A,0.5)], uniform split; bounds and "
                "split are recorded per plot",
        "plots": [],
        "errors": [],
    }

    for src in csvs:
        src = os.path.abspath(src)
        uid, ext = os.path.splitext(os.path.basename(src))
        try:
            if not (tree_list.is_tree_list_path(src) or tree_list.is_tree_list_csv(src)):
                raise RuntimeError(f"{src}: not a canonical plot table (legacy SST "
                                   f"parameter CSVs carry no DBH)")
            ba, area, geom = plot_basics(
                src, default_radius_m=default_radius_m, default_shape=default_shape,
                default_length_m=default_length_m, default_width_m=default_width_m)
            prng = np.random.default_rng(per_plot_seed(seed, uid))
            env, cover_rec = draw_plot_env(prng, ranges, ba / area, dtm_files)

            table = os.path.join(out_dir, f"{uid}{ext}")
            # CSV tables read their plot geometry from the report sidecar; copy
            # it too so the copy is self-contained
            report_src = os.path.splitext(src)[0] + "_report.json"
            report_dst = os.path.join(out_dir, f"{uid}_report.json")
            try:
                shutil.copyfile(src, table)
                if os.path.exists(report_src):
                    shutil.copyfile(report_src, report_dst)
            except OSError as e:
                raise RuntimeError(f"failed writing into output directory "
                                   f"{out_dir!r}: {e.strerror or e} (copied "
                                   f"{src!r} -> {table!r})") from None
            env_path = os.path.join(out_dir, f"{uid}.env.json")
            with open(env_path, "w") as f:
                json.dump({"plot_uid": uid, "source": src, "env": env}, f, indent=2)

            manifest["plots"].append({
                "plot_uid": uid,
                "source": src,
                "table": os.path.abspath(table),
                "report_copied": os.path.exists(report_src),
                "env_path": os.path.abspath(env_path),
                "basal_area_m2": round(ba, 6),
                "plot_area_m2": round(area, 6),
                "geometry": geom,
                **cover_rec,
                "env": env,
            })
            if verbose:
                print(f"  {uid}: BA/A={cover_rec['ba_frac']:.4f} "
                      f"T={cover_rec['cover_total']} in "
                      f"[{cover_rec['cover_bounds'][0]}, {cover_rec['cover_bounds'][1]}] "
                      f"-> shrub={env['shrub_cover']} grass={env['grass_cover']} "
                      f"bump={env['gen_bumpiness']} relief={env['gen_relief']} "
                      f"seed={env['seed']}")
        except (RuntimeError, ValueError, OSError) as e:
            # one bad table must not cancel the batch: record it, keep going
            manifest["errors"].append({"source": src, "error": str(e)})
            if verbose:
                print(f"  {uid}: skipped ({e})", file=sys.stderr)

    manifest_path = os.path.join(out_dir, "vary_env_manifest.json")
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)
    if verbose:
        print(f"\nwrote {len(manifest['plots'])} plot draw(s) from "
              f"{len(csvs)} table(s) -> {out_dir}")
        if manifest["errors"]:
            print(f"{len(manifest['errors'])} table(s) skipped; see 'errors' "
                  f"in the manifest", file=sys.stderr)
        print(f"manifest: {manifest_path}")
        print(f"next: python run_pipeline.py {out_dir} --seed 42 "
              f"(sidecars are applied per plot)")
    return manifest


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Draw one environment per plot table in a batch "
                    "(terrain ranges + joint tree-conditioned shrub/grass "
                    "cover), writing copied tables + env sidecars to out-dir")
    ap.add_argument("inputs", nargs="*",
                    help="plot table csv/json(s) and/or directories of them "
                         "(combined with -i/--input)")
    ap.add_argument("-i", "--input", dest="input", action="append", default=[],
                    help="input plot table or directory; repeatable. Use "
                         "-i INPUT -o OUT so the two paths cannot be "
                         "confused")
    ap.add_argument("-o", "--out-dir", default="plots_var",
                    help="directory for the copied tables + env sidecars + "
                         "manifest (default plots_var); must be empty "
                         "unless --force")
    ap.add_argument("--force", action="store_true",
                    help="overwrite a non-empty out-dir")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--dtm-dir", default=None,
                    help="directory of .asc DTMs; each plot draws one instead "
                         "of synthesizing terrain")
    g = ap.add_argument_group("plot geometry fallback",
                              "used only for tables that do not carry their "
                              "own geometry (hand-authored CSV without a "
                              "report sidecar); defaults match run_pipeline")
    g.add_argument("--plot-shape", default="circle",
                   choices=["circle", "rectangle", "square"])
    g.add_argument("--radius", type=float, default=7.315)
    g.add_argument("--plot-length", type=float, default=None)
    g.add_argument("--plot-width", type=float, default=None)
    r = ap.add_argument_group("draw ranges", "uniform [min:max] for the "
                                               "terrain/ground keys; a "
                                               "single value is fixed. "
                                               "shrub/grass cover is NOT "
                                               "drawn here -- it is the "
                                               "joint BA-conditioned draw")
    r.add_argument("--dtm-resolution", default=None, help="default 1.0:1.0")
    r.add_argument("--bumpiness", default=None,
                   help="DTM bumpiness, default 0.2:0.8")
    r.add_argument("--relief", default=None, help="DTM relief (m), default 3:15")
    r.add_argument("--ground-density", default=None,
                   help="ground points per DTM cell, default 0.5:2.0")
    r.add_argument("--ground-noise", default=None,
                   help="vertical ground noise half-range (m), default 0:0.03")
    r.add_argument("--shrub-cell-size", default=None, help="default 2.0:2.0")
    r.add_argument("--shrub-width", default=None,
                   help="range both width endpoints are drawn from; the pair "
                        "is sorted into min/max. default 0.8:1.6")
    r.add_argument("--grass-density", default=None,
                   help="grass candidate points per m2, default 100:300")
    r.add_argument("--attenuation", default=None,
                   help="grass attenuation under canopy, default 0.2:0.5")
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args(argv)

    inputs = list(args.inputs) + list(args.input)
    if not inputs:
        ap.error("no inputs: pass plot tables/directories positionally or "
                 "via -i/--input (output goes to -o/--out-dir)")

    range_args = {
        "dtm_resolution": args.dtm_resolution,
        "gen_bumpiness": args.bumpiness,
        "gen_relief": args.relief,
        "ground_density": args.ground_density,
        "ground_noise": args.ground_noise,
        # shrub_cover / grass_cover ranges are intentionally not settable:
        # the joint draw replaces them (they still get defaults below)
        "shrub_cover": None,
        "grass_cover": None,
        "shrub_cell_size": args.shrub_cell_size,
        "shrub_width": args.shrub_width,
        "grass_candidate_density": args.grass_density,
        "canopy_attenuation": args.attenuation,
    }
    ranges = {k: replicate.parse_range(v) for k, v in range_args.items()
              if v is not None}
    try:
        manifest = vary(inputs, args.out_dir, seed=args.seed, ranges=ranges,
                        dtm_dir=args.dtm_dir,
                        default_radius_m=args.radius, default_shape=args.plot_shape,
                        default_length_m=args.plot_length,
                        default_width_m=args.plot_width,
                        force=args.force, verbose=not args.quiet)
    except (RuntimeError, ValueError, OSError) as e:
        ap.exit(1, f"vary_env.py: error: {e}\n")
    if manifest.get("errors"):
        ap.exit(1, f"vary_env.py: error: {len(manifest['errors'])} input "
                   f"table(s) failed; the rest were written -- see 'errors' "
                   f"in {os.path.join(args.out_dir, 'vary_env_manifest.json')}\n")


if __name__ == "__main__":
    main()
