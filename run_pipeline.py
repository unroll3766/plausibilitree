"""
run_pipeline.py: run the PlausibiliTree pipeline over a directory of
already-generated per-plot tables.

fia_helper.py is now a separate, standalone step -- run it first to turn a raw
FIA TREE table into <plotdir>/ tree lists. This orchestrator only consumes
those tables and runs the processing stages, in order, for every one:

    generate_trees -> assemble_plot -> restructure -> terrain
                    (always; flat plane when no DTM is configured)
                    -> understory (opt.) -> build_info
    derived, opt-in (spec Section 5): stem_map (--stem-map),
                    canopy_map (--canopy-map), normalize (--normalize,
                    dataset-wide over the restructured plots)

Terrain always runs per spec Section 3.3: a --dtm file or --generate-dtm
selects the relief, and neither means the None topography (a flat plane at
z = 0; a warning applies if no ground density was configured).

Each stage is invoked as its own subprocess (same as running it by hand), so
each stage's own arguments, logging and error handling apply unchanged. All
plots use the same parameters.

Example (mirrors the documented per-stage commands):

    python fia_helper.py CA_TREE_SMALL.csv -o plots -r 15 --basis stem --seed 42

    python run_pipeline.py plots \\
        --radius 7.315 --seed 42 \\
        --points-per-metre 4096 --max-attempts 10 \\
        --stem-map --stem-resolution 0.05 --stem-section 0.25 \\
        --canopy-map --chm-resolution 0.05 \\
        --max-density 5000 --spatial --density cube --plot-l 30 \\
        --tree-l 30 --normalize

Parameters not given fall back to the defaults below, which match that example.
By default a failing plot is logged and skipped so the batch continues; use
--stop-on-error to abort instead. Use --limit to process only the first N of
the plot tables found in plotdir.

Replicated plots (see replicate.py) carry a plot_<UID>.env.json sidecar next
to their table; this orchestrator applies that draw's terrain and understory
settings to that plot only, while every other stage and every sidecar-free
plot keeps the global settings unchanged.

Requires the stage scripts to be importable/runnable from the same directory.
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import shutil
import subprocess
import sys
import time

import tree_list

HERE = os.path.dirname(os.path.abspath(__file__))
PY = sys.executable


def _run(script, args, dry_run=False):
    """Run `python <script> <args...>`; return True on success."""
    cmd = [PY, os.path.join(HERE, script), *[str(a) for a in args]]
    if dry_run:
        print("    " + " ".join(cmd))
        return True
    proc = subprocess.run(cmd, capture_output=True, text=True)
    out = (proc.stdout or "").strip()
    if out:
        # stage scripts print a one-line summary; indent it under the stage
        for line in out.splitlines():
            print(f"      {line}")
    if proc.returncode != 0:
        err = (proc.stderr or "").strip().splitlines()
        tail = err[-1] if err else f"exit code {proc.returncode}"
        print(f"      ERROR: {tail}")
        return False
    return True


# ---------------------------------------------------------------------------
# per-plot environment sidecars (written by replicate.py)
# ---------------------------------------------------------------------------
# A plot_<UID>.env.json next to a plot table overrides the global terrain and
# understory settings for that plot only, so replicated tables (same trees,
# different environment) can each draw their own ground, shrub and grass
# layer. Keys mirror run_pipeline's kwarg names; "dtm" (a path) replaces DTM
# synthesis for that plot. See replicate.py / README_replicate.md.

ENV_KEYS = {
    "seed", "dtm",
    "dtm_resolution", "gen_bumpiness", "gen_relief",
    "ground_density", "ground_noise",
    "shrub_cover", "shrub_cell_size", "shrub_width_min", "shrub_width_max",
    "grass_cover", "grass_candidate_density", "canopy_attenuation",
}


def _load_env(csv):
    """Load a plot's env sidecar, or None if absent/unreadable."""
    path = os.path.splitext(csv)[0] + ".env.json"
    if not os.path.exists(path):
        return None
    try:
        with open(path) as f:
            doc = json.load(f)
    except (OSError, ValueError) as e:
        print(f"    WARNING: unreadable env sidecar {path} ({e}); "
              f"using global settings")
        return None
    env = doc.get("env", doc)          # replicate.py wraps the draw in "env"
    unknown = sorted(k for k in env if k not in ENV_KEYS)
    if unknown:
        print(f"    WARNING: {path}: ignoring unknown key(s) "
              f"{', '.join(unknown)}")
    return {k: v for k, v in env.items() if k in ENV_KEYS}


def run_pipeline(plotdir, *,
                 radius=7.315, seed=42,
                 plot_shape="circle", plot_length=None, plot_width=None,
                 limit=None, generator="sst",
                 points_per_metre=4096, max_attempts=10,
                 branch_radius_frac=0.10, tolerance=0.0,
                 stem_resolution=0.05, stem_section=0.25, bbox_padding=1,
                 chm_resolution=0.05, canopy_buffer=2, canopy_source="trimmed",
                 dtm=None, generate_dtm=False, dtm_resolution=0.25,
                 gen_bumpiness=0.5, gen_relief=10.0,
                 ground_density=1.0, ground_noise=0.0,
                 terrain_source="trimmed",
                 understory=False, shrub_cover=0.15, shrub_cell_size=2.0,
                 shrub_width_min=0.8, shrub_width_max=1.6,
                 foliage_voxel=0.1, grass_voxel=0.1, voxel_threshold=1,
                 grass_cover=0.6, grass_candidate_density=200.0,
                 canopy_attenuation=0.3,
                 rotation_step=30.0,
                 max_density=5000, spatial=True, density="none",
                 plot_l=30, plot_z=None, tree_l=30, tree_z=None,
                 norm_source="understory",
                 stem_map=False, canopy_map=False, normalize=False,
                 stop_on_error=False, dry_run=False):
    t0 = time.time()
    if density == "none":          # accept the CLI-style sentinel programmatically
        density = None

    # ---- clear stale outputs ----------------------------------------------
    # fia_helper.py is a separate step now and owns plot_<id>.csv in plotdir;
    # this orchestrator only writes the per-plot plot_<id>/ stage directories
    # beside them. A previous run of this orchestrator can leave plot_*/
    # directories behind that this run would silently pick up via the glob
    # below, so clear those -- but never touch the plot_*.csv tables, which
    # are input here, not output.
    if os.path.isdir(plotdir) and not dry_run:
        stale_dir = [d for d in glob.glob(os.path.join(plotdir, "plot_*"))
                     if os.path.isdir(d)]
        if stale_dir:
            print(f"[cleanup] removing {len(stale_dir)} stale plot dir(s) in "
                  f"{plotdir}/")
            for d in stale_dir:
                shutil.rmtree(d)

    # ---- discover already-produced plot tables -----------------------------
    # fia_helper.py is run separately beforehand to produce these; see the
    # module docstring for the two-step invocation.
    plot_csvs = tree_list.find_plots(plotdir)
    if limit is not None:
        plot_csvs = plot_csvs[:limit]
    if not plot_csvs:
        print(f"no plot_*.json tree lists (or legacy plot_*.csv tables) "
              f"found in {plotdir}/; run fia_helper.py first, e.g.\n"
              f"  python fia_helper.py <fia_csv> -o {plotdir} -r {radius}")
        return 1
    print(f"[plots] found {len(plot_csvs)} plot table(s) in {plotdir}/\n")

    # ---- per-plot stages --------------------------------------------------
    # Stages 2-5 + restructure run per plot. Terrain (optional) runs after
    # restructure. Normalize is deferred and run once over all plots together,
    # so its branch-order scan spans the whole dataset and every _normalized.las
    # shares one column schema (uniform across the dataset, not just within a
    # plot).
    #
    # Env sidecars (replicate.py): a plot may carry its own terrain/understory
    # draw. Builders below merge the sidecar over the globals per plot and
    # return None when the stage isn't configured for that plot, so a mixed
    # directory (some replicated, some not) works unchanged.
    envs = {csv: _load_env(csv) for csv in plot_csvs}
    any_env = any(envs.values())

    stages = []
    stages.append(("generate_trees.py", lambda c: [
        c, "--seed", seed, "--points-per-metre", points_per_metre,
        "--generator", generator,
        "--max-attempts", max_attempts]))
    _assemble_shape = (["--plot-shape", plot_shape]
                       + (["--length", plot_length, "--width", plot_width]
                          if plot_shape == "rectangle" else []))
    stages.append(("assemble_plot.py", lambda c: [
        c, "--seed", seed, "--plot-radius", radius,
        "--branch-radius-frac", branch_radius_frac, "--tolerance", tolerance,
        "--rotation-step", rotation_step,
        *_assemble_shape]))
    # derived stages are opt-in (spec Section 5, item W6): stem/canopy
    # maps run only when requested; restructure tolerates their absence
    if stem_map:
        stages.append(("stem_map.py", lambda c: [
            c, "--stem-resolution", stem_resolution,
            "--stem-section", stem_section,
            "--bbox-padding", bbox_padding]))
    if canopy_map:
        stages.append(("canopy_map.py", lambda c: [
            c, "--chm-resolution", chm_resolution, "--buffer", canopy_buffer,
            "--source", canopy_source]))
    stages.append(("restructure.py", lambda c: [c]))

    # terrain is optional and runs after restructure, before normalize.
    # Sidecar "dtm" (if any) wins over synthesis; a sidecar without "dtm"
    # synthesizes with the sidecar's drawn gen_* parameters.
    def _terrain_args(csv):
        env = envs.get(csv) or {}
        dtm_p = env.get("dtm", dtm)
        common = [csv,
                  "--density", env.get("ground_density", ground_density),
                  "--noise", env.get("ground_noise", ground_noise),
                  "--source", terrain_source,
                  "--seed", env.get("seed", seed)]
        if dtm_p is not None:
            return common + ["--dtm", dtm_p]
        if generate_dtm or env:
            return common + ["--generate-dtm",
                             "--dtm-resolution", env.get("dtm_resolution", dtm_resolution),
                             "--gen-bumpiness", env.get("gen_bumpiness", gen_bumpiness),
                             "--gen-relief", env.get("gen_relief", gen_relief)]
        # no topography source: the None topography (spec Section 3.3) --
        # terrain.py drapes onto a flat plane at z = 0 and warns if no
        # density was configured
        return common

    # understory (shrubs + grass) is optional; terrain always runs
    # (flat plane when no topography is configured, spec Section 3.3)
    def _understory_args(csv):
        env = envs.get(csv) or {}
        terr = _terrain_args(csv)
        if not (understory or env):
            return None                 # not requested for this plot
        args = [csv,
                "--seed", env.get("seed", seed),
                "--shrub-cover", env.get("shrub_cover", shrub_cover),
                "--shrub-cell-size", env.get("shrub_cell_size", shrub_cell_size),
                "--shrub-width-min", env.get("shrub_width_min", shrub_width_min),
                "--shrub-width-max", env.get("shrub_width_max", shrub_width_max),
                "--grass-cover", env.get("grass_cover", grass_cover),
                "--grass-candidate-density",
                env.get("grass_candidate_density", grass_candidate_density),
                "--canopy-attenuation",
                env.get("canopy_attenuation", canopy_attenuation)]
        if "--dtm" in terr:             # terrain used a provided DTM
            args += ["--dtm", terr[terr.index("--dtm") + 1]]
        return args

    # terrain always runs (after restructure, before normalize): with a
    # provided DTM, a generated DTM, or -- when neither is configured -- the
    # flat plane (None topography, spec Section 3.3). Sidecar "dtm" (if any)
    # wins over synthesis; a sidecar without "dtm" synthesizes with the
    # sidecar's drawn gen_* parameters.
    stages.append(("terrain.py", _terrain_args))
    want_understory = (understory
                      or any("shrub_cover" in (e or {}) for e in envs.values()))
    if want_understory:
        stages.append(("understory.py", _understory_args))

    # info JSONs (plot + per-tree) are always built last, after every cloud,
    # skeleton and map exists. Replaces the old per-tree/plot volume JSONs.
    info_args = ["--foliage-voxel", foliage_voxel,
                 "--grass-voxel", grass_voxel,
                 "--voxel-threshold", voxel_threshold]
    stages.append(("build_info.py", lambda c: [c, *info_args]))

    ok, failed = [], []
    restructured = []          # plot CSVs that reached restructure, for normalize
    if any_env:
        n_env = sum(1 for e in envs.values() if e)
        print(f"[env] {n_env} of {len(plot_csvs)} plot(s) carry an env sidecar "
              f"(per-plot terrain/understory settings)\n")
    for i, csv in enumerate(plot_csvs, 1):
        name = os.path.splitext(os.path.basename(csv))[0]
        print(f"[{i}/{len(plot_csvs)}] {name}")
        plot_failed = False
        for script, build in stages:
            args = build(csv)
            if args is None:
                print(f"    {script}: skipped (not configured for this plot)")
                continue
            print(f"    {script}")
            if not _run(script, args, dry_run):
                plot_failed = True
                if stop_on_error:
                    print(f"\nAborting: {script} failed on {name} "
                          f"(--stop-on-error).")
                    return 2
                print(f"    skipping rest of {name}")
                break
        if plot_failed:
            failed.append(name)
        else:
            ok.append(name)
            restructured.append(csv)
        print()

    # ---- normalize: one call over every restructured plot -----------------
    # (opt-in derived stage, spec Section 5, item W6)
    if normalize and restructured:
        norm = ["--max-density", max_density] if max_density is not None else []
        if spatial:
            norm.append("--spatial")
        if density is not None:
            norm += ["--density", density, "--plot-l", plot_l, "--tree-l", tree_l]
            if density == "square":
                norm += ["--plot-z", plot_z, "--tree-z", tree_z]
        norm += ["--source", norm_source, "--seed", seed]
        print(f"[normalize] {len(restructured)} plot(s) together "
              f"(dataset-wide uniform schema)")
        if not _run("normalize.py", [*restructured, *norm], dry_run):
            print("      normalize failed on the batch")
            # move any ok plots that only failed at normalize into failed
            failed.extend(os.path.splitext(os.path.basename(c))[0]
                          for c in restructured)
            ok = [n for n in ok if n not in failed]
        print()

    dt = time.time() - t0
    print("=" * 60)
    print(f"done in {dt:.0f}s: {len(ok)} ok, {len(failed)} failed, "
          f"of {len(plot_csvs)} plots")
    if failed:
        print("failed plots: " + ", ".join(failed))
    return 0 if not failed else 3


def main():
    ap = argparse.ArgumentParser(
        description="Run the PlausibiliTree pipeline over a directory of "
                    "plot_*.csv tables produced by fia_helper.py")
    ap.add_argument("plotdir",
                    help="directory of plot_*.csv tables (run fia_helper.py "
                         "first to produce them)")
    # shared / assemble_plot geometry
    ap.add_argument("-r", "--radius", type=float, default=7.315)
    ap.add_argument("--plot-shape", choices=["circle", "rectangle"],
                    default="circle", help="plot geometry (default circle)")
    ap.add_argument("--length", type=float, default=None,
                    help="rectangle length in m (x), for --plot-shape rectangle")
    ap.add_argument("--width", type=float, default=None,
                    help="rectangle width in m (y), for --plot-shape rectangle")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--limit", type=int, default=None,
                    help="process only the first N plot tables found in plotdir")
    # generate_trees
    ap.add_argument("--points-per-metre", type=int, default=4096)
    ap.add_argument("--max-attempts", type=int, default=10)
    ap.add_argument("-g", "--generator", default="sst",
                    help="generator module (spec Section 7.1; registry in "
                         "generator_module.py)")
    # assemble_plot
    ap.add_argument("--branch-radius-frac", type=float, default=0.10)
    ap.add_argument("--tolerance", type=float, default=0.0)
    ap.add_argument("--rotation-step", type=float, default=30.0,
                    help="rotation quantization in degrees; a multiple of it "
                         "is drawn from [0, 360) per tree (default 30)")
    # stem_map
    ap.add_argument("--stem-resolution", type=float, default=0.05,
                    help="stem map cell size, m (spec STEM_RESOLUTION)")
    ap.add_argument("--stem-section", type=float, default=0.25,
                    help="DBH z-band half-thickness, m (spec STEM_SECTION)")
    ap.add_argument("--bbox-padding", type=int, default=1)
    # canopy_map
    ap.add_argument("--chm-resolution", type=float, default=0.05,
                    help="CHM cell size, m (spec CHM_RESOLUTION)")
    ap.add_argument("--canopy-buffer", type=int, default=2)
    ap.add_argument("--canopy-source", default="trimmed", choices=["raw", "trimmed"])
    # terrain (always runs; the None topography is a flat plane)
    ap.add_argument("--dtm", default=None,
                    help="DTM .asc (metres) to drape onto (wins over "
                         "--generate-dtm)")
    ap.add_argument("--generate-dtm", action="store_true",
                    help="synthesize a per-plot DTM instead of providing one; "
                         "enables the terrain step")
    ap.add_argument("--dtm-resolution", type=float, default=0.25,
                    help="generated/flat DTM cell size, m (spec DTM_RESOLUTION)")
    ap.add_argument("--gen-bumpiness", type=float, default=0.5)
    ap.add_argument("--gen-relief", type=float, default=10.0)
    ap.add_argument("--ground-density", type=float, default=1.0,
                    help="ground points sampled per DTM cell (default 1)")
    ap.add_argument("--ground-noise", type=float, default=0.0,
                    help="vertical noise half-range for ground points, m (default 0)")
    ap.add_argument("--terrain-source", default="trimmed", choices=["raw", "trimmed"])
    # understory (optional; terrain always runs -- flat plane when no
    # topography source is configured)
    ap.add_argument("--understory", action="store_true",
                    help="add shrubs + grass after terrain")
    ap.add_argument("--shrub-cover", type=float, default=0.15)
    ap.add_argument("--shrub-cell-size", type=float, default=2.0)
    ap.add_argument("--shrub-width-min", type=float, default=0.8)
    ap.add_argument("--shrub-width-max", type=float, default=1.6)
    ap.add_argument("--grass-cover", type=float, default=0.6)
    ap.add_argument("--grass-candidate-density", type=float, default=200.0)
    ap.add_argument("--canopy-attenuation", type=float, default=0.3)
    ap.add_argument("--foliage-voxel", type=float, default=0.1,
                    help="FOLIAGE_VOXEL: voxel edge (m) for foliage volume "
                         "accounting in build_info (default 0.1)")
    ap.add_argument("--grass-voxel", type=float, default=0.1,
                    help="GRASS_VOXEL: voxel edge (m) for grass volume "
                         "accounting in build_info (default 0.1)")
    ap.add_argument("--voxel-threshold", type=int, default=1,
                    help="min points for an occupied voxel (default 1)")
    # normalize
    ap.add_argument("--max-density", type=int, default=5000)
    ap.add_argument("--spatial", action="store_true", default=True)
    ap.add_argument("--no-spatial", dest="spatial", action="store_false")
    ap.add_argument("--density", default="none", choices=["square", "cube", "none"],
                    help="density-dependent voxelization: none (default; write "
                         "the subsampled/rescaled point cloud), cube (L^3 voxel "
                         "lattice), or square (LxLxZ column-profile lattice). "
                         "Voxelization is opt-in.")
    ap.add_argument("--plot-l", type=int, default=30)
    ap.add_argument("--plot-z", type=int, default=None)
    ap.add_argument("--tree-l", type=int, default=30)
    ap.add_argument("--tree-z", type=int, default=None)
    ap.add_argument("--norm-source", default="understory",
                    choices=["raw", "trimmed", "terrain", "understory"],
                    help="cloud normalize reads; default resolves to the "
                         "furthest-processed present (understory->terrain->"
                         "trimmed->raw)")
    # opt-in derived stages (spec Section 5, item W6)
    ap.add_argument("--stem-map", action="store_true",
                    help="run the stem map derived stage (2-band breast-"
                         "height raster + bounding boxes)")
    ap.add_argument("--canopy-map", action="store_true",
                    help="run the canopy map derived stage (2-band CHM + "
                         "crown bboxes/hulls)")
    ap.add_argument("--normalize", action="store_true",
                    help="run the normalization derived stage over the "
                         "restructured plots")
    # control
    ap.add_argument("--stop-on-error", action="store_true",
                    help="abort the batch on the first stage failure")
    ap.add_argument("--dry-run", action="store_true",
                    help="print the commands without running them")
    args = ap.parse_args()

    density = None if args.density == "none" else args.density
    if density == "square" and (args.plot_z is None or args.tree_z is None):
        ap.error("--density square requires --plot-z and --tree-z")
    if args.dtm and args.generate_dtm:
        print("warning: both --dtm and --generate-dtm given; using the "
              "provided DTM (filepath wins, generation ignored) per the "
              "spec's topography resolution order")
    if args.plot_shape == "rectangle" and (args.length is None or args.width is None):
        ap.error("--plot-shape rectangle requires --length and --width")

    rc = run_pipeline(
        args.plotdir,
        radius=args.radius,
        plot_shape=args.plot_shape, plot_length=args.length,
        plot_width=args.width,
        seed=args.seed, limit=args.limit, generator=args.generator,
        points_per_metre=args.points_per_metre, max_attempts=args.max_attempts,
        branch_radius_frac=args.branch_radius_frac, tolerance=args.tolerance,
        rotation_step=args.rotation_step,
        stem_resolution=args.stem_resolution, stem_section=args.stem_section,
        bbox_padding=args.bbox_padding,
        chm_resolution=args.chm_resolution, canopy_buffer=args.canopy_buffer,
        canopy_source=args.canopy_source,
        dtm=args.dtm, generate_dtm=args.generate_dtm,
        dtm_resolution=args.dtm_resolution, gen_bumpiness=args.gen_bumpiness,
        gen_relief=args.gen_relief, ground_density=args.ground_density,
        ground_noise=args.ground_noise, terrain_source=args.terrain_source,
        understory=args.understory, shrub_cover=args.shrub_cover,
        shrub_width_min=args.shrub_width_min, shrub_width_max=args.shrub_width_max,
        foliage_voxel=args.foliage_voxel,
        grass_voxel=args.grass_voxel,
        voxel_threshold=args.voxel_threshold,
        shrub_cell_size=args.shrub_cell_size, grass_cover=args.grass_cover,
        grass_candidate_density=args.grass_candidate_density,
        canopy_attenuation=args.canopy_attenuation,
        max_density=args.max_density, spatial=args.spatial, density=density,
        plot_l=args.plot_l, plot_z=args.plot_z, tree_l=args.tree_l,
        tree_z=args.tree_z, norm_source=args.norm_source,
        stem_map=args.stem_map, canopy_map=args.canopy_map,
        normalize=args.normalize,
        stop_on_error=args.stop_on_error, dry_run=args.dry_run)
    sys.exit(rc)


if __name__ == "__main__":
    main()
