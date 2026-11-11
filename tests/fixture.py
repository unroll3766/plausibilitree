"""Shared fixture builder for the PlausibiliTree regression suites.

Builds small, deterministic fixtures in a temporary directory using the
pipeline's programmatic entry points (no CLI subprocesses):

  - build_tree_list(): FIA placement (radius 20, seed 42) -> canonical
    6-tree JSON tree list copied from CA_TREE_SMALL.csv (the repo's sample
    inventory; 111 documented trees, TPA_UNADJ 1.0).
  - build_generated_plot(): run generate_trees over that list at 64
    points/metre (the same fixture the VERIFICATION_LOG W7 rounds used).
  - build_full_fixture(): a complete run_pipeline pass with a generated
    DTM and understory -- the base for the cross-artifact invariant suite.

Fixtures land in a fresh tempfile.mkdtemp() directory unless --outdir is
given; --keep on the test CLIs skips cleanup for debugging.
"""
from __future__ import annotations

import json
import os
import shutil
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO not in sys.path:
    sys.path.insert(0, REPO)

import pandas as pd

from fia_helper import fia_to_tree_list

FIA_SAMPLE = os.path.join(REPO, "CA_TREE_SMALL.csv")
PLOT_UID = "plot_11839"
SEED = 42
RADIUS_M = 20.0
N_TREES = 6


def workdir(keep: bool = False, tag: str = "plt") -> str:
    """A scratch dir; kept when `keep`, removed at exit otherwise.

    Returns (path, fuller) where fuller() removes it later (a no-op for
    kept dirs)."""
    path = tempfile.mkdtemp(prefix=f"{tag}_tests_")
    if keep:
        return path, lambda: None
    return path, lambda: shutil.rmtree(path, ignore_errors=True)


def build_tree_list(outdir: str, n_trees: int = N_TREES,
                    radius: float = RADIUS_M, seed: int = SEED) -> str:
    """FIA sample -> placed canonical tree list (n trees) at outdir."""
    df = pd.read_csv(FIA_SAMPLE, low_memory=False)
    plot = fia_to_tree_list(df, plot_radius_m=radius, seed=seed,
                            plot_uid=PLOT_UID)
    plot["trees"] = plot["trees"][:n_trees]
    path = os.path.join(outdir, f"{PLOT_UID}.json")
    with open(path, "w") as f:
        json.dump(plot, f, indent=2)
    return path


def build_generated_plot(table_path: str, ppm: int = 64,
                        seed: int = SEED) -> str:
    """generate_trees over the tree list; returns the plot dir."""
    import generate_trees
    generate_trees.generate_plot(table_path, out_dir=None, seed=seed,
                                 points_per_metre=ppm, overwrite=True,
                                 verbose=False)
    base = os.path.dirname(os.path.abspath(table_path))
    return os.path.join(base, PLOT_UID)


def build_full_fixture(outdir: str, ppm: int = 64, seed: int = SEED) -> str:
    """Full pipeline (generated DTM + understory) on the 6-tree list.

    run_pipeline() returns the batch exit code (0 = every plot ok)."""
    from run_pipeline import run_pipeline
    rc = run_pipeline(outdir, radius=RADIUS_M, seed=seed,
                     points_per_metre=ppm, limit=1,
                     generate_dtm=True, understory=True,
                     stem_map=True, canopy_map=True, normalize=True)
    if rc not in (0, None):
        raise RuntimeError(f"fixture pipeline failed with exit code {rc}")
    return os.path.join(outdir, PLOT_UID)


def build_core_fixture(outdir: str, ppm: int = 16,
                       seed: int = SEED) -> str:
    """Core-only pipeline (no derived stages, item W6): the tree list is
    generated, assembled, restructured, terrain-draped and info-built, with
    NO stem/canopy maps and NO normalization. Returns the plot dir."""
    from run_pipeline import run_pipeline
    build_tree_list(outdir, n_trees=2)
    rc = run_pipeline(outdir, radius=RADIUS_M, seed=seed,
                      points_per_metre=ppm, limit=1)
    if rc not in (0, None):
        raise RuntimeError(f"core fixture pipeline failed with exit code {rc}")
    return os.path.join(outdir, PLOT_UID)


GEO_CENTER = {"lon": -120.6, "lat": 38.9, "elev": 1500.0}


def build_geo_fixture(outdir: str, ppm: int = 32, seed: int = 7) -> str:
    """A small georeferenced plot (anchored tree list, spec Section 10)
    for the W5 geo-export checks. Returns the plot dir."""
    trees = [{"TREE_UID": "001", "X": -3.0, "Y": 2.0,
              "DBH": 0.45, "HT": 25.0, "TREE_CODE": "DF"},
             {"TREE_UID": "002", "X": 4.0, "Y": -1.5,
              "DBH": 0.35, "HT": 20.0, "TREE_CODE": "PP"},
             {"TREE_UID": "003", "X": 1.0, "Y": 5.0,
              "DBH": 0.30, "HT": 18.0, "TREE_CODE": "HW"}]
    path = os.path.join(outdir, "plot_geo.json")
    with open(path, "w") as f:
        json.dump({"plot_uid": "plot_geo", "plot_shape": "circle",
                   "plot_radius_m": 12.0, "plot_center": GEO_CENTER,
                   "trees": trees}, f, indent=2)
    from run_pipeline import run_pipeline
    # shrub_cover/width raised so the shrub raster exports reliably engage
    rc = run_pipeline(outdir, radius=12.0, seed=seed, points_per_metre=ppm,
                      limit=1, generate_dtm=True, understory=True,
                      stem_map=True, canopy_map=True, normalize=True,
                      shrub_cover=0.9, shrub_width_min=0.8,
                      shrub_width_max=1.6)
    if rc not in (0, None):
        raise RuntimeError(f"geo fixture pipeline failed (exit {rc})")
    return os.path.join(outdir, "plot_geo")
