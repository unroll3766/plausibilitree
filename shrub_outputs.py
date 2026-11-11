"""
shrub_outputs.py: per-shrub individual outputs, mirroring the per-tree ones.

For every shrub in a plot this writes, under individual_shrubs/<tag>/:
  las/<tag>_raw.las          the shrub's own points only
  las/<tag>_understory.las   the shrub + surrounding context (crop of the scene)
  deep_learning/canopy_map/...  crown CHM + bbox/hull (foliage-style, from shrub pts)
  deep_learning/stem_map/...    breast-height (DBH-band) cross-section (may be empty)

Each shrub's volume is NOT a per-object file (2026-11 review, item W24):
the sampled-shape volumes (primitive + mesh) plus the geometry are carried
by the plot-level shrub index
(plot_placement/textual/<plot>_shrub_index.json) and, from there, into each
shrub's info JSON (build_info.py).
"""

from __future__ import annotations

import json
import os

import numpy as np

import georef
import layout
import schema
from las import read_canonical_las, export_canonical_las
from stem_map import RasterGrid, _write_geotiff, BREAST_HEIGHT_M, STEM_SECTION
from canopy_map import canopy_hull


def _shrub_uid(sid):
    return f"{int(sid):03d}"


def build_shrub_outputs(plot_dir, plot_uid, shrub_geom,
                        canopy_cell_size=0.25, stem_cell_size=0.25,
                        dbh_height=BREAST_HEIGHT_M, dbh_threshold=STEM_SECTION,
                        verbose=False):
    """Write per-shrub outputs from the plot understory cloud.

    shrub_geom: {shrub_id: {height_m, width_m, center_xy, primitive_volume_m3,
                            mesh_volume_m3, ...}} from gen_shrub_layer.
    """
    understory_las = layout.plot_path(plot_dir, plot_uid, "understory_las",
                                      create=False)
    if not os.path.exists(understory_las):
        if verbose:
            print(f"{plot_uid}: no understory cloud; skipping shrub outputs")
        return 0

    xyz, cols = read_canonical_las(understory_las, level="plot")
    sem = cols["semantic_code"].astype(int)
    sid = cols["shrub_id"].astype(int)

    shrub_pt = sem == schema.SEM_SHRUB
    written = 0
    shrub_index = []

    for s in sorted(shrub_geom):
        geom = shrub_geom[s]
        uid = _shrub_uid(s)
        mine = shrub_pt & (sid == s)
        n_pts = int(mine.sum())

        # (per-shrub sub-clouds retired, item W6: the plot understory
        # cloud + shrub_id are the record; extraction crops them on demand)

        # ---- canopy map (crown envelope from shrub points) ----
        _shrub_canopy_map(plot_dir, plot_uid, uid, xyz[mine], canopy_cell_size)

        # ---- DBH/stem map (breast-height band cross-section; may be empty) ----
        _shrub_stem_map(plot_dir, plot_uid, uid, xyz[mine], stem_cell_size,
                        dbh_height, dbh_threshold)

        # full per-shrub record for the plot-level index: the geometry and
        # the sampled-shape volumes (primitive + mesh) that build_info
        # copies into the shrub info JSONs (volume JSONs are retired,
        # 2026-11 review item W24)
        shrub_index.append({
            "shrub_id": int(s), "shrub_uid": uid, "n_points": n_pts,
            "height_m": geom["height_m"],
            "width_m": geom["width_m"],
            "base_z_m": geom.get("base_z_m", 0.0),
            "center_xy": geom["center_xy"],
            "volume_m3": {
                "primitive": geom["primitive_volume_m3"],
                "mesh": geom["mesh_volume_m3"],
            },
        })
        written += 1

    # plot-level shrub index in the scene textual dir
    idx_dir = os.path.join(plot_dir, layout.PLOT_PLACEMENT, "textual")
    os.makedirs(idx_dir, exist_ok=True)
    with open(os.path.join(idx_dir, f"{plot_uid}_shrub_index.json"), "w") as f:
        json.dump({"plot_uid": plot_uid, "n_shrubs": written,
                   "shrubs": shrub_index}, f, indent=2)

    if verbose:
        print(f"{plot_uid}: wrote {written} individual shrub outputs")
    return written


def _shrub_canopy_map(plot_dir, plot_uid, uid, pts, cell_size):
    """Crown CHM + convex-hull/bbox JSON for a single shrub."""
    if len(pts) < 1:
        return
    anchor = georef.load_anchor(plot_dir, plot_uid)
    x, y, z = pts[:, 0], pts[:, 1], pts[:, 2]
    grid = RasterGrid(x.min(), y.min(), x.max(), y.max(), cell_size)
    # CHM feature: max height per cell
    col, row = grid.to_pixel(x, y)
    chm = np.zeros((grid.height, grid.width), dtype=np.float32)
    flat = row * grid.width + col
    order = np.argsort(z)
    chm.reshape(-1)[flat[order]] = z[order].astype(np.float32)
    x0, y0 = grid.x_min, grid.y_max
    _write_geotiff(layout.shrub_path(plot_dir, plot_uid, uid, "canopy_feature_raw"),
                   chm, x0, y0, grid.cell, anchor=anchor)
    # label mask: cells with any shrub point
    mask = (chm > 0).astype(np.uint8)
    _write_geotiff(layout.shrub_path(plot_dir, plot_uid, uid, "canopy_label_raw"),
                   mask, x0, y0, grid.cell, anchor=anchor)
    # bbox + hull JSON
    bbox = {"x_min": float(x.min()), "x_max": float(x.max()),
            "y_min": float(y.min()), "y_max": float(y.max()),
            "z_min": float(z.min()), "z_max": float(z.max())}
    with open(layout.shrub_path(plot_dir, plot_uid, uid, "canopy_bbox_json"), "w") as f:
        json.dump({"shrub_id": int(uid), "bbox_world": bbox}, f, indent=2)
    hull = {"shrub_id": int(uid)}
    if len(pts) >= 4:
        try:
            _, _, vol, area = canopy_hull(pts)
            hull.update({"convex_hull_volume_m3": float(vol),
                         "convex_hull_area_m2": float(area)})
        except Exception:
            hull.update({"convex_hull_volume_m3": None, "convex_hull_area_m2": None})
    with open(layout.shrub_path(plot_dir, plot_uid, uid, "canopy_hull_json"), "w") as f:
        json.dump(hull, f, indent=2)


def _shrub_stem_map(plot_dir, plot_uid, uid, pts, cell_size, dbh_height, threshold):
    """Breast-height band cross-section for a shrub (DBH-map analog).

    Uses the same DBH band parameters as trees. A shrub has no trunk, so this is
    often empty or a small patch near the shrub's mid-height -- that's expected
    and intentional (kept for parity with the tree DBH map).
    """
    if len(pts) < 1:
        return
    z = pts[:, 2]
    z0 = z.min()
    band = (z >= z0 + dbh_height - threshold) & (z <= z0 + dbh_height + threshold)
    x, y = pts[:, 0], pts[:, 1]
    grid = RasterGrid(x.min(), y.min(), x.max(), y.max(), cell_size)
    feat = grid.counts(x[band], y[band]) if band.any() else \
        np.zeros((grid.height, grid.width), dtype=np.uint16)
    x0, y0 = grid.x_min, grid.y_max
    anchor = georef.load_anchor(plot_dir, plot_uid)
    _write_geotiff(layout.shrub_path(plot_dir, plot_uid, uid, "stem_feature_raw"),
                   feat, x0, y0, grid.cell, anchor=anchor)
    mask = (feat > 0).astype(np.uint8)
    _write_geotiff(layout.shrub_path(plot_dir, plot_uid, uid, "stem_label_raw"),
                   mask, x0, y0, grid.cell, anchor=anchor)
