"""
shrub_outputs.py: per-shrub individual outputs, mirroring the per-tree ones.

For every shrub in a plot this writes, under individual_shrubs/<tag>/:
  las/<tag>_raw.las          the shrub's own points only
  las/<tag>_understory.las   the shrub + surrounding context (crop of the scene)
  deep_learning/volume/labels/<tag>_volume.json
        shrub_id, height/width, and three volumes:
          primitive  -- nominal ellipsoid (4/3)pi(w/2)^2(h/2)
          mesh       -- actual jittered primitive mesh volume
          voxel      -- occupied-voxel count * occupied_voxel_size^3
  deep_learning/canopy_map/...  crown CHM + bbox/hull (foliage-style, from shrub pts)
  deep_learning/stem_map/...    breast-height (DBH-band) cross-section (may be empty)

The (raw / understory / normalized) cloud triple matches trees so a normalize
run treats shrubs exactly like trees.
"""

from __future__ import annotations

import json
import os

import numpy as np

import layout
import schema
from las import read_canonical_las, export_canonical_las
from stem_map import RasterGrid, _write_geotiff, BREAST_HEIGHT_M, DEFAULT_THRESHOLD
from canopy_map import canopy_hull


def _shrub_uid(sid):
    return f"{int(sid):03d}"


def voxel_volume(xyz, voxel_size, threshold):
    """Occupied-voxel volume: count voxels holding >= threshold points.

    Returns (volume_m3, n_occupied, n_total_nonempty). A voxel is occupied when
    it contains at least `threshold` of the given points; volume is that count
    times voxel_size**3.
    """
    if len(xyz) == 0 or voxel_size <= 0:
        return 0.0, 0, 0
    keys = np.floor(xyz / voxel_size).astype(np.int64)
    # unique voxel -> point count
    _, counts = np.unique(keys, axis=0, return_counts=True)
    occ = int((counts >= threshold).sum())
    return occ * (voxel_size ** 3), occ, int(len(counts))


def build_shrub_outputs(plot_dir, plot_uid, shrub_geom,
                        occupied_voxel_size=0.1, occupied_voxel_threshold=1,
                        canopy_cell_size=0.25, stem_cell_size=0.25,
                        dbh_height=BREAST_HEIGHT_M, dbh_threshold=DEFAULT_THRESHOLD,
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

        # ---- raw cloud: this shrub's points only ----
        raw_out = layout.shrub_path(plot_dir, plot_uid, uid, "raw_las")
        export_canonical_las(
            raw_out, xyz[mine],
            {"shrub_id": cols["shrub_id"][mine],
             "semantic_code": cols["semantic_code"][mine],
             "classification": cols["classification"][mine],
             "height_above_ground": cols["height_above_ground"][mine]},
            level="tree")

        # ---- understory cloud: shrub + surrounding context (bbox crop) ----
        # context window = shrub footprint bbox padded by its width, all points
        # (any class) inside it, so the shrub sits in its neighbourhood.
        cx, cy = geom["center_xy"]
        half = max(geom["width_m"], 0.1)
        pad = half  # one shrub-width of context around the footprint
        in_box = ((xyz[:, 0] >= cx - half - pad) & (xyz[:, 0] <= cx + half + pad) &
                  (xyz[:, 1] >= cy - half - pad) & (xyz[:, 1] <= cy + half + pad))
        u_out = layout.shrub_path(plot_dir, plot_uid, uid, "understory_las")
        export_canonical_las(
            u_out, xyz[in_box],
            {"tree_id": cols["tree_id"][in_box],
             "shrub_id": cols["shrub_id"][in_box],
             "semantic_code": cols["semantic_code"][in_box],
             "classification": cols["classification"][in_box],
             "branch_code": cols["branch_code"][in_box],
             "branch_id": cols["branch_id"][in_box],
             "height_above_ground": cols["height_above_ground"][in_box]},
            level="tree")

        # ---- volume JSON: primitive + mesh + voxel ----
        vvol, vocc, vtot = voxel_volume(xyz[mine], occupied_voxel_size,
                                        occupied_voxel_threshold)
        vol_doc = {
            "shrub_id": int(s),
            "n_points": n_pts,
            "height_m": geom["height_m"],
            "width_m": geom["width_m"],
            "base_z_m": geom.get("base_z_m", 0.0),
            "center_xy": geom["center_xy"],
            "volume_m3": {
                "primitive": geom["primitive_volume_m3"],
                "mesh": geom["mesh_volume_m3"],
                "voxel": vvol,
            },
            "voxel_volume_params": {
                "occupied_voxel_size_m": occupied_voxel_size,
                "occupied_voxel_threshold": occupied_voxel_threshold,
                "n_occupied_voxels": vocc,
                "n_nonempty_voxels": vtot,
            },
        }
        vol_out = layout.shrub_path(plot_dir, plot_uid, uid, "volume_label")
        with open(vol_out, "w") as f:
            json.dump(vol_doc, f, indent=2)

        # ---- canopy map (crown envelope from shrub points) ----
        _shrub_canopy_map(plot_dir, plot_uid, uid, xyz[mine], canopy_cell_size)

        # ---- DBH/stem map (breast-height band cross-section; may be empty) ----
        _shrub_stem_map(plot_dir, plot_uid, uid, xyz[mine], stem_cell_size,
                        dbh_height, dbh_threshold)

        shrub_index.append({
            "shrub_id": int(s), "shrub_uid": uid, "n_points": n_pts,
            "volume_m3": vol_doc["volume_m3"],
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
                   chm, x0, y0, grid.cell)
    # label mask: cells with any shrub point
    mask = (chm > 0).astype(np.uint8)
    _write_geotiff(layout.shrub_path(plot_dir, plot_uid, uid, "canopy_label_raw"),
                   mask, x0, y0, grid.cell)
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
    _write_geotiff(layout.shrub_path(plot_dir, plot_uid, uid, "stem_feature_raw"),
                   feat, x0, y0, grid.cell)
    mask = (feat > 0).astype(np.uint8)
    _write_geotiff(layout.shrub_path(plot_dir, plot_uid, uid, "stem_label_raw"),
                   mask, x0, y0, grid.cell)
