"""
restructure.py: reorganise a plot's flat stage outputs into the canonical
layout (see layout.py) and generate the derived artifacts the layout implies
but the stages don't emit yet:

  - per-tree canopy bbox / hull JSON (split out of the plot-level combined files)
  - volume JSONs are RETIRED (per-tree and plot-level, raw + trimmed):
    superseded by build_info's per-tree _info.json; restructure only
    prints the canopy totals now, and the scene JSON carries no link
  - plot-level stem mask  (features/labels pair with the stem map)
  - plot-level canopy mask (features/labels pair with the canopy height map)

Run after stages 2-5:

    python restructure.py plot_11519.csv

Existing files are MOVED (not copied) into the new tree, and the scene JSON is
rewritten with relative links to every artifact's new home. The plot-level
stem bounding-box list, which has no dedicated file in the new layout, is folded
into the scene JSON. Normalized variants are part of the layout but are not
produced here.

The plot-level masks are INSTANCE labels keyed by tree_id (0 = background),
uint32, aligned pixel-for-pixel with their feature raster:
  - stem mask:   tree_id of the stem points in each breast-height-band cell.
  - canopy mask: tree_id of the highest point in each cell (matches the CHM
    feature, which stores that point's height). Note this differs from the
    per-tree canopy masks, which mark ALL of a tree's branch cells including
    understory ones hidden under a taller crown.
(mask > 0) recovers a binary semantic mask.

Requires: numpy, pandas, scipy, rasterio, laspy (via las.py)
"""

from __future__ import annotations

import argparse
import json
import os
import shutil

import numpy as np
import pandas as pd
import rasterio
from rasterio.transform import from_origin

import georef
import layout
from las import read_points_las, read_canonical_las
from stem_map import BREAST_HEIGHT_M, STEM_SECTION, STEM_CLASS
from canopy_map import BRANCH_CLASS, canopy_hull
import schema


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def _move(src, dst):
    """Move src -> dst if src exists, creating parents. Returns True if moved.

    A legacy .las source keeps its extension at a .laz destination (the
    bytes are uncompressed LAS; only the name would say otherwise), so
    pre-LAZ plots keep their on-disk format through the move.
    """
    if not os.path.exists(src):
        return False
    if (src.endswith(".las") and dst.endswith(".laz")
            and not os.path.exists(dst)):
        dst = dst[:-4] + ".las"
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    shutil.move(src, dst)
    return True


def _read_raster(path):
    """Read a single-band GeoTIFF's array with rasterio."""
    with rasterio.open(path) as ds:
        return ds.read(1)


def _grid_from_tif(path, anchor=None):
    """(height, width, cell, x_min, y_max) of a GeoTIFF's grid.

    Values are in the LOCAL METRE frame: an anchored raster's EPSG:4326
    transform is inverted through the anchor first, because the mask
    helpers overlay metre-space skeleton/cloud geometry."""
    with rasterio.open(path) as ds:
        h, w = ds.height, ds.width
        t = ds.transform
        if anchor is not None and ds.crs is not None:
            x_min, y_max, cell = georef.local_grid_from_geo(anchor, t)
        else:
            cell = float(t.a)                       # pixel width
            x_min, y_max = float(t.c), float(t.f)   # c=x of UL, f=y of UL
        return h, w, cell, x_min, y_max


def _to_pixel(x, y, x_min, y_max, cell, w, h):
    col = np.floor((np.asarray(x) - x_min) / cell).astype(int)
    row = np.floor((y_max - np.asarray(y)) / cell).astype(int)
    np.clip(col, 0, w - 1, out=col)
    np.clip(row, 0, h - 1, out=row)
    return col, row


def _tree_uids(plot_dir):
    """UIDs of generated trees, from the flat *_raw_skeleton.json filenames."""
    d = os.path.join(plot_dir, layout.INDIVIDUAL_TREES)
    uids = []
    for fn in os.listdir(d):
        if fn.endswith("_raw_skeleton.json") and "_tree_" not in fn:
            uids.append(fn[: -len("_raw_skeleton.json")])
    return sorted(uids)


# ---------------------------------------------------------------------------
# plot-level instance masks
# ---------------------------------------------------------------------------

def _write_instance_geotiff(path, array, x_min, y_max, cell, anchor=None):
    """uint16 single-band GeoTIFF, LZW, geotransform embedded (rasterio).

    EPSG:4326 through the anchor when the plot is georeferenced (same
    degree-cell transform as the aligned feature raster; spec Section 10),
    else the local metre frame with no CRS."""
    array = np.ascontiguousarray(array)
    if anchor is not None:
        transform, crs = georef.geo_transform(anchor, x_min, y_max, cell)
    else:
        transform, crs = from_origin(x_min, y_max, cell, cell), None
    with rasterio.open(
        path, "w", driver="GTiff",
        height=array.shape[0], width=array.shape[1],
        count=1, dtype=array.dtype, transform=transform, crs=crs,
        compress="lzw",
    ) as dst:
        dst.write(array, 1)


def _stem_instance_mask(feature_tif, pts, dbh_height, threshold,
                        trees_dir=None, plot_uid=None, anchor=None):
    """tree_id per cell of the filled stem cross-section, from skeleton geometry.

    The stem geometry is ground truth from each tree's skeleton (centerline +
    taper), so every cell inside a stem's disc (union across the band) is labelled
    with that tree_id -- a solid disc, not a hollow ring, and present even where
    the scan left no points. Where two trees' discs overlap a cell, the cell is
    assigned to the NEARER stem centre.
    """
    from stem_map import _load_skeletons, stem_discs_in_band, RasterGrid
    h, w, cell, x0, y1 = _grid_from_tif(feature_tif, anchor)
    feat = _read_raster(feature_tif)
    z_lo, z_hi = dbh_height - threshold, dbh_height + threshold

    lab = np.zeros((h, w), dtype=np.uint32)
    best_d2 = np.full((h, w), np.inf)                 # nearest-centre arbitration

    skeletons = _load_skeletons(trees_dir, plot_uid) if trees_dir else {}
    grid = RasterGrid(x0, y1 - h * cell, x0 + w * cell, y1, cell)
    for t, skel in skeletons.items():
        discs = stem_discs_in_band(skel, z_lo, z_hi)
        if not discs:
            continue
        tmask = grid.rasterize_discs(discs, value=1).astype(bool)
        if not tmask.any():
            continue
        rows, cols = np.where(tmask)
        # distance from each set cell centre to the nearest of this tree's discs
        cxw = x0 + (cols + 0.5) * cell
        cyw = y1 - (rows + 0.5) * cell
        d2 = np.full(len(rows), np.inf)
        for cx, cy, r in discs:
            d2 = np.minimum(d2, (cxw - cx) ** 2 + (cyw - cy) ** 2)
        take = d2 < best_d2[rows, cols]
        rr, cc = rows[take], cols[take]
        lab[rr, cc] = int(t)
        best_d2[rr, cc] = d2[take]

    # NOTE: the label is the full skeleton disc; it is deliberately NOT clipped
    # to the point-count feature raster, so the stem cross-section stays a solid
    # filled circle even where the scan produced no returns.
    return lab, (x0, y1, cell)


def _canopy_instance_mask(feature_tif, pts, anchor=None):
    """tree_id of the highest point in each cell (matches the CHM feature).

    The CHM feature is built with a zero-initialised max, so a cell whose only
    points have z <= 0 (possible from the near-ground foliage z-jitter) stays 0
    and counts as empty. Mirror that here: after labelling by highest point,
    clear any cell the feature treats as empty, so the label and feature agree
    cell-for-cell.
    """
    h, w, cell, x0, y1 = _grid_from_tif(feature_tif, anchor)
    feat = _read_raster(feature_tif)
    col, row = _to_pixel(pts[:, 0], pts[:, 1], x0, y1, cell, w, h)
    lab = np.zeros((h, w), dtype=np.uint32)
    order = np.argsort(pts[:, 2])                  # ascending z -> highest wins
    lab[row[order], col[order]] = pts[order, 4].astype(np.uint32)
    lab[feat <= 0] = 0                             # keep label == feature support
    return lab, (x0, y1, cell)


# ---------------------------------------------------------------------------
# per-tree volumes
# ---------------------------------------------------------------------------

def _per_tree_volumes(raw_pts, trimmed_pts):
    """{tree_id: canopy hull volume} of convex-hull canopy volume.

    Canopy = foliage points (semantic 3). The trimmed cloud is already the
    culled geometry (no flag, no record of removed points, spec 3.2), so
    it is preferred when present.
    """
    out = {}
    pts = trimmed_pts if trimmed_pts is not None else raw_pts
    sem = pts[:, 3].astype(int) == schema.SEM_FOLIAGE
    for t in np.unique(pts[sem, 4].astype(int)):
        m = sem & (pts[:, 4].astype(int) == t)
        _, _, vol, _ = canopy_hull(pts[m, :3])
        out[int(t)] = float(vol)
    return out


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def restructure_plot(csv_path, plot_dir=None, verbose=True):
    plot_uid = os.path.splitext(os.path.basename(csv_path))[0]
    if plot_dir is None:
        plot_dir = os.path.join(os.path.dirname(os.path.abspath(csv_path)), plot_uid)
    flat_it = os.path.join(plot_dir, layout.INDIVIDUAL_TREES)
    flat_pp = os.path.join(plot_dir, layout.PLOT_PLACEMENT)

    uids = _tree_uids(plot_dir)

    # If the flat placement cloud is gone, this plot was already restructured.
    raw_las = layout.cloud_file(flat_pp, f"{plot_uid}_raw")
    if not os.path.exists(raw_las):
        if os.path.exists(layout.plot_path(plot_dir, plot_uid, "raw_las", create=False)):
            if verbose:
                print(f"{plot_uid}: already restructured, skipping")
            return _load(layout.plot_path(plot_dir, plot_uid, "scene", create=False))
        raise FileNotFoundError(f"placement cloud not found: {raw_las}")

    # ---- read clouds & combined JSONs before moving anything --------------
    # canonical reads -> positional (x,y,z,semantic_code,tree_id) arrays for the
    # mask/volume helpers below (which use column indices 3=semantic, 4=tree_id).
    trimmed_las = layout.cloud_file(flat_pp, f"{plot_uid}_trimmed")
    rxyz, rcols = read_canonical_las(raw_las, level="plot")
    raw_pts = np.column_stack([rxyz, rcols["semantic_code"].astype(float),
                               rcols["tree_id"].astype(float)])
    if os.path.exists(trimmed_las):
        txyz, tcols = read_canonical_las(trimmed_las, level="plot")
        trimmed_pts = np.column_stack([txyz, tcols["semantic_code"].astype(float),
                                       tcols["tree_id"].astype(float)])
    else:
        trimmed_pts = None
    has_trim = trimmed_pts is not None

    stem_feat_flat = os.path.join(flat_pp, f"{plot_uid}_raw_placement_stem_map.tif")
    canopy_feat_flat = os.path.join(flat_pp, f"{plot_uid}_raw_placement_canopy_map.tif")
    stem_bbox_doc = _load(os.path.join(
        flat_pp, f"{plot_uid}_raw_placement_stem_map_bounding_boxes.json"))
    canopy_bbox_doc = _load(os.path.join(
        flat_pp, f"{plot_uid}_raw_placement_canopy_bounding_box.json"))
    canopy_hull_doc = _load(os.path.join(
        flat_pp, f"{plot_uid}_raw_placement_canopy_convex_hull.json"))
    scene = _load(os.path.join(flat_pp, f"{plot_uid}_scene.json")) or {}

    canopy_bbox_by_tree = {b["tree_id"]: b for b in
                           (canopy_bbox_doc or {}).get("boxes", [])}
    canopy_hull_by_tree = {h["tree_id"]: h for h in
                           (canopy_hull_doc or {}).get("hulls", [])}

    vols = _per_tree_volumes(raw_pts, trimmed_pts if has_trim else None)

    # ---- plot-level instance masks (aligned to their feature rasters) -----
    anchor = georef.anchor_from_scene(scene)
    if os.path.exists(stem_feat_flat):
        lab, (x0, y1, cell) = _stem_instance_mask(
            stem_feat_flat, raw_pts, BREAST_HEIGHT_M, STEM_SECTION,
            trees_dir=flat_it, plot_uid=plot_uid, anchor=anchor)
        _write_instance_geotiff(
            layout.plot_path(plot_dir, plot_uid, "stem_label_raw"),
            lab, x0, y1, cell, anchor=anchor)
    if os.path.exists(canopy_feat_flat):
        lab, (x0, y1, cell) = _canopy_instance_mask(canopy_feat_flat, raw_pts,
                                                     anchor=anchor)
        _write_instance_geotiff(
            layout.plot_path(plot_dir, plot_uid, "canopy_label_raw"),
            lab, x0, y1, cell, anchor=anchor)

    # ---- move plot-level files into place ---------------------------------
    _move(os.path.join(plot_dir, "manifest.csv"),
          layout.plot_path(plot_dir, plot_uid, "manifest"))
    _move(raw_las, layout.plot_path(plot_dir, plot_uid, "raw_las"))
    _move(trimmed_las, layout.plot_path(plot_dir, plot_uid, "trimmed_las"))
    _move(stem_feat_flat, layout.plot_path(plot_dir, plot_uid, "stem_feature_raw"))
    _move(canopy_feat_flat,
          layout.plot_path(plot_dir, plot_uid, "canopy_feature_raw"))

    # ---- plot-level canopy total ------------------------------------------
    canopy_total = round(sum(vols.values()), 4)
    plot_vol = {
        "plot_uid": plot_uid,
        "n_trees": len(vols),
        "method": "sum of per-tree convex hull volumes of branch points",
        "total_canopy_volume_m3": canopy_total,
    }
    # (plot volume JSON retired; superseded by build_info's _info.json)

    # ---- per-tree: move + derive -----------------------------------------
    for uid in uids:
        tid = _uid_to_id(uid)

        # textual
        _move(os.path.join(flat_it, f"{uid}_generation.json"),
              layout.tree_path(plot_dir, plot_uid, uid, "generation"))
        _move(os.path.join(flat_it, f"{uid}_raw_skeleton.json"),
              layout.tree_path(plot_dir, plot_uid, uid, "raw_skeleton"))
        _move(os.path.join(flat_it, f"{uid}_trimmed_skeleton.json"),
              layout.tree_path(plot_dir, plot_uid, uid, "trimmed_skeleton"))
        # las
        _move(layout.cloud_file(flat_it, f"{uid}_raw"),
              layout.tree_path(plot_dir, plot_uid, uid, "raw_las"))
        _move(layout.cloud_file(flat_it, f"{uid}_trimmed"),
              layout.tree_path(plot_dir, plot_uid, uid, "trimmed_las"))
        # stem map crop + mask
        _move(os.path.join(flat_it, f"{plot_uid}_tree_{uid}_stem_map_bounding_box.tif"),
              layout.tree_path(plot_dir, plot_uid, uid, "stem_feature_raw"))
        _move(os.path.join(flat_it, f"{plot_uid}_tree_{uid}_stem_map_semantic_mask.tif"),
              layout.tree_path(plot_dir, plot_uid, uid, "stem_label_raw"))
        # canopy map crop + mask
        _move(os.path.join(flat_it, f"{plot_uid}_tree_{uid}_canopy_bounding_box.tif"),
              layout.tree_path(plot_dir, plot_uid, uid, "canopy_feature_raw"))
        _move(os.path.join(flat_it, f"{plot_uid}_tree_{uid}_canopy_semantic_mask.tif"),
              layout.tree_path(plot_dir, plot_uid, uid, "canopy_label_raw"))

        # per-tree canopy bbox / hull JSONs (split from combined)
        if tid in canopy_bbox_by_tree:
            _dump(layout.tree_path(plot_dir, plot_uid, uid, "canopy_bbox_json"),
                  {"plot_uid": plot_uid, "tree_id": tid, "tree_uid": uid,
                   **canopy_bbox_by_tree[tid]})
        if tid in canopy_hull_by_tree:
            _dump(layout.tree_path(plot_dir, plot_uid, uid, "canopy_hull_json"),
                  {"plot_uid": plot_uid, "tree_id": tid, "tree_uid": uid,
                   "coordinate_frame": "plot metres, z = 0 at ground",
                   **canopy_hull_by_tree[tid]})

        # (per-tree volume JSON retired; superseded by build_info's _info.json)

    # ---- rewrite scene JSON with new links; fold in stem bboxes -----------
    scene["structure"] = "canonical (see layout.py)"
    scene["plot_files"] = {
        "manifest": layout.plot_rel(plot_dir, plot_uid, "manifest"),
        "raw_las": layout.plot_rel(plot_dir, plot_uid, "raw_las"),
        "trimmed_las": layout.plot_rel(plot_dir, plot_uid, "trimmed_las"),
        # (no "volume" entry: the plot-level volume JSON is retired,
        # superseded by build_info's per-tree _info.json; a link here
        # would dangle, violating the scene's every-link-resolves role)
        "stem_map_feature": layout.plot_rel(plot_dir, plot_uid, "stem_feature_raw"),
        "stem_map_label": layout.plot_rel(plot_dir, plot_uid, "stem_label_raw"),
        "canopy_map_feature": layout.plot_rel(plot_dir, plot_uid, "canopy_feature_raw"),
        "canopy_map_label": layout.plot_rel(plot_dir, plot_uid, "canopy_label_raw"),
    }
    # the stage records above point at the FLAT-layout clouds as they were
    # when the stage ran; restructure has just MOVED those files, so rewrite
    # the provenance to the canonical paths (same files, new location)
    if scene.get("stem_map") is not None and stem_bbox_doc is not None:
        scene["stem_map"]["bounding_boxes_inline"] = stem_bbox_doc.get("boxes", [])
        scene["stem_map"].pop("bounding_boxes", None)
        scene["stem_map"]["feature"] = scene["plot_files"]["stem_map_feature"]
        scene["stem_map"]["label"] = scene["plot_files"]["stem_map_label"]
        scene["stem_map"]["source_cloud"] = scene["plot_files"]["raw_las"]
        scene["stem_map"].pop("raster", None)
    if scene.get("canopy_map") is not None:
        scene["canopy_map"]["feature"] = scene["plot_files"]["canopy_map_feature"]
        scene["canopy_map"]["label"] = scene["plot_files"]["canopy_map_label"]
        scene["canopy_map"]["source_cloud"] = scene["plot_files"]["trimmed_las"]
        scene["canopy_map"].pop("raster", None)
        # per-tree geometry now lives beside each tree
        scene["canopy_map"].pop("bounding_box", None)
        scene["canopy_map"].pop("convex_hull", None)
    def _tree_links(uid):
        keys = [("raw_las", "raw_las"), ("trimmed_las", "trimmed_las"),
                ("volume", "volume_label"),
                ("stem_feature", "stem_feature_raw"),
                ("stem_label", "stem_label_raw"),
                ("canopy_feature", "canopy_feature_raw"),
                ("canopy_label", "canopy_label_raw"),
                ("canopy_bbox", "canopy_bbox_json"),
                ("canopy_hull", "canopy_hull_json")]
        out = {}
        for name, key in keys:
            p = layout.tree_path(plot_dir, plot_uid, uid, key, create=False)
            if os.path.exists(p):
                out[name] = layout.tree_rel(plot_dir, plot_uid, uid, key)
        return out

    scene["trees"] = {uid: _tree_links(uid) for uid in uids}
    _dump(layout.plot_path(plot_dir, plot_uid, "scene"), scene)

    # ---- clean up consumed flat files -------------------------------------
    for leftover in [
        os.path.join(flat_pp, f"{plot_uid}_scene.json"),
        os.path.join(flat_pp, f"{plot_uid}_raw_placement_stem_map_bounding_boxes.json"),
        os.path.join(flat_pp, f"{plot_uid}_raw_placement_canopy_bounding_box.json"),
        os.path.join(flat_pp, f"{plot_uid}_raw_placement_canopy_convex_hull.json"),
    ]:
        if os.path.exists(leftover):
            os.remove(leftover)

    if verbose:
        print(f"{plot_uid}: restructured {len(uids)} trees into the canonical "
              f"layout (canopy volume {canopy_total:.1f} m3)")
    return scene


def _load(path):
    if path and os.path.exists(path):
        with open(path) as f:
            return json.load(f)
    return None


def _dump(path, obj):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        json.dump(obj, f, indent=2)


def _uid_to_id(uid):
    try:
        return int(str(uid).lstrip("0") or 0)
    except ValueError:
        return uid


def main():
    ap = argparse.ArgumentParser(
        description="Reorganise flat plot outputs into the canonical layout")
    ap.add_argument("csv", nargs="+", help="plot table csv(s)")
    ap.add_argument("-d", "--plot-dir", default=None)
    args = ap.parse_args()
    for path in args.csv:
        restructure_plot(path, plot_dir=args.plot_dir)


if __name__ == "__main__":
    main()
