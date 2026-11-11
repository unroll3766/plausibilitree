"""
canopy_map.py: rasterise a canopy height model (CHM) of a placed plot, with
per-tree 3D canopy bounding boxes, convex hulls, crops, masks and a total
canopy volume.

    python canopy_map.py plot_11519.csv --cell-size 0.01

The CHM is a top-down (x-y) raster whose value in each cell is the height of the
highest point falling in it (0 where no point falls, i.e. no tree). Built from
the placed scene (plot_placement/<PLOT_UID>_raw.las).

Outputs, under the plot directory:

    plot_placement/<PLOT_UID>_raw_placement_canopy_map.tif
        single-band float32 height raster (GeoTIFF), whole plot
    plot_placement/<PLOT_UID>_raw_placement_canopy_bounding_box.json
        per-tree 3D canopy bounding box (rectangular prism of branch points)
    plot_placement/<PLOT_UID>_raw_placement_canopy_convex_hull.json
        per-tree 3D convex hull of the branch points
    individual_trees/<PLOT_UID>_tree_<UID>_canopy_bounding_box.tif
        the CHM cropped to that tree's canopy x-y bbox plus a buffer
    individual_trees/<PLOT_UID>_tree_<UID>_canopy_semantic_mask.tif
        binary mask, same size as the crop: 1 where a cell holds >=1 branch
        point of that tree (using ALL of that tree's branch points -- so an
        understory tree's mask covers cells whose CHM height belongs to a
        taller tree above it)

The plot scene JSON (plot_placement/<PLOT_UID>_scene.json) gains a "canopy_map"
block with relative links and the total canopy volume (raw, and trimmed if the
trimmed cloud from the assembly step is available).

Requires: numpy, pandas, scipy, rasterio, laspy (via las.py)
Shares the raster grid + GeoTIFF writer with stem_map.py.
"""

from __future__ import annotations

import argparse
import json
import os

import numpy as np
import pandas as pd

import georef
from las import read_points_las, read_canonical_las
import layout
import schema
from stem_map import RasterGrid, _write_geotiff, _tree_uid

BRANCH_CLASS = schema.SEM_FOLIAGE    # canopy = foliage points (semantic 3)
CHM_RESOLUTION = 0.05       # metres per pixel (spec CHM_RESOLUTION)
DEFAULT_BUFFER_PX = 2        # crop buffer around a canopy bbox, pixels


# ---------------------------------------------------------------------------
# Convex hull
# ---------------------------------------------------------------------------

def canopy_hull(points):
    """Convex hull of a set of 3D points.

    Returns (vertices, faces, volume, area) with faces indexing into vertices,
    or (None, None, 0.0, 0.0) if a hull can't be formed (fewer than 4 points or
    a degenerate/coplanar set).
    """
    from scipy.spatial import ConvexHull
    try:
        from scipy.spatial import QhullError
    except ImportError:                              # older scipy
        from scipy.spatial.qhull import QhullError

    pts = np.asarray(points, dtype=np.float64)
    if len(pts) < 4:
        return None, None, 0.0, 0.0
    try:
        hull = ConvexHull(pts)
    except QhullError:
        return None, None, 0.0, 0.0

    used = np.unique(hull.simplices)
    remap = {int(o): i for i, o in enumerate(used)}
    verts = pts[used]

    # scipy's simplices are not consistently wound; orient each triangle so its
    # normal points outward (agrees with hull.equations), so the stored mesh is
    # directly usable (e.g. signed-volume / rendering) without re-orienting.
    faces = []
    for tri, eq in zip(hull.simplices, hull.equations):
        a, b, c = pts[tri]
        n = np.cross(b - a, c - a)
        if np.dot(n, eq[:3]) < 0:
            tri = tri[::-1]
        faces.append([remap[int(v)] for v in tri])
    faces = np.array(faces, dtype=int)
    return verts, faces, float(hull.volume), float(hull.area)


# ---------------------------------------------------------------------------
# CHM rasterisation
# ---------------------------------------------------------------------------

def max_height_raster(grid: RasterGrid, x, y, z):
    """Highest z per cell into a (height, width) float32 raster, 0 elsewhere."""
    col, row = grid.to_pixel(x, y)
    flat_idx = row * grid.width + col
    out = np.zeros(grid.height * grid.width, dtype=np.float64)
    # np.maximum.at accumulates the per-cell maximum over duplicate indices
    np.maximum.at(out, flat_idx, np.asarray(z, dtype=np.float64))
    return out.reshape(grid.height, grid.width).astype(np.float32)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def build_canopy_map(csv_path, plot_dir=None,
                     cell_size=CHM_RESOLUTION,
                     buffer_px=DEFAULT_BUFFER_PX,
                     source="raw",
                     verbose=True):
    """Build the CHM, canopy bboxes, hulls, per-tree crops/masks and volume."""
    plot_uid = os.path.splitext(os.path.basename(csv_path))[0]
    if plot_dir is None:
        plot_dir = os.path.join(os.path.dirname(os.path.abspath(csv_path)), plot_uid)
    place_dir = os.path.join(plot_dir, "plot_placement")
    trees_dir = os.path.join(plot_dir, "individual_trees")
    anchor = georef.load_anchor(plot_dir, plot_uid)
    os.makedirs(trees_dir, exist_ok=True)

    las_path = layout.cloud_file(place_dir, f"{plot_uid}_{source}")
    if not os.path.exists(las_path):
        raise FileNotFoundError(f"placement cloud not found: {las_path}")

    xyz, cols = read_canonical_las(las_path, level="plot")
    sem = cols["semantic_code"].astype(int)
    tid = cols["tree_id"].astype(int)

    # CHM over the extent of all points; values from vegetation only
    # (codes 1, 2, 3, 8 -- grass 9 is ground cover and never canopy)
    grid = RasterGrid(xyz[:, 0].min(), xyz[:, 1].min(),
                      xyz[:, 0].max(), xyz[:, 1].max(), cell_size)
    veg = np.isin(sem, (schema.SEM_TREE, schema.SEM_BRANCH,
                        schema.SEM_FOLIAGE, schema.SEM_SHRUB))
    chm = max_height_raster(grid, xyz[veg, 0], xyz[veg, 1], xyz[veg, 2])
    # band 2: foliage-presence mask (1 where a code-3 point lies in the cell)
    fol = sem == schema.SEM_FOLIAGE
    foliage_mask = (grid.counts(xyz[fol, 0], xyz[fol, 1]) > 0).astype(np.float32)

    raster_name = f"{plot_uid}_raw_placement_canopy_map.tif"
    _write_geotiff(os.path.join(place_dir, raster_name),
                   [chm, foliage_mask],
                   grid.x_min, grid.y_max, grid.cell, anchor=anchor)

    # per-tree canopy geometry from foliage points (semantic 3) -- the leafy
    # envelope; branch stems (2) are inside it and excluded by design.
    branch = sem == schema.SEM_FOLIAGE
    buf = int(buffer_px)
    bbox_entries, hull_entries = [], []
    per_tree_px = {}          # tree_id -> buffered (col0,row0,col1,row1)
    raw_volume_total = 0.0

    for t in np.unique(tid[branch]):
        m = branch & (tid == t)
        cx, cy, cz = xyz[m, 0], xyz[m, 1], xyz[m, 2]

        # 3D bounding box (rectangular prism)
        x0, x1 = float(cx.min()), float(cx.max())
        y0, y1 = float(cy.min()), float(cy.max())
        z0, z1 = float(cz.min()), float(cz.max())

        col, row = grid.to_pixel(cx, cy)
        c0, c1 = int(col.min()), int(col.max())
        r0, r1 = int(row.min()), int(row.max())
        bc0 = max(0, c0 - buf)
        br0 = max(0, r0 - buf)
        bc1 = min(grid.width - 1, c1 + buf)
        br1 = min(grid.height - 1, r1 + buf)
        per_tree_px[int(t)] = (bc0, br0, bc1, br1)

        # round outward (floor mins, ceil maxes) so the stored box always
        # encloses every branch point despite 4-dp rounding
        def _fl(v):
            return float(np.floor(v * 1e4) / 1e4)

        def _ce(v):
            return float(np.ceil(v * 1e4) / 1e4)

        bbox_entries.append({
            "tree_id": int(t),
            "n_branch_points": int(m.sum()),
            "bbox_world": {
                "x_min": _fl(x0), "x_max": _ce(x1),
                "y_min": _fl(y0), "y_max": _ce(y1),
                "z_min": _fl(z0), "z_max": _ce(z1),
                "dx": round(x1 - x0, 4), "dy": round(y1 - y0, 4),
                "dz": round(z1 - z0, 4),
            },
            "bbox_volume_m3": round((x1 - x0) * (y1 - y0) * (z1 - z0), 6),
            # buffered x-y pixel bbox used for the crop
            "bbox_pixel_buffered": {"col_min": bc0, "row_min": br0,
                                    "col_max": bc1, "row_max": br1,
                                    "width": bc1 - bc0 + 1,
                                    "height": br1 - br0 + 1,
                                    "buffer_px": buf},
        })

        verts, faces, vol, area = canopy_hull(np.column_stack([cx, cy, cz]))
        raw_volume_total += vol
        hull_entries.append({
            "tree_id": int(t),
            "n_branch_points": int(m.sum()),
            "volume_m3": round(vol, 6),
            "surface_area_m2": round(area, 6),
            "n_vertices": 0 if verts is None else int(len(verts)),
            "n_faces": 0 if faces is None else int(len(faces)),
            "vertices": [] if verts is None else
                        [[round(float(v[0]), 4), round(float(v[1]), 4),
                          round(float(v[2]), 4)] for v in verts],
            "faces": [] if faces is None else faces.tolist(),
            "degenerate": verts is None,
        })

    bbox_entries.sort(key=lambda b: b["tree_id"])
    hull_entries.sort(key=lambda h: h["tree_id"])

    # write the two geometry JSONs
    bbox_name = f"{plot_uid}_raw_placement_canopy_bounding_box.json"
    hull_name = f"{plot_uid}_raw_placement_canopy_convex_hull.json"
    common = {
        "plot_uid": plot_uid,
        "raster": raster_name,
        "cell_size_m": float(cell_size),
        "n_trees_with_canopy": len(bbox_entries),
    }
    with open(os.path.join(place_dir, bbox_name), "w") as f:
        json.dump({**common, "buffer_px": buf, "boxes": bbox_entries}, f, indent=2)
    with open(os.path.join(place_dir, hull_name), "w") as f:
        json.dump({**common,
                   "coordinate_frame": "plot metres, z = 0 at ground",
                   "hulls": hull_entries}, f, indent=2)

    # per-tree crops + canopy masks
    n_written = 0
    for t, (bc0, br0, bc1, br1) in per_tree_px.items():
        uid = _tree_uid(t, csv_path)
        crop = chm[br0:br1 + 1, bc0:bc1 + 1]

        # mask: every cell that holds >=1 branch point of THIS tree (all of its
        # branch points, not just the topmost -- see the understory note)
        m = branch & (tid == t)
        col, row = grid.to_pixel(xyz[m, 0], xyz[m, 1])
        mask = np.zeros_like(crop, dtype=np.uint8)
        mask[row - br0, col - bc0] = 1

        x0 = grid.x_min + bc0 * grid.cell
        y0 = grid.y_max - br0 * grid.cell
        _write_geotiff(os.path.join(
            trees_dir, f"{plot_uid}_tree_{uid}_canopy_bounding_box.tif"),
            crop, x0, y0, grid.cell, anchor=anchor)
        _write_geotiff(os.path.join(
            trees_dir, f"{plot_uid}_tree_{uid}_canopy_semantic_mask.tif"),
            mask, x0, y0, grid.cell, anchor=anchor)
        n_written += 1

    # single post-trim canopy volume (spec 3.2: no record of trimmed
    # points is kept, so there is no raw/trimmed pair)
    volume = {"volume_m3": round(raw_volume_total, 4),
              "method": "sum of per-tree convex hull volumes"}

    # link from the plot scene JSON
    rel = f"{plot_uid}/plot_placement"
    canopy_meta = {
        "raster": f"{rel}/{raster_name}",
        "bounding_box": f"{rel}/{bbox_name}",
        "convex_hull": f"{rel}/{hull_name}",
        "cell_size_m": float(cell_size),
        "buffer_px": buf,
        "source_cloud": f"{rel}/{os.path.basename(las_path)}",
        "raster_size": {"width": grid.width, "height": grid.height},
        "n_trees_with_canopy": len(bbox_entries),
        "total_canopy_volume": volume,
    }
    scene_path = os.path.join(place_dir, f"{plot_uid}_scene.json")
    if os.path.exists(scene_path):
        scene = json.load(open(scene_path))
        scene["canopy_map"] = canopy_meta
        with open(scene_path, "w") as f:
            json.dump(scene, f, indent=2)

    if verbose:
        print(f"{plot_uid}: CHM {grid.width}x{grid.height} px @ "
              f"{cell_size*100:g} cm, max height {chm.max():.1f} m, "
              f"{len(bbox_entries)} canopies, {n_written} crops+masks, "
              f"{volume['volume_m3']:.1f} m3 canopy volume")
    return canopy_meta


def main():
    ap = argparse.ArgumentParser(description="Canopy height map + 3D canopy geometry")
    ap.add_argument("csv", nargs="+", help="plot table csv(s)")
    ap.add_argument("-d", "--plot-dir", default=None)
    ap.add_argument("--chm-resolution", type=float, default=CHM_RESOLUTION,
                    help=f"metres per pixel (default {CHM_RESOLUTION})")
    ap.add_argument("--buffer", type=int, default=DEFAULT_BUFFER_PX,
                    help="crop buffer around a canopy bbox, pixels (default 2)")
    ap.add_argument("--source", default="raw", choices=["raw", "trimmed"])
    args = ap.parse_args()

    for path in args.csv:
        build_canopy_map(path, plot_dir=args.plot_dir, cell_size=args.chm_resolution,
                         buffer_px=args.buffer, source=args.source)


if __name__ == "__main__":
    main()
