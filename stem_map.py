"""
stem_map.py: rasterise a breast-height cross-section of a placed plot into a
stem map, with per-tree bounding boxes and semantic masks for deep learning.

    python stem_map.py plot_11519.csv --cell-size 0.01

A stem map is a top-down (x-y) raster of the point cloud in a thin horizontal
band at DBH (breast height, 1.37 m) +/- a threshold. Each cell's value is the
number of points falling in it -- a stem shows up as a dense ring/disc of
returns. Produced from the placed scene (plot_placement/<PLOT_UID>_raw.las).

Outputs, under the plot directory:

    plot_placement/<PLOT_UID>_raw_placement_stem_map.tif
        single-band count raster (GeoTIFF), whole plot
    plot_placement/<PLOT_UID>_raw_placement_stem_map_bounding_boxes.json
        per-tree stem bounding boxes (pixel + world) and tree ids
    individual_trees/<PLOT_UID>_tree_<UID>_stem_map_bounding_box.tif
        the stem map cropped to that tree's stem bbox
    individual_trees/<PLOT_UID>_tree_<UID>_stem_map_semantic_mask.tif
        binary mask, same size as the crop: 1 where a cell holds >=1 stem
        point of that tree

Relative links to the raster and the bounding-box JSON are added to the plot
scene JSON (plot_placement/<PLOT_UID>_scene.json) under "stem_map".

GeoTIFFs are written with rasterio (embedding the geotransform, so pixel<->world
is recoverable).

Requires: numpy, pandas, rasterio, laspy (via las.py)
"""

from __future__ import annotations

import argparse
import json
import math
import os

import numpy as np
import pandas as pd

from las import read_points_las, read_canonical_las
import layout
import schema


def _taper_radius(base_radius, total_height, h):
    """Trunk radius at height h under the linear taper d(h)=d_b*(1-h/h_t).

    Clamped to >= 0 (above the modelled tip the radius is 0).
    """
    if total_height <= 0:
        return 0.0
    return max(0.0, base_radius * (1.0 - h / total_height))


def stem_discs_in_band(skel, z_lo, z_hi, n_samples=9):
    """Ground-truth stem cross-section discs across a breast-height band.

    Uses the tree skeleton (not the point cloud): for every trunk polyline that
    passes through [z_lo, z_hi], sample the centerline at heights spanning the
    band and, at each, take the taper radius. Returns a list of (cx, cy, r) in
    plot metres -- the union of these discs is the stem mask for the band.

    A band of zero width (z_lo == z_hi) yields exactly one disc per polyline at
    that height (the single-slice case). Forked trees contribute one disc per
    fork polyline present in the band, so the general multi-stem case is handled.
    Radii come from the taper model, so a stem with no scanned points in the band
    still produces a correct filled disc.
    """
    base_r = float(skel.get("base_radius_m", 0.0))
    h_t = float(skel.get("height_m", 0.0))
    discs = []
    heights = (np.array([z_lo]) if z_hi <= z_lo
               else np.linspace(z_lo, z_hi, n_samples))
    for pl in skel.get("trunk", {}).get("polylines", []):
        arr = np.asarray(pl, dtype=float)          # (N,3) x,y,z
        if len(arr) < 2:
            continue
        pz = arr[:, 2]
        lo, hi = pz.min(), pz.max()
        for h in heights:
            if h < lo - 1e-9 or h > hi + 1e-9:
                continue                            # this polyline not present here
            cx = float(np.interp(h, pz, arr[:, 0]))
            cy = float(np.interp(h, pz, arr[:, 1]))
            r = _taper_radius(base_r, h_t, h)
            if r > 0:
                discs.append((cx, cy, r))
    return discs

BREAST_HEIGHT_M = 1.37       # DBH reference height (matches generate_trees.py)
DEFAULT_CELL_SIZE = 0.01     # metres per pixel
DEFAULT_THRESHOLD = 0.10     # half-thickness of the z-band, metres
STEM_CLASS = schema.SEM_TREE       # trunk; stem map also includes SEM_BRANCH (2)
STEM_CLASSES = (schema.SEM_TREE, schema.SEM_BRANCH)   # all woody stem


# ---------------------------------------------------------------------------
# GeoTIFF I/O
# ---------------------------------------------------------------------------

def _write_geotiff(path, array, x_min, y_max, cell_size):
    """Write a single-band north-up GeoTIFF with rasterio.

    The transform places pixel (0, 0)'s upper-left corner at (x_min, y_max)
    with y decreasing down the rows.
    """
    import rasterio
    from rasterio.transform import from_origin

    array = np.ascontiguousarray(array)
    transform = from_origin(x_min, y_max, cell_size, cell_size)
    with rasterio.open(
        path, "w", driver="GTiff",
        height=array.shape[0], width=array.shape[1],
        count=1, dtype=array.dtype, transform=transform,
        compress="lzw",
    ) as dst:
        dst.write(array, 1)


# ---------------------------------------------------------------------------
# Rasterisation
# ---------------------------------------------------------------------------

class RasterGrid:
    """North-up pixel grid over a world x-y extent."""

    def __init__(self, x_min, y_min, x_max, y_max, cell_size):
        self.cell = float(cell_size)
        self.x_min = float(x_min)
        self.y_max = float(y_max)
        self.width = max(1, int(math.ceil((x_max - x_min) / self.cell)))
        self.height = max(1, int(math.ceil((y_max - y_min) / self.cell)))

    def to_pixel(self, x, y):
        """World -> (col, row) integer arrays, clipped to the grid."""
        col = np.floor((np.asarray(x) - self.x_min) / self.cell).astype(int)
        row = np.floor((self.y_max - np.asarray(y)) / self.cell).astype(int)
        np.clip(col, 0, self.width - 1, out=col)
        np.clip(row, 0, self.height - 1, out=row)
        return col, row

    def counts(self, x, y):
        """Count points per cell into a (height, width) uint16 raster."""
        col, row = self.to_pixel(x, y)
        flat = np.bincount(row * self.width + col,
                           minlength=self.height * self.width)
        return flat.reshape(self.height, self.width).astype(np.uint16)

    def rasterize_discs(self, discs, out=None, value=1):
        """Fill each (cx, cy, r) disc into a (height, width) mask.

        A cell is set if its centre falls within r of a disc centre, so the
        result is a filled circle (interior included), not a ring. Discs smaller
        than a cell still set at least their centre cell, so a thin stem never
        vanishes. Returns the mask; pass `out` to accumulate into an existing one.
        """
        if out is None:
            out = np.zeros((self.height, self.width), dtype=np.uint16)
        # cell-centre world coordinates
        for cx, cy, r in discs:
            # candidate pixel window around the disc (+1 cell margin)
            rc = r + self.cell
            c_lo = max(0, int(np.floor((cx - rc - self.x_min) / self.cell)))
            c_hi = min(self.width - 1, int(np.ceil((cx + rc - self.x_min) / self.cell)))
            r_lo = max(0, int(np.floor((self.y_max - (cy + rc)) / self.cell)))
            r_hi = min(self.height - 1, int(np.ceil((self.y_max - (cy - rc)) / self.cell)))
            if c_lo > c_hi or r_lo > r_hi:
                continue
            cols = np.arange(c_lo, c_hi + 1)
            rows = np.arange(r_lo, r_hi + 1)
            cxw = self.x_min + (cols + 0.5) * self.cell
            cyw = self.y_max - (rows + 0.5) * self.cell
            dx = cxw[None, :] - cx
            dy = cyw[:, None] - cy
            inside = (dx * dx + dy * dy) <= (r * r)
            # ensure a sub-cell disc still marks its centre cell
            ccol = int(np.clip(np.floor((cx - self.x_min) / self.cell), 0, self.width - 1))
            crow = int(np.clip(np.floor((self.y_max - cy) / self.cell), 0, self.height - 1))
            out[crow, ccol] = value
            block = out[r_lo:r_hi + 1, c_lo:c_hi + 1]
            block[inside] = value
        return out


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def build_stem_map(csv_path, plot_dir=None,
                   cell_size=DEFAULT_CELL_SIZE,
                   dbh_height=BREAST_HEIGHT_M,
                   threshold=DEFAULT_THRESHOLD,
                   bbox_padding=1,
                   source="raw",
                   verbose=True):
    """Build the plot stem map, bounding boxes and per-tree crops/masks.

    source selects the placement LAS: 'raw' -> <plot>_raw.las (default, matches
    the *_raw_placement_* output names), 'trimmed' -> <plot>_trimmed.las.

    bbox_padding grows each tree's stem bounding box by that many pixels on
    every side (default 1), clipped to the raster. The padded box is the one
    used for the crop and mask and stored as bbox_pixel; the tight box is kept
    as bbox_pixel_tight.
    """
    plot_uid = os.path.splitext(os.path.basename(csv_path))[0]
    if plot_dir is None:
        plot_dir = os.path.join(os.path.dirname(os.path.abspath(csv_path)), plot_uid)
    place_dir = os.path.join(plot_dir, "plot_placement")
    trees_dir = os.path.join(plot_dir, "individual_trees")
    os.makedirs(trees_dir, exist_ok=True)

    las_path = layout.cloud_file(place_dir, f"{plot_uid}_{source}")
    if not os.path.exists(las_path):
        raise FileNotFoundError(f"placement cloud not found: {las_path}")

    xyz, cols = read_canonical_las(las_path, level="plot")
    sem = cols["semantic_code"].astype(int)
    tid = cols["tree_id"].astype(int)

    # breast-height band. Band POINTS are only used for the secondary count
    # raster and an initial grid extent -- the stem MASK is skeleton-based, so an
    # empty band (e.g. threshold 0, or a sparsely-scanned slice) is fine; the
    # grid is then seeded from the whole cloud and widened to the discs below.
    z_lo, z_hi = dbh_height - threshold, dbh_height + threshold
    band = (xyz[:, 2] >= z_lo) & (xyz[:, 2] <= z_hi)
    if band.any():
        bx, by = xyz[band, 0], xyz[band, 1]
        bsem, btid = sem[band], tid[band]
        grid = RasterGrid(bx.min(), by.min(), bx.max(), by.max(), cell_size)
    else:
        # no points in the band: empty count raster over the plot's x-y extent,
        # discs will define the real footprint.
        bx = by = np.zeros(0)
        bsem = btid = np.zeros(0, dtype=int)
        grid = RasterGrid(xyz[:, 0].min(), xyz[:, 1].min(),
                          xyz[:, 0].max(), xyz[:, 1].max(), cell_size)
    stem_map = grid.counts(bx, by)                  # count of band points (may be 0)

    raster_name = f"{plot_uid}_raw_placement_stem_map.tif"
    raster_path = os.path.join(place_dir, raster_name)
    _write_geotiff(raster_path, stem_map, grid.x_min, grid.y_max, grid.cell)

    # per-tree stem geometry is ground truth from the skeleton (centerline +
    # taper), NOT the point cloud: the mask is the filled stem cross-section
    # across the band, so it's a solid disc (interior included) even where the
    # scan left a hollow ring or no points at all. Every tree with a skeleton
    # gets a stem, including sparsely-scanned ones.
    skeletons = _load_skeletons(trees_dir, plot_uid)
    stem_mask_band = np.isin(bsem, (schema.SEM_TREE, schema.SEM_BRANCH))
    pad = int(bbox_padding)
    boxes = []
    per_tree = {}       # tree_id -> padded (col_min, row_min, col_max, row_max)
    tree_discs = {}     # tree_id -> [(cx, cy, r), ...] across the band

    if not skeletons:
        raise RuntimeError(
            f"no tree skeletons found under {trees_dir}; the stem mask is built "
            f"from skeleton geometry, so assembly must have run first")

    # the grid must cover every disc's extent (discs can reach slightly past the
    # band points for sparse stems); widen it to include disc bounds.
    all_discs = {}
    gx_min, gy_min = grid.x_min, grid.y_max - grid.height * grid.cell
    gx_max, gy_max = grid.x_min + grid.width * grid.cell, grid.y_max
    for t, skel in skeletons.items():
        discs = stem_discs_in_band(skel, z_lo, z_hi)
        if not discs:
            continue
        all_discs[int(t)] = discs
        for cx, cy, r in discs:
            gx_min = min(gx_min, cx - r); gx_max = max(gx_max, cx + r)
            gy_min = min(gy_min, cy - r); gy_max = max(gy_max, cy + r)

    # rebuild the grid if any disc fell outside the band-point extent
    if (gx_min < grid.x_min or gy_max > grid.y_max
            or gx_max > grid.x_min + grid.width * grid.cell
            or gy_min < grid.y_max - grid.height * grid.cell):
        grid = RasterGrid(gx_min, gy_min, gx_max, gy_max, cell_size)
        stem_map = grid.counts(bx, by)
        _write_geotiff(raster_path, stem_map, grid.x_min, grid.y_max, grid.cell)

    for t in sorted(all_discs):
        discs = all_discs[int(t)]
        tree_discs[int(t)] = discs

        # tight pixel bbox = union of disc extents
        cxs = np.array([d[0] for d in discs]); cys = np.array([d[1] for d in discs])
        rs = np.array([d[2] for d in discs])
        col_lo, row_hi = grid.to_pixel(cxs - rs, cys - rs)  # min x -> min col; min y -> max row
        col_hi, row_lo = grid.to_pixel(cxs + rs, cys + rs)
        c0, c1 = int(col_lo.min()), int(col_hi.max())
        r0, r1 = int(row_lo.min()), int(row_hi.max())

        pc0 = max(0, c0 - pad)
        pr0 = max(0, r0 - pad)
        pc1 = min(grid.width - 1, c1 + pad)
        pr1 = min(grid.height - 1, r1 + pad)
        per_tree[int(t)] = (pc0, pr0, pc1, pr1)

        # band stem points for THIS tree, for reporting only (not the mask source)
        m_pts = stem_mask_band & (btid == t)

        boxes.append({
            "tree_id": int(t),
            "n_stem_points": int(m_pts.sum()),
            "n_discs": len(discs),
            "stem_source": "skeleton",
            "padding_px": pad,
            "bbox_pixel": {"col_min": pc0, "row_min": pr0,
                           "col_max": pc1, "row_max": pr1,
                           "width": pc1 - pc0 + 1, "height": pr1 - pr0 + 1},
            "bbox_pixel_tight": {"col_min": c0, "row_min": r0,
                                 "col_max": c1, "row_max": r1,
                                 "width": c1 - c0 + 1, "height": r1 - r0 + 1},
            "bbox_world": {
                "x_min": round(grid.x_min + pc0 * grid.cell, 4),
                "x_max": round(grid.x_min + (pc1 + 1) * grid.cell, 4),
                "y_min": round(grid.y_max - (pr1 + 1) * grid.cell, 4),
                "y_max": round(grid.y_max - pr0 * grid.cell, 4),
            },
        })
    boxes.sort(key=lambda b: b["tree_id"])


    bbox_name = f"{plot_uid}_raw_placement_stem_map_bounding_boxes.json"
    bbox_path = os.path.join(place_dir, bbox_name)
    bbox_doc = {
        "plot_uid": plot_uid,
        "raster": raster_name,
        "cell_size_m": float(cell_size),
        "dbh_height_m": float(dbh_height),
        "band_half_thickness_m": float(threshold),
        "bbox_padding_px": pad,
        "raster_size": {"width": grid.width, "height": grid.height},
        "raster_origin_world": {"x_min": round(grid.x_min, 4),
                                "y_max": round(grid.y_max, 4)},
        "raster_origin_world_exact": {"x_min": grid.x_min,
                                      "y_max": grid.y_max},
        "n_trees_with_stem": len(boxes),
        "boxes": boxes,
    }
    with open(bbox_path, "w") as f:
        json.dump(bbox_doc, f, indent=2)

    # per-tree crops + semantic masks
    n_written = 0
    for t, (c0, r0, c1, r1) in per_tree.items():
        uid = _tree_uid(t, csv_path)
        crop = stem_map[r0:r1 + 1, c0:c1 + 1]

        # mask = filled stem cross-section from the skeleton discs (union across
        # the band), rasterized on the full grid then cropped. This is a solid
        # disc, not a ring: interior cells are set even with no points present.
        full_mask = grid.rasterize_discs(tree_discs[int(t)], value=1)
        mask = full_mask[r0:r1 + 1, c0:c1 + 1].astype(np.uint8)

        x0 = grid.x_min + c0 * grid.cell
        y0 = grid.y_max - r0 * grid.cell
        _write_geotiff(os.path.join(
            trees_dir, f"{plot_uid}_tree_{uid}_stem_map_bounding_box.tif"),
            crop, x0, y0, grid.cell)
        _write_geotiff(os.path.join(
            trees_dir, f"{plot_uid}_tree_{uid}_stem_map_semantic_mask.tif"),
            mask, x0, y0, grid.cell)
        n_written += 1

    # link from the plot scene JSON
    rel = f"{plot_uid}/plot_placement"
    stem_map_meta = {
        "raster": f"{rel}/{raster_name}",
        "bounding_boxes": f"{rel}/{bbox_name}",
        "cell_size_m": float(cell_size),
        "dbh_height_m": float(dbh_height),
        "band_half_thickness_m": float(threshold),
        "bbox_padding_px": pad,
        "band_z_range_m": [round(z_lo, 3), round(z_hi, 3)],
        "source_cloud": f"{rel}/{os.path.basename(las_path)}",
        "raster_size": {"width": grid.width, "height": grid.height},
        "n_trees_with_stem": len(boxes),
    }
    scene_path = os.path.join(place_dir, f"{plot_uid}_scene.json")
    if os.path.exists(scene_path):
        scene = json.load(open(scene_path))
        scene["stem_map"] = stem_map_meta
        with open(scene_path, "w") as f:
            json.dump(scene, f, indent=2)
    else:
        stem_map_meta["_note"] = "scene JSON not found; stem_map metadata orphaned"

    if verbose:
        print(f"{plot_uid}: stem map {grid.width}x{grid.height} px @ "
              f"{cell_size*100:g} cm, band {z_lo:.2f}-{z_hi:.2f} m, "
              f"{band.sum():,} band points, {len(boxes)} stems, "
              f"{n_written} tree crops+masks")
    return stem_map_meta


def _load_skeletons(trees_dir, plot_uid):
    """Map tree_id -> trimmed skeleton dict for every tree in the plot.

    The skeleton is the ground truth for stem geometry (centerline + taper), so
    the stem mask is built from it, not from the point cloud. Trees are keyed by
    their integer tree_id (from the skeleton's own field).
    """
    out = {}
    if not os.path.isdir(trees_dir):
        return out
    for name in sorted(os.listdir(trees_dir)):
        tdir = os.path.join(trees_dir, name)
        if not os.path.isdir(tdir):
            continue
        # skeleton lives under textual/ (canonical) or the tree dir (flat)
        cand = [os.path.join(tdir, "textual", f"{name}_trimmed_skeleton.json"),
                os.path.join(tdir, f"{name}_trimmed_skeleton.json")]
        # flat layout: <trees_dir>/<uid>_trimmed_skeleton.json
        cand.append(os.path.join(trees_dir, f"{name}_trimmed_skeleton.json"))
        for p in cand:
            if os.path.exists(p):
                try:
                    d = json.load(open(p))
                    out[int(d.get("tree_id"))] = d
                except Exception:
                    pass
                break
    # flat layout fallback: skeleton files directly in trees_dir
    if not out:
        for fn in sorted(os.listdir(trees_dir)):
            if fn.endswith("_trimmed_skeleton.json"):
                try:
                    d = json.load(open(os.path.join(trees_dir, fn)))
                    out[int(d.get("tree_id"))] = d
                except Exception:
                    pass
    return out


def _tree_uid(tree_id, csv_path):
    """Best-effort zero-padded TREE_UID string for a numeric tree id."""
    try:
        df = pd.read_csv(csv_path, dtype={"TREE_UID": str})
        hit = df[df["TREE_UID"].astype(int) == int(tree_id)]
        if len(hit):
            return str(hit.iloc[0]["TREE_UID"])
    except Exception:
        pass
    return f"{int(tree_id):03d}"


def main():
    ap = argparse.ArgumentParser(description="Breast-height stem map raster for a plot")
    ap.add_argument("csv", nargs="+", help="plot table csv(s)")
    ap.add_argument("-d", "--plot-dir", default=None)
    ap.add_argument("-c", "--cell-size", type=float, default=DEFAULT_CELL_SIZE,
                    help="metres per pixel (default 0.01 = 1 cm)")
    ap.add_argument("--dbh-height", type=float, default=BREAST_HEIGHT_M)
    ap.add_argument("--threshold", type=float, default=DEFAULT_THRESHOLD,
                    help="half-thickness of the z-band in metres (default 0.10)")
    ap.add_argument("--bbox-padding", type=int, default=1,
                    help="pixels of padding on each side of a tree bbox (default 1)")
    ap.add_argument("--source", default="raw", choices=["raw", "trimmed"])
    args = ap.parse_args()

    for path in args.csv:
        build_stem_map(path, plot_dir=args.plot_dir, cell_size=args.cell_size,
                       dbh_height=args.dbh_height, threshold=args.threshold,
                       bbox_padding=args.bbox_padding, source=args.source)


if __name__ == "__main__":
    main()
