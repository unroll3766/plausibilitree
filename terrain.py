"""
terrain.py: drape a placed plot onto a DTM, add ground points and a
height-above-ground column, and classify points per the LAS spec.

Runs after restructure.py, before normalize.py. Reads the canonical placement
cloud, reads a DTM (ESRI ASCII .asc, metres, x/y/z), and writes:

    plot_placement/las/<PLOT_UID>_terrain.las
    individual_trees/<PLOT_UID_TREE_UID>/las/<..>_terrain.las

Steps:
  - The DTM is used as-is. Its horizontal centre is aligned to the plot's
    horizontal centre (the placement origin, 0,0). The DTM must cover the whole
    plot footprint; if it does not, an error is raised (no tiling/stretching).
  - Every tree is shifted vertically so its lowest point sits on the ground at
    the tree's x,y (rigid drop, geometry unchanged).
  - Ground points are sampled across the plot bounding box: `density` points per
    DTM cell, z = DTM height + vertical noise, classified 2 (bare earth).
  - Points are classified by height above ground: <=0.5 m -> 3 (low veg),
    0.5-2 m -> 4 (med veg), >2 m -> 5 (high veg); ground stays 2.
  - Every point gets a height_above_ground column = z - DTM(x,y).
  - Per tree, the tree plus surrounding points (all points within the tree's
    canopy x-y bounding box) are written to its terrain.las.

Requires: numpy, laspy (via las.py), layout.py
"""

from __future__ import annotations

import argparse
import json
import os

import numpy as np

import georef
import layout
from las import read_points_las, export_points_las, export_canonical_las, read_canonical_las
import schema

# LAS classification codes
CLASS_GROUND = 2
CLASS_LOW_VEG = 3      # <= 0.5 m
CLASS_MED_VEG = 4      # 0.5 - 2 m
CLASS_HIGH_VEG = 5     # > 2 m

# Default ground sampling density for the flat/None topography when the
# user provides none (spec Sections 3.3 / 8.7: warn and use a default).
DEFAULT_GROUND_DENSITY = 1.0


# ---------------------------------------------------------------------------
# DTM
# ---------------------------------------------------------------------------

class DTM:
    """An ESRI ASCII grid, sampled by bilinear interpolation in world x,y.

    Rows in the file run north-to-south; internally we store south-to-north so
    that row index increases with y.
    """

    def __init__(self, grid, xll, yll, cell, nodata):
        self.grid = grid                     # (nrows, ncols), row 0 = south
        self.xll = xll
        self.yll = yll
        self.cell = cell
        self.nrows, self.ncols = grid.shape
        self.nodata = nodata
        self.width = self.ncols * cell
        self.height = self.nrows * cell

    @classmethod
    def read_asc(cls, path):
        with open(path) as f:
            lines = f.read().splitlines()
        h = {}
        i = 0
        while i < len(lines) and len(lines[i].split()) == 2:
            k, v = lines[i].split()
            h[k.lower()] = float(v)
            i += 1
        nrows, ncols = int(h["nrows"]), int(h["ncols"])
        rows = [[float(x) for x in lines[i + r].split()] for r in range(nrows)]
        grid_ns = np.array(rows)             # row 0 = north edge
        grid = np.flipud(grid_ns)            # -> row 0 = south edge
        cell = h["cellsize"]
        xll = h.get("xllcorner", h.get("xllcenter", 0.0))
        yll = h.get("yllcorner", h.get("yllcenter", 0.0))
        return cls(grid, xll, yll, cell, h.get("nodata_value", -9999.0))

    @classmethod
    def flat(cls, width, height, cell_size=0.25):
        """A flat plane at z = 0 (spec Section 3.3, the None topography).

        Same DTM interface as generated/provided grids -- sampling returns
        0 everywhere -- so the whole terrain stage (tree draping, HAG,
        classification, ground sampling) runs unchanged over the plane."""
        grid = np.zeros((int(round(height / cell_size)),
                         int(round(width / cell_size))), dtype=float)
        return cls(grid, -width / 2.0, -height / 2.0, cell_size, -9999.0)

    @classmethod
    def generate(cls, width, height, cell_size, bumpiness, relief, seed=None):
        """Synthesize a DTM of the given footprint via generate_dtm.

        width/height are the world extent to cover (metres). The grid's lower-
        left is left at (0,0) here; the caller recenters it onto the plot.
        Uses the project's diamond-square generator so provided and generated
        DTMs share the same terrain model.
        """
        from generate_dtm import generate_dtm
        grid, xll, yll = generate_dtm(
            width=width, height=height, cell_size=cell_size,
            bumpiness=bumpiness, relief=relief, seed=seed)
        # generate_dtm returns grid with row 0 at the south edge already
        return cls(np.asarray(grid, dtype=float), xll, yll, cell_size, -9999.0)

    def recenter_to(self, cx, cy):
        """Shift the grid so its horizontal centre sits at world (cx, cy)."""
        self.xll = cx - self.width / 2.0
        self.yll = cy - self.height / 2.0

    def write_asc(self, path):
        """Write this DTM as an ESRI ASCII grid (metres, row 0 = north edge)."""
        os.makedirs(os.path.dirname(path), exist_ok=True)
        grid_ns = np.flipud(self.grid)          # store north-first on disk
        with open(path, "w") as f:
            f.write(f"ncols {self.ncols}\n")
            f.write(f"nrows {self.nrows}\n")
            f.write(f"xllcorner {self.xll}\n")
            f.write(f"yllcorner {self.yll}\n")
            f.write(f"cellsize {self.cell}\n")
            f.write(f"NODATA_value {self.nodata}\n")
            for row in grid_ns:
                f.write(" ".join(repr(float(v)) for v in row) + "\n")

    def write_geotiff(self, anchor, path):
        """Georeferenced DTM GeoTIFF for an anchored plot (spec Section 10).

        Compound CRS EPSG:9707 (EPSG:4326 horizontal + NAVD88 height
        EPSG:5773 vertical); the band carries absolute heights = DTM
        height + the anchor's plot-center elevation (0 when the input
        carried no elevation). Nodata cells stay nodata."""
        import rasterio
        from rasterio.crs import CRS

        os.makedirs(os.path.dirname(path), exist_ok=True)
        transform = georef.geo_transform(
            anchor, self.xll, self.yll + self.height, self.cell)[0]
        band = np.flipud(self.grid) + float(anchor.get("elev") or 0.0)
        band = np.where(self.grid == self.nodata, self.nodata, band)
        with rasterio.open(
            path, "w", driver="GTiff",
            height=self.nrows, width=self.ncols, count=1, dtype="float32",
            transform=transform, crs=CRS.from_epsg(georef.DTM_CRS_EPSG),
            nodata=self.nodata, compress="lzw",
        ) as dst:
            dst.write(band.astype(np.float32), 1)

    def crop(self, x_min, y_min, x_max, y_max):
        """Return a new DTM covering [x_min,x_max] x [y_min,y_max].

        Snaps outward to whole cells so the crop fully contains the box; clamps
        to this DTM's extent. Used to store the slice of a generated DTM under a
        tree's footprint (its canopy bounding box, with whatever buffer it has).
        """
        c0 = int(np.floor((x_min - self.xll) / self.cell))
        c1 = int(np.ceil((x_max - self.xll) / self.cell))
        r0 = int(np.floor((y_min - self.yll) / self.cell))
        r1 = int(np.ceil((y_max - self.yll) / self.cell))
        c0 = max(0, c0); r0 = max(0, r0)
        c1 = min(self.ncols, c1); r1 = min(self.nrows, r1)
        if c1 <= c0 or r1 <= r0:
            return None
        sub = self.grid[r0:r1, c0:c1]
        return DTM(sub.copy(), self.xll + c0 * self.cell,
                   self.yll + r0 * self.cell, self.cell, self.nodata)

    def covers(self, x_min, y_min, x_max, y_max, eps=1e-6):
        return (x_min >= self.xll - eps and y_min >= self.yll - eps and
                x_max <= self.xll + self.width + eps and
                y_max <= self.yll + self.height + eps)

    def sample(self, x, y):
        """Bilinearly interpolated height at world x,y (arrays)."""
        # continuous cell coordinate; cell centres are at (i+0.5)
        fx = (np.asarray(x) - self.xll) / self.cell - 0.5
        fy = (np.asarray(y) - self.yll) / self.cell - 0.5
        x0 = np.clip(np.floor(fx).astype(int), 0, self.ncols - 1)
        y0 = np.clip(np.floor(fy).astype(int), 0, self.nrows - 1)
        x1 = np.clip(x0 + 1, 0, self.ncols - 1)
        y1 = np.clip(y0 + 1, 0, self.nrows - 1)
        tx = np.clip(fx - x0, 0.0, 1.0)
        ty = np.clip(fy - y0, 0.0, 1.0)
        g = self.grid
        v00 = g[y0, x0]; v10 = g[y0, x1]
        v01 = g[y1, x0]; v11 = g[y1, x1]
        return (v00 * (1 - tx) * (1 - ty) + v10 * tx * (1 - ty) +
                v01 * (1 - tx) * ty + v11 * tx * ty)


# ---------------------------------------------------------------------------
# classification
# ---------------------------------------------------------------------------

def classify_by_hag(hag):
    """LAS vegetation class from height above ground (array)."""
    cls = np.full(len(hag), CLASS_HIGH_VEG, dtype=np.uint8)
    cls[hag <= 2.0] = CLASS_MED_VEG
    cls[hag <= 0.5] = CLASS_LOW_VEG
    return cls


# ---------------------------------------------------------------------------
# ground sampling
# ---------------------------------------------------------------------------

def sample_ground(dtm, x_min, y_min, x_max, y_max, density, noise, rng):
    """`density` points per DTM cell over the bbox; z = DTM + noise; class 2.

    Returns (N, 5+) array: x, y, z, semantic(class), tree_id(0), then the ground
    height_above_ground (~ the noise) is added by the caller for one HAG column.
    """
    area_cells = ((x_max - x_min) / dtm.cell) * ((y_max - y_min) / dtm.cell)
    n = max(0, int(round(density * area_cells)))
    if n == 0:
        return np.zeros((0, 5)), np.zeros(0)
    x = rng.uniform(x_min, x_max, n)
    y = rng.uniform(y_min, y_max, n)
    z_ground = dtm.sample(x, y)
    dz = rng.uniform(-noise, noise, n) if noise > 0 else np.zeros(n)
    z = z_ground + dz
    pts = np.column_stack([x, y, z,
                           np.full(n, CLASS_GROUND, dtype=float),
                           np.zeros(n)])
    hag = dz                                          # z - DTM(x,y)
    return pts, hag


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def build_terrain(csv_path, dtm_path=None, plot_dir=None, *, density=None,
                  noise=0.0,
                  source="trimmed", seed=0, generate=False, dtm_resolution=1.0,
                  gen_bumpiness=0.5, gen_relief=10.0, verbose=True):
    plot_uid = os.path.splitext(os.path.basename(csv_path))[0]
    if plot_dir is None:
        plot_dir = os.path.join(os.path.dirname(os.path.abspath(csv_path)), plot_uid)

    # --- topography resolution (spec Section 3.3) ---
    # None wins if present, else filepath, else random parameters; a warning
    # describes any choice that ignored an input. With neither a DTM path nor
    # generation requested, the topography IS None: a uniform flat plane at
    # z = 0.
    flat = False
    if dtm_path:
        if generate:
            print(f"  [terrain] {plot_uid}: warning - both a provided DTM and "
                  f"--generate-dtm were requested; using the provided DTM "
                  f"(filepath wins, generation ignored)")
        generate = False
    elif not generate:
        flat = True
    if flat and density is None:
        print(f"  [terrain] {plot_uid}: warning - no ground sampling density "
              f"provided for flat topography; using the default "
              f"{DEFAULT_GROUND_DENSITY} points per cell")
    if density is None:
        density = DEFAULT_GROUND_DENSITY

    # source cloud (canonical layout; restructure must have run)
    src = layout.plot_path(plot_dir, plot_uid, f"{source}_las", create=False)
    if not os.path.exists(src):
        alt = layout.plot_path(plot_dir, plot_uid, "raw_las", create=False)
        if os.path.exists(alt):
            src = alt
        else:
            raise FileNotFoundError(
                f"{plot_uid}: source cloud not found ({src}). Run restructure.py "
                f"before terrain.py.")

    src_xyz, src_cols = read_canonical_las(src, level="plot")
    xyz = src_xyz
    tid = src_cols["tree_id"].astype(int)
    src_sem = src_cols["semantic_code"].astype(np.int8)      # object identity
    src_code = src_cols["branch_code"].astype(np.int8)       # branch depth
    src_bid = src_cols["branch_id"].astype(np.int32)         # branch instance

    # plot horizontal centre = placement origin (0,0)
    cx, cy = 0.0, 0.0
    x_min, y_min = float(xyz[:, 0].min()), float(xyz[:, 1].min())
    x_max, y_max = float(xyz[:, 0].max()), float(xyz[:, 1].max())

    if generate:
        # bounds auto-sized to the plot point-cloud bbox (symmetric about the
        # centre so the recentred grid always covers it). One cell of margin
        # avoids edge points landing exactly on the grid boundary.
        half_w = max(abs(x_min - cx), abs(x_max - cx)) + dtm_resolution
        half_h = max(abs(y_min - cy), abs(y_max - cy)) + dtm_resolution
        width = 2.0 * half_w
        height = 2.0 * half_h
        dtm = DTM.generate(width, height, dtm_resolution, gen_bumpiness,
                           gen_relief, seed=seed)
        dtm.recenter_to(cx, cy)
        dtm_desc = (f"generated {dtm.width:.1f}x{dtm.height:.1f} m "
                    f"(cell {dtm_resolution}, bumpiness {gen_bumpiness}, "
                    f"relief {gen_relief})")
    elif flat:
        # the None topography: a uniform flat plane at z = 0, sized to the
        # same auto-bounds so recentring/coverage checks behave identically
        half_w = max(abs(x_min - cx), abs(x_max - cx)) + dtm_resolution
        half_h = max(abs(y_min - cy), abs(y_max - cy)) + dtm_resolution
        dtm = DTM.flat(2.0 * half_w, 2.0 * half_h, dtm_resolution)
        dtm.recenter_to(cx, cy)
        dtm_desc = "flat plane (topography = None, z = 0)"
    else:
        dtm = DTM.read_asc(dtm_path)
        dtm.recenter_to(cx, cy)
        if not dtm.covers(x_min, y_min, x_max, y_max):
            raise ValueError(
                f"{plot_uid}: DTM does not cover the plot footprint. "
                f"DTM {dtm.width:.1f}x{dtm.height:.1f} m centred at ({cx},{cy}) "
                f"spans x[{dtm.xll:.1f},{dtm.xll+dtm.width:.1f}] y[{dtm.yll:.1f},"
                f"{dtm.yll+dtm.height:.1f}]; plot needs x[{x_min:.1f},{x_max:.1f}] "
                f"y[{y_min:.1f},{y_max:.1f}]. Provide a larger DTM or use "
                f"--generate-dtm (no tiling/stretching is applied).")
        dtm_desc = os.path.basename(dtm_path)

    rng = np.random.default_rng(seed)

    # --- drape each tree per spec Section 3.3: drop it vertically so its
    # lowest point sits on the DTM at the tree's (X, Y) -- the BASE
    # position from the plot table, NOT whichever cloud point happens to
    # be lowest (a downhill branch tip could otherwise anchor the tree to
    # a gully elevation and bury the trunk into an upslope). ---
    base_xy = {}
    try:
        import tree_list
        fdf = tree_list.trees_frame(csv_path)
        base_xy = {int(str(u)): (float(x), float(y))
                   for u, x, y in zip(fdf["TREE_UID"], fdf["X"], fdf["Y"])}
    except Exception:
        base_xy = {}
    draped = xyz.copy()
    for t in np.unique(tid):
        m = tid == t
        sub = xyz[m]
        i_low = np.argmin(sub[:, 2])
        bx, by = base_xy.get(int(t), (sub[i_low, 0], sub[i_low, 1]))
        ground_z = float(dtm.sample(np.array([bx]), np.array([by]))[0])
        shift = ground_z - sub[i_low, 2]
        draped[m, 2] += shift

    # height above ground for tree points = z - DTM(x,y)
    veg_ground = dtm.sample(draped[:, 0], draped[:, 1])
    veg_hag = draped[:, 2] - veg_ground
    veg_class = classify_by_hag(veg_hag)          # LAS classification (height-based)

    # tree points: xyz + classification + tree_id. semantic_code and branch_code
    # are the object identity -- preserved unchanged from the source cloud
    # (orthogonal to the height-based LAS classification).
    veg = np.column_stack([draped[:, 0], draped[:, 1], draped[:, 2],
                           veg_class.astype(float), tid.astype(float)])
    veg_sem = src_sem
    veg_code = src_code
    veg_bid = src_bid

    # --- ground points ---
    ground, ground_hag = sample_ground(dtm, x_min, y_min, x_max, y_max,
                                        density, noise, rng)
    ground_sem = np.full(len(ground), schema.SEM_GROUND, dtype=np.int8)   # 6
    ground_code = np.zeros(len(ground), dtype=np.int8)                    # not-tree
    ground_bid = np.zeros(len(ground), dtype=np.int32)                    # not-branch

    all_pts = np.vstack([veg, ground]) if len(ground) else veg
    all_hag = np.concatenate([veg_hag, ground_hag]) if len(ground) else veg_hag
    all_sem = np.concatenate([veg_sem, ground_sem]) if len(ground) else veg_sem
    all_code = np.concatenate([veg_code, ground_code]) if len(ground) else veg_code
    all_bid = np.concatenate([veg_bid, ground_bid]) if len(ground) else veg_bid

    # write plot terrain cloud (canonical schema)
    out = layout.plot_path(plot_dir, plot_uid, "terrain_las")
    export_canonical_las(
        out, all_pts[:, :3],
        {"tree_id": all_pts[:, 4].astype(np.int32),
         "semantic_code": all_sem,
         "classification": all_pts[:, 3].astype(np.uint8),
         "branch_code": all_code,
         "branch_id": all_bid,
         "height_above_ground": all_hag.astype(np.float32)},
        level="plot")

    # save the generated DTM (plot level) as an .asc in textual/. Provided DTMs
    # are not re-saved (the user already has them).
    if generate:
        dtm_out = os.path.join(plot_dir, layout.PLOT_PLACEMENT, "textual",
                               f"{plot_uid}_dtm.asc")
        dtm.write_asc(dtm_out)
        # EPSG:9707 DTM GeoTIFF rides along when the plot is anchored (W5)
        anchor = georef.load_anchor(plot_dir, plot_uid)
        if anchor is not None:
            dtm.write_geotiff(
                anchor, layout.plot_path(plot_dir, plot_uid, "dtm_tif"))

    # --- per-tree terrain clouds: tree + surrounding points in canopy bbox ---
    n_tree_files = _write_tree_terrain(plot_dir, plot_uid, all_pts, all_hag,
                                       all_sem, all_code, all_bid,
                                       dtm=(dtm if generate else None))

    # record in scene
    scene_path = layout.plot_path(plot_dir, plot_uid, "scene", create=False)
    if os.path.exists(scene_path):
        scene = json.load(open(scene_path))
        scene["terrain"] = {
            "dtm": dtm_desc,
            "dtm_generated": bool(generate),
            "dtm_cell_m": dtm.cell,
            "dtm_extent_m": [dtm.width, dtm.height],
            "ground_density_per_cell": density,
            "vertical_noise_m": noise,
            "source": source,
            "n_ground_points": int(len(ground)),
            "n_veg_points": int(len(veg)),
            "classes": {"2": "ground", "3": "low_veg<=0.5m",
                        "4": "med_veg_0.5-2m", "5": "high_veg>2m"},
            "semantic_note": ("semantic_code carries object identity "
                              "(1 trunk, 2 branch, 3 foliage, 6 ground), "
                              "orthogonal to the height-based classification"),
            "plot_terrain_las": layout.plot_rel(plot_dir, plot_uid, "terrain_las"),
        }
        with open(scene_path, "w") as f:
            json.dump(scene, f, indent=2)

    if verbose:
        vc = np.unique(veg_class, return_counts=True)
        dist = {int(k): int(v) for k, v in zip(*vc)}
        print(f"{plot_uid}: terrain cloud {len(all_pts):,} points "
              f"({len(veg):,} veg + {len(ground):,} ground), "
              f"veg classes {dist}, {n_tree_files} tree files")
    return out


def _write_tree_terrain(plot_dir, plot_uid, all_pts, all_hag, all_sem, all_code,
                        all_bid, dtm=None):
    """Write each tree's terrain.las: all points inside its canopy x-y bbox.

    When `dtm` is given (a generated DTM), also writes the DTM slice under each
    tree's canopy bounding box as <tag>_dtm.asc in the tree's textual/ dir.
    """
    scene_path = layout.plot_path(plot_dir, plot_uid, "scene", create=False)
    tree_uids = []
    if os.path.exists(scene_path):
        tree_uids = list(json.load(open(scene_path)).get("trees", {}).keys())
    if not tree_uids:
        d = os.path.join(plot_dir, layout.INDIVIDUAL_TREES)
        tree_uids = [n.rsplit("_", 1)[-1] for n in sorted(os.listdir(d))
                     if os.path.isdir(os.path.join(d, n))]

    x, y = all_pts[:, 0], all_pts[:, 1]
    n = 0
    for uid in tree_uids:
        bbox_path = layout.tree_path(plot_dir, plot_uid, uid, "canopy_bbox_json",
                                     create=False)
        if not os.path.exists(bbox_path):
            continue
        bw = json.load(open(bbox_path)).get("bbox_world")
        if not bw:
            continue
        m = ((x >= bw["x_min"]) & (x <= bw["x_max"]) &
             (y >= bw["y_min"]) & (y <= bw["y_max"]))
        if not m.any():
            continue
        out = layout.tree_path(plot_dir, plot_uid, uid, "terrain_las")
        export_canonical_las(
            out, all_pts[m, :3],
            {"tree_id": all_pts[m, 4].astype(np.int32),
             "semantic_code": all_sem[m],
             "classification": all_pts[m, 3].astype(np.uint8),
             "branch_code": all_code[m],
             "branch_id": all_bid[m],
             "height_above_ground": all_hag[m].astype(np.float32)},
            level="tree")

        # DTM crop under the tree's canopy bbox (generated DTMs only)
        if dtm is not None:
            sub = dtm.crop(bw["x_min"], bw["y_min"], bw["x_max"], bw["y_max"])
            if sub is not None:
                dtm_out = os.path.join(layout.tree_dir(plot_dir, plot_uid, uid),
                                       "textual",
                                       f"{layout.tree_tag(plot_uid, uid)}_dtm.asc")
                sub.write_asc(dtm_out)
        n += 1
    return n


def main():
    ap = argparse.ArgumentParser(
        description="Drape a plot onto a DTM (provided, generated, or a "
                    "flat plane when neither is given), add ground + HAG")
    ap.add_argument("csv", nargs="+", help="plot table csv(s)")
    ap.add_argument("--dtm", help="provided DTM .asc (metres); wins over "
                                   "--generate-dtm (spec Section 3.3)")
    ap.add_argument("--generate-dtm", action="store_true",
                    help="synthesize a DTM sized to each plot's footprint")
    ap.add_argument("--dtm-resolution", type=float, default=0.25,
                    help="generated/flat DTM cell size, m (spec DTM_RESOLUTION)")
    ap.add_argument("--gen-bumpiness", type=float, default=0.5,
                    help="generated DTM bumpiness/roughness (with --generate-dtm)")
    ap.add_argument("--gen-relief", type=float, default=10.0,
                    help="generated DTM relief / vertical range, m (with --generate-dtm)")
    ap.add_argument("-d", "--plot-dir", default=None)
    ap.add_argument("--density", type=float, default=None,
                    help="ground points per DTM cell (default 1)")
    ap.add_argument("--noise", type=float, default=0.0,
                    help="vertical noise half-range in metres (default 0)")
    ap.add_argument("--source", default="trimmed", choices=["raw", "trimmed"])
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    for path in args.csv:
        build_terrain(path, args.dtm, plot_dir=args.plot_dir,
                      density=args.density, noise=args.noise,
                      source=args.source, seed=args.seed,
                      generate=args.generate_dtm,
                      dtm_resolution=args.dtm_resolution,
                      gen_bumpiness=args.gen_bumpiness,
                      gen_relief=args.gen_relief)


if __name__ == "__main__":
    main()
