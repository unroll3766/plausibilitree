"""
generate_dtm.py

Generate a synthetic Digital Terrain Model (DTM) using fractal
(diamond-square) terrain synthesis.

Usage (CLI):
    # width/height are in meters; the grid resolution (number of cells) is
    # derived from width/height divided by cell-size.
    python generate_dtm.py --width 200 --height 200 --cell-size 1.0 --bumpiness 0.6 \
        --relief 15 --seed 42 --output my_dtm

    # With real-world coordinates:
    python generate_dtm.py --width 400 --height 400 --cell-size 2.0 --bumpiness 0.3 \
        --origin-x 500000 --origin-y 6200000 --origin-z 120 --epsg 32755 \
        --output surveyed_area

Usage (as a library):
    from generate_dtm import generate_dtm, write_outputs

    grid, xll, yll = generate_dtm(width=256, height=256, cell_size=1.0, bumpiness=0.5)
    write_outputs("my_dtm", grid, cell_size=1.0, xllcorner=xll, yllcorner=yll)

Outputs:
    <output>.asc        ESRI ASCII Grid (readable by QGIS, ArcGIS, GDAL, etc.)
    <output>.prj        Projection file (only written when --epsg is recognised)
    <output>_meta.json  Metadata describing how the DTM was generated
    <output>_preview.png  Quick-look shaded image of the terrain

Author: generated for tree-simulation project
"""

import argparse
import json

import numpy as np
from scipy.ndimage import gaussian_filter


# ---------------------------------------------------------------------------
# Terrain synthesis (diamond-square fractal algorithm)
# ---------------------------------------------------------------------------

def _next_pow2_plus_one(n):
    """Smallest size = 2**k + 1 that is >= n."""
    size = 1
    while size + 1 < n:
        size *= 2
    return size + 1


def _diamond_square(size, hurst, rng):
    """
    Generate a fractal heightfield of shape (size, size), size = 2**k + 1,
    using the diamond-square algorithm. Output is unnormalised (roughly
    in the range [-1, 1] but not guaranteed).

    hurst: float in (0, 1]. Lower values produce rougher/bumpier terrain,
           higher values produce smoother, more rolling terrain.
    """
    grid = np.zeros((size, size), dtype=np.float64)

    grid[0, 0] = rng.uniform(-1, 1)
    grid[0, -1] = rng.uniform(-1, 1)
    grid[-1, 0] = rng.uniform(-1, 1)
    grid[-1, -1] = rng.uniform(-1, 1)

    step = size - 1
    scale = 1.0

    while step > 1:
        half = step // 2

        # --- diamond step: set center of each square from its 4 corners ---
        for y in range(half, size - 1, step):
            for x in range(half, size - 1, step):
                avg = (grid[y - half, x - half] + grid[y - half, x + half] +
                       grid[y + half, x - half] + grid[y + half, x + half]) / 4.0
                grid[y, x] = avg + rng.uniform(-scale, scale)

        # --- square step: set center of each diamond from its 4 neighbours ---
        for y in range(0, size, half):
            for x in range((y // half + 1) % 2 * half, size, step):
                pts = []
                if y - half >= 0:
                    pts.append(grid[y - half, x])
                if y + half < size:
                    pts.append(grid[y + half, x])
                if x - half >= 0:
                    pts.append(grid[y, x - half])
                if x + half < size:
                    pts.append(grid[y, x + half])
                grid[y, x] = float(np.mean(pts)) + rng.uniform(-scale, scale)

        step = half
        scale *= 2.0 ** (-hurst)

    return grid


def generate_dtm(width, height, cell_size=1.0, bumpiness=0.5, relief=10.0,
                  seed=None, smooth=True,
                  origin_x=0.0, origin_y=0.0, origin_z=0.0):
    """
    Generate a random DTM.

    Parameters
    ----------
    width, height : real-world size of the DTM in meters (x and y extent).
                  The number of grid cells (nx, ny) is derived as
                  round(width / cell_size) and round(height / cell_size).
    cell_size   : size of each grid cell, in meters
    bumpiness   : float in [0, 1]; 0 = smooth/rolling, 1 = very rough/bumpy
    relief      : total vertical relief of the terrain (max - min elevation), in meters
    seed        : optional int, for reproducible terrain
    smooth      : whether to apply a light gaussian smoothing pass (reduces
                  diamond-square grid artifacts without removing bumpiness)
    origin_x, origin_y, origin_z : real-world coordinates of the DTM's
                  lower-left corner (origin_z is added as a vertical offset).
                  Left at 0.0 the output is in local meters starting at (0,0,0).

    Returns
    -------
    grid : (ny, nx) numpy array of elevations
    xllcorner, yllcorner : coordinates of the lower-left cell corner
    """
    if not (0.0 <= bumpiness <= 1.0):
        raise ValueError("bumpiness must be between 0 and 1")
    if cell_size <= 0:
        raise ValueError("cell_size must be positive")

    # derive the pixel grid resolution from the requested real-world extent
    nx = max(2, int(round(width / cell_size)) + 1)
    ny = max(2, int(round(height / cell_size)) + 1)

    rng = np.random.default_rng(seed)

    # map bumpiness -> Hurst-like exponent (low exponent = rougher terrain)
    hurst = 1.0 - 0.9 * bumpiness  # bumpiness=0 -> 1.0 (smooth), bumpiness=1 -> 0.1 (rough)

    size = _next_pow2_plus_one(max(nx, ny))
    grid = _diamond_square(size, hurst, rng)
    grid = grid[:ny, :nx]

    if smooth:
        grid = gaussian_filter(grid, sigma=0.6)

    # normalise so the terrain spans exactly `relief` meters of elevation
    grid = grid - grid.min()
    peak = grid.max()
    if peak > 0:
        grid = grid / peak * relief

    grid = grid + origin_z

    xllcorner = origin_x
    yllcorner = origin_y

    return grid, xllcorner, yllcorner


# ---------------------------------------------------------------------------
# Projection helpers (no external CRS libraries required)
# ---------------------------------------------------------------------------

def _utm_wkt(zone, northern):
    hemi = "Northern" if northern else "Southern"
    false_northing = 0 if northern else 10000000
    central_meridian = -183 + 6 * zone
    return (
        f'PROJCS["WGS 84 / UTM zone {zone}{"N" if northern else "S"}",'
        f'GEOGCS["WGS 84",DATUM["WGS_1984",'
        f'SPHEROID["WGS 84",6378137,298.257223563]],'
        f'PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]],'
        f'PROJECTION["Transverse_Mercator"],'
        f'PARAMETER["latitude_of_origin",0],'
        f'PARAMETER["central_meridian",{central_meridian}],'
        f'PARAMETER["scale_factor",0.9996],'
        f'PARAMETER["false_easting",500000],'
        f'PARAMETER["false_northing",{false_northing}],'
        f'UNIT["metre",1]]'
    )


_KNOWN_WKT = {
    4326: 'GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563]],'
          'PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]]',
    3857: 'PROJCS["WGS 84 / Pseudo-Mercator",GEOGCS["WGS 84",DATUM["WGS_1984",'
          'SPHEROID["WGS 84",6378137,298.257223563]],PRIMEM["Greenwich",0],'
          'UNIT["degree",0.0174532925199433]],PROJECTION["Mercator_1SP"],'
          'PARAMETER["central_meridian",0],PARAMETER["scale_factor",1],'
          'PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["metre",1]]',
}


def get_wkt_for_epsg(epsg):
    """
    Return an ESRI-style WKT string for a handful of common EPSG codes
    (4326, 3857, and any WGS84 UTM zone: 32601-32660 north, 32701-32760 south).
    Returns None if the EPSG code isn't recognised (no pyproj/rasterio
    available in this environment to look up arbitrary codes).
    """
    if epsg in _KNOWN_WKT:
        return _KNOWN_WKT[epsg]
    if 32601 <= epsg <= 32660:
        return _utm_wkt(epsg - 32600, northern=True)
    if 32701 <= epsg <= 32760:
        return _utm_wkt(epsg - 32700, northern=False)
    return None


# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------

def write_outputs(output_prefix, grid, cell_size, xllcorner, yllcorner,
                   epsg=None, nodata=-9999.0, preview=True, meta_extra=None):
    """Write the .asc grid, optional .prj, a metadata JSON, and a PNG preview."""

    asc_path = f"{output_prefix}.asc"
    _write_esri_ascii(asc_path, grid, cell_size, xllcorner, yllcorner, nodata)

    prj_path = None
    if epsg is not None:
        wkt = get_wkt_for_epsg(epsg)
        if wkt is not None:
            prj_path = f"{output_prefix}.prj"
            with open(prj_path, "w") as f:
                f.write(wkt)
        else:
            prj_path = f"{output_prefix}_epsg_note.txt"
            with open(prj_path, "w") as f:
                f.write(
                    f"EPSG:{epsg} was requested but is not one of the codes this "
                    "script can generate a .prj for automatically (supported: "
                    "4326, 3857, and any WGS84 UTM zone 32601-32660 / 32701-32760).\n"
                    f"Assign EPSG:{epsg} to {output_prefix}.asc manually in your GIS "
                    "software (the grid coordinates written into the .asc file are "
                    "already correct for that CRS).\n"
                )

    meta = {
        "ncols": int(grid.shape[1]),
        "nrows": int(grid.shape[0]),
        "cell_size_m": cell_size,
        "xllcorner": xllcorner,
        "yllcorner": yllcorner,
        "elevation_min": float(grid.min()),
        "elevation_max": float(grid.max()),
        "epsg": epsg,
    }
    if meta_extra:
        meta.update(meta_extra)
    meta_path = f"{output_prefix}_meta.json"
    with open(meta_path, "w") as f:
        json.dump(meta, f, indent=2)

    preview_path = None
    if preview:
        preview_path = f"{output_prefix}_preview.png"
        _write_preview(preview_path, grid, cell_size)

    return {
        "asc": asc_path,
        "prj": prj_path,
        "meta": meta_path,
        "preview": preview_path,
    }


def _write_esri_ascii(path, grid, cell_size, xllcorner, yllcorner, nodata):
    nrows, ncols = grid.shape
    with open(path, "w") as f:
        f.write(f"ncols {ncols}\n")
        f.write(f"nrows {nrows}\n")
        f.write(f"xllcorner {xllcorner}\n")
        f.write(f"yllcorner {yllcorner}\n")
        f.write(f"cellsize {cell_size}\n")
        f.write(f"NODATA_value {nodata}\n")
        # ESRI ASCII grids list rows north-to-south; grid row 0 is the south
        # edge (yllcorner), so flip vertically when writing.
        np.savetxt(f, np.flipud(grid), fmt="%.4f")


def _write_preview(path, grid, cell_size):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(6, 5))
    im = ax.imshow(grid, origin="lower", cmap="terrain",
                    extent=[0, grid.shape[1] * cell_size, 0, grid.shape[0] * cell_size])
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_title("Synthetic DTM preview")
    fig.colorbar(im, ax=ax, label="elevation (m)")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(description="Generate a random synthetic DTM.")
    p.add_argument("--width", type=float, default=200.0,
                    help="width of the DTM in meters (x extent)")
    p.add_argument("--height", type=float, default=200.0,
                    help="height of the DTM in meters (y extent)")
    p.add_argument("--cell-size", type=float, default=1.0, help="cell size in meters")
    p.add_argument("--bumpiness", type=float, default=0.5,
                    help="terrain roughness, 0 (smooth) to 1 (very bumpy)")
    p.add_argument("--relief", type=float, default=10.0,
                    help="total elevation range of the terrain, in meters")
    p.add_argument("--seed", type=int, default=None, help="random seed")
    p.add_argument("--no-smooth", action="store_true",
                    help="disable light gaussian smoothing of the raw fractal grid")
    p.add_argument("--origin-x", type=float, default=0.0)
    p.add_argument("--origin-y", type=float, default=0.0)
    p.add_argument("--origin-z", type=float, default=0.0)
    p.add_argument("--epsg", type=int, default=None,
                    help="EPSG code to georeference the output (optional)")
    p.add_argument("--output", type=str, default="dtm_output",
                    help="output file prefix")
    p.add_argument("--no-preview", action="store_true",
                    help="skip generating the PNG preview image")
    args = p.parse_args()

    grid, xll, yll = generate_dtm(
        width=args.width, height=args.height, cell_size=args.cell_size,
        bumpiness=args.bumpiness, relief=args.relief, seed=args.seed,
        smooth=not args.no_smooth,
        origin_x=args.origin_x, origin_y=args.origin_y, origin_z=args.origin_z,
    )

    paths = write_outputs(
        args.output, grid, cell_size=args.cell_size,
        xllcorner=xll, yllcorner=yll, epsg=args.epsg,
        preview=not args.no_preview,
        meta_extra={"bumpiness": args.bumpiness, "relief": args.relief, "seed": args.seed},
    )

    print(f"Grid: {grid.shape[1]} x {grid.shape[0]} cells @ {args.cell_size} m")
    print(f"Elevation range: {grid.min():.2f} to {grid.max():.2f} m")
    for k, v in paths.items():
        if v:
            print(f"  {k}: {v}")


if __name__ == "__main__":
    main()
