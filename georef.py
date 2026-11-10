"""
georef.py: the plot-anchor CRS module (spec Section 10, item W5).

Geographic position is optional. When the plot center reaches the
pipeline (FIA PLOT join or a hand-authored tree list), the scene JSON
records a `geo` block holding the UTM anchor computed from it with the
`utm` library (zone from the plot center; EPSG 326zz north / 327zz
south). Internal coordinates stay local metres -- they are metre offsets
from the anchor's UTM position, and placement centers every shape on
the local origin `(0, 0)` (the anchor's local point, recorded in the
block). Exports convert through the anchor:

  - EPSG:4326 GeoTIFF transforms from local raster corners, with metre
    cells converted to degrees at the plot latitude (Section 10);
  - lon/lat/elev positions for the info JSONs (elev = plot-center
    EPSG:5773 height + local z).

Unanchored plots keep the local metre frame and gain no CRS metadata.
"""

from __future__ import annotations

import json
import math
import os

import layout

M_PER_DEG_LAT = 111320.0        # documented metre-per-degree factor (Section 10)
LAT_RANGE = (-90.0, 90.0)
LON_RANGE = (-180.0, 180.0)
ELEV_RANGE = (-500.0, 5000.0)   # EPSG:5773 limits of agreement (Section 10)

# horizontal WGS 84 + vertical NAVD88 height (EPSG:4326 + EPSG:5773)
DTM_CRS_EPSG = 9707
GEOGRAPHIC_CRS = "EPSG:4326"


def _checked(value, name, bounds):
    if value is None:
        return None
    lo, hi = bounds
    if not (lo <= float(value) <= hi):
        raise ValueError(
            f"plot center {name}={value} outside the limits of agreement "
            f"[{lo}, {hi}] (spec Section 10)")
    return float(value)


def anchor_from_latlon(lat, lon, elev=None) -> dict | None:
    """The UTM anchor for a plot center, or None when lat/lon are absent.

    elev is the EPSG:5773 height of the plot center; the anchor's local
    point is the placement origin (0, 0) by construction (spec Section 10).
    """
    lat = _checked(lat, "lat", LAT_RANGE)
    lon = _checked(lon, "lon", LON_RANGE)
    elev = _checked(elev, "elev", ELEV_RANGE)
    if lat is None or lon is None:
        return None
    import utm
    easting, northing, zone, _ = utm.from_latlon(lat, lon)
    northern = lat >= 0.0
    return {
        "lat": float(lat),
        "lon": float(lon),
        "elev": elev,
        "epsg": (32600 if northern else 32700) + int(zone),
        "utm_zone": int(zone),
        "utm_hemisphere": "N" if northern else "S",
        "easting": float(easting),
        "northing": float(northing),
        "local_xy": [0.0, 0.0],
    }


def scene_record(anchor: dict) -> dict:
    """The scene JSON's `geo` block for a computed anchor."""
    return {
        "epsg": anchor["epsg"],
        "utm_zone": anchor["utm_zone"],
        "utm_hemisphere": anchor["utm_hemisphere"],
        "easting_m": anchor["easting"],
        "northing_m": anchor["northing"],
        "elev_m": anchor["elev"],
        "local_xy": list(anchor["local_xy"]),
        "lon": anchor["lon"],
        "lat": anchor["lat"],
    }


def anchor_from_scene(scene: dict) -> dict | None:
    """Rebuild an anchor from a scene's `geo` block (None when absent)."""
    geo = (scene or {}).get("geo")
    if not geo or geo.get("easting_m") is None or geo.get("northing_m") is None:
        return None
    return {
        "lat": geo.get("lat"),
        "lon": geo.get("lon"),
        "elev": geo.get("elev_m"),
        "epsg": geo.get("epsg"),
        "utm_zone": geo.get("utm_zone"),
        "utm_hemisphere": geo.get("utm_hemisphere"),
        "easting": geo.get("easting_m"),
        "northing": geo.get("northing_m"),
        "local_xy": list(geo.get("local_xy") or [0.0, 0.0]),
    }


def load_anchor(plot_dir: str, plot_uid: str) -> dict | None:
    """The anchor of a built plot, read from its scene JSON.

    Both locations are checked: the canonical one
    (plot_placement/textual/) and the pre-restructure one
    (plot_placement/, before restructure files the scene into place),
    because the stem and canopy maps run before restructure."""
    scene_paths = (
        layout.plot_path(plot_dir, plot_uid, "scene", create=False),
        os.path.join(plot_dir, layout.PLOT_PLACEMENT,
                     f"{plot_uid}_scene.json"),
    )
    for scene_path in scene_paths:
        if not os.path.exists(scene_path):
            continue
        try:
            with open(scene_path) as fh:
                return anchor_from_scene(json.load(fh))
        except (OSError, json.JSONDecodeError, ValueError):
            continue
    return None


# ---------------------------------------------------------------------------
# conversions local -> real world
# ---------------------------------------------------------------------------
def local_to_utm(anchor: dict, x: float, y: float) -> tuple[float, float]:
    """UTM easting/northing of a local point (anchor + offset)."""
    ax, ay = anchor["local_xy"]
    return anchor["easting"] + (x - ax), anchor["northing"] + (y - ay)


def local_to_deg(anchor: dict, x: float, y: float) -> tuple[float, float]:
    """(lon, lat) of a local point via the anchor's UTM zone."""
    import utm
    easting, northing = local_to_utm(anchor, x, y)
    northern = anchor["utm_hemisphere"] != "S"
    lat, lon = utm.to_latlon(easting, northing, anchor["utm_zone"],
                             northern=northern)
    return float(lon), float(lat)


def deg_cell_size(cell_m: float, lat: float) -> tuple[float, float]:
    """(deg per cell in x, deg per cell in y) at the plot's latitude
    (Section 10: 1 deg latitude ~ 111320 m; longitude scales by 1/cos)."""
    return (cell_m / (M_PER_DEG_LAT * math.cos(math.radians(lat))),
            cell_m / M_PER_DEG_LAT)


def geo_transform(anchor: dict, x_min: float, y_max: float,
                  cell_m: float):
    """(rasterio transform, CRS) placing a local-metre raster grid in
    EPSG:4326: the upper-left corner converted to lon/lat, pixel sizes
    in degrees at the plot latitude. Pixel values are unchanged."""
    import rasterio
    from rasterio.transform import from_origin
    west, north = local_to_deg(anchor, x_min, y_max)
    cell_x, cell_y = deg_cell_size(cell_m, anchor["lat"])
    return from_origin(west, north, cell_x, cell_y), GEOGRAPHIC_CRS
