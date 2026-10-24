"""
build_info.py: emit the plot- and tree-level *_info.json files.

Runs after understory (reads clouds, skeletons, and canopy maps; order-free as
long as those exist). Replaces the old *_volume.json files. Produces:

  plot_placement/textual/<plot>_info.json
  individual_trees/<plot>_<uid>/textual/<plot>_<uid>_info.json   (one per tree)

Volumes:
  trunk   -- frustum sum along the skeleton centerline (taper radius per node)
  branch  -- cylinder per branch: pi * r(base_height)^2 * length
  foliage -- occupied-voxel volume of foliage points (semantic 3)
  canopy  -- convex-hull volume of foliage points
Each of trunk/branch/foliage has raw (all elements / all points) and trimmed
(non-trimmed elements / trimmed cloud) variants.

Crown major axis = the longest chord across the crown hull (x-y); minor axis =
the widest span perpendicular to it. Lengths in m, orientation in degrees CCW
from +x (local frame).
"""

from __future__ import annotations

import argparse
import json
import os

import numpy as np

import layout
import schema
from las import read_canonical_las

BREAST_HEIGHT_M = 1.37


# ---------------------------------------------------------------------------
# volume helpers
# ---------------------------------------------------------------------------
def _taper_radius(base_radius, total_height, h):
    if total_height <= 0:
        return 0.0
    return max(0.0, base_radius * (1.0 - h / total_height))


def trunk_volume(skel, trimmed_only):
    """Frustum sum along each trunk polyline using the taper radius per node.

    A frustum (truncated cone) between consecutive nodes at heights h0,h1 with
    radii r0,r1 and slant length L has volume (pi/3) L (r0^2 + r0 r1 + r1^2).
    Uses 3D segment length between nodes. When trimmed_only, segments whose
    end node is flagged node_trimmed are skipped.
    """
    base_r = float(skel.get("base_radius_m", 0.0))
    h_t = float(skel.get("height_m", 0.0))
    total = 0.0
    trunk = skel.get("trunk", {})
    polylines = trunk.get("polylines", [])
    segments = trunk.get("segments", [])
    for pi, pl in enumerate(polylines):
        arr = np.asarray(pl, dtype=float)
        if len(arr) < 2:
            continue
        node_trimmed = None
        if trimmed_only and pi < len(segments):
            nt = segments[pi].get("node_trimmed")
            if nt is not None and len(nt) == len(arr):
                node_trimmed = np.asarray(nt, dtype=bool)
        for k in range(len(arr) - 1):
            if node_trimmed is not None and (node_trimmed[k] or node_trimmed[k + 1]):
                continue
            p0, p1 = arr[k], arr[k + 1]
            L = float(np.linalg.norm(p1 - p0))
            r0 = _taper_radius(base_r, h_t, p0[2])
            r1 = _taper_radius(base_r, h_t, p1[2])
            total += (np.pi / 3.0) * L * (r0 * r0 + r0 * r1 + r1 * r1)
    return total


def branch_volume(skel, trimmed_only):
    """Cylinder per branch: pi * r^2 * length, r = taper radius at branch base."""
    base_r = float(skel.get("base_radius_m", 0.0))
    h_t = float(skel.get("height_m", 0.0))
    total = 0.0
    for b in skel.get("branches", []):
        if trimmed_only and b.get("trimmed"):
            continue
        base_h = float(b.get("base", [0, 0, 0])[2])
        r = _taper_radius(base_r, h_t, base_h)
        L = float(b.get("length_m", 0.0))
        total += np.pi * r * r * L
    return total


def voxel_volume(xyz, voxel_size, threshold):
    """Occupied-voxel volume: voxels holding >= threshold points, times size^3."""
    if len(xyz) == 0 or voxel_size <= 0:
        return 0.0
    keys = np.floor(np.asarray(xyz) / voxel_size).astype(np.int64)
    _, counts = np.unique(keys, axis=0, return_counts=True)
    return int((counts >= threshold).sum()) * (voxel_size ** 3)


# ---------------------------------------------------------------------------
# crown axes
# ---------------------------------------------------------------------------
def crown_axes(xy):
    """Major/minor axis of the crown polygon (convex hull of x-y points).

    Major = the longest chord between hull vertices. Minor = the widest extent
    of the hull measured perpendicular to the major direction. Returns
    {major:{length,orientation_deg}, minor:{length,orientation_deg}} or None.
    """
    xy = np.asarray(xy, dtype=float)
    if len(xy) < 3:
        return None
    try:
        from scipy.spatial import ConvexHull
        hull = ConvexHull(xy)
        pts = xy[hull.vertices]
    except Exception:
        pts = xy
    if len(pts) < 2:
        return None
    # major axis: max pairwise distance among hull vertices
    d2 = ((pts[:, None, :] - pts[None, :, :]) ** 2).sum(-1)
    i, j = np.unravel_index(np.argmax(d2), d2.shape)
    p_i, p_j = pts[i], pts[j]
    major_len = float(np.hypot(*(p_j - p_i)))
    major_dir = (p_j - p_i)
    major_ang = float(np.degrees(np.arctan2(major_dir[1], major_dir[0])))
    # minor axis: extent of hull perpendicular to major direction
    u = major_dir / (np.linalg.norm(major_dir) + 1e-12)
    perp = np.array([-u[1], u[0]])
    proj = pts @ perp
    minor_len = float(proj.max() - proj.min())
    minor_ang = float(np.degrees(np.arctan2(perp[1], perp[0])))
    # normalise angles to [0,180)
    major_ang %= 180.0
    minor_ang %= 180.0
    return {
        "major_axis": {"length_m": major_len, "orientation_deg": major_ang},
        "minor_axis": {"length_m": minor_len, "orientation_deg": minor_ang},
    }


# ---------------------------------------------------------------------------
# io helpers
# ---------------------------------------------------------------------------
def _load(path):
    return json.load(open(path)) if os.path.exists(path) else None


def _bbox(xyz):
    if len(xyz) == 0:
        return None
    lo = xyz.min(0); hi = xyz.max(0)
    return {"x_min": float(lo[0]), "y_min": float(lo[1]), "z_min": float(lo[2]),
            "x_max": float(hi[0]), "y_max": float(hi[1]), "z_max": float(hi[2])}


def _dbh(skel):
    """Diameter at breast height from the taper model (metres)."""
    base_r = float(skel.get("base_radius_m", 0.0))
    h_t = float(skel.get("height_m", 0.0))
    return 2.0 * _taper_radius(base_r, h_t, BREAST_HEIGHT_M)


def _skeletons(plot_dir, plot_uid):
    """tree_id -> (uid, skeleton dict) for every tree."""
    out = {}
    d = os.path.join(plot_dir, layout.INDIVIDUAL_TREES)
    if not os.path.isdir(d):
        return out
    for name in sorted(os.listdir(d)):
        tdir = os.path.join(d, name)
        if not os.path.isdir(tdir):
            continue
        uid = name.rsplit("_", 1)[-1]
        for cand in (os.path.join(tdir, "textual", f"{name}_trimmed_skeleton.json"),
                     os.path.join(tdir, f"{name}_trimmed_skeleton.json")):
            if os.path.exists(cand):
                sk = json.load(open(cand))
                # The skeleton's own tree_uid is authoritative; the dirname
                # suffix is only a fallback. Uids may contain underscores
                # (e.g. table TREE_UID values that are already full tags),
                # which the suffix parse would truncate.
                uid = str(sk.get("tree_uid") or name.rsplit("_", 1)[-1])
                out[int(sk.get("tree_id"))] = (uid, sk)
                break
    return out


# ---------------------------------------------------------------------------
# per-tree info
# ---------------------------------------------------------------------------
def _foliage_xyz(las_path):
    if not os.path.exists(las_path):
        return np.zeros((0, 3))
    xyz, cols = read_canonical_las(las_path, level="tree")
    return xyz[cols["semantic_code"].astype(int) == schema.SEM_FOLIAGE]


def build_tree_info(plot_dir, plot_uid, uid, skel, coord_system, units,
                    voxel_size, voxel_threshold):
    tag = layout.tree_tag(plot_uid, uid)
    raw_las = layout.tree_path(plot_dir, plot_uid, uid, "raw_las", create=False)
    trim_las = layout.tree_path(plot_dir, plot_uid, uid, "trimmed_las", create=False)

    # clouds / bbox
    xyz_all = np.zeros((0, 3))
    if os.path.exists(trim_las):
        xyz_all, _ = read_canonical_las(trim_las, level="tree")
    elif os.path.exists(raw_las):
        xyz_all, _ = read_canonical_las(raw_las, level="tree")

    fol_raw = _foliage_xyz(raw_las)
    fol_trim = _foliage_xyz(trim_las)

    # canopy hull volumes (raw/trimmed) from the hull json
    hull = _load(layout.tree_path(plot_dir, plot_uid, uid, "canopy_hull_json",
                                  create=False)) or {}
    hv = hull.get("volume_m3")
    if isinstance(hv, dict):
        hull_raw, hull_trim = hv.get("raw"), hv.get("trimmed")
    else:
        hull_raw = hull_trim = hv

    # crown axes from foliage x-y (prefer trimmed, else raw)
    crown_src = fol_trim if len(fol_trim) else fol_raw
    axes = crown_axes(crown_src[:, :2]) if len(crown_src) else None

    info = {
        "meta": {
            "tree_uid": uid,
            "plot_uid": plot_uid,
            "coordinate_system": coord_system,
            "units": units,
            "bounding_box": _bbox(xyz_all),
        },
        "tree_characteristics": {
            "center": [float(skel.get("base_xy", [None, None])[0]),
                       float(skel.get("base_xy", [None, None])[1])],
            "dbh_m": _dbh(skel),
            "height_m": float(skel.get("height_m", 0.0)),
            "canopy": {
                "convex_hull_volume_m3": {"raw": hull_raw, "trimmed": hull_trim},
                **(axes or {"major_axis": None, "minor_axis": None}),
            },
            "volume_m3": {
                "trunk": {"raw": trunk_volume(skel, False),
                          "trimmed": trunk_volume(skel, True)},
                "branch": {"raw": branch_volume(skel, False),
                           "trimmed": branch_volume(skel, True)},
                "foliage": {"raw": voxel_volume(fol_raw, voxel_size, voxel_threshold),
                            "trimmed": voxel_volume(fol_trim, voxel_size, voxel_threshold)},
            },
        },
    }
    return info, xyz_all


def _load_plot_cloud(plot_dir, plot_uid):
    """Decode the plot-level cloud once for the whole build: the understory
    cloud when the understory stage ran, else the raw scene.

    Returns (xyz, cols, is_understory), or (None, None, False) when no cloud
    exists. Every consumer below -- per-tree and per-shrub neighborhoods,
    plot aggregates, grass stats -- used to decode the understory cloud
    itself, once per tree and once per shrub, which dominated this stage's
    runtime on real plots; they now share this single decode.
    """
    understory_las = layout.plot_path(plot_dir, plot_uid, "understory_las",
                                      create=False)
    if os.path.exists(understory_las):
        xyz, cols = read_canonical_las(understory_las, level="plot")
        return xyz, cols, True
    raw_las = layout.plot_path(plot_dir, plot_uid, "raw_las", create=False)
    if os.path.exists(raw_las):
        xyz, cols = read_canonical_las(raw_las, level="plot")
        return xyz, cols, False
    return None, None, False


def _neighborhood(plot_dir, plot_uid, uid, tree_bbox, plot_cloud,
                  voxel_size, voxel_threshold, tree_id=None):
    """Counts + volume proportions of other veg within the tree's bbox.

    All volumes here (target tree and every neighbour class) are occupied-voxel
    volumes measured from the understory cloud with the same voxel params, so the
    proportions are on one consistent basis. The target is this tree's own points
    inside its bbox.

    `plot_cloud` is the shared decode from _load_plot_cloud. Understory only:
    like the old existence check on the understory file, a plot with no
    understory cloud gets no tree neighborhoods.
    """
    xyz, cols, is_under = plot_cloud if plot_cloud else (None, None, False)
    if tree_bbox is None or not is_under:
        return None
    sem = cols["semantic_code"].astype(int)
    tid = cols["tree_id"].astype(int)
    sid = cols["shrub_id"].astype(int)
    bx = ((xyz[:, 0] >= tree_bbox["x_min"]) & (xyz[:, 0] <= tree_bbox["x_max"]) &
          (xyz[:, 1] >= tree_bbox["y_min"]) & (xyz[:, 1] <= tree_bbox["y_max"]))
    this_tid = int(tree_id if tree_id is not None else _tree_id_of(uid))
    target = bx & np.isin(sem, (1, 2, 3)) & (tid == this_tid)
    other_tree = bx & np.isin(sem, (1, 2, 3)) & (tid != this_tid) & (tid > 0)
    shrub = bx & (sem == schema.SEM_SHRUB)
    grass = bx & (sem == schema.SEM_GRASS)
    n_other_trees = len(np.unique(tid[other_tree])) if other_tree.any() else 0
    n_shrubs = len(np.unique(sid[shrub])) if shrub.any() else 0

    def vv(mask):
        return voxel_volume(xyz[mask], voxel_size, voxel_threshold)
    tv = vv(target)
    v_other = vv(other_tree)
    v_shrub = vv(shrub)
    v_grass = vv(grass)

    def prop(v):
        denom = tv + v
        return (tv / denom) if denom > 0 else None
    return {
        "n_other_trees_in_bbox": int(n_other_trees),
        "n_shrubs_in_bbox": int(n_shrubs),
        "target_voxel_volume_m3": tv,
        "proportion_target_vs_other_trees": prop(v_other),
        "proportion_target_vs_shrubs": prop(v_shrub),
        "proportion_target_vs_grass": prop(v_grass),
        "proportion_target_vs_all_other_veg": prop(v_other + v_shrub + v_grass),
    }


def _tree_id_of(uid):
    try:
        return int(uid)
    except ValueError:
        return -1


# ---------------------------------------------------------------------------
# plot info + orchestration
# ---------------------------------------------------------------------------
def _plot_shape(plot_dir, plot_uid):
    """Read the plot shape descriptor from the scene JSON."""
    scene_path = layout.plot_path(plot_dir, plot_uid, "scene", create=False)
    if not os.path.exists(scene_path):
        return None
    try:
        s = json.load(open(scene_path))
    except Exception:
        return None
    shape = s.get("plot_shape", "circle")
    if shape == "rectangle":
        return {"shape": "rectangle",
                "length_m": s.get("plot_length_m"),
                "width_m": s.get("plot_width_m")}
    return {"shape": "circle", "radius_m": s.get("plot_radius_m")}


def _plot_center(plot_dir, plot_uid):
    """Plot-center lat/lon/elev from the scene JSON, passed through unchanged
    from the FIA input (see fia_helper.preprocess_fia / assemble_plot.py's
    scene.json write). None for any field the FIA input didn't carry."""
    scene_path = layout.plot_path(plot_dir, plot_uid, "scene", create=False)
    if not os.path.exists(scene_path):
        return {"lat": None, "lon": None, "elev_m": None}
    try:
        s = json.load(open(scene_path))
    except Exception:
        return {"lat": None, "lon": None, "elev_m": None}
    return {"lat": s.get("plot_lat"), "lon": s.get("plot_lon"),
            "elev_m": s.get("plot_elev")}


def build_plot_info(plot_dir, plot_uid, skels, tree_infos, coord_system, units,
                    voxel_size, voxel_threshold, plot_cloud=None):
    xyz_all = plot_cloud[0] if plot_cloud else np.zeros((0, 3))

    # tree aggregates from the per-tree infos (single source of truth)
    n_trees = len(tree_infos)
    basal_area = 0.0
    tv = {"trunk": {"raw": 0.0, "trimmed": 0.0},
          "branch": {"raw": 0.0, "trimmed": 0.0},
          "foliage": {"raw": 0.0, "trimmed": 0.0},
          "canopy": {"raw": 0.0, "trimmed": 0.0}}
    for info in tree_infos.values():
        tc = info["tree_characteristics"]
        basal_area += np.pi * (tc["dbh_m"] / 2.0) ** 2
        for k in ("trunk", "branch", "foliage"):
            for v in ("raw", "trimmed"):
                tv[k][v] += tc["volume_m3"][k][v] or 0.0
        ch = tc["canopy"]["convex_hull_volume_m3"]
        tv["canopy"]["raw"] += ch["raw"] or 0.0
        tv["canopy"]["trimmed"] += ch["trimmed"] or 0.0

    total_vol = {v: tv["trunk"][v] + tv["branch"][v] + tv["foliage"][v]
                 for v in ("raw", "trimmed")}

    # shrubs / grass from the understory cloud + shrub index
    shrub_total = _shrub_totals(plot_dir, plot_uid, voxel_size, voxel_threshold)
    grass_vol, grass_cover = _grass_stats(plot_cloud,
                                          voxel_size, voxel_threshold)

    return {
        "meta": {
            "plot_uid": plot_uid,
            "coordinate_system": coord_system,
            "units": units,
            "plot_shape": _plot_shape(plot_dir, plot_uid),
            "plot_center": _plot_center(plot_dir, plot_uid),
            "bounding_box": _bbox(xyz_all),
        },
        "vegetation_characteristics": {
            "trees": {
                "n_trees": n_trees,
                "basal_area_m2": float(basal_area),
                "total_volume_m3": total_vol,
                "trunk_volume_m3": tv["trunk"],
                "branch_volume_m3": tv["branch"],
                "foliage_volume_m3": tv["foliage"],
                "canopy_volume_m3": tv["canopy"],
            },
            "shrubs": {
                "n_shrubs": shrub_total["n"],
                "total_volume_m3": shrub_total["volume"],
            },
            "grass": {
                "percent_cover": grass_cover,
                "total_volume_m3": grass_vol,
            },
        },
    }


def _shrub_totals(plot_dir, plot_uid, voxel_size, voxel_threshold):
    """Sum shrub volumes (all three methods) from the shrub index / volume JSONs."""
    idx = _load(os.path.join(plot_dir, layout.PLOT_PLACEMENT, "textual",
                             f"{plot_uid}_shrub_index.json"))
    total = {"primitive": 0.0, "mesh": 0.0, "voxel": 0.0}
    n = 0
    if idx:
        for s in idx.get("shrubs", []):
            n += 1
            v = s.get("volume_m3", {})
            for k in total:
                val = v.get(k)
                if val is not None and val == val:      # skip NaN
                    total[k] += val
    return {"n": n, "volume": total}


def _grass_stats(plot_cloud, voxel_size, voxel_threshold):
    """Grass voxel volume + realized percent cover.

    percent_cover is the *realized* ground cover: the fraction of plot footprint
    cells (at COVER_CELL_M) that contain any grass point. This is measured from
    the cloud, not the input `grass_cover` generation parameter -- the generator
    treats cover as a thinning weight modulated by the vigor field and canopy
    attenuation, so realized cover differs from the input and the measured value
    is the honest one.
    """
    COVER_CELL_M = 0.5
    if not plot_cloud or not plot_cloud[2]:
        return 0.0, None
    xyz, cols = plot_cloud[0], plot_cloud[1]
    grass = cols["semantic_code"].astype(int) == schema.SEM_GRASS
    gvol = voxel_volume(xyz[grass], voxel_size, voxel_threshold)
    cover = None
    if len(xyz):
        def cells(pts):
            k = np.floor(pts[:, :2] / COVER_CELL_M).astype(np.int64)
            return set(map(tuple, k))
        plot_cells = cells(xyz)
        gcells = cells(xyz[grass]) if grass.any() else set()
        if plot_cells:
            cover = 100.0 * len(gcells) / len(plot_cells)
    return gvol, cover


def build_info(csv_path, plot_dir=None, *, coord_system="local", units="meters",
               occupied_voxel_size=0.1, occupied_voxel_threshold=1, verbose=True):
    # plot_uid is the plot dir's basename when given (matches the pipeline), else
    # the csv stem. Both yield e.g. "plot_1_1".
    if plot_dir is not None:
        plot_uid = os.path.basename(os.path.normpath(plot_dir))
    else:
        plot_uid = os.path.splitext(os.path.basename(csv_path))[0]
        plot_dir = os.path.join(os.path.dirname(csv_path), plot_uid)

    skels = _skeletons(plot_dir, plot_uid)
    # one decode of the plot cloud for the whole build (see _load_plot_cloud)
    plot_cloud = _load_plot_cloud(plot_dir, plot_uid)

    # per-tree info
    tree_infos = {}
    for tid, (uid, sk) in skels.items():
        info, xyz_all = build_tree_info(
            plot_dir, plot_uid, uid, sk, coord_system, units,
            occupied_voxel_size, occupied_voxel_threshold)
        nb = _neighborhood(plot_dir, plot_uid, uid,
                           info["meta"]["bounding_box"], plot_cloud,
                           occupied_voxel_size, occupied_voxel_threshold,
                           tree_id=tid)
        info["neighborhood"] = nb
        tree_infos[tid] = info

        out = os.path.join(_tree_textual(plot_dir, plot_uid, uid),
                           f"{layout.tree_tag(plot_uid, uid)}_info.json")
        os.makedirs(os.path.dirname(out), exist_ok=True)
        with open(out, "w") as f:
            json.dump(info, f, indent=2)

    # plot info
    plot_info = build_plot_info(plot_dir, plot_uid, skels, tree_infos,
                                coord_system, units,
                                occupied_voxel_size, occupied_voxel_threshold,
                                plot_cloud)
    pout = os.path.join(plot_dir, layout.PLOT_PLACEMENT, "textual",
                        f"{plot_uid}_info.json")
    os.makedirs(os.path.dirname(pout), exist_ok=True)
    with open(pout, "w") as f:
        json.dump(plot_info, f, indent=2)

    n_shrub = build_shrub_info(plot_dir, plot_uid, coord_system, units,
                               occupied_voxel_size, occupied_voxel_threshold,
                               plot_cloud)

    if verbose:
        print(f"{plot_uid}: wrote plot info + {len(tree_infos)} tree info files"
              + (f" + {n_shrub} shrub info files" if n_shrub else ""))
    return pout


def _tree_textual(plot_dir, plot_uid, uid):
    return os.path.join(layout.tree_dir(plot_dir, plot_uid, uid), "textual")


def _shrub_textual(plot_dir, plot_uid, uid):
    return os.path.join(layout.shrub_dir(plot_dir, plot_uid, uid), "textual")


def _ellipsoid_dbh(height, width, base_z, bh=BREAST_HEIGHT_M):
    """Diameter of the shrub ellipsoid's horizontal cross-section at breast
    height (bh above the shrub base). 0 if the shrub is shorter than bh.

    The ellipsoid is centered at height/2 with vertical semi-axis h/2 and
    horizontal semi-axis w/2. Width at local height z (measured from base):
      w(z) = width * sqrt(1 - ((z - h/2)/(h/2))^2)   for 0 <= z <= h, else 0.
    """
    h = float(height)
    if h <= 0 or bh > h:
        return 0.0
    c = h / 2.0
    frac = 1.0 - ((bh - c) / c) ** 2
    if frac <= 0:
        return 0.0
    return float(width) * (frac ** 0.5)


def _shrub_geoms(plot_dir, plot_uid):
    """shrub_uid -> geometry dict from each shrub's volume JSON."""
    out = {}
    d = os.path.join(plot_dir, layout.INDIVIDUAL_SHRUBS)
    if not os.path.isdir(d):
        return out
    for name in sorted(os.listdir(d)):
        sdir = os.path.join(d, name)
        if not os.path.isdir(sdir):
            continue
        uid = name.rsplit("_", 1)[-1]
        vpath = layout.shrub_path(plot_dir, plot_uid, uid, "volume_label",
                                  create=False)
        if os.path.exists(vpath):
            out[uid] = json.load(open(vpath))
    return out


def _shrub_neighborhood(plot_cloud, this_sid, bbox, voxel_size, voxel_threshold):
    """Counts + voxel-volume proportions of other veg within the shrub's bbox."""
    if bbox is None or not plot_cloud or not plot_cloud[2]:
        return None
    xyz, cols = plot_cloud[0], plot_cloud[1]
    sem = cols["semantic_code"].astype(int)
    tid = cols["tree_id"].astype(int)
    sid = cols["shrub_id"].astype(int)
    bx = ((xyz[:, 0] >= bbox["x_min"]) & (xyz[:, 0] <= bbox["x_max"]) &
          (xyz[:, 1] >= bbox["y_min"]) & (xyz[:, 1] <= bbox["y_max"]))
    target = bx & (sem == schema.SEM_SHRUB) & (sid == this_sid)
    trees = bx & np.isin(sem, (1, 2, 3)) & (tid > 0)
    other_shrub = bx & (sem == schema.SEM_SHRUB) & (sid != this_sid) & (sid > 0)
    grass = bx & (sem == schema.SEM_GRASS)
    n_trees = len(np.unique(tid[trees])) if trees.any() else 0
    n_other_shrubs = len(np.unique(sid[other_shrub])) if other_shrub.any() else 0

    def vv(mask):
        return voxel_volume(xyz[mask], voxel_size, voxel_threshold)
    tv = vv(target)
    v_trees = vv(trees)
    v_shrub = vv(other_shrub)
    v_grass = vv(grass)

    def prop(v):
        denom = tv + v
        return (tv / denom) if denom > 0 else None
    return {
        "n_trees_in_bbox": int(n_trees),
        "n_other_shrubs_in_bbox": int(n_other_shrubs),
        "target_voxel_volume_m3": tv,
        "proportion_target_vs_other_shrubs": prop(v_shrub),
        "proportion_target_vs_trees": prop(v_trees),
        "proportion_target_vs_grass": prop(v_grass),
        "proportion_target_vs_all_other_veg": prop(v_trees + v_shrub + v_grass),
    }


def build_shrub_info(plot_dir, plot_uid, coord_system, units,
                     voxel_size, voxel_threshold, plot_cloud=None):
    """Write per-shrub _info.json files (mirrors the tree info)."""
    geoms = _shrub_geoms(plot_dir, plot_uid)
    if not geoms:
        return 0
    # shrub points from the shared plot-cloud decode, for per-shrub bbox
    xyz = np.zeros((0, 3)); sid_all = np.zeros(0, dtype=int); sem_all = np.zeros(0, dtype=int)
    if plot_cloud and plot_cloud[2]:
        cols = plot_cloud[1]
        xyz = plot_cloud[0]
        sid_all = cols["shrub_id"].astype(int)
        sem_all = cols["semantic_code"].astype(int)

    written = 0
    for uid, g in geoms.items():
        s = int(g.get("shrub_id", _tree_id_of(uid)))
        mine = (sem_all == schema.SEM_SHRUB) & (sid_all == s)
        bbox = _bbox(xyz[mine]) if mine.any() else None

        height = float(g.get("height_m", 0.0))
        width = float(g.get("width_m", 0.0))
        base_z = float(g.get("base_z_m", 0.0))
        cx, cy = g.get("center_xy", [None, None])

        info = {
            "meta": {
                "shrub_uid": uid,
                "plot_uid": plot_uid,
                "coordinate_system": coord_system,
                "units": units,
                "bounding_box": bbox,
            },
            "shrub_characteristics": {
                "center": [cx, cy],
                "dbh_m": _ellipsoid_dbh(height, width, base_z),
                "max_width_m": width,
                "height_of_max_width_m": base_z + height / 2.0,
                "height_m": height,
                "volume_m3": {
                    "primitive": g.get("volume_m3", {}).get("primitive"),
                    "mesh": g.get("volume_m3", {}).get("mesh"),
                    "voxel": g.get("volume_m3", {}).get("voxel"),
                },
            },
            "neighborhood": _shrub_neighborhood(
                plot_cloud, s, bbox, voxel_size, voxel_threshold),
        }
        out = os.path.join(_shrub_textual(plot_dir, plot_uid, uid),
                           f"{layout.shrub_tag(plot_uid, uid)}_info.json")
        os.makedirs(os.path.dirname(out), exist_ok=True)
        with open(out, "w") as f:
            json.dump(info, f, indent=2)
        written += 1
    return written


def main():
    ap = argparse.ArgumentParser(description="Build plot/tree info JSONs")
    ap.add_argument("csv", nargs="+")
    ap.add_argument("-d", "--plot-dir", default=None)
    ap.add_argument("--coord-system", default="local")
    ap.add_argument("--units", default="meters")
    ap.add_argument("--occupied-voxel-size", type=float, default=0.1)
    ap.add_argument("--occupied-voxel-threshold", type=int, default=1)
    args = ap.parse_args()
    for c in args.csv:
        build_info(c, plot_dir=args.plot_dir, coord_system=args.coord_system,
                   units=args.units,
                   occupied_voxel_size=args.occupied_voxel_size,
                   occupied_voxel_threshold=args.occupied_voxel_threshold)


if __name__ == "__main__":
    main()
