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
Each volume is a single post-trim value (spec 3.2): trimmed elements are
physically absent from the skeleton, so sums over the skeleton / inner
geometry are already over the surviving structure; clouds are the
surviving points.

Fields follow notes.md (spec 8.4): plot info carries NUM_TREES,
BASAL_AREA, NUM_SHRUBS, VOL_TREES/TRUNK/BRANCH/FOLIAGE/SHRUB/GRASS,
PLOT_SHAPE, PLOT_AREA, BOUNDING_BOX (6 floats), CHM_PATH/STEM_PATH and
PLAUSIBILITREE_CONFIGURATION at the top level; tree info carries
TREE_UID/DBH/HT/X/Y/Z and VOL_TREE/TRUNK/BRANCH/FOLIAGE. Foliage volume
voxelizes at FOLIAGE_VOXEL, grass at GRASS_VOXEL (separate knobs,
shared occupancy threshold); extension fields live in meta/canopy/
grass/shrubs/neighborhood blocks.

Crown major axis = the longest chord across the crown hull (x-y); minor axis =
the widest span perpendicular to it. Lengths in m, orientation in degrees CCW
from +x (local frame).
"""

from __future__ import annotations

import argparse
import json
import os

import numpy as np

import georef
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


def trunk_volume(skel):
    """Frustum sum along each trunk polyline using the taper radius per node.

    A frustum (truncated cone) between consecutive nodes at heights h0,h1 with
    radii r0,r1 and slant length L has volume (pi/3) L (r0^2 + r0 r1 + r1^2).
    Uses 3D segment length between nodes.
    """
    base_r = float(skel.get("base_radius_m", 0.0))
    h_t = float(skel.get("height_m", 0.0))
    total = 0.0
    for pl in skel.get("trunk", {}).get("polylines", []):
        arr = np.asarray(pl, dtype=float)
        if len(arr) < 2:
            continue
        for k in range(len(arr) - 1):
            p0, p1 = arr[k], arr[k + 1]
            L = float(np.linalg.norm(p1 - p0))
            r0 = _taper_radius(base_r, h_t, p0[2])
            r1 = _taper_radius(base_r, h_t, p1[2])
            total += (np.pi / 3.0) * L * (r0 * r0 + r0 * r1 + r1 * r1)
    return total


def branch_volume(skel):
    """Cylinder per branch: pi * r^2 * length, r = taper radius at branch base."""
    base_r = float(skel.get("base_radius_m", 0.0))
    h_t = float(skel.get("height_m", 0.0))
    total = 0.0
    for b in skel.get("branches", []):
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
    """notes.md BOUNDING_BOX: [x_min, y_min, z_min, x_max, y_max, z_max]."""
    if len(xyz) == 0:
        return None
    lo = xyz.min(0)
    hi = xyz.max(0)
    return [float(lo[0]), float(lo[1]), float(lo[2]),
            float(hi[0]), float(hi[1]), float(hi[2])]


def _bbox_geo(bbox, anchor):
    """BOUNDING_BOX in EPSG:4326/5773 for an anchored plot (spec Section 10):
    [lon_min, lat_min, elev_min, lon_max, lat_max, elev_max], elev = the
    anchor's plot-center height + local z (the anchor height defaults to
    0 when the input carried no elevation)."""
    if bbox is None:
        return None
    lon0, lat0 = georef.local_to_deg(anchor, bbox[0], bbox[1])
    lon1, lat1 = georef.local_to_deg(anchor, bbox[3], bbox[4])
    ez = float(anchor.get("elev") or 0.0)
    return [lon0, lat0, ez + bbox[2], lon1, lat1, ez + bbox[5]]


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


def _tree_centroid_z(plot_cloud, tree_id, xyz_tree):
    """TREE_Z: centroid height above the tree's base (spec 8.4).

    Primary: mean height-above-ground of the tree's points in the
    understory cloud. Fallback (no understory cloud, or the tree has no
    surviving points there): centroid z relative to the lowest point of
    the tree's own cloud. None when both are unavailable.
    """
    if plot_cloud and plot_cloud[2]:
        cols = plot_cloud[1]
        m = cols["tree_id"].astype(int) == int(tree_id)
        if m.any():
            hag = cols["height_above_ground"].astype(float)
            return float(hag[m].mean())
    if len(xyz_tree):
        return float(xyz_tree[:, 2].mean() - xyz_tree[:, 2].min())
    return None


def build_tree_info(plot_dir, plot_uid, uid, skel, coord_system, units,
                    foliage_voxel, voxel_threshold, tree_id=None,
                    plot_cloud=None, anchor=None):
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
    # single post-trim foliage source (spec 3.2): the surviving points
    fol_src = fol_trim if len(fol_trim) else fol_raw

    # canopy hull volume from the hull json (post-trim single value)
    hull = _load(layout.tree_path(plot_dir, plot_uid, uid, "canopy_hull_json",
                                  create=False)) or {}
    hull_vol = hull.get("volume_m3")

    # crown axes from surviving foliage x-y
    crown_src = fol_src
    axes = crown_axes(crown_src[:, :2]) if len(crown_src) else None

    base_xy = [float(skel.get("base_xy", [None, None])[0]),
               float(skel.get("base_xy", [None, None])[1])]
    vol_trunk = trunk_volume(skel)
    vol_branch = branch_volume(skel)
    vol_foliage = voxel_volume(fol_src, foliage_voxel, voxel_threshold)

    bbox = _bbox(xyz_all)
    if anchor is not None:
        bbox = _bbox_geo(bbox, anchor) if bbox is not None else None
        base_z = float(xyz_all[:, 2].min()) if len(xyz_all) else 0.0
        lon, lat = georef.local_to_deg(anchor, base_xy[0], base_xy[1])
        tree_center_real = [lon, lat,
                            float(anchor.get("elev") or 0.0) + base_z]
    else:
        tree_center_real = None

    info = {
        # notes.md required fields (spec 8.4)
        "TREE_UID": uid,
        "TREE_DBH": _dbh(skel),
        "TREE_HT": float(skel.get("height_m", 0.0)),
        "TREE_X": base_xy[0],
        "TREE_Y": base_xy[1],
        "TREE_Z": _tree_centroid_z(plot_cloud, tree_id, xyz_all),
        "VOL_TREE": vol_trunk + vol_branch + vol_foliage,
        "VOL_TRUNK": vol_trunk,
        "VOL_BRANCH": vol_branch,
        "VOL_FOLIAGE": vol_foliage,
        "BOUNDING_BOX": bbox,
        # extensions
        "meta": {
            "plot_uid": plot_uid,
            "coordinate_system": coord_system,
            "units": units,
        },
        "canopy": {
            "convex_hull_volume_m3": hull_vol,
            **(axes or {"major_axis": None, "minor_axis": None}),
        },
    }
    if tree_center_real is not None:
        info["TREE_CENTER_REAL"] = tree_center_real
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
    volumes measured from the understory cloud with the same voxel params
    (FOLIAGE_VOXEL / shared threshold), so the proportions are on one
    consistent basis. The target is this tree's own points inside its bbox.

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
    x0, y0, _, x1, y1, _ = tree_bbox
    bx = ((xyz[:, 0] >= x0) & (xyz[:, 0] <= x1) &
          (xyz[:, 1] >= y0) & (xyz[:, 1] <= y1))
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
def _scene(plot_dir, plot_uid):
    """Scene JSON dict, or {} when missing/unparsable."""
    scene_path = layout.plot_path(plot_dir, plot_uid, "scene", create=False)
    if not os.path.exists(scene_path):
        return {}
    try:
        return json.load(open(scene_path))
    except Exception:
        return {}


def _plot_area(scene):
    """Plot footprint area from the scene descriptor, or None when the
    geometry is unknown."""
    if scene.get("plot_shape") == "rectangle":
        length = scene.get("plot_length_m")
        width = scene.get("plot_width_m")
        if length is not None and width is not None:
            return float(length) * float(width)
        return None
    radius = scene.get("plot_radius_m")
    if radius is not None:
        return float(np.pi * radius * radius)
    return None


def _relpath_if_exists(path, plot_dir):
    """Path relative to the plot dir when the artifact exists, else None."""
    if path and os.path.exists(path):
        return os.path.relpath(path, plot_dir)
    return None


def build_plot_info(plot_dir, plot_uid, skels, tree_infos, coord_system, units,
                    foliage_voxel, grass_voxel, voxel_threshold,
                    plot_cloud=None, anchor=None):
    xyz_all = plot_cloud[0] if plot_cloud else np.zeros((0, 3))
    scene = _scene(plot_dir, plot_uid)

    # tree aggregates from the per-tree infos (single source of truth)
    n_trees = len(tree_infos)
    basal_area = 0.0
    tv = {"trunk": 0.0, "branch": 0.0, "foliage": 0.0, "canopy": 0.0}
    for info in tree_infos.values():
        basal_area += np.pi * (info["TREE_DBH"] / 2.0) ** 2
        for k in ("trunk", "branch", "foliage"):
            tv[k] += info["VOL_" + k.upper()]
        tv["canopy"] += info["canopy"]["convex_hull_volume_m3"] or 0.0

    vol_trees = tv["trunk"] + tv["branch"] + tv["foliage"]

    # shrubs / grass from the understory cloud + shrub volume JSONs
    shrub_total = _shrub_totals(plot_dir, plot_uid)
    grass_vol, grass_cover = _grass_stats(plot_cloud, grass_voxel,
                                          voxel_threshold)

    plot_center = None
    lon, lat, elev = (scene.get("plot_lon"), scene.get("plot_lat"),
                      scene.get("plot_elev"))
    if any(v is not None for v in (lon, lat, elev)):
        plot_center = [lon, lat, elev]

    # artifact paths (relative to the plot dir; None when not generated)
    chm_path = _relpath_if_exists(
        layout.plot_path(plot_dir, plot_uid, "canopy_feature_raw",
                         create=False), plot_dir)
    stem_path = _relpath_if_exists(
        layout.plot_path(plot_dir, plot_uid, "stem_feature_raw",
                         create=False), plot_dir)
    config_path = _relpath_if_exists(
        layout.plot_path(plot_dir, plot_uid, "scene", create=False), plot_dir)
    dtm_path = _relpath_if_exists(
        os.path.join(plot_dir, layout.PLOT_PLACEMENT, "textual",
                     f"{plot_uid}_dtm.asc"), plot_dir)

    info = {
        # notes.md required fields (spec 8.4)
        "NUM_TREES": n_trees,
        "BASAL_AREA": float(basal_area),
        "NUM_SHRUBS": shrub_total["n"],
        "VOL_TREES": vol_trees,
        "VOL_TRUNK": tv["trunk"],
        "VOL_BRANCH": tv["branch"],
        "VOL_FOLIAGE": tv["foliage"],
        "VOL_SHRUB": shrub_total["volume"],
        "VOL_GRASS": grass_vol,
        "PLOT_SHAPE": scene.get("plot_shape", "circle"),
        "PLOT_AREA": _plot_area(scene),
        "BOUNDING_BOX": (_bbox_geo(_bbox(xyz_all), anchor)
                         if anchor is not None else _bbox(xyz_all)),
        "CHM_PATH": chm_path,
        "STEM_PATH": stem_path,
        "PLAUSIBILITREE_CONFIGURATION": config_path,
        # extensions
        "meta": {
            "plot_uid": plot_uid,
            "coordinate_system": coord_system,
            "units": units,
        },
        "canopy": {"convex_hull_volume_m3": tv["canopy"]},
        "shrubs": {"primitive_reference_m3": shrub_total["primitive"]},
        "grass": {"percent_cover": grass_cover},
    }
    # optional notes.md fields: present only when the data exists
    if plot_center is not None:
        info["PLOT_CENTER"] = plot_center
    if dtm_path is not None:
        info["DTM_PATH"] = dtm_path
    return info


def _shrub_totals(plot_dir, plot_uid):
    """Sum shrub volumes from the per-shrub volume JSONs (sampled shape:
    mesh is the shrub volume per the notes.md ruling; primitive is the
    analytic NOMINAL reference)."""
    import glob as _glob
    pattern = os.path.join(plot_dir, "individual_shrubs",
                          f"{plot_uid}_shrub_*", "deep_learning", "volume",
                          "labels", "*_volume.json")
    n, mesh, prim = 0, 0.0, 0.0
    for p in sorted(_glob.glob(pattern)):
        doc = _load(p) or {}
        v = doc.get("volume_m3", {})
        n += 1
        mesh += float(v.get("mesh") or 0.0)
        prim += float(v.get("primitive") or 0.0)
    return {"n": n, "volume": mesh, "primitive": prim}


def _grass_stats(plot_cloud, grass_voxel, voxel_threshold):
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
    gvol = voxel_volume(xyz[grass], grass_voxel, voxel_threshold)
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
               foliage_voxel=0.1, grass_voxel=0.1, voxel_threshold=1,
               verbose=True):
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
    anchor = georef.load_anchor(plot_dir, plot_uid)

    # per-tree info
    tree_infos = {}
    for tid, (uid, sk) in skels.items():
        info, xyz_all = build_tree_info(
            plot_dir, plot_uid, uid, sk, coord_system, units,
            foliage_voxel, voxel_threshold, tree_id=tid, plot_cloud=plot_cloud,
            anchor=anchor)
        nb = _neighborhood(plot_dir, plot_uid, uid,
                           info["BOUNDING_BOX"], plot_cloud,
                           foliage_voxel, voxel_threshold,
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
                                foliage_voxel, grass_voxel, voxel_threshold,
                                plot_cloud, anchor=anchor)
    pout = os.path.join(plot_dir, layout.PLOT_PLACEMENT, "textual",
                        f"{plot_uid}_info.json")
    os.makedirs(os.path.dirname(pout), exist_ok=True)
    with open(pout, "w") as f:
        json.dump(plot_info, f, indent=2)

    n_shrub = build_shrub_info(plot_dir, plot_uid, coord_system, units,
                               foliage_voxel, voxel_threshold,
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
    x0, y0, _, x1, y1, _ = bbox
    bx = ((xyz[:, 0] >= x0) & (xyz[:, 0] <= x1) &
          (xyz[:, 1] >= y0) & (xyz[:, 1] <= y1))
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
    ap.add_argument("--foliage-voxel", type=float, default=0.1,
                    help="FOLIAGE_VOXEL: m per side for foliage voxel volumes")
    ap.add_argument("--grass-voxel", type=float, default=0.1,
                    help="GRASS_VOXEL: m per side for grass voxel volumes")
    ap.add_argument("--voxel-threshold", type=int, default=1,
                    help="points per occupied voxel for all voxel volumes")
    args = ap.parse_args()
    for c in args.csv:
        build_info(c, plot_dir=args.plot_dir, coord_system=args.coord_system,
                   units=args.units,
                   foliage_voxel=args.foliage_voxel,
                   grass_voxel=args.grass_voxel,
                   voxel_threshold=args.voxel_threshold)


if __name__ == "__main__":
    main()
