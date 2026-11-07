"""
understory.py: add procedural shrubs and grass to a plot's terrain cloud.

Runs after terrain.py, before normalize.py. Reads the plot terrain cloud
(trees classified + ground, with height_above_ground / stem_foliage dims),
generates a shrub layer and a grass layer via vegetation.py, drapes them onto
the same DTM the terrain step used, and writes:

    plot_placement/las/<PLOT_UID>_understory.las
    individual_trees/<PLOT_UID_TREE_UID>/las/<..>_understory.las

Schema (all points): x, y, z, classification, tree_id, and extra dims
height_above_ground, stem_foliage (type: 1 stem/0 foliage/2 shrub/3 grass/255
ground), shrub_id. Trees keep tree_id, shrub_id=0; shrubs get tree_id=0,
shrub_id>=1; grass and ground get both=0. LAS classification stays authoritative
by height above ground (3 low <=0.5 m, 4 med 0.5-2 m, 5 high >2 m; ground 2).

DTM: reproduced per the terrain step's record in the scene JSON (generated with
the same params/seed, or re-read from the provided file) so the understory sits
on the identical surface.

Requires: numpy, scipy, laspy (via las.py), vegetation.py, terrain.py, layout.py
"""

from __future__ import annotations

import argparse
import json
import os

import numpy as np

import layout
from las import read_points_las, export_points_las, export_canonical_las, read_canonical_las
from terrain import DTM, classify_by_hag, CLASS_GROUND
import schema
import vegetation as veg


def _dtm_for_plot(scene, dtm_path, plot_bounds, seed):
    """Rebuild the DTM the terrain step used, for consistent draping."""
    terr = scene.get("terrain", {}) if scene else {}
    xmin, xmax, ymin, ymax = plot_bounds
    if terr.get("dtm_generated"):
        # bounds were auto-sized symmetric about (0,0) with one cell margin
        cell = terr.get("dtm_cell_m", 1.0)
        half_w = max(abs(xmin), abs(xmax)) + cell
        half_h = max(abs(ymin), abs(ymax)) + cell
        # relief/bumpiness parsed from the recorded description
        desc = terr.get("dtm", "")
        relief = _num_after(desc, "relief")
        bump = _num_after(desc, "bumpiness")
        d = DTM.generate(2 * half_w, 2 * half_h, cell,
                         bump if bump is not None else 0.5,
                         relief if relief is not None else 10.0, seed=seed)
        d.recenter_to(0.0, 0.0)
        return d
    if dtm_path and os.path.exists(dtm_path):
        d = DTM.read_asc(dtm_path)
        d.recenter_to(0.0, 0.0)
        return d
    return None


def _num_after(s, key):
    import re
    m = re.search(rf"{key}\s+([-\d.]+)", s)
    return float(m.group(1)) if m else None


def build_understory(csv_path, plot_dir=None, *, dtm_path=None, seed=0,
                     shrub_cover=0.15, shrub_cell_size=2.0,
                     shrub_height_range=(0.4, 1.8), shrub_width_range=(0.8, 1.6),
                     points_per_shrub=300, shrub_margin=0.5,
                     grass_cover=0.6, grass_height_range=(0.03, 0.35),
                     grass_patch_scale=5.0, grass_candidate_density=200.0,
                     canopy_attenuation=0.3,
                     shrub_outputs=True,
                     shrub_canopy_cell_size=0.25, shrub_stem_cell_size=0.25,
                     verbose=True):
    plot_uid = os.path.splitext(os.path.basename(csv_path))[0]
    if plot_dir is None:
        plot_dir = os.path.join(os.path.dirname(os.path.abspath(csv_path)), plot_uid)

    terr_las = layout.plot_path(plot_dir, plot_uid, "terrain_las", create=False)
    if not os.path.exists(terr_las):
        raise FileNotFoundError(
            f"{plot_uid}: terrain cloud not found ({terr_las}). Run terrain.py "
            f"before understory.py.")

    txyz, tcols = read_canonical_las(terr_las, level="plot")
    tcls = tcols["classification"].astype(np.uint8)
    ttid = tcols["tree_id"].astype(np.int32)
    thag = tcols["height_above_ground"].astype(np.float32)
    tsem = tcols["semantic_code"].astype(np.int8)       # object identity
    tcode = tcols["branch_code"].astype(np.int8)
    tbid = tcols["branch_id"].astype(np.int32)

    x_min, y_min = float(txyz[:, 0].min()), float(txyz[:, 1].min())
    x_max, y_max = float(txyz[:, 0].max()), float(txyz[:, 1].max())
    plot_bounds = (x_min, x_max, y_min, y_max)

    scene_path = layout.plot_path(plot_dir, plot_uid, "scene", create=False)
    scene = json.load(open(scene_path)) if os.path.exists(scene_path) else {}
    dtm = _dtm_for_plot(scene, dtm_path, plot_bounds, seed)

    # tree positions + canopy radii from the terrain cloud's tree points
    tree_positions, tree_canopy_r = _tree_footprints(txyz, ttid, tsem)

    # --- shrubs ---
    s_pts, s_ids, s_centers, s_radii, shrub_geom = veg.gen_shrub_layer(
        plot_bounds, shrub_cell_size, shrub_cover, shrub_height_range,
        shrub_width_range, tree_positions, tree_canopy_r,
        points_per_shrub=points_per_shrub, margin=shrub_margin,
        dtm=dtm, seed=seed)

    # --- grass ---
    g_pts = veg.gen_grass_layer(
        plot_bounds, grass_cover, grass_height_range, grass_patch_scale,
        tree_positions, tree_canopy_r,
        shrub_positions=s_centers if len(s_centers) else None,
        shrub_radii=s_radii if len(s_radii) else None,
        canopy_attenuation=canopy_attenuation,
        candidate_density=grass_candidate_density, dtm=dtm, seed=seed + 1)

    # --- clip understory vegetation to the plot boundary (circle only) ---
    # Placement already respects a rectangle, so this is a no-op there. For a
    # circular plot the shrub/grass fill covers the bounding square, so drop any
    # shrub or grass point outside the disc. Trees are already inside the plot.
    if scene.get("plot_shape", "circle") == "circle":
        r_plot = float(scene.get("plot_radius_m") or 0.0)
        if r_plot > 0:
            def _in_disc(pts):
                if not len(pts):
                    return np.zeros(0, dtype=bool)
                return (pts[:, 0] ** 2 + pts[:, 1] ** 2) <= r_plot ** 2
            if len(s_pts):
                keep = _in_disc(s_pts)
                s_pts, s_ids = s_pts[keep], s_ids[keep]
            if len(g_pts):
                g_pts = g_pts[_in_disc(g_pts)]

    # --- assemble new points with the canonical schema ---
    def _hag(pts):
        if dtm is None or not len(pts):
            return pts[:, 2] if len(pts) else np.zeros(0)
        return pts[:, 2] - dtm.sample(pts[:, 0], pts[:, 1])

    # new points as xyz + classification(height) + tree_id(0) + semantic + shrub_id
    new_xyz, new_hag, new_sem, new_cls, new_tid, new_sid = [], [], [], [], [], []
    if len(s_pts):
        hag = _hag(s_pts)
        new_xyz.append(s_pts)
        new_hag.append(hag)
        new_cls.append(classify_by_hag(hag))
        new_sem.append(np.full(len(s_pts), veg.TYPE_SHRUB, dtype=np.int8))   # 4
        new_tid.append(np.zeros(len(s_pts), dtype=np.int32))
        new_sid.append(s_ids.astype(np.int32))
    if len(g_pts):
        hag = _hag(g_pts)
        new_xyz.append(g_pts)
        new_hag.append(hag)
        new_cls.append(classify_by_hag(hag))
        new_sem.append(np.full(len(g_pts), veg.TYPE_GRASS, dtype=np.int8))   # 5
        new_tid.append(np.zeros(len(g_pts), dtype=np.int32))
        new_sid.append(np.zeros(len(g_pts), dtype=np.int32))

    # terrain points carry through unchanged (shrub_id defaults 0)
    terr_sid = np.zeros(len(txyz), dtype=np.int32)

    if new_xyz:
        add_xyz = np.vstack(new_xyz)
        all_xyz = np.vstack([txyz, add_xyz])
        all_hag = np.concatenate([thag, np.concatenate(new_hag)]).astype(np.float32)
        all_sem = np.concatenate([tsem, np.concatenate(new_sem)]).astype(np.int8)
        all_cls = np.concatenate([tcls, np.concatenate(new_cls)]).astype(np.uint8)
        all_tid = np.concatenate([ttid, np.concatenate(new_tid)]).astype(np.int32)
        all_sid = np.concatenate([terr_sid, np.concatenate(new_sid)]).astype(np.int32)
        # new understory points are not tree parts -> branch_code / branch_id 0
        all_code = np.concatenate(
            [tcode, np.zeros(len(add_xyz), dtype=np.int8)]).astype(np.int8)
        all_bid = np.concatenate(
            [tbid, np.zeros(len(add_xyz), dtype=np.int32)]).astype(np.int32)
    else:
        all_xyz, all_hag, all_sem = txyz, thag, tsem
        all_cls, all_tid, all_sid = tcls, ttid, terr_sid
        all_code, all_bid = tcode, tbid

    out = layout.plot_path(plot_dir, plot_uid, "understory_las")
    export_canonical_las(
        out, all_xyz,
        {"tree_id": all_tid, "shrub_id": all_sid, "semantic_code": all_sem,
         "classification": all_cls, "branch_code": all_code,
         "branch_id": all_bid,
         "height_above_ground": all_hag},
        level="plot")

    n_tree_files = _write_tree_understory(
        plot_dir, plot_uid, all_xyz, all_tid, all_sid, all_sem, all_cls,
        all_code, all_bid, all_hag)

    # per-shrub individual outputs (clouds, volume JSON, canopy + DBH maps)
    n_shrub_files = 0
    if shrub_outputs and shrub_geom:
        from shrub_outputs import build_shrub_outputs
        n_shrub_files = build_shrub_outputs(
            plot_dir, plot_uid, shrub_geom,
            canopy_cell_size=shrub_canopy_cell_size,
            stem_cell_size=shrub_stem_cell_size,
            verbose=verbose)

    # record in scene
    if os.path.exists(scene_path):
        scene["understory"] = {
            "shrub_cover": shrub_cover, "shrub_cell_size": shrub_cell_size,
            "shrub_height_range": list(shrub_height_range),
            "grass_cover": grass_cover,
            "grass_height_range": list(grass_height_range),
            "grass_patch_scale": grass_patch_scale,
            "canopy_attenuation": canopy_attenuation,
            "n_shrubs": int(len(s_centers)),
            "n_shrub_points": int(len(s_pts)),
            "n_grass_points": int(len(g_pts)),
            "semantic_note": ("shrubs -> semantic_code 4, grass -> 5, on top of "
                              "the existing tree(1/2/3)/ground(6) codes; shrub_id "
                              ">=1 on shrub points, 0 elsewhere"),
            "plot_understory_las": layout.plot_rel(plot_dir, plot_uid, "understory_las"),
        }
        with open(scene_path, "w") as f:
            json.dump(scene, f, indent=2)

    if verbose:
        print(f"{plot_uid}: understory +{len(s_pts):,} shrub ({len(s_centers)} "
              f"shrubs) +{len(g_pts):,} grass -> {len(all_xyz):,} total, "
              f"{n_tree_files} tree files")
    return out


def _tree_footprints(xyz, tid, tsem):
    """Per-tree (x,y) centre and canopy radius from the cloud's tree points."""
    is_tree = np.isin(tsem, (schema.SEM_TREE, schema.SEM_BRANCH, schema.SEM_FOLIAGE))
    positions, radii = [], []
    for t in np.unique(tid[is_tree & (tid > 0)]):
        m = is_tree & (tid == t)
        px, py = xyz[m, 0], xyz[m, 1]
        cxp, cyp = float(px.mean()), float(py.mean())
        r = float(np.percentile(np.hypot(px - cxp, py - cyp), 95))
        positions.append((cxp, cyp))
        radii.append(max(r, 0.1))
    return (np.asarray(positions) if positions else np.zeros((0, 2)),
            np.asarray(radii) if radii else np.zeros(0))


def _write_tree_understory(plot_dir, plot_uid, all_xyz, all_tid, all_sid,
                           all_sem, all_cls, all_code, all_bid, all_hag):
    """Per-tree understory crop: all points inside the tree's canopy x-y bbox."""
    scene_path = layout.plot_path(plot_dir, plot_uid, "scene", create=False)
    tree_uids = []
    if os.path.exists(scene_path):
        tree_uids = list(json.load(open(scene_path)).get("trees", {}).keys())
    if not tree_uids:
        d = os.path.join(plot_dir, layout.INDIVIDUAL_TREES)
        tree_uids = [n.rsplit("_", 1)[-1] for n in sorted(os.listdir(d))
                     if os.path.isdir(os.path.join(d, n))]

    x, y = all_xyz[:, 0], all_xyz[:, 1]
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
        out = layout.tree_path(plot_dir, plot_uid, uid, "understory_las")
        export_canonical_las(
            out, all_xyz[m],
            {"tree_id": all_tid[m], "shrub_id": all_sid[m],
             "semantic_code": all_sem[m], "classification": all_cls[m],
             "branch_code": all_code[m], "branch_id": all_bid[m],
             "height_above_ground": all_hag[m].astype(np.float32)},
            level="tree")
        n += 1
    return n


def main():
    ap = argparse.ArgumentParser(description="Add shrubs + grass to the terrain cloud")
    ap.add_argument("csv", nargs="+", help="plot table csv(s)")
    ap.add_argument("-d", "--plot-dir", default=None)
    ap.add_argument("--dtm", default=None,
                    help="provided DTM .asc (only if terrain used one)")
    ap.add_argument("--seed", type=int, default=0)
    # shrubs
    ap.add_argument("--shrub-cover", type=float, default=0.15)
    ap.add_argument("--shrub-cell-size", type=float, default=2.0)
    ap.add_argument("--shrub-height-min", type=float, default=0.4)
    ap.add_argument("--shrub-height-max", type=float, default=1.8)
    ap.add_argument("--shrub-width-min", type=float, default=0.8)
    ap.add_argument("--shrub-width-max", type=float, default=1.6)
    ap.add_argument("--points-per-shrub", type=int, default=300)
    ap.add_argument("--shrub-margin", type=float, default=0.5)
    # grass
    ap.add_argument("--grass-cover", type=float, default=0.6)
    ap.add_argument("--grass-height-min", type=float, default=0.03)
    ap.add_argument("--grass-height-max", type=float, default=0.35)
    ap.add_argument("--grass-patch-scale", type=float, default=5.0)
    ap.add_argument("--grass-candidate-density", type=float, default=200.0)
    ap.add_argument("--canopy-attenuation", type=float, default=0.3)
    # per-shrub individual outputs
    ap.add_argument("--no-shrub-outputs", action="store_true",
                    help="skip writing per-shrub individual outputs")
    ap.add_argument("--shrub-canopy-cell-size", type=float, default=0.25)
    ap.add_argument("--shrub-stem-cell-size", type=float, default=0.25)
    args = ap.parse_args()

    for path in args.csv:
        build_understory(
            path, plot_dir=args.plot_dir, dtm_path=args.dtm, seed=args.seed,
            shrub_cover=args.shrub_cover, shrub_cell_size=args.shrub_cell_size,
            shrub_height_range=(args.shrub_height_min, args.shrub_height_max),
            shrub_width_range=(args.shrub_width_min, args.shrub_width_max),
            points_per_shrub=args.points_per_shrub, shrub_margin=args.shrub_margin,
            grass_cover=args.grass_cover,
            grass_height_range=(args.grass_height_min, args.grass_height_max),
            grass_patch_scale=args.grass_patch_scale,
            grass_candidate_density=args.grass_candidate_density,
            canopy_attenuation=args.canopy_attenuation,
            shrub_outputs=not args.no_shrub_outputs,
            shrub_canopy_cell_size=args.shrub_canopy_cell_size,
            shrub_stem_cell_size=args.shrub_stem_cell_size)


if __name__ == "__main__":
    main()
