"""Cross-artifact invariant suite (whole-repo review round).

Run: conda run -n ssf python tests/test_cross_artifact_invariants.py [--keep]

Builds a full-pipeline fixture (6 trees, seed 42, 64 pts/m, generated
DTM + understory) in a temp dir, then recomputes every artifact's key
quantities from a DIFFERENT artifact than the one that produced it --
the pattern that caught the place() frame bug, the draping anchor bug
and the shrub center-anchoring bug (see VERIFICATION_LOG.md, W7 review
and whole-repo review entries):

  A. scene JSON records vs table (n_trees, rotations/reflections cover all)
  B. per-tree cloud <-> skeleton <-> info JSON (stem points hug the
     serialized centerline; info bbox == cloud bounds; trunk z-extent ==
     height_m; dbh == 2 x taper radius at 1.37 m; semantics domain)
  C. trunk frustum volume: independent recompute == info JSON
  D. branch_id encoding domain (spec Section 8.3/8.5)
  E. terrain/DTM invariants: spec 3.3 anchoring (lowest point == DTM at
     the tree's (X, Y)), shrub/grass HAG >= -1mm (per-point draping),
     ground on the DTM, zero noise
  F. understory domains and plot containment
  G. branch_code domain
  H. generation.json seed + realised records
  Hb. parameter dictionary: notes.md schema (A0/A1/TAPER + DBH/HT
     applicability bounds), withdrawn keys absent, bounds fall-through
     (W18)
  K. georeferencing (spec Section 10, W5): an anchored mini plot's UTM
     anchor, EPSG:4326 raster exports (features and instance masks),
     compound-CRS DTM GeoTIFF and geo info JSONs; the unanchored
     fixture stays CRS-less
  I. SCENE LINK AUDIT every relative path in the scene JSON resolves

Exit 0 = all checks passed.
"""
from __future__ import annotations

import argparse
import glob
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, REPO)

import numpy as np

import fixture
import layout
import terrain as tmod
from las import read_canonical_las
from tree_list import trees_frame

fails = []


def check(name, cond, msg=""):
    print(f"  {'PASS' if cond else 'FAIL'}  {name}"
          f"{('  -- ' + msg) if (msg and not cond) else ''}")
    if not cond:
        fails.append(name)


def frustum_volume(pl, base_r, h_t):
    """Independent frustum sum: pi/3 L (r0^2 + r0 r1 + r1^2) per segment."""
    tot = 0.0
    for k in range(len(pl) - 1):
        p0, p1 = pl[k], pl[k + 1]
        r0 = max(base_r * (1 - p0[2] / h_t), 1e-9) if h_t else 1e-9
        r1 = max(base_r * (1 - p1[2] / h_t), 1e-9) if h_t else 1e-9
        tot += np.pi / 3.0 * (r0 * r0 + r0 * r1 + r1 * r1) \
            * float(np.linalg.norm(p1 - p0))
    return tot


def leaf_paths(obj, out):
    if isinstance(obj, dict):
        for v in obj.values():
            leaf_paths(v, out)
    elif isinstance(obj, list):
        for v in obj:
            leaf_paths(v, out)
    elif (isinstance(obj, str) and "/" in obj
          and os.path.splitext(obj)[1] in (".laz", ".las", ".asc", ".json")):
        out.append(obj)


def main(keep: bool = False, ppm: int = 64):
    tmp, cleanup = fixture.workdir(keep=keep, tag="inv")
    table = fixture.build_tree_list(tmp)
    plot_dir = fixture.build_full_fixture(tmp, ppm=ppm)
    try:
        run_checks(tmp, plot_dir)
    finally:
        if not keep:
            cleanup()
    if fails:
        print(f"\nINVARIANT SUITE FAILED ({len(fails)}): {fails}")
        return 1
    print("\ncross-artifact invariant suite: all checks passed")
    return 0


def run_checks(fixroot, plot_dir):
    uid_tag = fixture.PLOT_UID
    place = os.path.join(plot_dir, "plot_placement")
    treed = lambda uid: os.path.join(plot_dir, "individual_trees",
                                     f"{uid_tag}_{uid}")
    scene = json.load(open(os.path.join(place, "textual",
                                        f"{uid_tag}_scene.json")))
    table_path = os.path.join(fixroot, f"{uid_tag}.json")
    table = trees_frame(table_path)

    print("== A. scene records vs table ==")
    check("scene n_trees == table rows", scene["n_trees"] == len(table))
    check("rotations and reflections cover every tree",
          set(scene["rotations_deg"]) == set(scene["reflections"])
          == set(table["TREE_UID"].astype(str)))
    for _, row in table.iterrows():
        uid = str(row["TREE_UID"])
        sk = json.load(open(os.path.join(treed(uid), "textual",
                                         f"{uid_tag}_{uid}_trimmed_skeleton.json")))
        check(f"{uid} skeleton base == table X/Y",
              np.allclose(sk["base_xy"], [row["X"], row["Y"]], atol=1e-9),
              f"{sk['base_xy']} vs ({row['X']}, {row['Y']})")

    print("== B. cloud <-> skeleton <-> info JSON ==")
    for _, row in table.iterrows():
        uid = str(row["TREE_UID"])
        d = treed(uid)
        xyz, cols = read_canonical_las(
            os.path.join(d, "las", f"{uid_tag}_{uid}_trimmed.laz"),
            level="tree")
        sk = json.load(open(os.path.join(d, "textual",
                                         f"{uid_tag}_{uid}_trimmed_skeleton.json")))
        ti = json.load(open(os.path.join(d, "textual",
                                         f"{uid_tag}_{uid}_info.json")))
        stem = cols["semantic_code"] == 1
        pts = xyz[stem][:, :2]
        dmin = np.full(len(pts), np.inf)
        for pl in sk["trunk"]["polylines"]:
            pl = np.asarray(pl, dtype=float)[:, :2]
            for k in range(len(pl) - 1):
                a_, b_ = pl[k], pl[k + 1]
                dv = b_ - a_
                dd = float(dv @ dv) or 1e-12
                t = np.clip(((pts - a_) @ dv) / dd, 0, 1)
                dmin = np.minimum(dmin,
                                  np.hypot(*(pts - (a_ + t[:, None] * dv)).T))
        check(f"{uid} stem points hug the skeleton (max {dmin.max():.3f} m)",
              dmin.max() < 1.0)
        bb = ti["BOUNDING_BOX"]
        check(f"{uid} info bbox == cloud bounds (1 mm)",
              np.allclose(bb[:3], xyz.min(axis=0), atol=1e-3)
              and np.allclose(bb[3:], xyz.max(axis=0), atol=1e-3))
        zmax_poly = max(float(np.asarray(pl, dtype=float)[:, 2].max())
                        for pl in sk["trunk"]["polylines"])
        check(f"{uid} trunk z-extent == height_m",
              abs(zmax_poly - sk["height_m"]) < 1e-9)
        check(f"{uid} cloud max z <= height + 1.5 m (foliage overhang)",
              float(xyz[:, 2].max()) <= sk["height_m"] + 1.5)
        r_bh = sk["base_radius_m"] * (1 - 1.37 / sk["height_m"])
        check(f"{uid} info dbh == 2 x taper radius at 1.37 m",
              abs(ti["TREE_DBH"] - 2 * max(r_bh, 1e-3)) < 1e-9)
        check(f"{uid} tree-cloud semantics in 1/2/3",
              set(np.unique(cols["semantic_code"])) <= {1, 2, 3})

    print("== C. trunk frustum volume vs independent recompute ==")
    for _, row in table.iterrows():
        uid = str(row["TREE_UID"])
        d = treed(uid)
        sk = json.load(open(os.path.join(d, "textual",
                                         f"{uid_tag}_{uid}_trimmed_skeleton.json")))
        ti = json.load(open(os.path.join(d, "textual",
                                         f"{uid_tag}_{uid}_info.json")))
        vol = sum(frustum_volume(np.asarray(pl, dtype=float),
                                 sk["base_radius_m"], sk["height_m"])
                  for pl in sk["trunk"]["polylines"])
        reported = ti["VOL_TRUNK"]
        check(f"{uid} trunk frustum volume {vol:.3f} == info {reported:.3f}",
              abs(vol - reported) < 1e-3)

    print("== D. branch_id encoding ==")
    xyz, cols = read_canonical_las(
        os.path.join(place, "las", f"{uid_tag}_trimmed.laz"), level="plot")
    bid = cols.get("branch_id")
    check("branch_id dim present", bid is not None)
    if bid is not None:
        check("branch_id fits int32", bid.max() < 2 ** 31)

    print("== E. terrain / DTM invariants (spec 3.3, zero-noise run) ==")
    uxyz, ucols = read_canonical_las(
        os.path.join(place, "las", f"{uid_tag}_understory.laz"), level="plot")
    sem = ucols["semantic_code"]
    hag = ucols["height_above_ground"].astype(float)
    ground = sem == 0
    check("understory cloud finite",
          np.all(np.isfinite(uxyz)) and np.all(np.isfinite(hag)))
    check("ground HAG == 0 exactly (this run draws no ground noise)",
          np.all(hag[ground] == 0.0))
    dtm_path = os.path.join(place, "textual", f"{uid_tag}_dtm.asc")
    check("plot-level DTM asc exists", os.path.exists(dtm_path))
    dtm = tmod.DTM.read_asc(dtm_path) if os.path.exists(dtm_path) else None
    tid = ucols["tree_id"].astype(int)
    for _, row in table.iterrows():
        uid = str(row["TREE_UID"])
        m = (tid == int(uid)) & (sem != 6)
        if not m.any() or dtm is None:
            continue
        zmin = float(uxyz[m][:, 2].min())
        gz = float(dtm.sample([row["X"]], [row["Y"]])[0])
        check(f"tree {uid} lowest point on DTM at its (X,Y)",
              abs(zmin - gz) < 1.5e-3, f"{zmin:.4f} vs {gz:.4f}")
    under = np.isin(sem, [8, 9])
    check("shrub+grass HAG >= -1.5 mm (per-point draping)",
          float(hag[under].min()) >= -1.5e-3,
          f"min {hag[under].min():.4f}")
    if dtm is not None:
        gz = dtm.sample(uxyz[ground][:, 0], uxyz[ground][:, 1])
        check("ground on the DTM (2.5 mm = xy-quant x slope + z scale)",
              np.abs(gz - uxyz[ground][:, 2]).max() < 2.5e-3)

    print("== F. understory domains (notes.md semantic codes, W14) ==")
    check("understory semantic domain is the notes.md set",
          set(np.unique(sem)).issubset({0, 1, 2, 3, 8, 9}),
          str(sorted(set(np.unique(sem)))))
    check("shrubs carry semantic 8", int((sem == 8).sum()) > 0)
    check("grass carries semantic 9", int((sem == 9).sum()) > 0)
    check("ground carries semantic 0", int((sem == 0).sum()) > 0)
    pi = json.load(open(os.path.join(
        place, "textual", f"{uid_tag}_info.json")))
    check("plot info carries the notes.md required keys (W16)",
          {"NUM_TREES", "BASAL_AREA", "NUM_SHRUBS", "VOL_TREES", "VOL_TRUNK",
           "VOL_BRANCH", "VOL_FOLIAGE", "VOL_SHRUB", "VOL_GRASS",
           "PLOT_SHAPE", "PLOT_AREA", "BOUNDING_BOX", "CHM_PATH",
           "STEM_PATH", "PLAUSIBILITREE_CONFIGURATION"} <= set(pi))
    check(f"plot info NUM_TREES == table rows ({pi['NUM_TREES']})",
          pi["NUM_TREES"] == len(table))
    tree_infos = [json.load(open(os.path.join(
        treed(str(r["TREE_UID"])), "textual",
        f"{uid_tag}_{r['TREE_UID']}_info.json")))
        for _, r in table.iterrows()]
    ba = sum(np.pi * (t["TREE_DBH"] / 2.0) ** 2 for t in tree_infos)
    check("plot info BASAL_AREA == sum of per-tree taper basal areas",
          abs(pi["BASAL_AREA"] - ba) < 1e-9)
    check("plot info VOL_TREES == sum of per-tree VOL_TREE",
          abs(pi["VOL_TREES"] - sum(t["VOL_TREE"] for t in tree_infos))
          < 1e-9)
    check("VOL_TREE == trunk + branch + foliage in every tree info",
          all(abs(t["VOL_TREE"] - t["VOL_TRUNK"] - t["VOL_BRANCH"]
                  - t["VOL_FOLIAGE"]) < 1e-12 for t in tree_infos))
    check(f"plot info counts its shrubs (n={pi['NUM_SHRUBS']}, "
          f"total {pi['VOL_SHRUB']:.2f} m3)",
          pi["NUM_SHRUBS"] > 0 and pi["VOL_SHRUB"] > 0.0)
    check("tree info carries the notes.md required keys (W16)",
          {"TREE_UID", "TREE_DBH", "TREE_HT", "TREE_X", "TREE_Y", "TREE_Z",
           "VOL_TREE", "VOL_TRUNK", "VOL_BRANCH", "VOL_FOLIAGE",
           "BOUNDING_BOX"}
          <= set(tree_infos[0]))
    shrub_mesh = sum(
        (json.load(open(fp)).get("volume_m3", {}) or {}).get("mesh") or 0.0
        for fp in glob.glob(os.path.join(
            plot_dir, "individual_shrubs", f"{uid_tag}_shrub_*",
            "deep_learning", "volume", "labels", "*_volume.json")))
    check("VOL_SHRUB == recomputed sum of per-shrub mesh volumes "
          "(sampled shape, never voxels)",
          abs(pi["VOL_SHRUB"] - shrub_mesh) < 1e-9)
    check("VOL_GRASS is a non-negative occupied-voxel volume",
          pi["VOL_GRASS"] >= 0.0)
    check("PLOT_AREA == pi r^2 for the circular fixture",
          abs(pi["PLOT_AREA"] - np.pi * scene["plot_radius_m"] ** 2) < 1e-9)
    for key in ("CHM_PATH", "STEM_PATH", "PLAUSIBILITREE_CONFIGURATION"):
        check(f"{key} points at an existing file",
              pi[key] is not None
              and os.path.exists(os.path.join(plot_dir, pi[key])))
    grass = sem == 9
    check("grass z above ground", float(uxyz[grass][:, 2].min()) > -0.01)
    check("all points within 1.5x plot radius",
          np.hypot(*uxyz[:, :2].T).max() <= 1.5 * scene["plot_radius_m"])

    print("== G. branch_code domain ==")
    for _, row in table.iterrows():
        uid = str(row["TREE_UID"])
        xyz, cols = read_canonical_las(
            os.path.join(treed(uid), "las", f"{uid_tag}_{uid}_trimmed.laz"),
            level="tree")
        bc = set(np.unique(cols.get("branch_code", np.array([]))))
        check(f"{uid} branch_code in 0..4", bc <= {0, 1, 2, 3, 4}, str(bc))

    print("== H. generation records ==")
    for _, row in table.iterrows():
        uid = str(row["TREE_UID"])
        gj = json.load(open(os.path.join(treed(uid), "textual",
                                         f"{uid_tag}_{uid}_generation.json")))
        ok = (gj.get("attempts")
              and all(a.get("seed") is not None for a in gj["attempts"])
              and "realised" in gj)
        check(f"{uid} generation.json carries per-attempt seeds + realised",
              ok)

    print("== Hb. parameter dictionary (notes.md schema, W18) ==")
    repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    fg = json.load(open(os.path.join(repo, "data", "functional_groups.json")))
    check("every dictionary entry carries A0/A1/TAPER + DBH/HT bounds",
          all({"A0", "A1", "TAPER", "DBH_MIN", "DBH_MAX",
               "HT_MIN", "HT_MAX"} <= set(e) for e in fg.values()))
    check("no withdrawn keys in the shipped dictionary (CW_*/TAPER_B/*_TOL)",
          not any({"CW_A0", "CW_A1", "TAPER_B", "DBH_TOL", "HT_TOL"}
                  & set(e) for e in fg.values()))
    check("every entry's bounds are ordered",
          all(e["DBH_MIN"] <= e["DBH_MAX"] and e["HT_MIN"] <= e["HT_MAX"]
              for e in fg.values()))
    import sst_generator as sst_mod
    sst = sst_mod.SSTGenerator()
    for _, row in table.iterrows():
        uid = str(row["TREE_UID"])
        entry = sst.entry_for(str(row["TREE_CODE"]), dbh=float(row["DBH"]),
                              ht=float(row["HT"]))
        check(f"{uid} resolved entry uses notes.md names",
              "A0" in entry and "TAPER" in entry
              and not {"CW_A0", "TAPER_B", "DBH_TOL", "HT_TOL"} & set(entry),
              str(row["TREE_CODE"]))
    check("in-bounds tree resolves to its own code entry",
          sst.entry_for("DF", dbh=0.5, ht=30.0)["DESCRIPTION"]
          == fg["DF"]["DESCRIPTION"])
    check("out-of-bounds tree falls through to the DEFAULT entry",
          sst.entry_for("DF", dbh=4.0, ht=30.0)["DESCRIPTION"]
          == fg["DEFAULT"]["DESCRIPTION"])
    check("DEFAULT always applies even outside its own bounds",
          sst.entry_for("999", dbh=9.9, ht=999.0)["DESCRIPTION"]
          == fg["DEFAULT"]["DESCRIPTION"])

    print("== I. scene link audit ==")
    links = []
    leaf_paths(scene, links)
    missing = [p for p in links
               if not os.path.exists(os.path.join(fixroot, p))]
    check(f"every scene-JSON path resolves ({len(links)} links)",
          not missing, str(missing[:5]))

    print("== J. element graph (notes.md Section 7.4, W15) ==")
    for _, row in table.iterrows():
        uid = str(row["TREE_UID"])
        sk = json.load(open(os.path.join(treed(uid), "textual",
                                         f"{uid_tag}_{uid}_raw_skeleton.json")))
        els = sk.get("elements")
        check(f"{uid} skeleton carries the element graph",
              isinstance(els, list) and len(els) > 0)
        if not els:
            continue
        uids = {e["BRANCH_UID"] for e in els}
        # every element has exactly one SOURCE (null only for the trunk
        # root) and DESTINATION is its exact inverse
        roots = [e for e in els if e["SOURCE"] is None]
        check(f"{uid} exactly one null SOURCE (the trunk root)",
              len(roots) == 1 and roots[0]["BRANCH_UID"] == "T00")
        check(f"{uid} every SOURCE resolves",
              all(e["SOURCE"] in uids for e in els if e["SOURCE"] is not None))
        by_uid = {e["BRANCH_UID"]: e for e in els}
        inv_ok = all(by_uid[c]["SOURCE"] == e["BRANCH_UID"]
                     for e in els for c in e["DESTINATION"])
        cov_ok = all(e["BRANCH_UID"] in by_uid[e["SOURCE"]]["DESTINATION"]
                     for e in els if e["SOURCE"] is not None)
        check(f"{uid} DESTINATION is the exact inverse of SOURCE",
              inv_ok and cov_ok)
        # partition: one order-0 element per trunk stem, one order-1 per
        # branch, lengths recomputed from the polyline geometry
        n_stems = len(sk["trunk"]["segments"])
        n_br = len(sk["branches"])
        check(f"{uid} element count == trunk segments + branches",
              len(els) == n_stems + n_br)
        check(f"{uid} trunk elements are order 0, branches order 1",
              all(e["BRANCH_ORDER"] == (0 if e["BRANCH_UID"].startswith("T")
                                        else 1) for e in els))
        seg_ok = vol_ok = rad_ok = True
        for i, seg in enumerate(sk["trunk"]["segments"]):
            e = by_uid[f"T{i:02d}"]
            poly = np.asarray(sk["trunk"]["polylines"][i], dtype=float)
            arc = float(np.linalg.norm(np.diff(poly, axis=0), axis=1).sum())
            seg_ok = seg_ok and abs(e["BRANCH_LENGTH"] - arc) < 1e-6
            r_exp = max(0.0, sk["base_radius_m"]
                        * (1.0 - seg["z_range"][0] / sk["height_m"]))
            rad_ok = rad_ok and abs(e["BRANCH_RADIUS"] - r_exp) < 1e-9
        for k, br in enumerate(sk["branches"]):
            e = by_uid[f"B{k:03d}"]
            d = float(np.linalg.norm(np.asarray(br["tip"])
                                     - np.asarray(br["base"])))
            seg_ok = seg_ok and abs(e["BRANCH_LENGTH"] - d) < 1e-6
        for e in els:
            v = math.pi * e["BRANCH_RADIUS"] ** 2 * e["BRANCH_LENGTH"]
            vol_ok = vol_ok and abs(e["BRANCH_VOLUME"] - v) <= 1e-9 * max(v, 1.0)
        check(f"{uid} element lengths recompute from the geometry", seg_ok)
        check(f"{uid} stem radii follow the taper profile", rad_ok)
        check(f"{uid} element volumes are pi r^2 L", vol_ok)
        # the trimmed variant carries the surviving subgraph (removed
        # elements are physically absent, spec 3.2): every element must
        # match its raw counterpart exactly and resolve within itself
        trm = json.load(open(os.path.join(
            treed(uid), "textual", f"{uid_tag}_{uid}_trimmed_skeleton.json")))
        tels = trm.get("elements")
        check(f"{uid} trimmed skeleton carries the surviving element graph",
              isinstance(tels, list) and len(tels) > 0)
        if tels:
            raw_map = {e["BRANCH_UID"]: e for e in els}
            check(f"{uid} trimmed graph matches the raw graph element-wise "
                  f"(DESTINATION may only shrink)",
                  all(raw_map.get(e["BRANCH_UID"]) and
                      all(raw_map[e["BRANCH_UID"]][k] == e[k] for k in
                          ("BRANCH_ORDER", "BRANCH_LENGTH", "BRANCH_RADIUS",
                           "BRANCH_VOLUME", "SOURCE")) and
                      set(e["DESTINATION"])
                      <= set(raw_map[e["BRANCH_UID"]]["DESTINATION"])
                      for e in tels))
            check(f"{uid} trimmed graph sources resolve internally",
                  all(e["SOURCE"] in {x["BRANCH_UID"] for x in tels}
                      or e["SOURCE"] is None for e in tels))
            check(f"{uid} trimmed graph destinations stay internal",
                  all(c in {x["BRANCH_UID"] for x in tels}
                      for e in tels for c in e["DESTINATION"]))

    print("== K. georeferencing (spec Section 10, W5) ==")
    import georef
    import rasterio
    # the unanchored main fixture stays CRS-less and unanchored
    with rasterio.open(os.path.join(place, "deep_learning", "stem_map",
                                    "features",
                                    f"{uid_tag}_stem_map_raw.tif")) as src:
        check("unanchored fixture raster carries no CRS", src.crs is None)
    check("unanchored fixture has no geo block",
          "geo" not in scene or scene["geo"] is None)
    check("unanchored fixture has no DTM GeoTIFF",
          not os.path.exists(layout.plot_path(fixroot, uid_tag,
                                              "dtm_tif", create=False)))
    # anchored mini plot through the same pipeline
    geowd, geo_cleanup = fixture.workdir(tag="inv_geo")
    gd = fixture.build_geo_fixture(geowd)
    try:
        gplace = os.path.join(gd, "plot_placement")
        gscene = json.load(open(os.path.join(gplace, "textual",
                                            "plot_geo_scene.json")))
        geo = gscene["geo"]
        check("scene geo block: epsg 32610 / zone 10N / hemisphere N",
              (geo["epsg"] == 32610 and geo["utm_zone"] == 10
               and geo["utm_hemisphere"] == "N"), str(geo))
        check("scene geo elev = plot-center EPSG:5773 height",
              geo["elev_m"] == 1500.0)
        anchor = georef.anchor_from_latlon(38.9, -120.6, 1500.0)
        e1, n1 = georef.local_to_utm(anchor, 0, 0)
        check("geo easting/northing == recomputed anchor",
              abs(geo["easting_m"] - e1) < 1e-6
              and abs(geo["northing_m"] - n1) < 1e-6)
        # raster exports carry EPSG:4326 with degree pixel sizes
        for key, cell in (("stem_feature_raw", 0.05),
                          ("canopy_feature_raw", 0.05)):
            with rasterio.open(layout.plot_path(gd, "plot_geo", key,
                                                create=False)) as src:
                check(f"{key} CRS is EPSG:4326", src.crs.to_epsg() == 4326)
                cx, cy = georef.deg_cell_size(cell, 38.9)
                check(f"{key} degree pixel sizes at plot latitude",
                      abs(src.transform.a - cx) < 1e-12
                      and abs(abs(src.transform.e) - cy) < 1e-12)
                west, north = src.transform.c, src.transform.f
                check(f"{key} corner is lon/lat near the plot center",
                      abs(west + 120.6) < 1e-3 and abs(north - 38.9) < 1e-3,
                      f"({west:.6f}, {north:.6f})")
        # plot-level instance masks: anchored like their feature rasters
        for lkey, fkey in (("stem_label_raw", "stem_feature_raw"),
                          ("canopy_label_raw", "canopy_feature_raw")):
            with (rasterio.open(layout.plot_path(gd, "plot_geo", lkey,
                                               create=False)) as lab,
                  rasterio.open(layout.plot_path(gd, "plot_geo", fkey,
                                               create=False)) as feat):
                arr = lab.read(1)
                check(f"{lkey} CRS is EPSG:4326", lab.crs.to_epsg() == 4326)
                drift = max(abs(lab.transform.a - feat.transform.a),
                            abs(lab.transform.e - feat.transform.e),
                            abs(lab.transform.c - feat.transform.c),
                            abs(lab.transform.f - feat.transform.f))
                check(f"{lkey} transform matches its feature raster (<1e-6 deg)",
                      drift < 1e-6, f"drift {drift:.2e}")
                check(f"{lkey} labels all three planted trees",
                      set(np.unique(arr).tolist()) >= {1, 2, 3},
                      str(np.unique(arr).tolist()))
        shrub_tifs = [rp for rp in glob.glob(
            os.path.join(gd, "individual_shrubs", "plot_geo_shrub_*",
                         "deep_learning", "**", "*.tif"), recursive=True)]
        check("anchored plot produced shrub raster exports",
              len(shrub_tifs) > 0, f"{len(shrub_tifs)} tifs")
        for rp in shrub_tifs[:4]:
            with rasterio.open(rp) as src:
                check(f"shrub raster CRS is EPSG:4326 ({os.path.basename(rp)})",
                      src.crs is not None and src.crs.to_epsg() == 4326)
        # DTM GeoTIFF: compound CRS + absolute heights
        dtif = layout.plot_path(gd, "plot_geo", "dtm_tif", create=False)
        check("DTM GeoTIFF written for the anchored plot",
              os.path.exists(dtif))
        with rasterio.open(dtif) as src:
            check("DTM CRS is compound EPSG:9707",
                  src.crs.to_epsg() == 9707, str(src.crs.to_epsg()))
            band = src.read(1)
            c = band[band != src.nodata]
            check("DTM band carries absolute heights near 1500 m",
                  1400 < float(c.min()) < 1600, f"min {float(c.min()):.1f}")
        # info JSONs: geo bbox, TREE_CENTER_REAL, PLOT_CENTER
        gpi = json.load(open(os.path.join(gplace, "textual",
                                         "plot_geo_info.json")))
        gbb = gpi["BOUNDING_BOX"]
        check("plot info bbox is lon/lat/elev",
              abs(gbb[0] + 120.6) < 1e-3 and abs(gbb[1] - 38.9) < 1e-3
              and gbb[2] > 1400, str(gbb[:3]))
        check("PLOT_CENTER echoes the tree list center",
              gpi["PLOT_CENTER"] == [-120.6, 38.9, 1500.0])
        gti = json.load(open(os.path.join(gd, "individual_trees",
                                          "plot_geo_001", "textual",
                                          "plot_geo_001_info.json")))
        tcr = gti["TREE_CENTER_REAL"]
        exp_lon, exp_lat = georef.local_to_deg(anchor, gti["TREE_X"],
                                              gti["TREE_Y"])
        check("tree info TREE_CENTER_REAL = plot center + base offset",
              abs(tcr[0] - exp_lon) < 1e-9 and abs(tcr[1] - exp_lat) < 1e-9
              and tcr[2] > 1400)
        # anchor math round-trips within 1 mm; limits of agreement raise
        for (x, y) in ((-12.0, -12.0), (12.0, 12.0), (0.0, 0.0)):
            lon, lat = georef.local_to_deg(anchor, x, y)
            a2 = georef.anchor_from_latlon(lat, lon)
            e2, n2 = georef.local_to_utm(a2, 0, 0)
            ee, nn = georef.local_to_utm(anchor, x, y)
            check(f"local->deg->utm round-trip {x:.0f},{y:.0f} (1 mm)",
                  abs(e2 - ee) < 1e-3 and abs(n2 - nn) < 1e-3)
        try:
            georef.anchor_from_latlon(95.0, 0.0)
            check("lat 95 raises (limits of agreement)", False)
        except ValueError:
            check("lat 95 raises (limits of agreement)", True)
    finally:
        geo_cleanup()


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--keep", action="store_true",
                    help="keep the fixture directory for debugging")
    a = ap.parse_args()
    sys.exit(main(keep=a.keep))
