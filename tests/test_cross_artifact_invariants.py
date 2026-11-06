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
  I. SCENE LINK AUDIT every relative path in the scene JSON resolves

Exit 0 = all checks passed.
"""
from __future__ import annotations

import argparse
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, REPO)

import numpy as np

import fixture
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
        bb = ti["meta"]["bounding_box"]
        check(f"{uid} info bbox == cloud bounds (1 mm)",
              np.allclose([bb["x_min"], bb["y_min"], bb["z_min"]],
                          xyz.min(axis=0), atol=1e-3)
              and np.allclose([bb["x_max"], bb["y_max"], bb["z_max"]],
                              xyz.max(axis=0), atol=1e-3))
        zmax_poly = max(float(np.asarray(pl, dtype=float)[:, 2].max())
                        for pl in sk["trunk"]["polylines"])
        check(f"{uid} trunk z-extent == height_m",
              abs(zmax_poly - sk["height_m"]) < 1e-9)
        check(f"{uid} cloud max z <= height + 1.5 m (foliage overhang)",
              float(xyz[:, 2].max()) <= sk["height_m"] + 1.5)
        r_bh = sk["base_radius_m"] * (1 - 1.37 / sk["height_m"])
        check(f"{uid} info dbh == 2 x taper radius at 1.37 m",
              abs(ti["tree_characteristics"]["dbh_m"]
                  - 2 * max(r_bh, 1e-3)) < 1e-9)
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
        reported = ti["tree_characteristics"]["volume_m3"]["trunk"]["trimmed"]
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
    ground = sem == 6
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
    under = np.isin(sem, [4, 5])
    check("shrub+grass HAG >= -1.5 mm (per-point draping)",
          float(hag[under].min()) >= -1.5e-3,
          f"min {hag[under].min():.4f}")
    if dtm is not None:
        gz = dtm.sample(uxyz[ground][:, 0], uxyz[ground][:, 1])
        check("ground on the DTM (2.5 mm = xy-quant x slope + z scale)",
              np.abs(gz - uxyz[ground][:, 2]).max() < 2.5e-3)

    print("== F. understory domains ==")
    grass = sem == 5
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

    print("== I. scene link audit ==")
    links = []
    leaf_paths(scene, links)
    missing = [p for p in links
               if not os.path.exists(os.path.join(fixroot, p))]
    check(f"every scene-JSON path resolves ({len(links)} links)",
          not missing, str(missing[:5]))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--keep", action="store_true",
                    help="keep the fixture directory for debugging")
    a = ap.parse_args()
    sys.exit(main(keep=a.keep))
