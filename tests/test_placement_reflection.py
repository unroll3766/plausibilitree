"""Placement & frame regression suite (W7 + the W7-review frame fixes).

Run: conda run -n ssf python tests/test_placement_reflection.py [--keep]

Builds its own fixture (FIA sample -> 6 trees, seed 42, 64 pts/m) in a
temp dir and checks, in order:

  0. FRAME CONSISTENCY (regression guard for the place() polyline bug
     found by the W7 review): raw skeleton polylines anchored at the
     global base_xy, trunk geometry within crown reach of the base.
  1. Same-seed determinism of assembly (bit-identical trimmed cloud,
     identical rotation/reflection records).
  2. Scene JSON records: reflections cover every placed tree; both
     reflection states occur in the fixture.
  3. CLOUD-VS-RECORDS parity: replaying the recorded rotation + mirror
     on each per-tree raw cloud reproduces the merged raw cloud (up to
     the merged cloud's 1 mm re-quantization at its own offset).
  4. SKELETON-VS-RECORDS parity: replaying the recorded 4dp rotation +
     mirror on the raw skeleton reproduces the trimmed skeleton's
     polyline and branch geometry (the residual is the 4dp angle
     record itself: 5e-5 deg ~ 8.8e-7 rad reaches ~1.6e-5 m at an 18 m
     crown; a real transform error is metres).
  5. Reflection draw distribution across seeds (both states, ~50%).
  6. CHIRALITY: signed area over ordered branch bases flips sign
     exactly on reflected trees and is preserved on the others -- the
     mirror is a true handedness flip, not a substitute rotation.

Exit 0 = all checks passed. See VERIFICATION_LOG.md (W7 + W7 review
entries) for the history behind every check.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, REPO)

import numpy as np

import assemble_plot as ap
import fixture
from las import read_canonical_las

fails = []


def check(name, cond, msg=""):
    print(f"  {'PASS' if cond else 'FAIL'}  {name}"
          f"{('  -- ' + msg) if (msg and not cond) else ''}")
    if not cond:
        fails.append(name)


def manual_transform(xy, angle_rad, reflected, cx, cy):
    out = ap.rotate_z(np.asarray(xy, dtype=float), angle_rad, cx, cy)
    if reflected:
        out[..., 0] = 2.0 * cx - out[..., 0]
    return out


def signed_area(pts):
    x, y = pts[:, 0], pts[:, 1]
    return 0.5 * float(np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y))


def main(keep: bool = False):
    tmp, cleanup = fixture.workdir(keep=keep, tag="place")
    try:
        table = fixture.build_tree_list(tmp)
        plot_dir = fixture.build_generated_plot(table)
    except Exception:
        cleanup()
        raise
    try:
        run_checks(tmp, table, plot_dir)
    finally:
        if not keep:
            cleanup()
    if fails:
        print(f"\nPLACEMENT/FRAME SUITE FAILED ({len(fails)}): {fails}")
        return 1
    print("\nplacement & frame suite: all checks passed")
    return 0


def run_checks(tmp, table, plot_dir):
    treedir = os.path.join(plot_dir, "individual_trees")

    print("== 0. frame consistency (regression guard) ==")
    for uid in ["001", "006"]:
        sk = json.load(open(os.path.join(treedir, f"{uid}_raw_skeleton.json")))
        p0 = np.asarray(sk["trunk"]["polylines"][0])[0]
        check(f"{uid} raw polyline0 anchored at base_xy",
              np.allclose(p0[:2], sk["base_xy"], atol=1e-12),
              f"{p0[:2]} vs {sk['base_xy']}")
        far = max(float(np.hypot(*(np.asarray(pl, dtype=float)[:, :2]
                                   - np.asarray(sk["base_xy"])).T).max())
                  for pl in sk["trunk"]["polylines"])
        check(f"{uid} trunk within 25 m of base", far < 25.0, f"{far:.1f} m")

    print("== 1. determinism: same seed twice ==")
    scene1 = ap.assemble_plot(table, seed=fixture.SEED, verbose=False)
    trim = os.path.join(plot_dir, "plot_placement",
                        f"{fixture.PLOT_UID}_trimmed.laz")
    shutil.copy(trim, os.path.join(tmp, "run1.laz"))
    scene2 = ap.assemble_plot(table, seed=fixture.SEED, verbose=False)
    shutil.copy(trim, os.path.join(tmp, "run2.laz"))
    xyz1, _ = read_canonical_las(os.path.join(tmp, "run1.laz"), level="plot")
    xyz2, _ = read_canonical_las(os.path.join(tmp, "run2.laz"), level="plot")
    check("trimmed cloud bit-identical across same-seed runs",
          xyz1.shape == xyz2.shape and np.array_equal(xyz1, xyz2))
    check("scene reflections identical across runs",
          scene1["reflections"] == scene2["reflections"])
    scene = ap.assemble_plot(table, seed=fixture.SEED, verbose=False)
    refl, rots = scene["reflections"], scene["rotations_deg"]
    uids = sorted(refl)

    print("== 2. scene records ==")
    check("reflections covers every placed tree", len(uids) == fixture.N_TREES)
    check("reflections values are bool",
          all(isinstance(v, bool) for v in refl.values()))
    check("both reflection states occur",
          True in refl.values() and False in refl.values(),
          str(refl))

    print("== 3. cloud-vs-records parity (tol = 1mm re-quantization) ==")
    rawm = os.path.join(plot_dir, "plot_placement",
                        f"{fixture.PLOT_UID}_raw.laz")
    raw_xyz, raw_cols = read_canonical_las(rawm, level="plot")
    for uid in uids:
        txyz, _ = read_canonical_las(os.path.join(treedir, f"{uid}_raw.laz"),
                                     level="tree")
        sk = json.load(open(os.path.join(treedir, f"{uid}_raw_skeleton.json")))
        cx, cy = float(sk["base_xy"][0]), float(sk["base_xy"][1])
        tid = int(sk["tree_id"])
        sel = raw_cols["tree_id"] == tid
        exp = manual_transform(txyz[:, :2], math.radians(rots[uid]),
                               refl[uid], cx, cy)
        dmax = float(np.abs(exp - raw_xyz[sel][:, :2]).max())
        check(f"{uid} ({int(sel.sum())} pts) manual transform matches "
              f"(max |d| {dmax:.2e})",
              exp.shape == raw_xyz[sel][:, :2].shape and dmax < 6e-4)

    print("== 4. skeleton-vs-records parity (tol = 4dp angle record) ==")
    for uid in uids:
        raw_sk = json.load(open(os.path.join(treedir,
                                             f"{uid}_raw_skeleton.json")))
        trm_sk = json.load(open(os.path.join(
            treedir, f"{uid}_trimmed_skeleton.json")))
        check(f"{uid} skeleton records reflected flag",
              isinstance(trm_sk.get("reflected"), bool))
        cx, cy = float(raw_sk["base_xy"][0]), float(raw_sk["base_xy"][1])
        ang = math.radians(trm_sk["rotation_deg"])
        ok = True
        for pi, pl in enumerate(raw_sk["trunk"]["polylines"]):
            exp = manual_transform(np.asarray(pl, dtype=float)[:, :2], ang,
                                   trm_sk["reflected"], cx, cy)
            got = np.asarray(trm_sk["trunk"]["polylines"][pi])[:, :2]
            ok = ok and np.allclose(exp, got, atol=2e-5, rtol=0)
        for br in raw_sk["branches"]:
            for key in ("base", "tip"):
                v = np.asarray(br[key], dtype=float)[:2]
                exp = manual_transform(v[None, :], ang, trm_sk["reflected"],
                                       cx, cy)[0]
                cands = [b for b in trm_sk["branches"]
                         if abs(b[key][2] - br[key][2]) < 1e-12]
                got_b = min(cands, key=lambda b: abs(b[key][0] - exp[0])
                            + abs(b[key][1] - exp[1]))
                ok = ok and np.allclose(
                    exp, np.asarray(got_b[key], dtype=float)[:2],
                    atol=2e-5, rtol=0)
        check(f"{uid} serialized skeleton matches manual transform", ok)

    print("== 5. reflection draws vary across seeds (~50%) ==")
    total = 0
    for s in range(8):
        sc = ap.assemble_plot(table, seed=1000 + s, verbose=False)
        total += sum(sc["reflections"].values())
    rate = total / (fixture.N_TREES * 8)
    check("both states occur across seeds", 0 < total < fixture.N_TREES * 8)
    check(f"rate roughly balanced (observed {rate:.2f})",
          0.25 <= rate <= 0.75)
    ap.assemble_plot(table, seed=fixture.SEED, verbose=False)  # restore

    print("== 6. chirality: reflection flips orientation exactly ==")
    for uid in uids:
        raw_sk = json.load(open(os.path.join(treedir,
                                             f"{uid}_raw_skeleton.json")))
        trm_sk = json.load(open(os.path.join(
            treedir, f"{uid}_trimmed_skeleton.json")))
        P_raw = np.array([b["base"][:2] for b in raw_sk["branches"]])
        P_out = np.array([b["base"][:2] for b in trm_sk["branches"]])
        s_raw, s_out = signed_area(P_raw), signed_area(P_out)
        keep = abs(s_out - s_raw) <= 1e-6 * max(abs(s_raw), 1.0)
        flip = abs(s_out + s_raw) <= 1e-6 * max(abs(s_raw), 1.0)
        ok = flip if trm_sk["reflected"] else keep
        check(f"{uid} signed area {'flips' if trm_sk['reflected'] else 'preserved'}"
              f" ({s_raw:.4f} -> {s_out:.4f})", ok)


if __name__ == "__main__":
    aparser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    aparser.add_argument("--keep", action="store_true",
                         help="keep the fixture directory for debugging")
    args = aparser.parse_args()
    sys.exit(main(keep=args.keep))
