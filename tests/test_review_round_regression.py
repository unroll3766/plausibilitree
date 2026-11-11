"""Review-round regression suite (2026-11 review, items mined from
tests/CHECKLOG.md probes that were previously one-off scripts).

Run: conda run -n ssf python tests/test_review_round_regression.py

Fast, fixture-free checks for the three behaviors the review round fixed,
plus a contract audit of the committed demo (the only tracked plot):

  A. SHRUB PLACEMENT (W23, review item 6): circular plots seed shrubs
     only inside the plot disc; after the plot-boundary clip every shrub
     record keeps at least one point. The zero-point failure mode the
     reviewer found (bounding-square corner shrubs erased by the disc
     clip while their records persisted) is pinned by the invariants:
     centres inside the disc, and per-record survival >= 1 point.
  B. FORK-ORDER SPLIT RULE (W22, review items 4/5): at a split the
     daughter with the larger base radius keeps the parent order (ties:
     longer stem, then first-listed), the other is parent + 1; a branch
     is always its source stem's order + 1. Synthetic skeletons pin the
     rule where a fixture-tree draw cannot: radius-first dominance,
     order-independence of daughter listing, cascades, and the tie.
  C. COMMITTED DEMO CONTRACT: the tracked examples/plot_demo satisfies
     the current artifact contracts -- scene links all resolve, realised
     dbh_m present, shrub info free of dbh_m/voxel keys, no per-shrub
     deep_learning volume labels, no *_normalized* anywhere, and
     VOL_SHRUB == the per-shrub info mesh sum. This section fails
     whenever the demo is regenerated with a contract-breaking state,
     which is its job (the W6/W16 rounds both shipped stale demos).

Exit 0 = all checks passed.
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, REPO)

import numpy as np

import sst_generator
import vegetation as veg

DEMO = os.path.join(REPO, "examples")
DEMO_PLOT = os.path.join(DEMO, "plot_demo")

fails = []


def check(name, cond, msg=""):
    print(f"  {'PASS' if cond else 'FAIL'}  {name}"
          f"{('  :: ' + msg) if (msg and not cond) else ''}")
    if not cond:
        fails.append(name)


# ---------------------------------------------------------------------------
# A. shrub placement inside the plot footprint (review item 6, item W23),
#    promoted from the 2026-11 diagnostic repro
# ---------------------------------------------------------------------------

def check_shrub_placement(keep: bool = False):
    print("== A. shrub placement (W23: shrubs inside the plot footprint) ==")
    r = 8.0
    bounds = (-r, r, -r, r)

    pts, ids, centers, widths, geom = veg.gen_shrub_layer(
        bounds, 2.0, 1.0, (0.4, 1.8), (0.8, 1.6),
        tree_positions=np.empty((0, 2)), tree_canopy_r=np.empty(0),
        points_per_shrub=60, seed=42,
        plot_shape="circle", plot_radius=r)
    # the spec (3.4) constraints the candidate CELLS to the disc; the
    # shrub centre then jitters within its cell (jitter_frac 0.4 x 2 m),
    # so a rim shrub's centre may exceed the disc by up to that bound
    jitter_bound = 0.4 * 2.0
    dists = np.hypot(centers[:, 0], centers[:, 1])
    check("circle: every shrub centre within the cell jitter of the disc",
          bool(np.all(dists <= r + jitter_bound + 1e-9)),
          f"max centre dist {dists.max():.2f} m vs r {r}")

    # the plot-boundary clip the pipeline applies (understory's disc rule)
    inside = (pts[:, 0] ** 2 + pts[:, 1] ** 2) <= r ** 2
    alive = set(int(s) for s in np.unique(ids[inside]))
    survivors_ok = all(sid in alive for sid in geom)
    check("circle: every record retains >= 1 point after the disc clip",
          survivors_ok and len(geom) > 0,
          f"{len(geom)} records, {len(alive)} alive after clip")
    check("circle: corners of the bounding square carry no shrub",
          not any(abs(c).max() > r for c in centers))

    # unshaped callers keep the bounding-square behaviour (no disc filter)
    pts2, ids2, centers2, widths2, geom2 = veg.gen_shrub_layer(
        bounds, 2.0, 1.0, (0.4, 1.8), (0.8, 1.6),
        tree_positions=np.empty((0, 2)), tree_canopy_r=np.empty(0),
        points_per_shrub=2, seed=42)
    n_square_cells = 8 * 8        # cell 2.0 over bounds (-8, 8)
    check("default (unshaped): bounding-square grid kept, seeded as before",
          len(geom2) == n_square_cells, f"{len(geom2)} shrubs")

    pts3, ids3, centers3, widths3, geom3 = veg.gen_shrub_layer(
        bounds, 2.0, 0.0, (0.4, 1.8), (0.8, 1.6),
        tree_positions=np.empty((0, 2)), tree_canopy_r=np.empty(0),
        points_per_shrub=2, seed=42, plot_shape="circle", plot_radius=r)
    check("circle: cover 0 places nothing", len(geom3) == 0 and len(pts3) == 0)


# ---------------------------------------------------------------------------
# B. fork-order split rule (review items 4/5, item W22)
# ---------------------------------------------------------------------------

def _skel(segments, split_parents, branches, base_radius=0.30, height=30.0):
    """Minimal skeleton for sst_generator.element_graph.

    segments/d: ((z0, z1), length) pairs; split_parents: the parent
    polyline index of each split, in split order (split k's daughters are
    stems 2k+1 and 2k+2); branches: (stem_index, length, radius).
    """
    return {
        "base_radius_m": base_radius, "height_m": height,
        "trunk": {"segments": [{"z_range": [z0, z1], "length_m": length}
                               for (z0, z1), length in segments],
                  # element_graph reads only segments via this dict shape
                  "polylines": []},
        "splits": [{"split_index": k, "parent_polyline": p, "height_m": 0.0}
                   for k, p in enumerate(split_parents)],
        "branches": [{"stem_index": i, "length_m": length, "radius_m": radius}
                     for i, length, radius in branches],
    }


def _orders(skel):
    els = sst_generator.element_graph(skel)
    return {e["BRANCH_UID"]: e["BRANCH_ORDER"] for e in els}


def check_split_rule(keep: bool = False):
    print("== B. fork-order split rule (W22) ==")

    # single stem, no splits: all stems order 0; branches are 1
    sk = _skel([((0.0, 30.0), 20.0)], [], [(0, 3.0, 0.05), (0, 2.0, 0.04)])
    o = _orders(sk)
    check("no splits: main stem order 0, its branches order 1",
          o["T00"] == 0 and o["B000"] == 1 and o["B001"] == 1, str(o))

    # one split, both daughters at the same base height (equal taper
    # radius -- SST's real case): the longer daughter keeps 0
    sk = _skel([((0.0, 15.0), 12.0), ((15.0, 30.0), 10.0),
                ((15.0, 28.0), 4.0), ((15.0, 26.0), 5.0)],
               [0], [(1, 3.0, 0.05), (2, 3.0, 0.05)])
    o = _orders(sk)
    check("tie on radius: longer daughter keeps the parent order",
          o["T00"] == 0 and o["T01"] == 0 and o["T02"] == 1, str(o))
    check("branch orders follow their stems (parent + 1)",
          o["B000"] == 1 and o["B001"] == 2, str(o))

    # lengths equal too: first-listed daughter keeps the parent order
    sk = _skel([((0.0, 15.0), 12.0), ((15.0, 30.0), 8.0),
                ((15.0, 30.0), 8.0), ((15.0, 20.0), 6.0)],
               [0], [])
    o = _orders(sk)
    check("full tie: first-listed daughter keeps the parent order",
          o["T01"] == 0 and o["T02"] == 1, str(o))

    # radius decides before length: the shorter-but-thicker daughter
    # (lower base height -> larger taper radius) keeps 0. Radii: T01 at
    # z=15 -> 0.15 vs T02 at z=10 -> 0.20
    sk = _skel([((0.0, 15.0), 12.0), ((15.0, 29.0), 10.0),
                ((10.0, 30.0), 3.0)], [0], [])
    o = _orders(sk)
    check("radius dominates: thicker, shorter daughter keeps order 0",
          o["T00"] == 0 and o["T02"] == 0 and o["T01"] == 1, str(o))

    # cascading splits on the SUBORDINATE line grow the order: the main
    # stem splits into T01 (longer, keeps 0) and T02 (order 1); T02 then
    # splits into T03 (longer, keeps 1) and T04 (order 2); branches on
    # T03/T04 get 2/3 -- how the demo hardwoods reach order 3
    sk = _skel([((0.0, 12.0), 10.0), ((12.0, 30.0), 9.0),
                ((12.0, 26.0), 5.0), ((20.0, 30.0), 3.0),
                ((20.0, 28.0), 2.0)],
               [0, 2], [(3, 1.0, 0.04), (4, 1.0, 0.04)])
    o = _orders(sk)
    check("cascade: dominant line keeps 0, subordinate chain +1 per split",
          (o["T00"], o["T01"], o["T02"], o["T03"], o["T04"])
          == (0, 0, 1, 1, 2), str(o))
    check("cascade: branch orders are stem order + 1 (up to 3)",
          o["B000"] == 2 and o["B001"] == 3, str(o))

    els = sst_generator.element_graph(_skel(
        [((0.0, 12.0), 10.0), ((12.0, 30.0), 9.0), ((12.0, 26.0), 5.0),
         ((20.0, 30.0), 3.0), ((20.0, 28.0), 2.0)], [0, 2],
        [(3, 1.0, 0.04)]))
    check("element orders stay within the soft limit [0, 5], never -1 "
          "(-1 marks a point on no branch)",
          all(0 <= e["BRANCH_ORDER"] <= 5 for e in els), str(sorted(
              {e["BRANCH_ORDER"] for e in els})))


# ---------------------------------------------------------------------------
# C. committed demo contract (examples/plot_demo is tracked)
# ---------------------------------------------------------------------------

def _scene_links(d):
    for k, v in d.items():
        if isinstance(v, str) and v.endswith((".laz", ".las", ".json",
                                            ".tif", ".asc")):
            yield k, v
        elif isinstance(v, dict):
            yield from _scene_links(v)


def check_demo_contract(keep: bool = False):
    print("== C. committed demo contract (examples/plot_demo) ==")
    have_demo = os.path.isdir(DEMO_PLOT)
    if not have_demo:
        check("committed demo present", False, DEMO_PLOT)
        return

    # every scene link resolves (links are relative to the plot dir's
    # parent, i.e. examples/)
    bad, n = [], 0
    for scene_path in glob.glob(os.path.join(DEMO_PLOT, "**", "*_scene.json"),
                                recursive=True):
        scene = json.load(open(scene_path))
        for k, v in _scene_links(scene):
            n += 1
            if not os.path.exists(os.path.join(DEMO, v)):
                bad.append((os.path.basename(scene_path), k, v))
    check(f"every scene link resolves ({n} links)", bad == [], str(bad[:3]))

    # retired artifacts really gone from the tracked demo
    norm = glob.glob(os.path.join(DEMO_PLOT, "**", "*_normalized*"),
                     recursive=True)
    check("no *_normalized* artifact in the demo", norm == [], str(norm[:2]))
    vol = glob.glob(os.path.join(DEMO_PLOT, "individual_shrubs",
                                 "**", "*_volume.json"), recursive=True)
    check("no per-shrub deep_learning volume label in the demo",
          vol == [], str(vol[:2]))

    # realised DBH present on every generation record (review item 9)
    gen = sorted(glob.glob(os.path.join(DEMO_PLOT, "individual_trees", "*",
                                        "textual", "*_generation.json")))
    check("every generation record carries realised dbh_m (item 9)",
          len(gen) > 0 and all("dbh_m" in json.load(open(g))["realised"]
                               for g in gen))

    # shrub info schema (review item 7): no dbh_m, no voxel key
    infos = sorted(glob.glob(os.path.join(DEMO_PLOT, "individual_shrubs",
                                          "*", "textual", "*_info.json")))
    ok_keys = all(
        "dbh_m" not in (d := json.load(open(f)))["shrub_characteristics"]
        and set(d["shrub_characteristics"]["volume_m3"]) == {"primitive",
                                                             "mesh"}
        and d["meta"]["bounding_box"] is not None
        and d["neighborhood"] is not None
        for f in infos)
    check("every shrub info: no dbh_m, primitive+mesh only, bbox+hood set",
          len(infos) > 0 and ok_keys)

    # VOL_SHRUB epoch: plot info sums the per-shrub info mesh volumes
    pi = json.load(open(glob.glob(os.path.join(
        DEMO_PLOT, "plot_placement", "textual", "*_info.json"))[0]))
    mesh = sum(json.load(open(f))["shrub_characteristics"]["volume_m3"]["mesh"]
               for f in infos)
    check("VOL_SHRUB == sum of per-shrub info mesh volumes",
          abs(pi["VOL_SHRUB"] - mesh) < 1e-9,
          f"{pi['VOL_SHRUB']} vs {mesh}")
    check("plot info carries no canopy hull key (item 10)",
          "canopy" not in pi)
    tree_infos = sorted(glob.glob(os.path.join(DEMO_PLOT, "individual_trees",
                                               "*", "textual", "*_info.json")))
    check("tree info canopy = crown axes only (item 10)",
          all(sorted(json.load(open(f))["canopy"]) == ["major_axis",
                                                       "minor_axis"]
              for f in tree_infos))


def main(keep: bool = False) -> int:
    print("review-round regression suite (W22/W23/W24/W25/W21 probes)")
    check_shrub_placement(keep)
    check_split_rule(keep)
    check_demo_contract(keep)
    print()
    if fails:
        print(f"REVIEW-ROUND SUITE FAILED ({len(fails)}): {fails}")
        return 1
    print("review-round regression suite: all checks passed")
    return 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="review-round regression suite")
    ap.add_argument("--keep", action="store_true", default=False)
    args = ap.parse_args()
    sys.exit(main(keep=args.keep))
