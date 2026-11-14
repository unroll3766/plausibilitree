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
  D. STEM-MAP MASK RULES (W28): only branch-order-0 points enter the stem
     map mask, and a tree that splits below the DBH range (fork height <
     DBH_HT - STEM_SECTION) is excluded from the mask entirely. Unit
     probes pin fork_below_band's boundary semantics, polyline_orders'
     segment/element fallback, and the order-0 disc filter; a synthetic
     3-tree plot run through stem_map.build_stem_map pins the end-to-end
     behaviour (band 2 cells, bounding-box JSON, per-tree masks).
  E. CHM MASK = FILLED CANOPY FOOTPRINT (W29): the canopy mask is generated
     from geometry like the stem map's — the filled 2D convex hull of each
     tree's foliage points, interior included — not from point occupancy.
     Unit probes pin rasterize_polygon (interior fill, degenerate fallback)
     and canopy_footprint_hull; a synthetic 3-tree plot run through
     canopy_map.build_canopy_map (+ restructure's instance mask) pins band 2,
     per-tree masks, hull JSON footprints and taller-crown-wins overlap.

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


# ---------------------------------------------------------------------------
# D. stem-map mask rules (W28): order-0 only + below-band-fork exclusion
# ---------------------------------------------------------------------------

def _w28_skel(tree_id, base_xy, polylines, orders, splits, h=10.0, r=0.15):
    """Minimal trimmed-style skeleton: polylines as (N,3) lists, per-segment
    branch_order, and a split log with height_m entries."""
    segs = []
    for pl, o in zip(polylines, orders):
        arr = np.asarray(pl, dtype=float)
        segs.append({"z_range": [float(arr[:, 2].min()), float(arr[:, 2].max())],
                     "length_m": float(np.linalg.norm(np.diff(arr, axis=0),
                                                      axis=1).sum()),
                     "branch_order": int(o)})
    return {
        "tree_id": int(tree_id),
        "base_xy": [float(base_xy[0]), float(base_xy[1])],
        "height_m": float(h), "base_radius_m": float(r),
        "fork": 2 if splits else 1, "n_splits": len(splits),
        "splits": splits,
        "trunk": {"polylines": [[list(p) for p in pl] for pl in polylines],
                  "segments": segs},
        "branches": [], "n_branches": 0,
        "elements": [],
    }


def check_stem_mask_rules(keep: bool = False):
    print("== D. stem-map mask rules (W28) ==")
    import shutil
    import tempfile
    import rasterio
    import stem_map as sm
    from las import export_canonical_las

    z_lo = sm.BREAST_HEIGHT_M - sm.STEM_SECTION      # 1.12
    z_hi = sm.BREAST_HEIGHT_M + sm.STEM_SECTION      # 1.62

    # --- unit: fork_below_band (a tree is excluded iff a fork lies strictly
    # below the bottom of the DBH band) ---
    def _sp(h):
        return [{"split_index": 0, "parent_polyline": 0, "height_m": float(h)}]
    check("fork_below_band: split at 0.5 m (< 1.12) is below the band",
          sm.fork_below_band(_w28_skel(1, (0, 0), [[[0, 0, 0], [0, 0, 10]]],
                                      [0], _sp(0.5)), z_lo) == 0.5)
    check("fork_below_band: split exactly at z_lo is NOT below (boundary)",
          sm.fork_below_band(_w28_skel(1, (0, 0), [[[0, 0, 0], [0, 0, 10]]],
                                      [0], _sp(z_lo)), z_lo) is None)
    check("fork_below_band: split within the band (1.3 m) is NOT below",
          sm.fork_below_band(_w28_skel(1, (0, 0), [[[0, 0, 0], [0, 0, 10]]],
                                      [0], _sp(1.3)), z_lo) is None)
    check("fork_below_band: no splits -> None",
          sm.fork_below_band(_w28_skel(1, (0, 0), [[[0, 0, 0], [0, 0, 10]]],
                                      [0], []), z_lo) is None)
    two = [{"split_index": 0, "parent_polyline": 0, "height_m": 5.0},
           {"split_index": 1, "parent_polyline": 2, "height_m": 0.8}]
    check("fork_below_band: lowest below-band split wins",
          sm.fork_below_band(_w28_skel(1, (0, 0), [[[0, 0, 0], [0, 0, 10]]],
                                      [0], two), z_lo) == 0.8)

    # --- unit: polyline_orders (segments first, element graph fallback) ---
    sk = _w28_skel(1, (0, 0), [[[0, 0, 0], [0, 0, 5]],
                               [[0, 0, 5], [0.1, 0, 10]],
                               [[0, 0, 5], [-0.1, 0, 8]]], [0, 0, 1],
                   _sp(5.0))
    check("polyline_orders reads per-segment branch_order",
          sm.polyline_orders(sk) == [0, 0, 1])
    sk_no_seg_order = json.loads(json.dumps(sk))
    for s in sk_no_seg_order["trunk"]["segments"]:
        s.pop("branch_order")
    sk_no_seg_order["elements"] = [
        {"BRANCH_UID": "T00", "BRANCH_ORDER": 0},
        {"BRANCH_UID": "T01", "BRANCH_ORDER": 0},
        {"BRANCH_UID": "T02", "BRANCH_ORDER": 1},
        {"BRANCH_UID": "B000", "BRANCH_ORDER": 1}]
    check("polyline_orders falls back to the T{i} element graph",
          sm.polyline_orders(sk_no_seg_order) == [0, 0, 1])
    sk_bare = json.loads(json.dumps(sk_no_seg_order))
    sk_bare["elements"] = []
    check("polyline_orders defaults to 0 with no order info anywhere",
          sm.polyline_orders(sk_bare) == [0, 0, 0])

    # --- unit: stem_discs_in_band honours the order-0 filter ---
    sk2 = _w28_skel(2, (0, 0), [[[0, 0, 0], [0, 0, 10]],
                                [[0, 0, 0.5], [1.0, 0, 10]]], [0, 1], _sp(0.5))
    n_all = len(sm.stem_discs_in_band(sk2, z_lo, z_hi))
    n_o0 = len(sm.stem_discs_in_band(sk2, z_lo, z_hi,
                                     orders=sm.polyline_orders(sk2)))
    check("stem_discs_in_band: order-0 filter drops the order-1 polyline's discs",
          n_all == 2 * n_o0 and n_o0 > 0, f"all={n_all} order0={n_o0}")

    # --- end-to-end on a synthetic 3-tree flat plot ---
    # tree 1: single stem (in mask); tree 2: fork at 0.5 m (below the band ->
    # excluded); tree 3: fork at 1.3 m (within the band -> kept, but only its
    # order-0 polylines may mask)
    tmp = tempfile.mkdtemp(prefix="w28_stem_")
    try:
        plot_uid = "plot_w28"
        csv_path = os.path.join(tmp, f"{plot_uid}.csv")
        with open(csv_path, "w") as f:
            f.write("TREE_UID,X,Y,DBH,HT,TREE_CODE\n"
                    "001,0.0,0.0,0.3,10.0,HW\n"
                    "002,5.0,0.0,0.3,10.0,HW\n"
                    "003,0.0,5.0,0.3,10.0,HW\n")
        plot_dir = os.path.join(tmp, plot_uid)
        place = os.path.join(plot_dir, "plot_placement")
        trees_dir = os.path.join(plot_dir, "individual_trees")
        os.makedirs(place, exist_ok=True)

        pts = np.array([
            [0.00, 0.0, 1.20],   # tree 1, order 0 -> masks
            [0.12, 0.0, 1.30],   # tree 1, branch_code 2 -> must NOT mask
            [5.00, 0.0, 1.30],   # tree 2 (excluded) -> must NOT mask
            [0.00, 5.0, 1.20],   # tree 3, order 0 -> masks
            [9.00, 9.0, 1.30],   # ground -> never masks
        ])
        export_canonical_las(
            os.path.join(place, f"{plot_uid}_raw.laz"), pts,
            {"tree_id": np.array([1, 1, 2, 3, 0], dtype=np.int32),
             "semantic_code": np.array([1, 2, 1, 1, 0], dtype=np.int8),
             "branch_code": np.array([1, 2, 1, 1, 0], dtype=np.int8)},
            level="plot")

        skels = {
            1: _w28_skel(1, (0, 0), [[[0, 0, 0], [0, 0, 10]]], [0], []),
            2: _w28_skel(2, (5, 0), [[[5, 0, 0], [5, 0, 0.5]],
                                     [[5, 0, 0.5], [5.1, 0, 10]],
                                     [[5, 0, 0.5], [4.9, 0, 8]]], [0, 0, 1],
                         _sp(0.5)),
            3: _w28_skel(3, (0, 5), [[[0, 5, 0], [0, 5, 1.3]],
                                     [[0, 5, 1.3], [0.1, 5, 10]],
                                     [[0, 5, 1.3], [-2.0, 5, 4.0]]], [0, 0, 1],
                         _sp(1.3)),
        }
        for tid, skel in skels.items():
            uid = f"{tid:03d}"
            d = os.path.join(trees_dir, f"{plot_uid}_{uid}", "textual")
            os.makedirs(d, exist_ok=True)
            with open(os.path.join(d, f"{plot_uid}_{uid}_trimmed_skeleton.json"),
                      "w") as f:
                json.dump(skel, f)

        meta = sm.build_stem_map(csv_path, plot_dir=plot_dir, verbose=False)

        bbox_doc = json.load(open(os.path.join(
            place, f"{plot_uid}_raw_placement_stem_map_bounding_boxes.json")))
        check("e2e: below-band-fork tree is excluded, others get boxes",
              [b["tree_id"] for b in bbox_doc["boxes"]] == [1, 3]
              and [e["tree_id"] for e in bbox_doc["excluded_trees"]] == [2]
              and bbox_doc["excluded_trees"][0]["lowest_fork_below_band_m"] == 0.5,
              str(bbox_doc.get("excluded_trees")))
        check("e2e: scene meta records the exclusion",
              meta["n_trees_excluded"] == 1 and meta["n_trees_with_stem"] == 2)

        # band 2 == cells holding a branch_code==1 point of a non-excluded tree
        with rasterio.open(os.path.join(
                place, f"{plot_uid}_raw_placement_stem_map.tif")) as src:
            b1, b2 = src.read(1), src.read(2)
            t = src.transform
        cell = float(t.a); x0 = float(t.c); y1 = float(t.f)
        grid = sm.RasterGrid(x0, y1 - b1.shape[0] * cell,
                             x0 + b1.shape[1] * cell, y1, cell)
        tid = np.array([1, 1, 2, 3, 0])
        code = np.array([1, 2, 1, 1, 0])
        sel = ((pts[:, 2] >= z_lo) & (pts[:, 2] <= z_hi)
               & (tid != 0) & (code == 1) & (tid != 2))
        exp = (grid.counts(pts[sel, 0], pts[sel, 1]) > 0).astype(np.uint16)
        check("e2e: band 2 masks exactly the order-0 points of non-excluded trees",
              np.array_equal(exp, b2) and int((b2 > 0).sum()) == 2,
              f"cells={int((b2 > 0).sum())} expected={int((exp > 0).sum())}")
        check("e2e: band 2 is still binary and implies band-1 support",
              set(np.unique(b2).tolist()) <= {0, 1}
              and not np.any((b2 > 0) & (b1 <= 0)))

        # per-tree crops/masks exist only for non-excluded trees
        for tid, want in ((1, True), (2, False), (3, True)):
            uid = f"{tid:03d}"
            p = os.path.join(trees_dir,
                             f"{plot_uid}_tree_{uid}_stem_map_semantic_mask.tif")
            check(f"e2e: per-tree stem mask written iff not excluded (tree {tid})",
                  os.path.exists(p) == want)

        # tree 3's mask == its order-0 discs only (the order-1 fork daughter
        # veers to x=-2 m and must not contribute)
        box3 = next(b for b in bbox_doc["boxes"] if b["tree_id"] == 3)
        bp = box3["bbox_pixel"]
        with rasterio.open(os.path.join(
                trees_dir,
                f"{plot_uid}_tree_003_stem_map_semantic_mask.tif")) as src:
            mask3 = src.read(1)
        discs_o0 = sm.stem_discs_in_band(skels[3], z_lo, z_hi,
                                        orders=sm.polyline_orders(skels[3]))
        discs_all = sm.stem_discs_in_band(skels[3], z_lo, z_hi)
        full_o0 = grid.rasterize_discs(discs_o0)
        full_all = grid.rasterize_discs(discs_all)
        crop_o0 = full_o0[bp["row_min"]:bp["row_max"] + 1,
                         bp["col_min"]:bp["col_max"] + 1]
        crop_all = full_all[bp["row_min"]:bp["row_max"] + 1,
                            bp["col_min"]:bp["col_max"] + 1]
        check("e2e: per-tree mask == order-0 skeleton discs",
              np.array_equal(mask3, crop_o0.astype(mask3.dtype)),
              f"mask cells={int((mask3 > 0).sum())} expected={int((crop_o0 > 0).sum())}")
        check("e2e: order-1 fork daughter would add cells (filter is load-bearing)",
              int((crop_all > 0).sum()) > int((crop_o0 > 0).sum())
              and not np.array_equal(mask3, crop_all.astype(mask3.dtype)))
        check("e2e: box n_discs counts order-0 discs only",
              box3["n_discs"] == len(discs_o0) < len(discs_all))
    finally:
        if not keep:
            shutil.rmtree(tmp, ignore_errors=True)


def check_chm_mask_rules(keep: bool = False):
    print("== E. CHM mask = filled canopy footprint (W29) ==")
    import shutil
    import tempfile
    import rasterio
    import canopy_map as cm
    import stem_map as sm
    from las import export_canonical_las

    # --- unit: RasterGrid.rasterize_polygon fills the interior ---
    g = sm.RasterGrid(0.0, 0.0, 1.0, 1.0, 0.05)     # 20x20 cells over [0,1]^2
    square = np.array([[0.1, 0.1], [0.9, 0.1], [0.9, 0.9], [0.1, 0.9]])
    m = g.rasterize_polygon(square)
    ccol, crow = g.to_pixel([0.5], [0.5])
    ocol, orow = g.to_pixel([0.02], [0.02])
    ccol, crow = int(ccol[0]), int(crow[0])
    ocol, orow = int(ocol[0]), int(orow[0])
    check("rasterize_polygon fills the interior cell (centre of the square)",
          m[crow, ccol] == 1)
    check("rasterize_polygon leaves outside cells 0", m[orow, ocol] == 0)
    deg = g.rasterize_polygon(np.array([[0.2, 0.2], [0.8, 0.8]]))
    d1 = g.to_pixel([0.2], [0.2]); d2 = g.to_pixel([0.8], [0.8])
    d1 = (int(d1[0][0]), int(d1[1][0])); d2 = (int(d2[0][0]), int(d2[1][0]))
    check("rasterize_polygon degenerate (<3 verts): vertex cells only",
          int(deg.sum()) == 2 and deg[d1[1], d1[0]] == 1 and deg[d2[1], d2[0]] == 1)

    # --- unit: canopy_footprint_hull ---
    check("canopy_footprint_hull: square -> 4 vertices",
          len(cm.canopy_footprint_hull(square)) == 4)
    check("canopy_footprint_hull: collinear -> None (fallback)",
          cm.canopy_footprint_hull(np.array([[0.0, 0.0], [1.0, 1.0],
                                             [2.0, 2.0]])) is None)
    check("canopy_footprint_hull: <3 points -> None",
          cm.canopy_footprint_hull(np.array([[0.0, 0.0], [1.0, 0.0]])) is None)

    # --- end-to-end on a synthetic 3-tree flat plot ---
    # tree 1: foliage on a 1 m square's corners at z=8 (interior empty);
    # tree 2: three collinear foliage points (degenerate hull -> point cells);
    # tree 3: foliage on a 0.8 m square at z=12, overlapping tree 1 (taller).
    tmp = tempfile.mkdtemp(prefix="w29_chm_")
    try:
        plot_uid = "plot_w29"
        csv_path = os.path.join(tmp, f"{plot_uid}.csv")
        with open(csv_path, "w") as f:
            f.write("TREE_UID,X,Y,DBH,HT,TREE_CODE\n"
                    "001,0.0,0.0,0.3,10.0,HW\n"
                    "002,6.0,0.0,0.2,6.0,HW\n"
                    "003,0.0,0.0,0.3,12.0,HW\n")
        plot_dir = os.path.join(tmp, plot_uid)
        place = os.path.join(plot_dir, "plot_placement")
        os.makedirs(place, exist_ok=True)

        pts = np.array([
            # tree 1 corners (z=8) + one trunk point away from the centre cell
            [-0.5, -0.5, 8.0], [0.5, -0.5, 8.0], [0.5, 0.5, 8.0], [-0.5, 0.5, 8.0],
            [0.3, 0.3, 1.0],
            # tree 2: collinear foliage
            [5.0, 0.0, 3.0], [6.0, 0.0, 3.0], [7.0, 0.0, 3.0],
            # tree 3 corners (z=12), inside tree 1's square
            [-0.4, -0.4, 12.0], [0.4, -0.4, 12.0], [0.4, 0.4, 12.0], [-0.4, 0.4, 12.0],
        ])
        sem = np.array([3, 3, 3, 3, 1, 3, 3, 3, 3, 3, 3, 3], dtype=np.int8)
        tid = np.array([1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3], dtype=np.int32)
        code = np.where(sem == 3, 2, 1).astype(np.int8)
        export_canonical_las(os.path.join(place, f"{plot_uid}_raw.laz"), pts,
                             {"tree_id": tid, "semantic_code": sem,
                              "branch_code": code}, level="plot")

        cm.build_canopy_map(csv_path, plot_dir=plot_dir, verbose=False)

        feat_path = os.path.join(place, f"{plot_uid}_raw_placement_canopy_map.tif")
        with rasterio.open(feat_path) as src:
            b1, b2 = src.read(1), src.read(2)
            tr = src.transform
        cell = float(tr.a); x0 = float(tr.c); y1 = float(tr.f)
        grid = sm.RasterGrid(x0, y1 - b1.shape[0] * cell,
                             x0 + b1.shape[1] * cell, y1, cell)

        # interior cell of tree 1's square with NO point in it: band 2 fills it
        icol, irow = grid.to_pixel([0.0], [0.0])
        icol, irow = int(icol[0]), int(irow[0])
        check("e2e: band 2 fills a canopy-interior cell that holds no point",
              int(b2[irow, icol]) == 1 and int(b1[irow, icol]) == 0,
              f"b2={b2[irow, icol]} b1={b1[irow, icol]}")
        for wx, wy in ((-0.5, -0.5), (0.5, 0.5), (-0.4, 0.4)):
            cc, rr = grid.to_pixel([wx], [wy])
            cc, rr = int(cc[0]), int(rr[0])
            check(f"e2e: band 2 covers foliage corner ({wx},{wy})",
                  int(b2[rr, cc]) == 1)
        # tree 2's midpoint is NOT filled (degenerate hull -> point cells only)
        mcol, mrow = grid.to_pixel([5.5], [0.0])
        mcol, mrow = int(mcol[0]), int(mrow[0])
        check("e2e: degenerate-hull tree fills point cells only (midpoint empty)",
              int(b2[mrow, mcol]) == 0)
        for wx in (5.0, 6.0, 7.0):
            cc, rr = grid.to_pixel([wx], [0.0])
            cc, rr = int(cc[0]), int(rr[0])
            check(f"e2e: tree 2 point cell ({wx},0) set", int(b2[rr, cc]) == 1)

        # per-tree masks: tree 1 filled (interior set), tree 2 point cells only
        trees_dir = os.path.join(plot_dir, "individual_trees")
        for uid, want_interior, npts in (("001", True, None), ("002", False, 3),
                                         ("003", True, None)):
            mp = os.path.join(trees_dir,
                              f"{plot_uid}_tree_{uid}_canopy_semantic_mask.tif")
            check(f"e2e: per-tree canopy mask written (tree {uid})",
                  os.path.exists(mp))
            with rasterio.open(mp) as src:
                mm = src.read(1)
                mtr = src.transform
            mcell = float(mtr.a); mx0 = float(mtr.c); my1 = float(mtr.f)
            mgrid = sm.RasterGrid(mx0, my1 - mm.shape[0] * mcell,
                                  mx0 + mm.shape[1] * mcell, my1, mcell)
            if want_interior:
                cc, rr = mgrid.to_pixel([0.0], [0.0])
                cc, rr = int(cc[0]), int(rr[0])
                check(f"e2e: tree {uid} mask fills the interior cell",
                      int(mm[rr, cc]) == 1)
            else:
                check("e2e: tree 002 mask has exactly its 3 point cells",
                      int((mm > 0).sum()) == npts, str(int((mm > 0).sum())))

        # hull JSON carries the footprint used for the mask
        hull_doc = json.load(open(os.path.join(
            place, f"{plot_uid}_raw_placement_canopy_convex_hull.json")))
        fp_by_tree = {h["tree_id"]: h for h in hull_doc["hulls"]}
        check("e2e: hull JSON records the 2D footprint per tree",
              len(fp_by_tree[1]["footprint_xy"]) == 4
              and fp_by_tree[1]["footprint_degenerate"] is False
              and fp_by_tree[2]["footprint_degenerate"] is True
              and fp_by_tree[2]["footprint_xy"] == [])

        # plot instance label: filled footprints, taller crown wins overlaps
        import restructure as rst
        pts5 = np.column_stack([pts, sem.astype(float), tid.astype(float)])
        lab, _ = rst._canopy_instance_mask(feat_path, pts5, anchor=None)
        cc, rr = grid.to_pixel([0.0], [0.0])
        cc, rr = int(cc[0]), int(rr[0])
        check("e2e: instance label fills the interior cell, taller tree 3 wins "
              "the overlap", int(lab[rr, cc]) == 3, str(int(lab[rr, cc])))
        oc, orr = grid.to_pixel([0.45], [0.45])   # inside tree 1, outside tree 3
        oc, orr = int(oc[0]), int(orr[0])
        check("e2e: instance label gives tree 1 the non-overlapping footprint",
              int(lab[orr, oc]) == 1, str(int(lab[orr, oc])))
        mc, mr = grid.to_pixel([6.0], [0.0])
        mc, mr = int(mc[0]), int(mr[0])
        check("e2e: instance label marks tree 2's point cell",
              int(lab[mr, mc]) == 2)
        mcol2, mrow2 = grid.to_pixel([5.5], [0.0])
        mcol2, mrow2 = int(mcol2[0]), int(mrow2[0])
        check("e2e: instance label leaves tree 2's midpoint empty",
              int(lab[mrow2, mcol2]) == 0)
    finally:
        if not keep:
            shutil.rmtree(tmp, ignore_errors=True)


def main(keep: bool = False) -> int:
    print("review-round regression suite (W22-W29 probes)")
    check_shrub_placement(keep)
    check_split_rule(keep)
    check_demo_contract(keep)
    check_stem_mask_rules(keep)
    check_chm_mask_rules(keep)
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
