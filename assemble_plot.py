"""
assemble_plot.py: assemble individually-generated trees into a plot scene and
resolve inter-tree collisions by trimming.

    python assemble_plot.py plot_11519.csv --seed 42

Reads the plot table and the per-tree outputs of generate_trees.py, applies a
random z-rotation to each tree about its own base, writes the merged scene,
then uses the skeletons to find intersections between trees and trims the
losing branch/stem according to the rules below.

Outputs, under the plot directory:

    plot_placement/<PLOT_UID>_raw.las        merged scene, rotated + placed
    plot_placement/<PLOT_UID>_scene.json     bounds, rotations, collision log
    plot_placement/<PLOT_UID>_trimmed.las    merged scene, surviving points,
                                             with branch_id + branch_order
    individual_trees/<TREE_UID>_trimmed.las  per-tree, surviving points, with
                                             branch_id + branch_order
    individual_trees/<TREE_UID>_trimmed_skeleton.json

Trimmed points are REMOVED, not flagged: the trimmed .las files contain only
the surviving points, each carrying its branch id. The scene JSON records how
many points trimming removed (under the legacy key name
"points_flagged_trimmed").

Trim rules (applied per tree in ascending TREE_UID order):
  branch x branch  -> shorter path-length-back-to-trunk is trimmed; ties broken
                      randomly (seeded)
  branch x trunk   -> the branch is trimmed
  trunk  x trunk   -> shorter stem is trimmed; ties broken by higher TREE_UID
When an element is trimmed, everything distal to its nearest joint goes with
it: the whole sub-tree hanging off that node, including child forks and any
branches attached above it.

Requires: numpy, pandas, laspy (via las.py)
"""

from __future__ import annotations

import argparse
import json
import math
import os

import numpy as np
import pandas as pd

from las import export_points_las, read_points_las, export_canonical_las, read_canonical_las
import layout
import schema
import tree_list

# Branch axes in the skeleton carry no radius. The simulator builds each branch
# as a ribbon of half-width 0.2 * branch length, so this fraction is the
# geometrically honest default -- but it is also the single most influential
# knob on how much of the canopy gets trimmed. Lower it for a sparser
# interpretation of "intersecting".
BRANCH_RADIUS_FRAC = 0.10
MIN_BRANCH_RADIUS = 0.02

# Per-point branch identity.
#
# Two forms are written, because one integer can't be both plot-global and fit
# a LAS field. FIA PLT_CN values are ~14 digits, so plot+tree+element does not
# fit in the int32 LAS extra dimension.
#
#   branch_id  (int32, in the LAS)  tree_id * ELEMENTS_PER_TREE + element_index
#              -- unique *within one plot*, which is the scope of a scene file.
#   branch_uid (string, in JSON)    "<plot_uid>_t<tree_id>_e<element>"
#              -- globally unique, no width limit.
#
# element_index is the skeleton element the point was assigned to (trunk
# sub-segment or branch), so points on the same branch share an id and the
# trunk is distinguishable from every branch. -1 elements (unassignable
# points, e.g. a class with no matching skeleton parts) get branch id 0.
ELEMENTS_PER_TREE = 100000  # ceiling on skeleton elements per tree; ~300 seen


def branch_id_int(tree_id: int, element_index: int) -> int:
    """Plot-unique integer branch id for the LAS branch_id dimension."""
    if element_index < 0:
        return 0
    return int(tree_id) * ELEMENTS_PER_TREE + int(element_index) + 1


def branch_uid_str(plot_uid: str, tree_id: int, element_index: int) -> str:
    """Globally-unique string branch id (plot + tree + element)."""
    if element_index < 0:
        return f"{plot_uid}_t{int(tree_id)}_eNA"
    return f"{plot_uid}_t{int(tree_id)}_e{int(element_index)}"


def _first_numeric(df, col):
    """First non-null value of a plot-center pass-through column (plot_lat/
    plot_lon/plot_elev), or None if the column isn't present in this plot's
    table (i.e. the FIA input had no LAT/LON/ELEV to pass through -- see
    fia_helper.preprocess_fia). All rows of one plot carry the same value."""
    if col not in df.columns:
        return None
    vals = pd.to_numeric(df[col], errors="coerce").dropna()
    return float(vals.iloc[0]) if len(vals) else None


# ---------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------

def rotate_z(xy: np.ndarray, angle_rad: float, cx: float, cy: float) -> np.ndarray:
    """Rotate (...,2) coordinates about the vertical axis through (cx, cy)."""
    c, s = math.cos(angle_rad), math.sin(angle_rad)
    dx = xy[..., 0] - cx
    dy = xy[..., 1] - cy
    out = np.empty_like(xy)
    out[..., 0] = cx + c * dx - s * dy
    out[..., 1] = cy + s * dx + c * dy
    return out


def segment_segment_distance(p0, p1, q0, q1):
    """Pairwise minimum distance between two sets of 3D segments.

    p0/p1 are (n, 3), q0/q1 are (m, 3); returns (n, m). Clamped
    parametric solution (Ericson, Real-Time Collision Detection), vectorised
    over all pairs and safe for degenerate (zero-length) segments.
    """
    d1 = (p1 - p0)[:, None, :]
    d2 = (q1 - q0)[None, :, :]
    r = p0[:, None, :] - q0[None, :, :]

    a = np.sum(d1 * d1, axis=-1)
    e = np.sum(d2 * d2, axis=-1)
    f = np.sum(d2 * r, axis=-1)
    b = np.sum(d1 * d2, axis=-1)
    c = np.sum(d1 * r, axis=-1)

    a = np.broadcast_to(a, f.shape)
    e = np.broadcast_to(e, f.shape)

    denom = a * e - b * b
    with np.errstate(divide="ignore", invalid="ignore"):
        s = np.where(denom > 1e-12, (b * f - c * e) / np.where(denom > 1e-12, denom, 1.0), 0.0)
        s = np.clip(s, 0.0, 1.0)
        t = np.where(e > 1e-12, (b * s + f) / np.where(e > 1e-12, e, 1.0), 0.0)
        t_c = np.clip(t, 0.0, 1.0)
        s = np.where(a > 1e-12, np.clip((b * t_c - c) / np.where(a > 1e-12, a, 1.0), 0.0, 1.0), 0.0)

    cp = p0[:, None, :] + d1 * s[..., None]
    cq = q0[None, :, :] + d2 * t_c[..., None]
    return np.linalg.norm(cp - cq, axis=-1)


def point_segment_distance(pts, p0, p1):
    """Distance from (n,3) points to one segment, as (n,)."""
    d = p1 - p0
    dd = float(d @ d)
    if dd < 1e-12:
        return np.linalg.norm(pts - p0, axis=1)
    t = np.clip(((pts - p0) @ d) / dd, 0.0, 1.0)
    return np.linalg.norm(pts - (p0 + t[:, None] * d), axis=1)


# ---------------------------------------------------------------------------
# Skeleton -> element graph
# ---------------------------------------------------------------------------

class TreeSkeleton:
    """Flattened, rotatable element graph built from a *_raw_skeleton.json.

    Elements are trunk sub-segments (one per polyline vertex pair) and
    branches (one segment each). Each element knows its parent, so trimming is
    just "mark this element and its descendants".
    """

    def __init__(self, skel: dict, tree_uid: str, tree_id: int,
                 branch_radius_frac: float = BRANCH_RADIUS_FRAC):
        self.raw = skel
        self.tree_uid = str(tree_uid)
        self.tree_id = int(tree_id)
        self.base_xy = np.asarray(skel["base_xy"], dtype=float)
        self.height = float(skel["height_m"])
        self.base_radius = float(skel["base_radius_m"])
        self.rotation_rad = 0.0

        self.kind = []          # 'trunk' | 'branch'
        self.poly = []          # polyline index (trunk) or -1
        self.node = []          # start-node index within polyline, or -1
        self.branch_idx = []    # index into skel['branches'], or -1
        self.parent = []
        p0, p1, r0, r1 = [], [], [], []

        polylines = [np.asarray(pl, dtype=float) for pl in skel["trunk"]["polylines"]]
        self.poly_lengths = [float(s["length_m"]) for s in skel["trunk"]["segments"]]

        for si, pl in enumerate(polylines):
            for i in range(len(pl) - 1):
                idx = len(self.kind)
                self.kind.append("trunk")
                self.poly.append(si)
                self.node.append(i)
                self.branch_idx.append(-1)
                p0.append(pl[i])
                p1.append(pl[i + 1])
                r0.append(self._trunk_radius(pl[i][2]))
                r1.append(self._trunk_radius(pl[i + 1][2]))
                # chain within the polyline; forks attach to the trunk
                # sub-segment their base sits on (nearest earlier element),
                # so cascaded forks hang off their true parent stem rather
                # than always off the main stem tip. For the original single
                # fork the base sits exactly on the main stem's tip, so this
                # resolves to the same parent as before.
                if i > 0:
                    self.parent.append(idx - 1)
                elif si == 0:
                    self.parent.append(None)
                else:
                    d = point_segment_distance(pl[0][None, :], p0[0], p1[0])
                    best, bestd = 0, d[0]
                    for k in range(1, idx):
                        dk = point_segment_distance(pl[0][None, :],
                                                     p0[k], p1[k])[0]
                        if dk < bestd:
                            best, bestd = k, dk
                    self.parent.append(best)

        self.n_trunk = len(self.kind)
        trunk_p0 = np.asarray(p0, dtype=float) if p0 else np.zeros((0, 3))
        trunk_p1 = np.asarray(p1, dtype=float) if p1 else np.zeros((0, 3))

        for bi, br in enumerate(skel.get("branches", [])):
            base = np.asarray(br["base"], dtype=float)
            tip = np.asarray(br["tip"], dtype=float)
            idx = len(self.kind)
            self.kind.append("branch")
            self.poly.append(-1)
            self.node.append(-1)
            self.branch_idx.append(bi)
            p0.append(base)
            p1.append(tip)
            rad = max(branch_radius_frac * float(br.get("length_m", 0.0)),
                      MIN_BRANCH_RADIUS)
            r0.append(rad)
            r1.append(rad)
            # attach to whichever trunk sub-segment the branch base sits on
            if self.n_trunk:
                d = point_segment_distance(base[None, :], trunk_p0[0], trunk_p1[0])
                best, bestd = 0, d[0]
                for k in range(1, self.n_trunk):
                    dk = point_segment_distance(base[None, :], trunk_p0[k], trunk_p1[k])[0]
                    if dk < bestd:
                        best, bestd = k, dk
                self.parent.append(best)
            else:
                self.parent.append(None)

        self.p0 = np.asarray(p0, dtype=float)
        self.p1 = np.asarray(p1, dtype=float)
        self.r0 = np.asarray(r0, dtype=float)
        self.r1 = np.asarray(r1, dtype=float)
        self.radius = np.maximum(self.r0, self.r1)
        self.kind = np.asarray(self.kind)
        self.poly = np.asarray(self.poly)
        self.node = np.asarray(self.node)
        self.branch_idx = np.asarray(self.branch_idx)
        self.trimmed = np.zeros(len(self.kind), dtype=bool)

        self.children = [[] for _ in range(len(self.parent))]
        for i, par in enumerate(self.parent):
            if par is not None:
                self.children[par].append(i)

        self._path_len = None
        self._order = None

    def _trunk_radius(self, z: float) -> float:
        return max(self.base_radius * (1.0 - z / self.height), 1e-3)

    # -- rotation ----------------------------------------------------------

    def rotate(self, angle_rad: float):
        cx, cy = self.base_xy
        self.p0[:, :2] = rotate_z(self.p0[:, :2], angle_rad, cx, cy)
        self.p1[:, :2] = rotate_z(self.p1[:, :2], angle_rad, cx, cy)
        self.rotation_rad = angle_rad

    # -- metrics -----------------------------------------------------------

    @property
    def lengths(self):
        return np.linalg.norm(self.p1 - self.p0, axis=1)

    def path_length_to_trunk(self, i: int) -> float:
        """Summed length from element i back to the first trunk element.

        For the current simulator every branch hangs directly off the stem, so
        this equals the branch length; it is written generically so that a
        higher-order branching model doesn't need this code changed.
        """
        if self._path_len is None:
            self._path_len = {}
        if i in self._path_len:
            return self._path_len[i]
        lens = self.lengths
        total, j = 0.0, i
        while j is not None and self.kind[j] == "branch":
            total += lens[j]
            j = self.parent[j]
        self._path_len[i] = total
        return total

    def stem_length(self, i: int) -> float:
        """Skeletal length of the stem an element belongs to."""
        si = int(self.poly[i])
        if 0 <= si < len(self.poly_lengths):
            return self.poly_lengths[si]
        return float(self.raw["trunk"]["total_length_m"])

    @property
    def order(self) -> np.ndarray:
        """Branching order per element.

        0 = the main stem. Order increments each time a new structure starts:
        a fork (the apical stems above a split) is 1st order, a branch on the
        main stem is 1st order, a branch on a fork is 2nd order, and so on.
        Continuing along the same stem or fork keeps the same order.

        Concretely, order = parent order + 1 for a branch or for the first
        element of a non-main trunk polyline (a fork), and = parent order for
        any element continuing along the polyline it started on. The main
        stem's root has order 0.

        Computed by walking the parent graph, so it is correct for arbitrary
        branching/forking depth. The current simulator forks the main stem at
        most once and hangs every branch directly off a stem, so in practice
        orders run 0 (main stem), 1 (forks and trunk branches) and 2 (branches
        on forks).
        """
        if self._order is not None:
            return self._order
        order = np.zeros(len(self.kind), dtype=np.int32)
        cache: dict[int, int] = {}

        def _o(i: int) -> int:
            if i in cache:
                return cache[i]
            par = self.parent[i]
            if par is None:                       # root of the main stem
                val = 0
            elif self.kind[i] == "branch":        # branch: new structure
                val = _o(par) + 1
            elif int(self.node[i]) == 0 and int(self.poly[i]) > 0:
                val = _o(par) + 1                 # first element of a fork
            else:                                 # continuing along a stem/fork
                val = _o(par)
            cache[i] = val
            return val

        for i in range(len(self.kind)):
            order[i] = _o(i)
        self._order = order
        return order

    def trim_subtree(self, i: int) -> int:
        """Trim element i and everything distal to it. Returns count trimmed."""
        stack, n = [i], 0
        while stack:
            j = stack.pop()
            if self.trimmed[j]:
                continue
            self.trimmed[j] = True
            n += 1
            stack.extend(self.children[j])
        return n

    def aabb(self, pad: float = 0.0):
        pts = np.vstack([self.p0, self.p1])
        rad = float(self.radius.max()) if len(self.radius) else 0.0
        lo = pts.min(axis=0) - rad - pad
        hi = pts.max(axis=0) + rad + pad
        return lo, hi


# ---------------------------------------------------------------------------
# Collision resolution
# ---------------------------------------------------------------------------

def _aabb_overlap(a, b) -> bool:
    (alo, ahi), (blo, bhi) = a, b
    return bool(np.all(alo <= bhi) and np.all(blo <= ahi))


def resolve_collisions(trees: list[TreeSkeleton], tolerance: float = 0.0,
                       rng: np.random.Generator | None = None) -> list[dict]:
    """Trim intersecting elements across trees. Returns a collision log.

    Trees are visited in ascending TREE_UID order; within a visit, every
    still-live element of that tree is tested against every still-live element
    of every other tree.
    """
    rng = rng or np.random.default_rng(0)
    order = sorted(range(len(trees)), key=lambda k: trees[k].tree_uid)
    boxes = [t.aabb() for t in trees]
    log = []

    for ai in order:
        A = trees[ai]
        for bi in range(len(trees)):
            if bi == ai or not _aabb_overlap(boxes[ai], boxes[bi]):
                continue
            B = trees[bi]

            live_a = np.flatnonzero(~A.trimmed)
            live_b = np.flatnonzero(~B.trimmed)
            if not len(live_a) or not len(live_b):
                continue

            d = segment_segment_distance(A.p0[live_a], A.p1[live_a],
                                         B.p0[live_b], B.p1[live_b])
            thresh = (A.radius[live_a][:, None] + B.radius[live_b][None, :]
                      + tolerance)
            hits = np.argwhere(d < thresh)
            if not len(hits):
                continue

            # deepest interpenetration first, so the decisive contact is
            # resolved before incidental ones nearby
            depth = thresh[hits[:, 0], hits[:, 1]] - d[hits[:, 0], hits[:, 1]]
            for h in hits[np.argsort(-depth)]:
                ea, eb = int(live_a[h[0]]), int(live_b[h[1]])
                if A.trimmed[ea] or B.trimmed[eb]:
                    continue
                loser, elem, reason = _pick_loser(A, ea, B, eb, rng)
                n = loser.trim_subtree(elem)
                log.append({
                    "tree_a": A.tree_uid, "element_a": ea, "kind_a": str(A.kind[ea]),
                    "tree_b": B.tree_uid, "element_b": eb, "kind_b": str(B.kind[eb]),
                    "trimmed_tree": loser.tree_uid, "trimmed_element": elem,
                    "elements_trimmed": int(n), "rule": reason,
                })
    return log


def _pick_loser(A: TreeSkeleton, ea: int, B: TreeSkeleton, eb: int, rng):
    """Apply the trim rules to one contact. Returns (tree, element, rule)."""
    ka, kb = A.kind[ea], B.kind[eb]

    if ka == "branch" and kb == "trunk":
        return A, ea, "branch_vs_trunk"
    if ka == "trunk" and kb == "branch":
        return B, eb, "branch_vs_trunk"

    if ka == "branch" and kb == "branch":
        la, lb = A.path_length_to_trunk(ea), B.path_length_to_trunk(eb)
        if abs(la - lb) < 1e-9:
            return (A, ea, "branch_tie_random") if rng.random() < 0.5 \
                else (B, eb, "branch_tie_random")
        return (A, ea, "shorter_branch") if la < lb else (B, eb, "shorter_branch")

    # trunk vs trunk
    la, lb = A.stem_length(ea), B.stem_length(eb)
    if abs(la - lb) < 1e-9:
        # tie -> higher tree id is cut
        return (A, ea, "trunk_tie_higher_uid") if A.tree_uid > B.tree_uid \
            else (B, eb, "trunk_tie_higher_uid")
    return (A, ea, "shorter_trunk") if la < lb else (B, eb, "shorter_trunk")


# ---------------------------------------------------------------------------
# Point <-> element association
# ---------------------------------------------------------------------------

def assign_points(points: np.ndarray, tree: TreeSkeleton,
                  chunk: int = 20000) -> np.ndarray:
    """Nearest skeleton element for every point, respecting the class split.

    Stem points (semantic_code 1 = tree, provisional) are matched against trunk
    elements only, and foliage points (semantic_code 3) against branch elements
    only, because that is how they were sampled -- letting foliage snap to the
    trunk would trim canopy whenever a stem was cut.
    """
    xyz = points[:, :3]
    sem = points[:, 3].astype(int)
    out = np.full(len(points), -1, dtype=int)

    trunk_idx = np.flatnonzero(tree.kind == "trunk")
    branch_idx = np.flatnonzero(tree.kind == "branch")

    for sem_val, idx in ((schema.SEM_TREE, trunk_idx),
                         (schema.SEM_FOLIAGE, branch_idx)):
        mask = sem == sem_val
        if not mask.any() or not len(idx):
            continue
        sel = np.flatnonzero(mask)
        for start in range(0, len(sel), chunk):
            part = sel[start:start + chunk]
            pts = xyz[part]
            best_d = np.full(len(part), np.inf)
            best_i = np.zeros(len(part), dtype=int)
            for k in idx:
                d = point_segment_distance(pts, tree.p0[k], tree.p1[k])
                better = d < best_d
                best_d[better] = d[better]
                best_i[better] = k
            out[part] = best_i
    return out


# ---------------------------------------------------------------------------
# Skeleton serialisation
# ---------------------------------------------------------------------------

def trimmed_skeleton_dict(tree: TreeSkeleton, n_points_trimmed: int) -> dict:
    """Copy of the raw skeleton with per-element trimmed flags added."""
    out = json.loads(json.dumps(tree.raw))
    out["rotation_deg"] = round(math.degrees(tree.rotation_rad), 4)
    out["trim_summary"] = {
        "elements_total": int(len(tree.trimmed)),
        "elements_trimmed": int(tree.trimmed.sum()),
        "points_trimmed": int(n_points_trimmed),
        "any_trimmed": bool(tree.trimmed.any()),
    }
    plot_uid = tree.raw.get("plot_uid", "plot")

    # rotated geometry, so the file still matches the point cloud beside it
    polys = out["trunk"]["polylines"]
    for si, pl in enumerate(polys):
        arr = np.asarray(pl, dtype=float)
        arr[:, :2] = rotate_z(arr[:, :2], tree.rotation_rad, *tree.base_xy)
        polys[si] = arr.tolist()
        node_trimmed = [False] * len(arr)
        # one skeleton element per vertex pair; store its branch id against the
        # sub-segment's start node so the JSON maps 1:1 onto the LAS branch_id.
        sub_branch_ids = [None] * (len(arr) - 1)
        sub_branch_uids = [None] * (len(arr) - 1)
        sel = np.flatnonzero((tree.kind == "trunk") & (tree.poly == si))
        for k in sel:
            ni = int(tree.node[k])
            sub_branch_ids[ni] = branch_id_int(tree.tree_id, k)
            sub_branch_uids[ni] = branch_uid_str(plot_uid, tree.tree_id, k)
            if tree.trimmed[k]:
                node_trimmed[ni] = True
                node_trimmed[ni + 1] = True
        seg = out["trunk"]["segments"][si]
        seg["node_trimmed"] = node_trimmed
        seg["trimmed"] = bool(any(node_trimmed))
        seg["subsegment_branch_ids"] = sub_branch_ids
        seg["subsegment_branch_uids"] = sub_branch_uids
        # every element of a polyline shares one order (0 for the main stem,
        # 1 for a fork); read it from the first element of this polyline
        seg["branch_order"] = int(tree.order[sel[0]]) if len(sel) else 0

    for k in np.flatnonzero(tree.kind == "branch"):
        bi = int(tree.branch_idx[k])
        br = out["branches"][bi]
        base = rotate_z(np.asarray(br["base"], dtype=float)[None, :2],
                        tree.rotation_rad, *tree.base_xy)[0]
        tip = rotate_z(np.asarray(br["tip"], dtype=float)[None, :2],
                       tree.rotation_rad, *tree.base_xy)[0]
        br["base"] = [float(base[0]), float(base[1]), float(br["base"][2])]
        br["tip"] = [float(tip[0]), float(tip[1]), float(br["tip"][2])]
        br["trimmed"] = bool(tree.trimmed[k])
        br["branch_id"] = branch_id_int(tree.tree_id, int(k))
        br["branch_uid"] = branch_uid_str(plot_uid, tree.tree_id, int(k))
        br["branch_order"] = int(tree.order[k])
    return out


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def assemble_plot(table_path: str, plot_dir: str | None = None,
                  seed: int = 0,
                  plot_radius_m: float | None = None,
                  branch_radius_frac: float = BRANCH_RADIUS_FRAC,
                  tolerance: float = 0.0,
                  plot_shape: str = "circle",
                  plot_length_m: float | None = None,
                  plot_width_m: float | None = None,
                  verbose: bool = True) -> dict:
    """Rotate, place, collide and trim every tree of one plot.

    Accepts a canonical tree-list JSON or a legacy SST parameter CSV."""
    plot_uid = os.path.splitext(os.path.basename(table_path))[0]
    if plot_dir is None:
        plot_dir = os.path.join(os.path.dirname(os.path.abspath(table_path)), plot_uid)
    trees_dir = os.path.join(plot_dir, "individual_trees")
    place_dir = os.path.join(plot_dir, "plot_placement")
    os.makedirs(place_dir, exist_ok=True)

    df = tree_list.trees_frame(table_path)
    if "successful_placement" in df.columns:
        df = df[df["successful_placement"].astype(str).str.lower()
                .isin(["true", "1", "yes", "t"])]
    df = df.sort_values("TREE_UID").reset_index(drop=True)

    rng = np.random.default_rng(seed)

    trees, clouds, rotations = [], [], {}
    for _, row in df.iterrows():
        uid = str(row["TREE_UID"])
        las_path = layout.cloud_file(trees_dir, f"{uid}_raw")
        skel_path = os.path.join(trees_dir, f"{uid}_raw_skeleton.json")
        if not (os.path.exists(las_path) and os.path.exists(skel_path)):
            if verbose:
                print(f"  {uid}: no generated output, skipping")
            continue

        xyz_t, cols_t = read_canonical_las(las_path, level="tree")
        pts = np.column_stack([xyz_t, cols_t["semantic_code"].astype(float),
                               cols_t["tree_id"].astype(float)])
        # pts columns: x y z semantic_code tree_id  (provisional semantic)
        skel = json.load(open(skel_path))
        tree = TreeSkeleton(skel, uid, int(skel.get("tree_id", 0)),
                            branch_radius_frac=branch_radius_frac)

        # rotate about the tree's own base, then re-anchor on the table X/Y
        angle = float(rng.uniform(0.0, 2 * math.pi))
        cx, cy = float(skel["base_xy"][0]), float(skel["base_xy"][1])
        pts[:, :2] = rotate_z(pts[:, :2], angle, cx, cy)
        tree.rotate(angle)
        dx, dy = float(row["X"]) - cx, float(row["Y"]) - cy
        if abs(dx) > 1e-9 or abs(dy) > 1e-9:
            pts[:, 0] += dx
            pts[:, 1] += dy
            tree.p0[:, 0] += dx; tree.p0[:, 1] += dy
            tree.p1[:, 0] += dx; tree.p1[:, 1] += dy
            tree.base_xy = np.array([row["X"], row["Y"]], dtype=float)

        trees.append(tree)
        clouds.append(pts)
        rotations[uid] = round(math.degrees(angle), 3)

    if not trees:
        raise RuntimeError(f"no generated trees found under {trees_dir}")

    # ---- merged raw scene ------------------------------------------------
    raw = np.vstack(clouds)                 # x y z semantic_code tree_id
    raw_path = layout.cloud_file(place_dir, f"{plot_uid}_raw")
    export_canonical_las(
        raw_path, raw[:, :3],
        {"tree_id": raw[:, 4].astype(np.int32),
         "semantic_code": raw[:, 3].astype(np.int8)},
        level="plot")

    if plot_radius_m is None:
        plot_radius_m = float(np.hypot(df["X"].astype(float),
                                       df["Y"].astype(float)).max())
    lo = raw[:, :3].min(axis=0)
    hi = raw[:, :3].max(axis=0)

    # ---- collisions ------------------------------------------------------
    log = resolve_collisions(trees, tolerance=tolerance, rng=rng)

    # ---- per-tree trimmed outputs ---------------------------------------
    # Points are flagged, never removed: every raw point survives into the
    # trimmed cloud carrying a trimmed flag (0/1) and a branch id. Consumers
    # decide what to do with trimmed points; this stage does not clip.
    total_trimmed = 0
    _branch_id_remap = {}          # full branch_id_int -> compact int32 id
    scene_xyz, scene_tid, scene_sem, scene_code, scene_bid = [], [], [], [], []
    for tree, pts in zip(trees, clouds):
        elem = assign_points(pts, tree)                 # skeleton element / point

        valid = elem >= 0
        keep = np.ones(len(pts), dtype=bool)
        # trimmed points are now REMOVED, not flagged
        keep[valid] = ~tree.trimmed[elem[valid]].astype(bool)
        total_trimmed += int((~keep).sum())

        e = elem[keep]
        ev = e >= 0
        kpts = pts[keep]

        # branch id per surviving point (0 for non-tree / unassigned).
        # The full branch_id_int() value can exceed the LAS dim's width, so
        # the point-cloud branch_id is a compact plot-unique sequential id
        # (0 reserved for non-branch). The skeleton JSON keeps the full ids.
        # The remap counter is cumulative across every tree in the plot, so
        # it needs the same int32 headroom as tree_id/shrub_id -- a plot with
        # many large, many-branched trees can assign more than int16's 32767
        # distinct branch ids well before running out of trees.
        bid_full = np.zeros(len(kpts), dtype=np.int64)
        bid_full[ev] = [branch_id_int(tree.tree_id, ee) for ee in e[ev]]
        bid = np.zeros(len(kpts), dtype=np.int32)
        for full in np.unique(bid_full[bid_full > 0]):
            if full not in _branch_id_remap:
                _branch_id_remap[full] = len(_branch_id_remap) + 1
            bid[bid_full == full] = _branch_id_remap[full]

        # branch_code = branch_order + 1 on any assigned tree point; 0 elsewhere.
        # foliage points sit on branch elements, so they carry that branch's code.
        code = np.zeros(len(kpts), dtype=np.int8)
        code[ev] = (tree.order[e[ev]] + 1).astype(np.int8)

        # semantic refinement:
        #   stem points (provisional 1): trunk element -> 1 (trunk stem);
        #     fork element (still trunk-kind but branch_order>=1) -> 2 (branch stem)
        #   foliage points (provisional 3): stay 3 (foliage)
        sem = kpts[:, 3].astype(np.int8)                # provisional (1 or 3)
        is_stem = sem == schema.SEM_TREE
        stem_order = np.zeros(len(kpts), dtype=int)
        stem_order[ev] = tree.order[e[ev]]
        sem[is_stem & (stem_order == 0)] = schema.SEM_TREE     # 1 main trunk
        sem[is_stem & (stem_order >= 1)] = schema.SEM_BRANCH   # 2 fork/branch stem

        tid = kpts[:, 4].astype(np.int32)

        export_canonical_las(
            layout.cloud_file(trees_dir, f"{tree.tree_uid}_trimmed"),
            kpts[:, :3],
            {"tree_id": tid, "semantic_code": sem,
             "branch_code": code, "branch_id": bid},
            level="tree")
        with open(os.path.join(trees_dir,
                               f"{tree.tree_uid}_trimmed_skeleton.json"), "w") as f:
            json.dump(trimmed_skeleton_dict(tree, int((~keep).sum())), f, indent=2)

        scene_xyz.append(kpts[:, :3])
        scene_tid.append(tid)
        scene_sem.append(sem)
        scene_code.append(code)
        scene_bid.append(bid)
        if verbose and (~keep).any():
            print(f"  {tree.tree_uid}: {int(tree.trimmed.sum())} elements, "
                  f"{int((~keep).sum())} points removed (trimmed)")

    # merged trimmed scene: collision-culled geometry, canonical schema.
    export_canonical_las(
        layout.cloud_file(place_dir, f"{plot_uid}_trimmed"),
        np.vstack(scene_xyz),
        {"tree_id": np.concatenate(scene_tid),
         "semantic_code": np.concatenate(scene_sem),
         "branch_code": np.concatenate(scene_code),
         "branch_id": np.concatenate(scene_bid)},
        level="plot")

    scene = {
        "plot_uid": plot_uid,
        "n_trees": len(trees),
        "plot_shape": plot_shape,
        "plot_radius_m": float(plot_radius_m),
        "plot_length_m": (float(plot_length_m) if plot_length_m is not None else None),
        "plot_width_m": (float(plot_width_m) if plot_width_m is not None else None),
        "plot_lat": _first_numeric(df, "plot_lat"),
        "plot_lon": _first_numeric(df, "plot_lon"),
        "plot_elev": _first_numeric(df, "plot_elev"),
        "bounds": {"min": [float(v) for v in lo], "max": [float(v) for v in hi]},
        "seed": int(seed),
        "branch_radius_frac": float(branch_radius_frac),
        "contact_tolerance_m": float(tolerance),
        "rotations_deg": rotations,
        "points_raw": int(len(raw)),
        "points_flagged_trimmed": int(total_trimmed),
        "branch_id_encoding": {
            "las_field": "branch_id",
            "formula": "tree_id * {} + element_index + 1".format(ELEMENTS_PER_TREE),
            "scope": "unique within this plot",
            "string_form": "<plot_uid>_t<tree_id>_e<element_index> (see skeletons)",
            "unassigned_value": 0,
        },
        "branch_order_encoding": {
            "las_field": "branch_order",
            "definition": "0 = main stem; a fork or a branch is parent order + "
                          "1 (fork = 1, branch on trunk = 1, branch on fork = 2)",
            "note": "current simulator produces orders 0, 1 and 2",
            "unassigned_value": 0,
        },
        "n_collisions": len(log),
        "collisions": log,
    }
    with open(os.path.join(place_dir, f"{plot_uid}_scene.json"), "w") as f:
        json.dump(scene, f, indent=2)

    if verbose:
        print(f"{plot_uid}: {len(trees)} trees, {len(raw):,} points, "
              f"{len(log)} collisions, {total_trimmed:,} points trimmed "
              f"({100 * total_trimmed / max(len(raw), 1):.1f}%)")
    return scene


def main():
    ap = argparse.ArgumentParser(description="Assemble and trim a plot scene")
    ap.add_argument("csv", nargs="+", help="plot table csv(s)")
    ap.add_argument("-d", "--plot-dir", default=None,
                    help="directory holding individual_trees/ (default: from csv name)")
    ap.add_argument("-s", "--seed", type=int, default=0)
    ap.add_argument("-r", "--plot-radius", type=float, default=None)
    ap.add_argument("--plot-shape", choices=["circle", "rectangle"],
                    default="circle")
    ap.add_argument("--length", type=float, default=None)
    ap.add_argument("--width", type=float, default=None)
    ap.add_argument("--branch-radius-frac", type=float, default=BRANCH_RADIUS_FRAC)
    ap.add_argument("--tolerance", type=float, default=0.0,
                    help="extra clearance required between elements, metres")
    args = ap.parse_args()

    for path in args.csv:
        assemble_plot(path, plot_dir=args.plot_dir, seed=args.seed,
                      plot_radius_m=args.plot_radius,
                      branch_radius_frac=args.branch_radius_frac,
                      tolerance=args.tolerance,
                      plot_shape=args.plot_shape,
                      plot_length_m=args.length,
                      plot_width_m=args.width)


if __name__ == "__main__":
    main()
