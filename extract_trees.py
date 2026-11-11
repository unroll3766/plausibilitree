"""
extract_trees.py: the opt-in individual-extraction stage (spec Section 5,
item W6).

The core outputs are one plot cloud plus per-object JSONs; nothing per-tree
or per-shrub is written by a core run. This convenience stage crops
per-object clouds back out of the plot cloud on demand:

  - trees: `tree_id` slices, keyed by the skeletons' TREE_UID;
  - shrubs: `shrub_id` slices of the shrub returns (understory clouds).

`--source` picks the plot cloud to slice (understory -> terrain -> trimmed
-> raw resolution); each crop is written to the canonical per-object `las/`
position under that source's name, so the scene links, normalize and the
editor read it unchanged. `--surroundings N` widens every crop to an
N-metre window around the object's own points so the crop carries
background; without it a crop is the object's own points only.

Usage:
    python extract_trees.py plots/plot_demo --source understory
    python extract_trees.py plots/plot_demo --surroundings 2
"""

from __future__ import annotations

import argparse
import json
import os

import numpy as np

import layout
from las import read_canonical_las, export_canonical_las

# resolution order for --source understory (the default): the
# furthest-processed plot cloud present
_RESOLUTION = ("understory", "terrain", "trimmed", "raw")


def _resolve_source(plot_dir: str, plot_uid: str, source: str) -> str | None:
    """The concrete plot-cloud variant to slice."""
    if source != "understory":
        return source
    for cand in _RESOLUTION:
        p = layout.plot_path(plot_dir, plot_uid, f"{cand}_las", create=False)
        if os.path.exists(p):
            return cand
    return None


def _skeleton_index(plot_dir: str, plot_uid: str) -> dict:
    """{tree_id: tree_uid} from the canonical trimmed skeletons."""
    out = {}
    it = os.path.join(plot_dir, layout.INDIVIDUAL_TREES)
    if not os.path.isdir(it):
        return out
    for name in sorted(os.listdir(it)):
        # dir names are the full <PLOT_UID>_<TREE_UID> tags; the skeleton's
        # own tree_uid is authoritative (uids may contain underscores)
        tdir = os.path.join(it, name)
        cand = os.path.join(tdir, "textual",
                            f"{name}_trimmed_skeleton.json")
        if not os.path.exists(cand):
            continue
        sk = json.load(open(cand))
        uid = str(sk.get("tree_uid") or name.rsplit("_", 1)[-1])
        if sk.get("tree_id") is not None:
            out[int(sk["tree_id"])] = uid
    return out


def _slice_mask(cols, xyz, *, tree_id=None, shrub_id=None):
    if tree_id is not None:
        m = cols["tree_id"].astype(int) == int(tree_id)
    elif shrub_id is not None:
        m = (cols["shrub_id"].astype(int) == int(shrub_id)) & \
            (cols["semantic_code"].astype(int) == 8)
    else:
        raise ValueError("one of tree_id / shrub_id is required")
    if not m.any():
        return m, None
    box = (float(xyz[m, 0].min()), float(xyz[m, 1].min()),
           float(xyz[m, 0].max()), float(xyz[m, 1].max()))
    return m, box


def _surroundings_mask(xyz, box, radius):
    x0, y0, x1, y1 = box
    return ((xyz[:, 0] >= x0 - radius) & (xyz[:, 0] <= x1 + radius) &
            (xyz[:, 1] >= y0 - radius) & (xyz[:, 1] <= y1 + radius))


def _export(out_path, xyz, cols, m):
    sub = {k: v[m] for k, v in cols.items()}
    export_canonical_las(out_path, xyz[m], sub, level="tree")


def extract_plot(plot_dir, plot_uid, source="understory", surroundings=0.0,
                 verbose=True) -> int:
    """Slice per-object clouds out of the resolved plot cloud. Returns the
    number of objects written (trees + shrubs)."""
    resolved = _resolve_source(plot_dir, plot_uid, source)
    if resolved is None:
        raise FileNotFoundError(
            f"{plot_uid}: no plot cloud to extract from (looked for "
            + ", ".join(f"{v}_las" for v in _RESOLUTION) + " in "
            f"{layout.PLOT_PLACEMENT}/las/)")
    las_path = layout.plot_path(plot_dir, plot_uid, f"{resolved}_las",
                                create=False)
    xyz, cols = read_canonical_las(las_path, level="plot")

    written = 0
    new_links: dict[str, dict] = {"trees": {}, "shrubs": {}}

    # ---- trees: tree_id slices keyed by UID ----
    for tid, uid in sorted(_skeleton_index(plot_dir, plot_uid).items()):
        m, box = _slice_mask(cols, xyz, tree_id=tid)
        if box is not None and surroundings > 0:
            m = _surroundings_mask(xyz, box, surroundings)
        out = layout.tree_path(plot_dir, plot_uid, uid, f"{resolved}_las")
        _export(out, xyz, cols, m)
        new_links["trees"][uid] = layout.tree_rel(plot_dir, plot_uid, uid,
                                                  f"{resolved}_las")
        written += 1
        if verbose:
            print(f"      tree {uid}: {int(m.sum())} points ({resolved})")

    # ---- shrubs: shrub_id slices of the shrub returns (understory only) ----
    if resolved == "understory":
        sem = cols["semantic_code"].astype(int)
        sid_all = cols["shrub_id"].astype(int)
        sids = np.unique(sid_all[sem == 8])
        for s in sorted(int(v) for v in sids if int(v) > 0):
            uid = f"{s:03d}"
            m, box = _slice_mask(cols, xyz, shrub_id=s)
            if m.any() and surroundings > 0:
                m = _surroundings_mask(xyz, box, surroundings)
            if not m.any():
                continue
            out = layout.shrub_path(plot_dir, plot_uid, uid, f"{resolved}_las")
            _export(out, xyz, cols, m)
            new_links["shrubs"][uid] = layout.shrub_rel(plot_dir, plot_uid,
                                                        uid, f"{resolved}_las")
            written += 1
            if verbose:
                print(f"      shrub {uid}: {int(m.sum())} points")

    if written:
        _update_scene(plot_dir, plot_uid, resolved, new_links)
    if verbose:
        print(f"  {plot_uid}: extracted {written} object cloud(s) "
              f"from {resolved}")
    return written


def _update_scene(plot_dir, plot_uid, resolved, new_links):
    """Add the extracted clouds to the scene's per-object links (the scene
    links every artifact that exists; nothing is overwritten)."""
    scene_path = layout.plot_path(plot_dir, plot_uid, "scene", create=False)
    if not os.path.exists(scene_path):
        return
    with open(scene_path) as fh:
        scene = json.load(fh)
    trees = scene.setdefault("trees", {})
    for uid, rel in new_links["trees"].items():
        trees.setdefault(uid, {})[f"{resolved}_las"] = rel
    shrubs = scene.setdefault("shrubs", {})
    for uid, rel in new_links["shrubs"].items():
        shrubs.setdefault(uid, {})[f"{resolved}_las"] = rel
    with open(scene_path, "w") as fh:
        json.dump(scene, fh, indent=2)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("csv", nargs="+", help="plot table csv(s)")
    ap.add_argument("-d", "--plot-dir", default=None)
    ap.add_argument("--source", default="understory",
                    choices=["understory", "terrain", "trimmed", "raw"],
                    help="plot cloud to slice; 'understory' resolves to the "
                         "furthest-processed variant present")
    ap.add_argument("--surroundings", type=float, default=0.0,
                    help="metres of surrounding-plot context to include in "
                         "each crop (default: the object's own points only)")
    args = ap.parse_args()
    n_total = 0
    for c in args.csv:
        plot_uid = os.path.splitext(os.path.basename(c))[0]
        if args.plot_dir is not None:
            plot_dir = args.plot_dir
        else:
            plot_dir = os.path.join(os.path.dirname(c), plot_uid)
        n_total += extract_plot(plot_dir, plot_uid, source=args.source,
                                surroundings=args.surroundings)
    print(f"extract_trees: {n_total} object cloud(s) written")


if __name__ == "__main__":
    main()
