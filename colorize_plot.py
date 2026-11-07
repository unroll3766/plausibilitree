"""Colorize a plot-level cloud by semantic rules and save a stripped RGB .las.

  ground (0)        -> black
  grass  (9)        -> light grey
  shrub  (8)        -> random colour per shrub_id (not black/brown/grey)
  trunk/branch(1,2) -> brown
  foliage(3)        -> random colour per tree_id (not black/brown/grey)

Usage: python colorize_plot.py <in.las> [out.las]
Keeps only xyz + RGB.
"""
import sys
import numpy as np
import laspy

import schema

BLACK = (0, 0, 0)
GREY = (200, 200, 200)
BROWN = (110, 70, 30)


def _rand_color(rng):
    """A random-ish colour that isn't black, brown, or grey (avoid dull/dark)."""
    while True:
        r, g, b = rng.integers(60, 256, size=3)
        mx, mn = int(max(r, g, b)), int(min(r, g, b))
        if mx < 80:                      # too dark / near black
            continue
        if mx - mn < 40:                 # too grey (low saturation)
            continue
        # brown ~ dark orange: r>g>b, moderate. Reject that region.
        if r > g > b and r < 170 and b < 90:
            continue
        return int(r), int(g), int(b)


def _palette(ids, rng):
    return {i: _rand_color(rng) for i in ids if i > 0}


def main():
    if len(sys.argv) < 2:
        sys.exit("usage: python colorize_plot.py <in.las> [out.las]")
    src = sys.argv[1]
    dst = sys.argv[2] if len(sys.argv) > 2 else src.rsplit(".", 1)[0] + "_rgb.las"
    rng = np.random.default_rng(0)

    las = laspy.read(src)
    xyz = np.column_stack([las.x, las.y, las.z])
    sem = np.asarray(las.semantic_code).astype(int)
    tid = np.asarray(las.tree_id).astype(int)
    sid = np.asarray(las.shrub_id).astype(int)

    rgb = np.zeros((len(xyz), 3), dtype=np.uint16)
    rgb[sem == schema.SEM_GROUND] = BLACK
    rgb[sem == schema.SEM_GRASS] = GREY
    rgb[np.isin(sem, (schema.SEM_TREE, schema.SEM_BRANCH))] = BROWN
    for i, col in _palette(np.unique(sid[sem == schema.SEM_SHRUB]), rng).items():
        rgb[(sem == schema.SEM_SHRUB) & (sid == i)] = col
    for i, col in _palette(np.unique(tid[sem == schema.SEM_FOLIAGE]), rng).items():
        rgb[(sem == schema.SEM_FOLIAGE) & (tid == i)] = col

    # 8-bit colours -> 16-bit LAS range
    rgb16 = (rgb.astype(np.uint16) * 257)

    header = laspy.LasHeader(point_format=2, version="1.2")   # fmt 2 has RGB
    out = laspy.LasData(header)
    out.x, out.y, out.z = xyz[:, 0], xyz[:, 1], xyz[:, 2]
    out.red, out.green, out.blue = rgb16[:, 0], rgb16[:, 1], rgb16[:, 2]
    out.write(dst)
    print(f"wrote {dst}  ({len(xyz):,} points)")


if __name__ == "__main__":
    main()
