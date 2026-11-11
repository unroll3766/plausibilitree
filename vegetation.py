"""
vegetation.py: procedural understory (shrubs + grass) for a placed plot.

Two ground-cover generators, run once per plot after trees are placed and the
terrain has been draped:

    gen_shrub_layer -- discrete jittered-ellipsoid shrubs, one mesh per shrub
    gen_grass_layer -- continuous procedural point field (Poisson thinning of a
                       diamond-square "vigor" field), no per-plant mesh

Both avoid/attenuate around trees (and, for grass, shrubs) using canopy radii.

Object-type convention lives in the `stem_foliage` ("type") dim, extending the
terrain scheme: 1 stem, 0 foliage, 2 shrub, 3 grass, 255 ground. The LAS
`classification` field stays authoritative and is set by height above ground
(3 low veg <=0.5 m, 4 med veg 0.5-2 m, 5 high veg >2 m) -- so a shrub is classed
by its height like any real return, not by a fixed code.

Requires: numpy, scipy, trimesh (via las stub in this sandbox), generate_dtm.
"""

from __future__ import annotations

import numpy as np
from scipy.spatial import cKDTree

# object-type codes for the stem_foliage / "type" dim
import schema

# object-type codes = canonical semantic_code values
TYPE_FOLIAGE = schema.SEM_FOLIAGE   # 3
TYPE_STEM = schema.SEM_TREE         # 1 (trunk); branch stems are SEM_BRANCH (2)
TYPE_BRANCH = schema.SEM_BRANCH     # 2
TYPE_SHRUB = schema.SEM_SHRUB       # 4
TYPE_GRASS = schema.SEM_GRASS       # 5
TYPE_GROUND = schema.SEM_GROUND     # 6


# ---------------------------------------------------------------------------
# shrub mesh primitive
# ---------------------------------------------------------------------------

def gen_shrub_ellipsoid_mesh(height, width, jitter_frac=0.4, subdivisions=2,
                             rng=None):
    """A single warped-ellipsoid shrub mesh, centred at the local origin.

    Icosphere -> anisotropic scale (xy by width/2, z by height) -> per-vertex
    radial jitter so shrubs aren't identical smooth ellipsoids. The mesh sits
    on z in [0, height] with its base at z=0 (caller drapes it onto the DTM).
    """
    import trimesh
    rng = rng or np.random.default_rng()
    ico = trimesh.creation.icosphere(subdivisions=subdivisions, radius=1.0)
    v = np.asarray(ico.vertices, dtype=float).copy()

    # radial jitter: scale each unit vertex radius by (1 + N(0, sigma))
    jitter = 1.0 + rng.normal(0.0, jitter_frac * 0.5, size=len(v))
    v *= np.clip(jitter, 0.3, 2.0)[:, None]

    # anisotropic scale to an ellipsoid, then lift so the base sits at z=0
    v[:, 0] *= width / 2.0
    v[:, 1] *= width / 2.0
    v[:, 2] = v[:, 2] * (height / 2.0) + height / 2.0
    v[:, 2] = np.clip(v[:, 2], 0.0, None)

    mesh = trimesh.Trimesh(vertices=v, faces=np.asarray(ico.faces), process=False)
    return mesh


def ellipsoid_primitive_volume(height, width):
    """Nominal ellipsoid volume (4/3)*pi*a*b*c for semi-axes (w/2, w/2, h/2).

    This is the un-jittered primitive the shrub was scaled from, before the
    per-vertex radial warp -- the analytic 'primitive' volume.
    """
    a = width / 2.0
    c = height / 2.0
    return (4.0 / 3.0) * np.pi * a * a * c


# ---------------------------------------------------------------------------
# shrub layer
# ---------------------------------------------------------------------------

def gen_shrub_layer(plot_bounds, cell_size, shrub_cover, height_range,
                    width_range, tree_positions, tree_canopy_r,
                    points_per_shrub=300, jitter_frac=0.4, margin=0.0,
                    dtm=None, seed=None,
                    plot_shape="rectangle", plot_radius=None):
    """Generate the shrub point layer.

    Each shrub's height and width are drawn independently from their ranges
    (width_range is a (min, max) in metres, like height_range).

    Candidate cells are constrained to the plot's own footprint: for a
    circle (plot_shape="circle", plot_radius in metres) only cells whose
    centre lies inside the disc are seeded. Placing shrubs over the full
    bounding square would put corner shrubs wholly outside the plot, and
    the point-level disc clip would strip every one of their points while
    their records persisted (the 2026-11 review's 0-point shrubs).

    Returns (points, centers, widths):
      points  : (N, 3) x, y, z  (draped onto the DTM)
      shrub_id: (N,)   1-based id per point
      centers : (Ns, 2) shrub centres  (for grass attenuation)
      widths  : (Ns,)  shrub radii     (for grass attenuation)
    """
    rng = np.random.default_rng(seed)
    xmin, xmax, ymin, ymax = plot_bounds

    # 1. cell-centre grid
    nx = max(1, int(np.floor((xmax - xmin) / cell_size)))
    ny = max(1, int(np.floor((ymax - ymin) / cell_size)))
    cx = xmin + (np.arange(nx) + 0.5) * cell_size
    cy = ymin + (np.arange(ny) + 0.5) * cell_size
    gx, gy = np.meshgrid(cx, cy)
    cells = np.column_stack([gx.ravel(), gy.ravel()])

    # 2. exclude cells within (canopy_r + margin) of any tree
    keep = np.ones(len(cells), dtype=bool)
    if len(tree_positions):
        tree = cKDTree(np.asarray(tree_positions, dtype=float))
        dist, idx = tree.query(cells, k=1)
        keep &= dist > (np.asarray(tree_canopy_r)[idx] + margin)

    # 2b. constrain candidate cells to the plot footprint (spec 3.4): a
    # circular plot seeds shrubs only inside the disc, never in the
    # bounding square's corners
    if plot_shape == "circle" and plot_radius:
        r = float(plot_radius)
        keep &= cells[:, 0] ** 2 + cells[:, 1] ** 2 <= r * r

    # 3. keep with probability shrub_cover
    keep &= rng.random(len(cells)) <= shrub_cover
    kept = cells[keep]
    if not len(kept):
        return (np.zeros((0, 3)), np.zeros(0, dtype=int),
                np.zeros((0, 2)), np.zeros(0), {})

    # 4-5. one shrub per kept cell
    all_pts, all_ids, centers, widths = [], [], [], []
    shrub_geom = {}          # shrub_id -> dict of geometry + volumes
    for sid, (ccx, ccy) in enumerate(kept, start=1):
        jx = rng.uniform(-jitter_frac, jitter_frac) * cell_size
        jy = rng.uniform(-jitter_frac, jitter_frac) * cell_size
        sx, sy = ccx + jx, ccy + jy
        h = rng.uniform(*height_range)
        w = rng.uniform(*width_range)
        mesh = gen_shrub_ellipsoid_mesh(h, w, jitter_frac=jitter_frac, rng=rng)
        # trimesh >= 5 ignores numpy seeds in Trimesh.sample, so the seed is
        # drawn explicitly from this function's seeded rng; without it the
        # same understory seed would scatter each shrub's points differently
        # on every run
        pts = mesh.sample(points_per_shrub,
                          seed=int(rng.integers(1, 2**31 - 1)))
        # drape per point onto the DTM at each point's own xy (matching
        # the grass layer): anchoring the whole shrub at its CENTER's
        # elevation buries the flanks into upslope terrain and yields
        # negative height-above-ground for spread-out shrubs
        pts = pts + np.array([sx, sy, 0.0])
        center_z = float(dtm.sample(np.array([sx]), np.array([sy]))[0]) \
            if dtm is not None else 0.0
        if dtm is not None:
            pts[:, 2] += dtm.sample(pts[:, 0], pts[:, 1])
        all_pts.append(pts)
        all_ids.append(np.full(len(pts), sid, dtype=int))
        centers.append((sx, sy))
        widths.append(w / 2.0)
        # mesh volume: the actual (jittered) primitive that produced the points.
        # Prefer trimesh's own volume; fall back to the convex-hull volume of the
        # mesh vertices if the mesh isn't watertight / volume is unavailable.
        mesh_vol = float("nan")
        try:
            mv = abs(float(mesh.volume))
            if np.isfinite(mv) and mv > 0:
                mesh_vol = mv
        except Exception:
            pass
        if not np.isfinite(mesh_vol) or mesh_vol <= 0:
            try:
                from scipy.spatial import ConvexHull
                mesh_vol = float(ConvexHull(np.asarray(mesh.vertices)).volume)
            except Exception:
                mesh_vol = float("nan")
        shrub_geom[sid] = {
            "shrub_id": int(sid),
            "center_xy": [float(sx), float(sy)],
            "base_z_m": center_z,   # ground at the shrub centre; points
                                    # drape per-point onto the DTM
            "height_m": float(h),
            "width_m": float(w),
            "primitive_volume_m3": float(ellipsoid_primitive_volume(h, w)),
            "mesh_volume_m3": mesh_vol,
        }

    pts = np.vstack(all_pts)
    ids = np.concatenate(all_ids)

    # 6. clip shrub points inside a tree's stem-taper volume (trunk cylinder
    # approximation: within the tree's own point footprint radius at that
    # height is overkill here; use a small trunk radius so shrubs may nestle
    # under canopy but not intersect the stem)
    # (kept simple: no hard clip needed since cells near trees were excluded;
    # attenuation/exclusion already handled placement)

    return pts, ids, np.asarray(centers), np.asarray(widths), shrub_geom


# ---------------------------------------------------------------------------
# grass vigor field
# ---------------------------------------------------------------------------

def gen_vigor_field(plot_bounds, patch_scale, seed=None):
    """Continuous [0,1] vigor field D(x,y) via the diamond-square generator.

    Returns a DTM-like object with a .sample(x, y) bilinear query, reusing the
    terrain DTM class so the interpolation matches the rest of the pipeline.
    patch_scale (m) sets the field's correlation length: larger -> coarser
    grid -> larger patches.
    """
    from terrain import DTM
    xmin, xmax, ymin, ymax = plot_bounds
    width = xmax - xmin
    height = ymax - ymin
    # coarse cell ~ patch_scale; relief 1.0 so the normalized field is ~[0,1]
    field = DTM.generate(width, height, cell_size=max(patch_scale, 1e-3),
                         bumpiness=0.6, relief=1.0, seed=seed)
    # normalize to exactly [0,1]
    g = field.grid
    lo, hi = float(g.min()), float(g.max())
    field.grid = (g - lo) / (hi - lo) if hi > lo else np.zeros_like(g)
    field.xll, field.yll = xmin, ymin
    return field


# ---------------------------------------------------------------------------
# grass layer
# ---------------------------------------------------------------------------

def gen_grass_layer(plot_bounds, cover, height_range, patch_scale,
                    tree_positions, tree_canopy_r, shrub_positions=None,
                    shrub_radii=None, canopy_attenuation=0.3,
                    candidate_density=200.0, dtm=None, seed=None):
    """Generate the grass point layer by Poisson thinning of a vigor field.

    Returns (N, 3) x, y, z.
    """
    rng = np.random.default_rng(seed)
    xmin, xmax, ymin, ymax = plot_bounds
    area = (xmax - xmin) * (ymax - ymin)

    field = gen_vigor_field(plot_bounds, patch_scale, seed=seed)

    # 2. candidates
    n = max(0, int(round(candidate_density * area)))
    if n == 0:
        return np.zeros((0, 3))
    x = rng.uniform(xmin, xmax, n)
    y = rng.uniform(ymin, ymax, n)

    # 3. local vigor
    d_local = np.clip(field.sample(x, y), 0.0, 1.0)

    # 4. attenuation: canopy_attenuation ** (num overlapping footprints)
    n_over = np.zeros(n, dtype=int)
    xy = np.column_stack([x, y])
    if len(tree_positions):
        tree = cKDTree(np.asarray(tree_positions, dtype=float))
        dist, idx = tree.query(xy, k=1)
        n_over += (dist <= np.asarray(tree_canopy_r)[idx]).astype(int)
    if shrub_positions is not None and len(shrub_positions):
        shrub = cKDTree(np.asarray(shrub_positions, dtype=float))
        dist, idx = shrub.query(xy, k=1)
        n_over += (dist <= np.asarray(shrub_radii)[idx]).astype(int)
    atten = np.power(canopy_attenuation, n_over)

    # 5. thinning
    accept = np.clip(cover * d_local * atten, 0.0, 1.0)
    keep = rng.random(n) < accept
    x, y, d_local = x[keep], y[keep], d_local[keep]
    if not len(x):
        return np.zeros((0, 3))

    # 6. height from vigor, u^2-biased blade sampling
    h = height_range[0] + (height_range[1] - height_range[0]) * d_local
    u = rng.random(len(x)) ** 2
    ground = dtm.sample(x, y) if dtm is not None else np.zeros(len(x))
    z = ground + h * u

    return np.column_stack([x, y, z])
