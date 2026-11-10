"""
fia_helper.py: Convert FIA (Forest Inventory and Analysis) tree records into a
per-plot parameter table that PlausibiliTree can consume directly.

Pipeline (see fia_to_sst_plot):

    preprocess -> expand by TPA -> crown width -> placement -> postprocess

Each stage is also usable on its own so that individual steps can be swapped
out.

If the input carries FIA PLOT fields LAT/LON/ELEV, they are passed through
unchanged as plot_lat/plot_lon/plot_elev columns on the output (plot-center
values, no offset) -- see preprocess_fia and PLOT_CENTER_COLUMNS.

Requires: numpy, pandas
Companion data files (same directory by default):
    data/functional_groups.json  parameter dictionary keyed by functional
                                 group (spec Section 8.6)
    data/spcd_func_group.csv     SPCD -> functional group lookup
"""

from __future__ import annotations

import json
import os
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Unit conversions
# ---------------------------------------------------------------------------

IN_TO_M = 0.0254        # FIA DIA is in inches
FT_TO_M = 0.3048        # FIA HT / ACTUALHT are in feet
SQM_PER_ACRE = 4046.86  # TPA_UNADJ is trees per acre

# NOTE: the spec for this module gave the height factor as 3.6576 (= 0.3048*12,
# i.e. a feet->meters conversion applied to a value already assumed to be in
# inches).  FIA stores HT in feet, so 3.6576 makes a 40 ft tree 146 m tall.
# The correct factor is used by default; override HT_TO_M below (or pass
# ht_factor= to preprocess_fia) if the 3.6576 behaviour is actually wanted.
HT_TO_M = FT_TO_M

REQUIRED_COLS = ["DIA", "HT", "TPA_UNADJ", "CCLCD"]

_HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_FG_DB = os.path.join(_HERE, "data", "functional_groups.json")
DEFAULT_SPCD_MAP = os.path.join(_HERE, "data", "spcd_func_group.csv")

# Columns written to the PlausibiliTree plot file, in order.
SST_COLUMNS = [
    "TREE_UID", "FUNC_GROUP", "X", "Y",
    "HT_LO", "HT_HI", "BASE_DIA_LO", "BASE_DIA_HI",
    "SPLIT_HT_LO", "SPLIT_HT_HI", "SPLIT_PROB",
    "SPLIT_NUM_LO", "SPLIT_NUM_HI",
    "NUM_BRANCH_LO", "NUM_BRANCH_HI",
    "MIN_CAN_HT_LO", "MIN_CAN_HT_HI",
    "MAX_CAN_W", "MAX_CAN_W_HT_LO", "MAX_CAN_W_HT_HI",
    "TOP_DIST", "MID_DIST", "FP_RAND_HT", "SPR_LO", "SPR_HI",
    "successful_placement",
]

# Plot-center pass-through columns: appended to SST_COLUMNS only when present
# in the source (see preprocess_fia). Unlike X/Y, these carry no offset.
PLOT_CENTER_COLUMNS = ["plot_lat", "plot_lon", "plot_elev"]


# ---------------------------------------------------------------------------
# Functional groups
# ---------------------------------------------------------------------------

# pre-W18 dictionary spellings -> the notes.md names (spec Sections 7.3,
# 8.6); legacy CSV/JSON files are accepted by renaming at load
LEGACY_PARAM_NAMES = {"CW_A0": "A0", "CW_A1": "A1", "TAPER_B": "TAPER"}


def normalize_param_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Rename legacy parameter-dictionary names to the notes.md names when
    only the legacy spelling is present (Section 8.6)."""
    ren = {old: new for old, new in LEGACY_PARAM_NAMES.items()
           if old in df.columns and new not in df.columns}
    return df.rename(columns=ren) if ren else df


def load_func_group_db(path: str = DEFAULT_FG_DB) -> pd.DataFrame:
    """Load the functional-group parameter table, indexed by FUNC_GROUP.

    Serialized as JSON per the spec (Section 8.6): one object keyed by
    FUNC_GROUP whose values are that code's parameters. Legacy CSV tables
    (one row per code) are still accepted; pre-W18 spellings rename at
    load (see normalize_param_columns)."""
    if str(path).endswith(".json"):
        with open(path) as fh:
            obj = json.load(fh)
        if "DEFAULT" not in obj:
            raise ValueError("functional group database must contain a DEFAULT key")
        df = pd.DataFrame.from_dict({k: dict(v) for k, v in obj.items()},
                                    orient="index")
        df.index.name = "FUNC_GROUP"
        return normalize_param_columns(df)
    db = normalize_param_columns(pd.read_csv(path).set_index("FUNC_GROUP"))
    if "DEFAULT" not in db.index:
        raise ValueError("functional group database must contain a DEFAULT row")
    return db


def load_spcd_map(path: str = DEFAULT_SPCD_MAP) -> dict:
    """Load the SPCD -> functional group lookup as a dict."""
    m = pd.read_csv(path)
    return dict(zip(m["SPCD"].astype(int), m["FUNC_GROUP"].astype(str)))


def assign_func_group(df: pd.DataFrame,
                      group_col: str = "SPCD",
                      spcd_map: dict | None = None,
                      fg_db: pd.DataFrame | None = None) -> pd.Series:
    """Return a FUNC_GROUP label for every row.

    Resolution order:
      1. explicit SPCD -> group lookup (species-level, most specific)
      2. SPGRPCD coarse split: 1-24 softwood (SW), 25-48 hardwood (HW)
      3. DEFAULT

    Any label that is not present in the functional-group database falls back
    to DEFAULT, so an incomplete lookup table is never fatal.
    """
    if spcd_map is None:
        spcd_map = load_spcd_map()
    if fg_db is None:
        fg_db = load_func_group_db()

    fg = pd.Series("DEFAULT", index=df.index, dtype=object)

    if "SPGRPCD" in df.columns:
        grp = pd.to_numeric(df["SPGRPCD"], errors="coerce")
        fg[grp.between(1, 24)] = "SW"
        fg[grp.between(25, 48)] = "HW"

    if group_col in df.columns and group_col == "SPCD":
        spcd = pd.to_numeric(df["SPCD"], errors="coerce")
        mapped = spcd.map(spcd_map)
        fg[mapped.notna()] = mapped[mapped.notna()]
    elif group_col not in ("SPCD", "SPGRPCD") and group_col in df.columns:
        # user-supplied functional group column
        fg = df[group_col].astype(str)

    known = set(fg_db.index)
    fg = fg.where(fg.isin(known), "DEFAULT")
    return fg


def _fg_param(fg: pd.Series, fg_db: pd.DataFrame, param: str) -> pd.Series:
    """Vectorised lookup of one functional-group parameter, with fallback."""
    lut = fg_db[param]
    out = fg.map(lut)
    return out.fillna(lut.loc["DEFAULT"]).astype(float)


# ---------------------------------------------------------------------------
# 1. Preprocessing
# ---------------------------------------------------------------------------

def plot_area_m2(shape="circle", *, radius_m=7.315, length_m=None, width_m=None,
                 side_m=None):
    """Plot area for the FIA per-acre -> count scaling."""
    if shape == "square":
        return float(side_m) ** 2
    if shape == "rectangle":
        return float(length_m) * float(width_m)
    return float(np.pi * radius_m ** 2)


def plot_bbox(shape="circle", *, radius_m=7.315, length_m=None, width_m=None,
             side_m=None):
    """Axis-aligned bounding box (x_min, x_max, y_min, y_max), origin-centred."""
    if shape == "square":
        h = float(side_m) / 2.0
        return (-h, h, -h, h)
    if shape == "rectangle":
        hx, hy = float(length_m) / 2.0, float(width_m) / 2.0
        return (-hx, hx, -hy, hy)
    return (-radius_m, radius_m, -radius_m, radius_m)


def _shape_attrs(shape, radius_m, length_m, width_m, side_m=None):
    """Normalised shape descriptor stored on df.attrs and passed downstream."""
    if shape == "square":
        return {"plot_shape": "square", "plot_side_m": float(side_m)}
    if shape == "rectangle":
        return {"plot_shape": "rectangle",
                "plot_length_m": float(length_m),
                "plot_width_m": float(width_m)}
    return {"plot_shape": "circle", "plot_radius_m": float(radius_m)}


def preprocess_fia(df: pd.DataFrame,
                   plot_radius_m: float = 7.315,
                   group_col: str = "SPCD",
                   ht_factor: float = HT_TO_M,
                   fg_db: pd.DataFrame | None = None,
                   spcd_map: dict | None = None,
                   plot_shape: str = "circle",
                   plot_length_m: float | None = None,
                   plot_width_m: float | None = None,
                   plot_side_m: float | None = None) -> pd.DataFrame:
    """Filter, unit-convert and annotate a raw FIA TREE table.

    - drops trees with DECAYCD >= 1 (snags are out of scope)
    - drops trees missing DIA, HT, TPA_UNADJ, CCLCD or the functional group col
    - DIA inches -> m, HT feet -> m
    - TPA_UNADJ trees/acre -> trees/m^2 -> expected trees in this plot
    - adds empty X, Y and successful_placement columns
    - if LAT/LON/ELEV are present in the input, passes them through unchanged
      as plot_lat/plot_lon/plot_elev (no unit conversion, no offset)

    plot_shape is 'circle' (use plot_radius_m), 'rectangle' (use plot_length_m
    x plot_width_m), or 'square' (use plot_side_m, a single side length). The
    plot area used to scale FIA per-acre density follows the chosen shape.
    """
    if fg_db is None:
        fg_db = load_func_group_db()

    out = df.copy()

    # snags out of scope; DECAYCD is only populated for standing dead trees
    if "DECAYCD" in out.columns:
        decay = pd.to_numeric(out["DECAYCD"], errors="coerce")
        out = out[~(decay >= 1)]

    needed = REQUIRED_COLS + [group_col]
    missing = [c for c in needed if c not in out.columns]
    if missing:
        raise KeyError(f"input is missing required column(s): {missing}")
    out = out.dropna(subset=needed)

    out["DIA"] = pd.to_numeric(out["DIA"], errors="coerce") * IN_TO_M
    out["HT"] = pd.to_numeric(out["HT"], errors="coerce") * ht_factor

    # Plot-center metadata, passed straight through unchanged (no unit
    # conversion, no offset) if present in the original input -- distinct
    # from X/Y below, which are local placement coordinates relative to plot
    # center. Only added when the source column is actually present, so an
    # input without FIA PLOT fields joined in is unaffected.
    for src_col, out_col in (("LAT", "plot_lat"), ("LON", "plot_lon"),
                             ("ELEV", "plot_elev")):
        if src_col in out.columns:
            out[out_col] = pd.to_numeric(out[src_col], errors="coerce")

    area = plot_area_m2(plot_shape, radius_m=plot_radius_m,
                        length_m=plot_length_m, width_m=plot_width_m,
                        side_m=plot_side_m)
    tpa = pd.to_numeric(out["TPA_UNADJ"], errors="coerce")
    out["TPA_M2"] = tpa / SQM_PER_ACRE          # trees per square metre
    out["TREES_IN_PLOT"] = out["TPA_M2"] * area

    out["FUNC_GROUP"] = assign_func_group(out, group_col=group_col,
                                          spcd_map=spcd_map, fg_db=fg_db)

    out["X"] = np.nan
    out["Y"] = np.nan
    out["successful_placement"] = pd.NA

    for k, v in _shape_attrs(plot_shape, plot_radius_m,
                             plot_length_m, plot_width_m,
                             side_m=plot_side_m).items():
        out.attrs[k] = v
    return out.reset_index(drop=True)


# ---------------------------------------------------------------------------
# 2. Variable plot sampling
# ---------------------------------------------------------------------------

def expand_by_tpa(df: pd.DataFrame, count_col: str = "TREES_IN_PLOT",
                  modifier: float = 1.0) -> pd.DataFrame:
    """Replicate each record to match its expansion factor.

    Per the design spec (Plausibilitree.md Section 6.1, step 5): for each tree
    the expected count is

        N = ceil(TPA [trees/m^2] * plot area [m^2] * modifier)

    with a floor of 1. Every FIA record is a tree that was physically
    tallied on the plot, so it places at least one tree regardless of its
    expansion factor; a modifier of 0 therefore generates exactly the
    documented trees (no expansion), values in (0, 1) thin globally, and
    values > 1 multiply the recorded stocking. When N > 1 the record is
    duplicated so the extra stems become independent trees.

    See Bechtold & Patterson (2005), GTR SRS-80, for the full FIA
    variable-radius / macroplot estimation design.
    """
    out = df.copy()
    raw = pd.to_numeric(out[count_col], errors="coerce").fillna(0.0) * modifier
    n = np.where(raw <= 0.0, 1, np.ceil(raw)).astype(int)
    n = np.maximum(n, 1)                    # every tallied tree is placed at least once
    out["N_REPS"] = n
    out = out.loc[out.index.repeat(n)].reset_index(drop=True)
    out["TREE_UID"] = [f"{i + 1:03d}" for i in range(len(out))]
    return out


# ---------------------------------------------------------------------------
# 3. Crown width allometry
# ---------------------------------------------------------------------------

def crown_width(dbh_m: pd.Series, func_group: pd.Series,
                fg_db: pd.DataFrame | None = None) -> pd.Series:
    """Maximum crown width (m) from CW = exp(a0 + a1 * ln(DBH)).

    a0/a1 are looked up per functional group.  DBH is passed in metres and
    converted to inches for the equation, which is the unit convention of the
    published FIA crown-width models; the result is converted back to metres.
    """
    if fg_db is None:
        fg_db = load_func_group_db()
    a0 = _fg_param(func_group, fg_db, "A0")
    a1 = _fg_param(func_group, fg_db, "A1")
    dbh_in = pd.to_numeric(dbh_m, errors="coerce") / IN_TO_M
    dbh_in = dbh_in.clip(lower=0.1)
    cw_ft = np.exp(a0 + a1 * np.log(dbh_in))
    return cw_ft * FT_TO_M


# ---------------------------------------------------------------------------
# 4. Tree placement
# ---------------------------------------------------------------------------

def _circle_intersection_area(d: float, r1: float, r2: float) -> float:
    """Area of intersection of two circles (radii r1, r2, centers d apart)."""
    if d >= r1 + r2:
        return 0.0
    if d <= abs(r1 - r2):
        return float(np.pi * min(r1, r2) ** 2)
    f1 = (d * d + r1 * r1 - r2 * r2) / (2.0 * d * r1)
    f2 = (d * d + r2 * r2 - r1 * r1) / (2.0 * d * r2)
    a1 = r1 * r1 * np.arccos(np.clip(f1, -1.0, 1.0))
    a2 = r2 * r2 * np.arccos(np.clip(f2, -1.0, 1.0))
    tri = 0.5 * np.sqrt(max((-d + r1 + r2) * (d + r1 - r2)
                            * (d - r1 + r2) * (d + r1 + r2), 0.0))
    return float(a1 + a2 - tri)


def _disc_outside_proportion(x: float, y: float, r: float, plot) -> float:
    """Proportion of a disc's area that lies outside the plot.

    plot is a shape dict from _shape_attrs (circle, rectangle or square).
    Analytic for circular plots; deterministic midpoint-grid quadrature
    (64 x 64) for rectangular/square plots.
    """
    if r <= 0.0:
        return 0.0
    if plot["plot_shape"] == "circle":
        R = plot["plot_radius_m"]
        d = float(np.hypot(x, y))
        if d + r <= R:
            return 0.0
        if d - r >= R:
            return 1.0
        disc = np.pi * r * r
        inter = _circle_intersection_area(d, R, r)
        return float(min(max((disc - inter) / disc, 0.0), 1.0))
    # rectangle / square
    if plot["plot_shape"] == "square":
        half = plot["plot_side_m"] / 2.0
        x_min, x_max, y_min, y_max = -half, half, -half, half
    else:
        hx, hy = plot["plot_length_m"] / 2.0, plot["plot_width_m"] / 2.0
        x_min, x_max, y_min, y_max = -hx, hx, -hy, hy
    n = 64
    step = 2.0 * r / n
    gx = x - r + (np.arange(n) + 0.5) * step
    gy = y - r + (np.arange(n) + 0.5) * step
    inside_disc = (gx[None, :] - x) ** 2 + (gy[:, None] - y) ** 2 <= r * r
    inside_plot = ((gx[None, :] >= x_min) & (gx[None, :] <= x_max)
                   & (gy[:, None] >= y_min) & (gy[:, None] <= y_max))
    total = int(inside_disc.sum())
    if total == 0:
        return 0.0
    return float((inside_disc & ~inside_plot).sum() / total)


def place_trees(df: pd.DataFrame,
                plot_radius_m: float | None = None,
                radius_col: str = "CROWN_W",
                overlap_basis: str = "stem",
                min_distance: float = 0.0,
                edge_effects: float = 0.0,
                max_tries: int = 100,
                seed: int | None = None,
                plot_shape: str | None = None,
                plot_length_m: float | None = None,
                plot_width_m: float | None = None,
                plot_side_m: float | None = None) -> pd.DataFrame:
    """Place trees inside the plot per the design spec (Section 6.1, step 7).

    Ordering: height descending, DBH as tiebreaker; exact ties broken by a
    seeded random draw.

    overlap_basis:
        'stem'  -- disc radius is half the stem diameter (the spec's trunk
                   disc; the default)
        'crown' -- disc radius is half the crown width (crown-cover packing)
    min_distance: required surface-to-surface separation between two trees'
    discs. Negative = discs may interpenetrate by up to |min_distance|;
    0 = discs cannot touch; positive = a buffer of at least min_distance.
    The protected disc of tree i therefore has radius r_i + min_distance / 2
    (for 'stem', diameter DBH + min_distance, exactly the spec's disc).

    edge_effects: 0.0-1.0, the allowed proportion of a tree's disc area
    outside the plot. 0.0 = the disc must fall entirely inside (the
    pre-existing behaviour); 1.0 = any position goes, except one that is
    100% outside, which cannot be drawn because positions are sampled
    inside the plot.

    A tree that fails max_tries attempts gets successful_placement = False
    and X/Y left as NaN. The per-tree counter resets after every tree.
    """
    if plot_radius_m is None:
        plot_radius_m = df.attrs.get("plot_radius_m", 7.315)
    if plot_shape is None:
        plot_shape = df.attrs.get("plot_shape", "circle")
    if plot_length_m is None:
        plot_length_m = df.attrs.get("plot_length_m")
    if plot_width_m is None:
        plot_width_m = df.attrs.get("plot_width_m")
    if plot_side_m is None:
        plot_side_m = df.attrs.get("plot_side_m")
    rng = np.random.default_rng(seed)
    plot = _shape_attrs(plot_shape, plot_radius_m,
                        plot_length_m, plot_width_m, side_m=plot_side_m)

    # Height descending, DBH tiebreaker, seeded random third key (spec)
    out = df.copy()
    out["_tiebreak"] = rng.random(len(out))
    out = (out.sort_values(["HT", "DIA", "_tiebreak"],
                           ascending=[False, False, True])
           .drop(columns="_tiebreak").reset_index(drop=True))

    if overlap_basis == "stem":
        radii = pd.to_numeric(out["DIA"], errors="coerce").fillna(0.0) / 2.0
    else:
        radii = pd.to_numeric(out[radius_col], errors="coerce").fillna(0.0) / 2.0

    # protected discs: surface-to-surface separation of min_distance
    disc_r = radii.to_numpy() + min_distance / 2.0

    placed_x, placed_y, placed_r = [], [], []
    xs = np.full(len(out), np.nan)
    ys = np.full(len(out), np.nan)
    ok = np.zeros(len(out), dtype=bool)

    is_rect = plot_shape in ("rectangle", "square")
    if plot_shape == "square":
        hx = hy = float(plot_side_m) / 2.0
    elif is_rect:
        hx, hy = float(plot_length_m) / 2.0, float(plot_width_m) / 2.0
    else:
        plot_r = float(plot_radius_m)

    for i in range(len(out)):
        r_i = disc_r[i]
        quota_r = max(r_i, 0.0)       # disc whose area the edge quota applies to
        for _attempt in range(max(1, int(max_tries))):
            if is_rect:
                x = rng.uniform(-hx, hx)
                y = rng.uniform(-hy, hy)
            else:
                # uniform over the full plot disc (positions are always
                # inside the plot; the edge quota decides how far out)
                rad = plot_r * np.sqrt(rng.random())
                th = 2 * np.pi * rng.random()
                x, y = rad * np.cos(th), rad * np.sin(th)
            # edge quota: allowed proportion of the disc outside the plot
            if _disc_outside_proportion(x, y, quota_r, plot) > edge_effects + 1e-9:
                continue
            # disc/disc separation against every previously placed tree
            if placed_x:
                dx = np.asarray(placed_x) - x
                dy = np.asarray(placed_y) - y
                clash = np.any(np.hypot(dx, dy)
                               < (np.maximum(np.asarray(placed_r), 0.0)
                                  + max(r_i, 0.0)))
            else:
                clash = False
            if clash:
                continue
            xs[i], ys[i], ok[i] = x, y, True
            placed_x.append(x)
            placed_y.append(y)
            placed_r.append(r_i)
            break

    out["X"] = xs
    out["Y"] = ys
    out["successful_placement"] = ok
    for k, v in plot.items():
        out.attrs[k] = v
    return out


# ---------------------------------------------------------------------------
# 5. Taper: DBH -> base diameter
# ---------------------------------------------------------------------------

BREAST_HEIGHT_M = 1.37


def dbh_to_base_diameter(dbh_m: pd.Series, ht_m: pd.Series,
                         func_group: pd.Series,
                         fg_db: pd.DataFrame | None = None) -> pd.Series:
    """Extrapolate DBH down to a ground-line (base) diameter.

    Simple single-parameter taper of the Kozak/Ormerod form

        d(h) = DBH * ((H - h) / (H - 1.37)) ** (1 / b)

    evaluated at h = 0, where b is a functional-group shape exponent
    (larger b = less butt swell).  Trees shorter than breast height fall back
    to DBH.
    """
    if fg_db is None:
        fg_db = load_func_group_db()
    b = _fg_param(func_group, fg_db, "TAPER")
    dbh = pd.to_numeric(dbh_m, errors="coerce")
    ht = pd.to_numeric(ht_m, errors="coerce")
    ratio = (ht / (ht - BREAST_HEIGHT_M)).where(ht > BREAST_HEIGHT_M + 0.1, 1.0)
    return dbh * ratio ** (1.0 / b)


# ---------------------------------------------------------------------------
# 6. Postprocessing
# ---------------------------------------------------------------------------

def build_sst_table(df: pd.DataFrame,
                    fg_db: pd.DataFrame | None = None,
                    dia_tol: float = 0.01,
                    ht_tol: float = 0.01,
                    drop_failed: bool = False) -> pd.DataFrame:
    """Reduce the working frame to the PlausibiliTree parameter table.

    Unplaced trees are retained by default with successful_placement=False so
    that the downstream generator can filter them (and so the plot table stays
    a complete record of the FIA input); set drop_failed=True to omit them.

    If plot_lat/plot_lon/plot_elev were added by preprocess_fia (i.e. LAT/LON/
    ELEV were present in the original input), they are carried through as
    extra trailing columns.
    """
    if fg_db is None:
        fg_db = load_func_group_db()

    work = df
    if drop_failed and "successful_placement" in work.columns:
        work = work[work["successful_placement"] == True]  # noqa: E712
    work = work.reset_index(drop=True)

    fg = work["FUNC_GROUP"].astype(str)
    base_dia = dbh_to_base_diameter(work["DIA"], work["HT"], fg, fg_db)
    ht = pd.to_numeric(work["HT"], errors="coerce")

    out = pd.DataFrame(index=work.index)
    out["TREE_UID"] = (work["TREE_UID"] if "TREE_UID" in work.columns
                       else [f"{i + 1:03d}" for i in range(len(work))])
    out["FUNC_GROUP"] = fg
    out["X"] = work["X"].round(4)
    out["Y"] = work["Y"].round(4)
    out["HT_LO"] = (ht * (1 - ht_tol)).round(4)
    out["HT_HI"] = (ht * (1 + ht_tol)).round(4)
    out["BASE_DIA_LO"] = (base_dia * (1 - dia_tol)).round(4)
    out["BASE_DIA_HI"] = (base_dia * (1 + dia_tol)).round(4)
    out["MAX_CAN_W"] = pd.to_numeric(work["CROWN_W"], errors="coerce").round(4)
    out["successful_placement"] = (work["successful_placement"]
                                   if "successful_placement" in work.columns else True)

    for p in ["SPLIT_HT_LO", "SPLIT_HT_HI", "SPLIT_PROB",
              "NUM_BRANCH_LO", "NUM_BRANCH_HI",
              "MIN_CAN_HT_LO", "MIN_CAN_HT_HI",
              "MAX_CAN_W_HT_LO", "MAX_CAN_W_HT_HI",
              "TOP_DIST", "MID_DIST", "FP_RAND_HT", "SPR_LO", "SPR_HI"]:
        out[p] = _fg_param(fg, fg_db, p)

    # SPLIT_NUM_* is an optional column pair: older functional-group tables
    # don't carry it, and its absence means the original single fork (1:1)
    if "SPLIT_NUM_LO" in fg_db.columns and "SPLIT_NUM_HI" in fg_db.columns:
        out["SPLIT_NUM_LO"] = _fg_param(fg, fg_db, "SPLIT_NUM_LO").astype(int)
        out["SPLIT_NUM_HI"] = _fg_param(fg, fg_db, "SPLIT_NUM_HI").astype(int)
    else:
        out["SPLIT_NUM_LO"] = 1
        out["SPLIT_NUM_HI"] = 1

    out["NUM_BRANCH_LO"] = out["NUM_BRANCH_LO"].astype(int)
    out["NUM_BRANCH_HI"] = out["NUM_BRANCH_HI"].astype(int)

    present = [c for c in PLOT_CENTER_COLUMNS if c in work.columns]
    for c in present:
        out[c] = work[c]
    return out[SST_COLUMNS + present]


# ---------------------------------------------------------------------------
# Pipeline
# ---------------------------------------------------------------------------

def fia_to_sst_plot(fia_df: pd.DataFrame,
                    plot_radius_m: float = 7.315,
                    group_col: str = "SPCD",
                    overlap_basis: str = "stem",
                    expansion_modifier: float = 1.0,
                    min_distance: float = 0.0,
                    edge_effects: float = 0.0,
                    max_tries: int = 100,
                    ht_factor: float = HT_TO_M,
                    dia_tol: float = 0.01,
                    ht_tol: float = 0.01,
                    seed: int | None = None,
                    out_csv: str | None = None,
                    plot_shape: str = "circle",
                    plot_length_m: float | None = None,
                    plot_width_m: float | None = None,
                    plot_side_m: float | None = None,
                    drop_failed: bool = True,
                    return_intermediate: bool = False):
    """Run the full FIA -> PlausibiliTree conversion for one plot.

    plot_shape selects 'circle' (plot_radius_m), 'rectangle' (plot_length_m x
    plot_width_m), or 'square' (plot_side_m, a single side length).

    By default only successfully placed trees are kept in the output table
    and TREE_UIDs are reassigned sequentially over the placed list (design
    Section 6.1, step 7); set drop_failed=False to keep unplaced trees with
    successful_placement = False (their UIDs are then the expansion-order
    ones).

    Returns the SST table, or (sst_table, working_frame) if
    return_intermediate is True.
    """
    fg_db = load_func_group_db()

    work = preprocess_fia(fia_df, plot_radius_m=plot_radius_m,
                          group_col=group_col, ht_factor=ht_factor, fg_db=fg_db,
                          plot_shape=plot_shape, plot_length_m=plot_length_m,
                          plot_width_m=plot_width_m, plot_side_m=plot_side_m)
    work = expand_by_tpa(work, modifier=expansion_modifier)
    work["CROWN_W"] = crown_width(work["DIA"], work["FUNC_GROUP"], fg_db)
    work = place_trees(work, plot_radius_m=plot_radius_m,
                       overlap_basis=overlap_basis,
                       min_distance=min_distance,
                       edge_effects=edge_effects,
                       max_tries=max_tries, seed=seed,
                       plot_shape=plot_shape, plot_length_m=plot_length_m,
                       plot_width_m=plot_width_m, plot_side_m=plot_side_m)
    sst = build_sst_table(work, fg_db=fg_db, dia_tol=dia_tol, ht_tol=ht_tol,
                          drop_failed=drop_failed)
    if drop_failed:
        sst = sst.reset_index(drop=True)
        sst["TREE_UID"] = [f"{i + 1:03d}" for i in range(len(sst))]

    if out_csv:
        sst.to_csv(out_csv, index=False)
    return (sst, work) if return_intermediate else sst


def placement_report(params: dict, work: pd.DataFrame) -> dict:
    """Report dict for one processed plot (design Section 6.1, step 8)."""
    placed = pd.Series(work["successful_placement"], dtype="boolean")
    n_placed = int(placed.fillna(False).sum())
    n_total = int(len(work))
    return {
        "parameters": params,
        "trees_expanded": n_total,
        "trees_placed": n_placed,
        "trees_not_placed": n_total - n_placed,
        "proportion_placed": (n_placed / n_total) if n_total else 1.0,
    }


def build_tree_list(work: pd.DataFrame) -> dict:
    """Placed trees (working frame) -> canonical tree list dict (Section 8.1).

    Inventory facts only: TREE_UID, X, Y, DBH, HT, TREE_CODE (+Z). The
    plot-center passthrough columns become the tree list's plot_center.
    Simulator parameters are NOT included -- generate_trees looks them up
    per TREE_CODE at generation time.
    """
    p = work[work["successful_placement"] == True]  # noqa: E712
    p = p.reset_index(drop=True)
    trees = []
    for i, row in p.iterrows():
        # DBH/HT carry full inventory precision: generation-time allometry
        # (crown width, taper) is computed FROM these values (Section 7.2),
        # so rounding here would shift every derived parameter.
        t = {"TREE_UID": f"{i + 1:03d}",
             "X": round(float(row["X"]), 4),
             "Y": round(float(row["Y"]), 4),
             "DBH": float(row["DIA"]),
             "HT": float(row["HT"]),
             "TREE_CODE": str(row["FUNC_GROUP"]),
             "Z": 0.0}
        trees.append(t)
    center = {}
    for col, key in (("plot_lat", "lat"), ("plot_lon", "lon"),
                     ("plot_elev", "elev")):
        if col in work.columns:
            vals = pd.to_numeric(work[col], errors="coerce").dropna()
            center[key] = (float(vals.iloc[0]) if len(vals) else None)
        else:
            center[key] = None
    return {
        "plot_uid": str(work.attrs.get("plot_uid", "")),
        "plot_shape": work.attrs.get("plot_shape", "circle"),
        "plot_radius_m": work.attrs.get("plot_radius_m"),
        "plot_length_m": work.attrs.get("plot_length_m"),
        "plot_width_m": work.attrs.get("plot_width_m"),
        "plot_side_m": work.attrs.get("plot_side_m"),
        "plot_center": center,
        "trees": trees,
    }


def fia_to_tree_list(fia_df: pd.DataFrame, *, return_intermediate: bool = False,
                     **kwargs) -> dict:
    """FIA TREE records -> one canonical tree list dict (Section 8.1).

    Runs the same placement pipeline as fia_to_sst_plot (placement lives in
    this helper, never the core), then emits the minimal inventory-fact
    contract instead of frozen SST parameters. kwargs are the placement
    controls of fia_to_sst_plot (expansion_modifier, min_distance,
    edge_effects, max_tries, seed, plot shape, plot_uid ...). Returns
    (tree_list, work_frame) if return_intermediate, else the tree list.
    """
    kw = dict(kwargs)
    kw.pop("out_csv", None)
    plot_uid = str(kw.pop("plot_uid", ""))
    kw["return_intermediate"] = True
    sst_unused, work = fia_to_sst_plot(fia_df, **kw)
    work = work.copy()
    work.attrs["plot_uid"] = plot_uid
    plot = build_tree_list(work)
    return (plot, work) if return_intermediate else plot


def split_plots(fia_df: pd.DataFrame, plot_col: str = "PLT_CN"):
    """Yield (plot_id, subframe) for every plot in a multi-plot FIA table."""
    for pid, sub in fia_df.groupby(plot_col):
        yield pid, sub.reset_index(drop=True)


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser(description="FIA -> PlausibiliTree plot tables")
    ap.add_argument("csv", help="FIA TREE csv")
    ap.add_argument("-o", "--outdir", default=".")
    ap.add_argument("-r", "--radius", type=float, default=None,
                    help="plot radius in metres, for --plot-shape circle "
                         "(default 7.315 = FIA subplot; mutually exclusive "
                         "with --plot-shape rectangle/square)")
    ap.add_argument("--plot-shape", choices=["circle", "rectangle", "square"],
                    default="circle",
                    help="plot geometry (default circle)")
    ap.add_argument("--length", type=float, default=None,
                    help="rectangle length in metres (x), for --plot-shape rectangle")
    ap.add_argument("--width", type=float, default=None,
                    help="rectangle width in metres (y), for --plot-shape rectangle")
    ap.add_argument("--side", type=float, default=None,
                    help="side length in metres, for --plot-shape square")
    ap.add_argument("--plot-col", default="PLT_CN")
    ap.add_argument("--group-col", default="SPCD")
    ap.add_argument("--basis", default="stem", choices=["stem", "crown"],
                    help="disc used for spacing: trunk diameter (spec default) "
                         "or crown width")
    ap.add_argument("--expansion-modifier", type=float,
                    dest="expansion_modifier", default=1.0,
                    help="multiplied into FIA-recorded TPA before rounding up. "
                         "0 = only the documented trees; (0,1) = global "
                         "reduction; >1 = more trees than the recorded TPA")
    ap.add_argument("--min-distance", type=float, dest="min_distance",
                    default=0.0,
                    help="required surface-to-surface separation between "
                         "trees' discs, metres. Negative = discs may overlap "
                         "by up to |d|; 0 = trees cannot touch; positive = "
                         "required buffer")
    ap.add_argument("--max-tries", type=int, dest="max_tries", default=100,
                    help="placement attempts per tree before it is dropped")
    ap.add_argument("--edge-effects", type=float, dest="edge_effects",
                    default=0.0,
                    help="allowed proportion of a tree's disc area outside "
                         "the plot, 0.0 (fully inside) to 1.0")
    ap.add_argument("--format", choices=["tree-list-csv", "tree-list-json",
                                        "sst-csv"],
                    default="tree-list-csv",
                    help="output contract: the minimal tree-list CSV "
                         "(design Section 8.1, notes.md default -- one row "
                         "per tree; plot geometry lives in the report JSON) "
                         "or the alternate canonical JSON; sst-csv emits the "
                         "legacy SST parameter table")
    ap.add_argument("--seed", type=int, default=None)
    ap.add_argument("--limit", type=int, default=None, help="max plots to write")
    args = ap.parse_args()

    # --radius only means anything for a circular plot (it sets the disc
    # radius used for both placement and the FIA per-acre -> count scaling);
    # for rectangle/square that scaling instead comes from --length/--width
    # or --side, so mixing in --radius as well is almost certainly a mistake.
    if args.plot_shape == "circle":
        if args.length is not None or args.width is not None or args.side is not None:
            ap.error("--length/--width/--side are not allowed with "
                     "--plot-shape circle")
        radius = 7.315 if args.radius is None else args.radius
    else:
        if args.radius is not None:
            ap.error(f"--radius is not allowed with --plot-shape "
                     f"{args.plot_shape} (mutually exclusive; use "
                     f"{'--length/--width' if args.plot_shape == 'rectangle' else '--side'} "
                     f"instead)")
        radius = 7.315  # unused for rectangle/square, kept only as a fallback
        if args.plot_shape == "rectangle":
            if args.length is None or args.width is None:
                ap.error("--plot-shape rectangle requires --length and --width")
            if args.side is not None:
                ap.error("--side is not allowed with --plot-shape rectangle")
        else:  # square
            if args.side is None:
                ap.error("--plot-shape square requires --side")
            if args.length is not None or args.width is not None:
                ap.error("--length/--width are not allowed with --plot-shape square")

    df = pd.read_csv(args.csv, low_memory=False)
    os.makedirs(args.outdir, exist_ok=True)
    import json
    for i, (pid, sub) in enumerate(split_plots(df, args.plot_col)):
        if args.limit is not None and i >= args.limit:
            break
        params = {"plot_uid": str(pid),
                  "plot_shape": args.plot_shape,
                  "plot_radius_m": radius if args.plot_shape == "circle" else None,
                  "plot_length_m": args.length,
                  "plot_width_m": args.width,
                  "plot_side_m": args.side,
                  "group_col": args.group_col,
                  "overlap_basis": args.basis,
                  "expansion_modifier": args.expansion_modifier,
                  "min_distance_m": args.min_distance,
                  "max_tries": args.max_tries,
                  "edge_effects": args.edge_effects,
                  "seed": args.seed,
                  "format": args.format}
        if args.format in ("tree-list-csv", "tree-list-json"):
            # canonical Section 8.1 tree list + report JSON; the CSV form
            # carries tree rows only, so the report doubles as the
            # plot-geometry sidecar (tree_list reads it back)
            if args.format == "tree-list-csv":
                path = os.path.join(args.outdir, f"plot_{pid}.csv")
            else:
                path = os.path.join(args.outdir, f"plot_{pid}.json")
            plot, work = fia_to_tree_list(sub, plot_radius_m=radius,
                                           group_col=args.group_col,
                                           overlap_basis=args.basis,
                                           expansion_modifier=args.expansion_modifier,
                                           min_distance=args.min_distance,
                                           edge_effects=args.edge_effects,
                                           max_tries=args.max_tries,
                                           seed=args.seed,
                                           plot_shape=args.plot_shape,
                                           plot_length_m=args.length,
                                           plot_width_m=args.width,
                                           plot_side_m=args.side,
                                           plot_uid=f"plot_{pid}",
                                           return_intermediate=True)
            if args.format == "tree-list-csv":
                import tree_list
                tree_list.save_plot_csv(plot, path)
            else:
                with open(path, "w") as fh:
                    json.dump(plot, fh, indent=2)
            report = placement_report(params, work)
            report["plot_center"] = plot["plot_center"]
            report_path = os.path.join(args.outdir, f"plot_{pid}_report.json")
            with open(report_path, "w") as fh:
                json.dump(report, fh, indent=2)
            print(f"{path}: {report['trees_placed']} trees placed "
                  f"({report['trees_not_placed']} not placed, "
                  f"{report['proportion_placed']:.3f}); report: {report_path}")
            continue
        path = os.path.join(args.outdir, f"plot_{pid}.csv")
        sst, work = fia_to_sst_plot(sub, plot_radius_m=radius,
                                    group_col=args.group_col,
                                    overlap_basis=args.basis,
                                    expansion_modifier=args.expansion_modifier,
                                    min_distance=args.min_distance,
                                    edge_effects=args.edge_effects,
                                    max_tries=args.max_tries,
                                    seed=args.seed, out_csv=path,
                                    plot_shape=args.plot_shape,
                                    plot_length_m=args.length,
                                    plot_width_m=args.width,
                                    plot_side_m=args.side,
                                    return_intermediate=True)
        report = placement_report(params, work)
        report_path = os.path.join(args.outdir, f"plot_{pid}_report.json")
        with open(report_path, "w") as fh:
            json.dump(report, fh, indent=2)
        print(f"{path}: {report['trees_placed']} trees placed "
              f"({report['trees_not_placed']} not placed, "
              f"{report['proportion_placed']:.3f}); report: {report_path}")
