"""
fia_helper.py: Convert FIA (Forest Inventory and Analysis) tree records into a
per-plot parameter table that PlausibiliTree can consume directly.

Pipeline (see fia_to_sst_plot):

    preprocess -> expand by TPA -> crown width -> placement -> postprocess

Each stage is also usable on its own so that individual steps can be swapped
out (in particular the TPA expansion, which is a placeholder -- see
expand_by_tpa).

If the input carries FIA PLOT fields LAT/LON/ELEV, they are passed through
unchanged as plot_lat/plot_lon/plot_elev columns on the output (plot-center
values, no offset) -- see preprocess_fia and PLOT_CENTER_COLUMNS.

Requires: numpy, pandas
Companion data files (same directory by default):
    data/functional_groups.csv   parameters keyed by functional group
    data/spcd_func_group.csv     SPCD -> functional group lookup
"""

from __future__ import annotations

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
DEFAULT_FG_DB = os.path.join(_HERE, "data", "functional_groups.csv")
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

def load_func_group_db(path: str = DEFAULT_FG_DB) -> pd.DataFrame:
    """Load the functional-group parameter table, indexed by FUNC_GROUP."""
    db = pd.read_csv(path).set_index("FUNC_GROUP")
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

def plot_area_m2(shape="circle", *, radius_m=7.32, length_m=None, width_m=None,
                 side_m=None):
    """Plot area for the FIA per-acre -> count scaling."""
    if shape == "square":
        return float(side_m) ** 2
    if shape == "rectangle":
        return float(length_m) * float(width_m)
    return float(np.pi * radius_m ** 2)


def plot_bbox(shape="circle", *, radius_m=7.32, length_m=None, width_m=None,
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
                   plot_radius_m: float = 7.32,
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
# 2. Variable plot sampling  (PLACEHOLDER)
# ---------------------------------------------------------------------------

def expand_by_tpa(df: pd.DataFrame, count_col: str = "TREES_IN_PLOT",
                  rng: np.random.Generator | None = None) -> pd.DataFrame:
    """Replicate each record to match its expansion factor.

    Each FIA record is a tree that was physically tallied on the plot, so every
    record places AT LEAST ONE tree regardless of its expansion factor: the
    count is rounded to the nearest integer but floored at 1, never 0. (An
    expansion factor below 1 scales the per-area estimate down, but it can never
    mean fewer trees than were actually observed in the field.) When the count
    rounds above 1 the record is duplicated to add the extra stems.

    Rounding is unbiased only in aggregate; with subplot-level TPA values near
    1.0 it collapses most of the between-tree variation, but the one observed
    tree is always placed. See Bechtold & Patterson (2005), GTR SRS-80, for the
    full FIA variable-radius / macroplot estimation design.
    """
    out = df.copy()
    raw = pd.to_numeric(out[count_col], errors="coerce").fillna(0)
    n = np.rint(raw).astype(int)
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
    a0 = _fg_param(func_group, fg_db, "CW_A0")
    a1 = _fg_param(func_group, fg_db, "CW_A1")
    dbh_in = pd.to_numeric(dbh_m, errors="coerce") / IN_TO_M
    dbh_in = dbh_in.clip(lower=0.1)
    cw_ft = np.exp(a0 + a1 * np.log(dbh_in))
    return cw_ft * FT_TO_M


# ---------------------------------------------------------------------------
# 4. Tree placement
# ---------------------------------------------------------------------------

def place_trees(df: pd.DataFrame,
                plot_radius_m: float | None = None,
                radius_col: str = "CROWN_W",
                overlap_basis: str = "crown",
                allowed_overlap: float = 0.0,
                max_tries: int = 100,
                seed: int | None = None,
                plot_shape: str | None = None,
                plot_length_m: float | None = None,
                plot_width_m: float | None = None,
                plot_side_m: float | None = None) -> pd.DataFrame:
    """Place non-overlapping circles inside the plot, largest/most dominant first.

    Ordering: CCLCD ascending (1 = open grown / dominant first), then DBH
    descending, then HT descending.

    overlap_basis:
        'crown' -- circle radius is half the crown width (crown-cover packing)
        'stem'  -- circle radius is half the stem diameter (stems only)
    allowed_overlap is the fraction of the circle radius that neighbouring
    circles are permitted to interpenetrate (0 = strictly tangent or apart).

    A tree that fails max_tries attempts gets successful_placement = False and
    X/Y left as NaN.  The per-tree counter resets after every tree, as
    specified.
    """
    if plot_radius_m is None:
        plot_radius_m = df.attrs.get("plot_radius_m", 7.32)
    if plot_shape is None:
        plot_shape = df.attrs.get("plot_shape", "circle")
    if plot_length_m is None:
        plot_length_m = df.attrs.get("plot_length_m")
    if plot_width_m is None:
        plot_width_m = df.attrs.get("plot_width_m")
    if plot_side_m is None:
        plot_side_m = df.attrs.get("plot_side_m")
    rng = np.random.default_rng(seed)

    out = df.sort_values(["CCLCD", "DIA", "HT"],
                         ascending=[True, False, False]).reset_index(drop=True)

    if overlap_basis == "stem":
        radii = pd.to_numeric(out["DIA"], errors="coerce").fillna(0.0) / 2.0
    else:
        radii = pd.to_numeric(out[radius_col], errors="coerce").fillna(0.0) / 2.0
    radii = radii.to_numpy() * (1.0 - allowed_overlap)

    placed_x, placed_y, placed_r = [], [], []
    xs = np.full(len(out), np.nan)
    ys = np.full(len(out), np.nan)
    ok = np.zeros(len(out), dtype=bool)

    is_rect = plot_shape in ("rectangle", "square")
    if plot_shape == "square":
        hx = hy = float(plot_side_m) / 2.0
    elif is_rect:
        hx, hy = float(plot_length_m) / 2.0, float(plot_width_m) / 2.0

    for i in range(len(out)):
        r_i = radii[i]
        # stem must fall inside the boundary, so the centre is inset by r_i
        if is_rect:
            fx = max(hx - r_i, 0.0)
            fy = max(hy - r_i, 0.0)
        else:
            r_free = max(plot_radius_m - r_i, 0.0)
        overlap_counter = 0
        while overlap_counter <= max_tries:
            if is_rect:
                # uniform over the inset rectangle (or inset square, same math)
                x = rng.uniform(-fx, fx)
                y = rng.uniform(-fy, fy)
            else:
                # uniform over the disc
                rad = r_free * np.sqrt(rng.random())
                th = 2 * np.pi * rng.random()
                x, y = rad * np.cos(th), rad * np.sin(th)
            if placed_x:
                dx = np.asarray(placed_x) - x
                dy = np.asarray(placed_y) - y
                clash = np.any(np.hypot(dx, dy) < (np.asarray(placed_r) + r_i))
            else:
                clash = False
            if clash:
                overlap_counter += 1
                continue
            xs[i], ys[i], ok[i] = x, y, True
            placed_x.append(x)
            placed_y.append(y)
            placed_r.append(r_i)
            break

    out["X"] = xs
    out["Y"] = ys
    out["successful_placement"] = ok
    for k, v in _shape_attrs(plot_shape, plot_radius_m,
                             plot_length_m, plot_width_m,
                             side_m=plot_side_m).items():
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
    b = _fg_param(func_group, fg_db, "TAPER_B")
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
                    plot_radius_m: float = 7.32,
                    group_col: str = "SPCD",
                    overlap_basis: str = "crown",
                    allowed_overlap: float = 0.0,
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
                    return_intermediate: bool = False):
    """Run the full FIA -> PlausibiliTree conversion for one plot.

    plot_shape selects 'circle' (plot_radius_m), 'rectangle' (plot_length_m x
    plot_width_m), or 'square' (plot_side_m, a single side length). Returns
    the SST table, or (sst_table, working_frame) if return_intermediate is
    True.
    """
    fg_db = load_func_group_db()

    work = preprocess_fia(fia_df, plot_radius_m=plot_radius_m,
                          group_col=group_col, ht_factor=ht_factor, fg_db=fg_db,
                          plot_shape=plot_shape, plot_length_m=plot_length_m,
                          plot_width_m=plot_width_m, plot_side_m=plot_side_m)
    work = expand_by_tpa(work)
    work["CROWN_W"] = crown_width(work["DIA"], work["FUNC_GROUP"], fg_db)
    work = place_trees(work, plot_radius_m=plot_radius_m,
                       overlap_basis=overlap_basis,
                       allowed_overlap=allowed_overlap,
                       max_tries=max_tries, seed=seed,
                       plot_shape=plot_shape, plot_length_m=plot_length_m,
                       plot_width_m=plot_width_m, plot_side_m=plot_side_m)
    sst = build_sst_table(work, fg_db=fg_db, dia_tol=dia_tol, ht_tol=ht_tol)

    if out_csv:
        sst.to_csv(out_csv, index=False)
    return (sst, work) if return_intermediate else sst


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
                         "(default 7.32 = FIA subplot; mutually exclusive "
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
    ap.add_argument("--basis", default="crown", choices=["crown", "stem"])
    ap.add_argument("--allowed-overlap", type=float, default=0.0)
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
        radius = 7.32 if args.radius is None else args.radius
    else:
        if args.radius is not None:
            ap.error(f"--radius is not allowed with --plot-shape "
                     f"{args.plot_shape} (mutually exclusive; use "
                     f"{'--length/--width' if args.plot_shape == 'rectangle' else '--side'} "
                     f"instead)")
        radius = 7.32  # unused for rectangle/square, kept only as a fallback
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
    for i, (pid, sub) in enumerate(split_plots(df, args.plot_col)):
        if args.limit is not None and i >= args.limit:
            break
        path = os.path.join(args.outdir, f"plot_{pid}.csv")
        sst = fia_to_sst_plot(sub, plot_radius_m=radius,
                              group_col=args.group_col,
                              overlap_basis=args.basis,
                              allowed_overlap=args.allowed_overlap,
                              seed=args.seed, out_csv=path,
                              plot_shape=args.plot_shape,
                              plot_length_m=args.length,
                              plot_width_m=args.width,
                              plot_side_m=args.side)
        print(f"{path}: {len(sst)} trees")
