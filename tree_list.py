"""
tree_list.py: the canonical JSON tree list (design spec Section 8.1) and the
one place in the codebase that understands both its serialization and the
legacy SST parameter CSV.

Canonical form (one JSON object per plot):

    {
      "plot_uid": "11839",
      "plot_shape": "circle",           # circle | rectangle | square
      "plot_radius_m": 20.0,            # only for circle
      "plot_length_m": 30.0,            # only for rectangle
      "plot_width_m": 30.0,             # only for rectangle
      "plot_side_m": 30.0,              # only for square
      "plot_center": {"lon": null, "lat": null, "elev": null },
      "trees": [
        {"TREE_UID": "001", "X": -3.1, "Y": 8.4, "DBH": 0.45,
         "HT": 44.6, "TREE_CODE": "RW", "Z": 0.0}
      ]
    }

Required tree attributes: TREE_UID, X, Y, DBH, HT, TREE_CODE.
Optional tree attributes: Z, LON, LAT, ELEV.

The tree list carries inventory facts only -- simulator parameters never
appear in it (design principle 2); generate_trees.py looks them up per
tree code at generation time.

Legacy SST parameter CSVs (TREE_UID, FUNC_GROUP, X, Y, HT_LO/HI,
BASE_DIA_LO/HI, ...) load through `legacy_trees_frame`, which supplies
everything downstream stages read from a table (TREE_UID, X, Y, the
plot-center passthrough columns) but cannot recover DBH/HT/TREE_CODE --
those lists precede the minimal contract. Generation from legacy CSVs is
still supported by generate_trees.py reading them directly.

Requires: pandas
"""

from __future__ import annotations

import json
import os

import pandas as pd

REQUIRED_TREE_ATTRS = ["TREE_UID", "X", "Y", "DBH", "HT", "TREE_CODE"]
OPTIONAL_TREE_ATTRS = ["Z", "LON", "LAT", "ELEV"]

# geo columns downstream stages have always read (assemble_plot's
# plot-center passthrough); supplied for JSON lists from plot_center.
GEO_COLS = ["plot_lat", "plot_lon", "plot_elev"]


# ---------------------------------------------------------------------------
# Canonical JSON
# ---------------------------------------------------------------------------

def save_plot(plot: dict, path: str) -> str:
    """Write one canonical tree-list JSON. `plot["trees"]` may be a DataFrame
    or a list of dicts. Returns the path written."""
    if hasattr(plot.get("trees"), "to_dict"):
        records = plot["trees"].to_dict(orient="records")
    else:
        records = list(plot.get("trees", []))
    obj = {k: v for k, v in plot.items() if k != "trees"}
    obj["trees"] = [{k: v for k, v in t.items() if v is not None or k in
                     REQUIRED_TREE_ATTRS} for t in records]
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=2)
    return path


def validate_plot(plot: dict) -> None:
    """Raise ValueError if the tree list violates the Section 8.1 contract."""
    trees = plot.get("trees")
    if trees is None:
        raise ValueError("tree list has no 'trees' array")
    if hasattr(trees, "iterrows"):
        records = trees.to_dict(orient="records")
    else:
        records = list(trees)
    if not records:
        raise ValueError("tree list 'trees' array is empty")
    for i, t in enumerate(records):
        missing = [a for a in REQUIRED_TREE_ATTRS
                   if t.get(a) is None or (isinstance(t.get(a), float)
                                           and pd.isna(t.get(a)))]
        if missing:
            raise ValueError(f"tree {i} is missing required attribute(s) "
                             f"{missing}")
        if float(t["DBH"]) <= 0 or not float(t["HT"]) > 0:
            raise ValueError(f"tree {t['TREE_UID']!r} has non-physical "
                             f"DBH={t['DBH']} / HT={t['HT']}")
    shape = plot.get("plot_shape", "circle")
    if shape not in ("circle", "rectangle", "square", "CIRCLE", "RECTANGLE"):
        raise ValueError(f"unknown plot_shape {shape!r}")


def load_plot(path: str) -> dict:
    """Load a canonical tree-list JSON -> plot dict with `trees` a DataFrame.

    plot_uid defaults from the filename when absent. Unknown/extra keys are
    preserved. plot_center defaults to all-null."""
    with open(path) as fh:
        obj = json.load(fh)
    if "trees" not in obj:
        raise ValueError(f"{path}: not a canonical tree list (no 'trees')")
    uid = obj.get("plot_uid") or os.path.splitext(os.path.basename(path))[0]
    df = pd.DataFrame(obj["trees"])
    if "TREE_UID" in df.columns:
        df["TREE_UID"] = df["TREE_UID"].astype(str)
    plot = dict(obj)
    plot["plot_uid"] = str(uid)
    plot["trees"] = df
    center = obj.get("plot_center") or {}
    plot["plot_center"] = {"lon": center.get("lon"),
                           "lat": center.get("lat"),
                           "elev": center.get("elev")}
    validate_plot(plot)
    return plot


# ---------------------------------------------------------------------------
# Legacy SST parameter CSV adapter
# ---------------------------------------------------------------------------

def is_tree_list_path(path: str) -> bool:
    """True for canonical tree-list JSON files (excludes *_report.json and
    *_env.json sidecars, which share the plot_*.json glob prefix)."""
    return (path.endswith(".json")
            and not path.endswith("_report.json")
            and not path.endswith("_env.json"))


def legacy_trees_frame(csv_path: str) -> pd.DataFrame:
    """Legacy SST parameter table -> downstream table frame.

    Supplies TREE_UID, FUNC_GROUP, X, Y and the plot_lat/plot_lon/plot_elev
    passthrough columns. DBH/HT/TREE_CODE cannot be recovered from frozen
    ranges -- they are omitted (generate_trees reads legacy CSVs directly).
    """
    df = pd.read_csv(csv_path, dtype={"TREE_UID": str})
    return df


def trees_frame(path: str) -> pd.DataFrame:
    """One loader for every downstream stage (assemble, stem/canopy maps).

    Accepts a canonical tree-list JSON or a legacy SST parameter CSV. Always
    returns columns: TREE_UID, X, Y (+ plot_lat/plot_lon/plot_elev, which
    are the per-plot center values the pipeline passes through; FIA lists
    whose input had no PLOT join get all-NaN geo columns).
    """
    if is_tree_list_path(path):
        plot = load_plot(path)
        df = plot["trees"].copy()
        center = plot["plot_center"]
        df["plot_lat"] = center.get("lat")
        df["plot_lon"] = center.get("lon")
        df["plot_elev"] = center.get("elev")
        return df
    return legacy_trees_frame(path)


def find_plots(plotdir: str) -> list:
    """Canonical plot tables in a directory: plot_*.json (tree lists,
    excluding report/env sidecars) plus legacy plot_*.csv. Sorted."""
    tables = []
    for name in sorted(os.listdir(plotdir)):
        full = os.path.join(plotdir, name)
        if not os.path.isfile(full):
            continue
        if is_tree_list_path(full) or name.endswith(".csv"):
            tables.append(full)
    return tables
