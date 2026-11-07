"""
tree_list.py: the canonical tree list (design spec Section 8.1) and the
one place in the codebase that understands its two serializations plus
the legacy SST parameter CSV.

Canonical form: the Section 8.1 attribute contract, serialized either as
one JSON object per plot,

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

or as one CSV row per tree (the notes.md default) carrying only the tree
rows -- never plot geometry: TREE_UID, X, Y, Z, DBH, HT, TREE_CODE,
optionally LON/LAT/ELEV. In CSV form the plot-level fields are supplied by
the sibling plot_<UID>_report.json the FIA helper writes (parameters used
plus plot_center); a hand-authored CSV without that sidecar defaults to an
unknown shape callers override (run_pipeline's --plot-shape etc.).

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

# Canonical CSV column order: the required attrs first, then the optional
# geo attrs; extra columns are preserved after these.
CSV_TREE_COLUMNS = ["TREE_UID", "X", "Y", "Z", "DBH", "HT", "TREE_CODE",
                    "LON", "LAT", "ELEV"]


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
    """Load a canonical tree list (JSON object or minimal CSV) -> plot dict
    with `trees` a DataFrame.

    plot_uid defaults from the filename when absent. Unknown/extra keys are
    preserved. plot_center defaults to all-null; for CSV lists, plot-level
    fields come from the sibling report (see _load_plot_csv)."""
    if path.endswith(".csv"):
        if not is_tree_list_csv(path):
            raise ValueError(f"{path}: not a canonical tree-list CSV "
                             f"(missing required attribute columns)")
        return _load_plot_csv(path)
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
# Canonical CSV (one row per tree; plot fields via the report sidecar)
# ---------------------------------------------------------------------------

def save_plot_csv(plot: dict, path: str) -> str:
    """Write one canonical tree-list CSV (Section 8.1): one row per tree,
    inventory facts only, in the CSV_TREE_COLUMNS order. Returns the path
    written."""
    if hasattr(plot.get("trees"), "to_dict"):
        records = plot["trees"].to_dict(orient="records")
    else:
        records = list(plot.get("trees", []))
    df = pd.DataFrame(records)
    cols = ([c for c in CSV_TREE_COLUMNS if c in df.columns]
            + [c for c in df.columns if c not in CSV_TREE_COLUMNS])
    df[cols].to_csv(path, index=False)
    return path


def is_tree_list_csv(path: str) -> bool:
    """True for canonical tree-list CSVs: a header carrying the required
    tree attributes (TREE_UID, X, Y, DBH, HT, TREE_CODE). Legacy SST
    parameter tables (FUNC_GROUP + HT_LO/BASE_DIA_LO ranges) return False."""
    if not path.endswith(".csv"):
        return False
    header = pd.read_csv(path, nrows=0)
    return set(REQUIRED_TREE_ATTRS).issubset(set(header.columns))


def _report_plot_fields(report_path: str) -> dict:
    """FIA report sidecar -> the plot-level fields of the Section 8.1
    contract. A minimal CSV tree list cannot carry plot geometry, so it is
    read against its sibling plot_<UID>_report.json: shape and dimensions
    come from the recorded placement parameters, the center from the
    report's plot_center when present."""
    with open(report_path) as fh:
        report = json.load(fh)
    params = report.get("parameters", {})
    return {"plot_shape": params.get("plot_shape"),
            "plot_radius_m": params.get("plot_radius_m"),
            "plot_length_m": params.get("plot_length_m"),
            "plot_width_m": params.get("plot_width_m"),
            "plot_side_m": params.get("plot_side_m"),
            "plot_center": (report.get("plot_center")
                            or {"lon": None, "lat": None, "elev": None})}


def _load_plot_csv(path: str) -> dict:
    """Canonical tree-list CSV -> plot dict. Plot-level fields come from the
    sibling plot_<UID>_report.json when it exists; a hand-authored CSV
    without one defaults to a circle with unknown dimensions, which callers
    override with their own geometry arguments."""
    df = pd.read_csv(path, dtype={"TREE_UID": str},
                     float_precision="round_trip")
    plot = {"plot_uid": os.path.splitext(os.path.basename(path))[0],
            "plot_shape": "circle", "plot_radius_m": None,
            "plot_length_m": None, "plot_width_m": None, "plot_side_m": None,
            "plot_center": {"lon": None, "lat": None, "elev": None}}
    report = os.path.splitext(path)[0] + "_report.json"
    if os.path.isfile(report):
        reported = _report_plot_fields(report)
        for key in ("plot_shape", "plot_radius_m", "plot_length_m",
                    "plot_width_m", "plot_side_m"):
            if reported[key] is not None:
                plot[key] = reported[key]
        plot["plot_center"] = reported["plot_center"]
    plot["trees"] = df
    validate_plot(plot)
    return plot


# ---------------------------------------------------------------------------
# Legacy SST parameter CSV adapter
# ---------------------------------------------------------------------------

def is_tree_list_path(path: str) -> bool:
    """True for canonical tree-list JSON files (excludes *_report.json
    reports and the replicate stage's *.env.json sidecars, which share the
    plot_*.json glob prefix)."""
    return (path.endswith(".json")
            and not path.endswith("_report.json")
            and not path.endswith(".env.json"))


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

    Accepts a canonical tree list (JSON or minimal CSV) or a legacy SST
    parameter CSV. Always
    returns columns: TREE_UID, X, Y (+ plot_lat/plot_lon/plot_elev, which
    are the per-plot center values the pipeline passes through; FIA lists
    whose input had no PLOT join get all-NaN geo columns).
    """
    if is_tree_list_path(path) or is_tree_list_csv(path):
        plot = load_plot(path)
        df = plot["trees"].copy()
        center = plot["plot_center"]
        df["plot_lat"] = center.get("lat")
        df["plot_lon"] = center.get("lon")
        df["plot_elev"] = center.get("elev")
        return df
    return legacy_trees_frame(path)


def find_plots(plotdir: str) -> list:
    """Canonical plot tables in a directory: plot_*.json tree lists
    (excluding *_report.json and *_env.json sidecars) plus legacy
    plot_*.csv tables. Only plot_-prefixed files qualify: directories can
    hold other json/csv artifacts (manifests, reports, env sidecars) that
    must never be treated as tables. Sorted."""
    tables = []
    for name in sorted(os.listdir(plotdir)):
        full = os.path.join(plotdir, name)
        if not os.path.isfile(full) or not name.startswith("plot_"):
            continue
        if is_tree_list_path(full) or name.endswith(".csv"):
            tables.append(full)
    return tables
