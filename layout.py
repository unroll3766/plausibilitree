"""
layout.py: the canonical output directory layout for a plot.

Every path in the restructured output is defined here exactly once, so the
structure can evolve (e.g. when the normalized variants are populated later)
without hunting through the stage scripts. Nothing writes files here; these
functions only build paths and, by default, create the parent directory.

Layout (per plot):

    PLOT_UID/
      individual_trees/
        PLOT_UID_TREE_UID/
          textual/          generation + skeleton JSONs
          las/              raw / trimmed / normalized clouds
          deep_learning/
            volume/labels/          per-tree volume JSON
            stem_map/features|labels
            canopy_map/features|labels   (+ per-tree bbox & hull JSONs)
      plot_placement/
        textual/            manifest.csv, scene JSON
        las/                raw / trimmed / normalized clouds
        deep_learning/
          volume/labels/            plot volume JSON
          stem_map/features|labels
          canopy_map/features|labels

"normalized" outputs are part of the layout but are not produced yet; their
path keys exist so callers can address them when that step is added.

Clouds are canonically LAZ: filenames end in `.laz` and laspy compresses /
decompresses transparently through the `lazrs` backend. The path helpers
fall back to a legacy `.las` file when only that exists, so plots generated
before the switch remain readable without regeneration.
"""

from __future__ import annotations

import os

INDIVIDUAL_TREES = "individual_trees"
INDIVIDUAL_SHRUBS = "individual_shrubs"
PLOT_PLACEMENT = "plot_placement"


# ---------------------------------------------------------------------------
# Tree-level artifacts:  key -> (subdir, filename template using {tag})
# ---------------------------------------------------------------------------

_TREE = {
    "generation":         ("textual", "{tag}_generation.json"),
    "raw_skeleton":       ("textual", "{tag}_raw_skeleton.json"),
    "trimmed_skeleton":   ("textual", "{tag}_trimmed_skeleton.json"),

    "raw_las":            ("las", "{tag}_raw.laz"),
    "trimmed_las":        ("las", "{tag}_trimmed.laz"),
    "terrain_las":        ("las", "{tag}_terrain.laz"),
    "understory_las":     ("las", "{tag}_understory.laz"),
    "normalized_las":     ("las", "{tag}_normalized.laz"),

    "volume_label":       ("deep_learning/volume/labels", "{tag}_volume.json"),

    "stem_feature_raw":   ("deep_learning/stem_map/features",
                           "{tag}_stem_bounding_box_raw.tif"),
    "stem_feature_norm":  ("deep_learning/stem_map/features",
                           "{tag}_stem_bounding_box_normalized.tif"),
    "stem_label_raw":     ("deep_learning/stem_map/labels",
                           "{tag}_stem_mask_raw.tif"),
    "stem_label_norm":    ("deep_learning/stem_map/labels",
                           "{tag}_stem_mask_normalized.tif"),

    "canopy_feature_raw": ("deep_learning/canopy_map/features",
                           "{tag}_canopy_bounding_box_raw.tif"),
    "canopy_feature_norm":("deep_learning/canopy_map/features",
                           "{tag}_canopy_bounding_box_normalized.tif"),
    "canopy_label_raw":   ("deep_learning/canopy_map/labels",
                           "{tag}_canopy_mask_raw.tif"),
    "canopy_label_norm":  ("deep_learning/canopy_map/labels",
                           "{tag}_canopy_mask_normalized.tif"),
    "canopy_bbox_json":   ("deep_learning/canopy_map/labels",
                           "{tag}_canopy_bounding_box.json"),
    "canopy_hull_json":   ("deep_learning/canopy_map/labels",
                           "{tag}_canopy_convex_hull.json"),
}


# ---------------------------------------------------------------------------
# Plot-level artifacts
# ---------------------------------------------------------------------------

_PLOT = {
    "manifest":           ("textual", "manifest.csv"),
    "scene":              ("textual", "{tag}_scene.json"),

    "raw_las":            ("las", "{tag}_raw.laz"),
    "trimmed_las":        ("las", "{tag}_trimmed.laz"),
    "terrain_las":        ("las", "{tag}_terrain.laz"),
    "understory_las":     ("las", "{tag}_understory.laz"),
    "normalized_las":     ("las", "{tag}_normalized.laz"),

    "volume_label":       ("deep_learning/volume/labels", "{tag}_volume.json"),

    "stem_feature_raw":   ("deep_learning/stem_map/features",
                           "{tag}_stem_map_raw.tif"),
    "stem_feature_norm":  ("deep_learning/stem_map/features",
                           "{tag}_stem_map_normalized.tif"),
    "stem_label_raw":     ("deep_learning/stem_map/labels",
                           "{tag}_stem_mask_raw.tif"),
    "stem_label_norm":    ("deep_learning/stem_map/labels",
                           "{tag}_stem_mask_normalized.tif"),

    "canopy_feature_raw": ("deep_learning/canopy_map/features",
                           "{tag}_canopy_map_raw.tif"),
    "canopy_feature_norm":("deep_learning/canopy_map/features",
                           "{tag}_canopy_map_normalized.tif"),
    "canopy_label_raw":   ("deep_learning/canopy_map/labels",
                           "{tag}_canopy_mask_raw.tif"),
    "canopy_label_norm":  ("deep_learning/canopy_map/labels",
                           "{tag}_canopy_mask_normalized.tif"),
}


# ---------------------------------------------------------------------------
# Shrub-level artifacts (mirrors trees; shrubs have raw / understory /
# normalized clouds, a volume JSON with primitive+mesh+voxel volumes, and
# canopy + DBH(stem) maps). Shrubs have no trimmed/skeleton stage.
# ---------------------------------------------------------------------------

_SHRUB = {
    "raw_las":            ("las", "{tag}_raw.laz"),
    "understory_las":     ("las", "{tag}_understory.laz"),
    "normalized_las":     ("las", "{tag}_normalized.laz"),

    "volume_label":       ("deep_learning/volume/labels", "{tag}_volume.json"),

    "stem_feature_raw":   ("deep_learning/stem_map/features",
                           "{tag}_stem_bounding_box_raw.tif"),
    "stem_feature_norm":  ("deep_learning/stem_map/features",
                           "{tag}_stem_bounding_box_normalized.tif"),
    "stem_label_raw":     ("deep_learning/stem_map/labels",
                           "{tag}_stem_mask_raw.tif"),
    "stem_label_norm":    ("deep_learning/stem_map/labels",
                           "{tag}_stem_mask_normalized.tif"),

    "canopy_feature_raw": ("deep_learning/canopy_map/features",
                           "{tag}_canopy_bounding_box_raw.tif"),
    "canopy_feature_norm":("deep_learning/canopy_map/features",
                           "{tag}_canopy_bounding_box_normalized.tif"),
    "canopy_label_raw":   ("deep_learning/canopy_map/labels",
                           "{tag}_canopy_mask_raw.tif"),
    "canopy_label_norm":  ("deep_learning/canopy_map/labels",
                           "{tag}_canopy_mask_normalized.tif"),
    "canopy_bbox_json":   ("deep_learning/canopy_map/labels",
                           "{tag}_canopy_bounding_box.json"),
    "canopy_hull_json":   ("deep_learning/canopy_map/labels",
                           "{tag}_canopy_convex_hull.json"),
}


def shrub_tag(plot_uid: str, shrub_uid: str) -> str:
    return f"{plot_uid}_shrub_{shrub_uid}"


def _resolve_cloud(path: str) -> str:
    """Cloud paths are canonically `.laz` (laspy writes LAZ-compressed data
    for a `.laz` target and decompresses transparently on read). A legacy
    `.las` file wins when the `.laz` one does not exist, so plots generated
    before the switch keep working without regeneration."""
    if path.endswith(".laz") and not os.path.exists(path):
        legacy = path[:-4] + ".las"
        if os.path.exists(legacy):
            return legacy
    return path


def cloud_file(directory: str, base: str) -> str:
    """Path of a point-cloud file whose name is `base` (no extension) in
    `directory`: the canonical `.laz` name, or the legacy `.las` file when
    only that exists. For flat-stage paths (generation/assembly) that do
    not go through the key-based path functions."""
    laz = os.path.join(directory, base + ".laz")
    legacy = os.path.join(directory, base + ".las")
    if not os.path.exists(laz) and os.path.exists(legacy):
        return legacy
    return laz


def shrub_dir(plot_dir: str, plot_uid: str, shrub_uid: str) -> str:
    return os.path.join(plot_dir, INDIVIDUAL_SHRUBS,
                        shrub_tag(plot_uid, shrub_uid))


def shrub_path(plot_dir: str, plot_uid: str, shrub_uid: str, key: str,
               create: bool = True) -> str:
    """Absolute path to a shrub-level artifact identified by `key`."""
    sub, tmpl = _SHRUB[key]
    tag = shrub_tag(plot_uid, shrub_uid)
    d = os.path.join(shrub_dir(plot_dir, plot_uid, shrub_uid), sub)
    if create:
        os.makedirs(d, exist_ok=True)
    return _resolve_cloud(os.path.join(d, tmpl.format(tag=tag)))


def shrub_rel(plot_dir: str, plot_uid: str, shrub_uid: str, key: str) -> str:
    return _rel(plot_dir, shrub_path(plot_dir, plot_uid, shrub_uid, key,
                                     create=False), plot_uid)


def tree_tag(plot_uid: str, tree_uid: str) -> str:
    return f"{plot_uid}_{tree_uid}"


def tree_dir(plot_dir: str, plot_uid: str, tree_uid: str) -> str:
    return os.path.join(plot_dir, INDIVIDUAL_TREES, tree_tag(plot_uid, tree_uid))


def tree_path(plot_dir: str, plot_uid: str, tree_uid: str, key: str,
              create: bool = True) -> str:
    """Absolute path to a tree-level artifact identified by `key`."""
    sub, tmpl = _TREE[key]
    tag = tree_tag(plot_uid, tree_uid)
    d = os.path.join(tree_dir(plot_dir, plot_uid, tree_uid), sub)
    if create:
        os.makedirs(d, exist_ok=True)
    return _resolve_cloud(os.path.join(d, tmpl.format(tag=tag)))


def plot_path(plot_dir: str, plot_uid: str, key: str,
              create: bool = True) -> str:
    """Absolute path to a plot-level artifact identified by `key`."""
    sub, tmpl = _PLOT[key]
    d = os.path.join(plot_dir, PLOT_PLACEMENT, sub)
    if create:
        os.makedirs(d, exist_ok=True)
    return _resolve_cloud(os.path.join(d, tmpl.format(tag=plot_uid)))


# ---------------------------------------------------------------------------
# Relative links for the scene JSON (from the collection root, i.e. the parent
# of the plot directory), using forward slashes for portability.
# ---------------------------------------------------------------------------

def _rel(plot_dir: str, abs_path: str, plot_uid: str) -> str:
    root = os.path.dirname(os.path.abspath(plot_dir))
    return os.path.relpath(abs_path, root).replace(os.sep, "/")


def tree_rel(plot_dir: str, plot_uid: str, tree_uid: str, key: str) -> str:
    return _rel(plot_dir, tree_path(plot_dir, plot_uid, tree_uid, key,
                                    create=False), plot_uid)


def plot_rel(plot_dir: str, plot_uid: str, key: str) -> str:
    return _rel(plot_dir, plot_path(plot_dir, plot_uid, key, create=False),
                plot_uid)


def tree_keys():
    return list(_TREE)


def plot_keys():
    return list(_PLOT)
