"""
sst_generator.py: the Simple Synth Tree generator module (spec Sections
7.2-7.5), implementing the Section 7.1 contract.

Everything SST-specific lives here behind the GeneratorModule interface:

- the parameter dictionary is the functional-group table (Section 7.3;
  crown width from the A0/A1 allometry, base diameter from the TAPER
  extrapolation; realized sizes drawn within the +/-1% generation
  constant around the tree list values; per-code DBH/HT applicability
  bounds gate entry resolution -- an out-of-bounds tree falls through
  the resolution order to DEFAULT);
- `build` translates an entry (or an explicit legacy SST row's frozen
  ranges) into treesim's model_params and calls gen_simtree at the tree
  base (tx = ty = 0) so the returned TreeModel is base-origin (Section 7.1);
- treesim draws from numpy's GLOBAL generator, so `build` bridges the rng
  the core hands it: one draw from rng seeds the global state (this is the
  documented, deterministic bridge -- the tree's shape derives entirely
  from the core's rng);
- SST's degenerate-build guards (non-finite coordinates, stem excursions
  past TOP_DIST + splits * MAX_FORK_OFFSET) raise ValueError so the core
  retries with a rehashed seed (MAX_ATTEMPTS);
- the skeleton is the module format documented in Section 7.4; `place`
  shifts base_xy, the trunk polylines and each branch's base/tip into plot
  coordinates.

Legacy support: build() also accepts an explicit params dict carrying
pre-expanded ranges (the frozen rows of the legacy SST parameter CSV),
through legacy_row_params() below -- both input contracts drive the same
build path.
"""

from __future__ import annotations

import math

import numpy as np

import fia_helper
from generator_module import GeneratorModule, TreeModel, register
from treesim import gen_simtree

# legacy plot-table column -> treesim model_params key
_RANGE_PARAMS = {
    "height_range": ("HT_LO", "HT_HI"),
    "diam_range": ("BASE_DIA_LO", "BASE_DIA_HI"),
    "split_height_range": ("SPLIT_HT_LO", "SPLIT_HT_HI"),
    "num_branches": ("NUM_BRANCH_LO", "NUM_BRANCH_HI"),
    "min_can_height": ("MIN_CAN_HT_LO", "MIN_CAN_HT_HI"),
    "max_can_width_height": ("MAX_CAN_W_HT_LO", "MAX_CAN_W_HT_HI"),
}
_SCALAR_PARAMS = {
    "split_prob": "SPLIT_PROB",
    "max_can_width": "MAX_CAN_W",
    "tree_top_dist": "TOP_DIST",
    "tree_mid_dist": "MID_DIST",
    "foliage_noise": "FP_RAND_HT",
}

DEFAULT_SIZE_TOL = 0.01        # +/-1% generation constant bounding realized
# sizes around the tree list values (spec Sections 7.3 / 8.6; W18)
DEFAULT_POINTS_PER_METRE = 1024   # overridden by the core's budget setting
# treesim offsets fork tops by up to 2.0*rand + 1.5 metres from the main stem
MAX_FORK_OFFSET = 3.5


def _scalar_taper_base_diameter(dbh_m: float, ht_m: float, b: float) -> float:
    """DBH extrapolated to ground line (fia_helper's taper, scalar form)."""
    if ht_m <= 1.37 + 0.1:
        return float(dbh_m)
    return float(dbh_m) * (ht_m / (ht_m - 1.37)) ** (1.0 / b)


def _scalar_crown_width(dbh_m: float, a0: float, a1: float) -> float:
    """CW (ft) = exp(a0 + a1 ln(DBH_in)), converted to metres."""
    dbh_in = max(float(dbh_m) / 0.0254, 0.1)
    return math.exp(a0 + a1 * math.log(dbh_in)) * 0.3048


def _max_stem_excursion(skeleton: dict) -> float:
    """Furthest any trunk vertex strays horizontally from the tree base."""
    bx, by = skeleton["base_xy"]
    worst = 0.0
    for poly in skeleton["trunk"]["polylines"]:
        arr = np.asarray(poly, dtype=float)
        worst = max(worst, float(np.hypot(arr[:, 0] - bx,
                                          arr[:, 1] - by).max()))
    return worst


def element_graph(skeleton: dict) -> list:
    """The notes.md element graph (spec Section 7.4) from the SST skeleton.

    One order-0 element per trunk stem (stem 0 is the root; each split k
    appends daughter stems 2k+1 and 2k+2, both children of the split's
    parent stem) and one order-1 element per branch, sourced to the stem
    element recorded at generation. Elements carry no coordinates -- only
    UIDs, parents, lengths, radii and cylinder volumes -- so they are
    invariant under place()'s rotation and reflection and can be copied
    verbatim into trimmed skeletons."""
    def stem_uid(i: int) -> str:
        return f"T{i:02d}"

    r0 = float(skeleton["base_radius_m"])
    h = float(skeleton["height_m"])
    segments = skeleton["trunk"]["segments"]
    branches = skeleton["branches"]

    # stem index -> parent stem index (root trunk has none)
    parent_of = {}
    for k, sp in enumerate(skeleton.get("splits", [])):
        parent_of[2 * k + 1] = int(sp["parent_polyline"])
        parent_of[2 * k + 2] = int(sp["parent_polyline"])

    elements = []
    children = {i: [] for i in range(len(segments) + len(branches))}
    for i, seg in enumerate(segments):
        z_base = float(seg["z_range"][0])
        radius = max(0.0, r0 * (1.0 - z_base / h)) if h > 0 else 0.0
        length = float(seg["length_m"])
        parent = parent_of.get(i)
        if parent is not None:
            children[parent].append(i)
        elements.append({
            "BRANCH_UID": stem_uid(i),
            "BRANCH_ORDER": 0,
            "BRANCH_LENGTH": length,
            "BRANCH_RADIUS": radius,
            "BRANCH_VOLUME": float(np.pi) * radius ** 2 * length,
            "SOURCE": (stem_uid(parent) if parent is not None else None),
            "DESTINATION": [],
        })
    for k, br in enumerate(branches):
        i = int(br.get("stem_index", 0))
        length = float(br["length_m"])
        radius = float(br.get("radius_m", 0.0))
        children[i].append(len(segments) + k)
        elements.append({
            "BRANCH_UID": f"B{k:03d}",
            "BRANCH_ORDER": 1,
            "BRANCH_LENGTH": length,
            "BRANCH_RADIUS": radius,
            "BRANCH_VOLUME": float(np.pi) * radius ** 2 * length,
            "SOURCE": stem_uid(i),
            "DESTINATION": [],
        })

    for j, e in enumerate(elements):
        e["DESTINATION"] = [elements[c]["BRANCH_UID"] for c in children[j]]
    uids = {e["BRANCH_UID"] for e in elements}
    assert all(e["SOURCE"] in uids or e["SOURCE"] is None
               for e in elements)
    return elements


@register
class SSTGenerator(GeneratorModule):
    """Simple Synth Tree behind the Section 7.1 contract."""

    name = "sst"

    def __init__(self, param_dict_path: str | None = None,
                 spcd_map_path: str | None = None):
        self._param_dict_path = param_dict_path
        self._spcd_map_path = spcd_map_path
        self._fg_db = None
        self._spcd_map = None

    # -- parameter dictionary ------------------------------------------

    def load_param_dict(self, path: str | None = None):
        """The functional-group table (CSV until W13), indexed by key."""
        explicit = path or self._param_dict_path
        self._fg_db = (fia_helper.load_func_group_db(explicit)
                       if explicit else fia_helper.load_func_group_db())
        return self._fg_db

    def _spcd(self):
        if self._spcd_map is None:
            self._spcd_map = (fia_helper.load_spcd_map(self._spcd_map_path)
                              if self._spcd_map_path
                              else fia_helper.load_spcd_map())
        return self._spcd_map

    @property
    def fg_db(self):
        if self._fg_db is None:
            self.load_param_dict()
        return self._fg_db

    def resolve_tree_code(self, code: str) -> str:
        """Spec Section 8.1 resolution order: species code -> functional
        group label -> DEFAULT."""
        code = str(code)
        known = set(self.fg_db.index)
        if code in known:
            return code
        mapped = self._spcd().get(code)
        if mapped is not None and mapped in known:
            return mapped
        try:
            mapped = self._spcd().get(int(float(code)))
            if mapped in known:
                return mapped
        except (ValueError, TypeError):
            pass
        return "DEFAULT"

    def _spcd_mapped(self, code: str) -> str | None:
        """The SPCD map hit for a numeric code string, or None."""
        try:
            return self._spcd().get(int(float(code)))
        except (ValueError, TypeError):
            return None

    @staticmethod
    def _entry_applicable(row, dbh: float | None, ht: float | None) -> bool:
        """Applicability bounds check (spec Section 7.3): both given sizes
        must lie inside the entry's DBH/HT bounds. A missing or NaN bound
        is unlimited; missing sizes are not checked (key-only resolution)."""
        for value, lo_col, hi_col in ((dbh, "DBH_MIN", "DBH_MAX"),
                                      (ht, "HT_MIN", "HT_MAX")):
            if value is None:
                continue
            lo, hi = row.get(lo_col), row.get(hi_col)
            try:
                if lo is not None and np.isfinite(lo) and value < float(lo):
                    return False
                if hi is not None and np.isfinite(hi) and value > float(hi):
                    return False
            except TypeError:
                continue
        return True

    def entry_for(self, tree_code: str, dbh: float | None = None,
                  ht: float | None = None) -> dict:
        """The resolved dictionary entry for a tree code (spec Section 7.3).

        Candidate keys in resolution order: the code itself, then its
        SPCD-mapped group, then DEFAULT. With sizes given, a candidate only
        applies when the tree's DBH/HT lie inside its applicability bounds
        (DBH_MIN/DBH_MAX, HT_MIN/HT_MAX; a missing bound is unlimited);
        the first applicable candidate wins, and DEFAULT is the
        unconditional last resort.
        """
        code = str(tree_code)
        known = set(self.fg_db.index)
        candidates = []
        if code in known:
            candidates.append(code)
        for mapped in (self._spcd().get(code), self._spcd_mapped(code)):
            if mapped is not None and mapped in known and mapped not in candidates:
                candidates.append(mapped)
        row = None
        for key in candidates:
            if self._entry_applicable(self.fg_db.loc[key], dbh, ht):
                row = self.fg_db.loc[key]
                break
        if row is None:
            row = self.fg_db.loc["DEFAULT"]
        out = {}
        for col, v in row.items():
            try:
                out[col] = None if np.isnan(v) else float(v)
            except TypeError:      # non-numeric column (e.g. DESCRIPTION)
                out[col] = v
        return out

    # -- params expansion ----------------------------------------------

    def _expand_entry(self, DBH: float, HT: float, entry: dict,
                       points_per_metre: int) -> dict:
        """A dictionary entry + tree list sizes -> explicit model params.

        The realized-size ranges are drawn within the +/-1% generation
        constant around the list values (the dictionary's tolerance keys
        are gone, spec Section 8.6); crown width and base diameter come
        from the entry's allometry (Section 7.2)."""
        ht_lo, ht_hi = HT * (1 - DEFAULT_SIZE_TOL), HT * (1 + DEFAULT_SIZE_TOL)
        base_dia = _scalar_taper_base_diameter(DBH, HT, entry["TAPER"])
        crown_w = _scalar_crown_width(DBH, entry["A0"], entry["A1"])

        model_params = {
            "height_range": [ht_lo, ht_hi],
            "diam_range": [base_dia * (1 - DEFAULT_SIZE_TOL),
                           base_dia * (1 + DEFAULT_SIZE_TOL)],
            "split_height_range": [entry["SPLIT_HT_LO"], entry["SPLIT_HT_HI"]],
            "split_prob": entry["SPLIT_PROB"],
            "num_branches": [int(entry["NUM_BRANCH_LO"]),
                             int(entry["NUM_BRANCH_HI"])],
            "min_can_height": [entry["MIN_CAN_HT_LO"], entry["MIN_CAN_HT_HI"]],
            "max_can_width_height": [entry["MAX_CAN_W_HT_LO"],
                                     entry["MAX_CAN_W_HT_HI"]],
            "max_can_width": crown_w,
            "tree_top_dist": entry["TOP_DIST"],
            "tree_mid_dist": entry["MID_DIST"],
            "foliage_noise": entry["FP_RAND_HT"],
        }
        if entry.get("SPLIT_NUM_LO") is not None:
            model_params["num_splits"] = [int(entry["SPLIT_NUM_LO"]),
                                          int(entry["SPLIT_NUM_HI"])]
        return {
            "model_params": model_params,
            "stem_frac_range": [entry["SPR_LO"], entry["SPR_HI"]],
            "points_per_metre": int(points_per_metre),
        }

    def legacy_row_params(self, row, points_per_metre: int) -> dict:
        """A legacy SST parameter-table row -> explicit build params."""
        p = {}
        for key, (lo, hi) in _RANGE_PARAMS.items():
            p[key] = [float(row[lo]), float(row[hi])]
        for key, col in _SCALAR_PARAMS.items():
            p[key] = float(row[col])
        p["num_branches"] = [int(round(v)) for v in p["num_branches"]]
        if "SPLIT_NUM_LO" in row.index and "SPLIT_NUM_HI" in row.index:
            p["num_splits"] = [int(round(float(row["SPLIT_NUM_LO"]))),
                               int(round(float(row["SPLIT_NUM_HI"])))]
        return {
            "model_params": p,
            "stem_frac_range": [float(row["SPR_LO"]), float(row["SPR_HI"])],
            "points_per_metre": int(points_per_metre),
        }

    # -- build -----------------------------------------------------------

    def build(self, rng: np.random.Generator, DBH: float, HT: float,
              tree_code: str, params: dict) -> TreeModel:
        if "model_params" not in params:
            params = self._expand_entry(DBH, HT, params,
                                        params.get("points_per_metre",
                                                   DEFAULT_POINTS_PER_METRE))
        model_params = params["model_params"]
        ppm = int(params["points_per_metre"])
        lo, hi = params["stem_frac_range"]
        stem_frac = [lo, hi] if hi > lo else lo

        ht_hi = float(model_params["height_range"][1])
        base_dia_hi = float(model_params["diam_range"][1])
        for v in [ht_hi, base_dia_hi, float(model_params["max_can_width"])]:
            if not (np.isfinite(v) and v > 0):
                raise ValueError(f"non-physical model parameters ({v})")
        np_target = max(1, int(round(ppm * math.ceil(ht_hi) * base_dia_hi)))

        # treesim draws from numpy's global generator: bridge the core's
        # rng into it with one deterministic draw (all shape entropy
        # derives from the core-provided rng; spec reproducibility holds)
        np.random.seed(int(rng.integers(0, 2 ** 31 - 1)))

        points, skeleton = gen_simtree(
            Np=np_target, stem_frac_lambda=stem_frac,
            model_params=model_params, tx=0.0, ty=0.0,
            return_skeleton=True)

        if not np.all(np.isfinite(points[:, :3])):
            raise ValueError("simulator produced non-finite coordinates")

        # divergence guard: a stem should never wander further from its
        # base than its lean plus the fork offsets allow
        split_hi = 1
        if "num_splits" in model_params:
            split_hi = max(1, int(round(model_params["num_splits"][1])))
        limit = (float(model_params["tree_top_dist"])
                 + split_hi * MAX_FORK_OFFSET + 1.0)
        excursion = _max_stem_excursion(skeleton)
        if excursion > limit:
            raise ValueError(f"stem excursion {excursion:.1f} m exceeds "
                             f"{limit:.1f} m")

        skeleton["elements"] = element_graph(skeleton)

        import schema
        old_sem = points[:, 3].astype(int)        # 1 stem, 0 foliage
        labels = np.where(old_sem == 1, schema.SEM_TREE,
                          schema.SEM_FOLIAGE).astype(np.int8)
        pts = np.column_stack([points[:, :3], labels]).astype(float)

        return TreeModel(
            points=pts,
            labels=labels,
            skeleton=skeleton,
            realized={"height_m": float(skeleton["height_m"]),
                     "base_radius_m": float(skeleton["base_radius_m"]),
                     "fork": int(skeleton["fork"]),
                     "n_points": int(len(pts)),
                     "n_stem_points": int((old_sem == 1).sum()),
                     "n_foliage_points": int((old_sem == 0).sum())},
            provenance={"model_params": model_params,
                        "stem_frac_lambda": stem_frac,
                        "target_points": np_target,
                        "points_per_metre": ppm,
                        "ht_hi_used_for_budget": ht_hi,
                        "base_dia_hi_used_for_budget": base_dia_hi,
                        "excursion_m": excursion,
                        "generator": "SimpleSynthTree / treesim.gen_simtree "
                                     "(skeleton patch)",
                        },
        )

    # -- placement -------------------------------------------------------

    def place(self, model: TreeModel, x: float, y: float) -> TreeModel:
        """Shift the base-origin model into plot coordinates (positions
        were measured to be translation-invariant in treesim, so this
        equals having generated at (x, y))."""
        model.points[:, 0] += float(x)
        model.points[:, 1] += float(y)
        sk = model.skeleton
        sk["base_xy"] = [sk["base_xy"][0] + float(x),
                         sk["base_xy"][1] + float(y)]
        polys = sk["trunk"]["polylines"]
        for i, poly in enumerate(polys):
            # polylines are plain lists (gen_simtree emits .tolist()): the
            # shift must be written back -- np.asarray on a list COPIES, so
            # mutating the array alone silently drops the shift (the trunk
            # polylines would stay base-local while points, base_xy and the
            # branches move to plot coords)
            arr = np.asarray(poly, dtype=float)
            arr[:, 0] += float(x)
            arr[:, 1] += float(y)
            polys[i] = arr.tolist()
        for br in sk.get("branches", []):
            br["base"][0] += float(x)
            br["base"][1] += float(y)
            br["tip"][0] += float(x)
            br["tip"][1] += float(y)
        return model
