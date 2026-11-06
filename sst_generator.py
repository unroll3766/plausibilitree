"""
sst_generator.py: the Simple Synth Tree generator module (spec Sections
7.2-7.5), implementing the Section 7.1 contract.

Everything SST-specific lives here behind the GeneratorModule interface:

- the parameter dictionary is the functional-group table (Section 7.3;
  per-tree DBH/HT tolerances HT_TOL/DBH_TOL, default +/-1%, bound the
  realized size around the tree list values; crown width from the CW
  allometry, base diameter from the TAPER_B extrapolation);
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

DEFAULT_SIZE_TOL = 0.01        # +/- 1% (spec Sections 7.3 / 8.6)
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

    def entry_for(self, tree_code: str) -> dict:
        """The resolved dictionary entry for a tree code."""
        key = self.resolve_tree_code(tree_code)
        if key not in self.fg_db.index:
            key = "DEFAULT"
        row = self.fg_db.loc[key]
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

        Tolerances bound the realized size around the list values; crown
        width and base diameter come from the entry's allometry
        (Section 7.2)."""
        ht_tol = entry.get("HT_TOL", DEFAULT_SIZE_TOL) or DEFAULT_SIZE_TOL
        dbh_tol = entry.get("DBH_TOL", DEFAULT_SIZE_TOL) or DEFAULT_SIZE_TOL
        ht_lo, ht_hi = HT * (1 - ht_tol), HT * (1 + ht_tol)
        base_dia = _scalar_taper_base_diameter(DBH, HT, entry["TAPER_B"])
        crown_w = _scalar_crown_width(DBH, entry["CW_A0"], entry["CW_A1"])

        model_params = {
            "height_range": [ht_lo, ht_hi],
            "diam_range": [base_dia * (1 - dbh_tol),
                           base_dia * (1 + dbh_tol)],
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
        for poly in sk["trunk"]["polylines"]:
            arr = np.asarray(poly, dtype=float)
            arr[:, 0] += float(x)
            arr[:, 1] += float(y)
        for br in sk.get("branches", []):
            br["base"][0] += float(x)
            br["base"][1] += float(y)
            br["tip"][0] += float(x)
            br["tip"][1] += float(y)
        return model
