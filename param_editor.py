"""
param_editor.py: interactively tune tree parameters and preview sample trees.

Load the parameter dictionary ``functional_groups.json`` (spec Section 8.6;
legacy CSVs still open), edit any
cell directly, generate N sample trees from the selected row with
SimpleSynthTree, and page through them in a 3D viewer. The table is shown
transposed: one row per parameter, one column per group, with the parameter
names frozen in the leftmost column. Move between trees with the Left/Right
keys or the Prev/Next
buttons; drag to rotate a tree, press R (or the Reset view button) to restore
the default camera. The viewer can show or hide each point class -- trunk,
branch stems, foliage -- with the checkboxes or the 1/2/3 keys. The plot title
carries the label of the row that generated the current samples, so they stay
traceable after you move the selection elsewhere. Happy with the results? Save
the edited dictionary as JSON (or legacy CSV).

    python param_editor.py data/functional_groups.json
    python param_editor.py                     # opens a file picker

How a sample tree is drawn from a row:

  - Columns that map to treesim model params (SPLIT_HT_*, NUM_BRANCH_*,
    MIN_CAN_HT_*, MAX_CAN_W_HT_*, SPLIT_PROB, TOP_DIST, MID_DIST,
    FP_RAND_HT, SPR_*) are read from the selected row.
  - A functional-group table has no heights or diameters (those come from
    inventory rows at fia_helper time), so the GUI samples HT and DBH from
    the ranges in the sampling fields and derives the rest exactly the way
    fia_helper does: crown width from CW = exp(A0 + A1 * ln(DBH_in))
    in feet, and base diameter from the TAPER extrapolation
    d(0) = DBH * (H / (H - 1.37)) ** (1/b).
  - If the CSV *does* carry HT_LO/HT_HI, BASE_DIA_LO/HI or MAX_CAN_W columns
    (e.g. a hand-built plot-style table), those columns win over the sampling
    fields.

Generation is deterministic for a given seed: tree i reseeds numpy's global
RNG at ``seed + i * 7919`` (the same per-attempt reseeding pattern as
generate_trees.py) and retries degenerate draws a few times, so the same
seed always shows the same N trees.

Requires: numpy, pandas, matplotlib, tkinter (stdlib), trimesh (via treesim)
"""

from __future__ import annotations

import argparse
import json
import os
import time

import fia_helper

import numpy as np
import pandas as pd

import tkinter as tk
from tkinter import filedialog, messagebox, ttk

import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402

from treesim import gen_simtree

import trimesh

# Determinism of the sampled surface points comes from treesim's
# _seed_trimesh_sampling (installed at import above): trimesh >= 5 ignores
# np.random.seed without it, so a (row, seed, i) triple would build the same
# skeleton but scatter different points on every run.

# ---------------------------------------------------------------------------
# parameter mapping (mirrors generate_trees.py / fia_helper.py)
# ---------------------------------------------------------------------------

IN_TO_M = 0.0254
FT_TO_M = 0.3048
BREAST_HEIGHT_M = 1.37

RANGE_PARAMS = {
    "height_range": ("HT_LO", "HT_HI"),
    "diam_range": ("BASE_DIA_LO", "BASE_DIA_HI"),
    "split_height_range": ("SPLIT_HT_LO", "SPLIT_HT_HI"),
    "num_splits": ("SPLIT_NUM_LO", "SPLIT_NUM_HI"),
    "num_branches": ("NUM_BRANCH_LO", "NUM_BRANCH_HI"),
    "min_can_height": ("MIN_CAN_HT_LO", "MIN_CAN_HT_HI"),
    "max_can_width_height": ("MAX_CAN_W_HT_LO", "MAX_CAN_W_HT_HI"),
}
SCALAR_PARAMS = {
    "split_prob": "SPLIT_PROB",
    "max_can_width": "MAX_CAN_W",
    "tree_top_dist": "TOP_DIST",
    "tree_mid_dist": "MID_DIST",
    "foliage_noise": "FP_RAND_HT",
}
# treesim's own built-in example parameters: the fallback for any column a
# CSV happens to lack, so any "similar CSV" still previews
TREESIM_DEFAULTS = {
    "height_range": [30, 50], "diam_range": [0.5, 1.0],
    "split_height_range": [0.15, 0.5], "split_prob": 0.5,
    "num_splits": [1, 1],
    "num_branches": [60, 100], "min_can_height": [0.2, 0.5],
    "max_can_width": 7.0, "max_can_width_height": [0.4, 0.8],
    "tree_top_dist": 2.5, "tree_mid_dist": 0.5, "foliage_noise": 0.5,
}
DEFAULT_SAMPLE_HT = (25.0, 35.0)
DEFAULT_SAMPLE_DBH = (0.30, 0.60)
DEFAULT_STEM_FRAC = 0.3

TRUNK_COLOR = "#8B5A2B"      # trunk stems (skeleton polyline 0)
BRANCH_COLOR = "#A0522D"     # fork / branch stems (skeleton polylines 1+)
FOLIAGE_COLOR = "#3A8C4A"    # foliage


# ---------------------------------------------------------------------------
# headless core (usable / testable without the GUI)
# ---------------------------------------------------------------------------

def load_table(path):
    """Load a parameter dictionary (JSON, spec Section 8.6), or a legacy CSV.

    JSON: one object keyed by FUNC_GROUP whose values are that code's
    parameters (what the pipeline consumes). CSV: one row per code."""
    if str(path).endswith(".json"):
        with open(path) as fh:
            obj = json.load(fh)
        df = pd.DataFrame.from_dict({k: dict(v) for k, v in obj.items()},
                                    orient="index")
        df.index.name = "FUNC_GROUP"
        df = fia_helper.normalize_param_columns(df).reset_index()
    else:
        df = fia_helper.normalize_param_columns(pd.read_csv(path))
    if df.empty:
        raise RuntimeError(f"{path}: table has no rows")
    return df


def save_table(df, path):
    """Write the parameter dictionary as JSON (.json) or legacy CSV."""
    if str(path).endswith(".json"):
        d = df.set_index("FUNC_GROUP")
        obj = {code: {k: (None if pd.isna(v) else
                          (v.item() if hasattr(v, "item") else v))
                      for k, v in row.items()}
               for code, row in d.iterrows()}
        with open(path, "w") as fh:
            json.dump(obj, fh, indent=2)
            fh.write("\n")
        return
    df.to_csv(path, index=False)


def _f(row, col, default=None):
    """Numeric cell value, or `default` when missing/unparseable."""
    if col not in row.index:
        return default
    v = pd.to_numeric(pd.Series([row[col]]), errors="coerce").iloc[0]
    return default if not np.isfinite(v) else float(v)


def crown_width_m(dbh_m, a0, a1):
    """fia_helper.crown_width's equation: CW = exp(a0 + a1*ln(DBH_in)) ft."""
    dbh_in = max(dbh_m / IN_TO_M, 0.1)
    return float(np.exp(a0 + a1 * np.log(dbh_in)) * FT_TO_M)


def base_diameter_m(dbh_m, ht_m, taper_b):
    """fia_helper.dbh_to_base_diameter's Kozak/Ormerod-style extrapolation."""
    if ht_m > BREAST_HEIGHT_M + 0.1:
        return float(dbh_m * (ht_m / (ht_m - BREAST_HEIGHT_M)) ** (1.0 / taper_b))
    return float(dbh_m)


def sample_tree_params(row, rng, sample_ht, sample_dbh):
    """One sample tree's (model_params, info) drawn from one CSV row.

    Table columns that carry inventory-derived values (HT_*, BASE_DIA_*,
    MAX_CAN_W) win over the GUI sampling ranges; everything else comes from
    the row, falling back to treesim's example parameters when a column is
    absent. Heights and diameters are pinned to their drawn values (passed
    as [v, v] ranges) so the drawn DBH/HT and the generated tree agree.
    """
    p, info = {}, {}

    if "HT_LO" in row.index and "HT_HI" in row.index:
        ht = float(rng.uniform(_f(row, "HT_LO"), _f(row, "HT_HI")))
    else:
        ht = float(rng.uniform(*sample_ht))
    p["height_range"] = [ht, ht]
    info["height_m"] = ht

    if "BASE_DIA_LO" in row.index and "BASE_DIA_HI" in row.index:
        p["diam_range"] = [_f(row, "BASE_DIA_LO"), _f(row, "BASE_DIA_HI")]
        info["base_diam_m"] = p["diam_range"][0]
    else:
        dbh = float(rng.uniform(*sample_dbh))
        taper = _f(row, "TAPER", 1.6)
        base = base_diameter_m(dbh, ht, taper)
        p["diam_range"] = [base, base]
        info["dbh_m"], info["base_diam_m"] = dbh, base

    if "MAX_CAN_W" in row.index:
        p["max_can_width"] = _f(row, "MAX_CAN_W")
    elif "A0" in row.index and "A1" in row.index:
        # representative dbh for the allometry: the drawn one, or the row's
        # own base diameter when the table supplies diameters directly
        dbh_for_cw = info.get("dbh_m", p["diam_range"][0])
        p["max_can_width"] = crown_width_m(dbh_for_cw, _f(row, "A0"),
                                           _f(row, "A1"))
    else:
        p["max_can_width"] = TREESIM_DEFAULTS["max_can_width"]
    info["crown_w_m"] = p["max_can_width"]

    for key, (lo, hi) in RANGE_PARAMS.items():
        if key in ("height_range", "diam_range"):
            continue
        vlo, vhi = _f(row, lo), _f(row, hi)
        if vlo is not None and vhi is not None:
            v = [vlo, vhi] if vhi >= vlo else [vhi, vlo]
            p[key] = [int(round(x)) for x in v] \
                if key in ("num_branches", "num_splits") else v
        else:
            p[key] = TREESIM_DEFAULTS[key]
    for key, col in SCALAR_PARAMS.items():
        if key == "max_can_width":
            continue
        v = _f(row, col)
        p[key] = v if v is not None else TREESIM_DEFAULTS[key]
    return p, info


def stem_frac_from_row(row):
    """Stem-point fraction from SPR_LO/SPR_HI (list = sampled by treesim)."""
    lo, hi = _f(row, "SPR_LO"), _f(row, "SPR_HI")
    if lo is None:
        return DEFAULT_STEM_FRAC
    if hi is None or hi <= lo:
        return lo
    return [lo, hi]


def _polyline_distance(xyz, poly):
    """Distance from (n,3) points to one polyline, as (n,)."""
    d = np.full(len(xyz), np.inf)
    for s in range(len(poly) - 1):
        p0, p1 = poly[s], poly[s + 1]
        v = p1 - p0
        l2 = float((v * v).sum())
        if l2 <= 0.0:
            proj = np.broadcast_to(p0, xyz.shape)
        else:
            t = np.clip(((xyz - p0) @ v) / l2, 0.0, 1.0)
            proj = p0 + t[:, None] * v
        d = np.minimum(d, np.linalg.norm(xyz - proj, axis=1))
    return d


def point_classes(points, skeleton):
    """Per-point display class for one sample cloud: 1 trunk, 2 branch
    stems, 3 foliage.

    Foliage is the treesim class-0 cloud. Stem points (class 1) are
    assigned to their nearest trunk polyline: polyline 0 is the base stem
    (the trunk), polylines 1.. are forks -- the same main-stem/fork split
    the pipeline's semantic codes 1/2 use, so what the editor shows as
    "trunk" matches what assembly labels trunk.
    """
    cls = np.where(points[:, 3] == 1, 1, 3)
    if skeleton is None:
        return cls
    polys = [np.asarray(p, dtype=float) for p in skeleton["trunk"]["polylines"]]
    stem = points[:, 3] == 1
    if len(polys) <= 1 or not stem.any():
        return cls
    xyz = points[stem, :3]
    best = np.full(len(xyz), np.inf)
    best_i = np.zeros(len(xyz), dtype=int)
    for pi, poly in enumerate(polys):
        d = _polyline_distance(xyz, poly)
        closer = d < best
        best[closer] = d[closer]
        best_i[closer] = pi
    cls[stem] = np.where(best_i == 0, 1, 2)
    return cls


def generate_samples(row, n, *, seed=0, points=8192,
                     sample_ht=DEFAULT_SAMPLE_HT,
                     sample_dbh=DEFAULT_SAMPLE_DBH):
    """Generate `n` deterministic sample trees from one CSV row.

    Returns a list of (points, info, classes) with points (N, 4): x, y, z,
    treesim class (1 = stem, 0 = foliage), and classes the display class
    per point (1 trunk, 2 branch stems, 3 foliage). Tree i reseeds numpy's
    global RNG and retries degenerate draws, mirroring generate_trees.py's
    attempt loop.
    """
    stem_frac = stem_frac_from_row(row)
    samples = []
    for i in range(n):
        tree_seed = seed + i * 7919
        last_err = None
        for attempt in range(5):
            rng = np.random.default_rng(tree_seed + attempt * 104729 + 1)
            p, info = sample_tree_params(row, rng, sample_ht, sample_dbh)
            try:
                np.random.seed(tree_seed + attempt * 104729)
                pts, skel = gen_simtree(Np=points, stem_frac_lambda=stem_frac,
                                       model_params=p, return_skeleton=True)
                if len(pts) < 4 or not np.all(np.isfinite(pts[:, :3])):
                    raise ValueError("degenerate cloud")
                samples.append((pts, info, point_classes(pts, skel)))
                last_err = None
                break
            except Exception as e:               # noqa: BLE001
                last_err = e
        if last_err is not None:
            raise RuntimeError(
                f"sample {i + 1} of {n} failed after 5 attempts: {last_err}")
    return samples


# ---------------------------------------------------------------------------
# GUI
# ---------------------------------------------------------------------------

class ParamEditor(tk.Tk):
    """Parameter table (left) + sample-tree viewer (right)."""

    def __init__(self, csv_path):
        super().__init__()
        self.df = load_table(csv_path)
        self.path = csv_path
        self.samples = []           # list of (points, info, classes)
        self.samples_group = None  # label of the row the samples came from
        self.index = 0
        self.selected = 0
        # point-class visibility (checkboxes; keys 1/2/3)
        self.show_trunk = tk.BooleanVar(value=True)
        self.show_branches = tk.BooleanVar(value=True)
        self.show_foliage = tk.BooleanVar(value=True)
        self.title(self._title_text())
        self._build()
        self._select_row(0)

    @staticmethod
    def _treesim_stamp():
        """Marker for which simulator build is loaded -- the editor keeps its
        treesim module from launch, so a session opened before an edit keeps
        rendering the old geometry; this makes that visible at a glance."""
        src = os.path.join(os.path.dirname(os.path.abspath(__file__)), "treesim.py")
        try:
            t = time.localtime(os.path.getmtime(src))
            return time.strftime("treesim %Y-%m-%d %H:%M", t)
        except OSError:
            return "treesim ?"

    def _title_text(self):
        return (f"PlausibiliTree parameter editor - {os.path.basename(self.path)} "
                f"[{self._treesim_stamp()}]")

    # -- construction ------------------------------------------------------

    def _build(self):
        # a sane default window: without this the 19-column table's requested
        # width exceeds the screen and pushes the viewer off the right edge
        self.geometry("1500x880")
        self.minsize(980, 600)

        toolbar = ttk.Frame(self)
        toolbar.pack(side=tk.TOP, fill=tk.X, padx=6, pady=4)
        ttk.Button(toolbar, text="Open...", command=self._open).pack(side=tk.LEFT)
        ttk.Button(toolbar, text="Save As...",
                   command=self._save).pack(side=tk.LEFT, padx=(8, 0))
        self.status = ttk.Label(toolbar, text=f"{len(self.df)} row(s), "
                                              f"{len(self.df.columns)} column(s)")
        self.status.pack(side=tk.RIGHT)

        hint = ttk.Label(
            self, foreground="#555",
            text="the table is transposed: one row per parameter, one column "
                 "per group -- click a cell in a column to select that group  "
                 "->  set N trees / seed / points (and the sample HT + DBH "
                 "ranges)  ->  Generate  ->  page with Left/Right, R resets "
                 "the view, 1/2/3 toggle trunk / branch stems / foliage; "
                 "wheel scrolls the parameters, Shift+wheel the groups  ->  "
                 "Save As... (JSON)")
        hint.pack(side=tk.TOP, fill=tk.X, padx=8)

        body = ttk.Frame(self)
        body.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        # viewer first / from the right edge: under space pressure the table
        # (which scrolls) gives way, never the Generate controls
        self._build_viewer(body)
        self._build_table(body)

        self.bind("<Left>", lambda e: self._step(-1))
        self.bind("<Right>", lambda e: self._step(+1))
        self.bind("r", lambda e: self._reset_view())
        self.bind("R", lambda e: self._reset_view())
        self.bind("1", lambda e: self._toggle_class(self.show_trunk))
        self.bind("2", lambda e: self._toggle_class(self.show_branches))
        self.bind("3", lambda e: self._toggle_class(self.show_foliage))

    def _build_table(self, parent):
        frame = ttk.Frame(parent)
        frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=(6, 6), pady=4)
        nrows = len(self.df)
        ncols = len(self.df.columns)

        # -- transposed table: one row per parameter, one column per group --
        # the frozen left column keeps the parameter names visible while the
        # group columns scroll; entries keep their (group-row, param-column)
        # keys, so readback and save are unchanged
        key_holder = tk.Canvas(frame, highlightthickness=0)
        key_grid = ttk.Frame(key_holder)
        key_holder.create_window((0, 0), window=key_grid, anchor="nw")
        key_grid.bind("<Configure>",
                      lambda e: key_holder.configure(scrollregion=key_holder.bbox("all")))
        ttk.Label(key_grid, text="parameter", relief="groove",
                  padding=(4, 2)).grid(row=0, column=0, sticky="nsew")

        # -- scrollable remainder: one column per group (rows of the CSV) ---
        # the body canvas is the single source of truth for vertical
        # position: every scroll path (wheel, scrollbar, programmatic)
        # updates its yscrollcommand, which mirrors the frozen column along
        holder = tk.Canvas(frame, highlightthickness=0)
        ybar = ttk.Scrollbar(frame, orient=tk.VERTICAL, command=holder.yview)
        xbar = ttk.Scrollbar(frame, orient=tk.HORIZONTAL, command=holder.xview)

        def _body_yview_changed(first, last):
            ybar.set(first, last)
            key_holder.yview_moveto(float(first))

        holder.configure(yscrollcommand=_body_yview_changed,
                         xscrollcommand=xbar.set)
        grid = ttk.Frame(holder)
        holder.create_window((0, 0), window=grid, anchor="nw")
        grid.bind("<Configure>",
                  lambda e: holder.configure(scrollregion=holder.bbox("all")))

        # -- cells -----------------------------------------------------------
        # the CSV's first column (the group labels) becomes the header row
        # of the transposed view and stays editable, like every other cell;
        # clicking any cell in a column selects that group
        self.entries = {}
        for r in range(nrows):
            val = self.df.iat[r, 0]
            e0 = ttk.Entry(grid, width=12, font=("TkFixedFont", 9))
            e0.insert(0, "" if pd.isna(val) else str(val))
            e0.grid(row=0, column=r, sticky="nsew")
            e0.bind("<Button-1>", lambda ev, r=r: self._select_row(r), add="+")
            self.entries[(r, 0)] = e0
        for c in range(1, ncols):
            ttk.Label(key_grid, text=str(self.df.columns[c]), relief="groove",
                      padding=(4, 2)).grid(row=c, column=0, sticky="nsew")
            for r in range(nrows):
                val = self.df.iat[r, c]
                e = ttk.Entry(grid, width=9, font=("TkFixedFont", 9))
                e.insert(0, "" if pd.isna(val) else str(val))
                e.grid(row=c, column=r, sticky="nsew")
                e.bind("<Button-1>", lambda ev, r=r: self._select_row(r), add="+")
                self.entries[(r, c)] = e

        # pin every row (both grids) to the same pixel height: the frozen
        # labels and the body entries live in two canvases that scroll in
        # lockstep, so each canvas lays out its own rows -- and a label row
        # is naturally taller than an entry row, letting the parameter
        # names drift past their value rows as the offset accumulates
        row_h = max(w.winfo_reqheight()
                    for w in (*key_grid.winfo_children(),
                              *grid.winfo_children()))
        for i in range(ncols):
            key_grid.rowconfigure(i, minsize=row_h, uniform="krows")
            grid.rowconfigure(i, minsize=row_h, uniform="brows")

        # size the frozen column to its content (pixels), capped
        import tkinter.font as tkfont
        base = tkfont.nametofont("TkFixedFont")
        texts = ["parameter"] + [str(c) for c in self.df.columns]
        key_holder.configure(width=min(240, max(70, max(base.measure(t)
                                                        for t in texts) + 18)))

        key_holder.grid(row=1, column=0, sticky="ns")
        holder.grid(row=1, column=1, sticky="nsew")
        ybar.grid(row=1, column=2, sticky="ns")
        xbar.grid(row=2, column=1, sticky="ew")
        frame.rowconfigure(1, weight=1)
        frame.columnconfigure(1, weight=1)

        # wheel scrolling; the body canvas mirrors the frozen column itself.
        # X11 sends Button-4 (up) / Button-5 (down); wayland/win/mac send
        # <MouseWheel> with delta +/-120. Tk does not propagate events from
        # a child widget to its parent, and the table's visible area is
        # mostly Entry/Label cells, so the binding has to go on every cell
        # too -- not just the canvases and frames underneath them.
        # Plain wheel scrolls the parameters (down); Shift+wheel scrolls
        # the group columns sideways -- without it, the cells are
        # unreachable whenever the table is wider than the pane
        def _wheel(ev):
            if getattr(ev, "num", 0) == 4 or getattr(ev, "delta", 0) > 0:
                step = -1                      # up
            else:
                step = 1                       # down
            if getattr(ev, "state", 0) & 0x0001:          # Shift held
                holder.xview_scroll(step, "units")
            else:
                holder.yview_scroll(step, "units")

        for w in (key_holder, key_grid, holder, grid,
                  *key_grid.winfo_children(), *grid.winfo_children()):
            w.bind("<Button-4>", _wheel, add="+")
            w.bind("<Button-5>", _wheel, add="+")
            w.bind("<MouseWheel>", _wheel, add="+")

    def _build_viewer(self, parent):
        right = ttk.Frame(parent)
        right.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True, pady=4)

        # generation controls
        ctl = ttk.Frame(right)
        ctl.pack(side=tk.TOP, fill=tk.X)
        self._spin(ctl, "N trees:", "n", 8, 1, 200)
        self._spin(ctl, "Seed:", "seed", 0, 0, 10**9)
        self._spin(ctl, "Points:", "points", 8192, 256, 200000, inc=256)
        ttk.Label(ctl, text="Sample HT (m):").pack(side=tk.LEFT, padx=(12, 2))
        self.ht_lo = self._entry(ctl, DEFAULT_SAMPLE_HT[0])
        self.ht_hi = self._entry(ctl, DEFAULT_SAMPLE_HT[1])
        ttk.Label(ctl, text="DBH (m):").pack(side=tk.LEFT, padx=(12, 2))
        self.dbh_lo = self._entry(ctl, DEFAULT_SAMPLE_DBH[0])
        self.dbh_hi = self._entry(ctl, DEFAULT_SAMPLE_DBH[1])
        self.gen_btn = ttk.Button(ctl, text="Generate", command=self._generate)
        self.gen_btn.pack(side=tk.LEFT, padx=(12, 0))

        # 3D view
        self.fig = Figure(figsize=(5.2, 6.2), dpi=100)
        self.ax = self.fig.add_subplot(111, projection="3d")
        # the camera the viewer opens with; R (Reset view) restores it
        self.home_view = (self.ax.elev, self.ax.azim)
        self.zoom = 1.0        # 3D axes ignore the wheel by default; wired below
        self.span = None
        self.canvas = FigureCanvasTkAgg(self.fig, master=right)
        self.canvas.mpl_connect("scroll_event", self._on_scroll)
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        # pager
        pager = ttk.Frame(right)
        pager.pack(side=tk.TOP, fill=tk.X)
        ttk.Button(pager, text="< Prev", command=lambda: self._step(-1)).pack(side=tk.LEFT)
        ttk.Button(pager, text="Reset view",
                  command=self._reset_view).pack(side=tk.LEFT, padx=(8, 0))
        ttk.Button(pager, text="Next >", command=lambda: self._step(+1)).pack(side=tk.RIGHT)
        self.nav = ttk.Label(pager, text="no samples yet - pick a row and Generate")
        self.nav.pack(side=tk.TOP)
        self.info = ttk.Label(pager, text="", justify=tk.CENTER)
        self.info.pack(side=tk.TOP)
        ttk.Label(pager, text="Left/Right keys page; R resets the view; "
                             "drag to rotate; scroll to zoom").pack(side=tk.TOP)

        # point-class visibility (keys 1/2/3 do the same)
        toggles = ttk.Frame(right)
        toggles.pack(side=tk.TOP, fill=tk.X)
        for text, var in (("trunk", self.show_trunk),
                          ("branch stems", self.show_branches),
                          ("foliage", self.show_foliage)):
            ttk.Checkbutton(toggles, text=text, variable=var,
                            command=self._redraw).pack(side=tk.LEFT,
                                                        padx=(12, 0))

    def _spin(self, parent, label, attr, default, lo, hi, inc=1):
        ttk.Label(parent, text=label).pack(side=tk.LEFT, padx=(12, 2))
        sv = tk.StringVar(value=str(default))
        ttk.Spinbox(parent, from_=lo, to=hi, increment=inc, textvariable=sv,
                    width=7).pack(side=tk.LEFT)
        setattr(self, f"ctl_{attr}", sv)

    def _entry(self, parent, default):
        e = ttk.Entry(parent, width=6)
        e.insert(0, f"{default:g}")
        e.pack(side=tk.LEFT)
        return e

    # -- table <-> df ------------------------------------------------------

    def _read_table(self):
        """Pull the edited cell values back into self.df (numeric where the
        cell parses as a number, string otherwise)."""
        for (r, c), e in self.entries.items():
            txt = e.get()
            try:
                self.df.iat[r, c] = float(txt)
            except ValueError:
                self.df.iat[r, c] = txt

    def _row(self, r):
        return self.df.iloc[r]

    def _select_row(self, r):
        self.selected = r
        # ttk.Entry has no per-widget background, so the active row is shown
        # in the status bar instead
        try:
            key = str(self.df.iloc[r, 0])
        except Exception:
            key = f"row {r}"
        self.status.configure(text=f"selected row: {key}  |  "
                              f"{len(self.df)} row(s), "
                              f"{len(self.df.columns)} column(s)")

    # -- actions -----------------------------------------------------------

    def _generate(self):
        self._read_table()
        row = self._row(self.selected)
        n = max(1, int(self._float(self.ctl_n.get(), 8)))
        seed = int(self._float(self.ctl_seed.get(), 0))
        pts_n = max(64, int(self._float(self.ctl_points.get(), 8192)))
        ht = (self._float(self.ht_lo.get(), DEFAULT_SAMPLE_HT[0]),
              self._float(self.ht_hi.get(), DEFAULT_SAMPLE_HT[1]))
        dbh = (self._float(self.dbh_lo.get(), DEFAULT_SAMPLE_DBH[0]),
               self._float(self.dbh_hi.get(), DEFAULT_SAMPLE_DBH[1]))
        self.gen_btn.configure(state="disabled", text="Generating...")
        self.update_idletasks()
        try:
            self.samples = generate_samples(row, n, seed=seed, points=pts_n,
                                            sample_ht=ht, sample_dbh=dbh)
            self.samples_group = self._row_label(self.selected)
            self.index = 0
            self._show()
        except Exception as e:                     # noqa: BLE001
            messagebox.showerror("generation failed", str(e))
        finally:
            self.gen_btn.configure(state="normal", text="Generate")

    def _row_label(self, r):
        """Short label for a table row: its first-column value (e.g. the
        functional-group code), so generated samples stay traceable after
        the user moves the selection elsewhere."""
        try:
            v = self.df.iloc[r, 0]
            txt = "" if pd.isna(v) else str(v)
        except Exception:
            txt = ""
        return txt or f"row {r}"

    @staticmethod
    def _float(txt, default):
        try:
            return float(txt)
        except (TypeError, ValueError):
            return default

    def _step(self, d):
        if self.samples:
            self.index = (self.index + d) % len(self.samples)
            self._show()

    def _toggle_class(self, var):
        var.set(not var.get())
        self._redraw()

    def _redraw(self):
        """Redraw the current sample after a point-class visibility change
        (a checkbox or the 1/2/3 keys)."""
        if self.samples:
            self._show()

    def _reset_view(self):
        """Restore the default camera (elevation, azimuth and zoom); bound
        to the R key and the Reset view button. Axes3D.clear() does not touch
        the camera, so the restore has to set it explicitly."""
        if not self.samples:
            return
        self.zoom = 1.0
        self.ax.view_init(*self.home_view)
        self._show()

    def _on_scroll(self, event):
        """Wheel = zoom. 3D axes don't wire the mouse wheel by default, so
        this applies the wheel to set_box_aspect's zoom factor directly."""
        if event.inaxes is not self.ax or self.span is None:
            return
        self.zoom = float(np.clip(self.zoom * 1.15 ** event.step, 0.3, 12.0))
        self.ax.set_box_aspect(self.span, zoom=self.zoom)
        self.canvas.draw_idle()

    def _show(self):
        pts, info, cls = self.samples[self.index]
        ax = self.ax
        ax.clear()
        xyz = pts[:, :3]
        # frame on the full cloud, not the visible classes, so toggling a
        # class on or off never moves the camera
        lo, hi = xyz.min(axis=0), xyz.max(axis=0)
        span = np.maximum(hi - lo, 1e-6)
        self.span = tuple(float(s) for s in span)
        ax.set_box_aspect(self.span, zoom=self.zoom)
        shown = []
        if self.show_trunk.get():
            m = cls == 1
            shown.append(ax.scatter(xyz[m, 0], xyz[m, 1], xyz[m, 2],
                                    s=2, c=TRUNK_COLOR, label="trunk",
                                    depthshade=True))
        if self.show_branches.get():
            m = cls == 2
            shown.append(ax.scatter(xyz[m, 0], xyz[m, 1], xyz[m, 2],
                                    s=2, c=BRANCH_COLOR, label="branch stems",
                                    depthshade=True))
        if self.show_foliage.get():
            m = cls == 3
            shown.append(ax.scatter(xyz[m, 0], xyz[m, 1], xyz[m, 2],
                                    s=2, c=FOLIAGE_COLOR, label="foliage",
                                    depthshade=True))
        ax.set_xlabel("x (m)"); ax.set_ylabel("y (m)"); ax.set_zlabel("z (m)")
        if shown:
            ax.legend(loc="upper left", fontsize=8)
        group = self.samples_group or f"row {self.selected}"
        ax.set_title(f"{group} - tree {self.index + 1} of {len(self.samples)}")
        self.nav.configure(text=f"{self.index + 1} / {len(self.samples)}   "
                                f"of {group}   (Left/Right to page, R to reset)")
        bits = [f"HT {info.get('height_m', float('nan')):.1f} m"]
        if "dbh_m" in info:
            bits.append(f"DBH {info['dbh_m'] * 100:.1f} cm")
        bits.append(f"crown width {info.get('crown_w_m', float('nan')):.1f} m")
        bits.append(f"{len(pts):,} points")
        self.info.configure(text="   |   ".join(bits))
        self.canvas.draw_idle()

    def _open(self):
        path = filedialog.askopenfilename(
            title="Open parameter dictionary",
            filetypes=[("JSON dictionary", "*.json"),
                       ("Legacy CSV", "*.csv"), ("All files", "*.*")])
        if not path:
            return
        self._read_table()          # don't silently lose current edits
        self.df = load_table(path)
        self.path = path
        self.samples = []
        self.index = 0
        self._rebuild()

    def _save(self):
        path = filedialog.asksaveasfilename(
            title="Save edited dictionary as",
            defaultextension=".json", initialfile=os.path.basename(self.path),
            filetypes=[("JSON dictionary", "*.json"),
                       ("Legacy CSV", "*.csv")])
        if not path:
            return
        self._read_table()
        save_table(self.df, path)
        self.status.configure(text=f"saved {len(self.df)} row(s) -> {path}")

    def _show_empty(self):
        self.ax.clear()
        self.nav.configure(text="no samples yet - pick a row and Generate")
        self.info.configure(text="")
        self.canvas.draw_idle()

    def _rebuild(self):
        """Tear down the widgets and rebuild for the (new) table."""
        for child in list(self.children.values()):
            child.destroy()
        self.title(self._title_text())
        self._build()
        self._select_row(0)
        self._show_empty()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Interactively edit the tree parameter dictionary "
                    "(JSON, spec Section 8.6; legacy CSV) and preview "
                    "sample trees generated with SimpleSynthTree")
    ap.add_argument("dictionary", nargs="?", default=None,
                    help="parameter dictionary JSON (default: opens a file picker)")
    args = ap.parse_args(argv)

    path = args.dictionary
    if path is None:
        root = tk.Tk()
        root.withdraw()
        path = filedialog.askopenfilename(
            title="Open parameter dictionary",
            filetypes=[("JSON dictionary", "*.json"),
                       ("Legacy CSV", "*.csv"), ("All files", "*.*")])
        root.destroy()
        if not path:
            return 0
    app = ParamEditor(path)
    app.mainloop()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
