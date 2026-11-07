"""
scene_editor.py: interactively tune the environment of an already-generated
plot -- topography, ground, shrubs and grass -- with immediate visual
feedback in a 3D viewer.

Point the editor at a run directory (or one plot table) whose plots have been
generated and restructured (the canonical layout with `*_trimmed.las` exists
in `plot_placement/las/`). Pick a plot, edit the terrain and understory
settings, press Generate: the editor runs `terrain.py` and `understory.py` on
that plot as subprocesses (exactly like `run_pipeline.py` does), then shows
the resulting cloud colored by semantic code. Trees are never touched --
Generate only rewrites the plot's derived terrain / understory outputs, which
are regenerable by definition.

    python scene_editor.py plots                  # run directory of plots
    python scene_editor.py plots/plot_11574.csv   # a single plot table

Workflow mirrors the parameter editor: edit -> Generate -> look -> repeat.
"Export settings" writes a `plot_<UID>.env.json` sidecar in the table's
directory using the same schema as `replicate.py`, so a good environment can
be re-applied by future pipeline runs of that directory.

Controls: drag to rotate, scroll to zoom, R (or the Reset view button)
restores the camera, Left/Right keys page through the plots in the directory.

Requires: numpy, pandas, matplotlib, tkinter (stdlib), laspy (via las.py)
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import subprocess
import sys
import time

import numpy as np
import pandas as pd

import tkinter as tk
from tkinter import filedialog, messagebox, ttk

import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402

import layout
import tree_list
import schema
from las import read_canonical_las

HERE = os.path.dirname(os.path.abspath(__file__))
PY = sys.executable

# semantic code -> display color / name (README_schema)
SEM_COLORS = {
    schema.SEM_TREE: "#8B5A2B",   # trunk stem
    schema.SEM_BRANCH: "#B08968",  # branch stem
    schema.SEM_FOLIAGE: "#3A8C4A",  # foliage
    schema.SEM_SHRUB: "#6B8E23",   # shrub
    schema.SEM_GRASS: "#C9B037",   # grass
    schema.SEM_GROUND: "#8A8A8A",  # ground
}
SEM_NAMES = {schema.SEM_TREE: "trunk", schema.SEM_BRANCH: "branch stem",
             schema.SEM_FOLIAGE: "foliage", schema.SEM_SHRUB: "shrub",
             schema.SEM_GRASS: "grass", schema.SEM_GROUND: "ground"}

DISPLAY_CAP = 120000        # scatter subsample cap; full clouds reach millions

# control defaults, in panel order. Keys mirror the stage CLIs.
DEFAULTS = {
    # terrain
    "dtm_resolution": "1.0", "gen_bumpiness": "0.5", "gen_relief": "10.0",
    "ground_density": "1.0", "ground_noise": "0.0",
    # understory: shrubs
    "shrub_cover": "0.15", "shrub_cell_size": "2.0",
    "shrub_width_min": "0.8", "shrub_width_max": "1.6",
    "shrub_height_min": "0.4", "shrub_height_max": "1.8",
    "points_per_shrub": "300", "shrub_margin": "0.5",
    # understory: grass
    "grass_cover": "0.6", "grass_height_min": "0.03", "grass_height_max": "0.35",
    "grass_patch_scale": "5.0", "grass_candidate_density": "200",
    "canopy_attenuation": "0.3",
    # shared
    "seed": "0",
}


# ---------------------------------------------------------------------------
# plot discovery / loading (headless-testable)
# ---------------------------------------------------------------------------

def find_plots(path):
    """Plots under a run directory or a single table path.

    Returns sorted [(plot_uid, csv_path, plot_dir, restructured)] where
    `restructured` means the canonical trimmed cloud exists (terrain reads
    from the canonical layout, so unrestructured plots cannot be edited).
    """
    if os.path.isdir(path):
        csvs = tree_list.find_plots(path)
    elif os.path.isfile(path):
        csvs = [path]
    else:
        raise FileNotFoundError(path)
    if not csvs:
        raise RuntimeError(f"no plot tables (plot_*.json / plot_*.csv) "
                           f"found under {path}")

    out = []
    for csv in csvs:
        uid = os.path.splitext(os.path.basename(csv))[0]
        plot_dir = os.path.join(os.path.dirname(os.path.abspath(csv)), uid)
        trimmed = layout.plot_path(plot_dir, uid, "trimmed_las", create=False)
        out.append((uid, csv, plot_dir, os.path.exists(trimmed)))
    return out


def load_plot_cloud(plot_dir, plot_uid, which):
    """(xyz, semantic_codes) of a plot-level canonical cloud.

    `which` is a layout key: 'trimmed_las', 'terrain_las' or 'understory_las'.
    """
    path = layout.plot_path(plot_dir, plot_uid, which, create=False)
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    xyz, cols = read_canonical_las(path, level="plot")
    return xyz, np.asarray(cols["semantic_code"])


# ---------------------------------------------------------------------------
# stage invocation (subprocesses, same as run_pipeline)
# ---------------------------------------------------------------------------

def terrain_args(values, dtm_mode, dtm_path):
    """CLI args (list) for terrain.py from the GUI's values dict."""
    args = ["--source", values.get("terrain_source", "trimmed"),
            "--seed", values["seed"],
            "--density", values["ground_density"],
            "--noise", values["ground_noise"]]
    if dtm_mode == "provided":
        if not dtm_path:
            raise ValueError("provided-DTM mode needs a .asc path")
        return args + ["--dtm", dtm_path]
    return args + ["--generate-dtm",
                   "--dtm-resolution", values["dtm_resolution"],
                   "--gen-bumpiness", values["gen_bumpiness"],
                   "--gen-relief", values["gen_relief"]]


def understory_args(values, dtm_mode, dtm_path, no_shrub_outputs):
    args = ["--seed", values["seed"],
            "--shrub-cover", values["shrub_cover"],
            "--shrub-cell-size", values["shrub_cell_size"],
            "--shrub-width-min", values["shrub_width_min"],
            "--shrub-width-max", values["shrub_width_max"],
            "--shrub-height-min", values["shrub_height_min"],
            "--shrub-height-max", values["shrub_height_max"],
            "--points-per-shrub", values["points_per_shrub"],
            "--shrub-margin", values["shrub_margin"],
            "--grass-cover", values["grass_cover"],
            "--grass-height-min", values["grass_height_min"],
            "--grass-height-max", values["grass_height_max"],
            "--grass-patch-scale", values["grass_patch_scale"],
            "--grass-candidate-density", values["grass_candidate_density"],
            "--canopy-attenuation", values["canopy_attenuation"]]
    if dtm_mode == "provided" and dtm_path:
        # a provided DTM is not stored in the plot, so understory needs it too
        args += ["--dtm", dtm_path]
    if no_shrub_outputs:
        args.append("--no-shrub-outputs")
    return args


def run_stage(script, args, log):
    """Run one stage script; append its output to `log`. True on success."""
    cmd = [PY, os.path.join(HERE, script), *args]
    log(f"$ {' '.join(os.path.basename(c) if i == 1 else c for i, c in enumerate(cmd))}")
    proc = subprocess.run(cmd, capture_output=True, text=True)
    for line in (proc.stdout or "").strip().splitlines():
        log(f"  {line}")
    if proc.returncode != 0:
        err = (proc.stderr or "").strip().splitlines()
        log(f"  ERROR: {err[-1] if err else 'exit %d' % proc.returncode}")
        return False
    return True


# ---------------------------------------------------------------------------
# env-sidecar export (schema shared with replicate.py)
# ---------------------------------------------------------------------------

ENV_KEYS = ["seed", "dtm_resolution", "gen_bumpiness", "gen_relief",
            "ground_density", "ground_noise",
            "shrub_cover", "shrub_cell_size", "shrub_width_min",
            "shrub_width_max", "grass_cover", "grass_candidate_density",
            "canopy_attenuation"]


def env_from_values(values, dtm_mode, dtm_path):
    """Flat env dict in replicate.py's sidecar schema."""
    env = {}
    for k in ENV_KEYS:
        v = float(values[k])
        env[k] = int(v) if k == "seed" else v
    if dtm_mode == "provided" and dtm_path:
        env["dtm"] = str(dtm_path)
    return env


# ---------------------------------------------------------------------------
# GUI
# ---------------------------------------------------------------------------

class SceneEditor(tk.Tk):
    """Settings panel (left) + semantic-colored 3D viewer (right)."""

    def __init__(self, path):
        super().__init__()
        self.plots = find_plots(path)
        self.index = 0
        self.zoom = 1.0
        self.span = None
        self.values = {}              # StringVar per DEFAULTS key
        self.dtm_mode = tk.StringVar(value="generate")
        self.dtm_path = tk.StringVar(value="")
        self.terrain_source = tk.StringVar(value="trimmed")
        self.no_shrub_outputs = tk.BooleanVar(value=False)
        self.cloud_choice = tk.StringVar(value="understory_las")
        # per-class visibility toggles (semantic codes), all on by default
        self.class_on = {code: tk.BooleanVar(value=True)
                         for code in SEM_COLORS}
        self._cloud_cache = None      # (key, xyz, sem) of the last load
        self.title(self._title_text())
        self._build()
        self._select_plot(0)

    @staticmethod
    def _treesim_stamp():
        return time.strftime("%Y-%m-%d %H:%M",
                             time.localtime(os.path.getmtime(
                                 os.path.join(HERE, "terrain.py"))))

    def _title_text(self):
        n_ok = sum(1 for p in self.plots if p[3])
        return (f"PlausibiliTree scene editor [{n_ok}/{len(self.plots)} plots ready] "
                f"[terrain.py {self._treesim_stamp()}]")

    # -- construction ------------------------------------------------------

    def _build(self):
        toolbar = ttk.Frame(self)
        toolbar.pack(side=tk.TOP, fill=tk.X, padx=6, pady=4)
        ttk.Button(toolbar, text="Open...",
                   command=self._open).pack(side=tk.LEFT)
        ttk.Label(toolbar, text="plot:").pack(side=tk.LEFT, padx=(14, 2))
        self.plot_combo = ttk.Combobox(toolbar, state="readonly", width=22,
                                       values=[p[0] for p in self.plots])
        self.plot_combo.pack(side=tk.LEFT)
        self.plot_combo.bind("<<ComboboxSelected>>",
                             lambda e: self._select_plot(self.plot_combo.current()))
        self.status = ttk.Label(toolbar, text="")
        self.status.pack(side=tk.RIGHT)

        body = ttk.Frame(self)
        body.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        self._build_controls(body)
        self._build_viewer(body)

        self.bind("<Left>", lambda e: self._step(-1))
        self.bind("<Right>", lambda e: self._step(+1))
        self.bind("r", lambda e: self._reset_view())
        self.bind("R", lambda e: self._reset_view())

    def _field(self, parent, key, label, width=7):
        ttk.Label(parent, text=label).pack(side=tk.LEFT, padx=(8, 2))
        var = tk.StringVar(value=DEFAULTS[key])
        ttk.Entry(parent, width=width, textvariable=var).pack(side=tk.LEFT)
        self.values[key] = var

    def _build_controls(self, parent):
        left = ttk.Frame(parent)
        left.pack(side=tk.LEFT, fill=tk.Y, padx=(6, 6), pady=4)

        # ---- terrain panel ----
        ter = ttk.LabelFrame(left, text="Topography (terrain.py)")
        ter.pack(side=tk.TOP, fill=tk.X, pady=(0, 6))
        row = ttk.Frame(ter); row.pack(fill=tk.X, pady=2)
        ttk.Radiobutton(row, text="generate DTM", variable=self.dtm_mode,
                        value="generate").pack(side=tk.LEFT)
        ttk.Radiobutton(row, text="provided .asc", variable=self.dtm_mode,
                        value="provided").pack(side=tk.LEFT, padx=(10, 0))
        row = ttk.Frame(ter); row.pack(fill=tk.X, pady=2)
        self._field(row, "dtm_resolution", "cell (m)")
        self._field(row, "gen_bumpiness", "bump")
        self._field(row, "gen_relief", "relief (m)")
        row = ttk.Frame(ter); row.pack(fill=tk.X, pady=2)
        ttk.Label(row, text=".asc:").pack(side=tk.LEFT, padx=(8, 2))
        ttk.Entry(row, textvariable=self.dtm_path, width=22).pack(side=tk.LEFT)
        ttk.Button(row, text="...", width=3,
                   command=self._browse_dtm).pack(side=tk.LEFT, padx=2)
        row = ttk.Frame(ter); row.pack(fill=tk.X, pady=2)
        self._field(row, "ground_density", "ground pts/cell")
        self._field(row, "ground_noise", "noise (m)")
        row = ttk.Frame(ter); row.pack(fill=tk.X, pady=2)
        ttk.Label(row, text="source:").pack(side=tk.LEFT, padx=(8, 2))
        ttk.Combobox(row, textvariable=self.terrain_source, width=8,
                     state="readonly", values=["trimmed", "raw"]).pack(side=tk.LEFT)

        # ---- understory panel ----
        und = ttk.LabelFrame(left, text="Shrubs + grass (understory.py)")
        und.pack(side=tk.TOP, fill=tk.X, pady=(0, 6))
        row = ttk.Frame(und); row.pack(fill=tk.X, pady=2)
        self._field(row, "shrub_cover", "shrub cover")
        self._field(row, "shrub_cell_size", "cell (m)")
        row = ttk.Frame(und); row.pack(fill=tk.X, pady=2)
        self._field(row, "shrub_width_min", "width min")
        self._field(row, "shrub_width_max", "max")
        row = ttk.Frame(und); row.pack(fill=tk.X, pady=2)
        self._field(row, "shrub_height_min", "height min")
        self._field(row, "shrub_height_max", "max")
        row = ttk.Frame(und); row.pack(fill=tk.X, pady=2)
        self._field(row, "points_per_shrub", "pts/shrub")
        self._field(row, "shrub_margin", "margin")
        row = ttk.Frame(und); row.pack(fill=tk.X, pady=2)
        self._field(row, "grass_cover", "grass cover")
        self._field(row, "grass_patch_scale", "patch (m)")
        row = ttk.Frame(und); row.pack(fill=tk.X, pady=2)
        self._field(row, "grass_height_min", "height min")
        self._field(row, "grass_height_max", "max")
        row = ttk.Frame(und); row.pack(fill=tk.X, pady=2)
        self._field(row, "grass_candidate_density", "cand/m2")
        self._field(row, "canopy_attenuation", "attenuation")
        row = ttk.Frame(und); row.pack(fill=tk.X, pady=2)
        self._field(row, "seed", "seed")
        ttk.Checkbutton(row, text="skip per-shrub outputs",
                        variable=self.no_shrub_outputs).pack(side=tk.LEFT,
                                                             padx=(10, 0))

        # ---- actions ----
        act = ttk.Frame(left)
        act.pack(side=tk.TOP, fill=tk.X)
        self.gen_btn = ttk.Button(act, text="Generate", command=self._generate)
        self.gen_btn.pack(side=tk.LEFT, fill=tk.X, expand=True)
        ttk.Button(act, text="Export settings",
                  command=self._export_env).pack(side=tk.LEFT, padx=(6, 0))

        note = ttk.Label(left, foreground="#555", justify=tk.LEFT, wraplength=280,
                         text="Generate re-runs terrain + understory on the "
                              "selected plot, overwriting its terrain / "
                              "understory outputs. Tree outputs are never "
                              "modified.")
        note.pack(side=tk.TOP, fill=tk.X, pady=(6, 0))

    def _build_viewer(self, parent):
        right = ttk.Frame(parent)
        right.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, pady=4)

        top = ttk.Frame(right)
        top.pack(side=tk.TOP, fill=tk.X)
        ttk.Label(top, text="cloud:").pack(side=tk.LEFT)
        self.cloud_combo = ttk.Combobox(top, textvariable=self.cloud_choice,
                                        state="readonly", width=12,
                                        values=["understory_las", "terrain_las",
                                                "trimmed_las"])
        self.cloud_combo.pack(side=tk.LEFT)
        self.cloud_combo.bind("<<ComboboxSelected>>",
                              lambda e: self._show_cloud())

        # per-class show/hide toggles; a hidden class is left out of the
        # scatter and of the display subsample, so the rest stay dense
        layers = ttk.Frame(right)
        layers.pack(side=tk.TOP, fill=tk.X)
        ttk.Label(layers, text="layers:").pack(side=tk.LEFT)
        for code, name in SEM_NAMES.items():
            tk.Checkbutton(layers, text=name, fg=SEM_COLORS[code],
                           variable=self.class_on[code],
                           command=self._show_cloud).pack(side=tk.LEFT,
                                                          padx=(6, 0))

        self.fig = Figure(figsize=(6.0, 5.6), dpi=100)
        self.ax = self.fig.add_subplot(111, projection="3d")
        self.home_view = (self.ax.elev, self.ax.azim)
        self.canvas = FigureCanvasTkAgg(self.fig, master=right)
        self.canvas.mpl_connect("scroll_event", self._on_scroll)
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        bar = ttk.Frame(right)
        bar.pack(side=tk.TOP, fill=tk.X)
        ttk.Button(bar, text="Reset view",
                   command=self._reset_view).pack(side=tk.LEFT)
        self.info = ttk.Label(bar, text="no cloud loaded")
        self.info.pack(side=tk.TOP)

        self.log_text = tk.Text(right, height=7, wrap="none",
                                font=("TkFixedFont", 8))
        self.log_text.pack(side=tk.TOP, fill=tk.X, pady=(4, 0))
        self.log_text.configure(state="disabled")

    def _browse_dtm(self):
        p = filedialog.askopenfilename(filetypes=[("ASC grid", "*.asc"),
                                                  ("All files", "*.*")])
        if p:
            self.dtm_path.set(p)
            self.dtm_mode.set("provided")

    # -- helpers -----------------------------------------------------------

    def _log(self, msg):
        self.log_text.configure(state="normal")
        self.log_text.insert("end", msg + "\n")
        self.log_text.see("end")
        self.log_text.configure(state="disabled")
        self.update_idletasks()

    def _select_plot(self, i):
        self.index = i
        uid, csv, plot_dir, ready = self.plots[i]
        self.plot_combo.set(uid)
        if not ready:
            self.status.configure(text=f"{uid}: not restructured -- run the "
                                      f"pipeline through restructure first")
            return
        self.status.configure(text=f"selected: {uid}")
        self._show_cloud(refresh=False)

    def _step(self, d):
        if len(self.plots) > 1:
            self._select_plot((self.index + d) % len(self.plots))

    def _values(self):
        return {k: v.get().strip() for k, v in self.values.items()}

    # -- viewer ------------------------------------------------------------

    def _show_cloud(self, refresh=True):
        uid, csv, plot_dir, ready = self.plots[self.index]
        which = self.cloud_choice.get()
        key = (plot_dir, uid, which)
        if self._cloud_cache is not None and self._cloud_cache[0] == key:
            xyz, sem = self._cloud_cache[1], self._cloud_cache[2]
        else:
            try:
                xyz, sem = load_plot_cloud(plot_dir, uid, which)
            except FileNotFoundError as e:
                self.status.configure(
                    text=f"{uid}: {os.path.basename(str(e))} not generated yet "
                         f"(run Generate)")
                return
            self._cloud_cache = (key, xyz, sem)

        visible = [c for c, on in self.class_on.items() if on.get()]
        vis_mask = np.isin(sem, visible) if visible else np.zeros(len(sem), bool)
        n_vis = int(vis_mask.sum())

        ax = self.ax
        ax.clear()
        if n_vis:
            # seeded display subsample over the VISIBLE classes only, so
            # hiding a dominant class (e.g. foliage) keeps the rest dense
            rng = np.random.default_rng(0)
            if n_vis > DISPLAY_CAP:
                sel = rng.choice(n_vis, DISPLAY_CAP, replace=False)
            else:
                sel = np.arange(n_vis)
            vis_xyz = xyz[vis_mask][sel]
            vis_sem = sem[vis_mask][sel]
            shown = len(sel)
            for code, color in SEM_COLORS.items():
                if code not in visible:
                    continue
                m = vis_sem == code
                if not m.any():
                    continue
                ax.scatter(vis_xyz[m, 0], vis_xyz[m, 1], vis_xyz[m, 2], s=2,
                           c=color, depthshade=True,
                           label=f"{SEM_NAMES[code]} "
                                 f"({int((sem == code).sum()):,})")
        else:
            shown = 0
        lo, hi = xyz.min(axis=0), xyz.max(axis=0)
        self.span = tuple(float(s) for s in np.maximum(hi - lo, 1e-6))
        ax.set_box_aspect(self.span, zoom=self.zoom)
        ax.set_xlabel("x (m)"); ax.set_ylabel("y (m)"); ax.set_zlabel("z (m)")
        if n_vis:
            ax.legend(loc="upper left", fontsize=8)
        name = {"understory_las": "understory", "terrain_las": "terrain",
                "trimmed_las": "trees (trimmed)"}[which]
        suffix = "" if len(visible) == len(SEM_COLORS) else \
            f" ({len(SEM_COLORS) - len(visible)} class(es) hidden)"
        ax.set_title(f"{uid} - {name} - {shown:,} shown of {n_vis:,} visible"
                     f" of {len(xyz):,} points{suffix}")
        self.canvas.draw_idle()

    def _on_scroll(self, event):
        if event.inaxes is not self.ax or self.span is None:
            return
        self.zoom = float(np.clip(self.zoom * 1.15 ** event.step, 0.3, 12.0))
        self.ax.set_box_aspect(self.span, zoom=self.zoom)
        self.canvas.draw_idle()

    def _reset_view(self):
        self.zoom = 1.0
        self.ax.view_init(*self.home_view)
        if self.plots[self.index][3]:
            self._show_cloud(refresh=False)

    # -- actions -----------------------------------------------------------

    def _generate(self):
        uid, csv, plot_dir, ready = self.plots[self.index]
        if not ready:
            messagebox.showwarning(
                "plot not restructured",
                f"{uid} has no canonical trimmed cloud. Run the pipeline "
                f"(through restructure) on it first; terrain reads from the "
                f"canonical layout.")
            return
        values = self._values()
        values["terrain_source"] = self.terrain_source.get()
        dtm_mode = self.dtm_mode.get()
        dtm_path = self.dtm_path.get().strip()
        try:
            t_args = terrain_args(values, dtm_mode, dtm_path)
            u_args = understory_args(values, dtm_mode, dtm_path,
                                     self.no_shrub_outputs.get())
        except ValueError as e:
            messagebox.showerror("invalid settings", str(e))
            return

        self.gen_btn.configure(state="disabled", text="Generating...")
        self._log(f"--- {uid} ---")
        self.update_idletasks()
        try:
            ok = run_stage("terrain.py", [csv, *t_args], self._log)
            if ok:
                ok = run_stage("understory.py", [csv, *u_args], self._log)
            if ok:
                self._log(f"{uid}: done")
                self.cloud_choice.set("understory_las")
                self._cloud_cache = None      # the clouds changed on disk
                self._show_cloud(refresh=False)
                self.status.configure(text=f"{uid}: environment regenerated")
        finally:
            self.gen_btn.configure(state="normal", text="Generate")

    def _export_env(self):
        uid, csv, plot_dir, ready = self.plots[self.index]
        values = self._values()
        env = env_from_values(values, self.dtm_mode.get(),
                              self.dtm_path.get().strip())
        out = os.path.join(os.path.dirname(os.path.abspath(csv)),
                           f"{uid}.env.json")
        with open(out, "w") as f:
            json.dump({"plot_uid": uid, "source": os.path.abspath(csv),
                       "env": env}, f, indent=2)
        self._log(f"exported settings -> {out} "
                  f"(run_pipeline.py applies this sidecar on future runs)")
        self.status.configure(text=f"settings exported: {os.path.basename(out)}")

    # -- open --------------------------------------------------------------

    def _open(self):
        p = filedialog.askopenfilename(title="plot table csv",
                                       filetypes=[("CSV", "*.csv")])
        if not p:
            # fall back to a directory picker
            p = filedialog.askdirectory(title="...or pick a run directory")
        if not p:
            return
        try:
            self.plots = find_plots(p)
        except (RuntimeError, FileNotFoundError) as e:
            messagebox.showerror("open", str(e))
            return
        self.index = 0
        self.title(self._title_text())
        self.plot_combo["values"] = [q[0] for q in self.plots]
        self._select_plot(0)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Interactively edit terrain / shrub / grass settings for "
                    "already-generated plots and view the result")
    ap.add_argument("path", nargs="?", default=None,
                    help="run directory of plot_*.csv tables, or one csv")
    args = ap.parse_args(argv)

    path = args.path
    if path is None:
        root = tk.Tk()
        root.withdraw()
        p = filedialog.askdirectory(title="run directory of plots")
        root.destroy()
        if not p:
            return 0
        path = p
    app = SceneEditor(path)
    app.mainloop()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
