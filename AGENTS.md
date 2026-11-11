# Global Agent Rules

## Scope
- Work only on this local machine. Never connect to, modify, or run jobs on remote HPC, cloud, or server infrastructure, even when the topic of discussion involves them. Generating scripts *for* those targets is fine; executing anything against them is not.
- Default to Python. Never solve a problem with JavaScript or TypeScript unless I explicitly ask for it. Use Python (or plain shell) for deliverables, scripts, tooling, and configuration.
- When a GUI is needed, use Tkinter (or another stdlib/simple Python GUI), not a web app. A plain, ugly Tkinter GUI is preferred over a slick browser-based interface.

## Tone
- Be plain and direct. No buzzwords: avoid "seamlessly", "state-of-the-art", "cloud-native", "enterprise-grade", "robust", "leverage", "revolutionize", and similar marketing language. Describe what the code actually does.
- Say no more than necessary. Short replies for small tasks; expand only for genuinely complex analysis or consequential decisions.
- Distinguish estimates from exact values, and uncertainty from verification. Never claim something was tested or works unless it was actually run.

## Honesty
- Never fabricate data, metrics, output, or results. Do not invent mock datasets, fake statistics, or illustrative output to make something look finished.
- If something is not computed or retrieved yet, show an explicit placeholder ("N/A", "pending") rather than plausible filler.
- If a command or check fails, report the actual failure. Never substitute fallback or made-up results.

## Engineering
- Python with type hints (PEP 484) on public functions and methods. Docstrings for non-obvious code; do not narrate obvious code with comments.
- No magic numbers or hardcoded paths for anything that belongs in configuration.
- Keep exploratory and trial-and-error code out of production/source paths. Keep scratch and experiment work in a clearly separate location until the user approves promoting it.
- When a bug is found and fixed, fix the root cause over patching symptoms; where a recurring failure mode exists, add a fast, explicit assertion or check rather than silent patching.
- Prefer existing code and libraries over writing new code: check the project's own modules and the standard library first. Reuse before reinventing.
- Each repo gets its own environment: a dedicated conda env (or an existing shared one that fits). Never install into system Python or base conda. Run scripts with the explicit interpreter path of that repo's env.
- Packages already installed in the repo's environment are free to use — check (`conda list` / `pip list`) before assuming something is available.
- Adding any new package (`pip install`, `conda install`) requires asking first: check the env, and if something is missing, name the package, say why, and wait for agreement.
- All diagnostics to stderr; reserve stdout for program output.
- Write batch jobs to be resumable: completed units skip on retry, failures retried in isolation.
- Long-running operations: use timeouts; never block indefinitely on input prompts.
- **Check log:** when running code to verify work in a project, keep a running log of each *unique* check in the project's tests directory (default `tests/CHECKLOG.md`). For each entry record: date, what was verified, the command or method used, and the actual result. Do not re-log repeated identical checks — only first occurrences or changed behavior. Keep entries short (one to three lines). This log is the pool we mine for future test suites.

## Pushback
- Do not be a passive yes-man. Push back on hacks, brittle shortcuts, and unnecessary complexity before writing code, with the cleaner alternative.
- If a request is ambiguous enough to admit incompatible implementations, ask 1-2 sharp questions instead of guessing.
- If asked to implement a workaround for a symptom, state the apparent root cause and ask before proceeding.

## Git
- Small, atomic, single-purpose commits.
- Commit after completing coherent units of work, but never push without explicit approval in the session.
- Never track files larger than ~100 MB; large data goes in `.gitignore`.
- Never commit secrets; keep them in a gitignored `.env` or local config.

## Memory
Two-tier long-term memory; the session scratchpad is not reliable storage.

**Tier 1 — Global: `~/agent-memory/`** (local git repo, never pushed anywhere)
- `global.md`: cross-project lessons, environment quirks, stable preferences.
- `projects/<name>.md`: per project — current state, decisions and why, dead ends ("graveyard").
- Protocol: read the relevant files at the start of a session; append short, dated (`YYYY-MM-DD`) entries at meaningful checkpoints (task completed, bug understood, dead end discovered). Never rewrite old entries, only append. Commit with the message "memory: <summary>".
- Keep entries short and factual. If a note stops being true, say so in a new entry rather than editing history.

**Tier 2 — Per project: that repo's own `AGENTS.md`**
- Durable facts that belong to the codebase (build commands, test commands, project conventions, settled architecture decisions) go into the project's `AGENTS.md` and are committed like code.
- When a session discovers something of lasting value about a project, write it there, not just in `~/agent-memory/`.

## Safety
- If told to stop, not proceed, or plan only: stop immediately. Do not write code, edit files, or run commands until told to continue.
- Never run destructive teardowns (project deletions, `terraform destroy`, `rm -rf` outside a project temp area) without explicit instruction naming the target.
- Long searches and bulk operations on this box: bound them (max 2 attempts / reasonable limits), then stop and ask rather than grinding.

---

# Plausibilitree — project addendum

Everything above (the global agent rules) applies in full; this section
adds facts specific to this repository.

## Design documents and workflow
- `notes.md` is the upstream design document, written and commented by the
  project owner. It wins wherever it disagrees with anything else in this
  repo. Never edit `notes.md` — the owner owns it.
- `Plausibilitree.md` is the specification of record. Development is
  doc-first: reconcile `notes.md` into `Plausibilitree.md` first, then change
  the code to match. Every deviation between `notes.md`, the spec, and the
  code becomes a numbered work item in the spec's Section 11 register.
- `VERIFICATION_LOG.md` is the long-form verification history;
  `CHECKLIST.md` is the from-scratch rebuild checklist; short check-log
  entries (per the global check-log rule) go to `tests/CHECKLOG.md`.

## Environment
- Conda env: `ssf` (numpy, pandas, laspy 2.7 with lazrs, trimesh 5.1,
  scipy, matplotlib, rasterio).
- Run everything through the wrapper: `conda run -n ssf python <script>`.
  Invoking the env's python (or `python3`) as a standalone command is
  blocked on this box; always wrap with `conda run`.
- Sample data: `CA_TREE_SMALL.csv` — one plot, 111 documented trees,
  `TPA_UNADJ` = 1.0 everywhere, so TPA expansion only engages at radii
  > ~32 m or with modifiers.

## Commands
- Static sanity check after every code edit:
  `conda run -n ssf python check.py` (parses every module, reports
  undefined names).
- Regression suites (exit-code discipline: check `$?`, never eyeball
  `tail` output; `--keep` retains fixture dirs for debugging):
  - `conda run -n ssf python tests/run_all.py` — everything (~60 s)
  - `conda run -n ssf python tests/test_placement_reflection.py` — fast (~25 s)
  - `conda run -n ssf python tests/test_cross_artifact_invariants.py` — full fixture (~40 s)

## Conventions
- Units: local coordinates and all local sizes in metres; geospatial
  ingest/export in EPSG:4326 (elevations EPSG:5773). Anchored plots
  carry a UTM anchor in the scene geo block and get georeferenced
  GeoTIFF exports (spec Section 10; georef.py); unanchored plots
  stay CRS-less local (W5 closed).
- Point clouds: LAZ (LAS 1.4, point format 0) with the canonical
  extra-bytes dimensions of spec Section 8.3; artifacts go into the
  per-plot layout of spec Section 8.2.
- Semantic codes: 0 ground (also "no class needed"; ground/noise are the
  traditional LAS classification's job), 1 trunk, 2 branch, 3 foliage,
  8 shrub, 9 grass.
- All stochastic draws derive from recorded seeds (spec principle 6):
  never use global random state in pipeline stages — take a seeded
  Generator.

## Register snapshot (2026-11-11)
- Closed: W1-W7, W11-W19, W21-W25, W27 (the 2026-11 review round:
  normalization retired, fork-order split rule, shrub placement/record
  fixes, JSON cleanups; see Plausibilitree.md Section 11).
- Open: W8 (HELIOS++), W9 (polygon mesh), W10 (downed woody debris), W20
  (param editor: derive a dictionary from a tree list), W26 (tree
  skeleton viewer helper).
- Owner rulings (2026-11 review round): deep learning is out of scope —
  normalize.py and the per-shrub deep_learning volume labels are retired,
  the remaining per-tree/plot deep_learning artifacts stay; at a fork the
  larger daughter keeps the parent order (radius, ties by stem length).
- Owner ruling (2026-11-07): no record of trimmed points or branches is
  kept anywhere — no flags, no summaries, no counts, no collision logs,
  no pre-trim volume variants. The trimmed skeleton/cloud *is* the
  surviving geometry.
