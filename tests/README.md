# Regression suite

Promoted from the session-ephemeral `/tmp` verification suites referenced
throughout `VERIFICATION_LOG.md` (W7 round, W7 review and the whole-repo
review). Same invariants, durable home.

Run from the repo root, through the project's conda env (`ssf`):

```
conda run -n ssf python tests/run_all.py              # everything (~60 s)
conda run -n ssf python tests/test_placement_reflection.py   # fast (~25 s)
conda run -n ssf python tests/test_cross_artifact_invariants.py  # full fixture (~40 s)
```

Every suite exits 0 only if every check passes (per the verification
log's exit-code discipline: check `$?`, never eyeball `tail` output).
`--keep` leaves the fixture directory in place for debugging.

## What each suite guards

- **`test_placement_reflection.py`** — the placement draws and the
  generator's coordinate frames: raw skeletons anchored at their global
  base (the regression guard for the `place()` trunk-polyline bug found
  by the W7 review), same-seed determinism, scene rotation/reflection
  records replaying exactly onto cloud and skeleton (1 mm LAS
  re-quantization and the 4 dp angle record are the documented
  tolerances), the ~50 % mirror draw distribution, and chirality — the
  mirror is a true handedness flip (signed area over ordered branch
  bases) rather than a substituted rotation, and the W19 rotation grid —
  every recorded rotation is a multiple of ROTATION_STEP (default 30
  degrees) from [0, 360), all orientations occurring across seeds.
- **`test_cross_artifact_invariants.py`** — recomputes every stage's
  key quantities from a *different* artifact than the one that wrote
  it: skeleton-vs-cloud stem coax, info-JSON bbox/DBH/frustum volumes,
  spec 3.3 draping rules (lowest point on the DTM at the tree's
  (X, Y); per-point shrub drape), scene-record coverage, and a full
  scene-JSON link audit. New in the W15 round: section J re-derives the
  notes.md element graph from the polyline skeleton independently of the
  writer (single null-SOURCE root, SOURCE/DESTINATION reciprocity,
  stem/branch partition, lengths recomputed from polylines and
  base-tip pairs, taper-profile radii, cylinder volumes), and checks
  the trimmed variant carries the identical, place()-invariant graph.
  This is the pattern that caught four
  long-hidden bugs (see the verification log), so the fixtures are
  deliberately *not* flat: a generated DTM (rugged terrain) plus
  understory exercises the paths flat runs never touch.

## Conventions

- `fixture.py` builds everything via the pipeline's programmatic entry
  points (`fia_helper.fia_to_tree_list`, `generate_trees.generate_plot`,
  `run_pipeline.run_pipeline`) — no CLI subprocesses, no dependence on
  files outside the repo and stdlib/tmp.
- Fixtures are deterministic (seed 42, radius 20 m, 6 trees at 64
  pts/m) so failures are reproducible; md5s recorded in
  `VERIFICATION_LOG.md` cross-reference the same fixtures.
- New rounds should add checks here in the same style as the log
  entries: one invariant per check, a `PASS`/`FAIL` line each, a
  message explaining any tolerance.
