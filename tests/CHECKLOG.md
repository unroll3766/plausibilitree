# Check Log

Running log of unique checks, per the global agent rules. One to three
lines per entry: date, what was verified, command/method, actual result.
Re-entries only when behavior changes. Long-form narrative lives in
`VERIFICATION_LOG.md`; this file is the short pool.

- 2026-11-07 — FIA tree-list CSV/JSON equivalence (W17): `fia_helper.py`
  twice on CA_TREE_SMALL.csv (seed 42, r 20) in both formats, then
  `assert_frame_equal` on loaded trees + plot fields via report sidecar —
  exit 0, 111 trees, geometry+center identical; CSV column order differs
  (Z canonical position), content identical.
- 2026-11-07 — CSV float round-trip must be exact: pandas' default fast
  parser was 1 ULP off on 17-sig-digit values (DBH 46 in -> 1.1684 vs
  1.1683999999999999); loader fixed with `float_precision="round_trip"`,
  proven by exact `==` of all 111 `tree_requests` from CSV vs JSON — pass.
- 2026-11-07 — W17 regression run: `conda run -n ssf python
  tests/run_all.py` — exit 0, all suites passed.
