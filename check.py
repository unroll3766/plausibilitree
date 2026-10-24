"""
check.py: fast static checks for the synthesis pipeline, so breakage is caught
here and not twenty minutes into a dataset run.

    python check.py

What it checks:

  1. **Syntax**, by compiling every script in this directory.
  2. **Undefined names**: every name used but never defined or imported, per
     file -- the kind of missing-import error that only surfaces at runtime,
     deep in a batch run.

The deep-learning project (training, tuning, diagnostics) is a separate
project now and carries its own checker.

Run the full pipeline on a small plot after structural changes:
    python run_pipeline.py <dir-with-plot-tables> --radius 15 --seed 42
"""

from __future__ import annotations

import ast
import builtins
import glob
import os
import sys


def undefined_names(path):
    """Names loaded but never bound anywhere in the file."""
    tree = ast.parse(open(path).read(), path)
    bound = set(dir(builtins)) | {"__name__", "__file__", "__doc__", "__spec__"}
    for n in ast.walk(tree):
        if isinstance(n, (ast.Import, ast.ImportFrom)):
            for a in n.names:
                bound.add((a.asname or a.name).split(".")[0])
        elif isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            bound.add(n.name)
        elif isinstance(n, ast.arg):
            bound.add(n.arg)
        elif isinstance(n, ast.Name) and isinstance(n.ctx, (ast.Store, ast.Del)):
            bound.add(n.id)
        elif isinstance(n, ast.ExceptHandler) and n.name:
            bound.add(n.name)
        elif isinstance(n, ast.Global):
            bound.update(n.names)
    used = {n.id for n in ast.walk(tree)
            if isinstance(n, ast.Name) and isinstance(n.ctx, ast.Load)}
    return sorted(used - bound)


def main():
    here = os.path.dirname(os.path.abspath(__file__)) or "."
    files = sorted(glob.glob(os.path.join(here, "*.py")))
    failures = []

    print("syntax")
    for f in files:
        try:
            compile(open(f).read(), f, "exec")
        except SyntaxError as e:
            failures.append(f"{f}: syntax error line {e.lineno}: {e.msg}")
            print(f"  FAIL {os.path.basename(f)}: line {e.lineno}: {e.msg}")
    print(f"  {len(files)} file(s) parsed")

    print("undefined names")
    found = False
    for f in files:
        missing = undefined_names(f)
        if missing:
            found = True
            failures.append(f"{f}: undefined {missing}")
            print(f"  FAIL {os.path.basename(f)}: {missing}")
    if not found:
        print("  none")

    if failures:
        print(f"\n{len(failures)} failure(s)")
        return 1
    print("\nall checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
