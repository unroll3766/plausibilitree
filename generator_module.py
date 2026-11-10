"""
generator_module.py: the generator module contract (design spec Section 7.1).

A generator module is the pluggable tree simulator behind the core
pipeline. The core builds every tree through this contract and nothing
else; simulator specifics (parameter formats, skeleton formats, draw
mechanics) live behind the implementing module (Section 7 quarantines
them there).

Contract:

    module.load_param_dict(path)   -> tree code -> parameter dictionary
    module.resolve_tree_code(code) -> the dictionary key for a tree code
    module.entry_for(code, dbh, ht) -> the resolved dictionary entry, gated
        by the entry's applicability bounds when sizes are given (7.3)
    module.build(rng, DBH, HT, tree_code, params) -> TreeModel
    module.place(model, x, y)      -> TreeModel positioned in plot coords

TreeModel holds:

    points      (N, 4) x, y, z in metres and one semantic label per point,
                origin at the tree base (z = 0 at the base; place() moves
                the model to its plot position)
    labels      (N,) semantic codes; at build time every point carries its
                coarse class (trunk stem or foliage by schema codes) --
                the trunk/branch split is resolved later at assembly
    skeleton    the module's documented skeleton format (for SST see
                Section 7.4); consumed by the trimming and volume stages,
                which therefore follow the configured module's format
    realized    realized dimensions from the skeleton (height, base
                radius, fork, ...)
    provenance  everything needed to reproduce/record the build
                (model_params, stem fraction, point budget, ...)

Reproducibility: the core seeds each tree from a stable hash of the plot
seed + TREE_UID, hands the module a fresh numpy Generator, and retries a
failing build with a rehashed seed (MAX_ATTEMPTS, spec Section 9). The
module must draw all of its randomness from that rng.

Module selection: generate_trees.py (-g/--generator, default "sst") and
run_pipeline.py select from the REGISTRY below.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass
class TreeModel:
    """One built tree, base-origin, plus its module-format skeleton."""
    points: np.ndarray              # (N, 4): x, y, z, semantic label
    labels: np.ndarray              # (N,) semantic codes (schema)
    skeleton: dict                  # module-documented format
    realized: dict = field(default_factory=dict)
    provenance: dict = field(default_factory=dict)


class GeneratorModule:
    """Abstract base class for generator modules (spec Section 7.1)."""

    name = "abstract"

    def load_param_dict(self, path: str | None = None):
        """Load the tree-code -> parameters dictionary (Section 8.6)."""
        raise NotImplementedError

    def resolve_tree_code(self, code: str) -> str:
        """A tree code -> the dictionary key (resolution order is the
        module's, documented in its section of the spec)."""
        raise NotImplementedError

    def entry_for(self, tree_code: str, dbh: float | None = None,
                  ht: float | None = None) -> dict:
        """The resolved dictionary entry for a tree code (Section 7.3).

        With sizes given, a candidate entry only applies when the tree's
        DBH/HT lie inside its applicability bounds (resolution order and
        DEFAULT fallback are the module's); the core passes sizes so the
        entry it records for a tree is the one the build will use."""
        raise NotImplementedError

    def build(self, rng: np.random.Generator, DBH: float, HT: float,
              tree_code: str, params: dict) -> TreeModel:
        """Build one tree. `params` is the resolved dictionary entry
        (often with the core's point-budget settings injected under
        'points_per_metre'); every stochastic draw MUST come from rng.
        Raises ValueError on a degenerate build so the core can retry
        with a rehashed seed."""
        raise NotImplementedError

    def place(self, model: TreeModel, x: float, y: float) -> TreeModel:
        """Position the base-origin model at plot coords (x, y): shift
        points and every positional field of the module's skeleton."""
        raise NotImplementedError


# ---------------------------------------------------------------------------
# Registry
# ---------------------------------------------------------------------------

REGISTRY: dict = {}


def register(cls):
    """Class decorator: add a generator module implementation."""
    REGISTRY[cls.name] = cls
    return cls


def get_generator(name: str) -> GeneratorModule:
    """Instantiate a registered generator module by name."""
    if name not in REGISTRY:
        # the SST module registers on import; keep the registry lazy
        import sst_generator                       # noqa: F401
    if name not in REGISTRY:
        raise KeyError(f"unknown generator module {name!r}; "
                       f"available: {sorted(REGISTRY)}")
    return REGISTRY[name]()
