# NOTICE

This repository contains modified third-party code.

## Simple Synth Tree (`treesim.py`, `generate_example_tree.py`)

The tree simulator `treesim.py`, and the demo entry point
`generate_example_tree.py`, are modified versions of code from **Simple
Synth Tree**, the library accompanying:

> M. Bryson, F. Wang and J. Allworth, "Using synthetic tree data in deep
> learning-based individual tree segmentation using LiDAR point clouds",
> *Remote Sensing* 2023, 15(9), 2380.
> https://www.mdpi.com/2072-4292/15/9/2380

The original code was authored by Mitch Bryson. Local modifications include
cascading stem splits with leader-cascade targeting, slice-proportional
foliage point allocation, skeleton extraction, seeded deterministic sampling,
and the branching/geometry changes documented in
[README_generation.md](README_generation.md) and
[README_assembly.md](README_assembly.md). Upstream behaviour that must stay
comparable (the stem spline and its published-equation basis) is documented
in README_assembly.md.

The copies this repository was built from carried no license notice beyond
the author attribution. The copyright holder of this repository should
confirm redistribution rights with the upstream authors (or their chosen
license) before publishing, or replace these files with a clean-room
implementation.

Everything else in this repository is original work under the repository
LICENSE.
