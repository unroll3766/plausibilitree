# Parameter editor (`param_editor.py`)

An interactive tool for tuning the tree parameters in `functional_groups.json`
(the parameter dictionary, spec Section 8.6; legacy CSVs open too)
(or any similar parameter CSV) with immediate visual feedback: edit any cell,
generate N sample trees from a row with SimpleSynthTree, page through them in
a 3D view, and save the edited table as a new CSV when you like the results.

```
python param_editor.py data/functional_groups.json
python param_editor.py            # no argument: opens a file picker
```

## The window

- **Left: the table, transposed.** Every CSV row and column, one editable
  cell each — but displayed as one row per *parameter* and one column per
  *group*, with the parameter names frozen in the leftmost column so they
  stay visible while the group columns scroll; both halves scroll
  vertically together (wheel over the table; **Shift+wheel** scrolls the
  group columns sideways). The CSV's first column (the group labels, e.g.
  `FUNC_GROUP`) becomes the editable header row of the group columns. Click
  any cell in a column to select that group (the active group shows in the
  status bar); the selected group is what Generate samples from. Saving and
  reading the table are unchanged by the transpose — the CSV keeps its
  original groups-as-rows orientation.
- **Right: controls + viewer.**
  - *N trees* (how many samples), *Seed*, *Points* per tree.
  - *Sample HT / DBH ranges* (m): what to draw heights and diameters from
    (functional-group tables carry no sizes — those normally come from
    inventory rows).
  - **Generate** regenerates the sample set from the current cell values —
    edit, press Generate, repeat.
  - The 3D view shows one tree at a time. Each point class can be shown or
    hidden: **trunk** (the main stem, brown), **branch stems** (forks,
    sienna) and **foliage** (green), via the checkboxes under the view or
    the **1/2/3** keys. Toggling never re-frames the camera — the view is
    always sized to the full cloud. Stem points are split trunk/branch by
    their nearest skeleton stem: polyline 0 is the trunk, the rest are
    forks, matching the pipeline's trunk/branch semantic split.
    **Left/Right keys** or the **Prev/Next buttons** page through the N
    samples; **drag** to rotate, **scroll** to zoom; **R** or the **Reset
    view button** restores the default camera. The plot title carries the
    label of the row that generated the current samples (e.g. `DF`), so
    they stay traceable even after you move the selection elsewhere. The
    line under the view shows the drawn height, DBH, crown width and point
    count of the current tree.
- **Toolbar:** Open... (load another dictionary; JSON or legacy CSV),
  Save As... (write the edited dictionary as JSON — or legacy CSV; nothing
  is written until you save).

## How a sample tree is drawn from a row

- Columns that map directly to treesim parameters are read from the row:
  `SPLIT_HT_LO/HI`, `SPLIT_PROB`, `SPLIT_NUM_LO/HI`, `NUM_BRANCH_LO/HI`,
  `MIN_CAN_HT_LO/HI`, `MAX_CAN_W_HT_LO/HI`, `TOP_DIST`, `MID_DIST`,
  `FP_RAND_HT`, `SPR_LO/HI`.
- The GUI samples HT and DBH from the sampling ranges and derives the rest
  exactly the way `fia_helper.py` does at pipeline time:
  - crown width `MAX_CAN_W` from `CW = exp(A0 + A1 * ln(DBH_in))` ft,
  - base diameter from the `TAPER` extrapolation
    `d(0) = DBH * (H / (H - 1.37)) ** (1/b)`.
- If the CSV carries its own `HT_LO/HI`, `BASE_DIA_LO/HI` or `MAX_CAN_W`
  columns (a hand-built, plot-style table), those columns win over the
  sampling fields.
- Missing columns fall back to treesim's built-in example parameters, so
  "similar" CSVs still preview.

## Parameter reference

What each column of a functional-group table controls, and how to read it.
"Fraction" values are multiplied by the drawn tree height. Heights and
diameters come from the inventory rows (or the editor's sampling fields); the
`A0/A1` and `TAPER` columns scale them, and everything else is consumed
directly by the simulator (`treesim.py`).

| Column(s) | What it is | Raising it... |
|---|---|---|
| `FUNC_GROUP` | The key `spcd_func_group.csv` maps species codes to; `DEFAULT` is the fallback when nothing matches | — |
| `DESCRIPTION` | Free text | — |
| `DBH_MIN`, `DBH_MAX` | Applicability bound (m) on the breast-height diameter of trees this entry may generate. A tree list DBH outside them falls through the resolution order (species -> functional group -> `DEFAULT`), so narrower bounds admit fewer trees to the entry | ...admits a wider size range to this entry |
| `HT_MIN`, `HT_MAX` | Same, for tree height (m) | ...admits a taller size range to this entry |
| `A0`, `A1` | Crown-width allometry (notes.md's a0/a1) `CW (m) = exp(A0 + A1 * ln(DBH_in))`, evaluated per tree from its DBH; the result becomes the tree's `MAX_CAN_W`. `A0` sets crown width at any given size, `A1` how strongly width grows with DBH | ...widens the crown at a given DBH (A0) / makes big trees disproportionately wide (A1); wider crowns also mean more collision trimming when plots are dense |
| `TAPER` | Butt-swell exponent of the taper function (notes.md: "taper function (base diameter)"): base diameter = `DBH * (H / (H - 1.37))^(1/b)`. The simulator's trunk tapers linearly from this base radius to zero at the tip | ...reduces butt swell (base diameter closer to DBH); lowering it thickens the trunk base |
| `SPLIT_PROB` | Probability the tree forks at all (0 = single straight stem) | ...more forked trees |
| `SPLIT_NUM_LO/HI` | How many cascading trunk splits a forked tree gets, drawn per tree. Each split targets the tallest remaining stem, so the forks climb the dominant leader and every new crown forms in the upper canopy, interlocked with its neighbours; subordinate forks also diverge at a damped rate. Default `1:1` = the original single fork; `1:5` gives hardwoods with multiple stacked leaders | ...more levels of forking and bushier multi-leader crowns; each leader's crown is proportionally airier at the same point count, so raise Points for denser multi-leader crowns |
| `SPLIT_HT_LO/HI` | Where each split happens, as a fraction of the stem being split (drawn in range) | ...pushes splits higher up the stems |
| `NUM_BRANCH_LO/HI` | Nominal branch count, drawn per tree and split between the lower crown band and the upper crown/fork band. The draw can land well below `HI` (the lower band especially can come out near-empty), so a lower `LO` shows up as a sparser, more open crown. The foliage point budget is allocated across 1 m height slices in proportion to the number of branches rooted in each slice (within a slice, by mesh area), so sections of crown carrying more branches are sampled proportionally more heavily | ...denser canopy and more foliage points |
| `MIN_CAN_HT_LO/HI` | Crown base: height of the lowest branch, as a fraction of tree height. Branches below this are zero-length | ...raises the crown base (the pruned-plantation look); more bare, visible stem |
| `MAX_CAN_W_HT_LO/HI` | Height (fraction) where the crown is widest: branch lengths ramp up from the crown base to here, then taper off toward the tip | ...moves the widest part of the crown upward |
| `TOP_DIST` | Maximum horizontal wander of the stem tip (m) - the lean of the whole tree | ...more leaning trees |
| `MID_DIST` | Horizontal wander of the mid-stem control point (m) - stem bow | ...more bowed / curved stems |
| `FP_RAND_HT` | Vertical jitter added to foliage points (m, upward only): softens the sharp geometric canopy surface | ...fuzzier, less clearly-bounded canopy |
| `SPR_LO/HI` | Stem-point fraction (drawn per tree): the share of the point budget sampled on the stem mesh rather than the foliage mesh | ...puts more of the point budget on the trunk, at the cost of canopy detail |

Useful combinations, in silhouette terms: a conifer look is low `SPLIT_PROB`,
low `MIN_CAN_HT` crown base variation, small `A0/A1` and high
`NUM_BRANCHES`; a hardwood look is high `SPLIT_PROB`, high `A0/A1`, a
low crown base, and `SPLIT_NUM_LO/HI` above `1:1` for cascading multi-leader
crowns. `TOP_DIST`/`MID_DIST` only change stem geometry, so they are
safe knobs for adding variety without changing the crown.

The editor's sampling fields complete the same chain: **Sample HT** draws the
tree height, **Sample DBH** draws the breast-height diameter that drives the
`A0/A1` allometry and the `TAPER` extrapolation. When the loaded CSV carries
`HT_LO/HI`, `BASE_DIA_LO/HI` or `MAX_CAN_W` columns of its own, those win over
the sampling fields.

## Reproducibility

A `(row, seed, sample index)` triple always generates the same tree: tree *i*
reseeds numpy's RNG at `seed + i * 7919`, and trimesh's surface sampling is
routed through that seeded RNG (trimesh >= 5 otherwise samples with its own
unseeded generator). Change the seed field to see different draws from
the same edited parameters.

The same boundary applies to the pipeline: tree *structure* is reproducible
across runs, while exact point placements under trimesh 5.x are not.

Requires: numpy, pandas, matplotlib, tkinter (stdlib), trimesh.
