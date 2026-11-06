import numpy as np
import laspy

# Point clouds are written as LAZ wherever the target filename ends in `.laz`
# (laspy compresses by extension via the lazrs backend and decompresses
# transparently on read, for either extension). The layout names its clouds
# `.laz`; `.las` targets stay plain LAS.


def export_canonical_las(filename, xyz, extra, level="plot", scale=0.001):
    """Write x/y/z plus the canonical extra-dim schema for `level`.

    xyz   : (N, 3) coordinates.
    extra : {name: array} covering the canonical fields for the level. Missing
            fields are filled with their schema default; unknown names error.
    level : 'plot' or 'tree' -- selects schema.PLOT_FIELDS / TREE_FIELDS.

    The LAS `classification` field is set from the canonical 'classification'
    dim. Every extra dim is written at its schema dtype (int16/int8/uint8/
    float32), so plot and tree clouds share one on-disk schema across the chain.
    """
    import schema
    fields = schema.PLOT_FIELDS if level == "plot" else schema.TREE_FIELDS
    xyz = np.asarray(xyz, dtype=np.float64)
    n = len(xyz)
    cols = schema.build_extra(fields, n, **{k: v for k, v in (extra or {}).items()})

    # point format 0: no RGB (nothing in the pipeline writes colour; the
    # unused 6 RGB bytes of format 2 were pure overhead). LAS 1.4 is the
    # target version per the design spec (Section 7.5): point formats 0-5
    # remain valid in 1.4, so only the version tag changes.
    header = laspy.LasHeader(point_format=0, version="1.4")
    header.scales = [scale, scale, scale]
    header.offsets = xyz[:, :3].min(axis=0) if n else [0.0, 0.0, 0.0]
    dims = [laspy.ExtraBytesParams(name=name, type=schema.dtype_of(name))
            for name in fields]
    header.add_extra_dims(dims)

    las = laspy.LasData(header)
    las.x = xyz[:, 0]; las.y = xyz[:, 1]; las.z = xyz[:, 2]
    las.classification = cols["classification"].astype(np.uint8)
    for name in fields:
        setattr(las, name, cols[name])
    las.write(filename)


def read_canonical_las(filename, level="plot"):
    """Read a canonical cloud into (xyz, {name: array}) for `level`'s fields."""
    import schema
    fields = schema.PLOT_FIELDS if level == "plot" else schema.TREE_FIELDS
    las = laspy.read(filename)
    xyz = np.column_stack([np.asarray(las.x), np.asarray(las.y), np.asarray(las.z)])
    out = {}
    for name in fields:
        if hasattr(las, name):
            out[name] = np.asarray(getattr(las, name))
        else:
            out[name] = np.full(len(xyz), schema.default_of(name),
                                dtype=schema.dtype_of(name))
    return xyz, out


def export_points_las(filename, points, scale=0.001, extra_fields=None):
    """
    Export points to a LAS file.

    points: (N, 3+) array, columns:
        0, 1, 2 = x, y, z  (real-world coordinates)
        3       = semantic class: 1 = trunk/stem, 0 = branch/foliage  [optional]
        4       = tree id: integer per placed tree                    [optional]

    extra_fields: optional {name: array} of further per-point columns, each
    registered as a LAS extra dimension. dtype is inferred: bool/uint8 -> uint8,
    floating -> float32, otherwise int32. (float32 is needed for the density /
    branch-order ratio columns written by the normalization step.)

    Coordinates go through laspy's scaled fields (lowercase las.x/y/z) so the
    real-world size is preserved. Both extra dimensions are registered in a
    single add_extra_dims() call: registering them with two separate
    add_extra_dim() calls left the second one (tree_id) unbound, so it was
    written as all zeros while the first (semantic) came through fine.
    """
    points = np.asarray(points, dtype=np.float64)
    n = points.shape[0]
    xyz = points[:, :3]

    semantic = (points[:, 3].astype(np.uint8) if points.shape[1] > 3
                else np.zeros(n, dtype=np.uint8))
    tree_id = (points[:, 4].astype(np.int32) if points.shape[1] > 4
               else np.zeros(n, dtype=np.int32))

    # point format 0: no RGB (see export_canonical_las)
    header = laspy.LasHeader(point_format=0, version="1.2")
    header.scales = [scale, scale, scale]
    header.offsets = xyz.min(axis=0) if n else [0.0, 0.0, 0.0]
    extra_fields = dict(extra_fields or {})
    extra_arrays = {}
    dims = [
        laspy.ExtraBytesParams(name="semantic", type=np.uint8),
        laspy.ExtraBytesParams(name="tree_id", type=np.int32),
    ]
    for name, values in extra_fields.items():
        arr = np.asarray(values)
        if arr.dtype == bool:
            dtype = np.uint8
        elif np.issubdtype(arr.dtype, np.floating):
            dtype = np.float32
        elif np.issubdtype(arr.dtype, np.integer):
            # honour the caller's integer width (int8/int16/int32, unsigned too)
            # so the canonical schema's int8/int16 dims land at their real size
            # instead of being widened to int32.
            dtype = arr.dtype.type
        else:
            dtype = np.int32
        extra_arrays[name] = arr.astype(dtype)
        dims.append(laspy.ExtraBytesParams(name=name, type=dtype))
    header.add_extra_dims(dims)

    las = laspy.LasData(header)          # schema (incl. both extra dims) fixed here
    las.x = xyz[:, 0]                    # then assign every field once
    las.y = xyz[:, 1]
    las.z = xyz[:, 2]
    las.classification = semantic
    las.semantic = semantic
    las.tree_id = tree_id
    for name, arr in extra_arrays.items():
        setattr(las, name, arr)

    las.write(filename)


def read_points_las(filename, extra_fields=None):
    """Read a LAS written by export_points_las back into an (N, 5[+]) array.

    Columns are x, y, z, semantic, tree_id, then any names given in
    extra_fields, in order. Missing extra dimensions come back as zeros so a
    raw cloud can be read with the same call as a trimmed one.
    """
    las = laspy.read(filename)
    cols = [np.asarray(las.x, dtype=np.float64),
            np.asarray(las.y, dtype=np.float64),
            np.asarray(las.z, dtype=np.float64)]
    n = len(cols[0])

    def _get(name):
        try:
            return np.asarray(getattr(las, name), dtype=np.float64)
        except (AttributeError, ValueError):
            return np.zeros(n, dtype=np.float64)

    cols.append(_get("semantic"))
    cols.append(_get("tree_id"))
    for name in (extra_fields or []):
        cols.append(_get(name))
    return np.column_stack(cols)