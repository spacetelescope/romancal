"""Helpers for propagating data quality arrays.

Roman carries two pixel-level data quality arrays:

``dq`` (``pixeldq`` on ramps)
    The actionable data quality array.  Pipeline behavior is driven from
    these bits, chiefly ``DO_NOT_USE``.

``dq2`` (``pixeldq2`` on ramps)
    Static, informational flags.  These are propagated but never consulted
    by the pipeline, by ``good_bits`` handling, or by ``stcal``/``drizzle``.
    Anything that must change pipeline behavior has to also set
    ``DO_NOT_USE`` in ``dq``.
"""

import numpy as np

__all__ = ["DQ2_DTYPE", "dq_names", "ensure_dq2", "update_dq"]

#: dtype used for newly created ``dq2`` arrays.  Matches ``dq``.
DQ2_DTYPE = np.uint32

#: Width of the reference pixel border trimmed from the science array
#: during ramp fitting.
BORDER = 4


def dq_names(model):
    """Return the data quality attribute names appropriate for ``model``.

    Ramps store the pixel-level arrays as ``pixeldq``/``pixeldq2`` while
    images use ``dq``/``dq2``.

    Parameters
    ----------
    model : `roman_datamodels.datamodels.DataModel`
        Model to inspect.

    Returns
    -------
    tuple of str
        The ``(dq, dq2)`` attribute names.
    """
    if "pixeldq" in model:
        return "pixeldq", "pixeldq2"
    return "dq", "dq2"


def ensure_dq2(model):
    """Return ``model``'s ``dq2`` array, creating an empty one if needed.

    Files written before ``dq2`` existed, and reference files that do not
    carry one, simply lack the array.  Rather than requiring every consumer
    to guard against that, materialize an all-zero array so that "no flags
    recorded" and "array absent" look the same downstream.

    Parameters
    ----------
    model : `roman_datamodels.datamodels.DataModel`
        Model to update in place.

    Returns
    -------
    `numpy.ndarray` or None
        The existing or newly created ``dq2`` array, or None if ``model``
        has no data quality array to match.
    """
    dq_name, dq2_name = dq_names(model)
    if dq2_name in model:
        return model[dq2_name]
    if dq_name not in model:
        return None
    model[dq2_name] = np.zeros(model[dq_name].shape, dtype=DQ2_DTYPE)
    return model[dq2_name]


def _align(ref_array, shape, name):
    """Match a reference array to the science array shape.

    Reference files are full frame while post-ramp-fitting science arrays
    have the reference pixel border removed.  Whether a given reference
    needs trimming depends on the reference type rather than on the step,
    so infer it from the shapes instead of requiring callers to know.
    """
    if ref_array.shape == shape:
        return ref_array

    untrimmed = tuple(size + 2 * BORDER for size in shape)
    if ref_array.shape == untrimmed:
        return ref_array[BORDER:-BORDER, BORDER:-BORDER]

    raise ValueError(
        f"reference {name} has shape {ref_array.shape}, which matches "
        f"neither the science shape {shape} nor the untrimmed shape "
        f"{untrimmed}"
    )


def update_dq(model, reference):
    """Combine a reference file's data quality arrays into ``model``.

    Both ``dq`` and ``dq2`` are or-ed in when the reference provides them;
    references that carry neither leave ``model`` unchanged apart from
    ensuring that ``dq2`` exists.  ``model`` is updated in place.

    The reference pixel border is trimmed from the reference arrays when
    needed, inferred from the array shapes.  This is inferred rather than
    passed in because the convention differs between reference types:
    dark references are full frame while flat references are already
    trimmed.  A shape matching neither convention raises `ValueError`.

    The operation is a bitwise or and so is idempotent, which makes it safe
    to call for steps where ``stcal`` has already folded the same reference
    ``dq`` into the science array.

    Parameters
    ----------
    model : `roman_datamodels.datamodels.DataModel`
        Science model to update in place.

    reference : `roman_datamodels.datamodels.DataModel`
        Reference model supplying the flags.

    Returns
    -------
    `roman_datamodels.datamodels.DataModel`
        The updated ``model``, for convenience.
    """
    dq_name, dq2_name = dq_names(model)
    ensure_dq2(model)

    for target_name, ref_name in ((dq_name, "dq"), (dq2_name, "dq2")):
        if ref_name not in reference:
            continue
        target = model[target_name]
        ref_array = _align(reference[ref_name], target.shape, ref_name)
        model[target_name] = np.bitwise_or(target, ref_array).astype(target.dtype)

    return model
