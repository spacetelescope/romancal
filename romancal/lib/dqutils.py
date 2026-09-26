"""Helpers for propagating data quality arrays.

Frequently steps update model data quality flags from reference files.
There may be multiple dq arrays (dq, dq2), in which case all of these
should be ORed into the models.  This module manages the handling of
these arrays.

dq2 arrays are always created, zero-filled when absent, but the pipeline
never sets any of their bits; those only arrive from reference files.
"""

import numpy as np

__all__ = ["update_dq"]

#: Width of the reference pixel border trimmed from the science array
#: during ramp fitting.
BORDER = 4


def _dq_names(model):
    """Return the data quality attribute names appropriate for ``model``.

    Ramps store the pixel-level arrays as ``pixeldq``/``pixeldq2`` while
    images use ``dq``/``dq2``.
    """
    if "pixeldq" in model:
        return "pixeldq", "pixeldq2"
    return "dq", "dq2"


def _ensure_dq2(model):
    """Return ``model``'s dq2 array, creating an empty one if needed.

    Files written before dq2 existed simply lack the array.  Materialize
    an all-zero one so that "no flags recorded" and "array absent" look
    the same downstream, rather than making every consumer guard.

    Models with no pixel-level data quality array at all, such as the L1
    ``ScienceRawModel``, have nothing for dq2 to accompany and raise
    `TypeError`.
    """
    dq_name, dq2_name = _dq_names(model)
    if dq2_name in model:
        return model[dq2_name]
    if dq_name not in model:
        raise TypeError(
            f"{type(model).__name__} has no {dq_name} array for {dq2_name} to accompany"
        )
    model[dq2_name] = np.zeros(model[dq_name].shape, dtype=model[dq_name].dtype)
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
    """Combine a reference file's data quality arrays into ``model`` with OR.

    Both dq and dq2 are or-ed in when the reference provides them;
    references that carry neither leave ``model`` unchanged apart from
    ensuring that dq2 exists.  ``model`` is updated in place.

    The reference pixel border is trimmed from the reference arrays when
    needed, inferred from the array shapes.

    Parameters
    ----------
    model : `roman_datamodels.datamodels.DataModel`
        Science model to update in place.

    reference : `roman_datamodels.datamodels.DataModel`
        Reference model supplying the flags.
    """
    dq_name, dq2_name = _dq_names(model)
    _ensure_dq2(model)

    for target_name, ref_name in ((dq_name, "dq"), (dq2_name, "dq2")):
        if ref_name not in reference:
            continue
        target = model[target_name]
        ref_array = _align(reference[ref_name], target.shape, ref_name)
        model[target_name] = np.bitwise_or(target, ref_array).astype(target.dtype)
