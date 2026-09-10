"""
Unit validation and in-place conversion for the Roman source catalog.

The functions here scale ``model.data``, ``model.err``, and a
convolved-data array between the various unit conventions used by the
source catalog pipeline (DN/s, MJy/sr, and a configurable flux density
unit).
"""

import astropy.units as u
from roman_datamodels.datamodels import ImageModel


def get_compatible_unit(*arrays):
    """
    Check if multiple arrays have compatible units and return the common
    unit.

    This function verifies that either all arrays are plain NumPy
    ndarrays (no units) or that they are all Astropy Quantity objects
    with the same units.

    Parameters
    ----------
    *arrays : `numpy.ndarray` or `astropy.units.Quantity`
        The data arrays to check. ``None`` entries are ignored.

    Returns
    -------
    result : `astropy.units.Unit` or None
        The common `astropy.units.Unit` object if all inputs are
        Quantity arrays with the same unit, otherwise `None` (when no
        inputs carry units, or when no inputs were supplied).

    Raises
    ------
    ValueError
        If one input is a Quantity and another is not, or if all are
        Quantities but with different units.
    """
    if len(arrays) == 0:
        return None

    # Filter out None values
    arrays = [arr for arr in arrays if arr is not None]
    if len(arrays) == 0:
        return None

    # Check if first array is a Quantity
    is_quantity = [isinstance(arr, u.Quantity) for arr in arrays]

    # All must be quantities or all must not be quantities
    if all(is_quantity):
        # Check that all have the same unit
        first_unit = arrays[0].unit
        for i, arr in enumerate(arrays[1:], start=1):
            if arr.unit != first_unit:
                raise ValueError(
                    f"Incompatible units: array 0 has unit '{first_unit}' "
                    f"but array {i} has unit '{arr.unit}'."
                )
        return first_unit

    if not any(is_quantity):
        return None

    # Mixed types
    raise ValueError("Incompatible types: some arrays have units while others do not.")


def scale_model_arrays(model, convolved_data, factor, *, unit=None):
    """
    Multiply ``model.data``, ``model.err``, and ``convolved_data`` by
    ``factor`` in place.

    Parameters
    ----------
    model : data model
        Mutable input model with ``data`` and ``err`` array fields.

    convolved_data : array-like or None
        Optional convolved-data array. Scaled in place when possible.
        When a ``unit`` is supplied, the ``<<=`` operator creates a new
        Quantity array if ``convolved_data`` is a plain ndarray, so the
        caller must use the return value of this function rather than
        relying on in-place mutation.

    factor : float or `~astropy.units.Quantity`
        Multiplicative factor.

    unit : `~astropy.units.Unit`, optional
        If given, attach this unit to each scaled array using the
        in-place ``<<=`` operator after the multiplication.

    Returns
    -------
    convolved_data : array-like or None
        The convolved-data array.

    Notes
    -----
    Dictionary syntax is used to assign back into ``model`` so that
    on-the-fly schema validation is not triggered for each intermediate
    state of the array.
    """
    for attr in ("data", "err"):
        model[attr] *= factor
        if unit is not None:
            model[attr] <<= unit

    if convolved_data is not None:
        convolved_data *= factor
        if unit is not None:
            # ``<<=`` on a plain ndarray creates a Quantity rather than
            # mutating in place, so we need to return it (if needed)
            convolved_data = convolved_data << unit

    return convolved_data


def validate_and_convert_to_flux_density(
    model, convolved_data, *, flux_unit, l2_to_sb, sb_to_flux
):
    """
    Validate and convert model data arrays to flux-density units.

    Checks that ``model.data`` and ``model.err`` carry compatible units,
    then converts both them and ``convolved_data`` (if not `None`) to
    ``flux_unit``. ``convolved_data`` is allowed to carry different
    units from ``model.data`` because it may originate from a separate
    source (e.g., a detection image from forced photometry).

    Conversion rules
    ----------------
    For arrays without units:

    * Level-2 (`ImageModel`):  ``DN/s -> MJy/sr -> flux_unit``.
    * Level-3 (mosaic):        ``MJy/sr -> flux_unit``.

    For arrays with units:

    * Convert to ``flux_unit`` if the existing unit is equivalent;
      otherwise raise.

    Parameters
    ----------
    model : `ImageModel` or mosaic model
        Data model with mutable ``data`` / ``err`` array fields
        (modified in place).

    convolved_data : array-like or None
        Convolved-data array. May be modified in place; the (possibly
        unit-attached) value is returned.

    flux_unit : `~astropy.units.Unit`
        Target flux-density unit.

    l2_to_sb : float
        Multiplicative factor converting L2 ``DN/s`` to ``MJy/sr``.

    sb_to_flux : `~astropy.units.Quantity`
        Multiplicative factor converting ``MJy/sr`` to ``flux_unit``.

    Returns
    -------
    convolved_data : array-like or None
        The (possibly converted) convolved-data array.

    Raises
    ------
    ValueError
        If ``model.data`` and ``model.err`` have incompatible units, or
        if either ``model`` or ``convolved_data`` carries a unit that is
        not equivalent to ``flux_unit``.
    """
    # Check that model.data and model.err have compatible units
    unit = get_compatible_unit(model.data, model.err)

    if unit is None:
        # No units present - convert to flux density units
        if isinstance(model, ImageModel):
            # Level-2: DN/s -> MJy/sr
            convolved_data = scale_model_arrays(model, convolved_data, l2_to_sb)

        # Level-2 or Level-3: MJy/sr -> flux density
        convolved_data = scale_model_arrays(
            model, convolved_data, sb_to_flux.value, unit=sb_to_flux.unit
        )
    else:
        # Units present - check compatibility and convert
        if not unit.is_equivalent(flux_unit):
            raise ValueError(
                f"Incompatible units: model data has unit '{unit}' "
                f"which is not equivalent to the desired flux unit "
                f"'{flux_unit}'."
            )

        # Convert to desired flux unit
        model["data"] = model["data"].to(flux_unit)
        model["err"] = model["err"].to(flux_unit)

    # Handle convolved_data separately as it may have different units
    if convolved_data is not None:
        conv_unit = get_compatible_unit(convolved_data)
        if conv_unit is None:
            # No units - apply same conversion as model data. When
            # ``unit`` was None above we already scaled convolved_data
            # via scale_model_arrays; otherwise we apply the same
            # L2->SB->flux chain to bring it into flux_unit here.
            if unit is not None:
                if isinstance(model, ImageModel):
                    convolved_data *= l2_to_sb
                convolved_data *= sb_to_flux.value
                convolved_data <<= sb_to_flux.unit
        else:
            if not conv_unit.is_equivalent(flux_unit):
                raise ValueError(
                    f"Incompatible units: convolved_data has unit "
                    f"'{conv_unit}' which is not equivalent to the "
                    f"desired flux unit '{flux_unit}'."
                )

            # Convert to desired flux unit
            convolved_data = convolved_data.to(flux_unit)

    return convolved_data
