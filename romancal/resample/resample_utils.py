import logging
import warnings
from typing import Tuple

import gwcs
import numpy as np
from astropy import wcs as fitswcs
from astropy.modeling import Model
from astropy.nddata.bitmask import interpret_bit_flags
from stcal.alignment.util import wcs_from_footprints

from romancal.assign_wcs.utils import wcs_bbox_from_shape
from romancal.lib.dqflags import pixel

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def make_output_wcs(
    input_models,
    pscale_ratio=None,
    pscale=None,
    rotation=None,
    shape=None,
    crpix: Tuple[float, float] = None,
    crval: Tuple[float, float] = None,
):
    """
    Generate output WCS here based on footprints of all input WCS objects

    Parameters
    ----------
    input_models : list of `roman_datamodels.datamodels.DataModel`
        Each datamodel must have a `gwcs.wcs.WCS` object.

    pscale_ratio : float, optional
        Ratio of input to output pixel scale. Ignored when ``pscale`` is provided.

    pscale : float, None, optional
        Absolute pixel scale in degrees. When provided, overrides ``pscale_ratio``.

    rotation : float, None, optional
        Position angle (in degrees) of output image's Y-axis relative to North.
        A value of 0.0 would orient the final output image to be North up.
        The default of `None` specifies that the images will not be rotated,
        but will instead be resampled in the default orientation for the camera
        with the x and y axes of the resampled image corresponding
        approximately to the detector axes.

    shape : tuple of int, None, optional
        Shape of the image (data array) using ``numpy.ndarray`` convention
        (``ny`` first and ``nx`` second). This value will be assigned
        to ``pixel_shape`` and ``array_shape`` properties of the returned
        WCS object.

    crpix : tuple of float, None, optional
        Position of the reference pixel in the image array.  If ``ref_pixel`` is not
        specified, it will be set to the center of the bounding box of the
        returned WCS object.

    crval : tuple of float, None, optional
        Right ascension and declination of the reference pixel. Automatically
        computed if not provided.

    Returns
    -------
    output_wcs : object
        WCS object, with defined domain, covering entire set of input frames
    """

    wcslist = [i.meta.wcs for i in input_models]
    for w, i in zip(wcslist, input_models):
        if w.bounding_box is None:
            w.bounding_box = wcs_bbox_from_shape(i.data.shape)
    naxes = wcslist[0].output_frame.naxes

    if naxes != 2:
        raise RuntimeError(f"Output WCS needs 2 axes.{wcslist[0]} has {naxes}.")

    output_wcs = wcs_from_footprints(
        input_models,
        pscale_ratio=pscale_ratio,
        pscale=pscale,
        rotation=rotation,
        shape=shape,
        crpix=crpix,
        crval=crval,
    )

    # Check that the output data shape has no zero length dimensions
    if not np.prod(output_wcs.array_shape):
        raise ValueError(
            f"Invalid output frame shape: {tuple(output_wcs.array_shape)}; dimensions"
            " must have length > 0."
        )

    return output_wcs


def build_driz_weight(model, weight_type=None, good_bits=None):
    """
    Builds the drizzle weight map for resampling.

    Parameters
    ----------
    model : object
        The input model.
    weight_type : str, optional
        The type of weight to use. Allowed values are 'ivm' or 'exptime'.
        Defaults to None.
    good_bits : int or list of int, optional
        The good bits to use for building the mask. Defaults to None.

    Returns
    -------
    numpy.ndarray
        The drizzle weight map.

    Raises
    ------
    ValueError
        If an invalid weight type is provided.

    Examples
    --------
    .. code-block:: none

        model = get_input_model()
        weight_map = build_driz_weight(model, weight_type='ivm', good_bits=[1, 2, 3])
        print(weight_map)
    """

    dqmask = build_mask(model.dq, good_bits)

    if weight_type == "ivm":
        if (
            hasattr(model, "var_rnoise")
            and model.var_rnoise is not None
            and model.var_rnoise.shape == model.data.shape
        ):
            with np.errstate(divide="ignore", invalid="ignore"):
                inv_variance = model.var_rnoise.value**-1
            inv_variance[~np.isfinite(inv_variance)] = 1
        else:
            warnings.warn(
                "var_rnoise array not available. Setting drizzle weight map to 1",
                RuntimeWarning,
            )
            inv_variance = 1.0
        result = inv_variance * dqmask
    elif weight_type == "exptime":
        exptime = model.meta.exposure.exposure_time
        result = exptime * dqmask
    elif weight_type is None:
        result = np.ones(model.data.shape, dtype=model.data.dtype) * dqmask
    else:
        raise ValueError(
            f"Invalid weight type: {weight_type}."
            "Allowed weight types are 'ivm' or 'exptime'."
        )

    return result.astype(np.float32)


def build_mask(dqarr, bitvalue):
    """
    Build a bit mask from an input DQ array and a bitvalue flag.

    Parameters
    ----------
    dqarr : numpy.ndarray
        Input DQ array.
    bitvalue : int
        Bitvalue flag.

    Returns
    -------
    ndarray
        Bit mask where 1 represents good and 0 represents bad.

    Notes
    -----
    - The function interprets the bitvalue flag using the
      `astropy.nddata.bitmask.interpret_bit_flags` function.
    - If the bitvalue is None, the function returns a bit mask with all elements
      set to 1.
    - Otherwise, the function performs a bitwise AND operation between the dqarr and
      the complement of the bitvalue, and then applies a logical NOT operation to
      obtain the bit mask.
    - The resulting bit mask is returned as an ndarray of dtype `numpy.uint8`.
    """
    bitvalue = interpret_bit_flags(bitvalue, flag_name_map=pixel)

    if bitvalue is None:
        return np.ones(dqarr.shape, dtype=np.uint8)
    return np.logical_not(np.bitwise_and(dqarr, ~bitvalue)).astype(np.uint8)


def calc_gwcs_pixmap(in_wcs, out_wcs, shape=None):
    """
    Generate a pixel map grid using the input and output WCS.

    Parameters
    ----------
    in_wcs : `~astropy.wcs.WCS`
        Input WCS.
    out_wcs : `~astropy.wcs.WCS`
        Output WCS.
    shape : tuple, optional
        Shape of the data. If provided, the bounding box will be calculated
        from the shape. If not provided, the bounding box will be calculated
        from the input WCS.

    Returns
    -------
    pixmap : `~numpy.ndarray`
        The calculated pixel map grid.

    """
    if shape:
        bb = wcs_bbox_from_shape(shape)
        log.debug(f"Bounding box from data shape: {bb}")
    else:
        bb = in_wcs.bounding_box
        log.debug(f"Bounding box from WCS: {in_wcs.bounding_box}")

    grid = gwcs.wcstools.grid_from_bounding_box(bb)
    return np.dstack(reproject(in_wcs, out_wcs)(grid[0], grid[1]))


def reproject(wcs1, wcs2):
    """
    Given two WCSs or transforms return a function which takes pixel
    coordinates in the first WCS or transform and computes them in the second
    one. It performs the forward transformation of ``wcs1`` followed by the
    inverse of ``wcs2``.

    Parameters
    ----------
    wcs1, wcs2 : `astropy.wcs.WCS` or `gwcs.wcs.WCS` or `astropy.modeling.Model`
        WCS objects.

    Returns
    -------
    : func
        Function to compute the transformations.  It takes `(x, y)`
        positions in ``wcs1`` and returns `(x, y)` positions in ``wcs2``.
    """

    if isinstance(wcs1, fitswcs.WCS):
        forward_transform = wcs1.all_pix2world
    elif isinstance(wcs1, gwcs.WCS):
        forward_transform = wcs1.forward_transform
    elif issubclass(wcs1, Model):
        forward_transform = wcs1
    else:
        raise TypeError(
            "Expected input to be astropy.wcs.WCS or gwcs.wcs.WCS "
            "object or astropy.modeling.Model subclass"
        )

    if isinstance(wcs2, fitswcs.WCS):
        backward_transform = wcs2.all_world2pix
    elif isinstance(wcs2, gwcs.WCS):
        backward_transform = wcs2.backward_transform
    elif issubclass(wcs2, Model):
        backward_transform = wcs2.inverse
    else:
        raise TypeError(
            "Expected input to be astropy.wcs.WCS or gwcs.wcs.WCS "
            "object or astropy.modeling.Model subclass"
        )

    def _reproject(x, y):
        sky = forward_transform(x, y)
        flat_sky = []
        for axis in sky:
            flat_sky.append(axis.flatten())
        # Filter out RuntimeWarnings due to computed NaNs in the WCS
        warnings.simplefilter("ignore")
        det = backward_transform(*tuple(flat_sky))
        warnings.resetwarnings()
        det_reshaped = []
        for axis in det:
            det_reshaped.append(axis.reshape(x.shape))
        return tuple(det_reshaped)

    return _reproject


def decode_context(context, x, y):
    """
    Get 0-based indices of input images that contributed to (resampled)
    output pixel with coordinates ``x`` and ``y``.

    Parameters
    ----------
    context: numpy.ndarray
        A 3D `~numpy.ndarray` of integral data type.

    x: int, list of integers, numpy.ndarray of integers
        X-coordinate of pixels to decode (3rd index into the ``context`` array)

    y: int, list of integers, numpy.ndarray of integers
        Y-coordinate of pixels to decode (2nd index into the ``context`` array)

    Returns
    -------

    A list of `numpy.ndarray` objects each containing indices of input images
    that have contributed to an output pixel with coordinates ``x`` and ``y``.
    The length of returned list is equal to the number of input coordinate
    arrays ``x`` and ``y``.

    Examples
    --------

    An example context array for an output image of array shape ``(5, 6)``
    obtained by resampling 80 input images.

    .. code-block:: none

        con = np.array(
            [[[0, 0, 0, 0, 0, 0],
              [0, 0, 0, 36196864, 0, 0],
              [0, 0, 0, 0, 0, 0],
              [0, 0, 0, 0, 0, 0],
              [0, 0, 537920000, 0, 0, 0]],
             [[0, 0, 0, 0, 0, 0,],
              [0, 0, 0, 67125536, 0, 0],
              [0, 0, 0, 0, 0, 0],
              [0, 0, 0, 0, 0, 0],
              [0, 0, 163856, 0, 0, 0]],
             [[0, 0, 0, 0, 0, 0],
              [0, 0, 0, 8203, 0, 0],
              [0, 0, 0, 0, 0, 0],
              [0, 0, 0, 0, 0, 0],
              [0, 0, 32865, 0, 0, 0]]],
            dtype=np.int32
        )
        decode_context(con, [3, 2], [1, 4])
        [array([ 9, 12, 14, 19, 21, 25, 37, 40, 46, 58, 64, 65, 67, 77]),
        array([ 9, 20, 29, 36, 47, 49, 64, 69, 70, 79])]

    """
    if context.ndim != 3:
        raise ValueError("'context' must be a 3D array.")

    x = np.atleast_1d(x)
    y = np.atleast_1d(y)

    if x.size != y.size:
        raise ValueError("Coordinate arrays must have equal length.")

    if x.ndim != 1:
        raise ValueError("Coordinates must be scalars or 1D arrays.")

    if not (np.issubdtype(x.dtype, np.integer) and np.issubdtype(y.dtype, np.integer)):
        raise ValueError("Pixel coordinates must be integer values")

    nbits = 8 * context.dtype.itemsize

    return [
        np.flatnonzero([v & (1 << k) for v in context[:, yi, xi] for k in range(nbits)])
        for xi, yi in zip(x, y)
    ]


def resample_range(data_shape, bbox=None):
    # Find range of input pixels to resample:
    if bbox is None:
        xmin = ymin = 0
        xmax = data_shape[1] - 1
        ymax = data_shape[0] - 1
    else:
        ((x1, x2), (y1, y2)) = bbox
        xmin = max(0, int(x1 + 0.5))
        ymin = max(0, int(y1 + 0.5))
        xmax = min(data_shape[1] - 1, int(x2 + 0.5))
        ymax = min(data_shape[0] - 1, int(y2 + 0.5))

    return xmin, xmax, ymin, ymax
