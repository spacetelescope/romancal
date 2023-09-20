import logging
import warnings
from typing import Tuple

import gwcs
import numpy as np
from astropy import units as u
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
    """Generate output WCS here based on footprints of all input WCS objects
    Parameters
    ----------
    input_models : list of `~roman_datamodels.datamodels.DataModel`
        Each datamodel must have a ~gwcs.WCS object.

    pscale_ratio : float, optional
        Ratio of input to output pixel scale. Ignored when ``pscale`` is provided.

    pscale : float, None, optional
        Absolute pixel scale in degrees. When provided, overrides
        ``pscale_ratio``.

    rotation : float, None, optional
        Position angle (in degrees) of output image's Y-axis relative to North.
        A value of 0.0 would orient the final output image to be North up.
        The default of `None` specifies that the images will not be rotated,
        but will instead be resampled in the default orientation for the camera
        with the x and y axes of the resampled image corresponding
        approximately to the detector axes.

    shape : tuple of int, None, optional
        Shape of the image (data array) using ``numpy.ndarray`` convention
        (``ny`` first and ``nx`` second). This value will be assigned to
        ``pixel_shape`` and ``array_shape`` properties of the returned
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
    """Create a weight map for use by drizzle"""
    dqmask = build_mask(model.dq, good_bits)

    if weight_type == "ivm":
        if (
            hasattr(model, "var_rnoise")
            and model.var_rnoise is not None
            and model.var_rnoise.shape == model.data.shape
        ):
            with np.errstate(divide="ignore", invalid="ignore"):
                inv_variance = model.var_rnoise**-1
            inv_variance[~np.isfinite(inv_variance)] = 1 * u.s**2 / u.electron**2
        else:
            warnings.warn(
                "var_rnoise array not available. Setting drizzle weight map to 1",
                RuntimeWarning,
            )
            inv_variance = 1.0 * u.s**2 / u.electron**2
        result = inv_variance * dqmask
    elif weight_type == "exptime":
        exptime = model.meta.exposure.exposure_time
        result = exptime * dqmask
    else:
        result = np.ones(model.data.shape, dtype=model.data.dtype) * dqmask

    return result.astype(np.float32)


def build_mask(dqarr, bitvalue):
    """Build a bit mask from an input DQ array and a bitvalue flag

    In the returned bit mask, 1 is good, 0 is bad
    """
    bitvalue = interpret_bit_flags(bitvalue, flag_name_map=pixel)

    if bitvalue is None:
        return np.ones(dqarr.shape, dtype=np.uint8)
    return np.logical_not(np.bitwise_and(dqarr, ~bitvalue)).astype(np.uint8)


def calc_gwcs_pixmap(in_wcs, out_wcs, shape=None):
    """Return a pixel grid map from input frame to output frame."""
    if shape:
        bb = wcs_bbox_from_shape(shape)
        log.debug(f"Bounding box from data shape: {bb}")
    else:
        bb = in_wcs.bounding_box
        log.debug(f"Bounding box from WCS: {in_wcs.bounding_box}")

    grid = gwcs.wcstools.grid_from_bounding_box(bb)
    pixmap = np.dstack(reproject(in_wcs, out_wcs)(grid[0], grid[1]))

    return pixmap


def reproject(wcs1, wcs2):
    """
    Given two WCSs or transforms return a function which takes pixel
    coordinates in the first WCS or transform and computes them in the second
    one. It performs the forward transformation of ``wcs1`` followed by the
    inverse of ``wcs2``.

    Parameters
    ----------
    wcs1, wcs2 : `~astropy.wcs.WCS` or `~gwcs.wcs.WCS` or `~astropy.modeling.Model`
        WCS objects.

    Returns
    -------
    _reproject : func
        Function to compute the transformations.  It takes x, y
        positions in ``wcs1`` and returns x, y positions in ``wcs2``.
    """

    if isinstance(wcs1, fitswcs.WCS):
        forward_transform = wcs1.all_pix2world
    elif isinstance(wcs1, gwcs.WCS):
        forward_transform = wcs1.forward_transform
    elif issubclass(wcs1, Model):
        forward_transform = wcs1
    else:
        raise TypeError(
            "Expected input to be astropy.wcs.WCS or gwcs.WCS "
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
            "Expected input to be astropy.wcs.WCS or gwcs.WCS "
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
