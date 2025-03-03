import logging
import math

import numpy as np
from stcal.alignment.util import compute_scale, wcs_from_sregions
from stcal.alignment.util import wcs_from_footprints as stcal_wcs_from_footprints
from stcal.resample.utils import compute_mean_pixel_area

from romancal.assign_wcs.utils import wcs_bbox_from_shape

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def wcs_from_footprints(
    input_models,
    pscale=None,
    pscale_ratio=1.0,
    shape=None,
    rotation=None,
    crpix=None,
    crval=None,
):
    sregions = []
    ref_wcs = None
    ref_wcsinfo = None
    ref_shape = None
    with input_models:
        for model in input_models:
            if ref_wcs is None:
                ref_wcs = model.meta.wcs
                ref_wcsinfo = model.meta.wcsinfo
                ref_shape = model.data.shape
            sregions.append(model.meta.wcs.footprint())
            input_models.shelve(model, modify=False)

    if pscale is None:
        pscale = (
            compute_scale(
                ref_wcs,
                fiducial=np.array([ref_wcsinfo["ra_ref"], ref_wcsinfo["dec_ref"]]),
            )
            * pscale_ratio
        )
    else:
        pscale_ratio = pscale / np.rad2deg(
            math.sqrt(compute_mean_pixel_area(ref_wcs, shape=ref_shape))
        )

    wcs = wcs_from_sregions(
        sregions,
        ref_wcs=ref_wcs,
        ref_wcsinfo=ref_wcsinfo,
        pscale_ratio=pscale_ratio,
        pscale=pscale,
        shape=shape,
        rotation=rotation,
        crpix=crpix,
        crval=crval,
    )

    return wcs, pscale * 3600.0, pscale_ratio


def make_output_wcs(
    input_models,
    pscale_ratio=None,
    pscale=None,
    rotation=None,
    shape=None,
    crpix: tuple[float, float] | None = None,
    crval: tuple[float, float] | None = None,
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
    for w, i in zip(wcslist, input_models, strict=False):
        if w.bounding_box is None:
            w.bounding_box = wcs_bbox_from_shape(i.data.shape)
    naxes = wcslist[0].output_frame.naxes

    if naxes != 2:
        raise RuntimeError(f"Output WCS needs 2 axes.{wcslist[0]} has {naxes}.")

    output_wcs = stcal_wcs_from_footprints(
        wcslist,
        None,
        dict(input_models[0].meta.wcsinfo),
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
