import logging

import numpy as np

from romancal.assign_wcs.utils import wcs_from_footprints, wcs_bbox_from_shape


log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def make_output_wcs(input_models, ref_wcs=None,
                    pscale_ratio=None, pscale=None, rotation=None, shape=None,
                    crpix=None, crval=None):
    """ Generate output WCS here based on footprints of all input WCS objects
    Parameters
    ----------
    input_models : list of `~jwst.datamodel.DataModel`
        Each datamodel must have a ~gwcs.WCS object.

    pscale_ratio : float, optional
        Ratio of input to output pixel scale. Ignored when ``pscale`` is provided.

    pscale : float, None, optional
        Absolute pixel scale in degrees. When provided, overrides
        ``pscale_ratio``.

    rotation : float, None, optional
        Position angle of output image’s Y-axis relative to North.
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
        Position of the reference pixel in the image array.  If ``crpix`` is not
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
        raise RuntimeError("Output WCS needs 2 spatial axes. "
                           f"{wcslist[0]} has {naxes}.")

    output_wcs = wcs_from_footprints(
        input_models,
        pscale_ratio=pscale_ratio,
        pscale=pscale,
        rotation=rotation,
        shape=shape,
        crpix=crpix,
        crval=crval
    )

    # Check that the output data shape has no zero length dimensions
    if not np.product(output_wcs.array_shape):
        raise ValueError(f"Invalid output frame shape: {tuple(output_wcs.array_shape)}")

    return output_wcs
