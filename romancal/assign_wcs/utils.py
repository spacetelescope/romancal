import logging

from stcal.alignment.util import compute_s_region_keyword, wcs_bbox_from_shape

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def create_footprint(wcs, shape=None, center=False):
    """Calculate sky footprint

    Parameters
    ----------
    wcs : `gwcs.WCS`
        The WCS information to get the footprint from

    shape : n-tuple or None
       Shape to use if wcs has no defined shape.

    center : bool
        If True use the center of the pixel, otherwise use the corner.

    Returns
    -------
    footprint : `numpy.ndarray`
        The footprint.
    """
    bbox = wcs.bounding_box

    if bbox is None:
        bbox = wcs_bbox_from_shape(shape)

    # footprint is an array of shape (2, 4) - i.e. 4 values for RA and 4 values for
    # Dec - as we are interested only in the footprint on the sky
    footprint = wcs.footprint(bbox, center=center, axis_type="spatial").T
    # take only imaging footprint
    footprint = footprint[:2, :]

    # Make sure RA values are all positive
    negative_ind = footprint[0] < 0
    if negative_ind.any():
        footprint[0][negative_ind] = 360 + footprint[0][negative_ind]

    footprint = footprint.T
    return footprint


def add_s_region(model):
    """
    Calculate the detector's footprint using ``WCS.footprint`` and save it in the
    ``S_REGION`` keyword

    Parameters
    ----------
    model : `~roman_datamodels.datamodels.ImageModel`
        The data model for processing

    Returns
    -------
    A formatted string representing the detector's footprint
    """
    _update_s_region_keyword(
        model, create_footprint(model.meta.wcs, shape=model.shape, center=False)
    )


def _update_s_region_keyword(model, footprint):
    s_region = compute_s_region_keyword(footprint)
    log.info(f"S_REGION VALUES: {s_region}")
    if "nan" in s_region:
        # do not update s_region if there are NaNs.
        log.info("There are NaNs in s_region, S_REGION not updated.")
    else:
        model.meta.wcsinfo.s_region = s_region
        log.info(f"Update S_REGION to {model.meta.wcsinfo.s_region}")
