"""
Image preparation shared by the catalog steps.

These routines mask, background subtract, and convert the units of an
image before sources are measured on it, so that every catalog made
from an image sees the same calibrated data.
"""

import logging
from contextlib import contextmanager

import astropy.units as u
import numpy as np
from roman_datamodels import datamodels
from roman_datamodels.datamodels import ImageModel
from roman_datamodels.dqflags import pixel

from romancal.source_catalog._background import RomanBackground
from romancal.source_catalog._unit_conversion import (
    validate_and_convert_to_flux_density,
)
from romancal.source_catalog._utils import make_model_mask
from romancal.source_catalog._wcs_utils import pixel_area_from_wcs, pixel_area_map
from romancal.source_catalog.psf import add_jitter
from romancal.stpipe.cal_logs import (
    copy_meta_without_cal_logs,
    strip_inherited_cal_logs,
)

log = logging.getLogger(__name__)

# The flux density unit of the catalog and of the saved detection image
DEFAULT_FLUX_UNIT = "nJy"


def make_image_mask(model):
    """
    Create the bad-pixel mask used for cataloging an image.

    Pixels are masked if their data or error is non-finite, if their
    error is not positive, or, for an `ImageModel`, if their DQ has
    DO_NOT_USE set.

    Parameters
    ----------
    model : `ImageModel` or `MosaicModel`
        The input data model.

    Returns
    -------
    mask : 2D `~numpy.ndarray` of bool
        True marks pixels to exclude.
    """
    mask = make_model_mask(model)
    if isinstance(model, ImageModel):
        if model.dq.shape != model.data.shape:
            msg = (
                f"model.dq shape {model.dq.shape} does not match "
                f"model.data shape {model.data.shape}; expected a 2D "
                "DQ array."
            )
            raise ValueError(msg)
        mask |= (model.dq & pixel.DO_NOT_USE) != 0
    return mask


def pixel_scale_arcsec(model):
    """
    Pixel scale in arcsec per pixel, measured at the center of ``model``.

    This is the area-preserving scale, sqrt of the pixel solid angle, to
    match the convention the rest of the catalog uses.
    """
    ysize, xsize = model.data.shape
    area = pixel_area_from_wcs(
        model.meta.wcs, np.array((xsize - 1) / 2.0), np.array((ysize - 1) / 2.0)
    )
    return float(np.sqrt(area.to_value(u.arcsec**2)))


def subtract_background(model, mask, box_size):
    """
    Estimate the 2D background and subtract it from ``model.data``, in place.

    Parameters
    ----------
    model : `ImageModel` or `MosaicModel`
        The data model; ``model.data`` is modified in place.

    mask : 2D `~numpy.ndarray` of bool
        Bad-pixel mask; True pixels are excluded from the estimate.

    box_size : int
        Edge length, in pixels, of the background mesh boxes.
    """
    log.info("Calculating and subtracting background")
    # bad pixel mask rather than coverage_mask to keep finite RMS
    # estimates in holes
    bkg = RomanBackground(model.data, box_size=box_size, mask=mask)
    model.data -= bkg.background


def convert_to_flux_density(
    model, convolved_data=None, *, flux_unit=DEFAULT_FLUX_UNIT, area_map=None
):
    """
    Convert ``model.data`` and ``model.err`` to flux density per pixel.

    Parameters
    ----------
    model : `ImageModel` or `MosaicModel`
        The data model; ``data`` and ``err`` are converted in place.

    convolved_data : array-like or `None`, optional
        A convolved image to convert along with the model arrays.

    flux_unit : str or `~astropy.units.Unit`, optional
        The target flux density unit.

    area_map : `~astropy.units.Quantity`, optional
        The per-pixel solid angle map. Computed from the WCS if not
        given.

    Returns
    -------
    convolved_data : array-like or `None`
        The converted ``convolved_data``.
    """
    flux_unit = u.Unit(flux_unit)
    if area_map is None:
        area_map = pixel_area_map(model.meta.wcs, model.data.shape)
    if "photometry" in model.meta:
        l2_to_sb = model.meta.photometry.conversion_megajanskys
    else:
        l2_to_sb = 1.0
    sb_to_flux = (area_map * (u.MJy / u.sr)).to(flux_unit)
    return validate_and_convert_to_flux_density(
        model,
        convolved_data,
        flux_unit=flux_unit,
        l2_to_sb=l2_to_sb,
        sb_to_flux=sb_to_flux,
    )


@contextmanager
def open_psf_model(ref_file, image_model, match_jitter=True):
    """
    Open an ePSF reference file and add the image's jitter to its PSF.

    Parameters
    ----------
    ref_file : str
        Path to the ePSF reference file.

    image_model : `ImageModel` or `MosaicModel`
        The image whose jitter is applied.

    match_jitter : bool, optional
        If `True`, adjust the PSF for the difference between the image's
        jitter and the reference file's jitter.  If `False`, use the
        reference PSF unchanged.

    Yields
    ------
    psf_model : ePSF reference model
    """
    log.info("Using ePSF reference file: %s", ref_file)
    with datamodels.open(ref_file) as psf_model:
        if match_jitter:
            psf_model.psf = add_jitter(psf_model, image_model)
        else:
            log.info("Skipping jitter matching; using reference PSF as is.")
        yield psf_model


def make_catalog_model(model, cat_model_cls):
    """
    Create a catalog model with metadata copied from an image.

    Upstream ``cal_logs`` are not copied; those remain on the parent
    image product (see issue #2510).

    Parameters
    ----------
    model : `ImageModel` or `MosaicModel`
        The image being cataloged.

    cat_model_cls : type
        The catalog data model class.

    Returns
    -------
    cat_model : data model
        The new catalog model.
    """
    cat_model = cat_model_cls.create_minimal(
        {"meta": copy_meta_without_cal_logs(model.meta)}
    )
    strip_inherited_cal_logs(cat_model)
    cat_model.meta["image"] = {
        "filename": model.meta.filename,
        "file_date": model.meta.file_date,
    }
    # If a data_release_id is present, copy it over.
    # Note that L2 products do not have a data_release_id and this
    # must be populated elsewhere.
    if "data_release_id" in model.meta:
        cat_model.meta.data_release_id = model.meta.data_release_id
    return cat_model
