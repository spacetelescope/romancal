"""
Module for PSF matching of images.
"""

import logging

import numpy as np
from astropy.convolution import convolve_fft
from roman_datamodels.datamodels import ImageModel, MosaicModel

from romancal.multiband_catalog.utils import add_filter_to_colnames
from romancal.source_catalog.psf import (
    create_convolution_kernel,
    create_l3_psf_model,
)
from romancal.source_catalog.utils import copy_model_arrays

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def create_psf_matched_image(
    model,
    psf_model,
    target_psf_model,
    min_fft_power_ratio=1e-5,
):
    """
    Create a PSF-matched version of the input model.

    This function convolves the input image with a matching kernel to
    transform its PSF to match a target PSF.

    Parameters
    ----------
    model : `ImageModel` or `MosaicModel`
        The input image to PSF-match. The image is assumed to be
        background subtracted. If the input model data and error arrays
        have units, they will be preserved in the output.

    psf_model : EpsfRefModel
        The PSF model for the input ``model`` image.

    target_psf_model : EpsfRefModel
        The PSF model for the reference/target PSF (typically the
        broadest PSF in a multiband set).

    min_fft_power_ratio : float, optional
        Regularization parameter for the matching kernel calculation.
        Controls the scale of regularization in terms of the peak power
        of the input PSF's FFT. Larger values correspond to stronger
        regularization. Default is 1e-5.

    Returns
    -------
    matched_model : `ImageModel` or `MosaicModel`
        A copy of the input model with PSF-matched data and error
        arrays.
    """
    if not isinstance(model, (ImageModel, MosaicModel)):
        raise ValueError("model must be an ImageModel or MosaicModel")

    # Get filter information
    input_filter = model.meta.instrument.optical_element
    target_filter = target_psf_model.meta.instrument.optical_element

    log.info(f"Creating PSF-matched image: {input_filter} -> {target_filter}")

    # Check if input PSF is the same or broader than target PSF by
    # comparing filter wavelengths (longer wavelength -> broader PSF)
    input_wavelength = get_filter_wavelength(input_filter)
    target_wavelength = get_filter_wavelength(target_filter)
    if input_wavelength >= target_wavelength:
        log.warning(
            f"Input filter ({input_filter}, {input_wavelength} μm) has the "
            f"same or longer wavelength than target filter ({target_filter}, "
            f"{target_wavelength} μm). No PSF matching performed."
        )

        # Return the input model
        return model

    # Create the matching kernel
    log.info("Computing PSF matching kernel")

    # Create L3 PSF models to access PSF data arrays
    input_l3_psf_model = create_l3_psf_model(psf_model)
    target_l3_psf_model = create_l3_psf_model(target_psf_model)

    # Check oversampling factors -- create_convolution_kernel requires
    # a single oversampling factor.
    input_oversampling = input_l3_psf_model.oversampling
    if input_oversampling[0] != input_oversampling[1]:
        msg = (
            "Input PSF model has different oversampling factors in x and y directions."
        )
        raise ValueError(msg)
    target_oversampling = target_l3_psf_model.oversampling
    if np.all(input_oversampling != target_oversampling):
        msg = "Input and target PSF models have different oversampling factors."
        raise ValueError(msg)
    oversampling = input_oversampling[0]

    matching_kernel = create_convolution_kernel(
        input_l3_psf_model.data,
        target_l3_psf_model.data,
        downsample=oversampling,
    )

    # Convolve the data
    # np.asarray() is needed to convert to an ndarray from
    # asdf.tags.core.ndarray.NDArrayType, which can cause issues with
    # convolve_fft
    log.info("Convolving image data with matching kernel")
    matched_data = convolve_fft(
        np.asanyarray(model.data),  # use np.asanyarray to preserve units
        matching_kernel,
        preserve_nan=True,
        normalize_kernel=True,
    )

    # Propagate errors - convolve variance with kernel^2
    log.info("Propagating errors")
    if model.err is not None:
        matched_variance = convolve_fft(
            model.err**2,
            matching_kernel**2,
            preserve_nan=True,
            normalize_kernel=False,
        )
        matched_err = np.sqrt(np.abs(matched_variance))
    else:
        matched_err = None

    # Create output model with copied data and err arrays
    matched_model = copy_model_arrays(model)
    matched_model.data = matched_data
    if matched_err is not None:
        matched_model.err = matched_err

    log.info(f"Created PSF-matched image: {input_filter} -> {target_filter}")

    return matched_model


def get_filter_wavelength(filter_name):
    """
    Extract approximate wavelength from filter name.

    Parameters
    ----------
    filter_name : str
        Filter name (e.g., 'F158', 'F184', 'f158', 'f158m')

    Returns
    -------
    wavelength : int
        Approximate wavelength in microns (e.g., 1.58, 1.84), or 0 if
        cannot parse.
    """
    try:
        # Remove 'm' suffix if present (for PSF-matched column names)
        name = filter_name.rstrip("m")
        return int(name[1:]) / 100.0  # Convert to microns
    except (ValueError, IndexError):
        return 0


def get_reddest_filter(library):
    """
    Get the reddest filter in the model library.

    The reddest filter is typically the one with the broadest PSF.

    Parameters
    ----------
    library : `~romancal.datamodels.ModelLibrary`
        The model library containing images from different filters.

    Returns
    -------
    reddest_filter : str
        The name of the reddest filter.

    Notes
    -----
    This function uses a simple heuristic that extracts filter
    wavelengths from the filter names.
    """
    # Get list of filters in the library
    filters = []
    with library:
        for model in library:
            filter_name = model.meta.instrument.optical_element
            if filter_name not in filters:
                filters.append(filter_name)
            library.shelve(model, modify=False)

    # Simple heuristic to extract wavelength from filter name
    # Roman WFI filters: F062, F087, F106, F129, F146, F158, F184, F213
    filter_wavelengths = {filt: get_filter_wavelength(filt) for filt in filters}

    # Return the longest wavelength filter
    reference_filter = max(filter_wavelengths, key=filter_wavelengths.get)

    log.info(
        f"Selected {reference_filter} as reference filter for PSF matching "
        f"(longest wavelength among {filters})"
    )

    return reference_filter


def compute_psf_correction_factors(
    ref_model,
    ref_psf_model,
    ref_catalog,
    target_model,
    target_psf_model,
    segment_img,
    star_kernel_fwhm,
    detection_cat,
    mask,
    ee_spline,
):
    """
    Compute correction factors for PSF-matched photometry on the
    reference filter.

    When the PSF match reference filter is not the reddest filter, we
    need to create synthetic PSF-matched photometry by:

    1. PSF-matching the reference image to the target (redder) filter's PSF.
    2. Measuring photometry on the PSF-matched reference image.
    3. Computing correction factors C = flux_original / flux_matched.
       These factors are typically larger than 1.0 because the
       PSF-matched image has broader PSF and thus lower flux in a fixed
       aperture. These factors are multiplied to the original measured
       fluxes in the target (redder) filter.

    Parameters
    ----------
    ref_model : `ImageModel` or `MosaicModel`
        The reference filter image (original, not PSF-matched).

    ref_psf_model : EpsfRefModel
        PSF model for the reference filter.

    ref_catalog : `~astropy.table.Table`
        Catalog from the original reference image.

    target_model : `ImageModel` or `MosaicModel`
        The target (redder) filter image that will be matched.

    target_psf_model : EpsfRefModel
        PSF model for the target filter.

    segment_img : `~photutils.segmentation.SegmentationImage`
        Segmentation map for source extraction.

    star_kernel_fwhm : float
        FWHM for star detection kernel.

    detection_cat : RomanSourceCatalog
        Detection catalog.

    mask : array-like
        Boolean mask array where True indicates valid pixels.

    ee_spline : callable
        The encircled energy spline function.

    Returns
    -------
    correction_factors : dict
        Dictionary mapping flux column names to arrays of correction
        factors. Values are arrays with one correction factor per
        source.
    """
    from romancal.source_catalog.source_catalog import RomanSourceCatalog

    ref_filter = ref_model.meta.instrument.optical_element
    target_filter = target_model.meta.instrument.optical_element
    log.info(
        f"Computing PSF correction factors for {ref_filter} matched to {target_filter}"
    )

    # PSF-match the reference image to the target filter
    psf_matched_ref = create_psf_matched_image(
        ref_model,
        ref_psf_model,
        target_psf_model,
    )

    # Measure photometry on the PSF-matched reference image
    catobj_matched = RomanSourceCatalog(
        psf_matched_ref,
        segment_img,
        None,
        star_kernel_fwhm,
        fit_psf=False,
        psf_model=None,
        mask=mask,
        detection_cat=detection_cat,
        cat_type="psf_matched",
        ee_spline=ee_spline,
    )

    # add filter name to catalog column names
    cat_matched = add_filter_to_colnames(catobj_matched.catalog, ref_filter)

    # Compute correction factors C = flux_original / flux_matched
    # for each flux type (segment, kron, aper*)
    correction_factors = {}

    # Flux column patterns to compute corrections for
    flux_columns = [col for col in ref_catalog.colnames if col.endswith("_flux")]

    for flux_col in flux_columns:
        if flux_col not in ref_catalog.colnames or flux_col not in cat_matched.colnames:
            continue

        flux_original = ref_catalog[flux_col]
        flux_matched = cat_matched[flux_col]

        # Compute correction factor, avoiding division by zero
        with np.errstate(divide="ignore", invalid="ignore"):
            correction = flux_original / flux_matched
            # Set correction to 1.0 where original flux is zero or
            # invalid
            correction = np.where(
                (flux_matched == 0) | ~np.isfinite(correction), 1.0, correction
            )

        correction_factors[flux_col] = correction

    return correction_factors
