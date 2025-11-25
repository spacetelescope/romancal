"""
Helper function for generating per-filter catalogs in multiband catalog
creation.
"""

import logging

import numpy as np
from astropy.table import join
from roman_datamodels import datamodels

from romancal.multiband_catalog.utils import add_filter_to_colnames
from romancal.source_catalog.psf_matching import (
    compute_psf_correction_factors,
    create_psf_matched_image,
    get_filter_wavelength,
)
from romancal.source_catalog.source_catalog import RomanSourceCatalog
from romancal.source_catalog.utils import get_ee_spline

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def create_filter_catalog(
    model,
    filter_name,
    ref_filter,
    ref_wavelength,
    segment_img,
    star_kernel_fwhm,
    detection_catobj,
    ref_model,
    ref_filter_catalog,
    ref_psf_model,
    fit_psf,
    get_reference_file_func,
):
    """
    Create catalog(s) for a single filter, including PSF matching if needed.

    This function handles three cases:
    1. Reference filter: Only original measurements (no PSF-matched catalog)
    2. Bluer filter: Normal PSF matching (convolve to reference filter)
    3. Redder filter: Synthetic PSF matching via correction factors

    Parameters
    ----------
    model : ImageModel or MosaicModel
        The input image model for this filter.

    filter_name : str
        Name of the filter (e.g., 'F158').

    ref_filter : str
        Name of the reference filter for PSF matching.

    ref_wavelength : int
        Wavelength of the reference filter in microns.

    segment_img : ndarray
        Segmentation image from detection.

    star_kernel_fwhm : float
        FWHM of star kernel for source catalog.

    detection_catobj : RomanSourceCatalog
        Detection catalog object.

    ref_model : ImageModel or MosaicModel
        The reference filter image model (for computing correction factors).

    ref_filter_catalog : Table or None
        The reference filter's catalog (needed for redder filters).

    ref_psf_model : EpsfRefModel
        PSF reference model for the reference filter.

    fit_psf : bool
        Whether to fit PSFs in the catalog.

    get_reference_file_func : callable
        Function to get reference files (self.get_reference_file from step).

    Returns
    -------
    result : dict
        Dictionary containing:
        - 'catalog': Combined catalog table for this filter
        - 'ee_fractions': Dictionary of ee_fractions for this filter

    Raises
    ------
    ValueError
        If trying to process a redder filter before the reference filter.
    """
    # Create mask
    mask = ~np.isfinite(model.data) | ~np.isfinite(model.err) | (model.err <= 0)

    # Load PSF reference model (needed for PSF matching and PSF fitting)
    log.info(f"Creating catalog for {filter_name} image")
    ref_file = get_reference_file_func(model, "epsf")
    log.info("Using ePSF reference file: %s", ref_file)
    psf_model = datamodels.open(ref_file)

    apcorr_ref = get_reference_file_func(model, "apcorr")
    ee_spline = get_ee_spline(model, apcorr_ref)

    # Create catalog with original (non-PSF-matched) data
    log.info(f"Creating catalog for original {filter_name} image")
    catobj_original = RomanSourceCatalog(
        model,
        segment_img,
        None,
        star_kernel_fwhm,
        fit_psf=fit_psf,
        psf_model=psf_model if fit_psf else None,
        mask=mask,
        detection_cat=detection_catobj,
        cat_type="dr_band",
        ee_spline=ee_spline,
    )

    # Store reference filter catalog for later use with redder filters
    updated_ref_filter_catalog = ref_filter_catalog
    if filter_name == ref_filter:
        updated_ref_filter_catalog = catobj_original.catalog

    # Add the filter name to the column names
    cat_original = add_filter_to_colnames(catobj_original.catalog, filter_name)

    # Store ee_fractions for this filter
    ee_fractions = {}
    ee_fractions[filter_name.lower()] = cat_original.meta["ee_fractions"]

    # Clear filter catalog metadata
    cat_original.meta = None

    # Determine if this filter is bluer or redder than reference
    filter_wavelength = get_filter_wavelength(filter_name)

    # Create PSF-matched catalog based on filter position
    if filter_name == ref_filter:
        # Reference filter case - only include original measurements
        log.info(
            f"Reference filter {filter_name}: including only original measurements"
        )
        cat = cat_original

    elif filter_wavelength < ref_wavelength:
        # Bluer filter - normal PSF matching (convolve to the reference
        # image's PSF)
        log.info(f"Creating PSF-matched image for {filter_name}")
        psf_matched_model = create_psf_matched_image(
            model,
            psf_model,
            ref_psf_model,
        )

        log.info(f"Creating catalog for PSF-matched {filter_name} image")
        catobj_matched = RomanSourceCatalog(
            psf_matched_model,
            segment_img,
            None,
            star_kernel_fwhm,
            fit_psf=False,  # No PSF fitting on matched images
            psf_model=None,
            mask=mask,
            detection_cat=detection_catobj,
            cat_type="psf_matched",
            ee_spline=ee_spline,
        )

        # Add filter name with "m" suffix for PSF-matched columns
        filter_name_matched = f"{filter_name}m"
        cat_matched = add_filter_to_colnames(
            catobj_matched.catalog, filter_name_matched
        )
        cat_matched.meta = None

        # Merge the original and PSF-matched catalogs
        cat = join(cat_original, cat_matched, keys="label", join_type="outer")

    else:
        # Redder filter - synthetic PSF matching using correction factors
        log.info(
            f"Creating synthetic PSF-matched catalog for "
            f"{filter_name} (redder than reference {ref_filter})"
        )

        if ref_filter_catalog is None:
            raise ValueError(
                f"Reference filter {ref_filter} catalog not yet created. "
                "Cannot compute correction factors for redder filter "
                f"{filter_name}."
            )

        # Compute correction factors by matching reference image to
        # the redder filter
        correction_factors = compute_psf_correction_factors(
            ref_model=ref_model,
            ref_psf_model=ref_psf_model,
            ref_catalog=ref_filter_catalog,
            target_model=model,
            target_psf_model=psf_model,
            segment_img=segment_img,
            star_kernel_fwhm=star_kernel_fwhm,
            detection_cat=detection_catobj,
            mask=mask,
            ee_spline=ee_spline,
        )

        # Create synthetic PSF-matched catalog by applying
        # correction factors to the original catalog
        cat_synthetic = cat_original.copy()

        # Apply correction factors to flux columns
        for flux_col, correction in correction_factors.items():
            # flux_col names contain the reference filter name, so
            # we need to replace it with the current filter name
            flux_col_syn = flux_col.replace(
                f"_{ref_filter.lower()}_", f"_{filter_name.lower()}_"
            )

            # Apply correction to flux column
            if flux_col_syn in cat_synthetic.colnames:
                cat_synthetic[flux_col_syn] *= correction

            # Also apply to error columns
            err_col_syn = flux_col_syn.replace("_flux", "_flux_err")
            if err_col_syn in cat_synthetic.colnames:
                cat_synthetic[err_col_syn] *= correction

            # Handle magnitude columns (subtract 2.5 * log10(C))
            mag_col_syn = flux_col_syn.replace("_flux", "_abmag")
            if mag_col_syn in cat_synthetic.colnames:
                with np.errstate(divide="ignore", invalid="ignore"):
                    mag_correction = -2.5 * np.log10(correction)
                    mag_correction = np.where(
                        ~np.isfinite(mag_correction),
                        0.0,
                        mag_correction,
                    )
                cat_synthetic[mag_col_syn] += mag_correction

        # Rename columns to add 'm' suffix and keep only flux columns
        filter_name_matched = f"{filter_name}m"
        cols_to_keep = ["label"]  # Always keep label for joining
        for colname in list(cat_synthetic.colnames):
            # Replace filter name with filter name + 'm'
            # e.g., kron_f158_flux -> kron_f158m_flux
            # We not keep columns that *end* with the filter name, e.g.,
            # sharpness_f158.
            if f"_{filter_name.lower()}_" in colname:
                new_colname = colname.replace(
                    f"_{filter_name.lower()}_",
                    f"_{filter_name_matched.lower()}_",
                )
                cat_synthetic.rename_column(colname, new_colname)
                cols_to_keep.append(new_colname)

        # Keep only the PSF-matched columns (and label)
        cat_synthetic = cat_synthetic[cols_to_keep]

        # Merge the two catalogs
        cat = join(
            cat_original,
            cat_synthetic,
            keys="label",
            join_type="outer",
        )
        log.info(
            f"Synthetic PSF-matched catalog created for {filter_name} "
            f"using correction factors"
        )

    log.info(f"Completed catalog for {filter_name} filter")

    return {
        "catalog": cat,
        "ee_fractions": ee_fractions,
        "ref_filter_catalog": updated_ref_filter_catalog,
    }
