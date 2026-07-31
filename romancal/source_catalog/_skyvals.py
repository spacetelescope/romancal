"""
Helpers for computing HEALPix-based sky summary arrays for segmentation products.
"""

import logging

import numpy as np
from astropy import units as u
from astropy_healpix import HEALPix

from romancal.source_catalog._utils import get_pixel_area_sr

log = logging.getLogger(__name__)

HEALPIX17_NSIDE = 2**17
HEALPIX11_DIVISOR = 4**6
SKYVALS_DTYPE = np.dtype(
    [
        ("healpix17", np.int64),
        ("data", np.float32),
        ("err", np.float32),
        ("covfrac", np.float32),
    ]
)


def compute_skyvals(input_model, segmentation, bad_pixel_mask):
    """
    Compute HEALPix-aggregated sky statistics from masked image pixels.

    Parameters
    ----------
    input_model
        Roman image model containing data/err arrays and metadata.
    segmentation : np.ndarray
        2D segmentation labels; source pixels are non-zero.
    bad_pixel_mask : np.ndarray
        Boolean 2D mask where True means "exclude pixel from sky calculation".
    """
    data = input_model.data
    err = input_model.err
    wcs = input_model.meta.wcs
    pixel_area_sr = get_pixel_area_sr(input_model)

    # Valid sky pixels are those not rejected by upstream masks and not part
    # of detected sources.
    valid = (~bad_pixel_mask) & (segmentation == 0)
    if not np.any(valid):
        # Fallback: if every source-masked pixel is rejected, use all unmasked
        # pixels so we still produce per-exposure sky estimates.
        valid = ~bad_pixel_mask
        if not np.any(valid):
            log.warning(
                "No usable pixels for skyvals after bad-pixel masking; "
                "returning empty skyvals/healpix11_cov arrays."
            )
            return np.empty((0,), dtype=SKYVALS_DTYPE), np.empty((0,), dtype=np.int64)

    yy, xx = np.nonzero(valid)
    ra, dec = wcs.pixel_to_world_values(xx, yy)

    # Convert sky coordinates to NEST-ordered HEALPix indices at nside=2**17.
    healpix = HEALPix(nside=HEALPIX17_NSIDE, order="nested", frame="icrs")
    healpix17 = healpix.lonlat_to_healpix(ra * u.deg, dec * u.deg).astype(np.int64)

    vals_data = data[yy, xx]
    vals_err = err[yy, xx]

    # Sort/group by healpixel so per-healpixel reductions are contiguous.
    order = np.argsort(healpix17)
    hp_sorted = healpix17[order]
    data_sorted = vals_data[order]
    err_sorted = vals_err[order]

    unique_hp, starts, counts = np.unique(
        hp_sorted, return_index=True, return_counts=True
    )

    med_data = np.empty(unique_hp.size, dtype=np.float32)
    med_err = np.empty(unique_hp.size, dtype=np.float32)
    # Use robust medians for sky and error estimates within each healpixel.
    for i, (start, count) in enumerate(zip(starts, counts, strict=True)):
        section = slice(start, start + count)
        med_data[i] = np.nanmedian(data_sorted[section]).astype(np.float32)
        med_err[i] = np.nanmedian(err_sorted[section]).astype(np.float32)

    # Normalize contributing-pixel counts to fractional HEALPix coverage.
    healpix_area_sr = 4.0 * np.pi / (12.0 * HEALPIX17_NSIDE * HEALPIX17_NSIDE)
    expected_pixels = healpix_area_sr / pixel_area_sr
    covfrac = np.clip(counts.astype(np.float64) / expected_pixels, 0.0, 1.0).astype(
        np.float32
    )

    skyvals = np.empty(unique_hp.size, dtype=SKYVALS_DTYPE)
    skyvals["healpix17"] = unique_hp
    skyvals["data"] = med_data
    skyvals["err"] = med_err
    skyvals["covfrac"] = covfrac

    # Coarse (nside=2**11) coverage index list for quick downstream footprint use.
    healpix11_cov = np.unique(unique_hp // HEALPIX11_DIVISOR).astype(np.int64)

    return skyvals, healpix11_cov
