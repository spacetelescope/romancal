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
# Process this many image rows at a time through WCS/HEALPix. Full-image
# conversion of a Roman SCA (~16M pixels) peaks several GB above baseline;
# chunking keeps temporary coordinate arrays small while preserving results.
_SKYVALS_CHUNK_ROWS = 256
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

    Returns
    -------
    skyvals : np.ndarray
        Structured array with one row per NEST-ordered HEALPix cell at
        ``nside=2**17``. Its ``healpix17``, ``data``, ``err``, and
        ``covfrac`` fields contain the cell index, median image value, median
        error value, and image-footprint coverage fraction, respectively.
    healpix11_cov : np.ndarray
        Unique NEST-ordered HEALPix cell indices at ``nside=2**11`` that
        contain at least one pixel not rejected by the bad-pixel mask,
        regardless of its segmentation label.
    """
    data = input_model.data
    err = input_model.err
    wcs = input_model.meta.wcs
    pixel_area_sr = get_pixel_area_sr(input_model)

    # Footprint coverage includes all usable pixels, including source pixels.
    coverage = ~bad_pixel_mask
    if not np.any(coverage):
        log.warning(
            "No usable pixels for skyvals after bad-pixel masking; "
            "returning empty skyvals/healpix11_cov arrays."
        )
        return np.empty((0,), dtype=SKYVALS_DTYPE), np.empty((0,), dtype=np.int64)

    # Valid sky pixels are those not rejected by upstream masks and not part
    # of detected sources.
    valid = coverage & (segmentation == 0)
    if not np.any(valid):
        # Fallback: if every source-masked pixel is rejected, use all unmasked
        # pixels so we still produce per-exposure sky estimates.
        valid = coverage

    ny = data.shape[0]
    n_sky = int(np.count_nonzero(valid))
    # Preallocate flat sky-pixel buffers. Peak memory is dominated by these
    # plus per-chunk WCS temporaries, not by a second full-image copy of RA/Dec.
    healpix17 = np.empty(n_sky, dtype=np.int64)
    vals_data = np.empty(n_sky, dtype=np.float32)
    vals_err = np.empty(n_sky, dtype=np.float32)
    n_filled = 0
    hp11_chunks: list[np.ndarray] = []

    healpix = HEALPix(nside=HEALPIX17_NSIDE, order="nested", frame="icrs")

    for y0 in range(0, ny, _SKYVALS_CHUNK_ROWS):
        y1 = min(y0 + _SKYVALS_CHUNK_ROWS, ny)
        cov_rows = coverage[y0:y1]
        if not np.any(cov_rows):
            continue

        yy_local, xx = np.nonzero(cov_rows)
        yy = yy_local + y0

        # Convert this row chunk only — full-image WCS was the main memory spike.
        ra, dec = wcs.pixel_to_world_values(xx, yy)
        chunk_healpix17 = np.asarray(
            healpix.lonlat_to_healpix(ra * u.deg, dec * u.deg),
            dtype=np.int64,
        )
        hp11_chunks.append(np.unique(chunk_healpix17 // HEALPIX11_DIVISOR))

        is_sky_pixel = valid[yy, xx]
        if np.any(is_sky_pixel):
            n = int(np.count_nonzero(is_sky_pixel))
            sl = slice(n_filled, n_filled + n)
            healpix17[sl] = chunk_healpix17[is_sky_pixel]
            yy_sky = yy[is_sky_pixel]
            xx_sky = xx[is_sky_pixel]
            vals_data[sl] = data[yy_sky, xx_sky]
            vals_err[sl] = err[yy_sky, xx_sky]
            n_filled += n

    if n_filled != n_sky:
        # Defensive trim if counts drift (should not happen with consistent masks).
        healpix17 = healpix17[:n_filled]
        vals_data = vals_data[:n_filled]
        vals_err = vals_err[:n_filled]

    if hp11_chunks:
        healpix11_cov = np.unique(np.concatenate(hp11_chunks)).astype(np.int64)
    else:
        healpix11_cov = np.empty((0,), dtype=np.int64)

    if n_filled == 0:
        return np.empty((0,), dtype=SKYVALS_DTYPE), healpix11_cov

    # Sort/group by healpixel so per-healpixel reductions are contiguous.
    order = np.argsort(healpix17)
    hp_sorted = healpix17[order]
    data_sorted = vals_data[order]
    err_sorted = vals_err[order]
    del healpix17, vals_data, vals_err, order

    unique_hp, starts, counts = np.unique(
        hp_sorted, return_index=True, return_counts=True
    )

    med_data = np.empty(unique_hp.size, dtype=np.float32)
    med_err = np.empty(unique_hp.size, dtype=np.float32)
    # Pixels entering this path are already finite (bad_pixel_mask excludes
    # non-finite data/err), so plain median is sufficient and cheaper than
    # nanmedian.
    for i, (start, count) in enumerate(zip(starts, counts, strict=True)):
        section = slice(start, start + count)
        med_data[i] = np.median(data_sorted[section])
        med_err[i] = np.median(err_sorted[section])

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

    return skyvals, healpix11_cov
