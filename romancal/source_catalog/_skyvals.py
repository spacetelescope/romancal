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
# Default row-chunk size for WCS/HEALPix conversion (~1M pixels on a 4k-wide SCA).
_DEFAULT_CHUNK_ROWS = 256
SKYVALS_DTYPE = np.dtype(
    [
        ("healpix17", np.int64),
        ("data", np.float32),
        ("err", np.float32),
        ("covfrac", np.float32),
    ]
)


def compute_skyvals(input_model, segmentation, bad_pixel_mask, chunk_rows=_DEFAULT_CHUNK_ROWS):
    """
    Compute HEALPix-aggregated sky statistics from masked image pixels.

    Coordinates are converted to HEALPix in row chunks of ``chunk_rows`` to
    bound temporary WCS/HEALPix arrays on large images.

    Parameters
    ----------
    input_model
        Roman image model containing data/err arrays and metadata.
    segmentation : np.ndarray
        2D segmentation labels; source pixels are non-zero.
    bad_pixel_mask : np.ndarray
        Boolean 2D mask where True means "exclude pixel from sky calculation".
    chunk_rows : int, optional
        Number of image rows to convert through WCS/HEALPix at a time.
        Must be >= 1. Default is 256.

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
    if chunk_rows < 1:
        raise ValueError(f"chunk_rows must be >= 1, got {chunk_rows}")

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

    # Valid sky pixels: unmasked and not part of a detected source.
    valid = coverage & (segmentation == 0)
    if not np.any(valid):
        # If every unmasked pixel is in a source, fall back to all unmasked
        # pixels so per-exposure sky estimates are still produced.
        valid = coverage

    ny = data.shape[0]
    n_sky = int(np.count_nonzero(valid))
    healpix17 = np.empty(n_sky, dtype=np.int64)
    vals_data = np.empty(n_sky, dtype=np.float32)
    vals_err = np.empty(n_sky, dtype=np.float32)
    n_filled = 0
    hp11_chunks: list[np.ndarray] = []

    healpix = HEALPix(nside=HEALPIX17_NSIDE, order="nested", frame="icrs")

    for y0 in range(0, ny, chunk_rows):
        y1 = min(y0 + chunk_rows, ny)
        cov_rows = coverage[y0:y1]
        if not np.any(cov_rows):
            continue

        yy_local, xx = np.nonzero(cov_rows)
        yy = yy_local + y0

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
            yy_sky = yy[is_sky_pixel]
            xx_sky = xx[is_sky_pixel]
            healpix17[sl] = chunk_healpix17[is_sky_pixel]
            vals_data[sl] = data[yy_sky, xx_sky]
            vals_err[sl] = err[yy_sky, xx_sky]
            n_filled += n

    if n_filled != n_sky:
        raise RuntimeError(
            f"skyvals fill count mismatch: filled {n_filled} sky pixels, "
            f"expected {n_sky}"
        )

    if hp11_chunks:
        healpix11_cov = np.unique(np.concatenate(hp11_chunks)).astype(np.int64)
    else:
        healpix11_cov = np.empty((0,), dtype=np.int64)

    # argsort once, then gather values into healpix order. Avoids keeping both
    # unsorted and fully duplicated sorted buffers longer than needed.
    order = np.argsort(healpix17, kind="mergesort")
    hp_sorted = np.take(healpix17, order)
    del healpix17
    unique_hp, starts, counts = np.unique(
        hp_sorted, return_index=True, return_counts=True
    )
    del hp_sorted

    vals_data = np.take(vals_data, order)
    vals_err = np.take(vals_err, order)
    del order

    med_data = np.empty(unique_hp.size, dtype=np.float32)
    med_err = np.empty(unique_hp.size, dtype=np.float32)
    # bad_pixel_mask already excludes non-finite data/err.
    for i, (start, count) in enumerate(zip(starts, counts, strict=True)):
        section = slice(start, start + count)
        med_data[i] = np.median(vals_data[section])
        med_err[i] = np.median(vals_err[section])

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
