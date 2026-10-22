"""Effective Point Spread Function utilities.

Utilities to create ePSF models to be used where PSFs are needed in source detection
and flux/position analysis
"""
from astropy.nddata import NDData
from astropy.stats import sigma_clipped_stats
from astropy.table import Table
import numpy as np
from photutils.psf import EPSFBuilder, EPSFStars, extract_stars

from roman_datamodels import dqflags


def build_epsf(stars, max_psferr=None, **kwargs_epsfbuilder):
    """Create the ePSF from a list of stars.

    Parameters
    ----------
    stars : EPSFStars
        The stars to base the ePSF on

    kwargs_epsfbuilder : dict
        Arguments to pass to EPSFBuilder. If the following arguments
        are left unspecified, these defaults will be used:
        - oversample=4
        - maxiters=30
        - central_accuracy: Value of argument `max_psferr` if max_psferr is not None.
        - progress_bar=false

    Returns
    -------
    epsf_result : EPFSBuildResults
        The created ePSF.
    """
    epsfbuilder_final_kwargs = {
            'oversampling': 4, 'maxiters': 30, 'progress_bar': False
        }
    if kwargs_epsfbuilder is None:
        kwargs_epsfbuilder = dict()
    epsfbuilder_final_kwargs.update(kwargs_epsfbuilder)
    if max_psferr is not None:
        epsfbuilder_final_kwargs['central_accuracy'] = epsfbuilder_final_kwargs.get('central_accuracy', max_psferr)

    epsfbuilder = EPSFBuilder(**epsfbuilder_final_kwargs)
    epsf_result = epsfbuilder(stars)

    return epsf_result


def pick_stars(model, catalog, stamp_size=9, mask_nans=True,
               mask_dq=(dqflags.pixel.DO_NOT_USE + dqflags.pixel.SATURATED),
               min_flux=500., min_stdev=None, max_psferr=None,
               mask_extended=True, mask_failed_psffit=True,
              ):
    """Pick stars out of an image based on various criteria.

    Parameters
    ----------
    model : roman_datamodels.datamodels.ImageModel
        L2 model

    catalog : astropy.table.Table
        catalog based on model

    stamp_size : int
        Size of the cutouts

    mask_dq : int or None
        Pixel flags to count as bad

    min_flux : float or None
        Minimum flux.

    min_stdev : float or None
        Minimum standard deviation with a stamp to be valid.

    max_psferr : float or None
        Maximum  sqrt(x_psf_err**2 + y_psf_err**2) allowed.
        If None, sigma-clipped stats to remove anything larger that
        two times the deviation.

    mask_extended : bool
        Remove extended sources

    mask_failed_psffit : bool
        Removed objects that have a bad PSF fit. Indicated by the x/y_psf_err being 0.

    Returns
    -------
    stars : EPSFStars
        The stars
    """
    # Massage the data
    mean_val, median_val, std_val = sigma_clipped_stats(model.data, sigma=2.0)
    data = model.data - median_val

    # Create simplified table
    hsize = (stamp_size - 1) / 2
    x = catalog['x_centroid']
    y = catalog['y_centroid']
    is_extended = catalog['is_extended']
    label = catalog['label']
    psf_err = np.sqrt(catalog['x_psf_err']**2 + catalog['y_psf_err']**2)
    mask = ((x > hsize) & (x < (data.shape[1] - 1 - hsize)) &
            (y > hsize) & (y < (data.shape[0] - 1 - hsize)))
    if mask_extended:
        mask = [m & (not ie) for m, ie in zip(mask, is_extended, strict=True)]
    if mask_failed_psffit:
        mask = [m & (psferr != 0.) for m, psferr in zip(mask, psf_err, strict=True)]

    stars_tbl = Table()
    stars_tbl['x'] = x[mask]
    stars_tbl['y'] = y[mask]
    stars_tbl['label'] = label[mask]
    stars_tbl['psf_err'] = psf_err[mask]

    # Extract stars
    stars = extract_stars(NDData(data), stars_tbl, size=stamp_size)
    for star, label, psf_err in zip(stars, stars_tbl['label'], stars_tbl['psf_err'], strict=True):
        star.label = label
        star.psf_err = psf_err
        star.stats = sigma_clipped_stats(star.data, sigma=2.0)

    # Throw out stars with bad values
    if mask_nans:
        ignore_mask = [np.any(np.isnan(star))
                       for star in stars
                      ]
    else:
        ignore_mask = [False] * len(stars)

    # Throw out stars with bad dq flags
    if mask_dq is not None:
        def dq_mask(dq, star, mask_dq):
            bb = star.bbox
            dq = dq[bb.iymin:bb.iymax, bb.ixmin:bb.ixmax]
            mask = (dq & mask_dq).any()
            return mask

        ignore_mask = [mask or dq_mask(model.dq, star, mask_dq)
                       for mask, star in zip(ignore_mask, stars, strict=True)
                      ]

    # Only bright
    if min_flux is not None:
        ignore_mask = [mask or star.estimate_flux() < min_flux
                       for mask, star in zip(ignore_mask, stars, strict=True)
                      ]

    # Only significant detections
    if min_stdev is not None:
        ignore_mask = [mask or star.stats[2] < min_stdev
                       for mask, star in zip(ignore_mask, stars, strict=True)
                      ]

    # With all gross undesirable detections removed,
    # make a final sweep using the psf error.
    if max_psferr is None:
        all_psf_errs = [star.psf_err
                        for mask, star in zip(ignore_mask, stars, strict=True)
                        if not mask]
        _, median, stddev = sigma_clipped_stats(all_psf_errs)
        max_psferr = median + (2 * stddev)
    ignore_mask = [mask or star.psf_err > max_psferr
                   for mask, star in zip(ignore_mask, stars, strict=True)]

    # Filter the stars based on the ignore mask
    stars = EPSFStars([star
                       for mask, star in zip(ignore_mask, stars, strict=True)
                       if not mask
                      ])

    # That's all folks
    return stars, max_psferr
