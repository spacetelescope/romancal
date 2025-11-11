"""
Module to estimate a 2D background and background RMS noise.
"""

import logging

import numpy as np
from astropy.stats import SigmaClip
from astropy.utils import lazyproperty
from photutils.background import Background2D, MedianBackground

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


class RomanBackground:
    """
    Class to estimate a 2D background and background RMS noise in an
    image.

    Parameters
    ----------
    data : 2D `~numpy.ndarray`
        The input 2D array.

    box_size : int or array_like (int)
        The box size along each axis.  If ``box_size`` is a scalar then
        a square box of size ``box_size`` will be used.  If ``box_size``
        has two elements, they should be in ``(ny, nx)`` order.

    mask : array_like (bool), optional
        A boolean mask, with the same shape as ``data``, where a `True`
        value indicates the corresponding element of ``data`` is
        masked. Masked data are excluded from calculations. ``mask``
        and ``coverage_mask`` differ only in that ``coverage_mask`` is
        applied to the output background and background RMS maps.

    coverage_mask : array_like (bool), optional
        A boolean mask, with the same shape as ``data``, where a `True`
        value indicates the corresponding element of ``data`` is masked.
        Masked data are excluded from calculations. ``coverage_mask``
        should be `True` where there is no coverage (i.e., no data) for
        a given pixel (e.g., blank areas in a mosaic image). It should
        not be used for bad pixels.

    Attributes
    ----------
    background : 2D `~numpy.ndimage`
        The estimated 2D background image.

    background_rms : 2D `~numpy.ndimage`
        The estimated 2D background RMS image.
    """

    def __init__(self, data, box_size=100, mask=None, coverage_mask=None):
        self.data = data
        self.box_size = np.asarray(box_size).astype(int)  # must be integer
        self.mask = mask
        self.coverage_mask = coverage_mask

    @lazyproperty
    def _background2d(self):
        """
        Estimate the 2D background and background RMS noise in an image.

        Returns
        -------
        background : `photutils.background.Background2D`
            A Background2D object containing the 2D background and
            background RMS noise estimates.
        """
        sigma_clip = SigmaClip(sigma=3.0)
        bkg_estimator = MedianBackground()
        filter_size = (3, 3)

        try:
            bkg = Background2D(
                self.data,
                self.box_size,
                filter_size=filter_size,
                mask=self.mask,
                coverage_mask=self.coverage_mask,
                sigma_clip=sigma_clip,
                bkg_estimator=bkg_estimator,
            )
        except ValueError:
            # use the entire unmasked array
            bkg = Background2D(
                self.data,
                self.data.shape,
                filter_size=filter_size,
                mask=self.mask,
                coverage_mask=self.coverage_mask,
                sigma_clip=sigma_clip,
                bkg_estimator=bkg_estimator,
                exclude_percentile=100.0,
            )
            log.info(
                "Background could not be estimated in meshes. "
                "Using the entire unmasked array for background "
                f"estimation: bkg_boxsize={self.data.shape}."
            )

        return bkg

    @lazyproperty
    def background(self):
        """
        The 2D background image.
        """
        return self._background2d.background

    @lazyproperty
    def background_rms(self):
        """
        The 2D background RMS image.
        """
        return self._background2d.background_rms
