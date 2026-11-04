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

    exclude_margin : float or `None`, optional
        Boxes are excluded from the background mesh, and interpolated
        over, when more than ``exclude_percentile`` percent of their
        pixels are masked or sigma clipped.  ``exclude_percentile`` is
        the percentage of unmasked pixels sigma clipped in the median
        box plus ``exclude_margin``, and at least 10 (the photutils
        default).  In crowded fields, where stars are clipped from every
        box, this keeps the boxes that are typical of the field while
        still excluding those contaminated by a large source.  If
        `None`, ``exclude_percentile`` is 10.

    Attributes
    ----------
    background : 2D `~numpy.ndarray`
        The estimated 2D background image.

    background_rms : 2D `~numpy.ndarray`
        The estimated 2D background RMS image.
    """

    def __init__(
        self, data, box_size=100, mask=None, coverage_mask=None, exclude_margin=10.0
    ):
        self.data = data
        self.box_size = np.asarray(box_size).astype(int)  # must be integer
        self.mask = mask
        self.coverage_mask = coverage_mask
        self.exclude_margin = exclude_margin

    def _exclude_percentile(self, **kwargs):
        """
        The percentage of masked or clipped pixels above which a box is
        excluded: the clipped percentage of the median box plus
        ``exclude_margin``, and at least 10.
        """
        min_percentile = 10.0
        if self.exclude_margin is None:
            return min_percentile

        try:
            bkg = Background2D(
                self.data, self.box_size, exclude_percentile=100.0, **kwargs
            )
        except ValueError:
            return min_percentile

        # unmasked pixels per box, on the same boxes as Background2D
        masked = ~np.isfinite(self.data)
        for mask in (self.mask, self.coverage_mask):
            if mask is not None:
                masked |= mask
        unmasked = (~masked).astype(float)
        # reduceat sums between consecutive start indices (the last runs to
        # the end), giving per-box counts on Background2D's mesh, including
        # its partial boxes at the top and right edges
        box_ny, box_nx = np.broadcast_to(self.box_size, 2)
        edges_y = np.arange(0, self.data.shape[0], box_ny)
        edges_x = np.arange(0, self.data.shape[1], box_nx)
        nunmasked = np.add.reduceat(
            np.add.reduceat(unmasked, edges_y, axis=0), edges_x, axis=1
        )

        nclipped = nunmasked - bkg.n_pixels_mesh
        good = nunmasked > 0
        clipped_percentile = 100.0 * np.median(nclipped[good] / nunmasked[good])
        percentile = min(
            max(min_percentile, clipped_percentile + self.exclude_margin), 100.0
        )
        log.info(
            f"Sigma clipping removes {clipped_percentile:.1f}% of the median "
            f"background box; using exclude_percentile={percentile:.1f}."
        )
        return percentile

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
        kwargs = {
            "filter_size": (3, 3),
            "mask": self.mask,
            "coverage_mask": self.coverage_mask,
            "sigma_clip": SigmaClip(sigma=3.0, stdfunc="mad_std"),
            "bkg_estimator": MedianBackground(),
        }

        try:
            bkg = Background2D(
                self.data,
                self.box_size,
                exclude_percentile=self._exclude_percentile(**kwargs),
                **kwargs,
            )
        except ValueError:
            # Use the entire unmasked array
            bkg = Background2D(
                self.data, self.data.shape, exclude_percentile=100.0, **kwargs
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
