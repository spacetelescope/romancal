"""
Module to calculate DAOFind sharpness and roundness1 properties.
"""

import warnings

import numpy as np
from astropy.convolution import Gaussian2DKernel
from astropy.nddata.utils import NoOverlapError, extract_array
from astropy.stats import gaussian_fwhm_to_sigma
from astropy.utils.decorators import lazyproperty
from scipy import ndimage


class DAOFindCatalog:
    """
    Class to calculate DAOFind ``sharpness`` and ``roundness1``
    statistics for a set of source positions.

    Parameters
    ----------
    data : 2D `~numpy.ndarray`
        The input 2D image data (already background-subtracted and on
        the same pixel grid as ``xypos``).

    xypos : 2D `~numpy.ndarray`
        Source positions with shape ``(N, 2)`` of ``(x, y)`` pixel
        coordinates.

    kernel_fwhm : float
        The full-width at half-maximum (FWHM) of the Gaussian kernel
        used for convolution, in pixels. The kernel size is derived from
        this value.

    requested_properties : iterable of str, optional
        If given, restrict ``self.properties`` to this subset of
        `available_properties` (``"sharpness"`` and/or ``"roundness1"``).

    Raises
    ------
    ValueError
        If ``data`` is not 2-D, ``xypos`` does not have shape ``(N, 2)``,
        or ``kernel_fwhm`` is not positive.
    """

    # Source-property column names this catalog can produce
    available_properties = ("sharpness", "roundness1")

    def __init__(self, data, xypos, kernel_fwhm, *, requested_properties=None):
        data = np.asarray(data)
        if data.ndim != 2:
            msg = f"data must be a 2-D array; got shape {data.shape}."
            raise ValueError(msg)

        xypos = np.asarray(xypos)
        if xypos.ndim != 2 or xypos.shape[1] != 2:
            msg = (
                "xypos must have shape (N, 2) of (x, y) pixel "
                f"positions; got shape {xypos.shape}."
            )
            raise ValueError(msg)

        if kernel_fwhm <= 0:
            msg = f"kernel_fwhm must be positive; got {kernel_fwhm}."
            raise ValueError(msg)

        self.data = data
        self.xypos = xypos
        self.kernel_sigma = kernel_fwhm * gaussian_fwhm_to_sigma

        if requested_properties is None:
            self.properties = list(self.available_properties)
        else:
            requested = set(requested_properties)
            self.properties = [
                prop for prop in self.available_properties if prop in requested
            ]

    @lazyproperty
    def kernel_size(self):
        """
        The DAOFind kernel size (in both x and y dimensions).
        """
        # Always odd
        return 2 * int(max(2.0, 1.5 * self.kernel_sigma)) + 1

    @lazyproperty
    def kernel_center(self):
        """
        The DAOFind kernel x/y center.
        """
        return (self.kernel_size - 1) // 2

    @lazyproperty
    def kernel_mask(self):
        """
        The DAOFind kernel circular mask.

        NOTE: 1=good pixels, 0=masked pixels
        """
        yy, xx = np.mgrid[0 : self.kernel_size, 0 : self.kernel_size]
        radius = np.hypot(xx - self.kernel_center, yy - self.kernel_center)
        return (radius <= max(2.0, 1.5 * self.kernel_sigma)).astype(int)

    @lazyproperty
    def kernel(self):
        """
        The DAOFind kernel, a 2D circular Gaussian normalized to have
        zero sum.
        """
        size = self.kernel_size
        kernel = Gaussian2DKernel(self.kernel_sigma, x_size=size, y_size=size).array
        kernel *= self.kernel_mask
        kernel /= np.max(kernel)

        # Normalize the kernel to zero sum
        npixels = self.kernel_mask.sum()
        denom = np.sum(kernel**2) - (np.sum(kernel) ** 2 / npixels)
        return ((kernel - (kernel.sum() / npixels)) / denom) * self.kernel_mask

    @lazyproperty
    def convolved_data(self):
        """
        The input data convolved with the DAOFind kernel.
        """
        return ndimage.convolve(self.data, self.kernel, mode="constant", cval=0.0)

    def _make_cutouts(self, array):
        """
        Make 2D cutouts centered on each source from the input data.

        Parameters
        ----------
        array : 2D `~numpy.ndarray`
            The input 2D image from which to make cutouts.

        Returns
        -------
        cutout : 3D `~numpy.ndarray`
            A 3D array containing 2D cutouts centered on each source from the
            input data. The cutout size always matches the size of the
            DAOFind kernel, which has odd dimensions.
        """
        cutouts = []
        for xcen, ycen in zip(*np.transpose(self.xypos), strict=True):
            try:
                cutout = extract_array(
                    array,
                    self.kernel.shape,
                    (ycen, xcen),
                    fill_value=0.0,
                )
            except NoOverlapError:
                cutout = np.zeros(self.kernel.shape)
            cutouts.append(cutout)

        return np.array(cutouts)  # all cutouts are the same size

    @lazyproperty
    def cutouts(self):
        """
        3D array containing 2D cutouts centered on each source from the
        input data.

        The cutout size always matches the size of the DAOFind kernel,
        which has odd dimensions.
        """
        return self._make_cutouts(self.data)

    @lazyproperty
    def cutouts_conv(self):
        """
        3D array containing 2D cutouts centered on each source from the
        input data convolved with the DAOFind kernel.

        The cutout size always matches the size of the DAOFind kernel,
        which has odd dimensions.
        """
        return self._make_cutouts(self.convolved_data)

    @lazyproperty
    def sharpness(self):
        """
        The DAOFind source sharpness statistic.

        The sharpness statistic measures the ratio of the difference
        between the height of the central pixel and the mean of the
        surrounding non-bad pixels to the height of the best fitting
        Gaussian function at that point (based on the convolved data).

        Stars generally have a "sharpness" between 0.2 and 1.0.
        """
        npixels = self.kernel_mask.sum() - 1  # exclude the peak pixel
        data_masked = self.cutouts * self.kernel_mask
        data_peak = self.cutouts[:, self.kernel_center, self.kernel_center]
        conv_peak = self.cutouts_conv[:, self.kernel_center, self.kernel_center]
        data_mean = (np.sum(data_masked, axis=(1, 2)) - data_peak) / npixels

        with warnings.catch_warnings():
            # Ignore 0 / 0 for non-finite xypos
            warnings.simplefilter("ignore", category=RuntimeWarning)
            value = (data_peak - data_mean) / conv_peak
            return value.astype(np.float32)

    @lazyproperty
    def roundness1(self):
        """
        The DAOFind source roundness1 statistic based on symmetry.

        The roundness1 statistic computes the ratio of a measure of the
        bilateral symmetry of the object to a measure of the four-fold
        symmetry of the object.

        "Round" objects have a "roundness1" close to 0, generally
        between -1 and 1.
        """
        # Set the central (peak) pixel to zero
        data = self.cutouts_conv.copy()
        data[:, self.kernel_center, self.kernel_center] = 0.0

        # Calculate the four roundness quadrants
        quad1 = data[:, 0 : self.kernel_center + 1, self.kernel_center + 1 :]
        quad2 = data[:, 0 : self.kernel_center, 0 : self.kernel_center + 1]
        quad3 = data[:, self.kernel_center :, 0 : self.kernel_center]
        quad4 = data[:, self.kernel_center + 1 :, self.kernel_center :]

        axis = (1, 2)
        sum2 = (
            -quad1.sum(axis=axis)
            + quad2.sum(axis=axis)
            - quad3.sum(axis=axis)
            + quad4.sum(axis=axis)
        )
        sum2[sum2 == 0] = 0.0

        sum4 = np.abs(data).sum(axis=axis)
        sum4[sum4 == 0] = np.nan

        with warnings.catch_warnings():
            # Ignore 0 / 0 for non-finite xypos
            warnings.simplefilter("ignore", category=RuntimeWarning)
            value = 2.0 * sum2 / sum4
            return value.astype(np.float32)
