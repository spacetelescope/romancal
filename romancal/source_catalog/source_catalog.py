"""
Module to calculate the source catalog.
"""

import logging
import warnings

import astropy.units as u
import numpy as np
from astropy.convolution import Gaussian2DKernel
from astropy.nddata.utils import NoOverlapError, extract_array
from astropy.stats import SigmaClip, gaussian_fwhm_to_sigma
from astropy.table import QTable, Table
from astropy.utils import lazyproperty
from astropy.utils.exceptions import AstropyUserWarning
from photutils.aperture import CircularAnnulus, CircularAperture, aperture_photometry
from photutils.segmentation import SourceCatalog
from roman_datamodels.datamodels import ImageModel, MosaicModel
from roman_datamodels.dqflags import pixel
from scipy import ndimage
from scipy.spatial import KDTree
from stpsf import __version__ as stpsf_version

from romancal import __version__ as romancal_version
from romancal.source_catalog import psf

from ._wcs_helpers import pixel_scale_angle_at_skycoord

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


class RomanSourceCatalog:
    """
    Class for the Roman source catalog.

    Parameters
    ----------
    model : `ImageModel` or `MosaicModel`
        The input data model. The image data is assumed to be background
        subtracted.

    segment_image : `~photutils.segmentation.SegmentationImage`
        A 2D segmentation image, with the same shape as the input data,
        where sources are marked by different positive integer values. A
        value of zero is reserved for the background.

    convolved_data : 2D `~numpy.ndarray` or `None`
        The 2D array used to calculate the source centroid and shape
        measurements. The image is assumed to be background subtracted.

    kernel_fwhm : float
        The full-width at half-maximum (FWHM) of the 2D Gaussian kernel.
        This is needed to calculate the DAOFind sharpness and roundness
        properties (DAOFind uses a special kernel that sums to zero).

    fit_psf : bool, optional
        Whether to fit a PSF model to the sources.

    mask : 2D `~numpy.ndarray` or `None`, optional
        A 2D boolean mask image with the same shape as the input data.
        This mask is used for PSF photometry. The mask should be the
        same one used to create the segmentation image.

    detection_cat : `None` or `~photutils.segmentation.SourceCatalog`, optional
        A `~photutils.segmentation.SourceCatalog` object for the
        detection image. The segmentation image used to create
        the detection catalog must be the same one input to
        ``segment_image``. If input, then the detection catalog source
        centroids and morphological/shape properties will be returned
        instead of calculating them from the input ``data``. The
        detection catalog centroids and shape properties will also be
        used to perform aperture photometry (i.e., circular and Kron).

    flux_unit : str, optional
        The unit of the flux density. Default is 'nJy'.

    Notes
    -----
    ``model.err`` is assumed to be the total error array corresponding
    to the input science ``model.data`` array. It is assumed to include
    *all* sources of error, including the Poisson error of the sources,
    and have the same shape and units as the science data array.
    """

    def __init__(
        self,
        model,
        segment_img,
        convolved_data,
        kernel_fwhm,
        *,
        fit_psf=True,
        mask=None,
        detection_cat=None,
        flux_unit="nJy",
    ):
        if not isinstance(model, ImageModel | MosaicModel):
            raise ValueError("The input model must be an ImageModel or MosaicModel.")

        self.model = model  # input model is background-subtracted
        self.segment_img = segment_img
        self.convolved_data = convolved_data
        self.kernel_sigma = kernel_fwhm * gaussian_fwhm_to_sigma
        self.fit_psf = fit_psf
        self.mask = mask
        self.detection_cat = detection_cat
        self.flux_unit = flux_unit

        self.n_sources = len(segment_img.labels)
        self.wcs = self.model.meta.wcs
        self.meta = {}

        # define flux unit conversion factors
        self.l2_to_sb = self.model.meta.photometry.conversion_megajanskys
        self.sb_to_flux = (1.0 * (u.MJy / u.sr) * self.pixel_area).to(
            u.Unit(self.flux_unit)
        )

        # needed for Kron photometry
        self.segm_sourcecat = None

    def __len__(self):
        return self.n_sources

    @lazyproperty
    def _pixscale_angle(self):
        """
        The pixel scale in arcseconds and the angle in degrees measured
        counterclockwise from the positive x axis to the "North" axis of
        the celestial coordinate system.

        Both are measured at the center of the image.
        """
        ysize, xsize = self.model.data.shape
        ycen = (ysize - 1) / 2.0
        xcen = (xsize - 1) / 2.0
        skycoord = self.wcs.pixel_to_world(xcen, ycen)
        _, pixscale, angle = pixel_scale_angle_at_skycoord(skycoord, self.wcs)
        return pixscale, angle

    @lazyproperty
    def _pixel_scale(self):
        """
        The pixel scale (as a Quantity in arcseconds) at the center of
        the image.
        """
        return self._pixscale_angle[0]

    @lazyproperty
    def _wcs_angle(self):
        """
        The angle (as a Quantity in degrees) measured counterclockwise
        from the positive x axis to the "North" axis of the celestial
        coordinate system.

        Measured at the center of the image.
        """
        return self._pixscale_angle[1]

    @lazyproperty
    def pixel_area(self):
        """
        The pixel area (as a Quantity in steradians).

        If the meta.photometry.pixel_area value is a negative
        placeholder (e.g., -999999 sr), the value is calculated from the
        WCS at the center of the image.
        """
        pixel_area = self.model.meta.photometry.pixel_area * u.sr
        if pixel_area < 0:
            pixel_area = (self._pixel_scale**2).to(u.sr)
        return pixel_area

    def convert_l2_to_sb(self):
        """
        Convert level-2 data from units of DN/s to MJy/sr (surface
        brightness).
        """
        # the conversion in done in-place to avoid making copies of the data;
        # use a dictionary to set the value to avoid on-the-fly validation
        self.model["data"] *= self.l2_to_sb
        self.model["err"] *= self.l2_to_sb
        if self.convolved_data is not None:
            self.convolved_data *= self.l2_to_sb

    def convert_sb_to_flux_density(self):
        """
        Convert level-3 data from units of MJy/sr (surface brightness)
        to flux density units.

        The flux density unit is defined by self.flux_unit.
        """
        # the conversion in done in-place to avoid making copies of the data;
        # use a dictionary to set the value to avoid on-the-fly validation
        self.model["data"] *= self.sb_to_flux.value
        self.model["data"] <<= self.sb_to_flux.unit
        self.model["err"] *= self.sb_to_flux.value
        self.model["err"] <<= self.sb_to_flux.unit
        if self.convolved_data is not None:
            self.convolved_data *= self.sb_to_flux.value
            self.convolved_data <<= self.sb_to_flux.unit

    def convert_flux_to_abmag(self, flux, flux_err):
        """
        Convert flux (and error) to AB magnitude (and error).

        Parameters
        ----------
        flux, flux_err : `~numpy.ndarray`
            The input flux and error arrays.

        Returns
        -------
        abmag, abmag_err : `~astropy.ndarray`
            The output AB magnitude and error arrays.
        """
        # ignore RunTimeWarning if flux or flux_err contains NaNs
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)

            abmag = flux.to(u.ABmag).value
            abmag_err = 2.5 * np.log10(1.0 + (flux_err / flux))

            # handle negative fluxes
            idx = flux.value < 0
            abmag[idx] = np.nan
            abmag_err[idx] = np.nan

        return abmag, abmag_err

    def calc_segment_properties(self):
        """
        Calculate the segment-based properties provided by
        `~photutils.segmentation.SourceCatalog`.

        The results are set as dynamic attributes on the class instance.
        """
        if self.detection_cat is not None:
            detection_cat = self.detection_cat.segm_sourcecat
        else:
            detection_cat = None

        segm_cat = SourceCatalog(
            self.model.data,
            self.segment_img,
            convolved_data=self.convolved_data,
            error=self.model.err,
            wcs=self.wcs,
            detection_cat=detection_cat,
        )
        self.segm_sourcecat = segm_cat
        self.meta.update(segm_cat.meta)

        # Extract the properties from the segment catalog. These
        # names are the SourceCatalog property names and the order
        # is not important.
        photutils_names = (
            "label",
            "xcentroid",
            "ycentroid",
            "sky_centroid",
            "bbox_xmin",
            "bbox_xmax",
            "bbox_ymin",
            "bbox_ymax",
            "area",
            "semimajor_sigma",
            "semiminor_sigma",
            "orientation",
            "ellipticity",
            "kron_radius",
            "segment_flux",
            "segment_fluxerr",
            "kron_flux",
            "kron_fluxerr",
        )

        # if needed, map names from photutils to the output catalog names
        name_map = {}
        name_map["xcentroid"] = "x_centroid"
        name_map["ycentroid"] = "y_centroid"
        name_map["area"] = "segment_area"
        name_map["orientation"] = "orientation_pix"
        name_map["segment_fluxerr"] = "segment_flux_err"
        name_map["kron_fluxerr"] = "kron_flux_err"

        # set the source properties as attributes of this instance
        for name in photutils_names:
            new_name = name_map.get(name, name)
            value = getattr(segm_cat, name)

            # change the photutils dtypes
            if new_name != "sky_centroid":
                if np.issubdtype(value.dtype, np.integer):
                    value = value.astype(np.int32)
                elif np.issubdtype(value.dtype, np.floating):
                    value = value.astype(np.float32)

            # handle any unit conversions
            if new_name in ("x_centroid", "y_centroid"):
                value *= u.pix
            if new_name == "segment_area":
                value = (value.value * self.pixel_area.to(u.arcsec**2)).astype(
                    np.float32
                )

            # split the sky_centroid values into separate RA and Dec
            # values
            if new_name == "sky_centroid":
                self.ra_centroid = value.ra
                self.dec_centroid = value.dec
            else:
                setattr(self, new_name, value)

    @lazyproperty
    def orientation_sky(self):
        """
        The position angle of the source major axis in degrees measured
        East of North.
        """
        return ((180.0 * u.deg) - self._wcs_angle + self.orientation_pix).astype(
            np.float32
        )

    @lazyproperty
    def _xypos(self):
        """
        The (x, y) source positions.

        If a detection catalog is input, then the detection catalog
        centroids are used. Otherwise, the segment catalog centroids are
        used.
        """
        if self.detection_cat is None:
            xycen = (self.x_centroid, self.y_centroid)
        else:
            xycen = (self.detection_cat.xcentroid, self.detection_cat.ycentroid)
        return np.transpose(xycen)

    @lazyproperty
    def _xypos_finite(self):
        """
        The (x, y) source positions where non-finite values have been
        replaced by to a large negative value.

        This is used with functions that fail with non-finite centroid
        values.

        For aperture photometry, at this position the aperture will not
        overlap the data, thus returning NaN fluxes and errors.
        """
        xypos = self._xypos.copy()
        nanmask = ~np.isfinite(xypos)
        xypos[nanmask] = -1000.0
        return xypos

    @lazyproperty
    def _xypos_nonfinite_mask(self):
        """
        A 1D boolean mask where `True` values denote sources where
        either the x_centroid or the y_centroid is not finite.
        """
        return ~np.isfinite(self._xypos).all(axis=1)

    @lazyproperty
    def aperture_radii(self):
        """
        A dictionary of the aperture radii used for aperture photometry.

        The radii are floats in units of pixels.
        """
        params = {}
        radii = np.array((0.1, 0.2, 0.4, 0.8, 1.6)) << u.arcsec
        annulus_radii = np.array((2.4, 2.8)) << u.arcsec
        params["circle"] = radii.copy()
        params["annulus"] = annulus_radii.copy()

        radii /= self._pixel_scale
        annulus_radii /= self._pixel_scale
        radii = radii.to(u.dimensionless_unscaled).value
        annulus_radii = annulus_radii.to(u.dimensionless_unscaled).value
        params["circle_pix"] = radii
        params["annulus_pix"] = annulus_radii

        return params

    @lazyproperty
    def _aperture_background(self):
        """
        The local background and error estimated using a circular
        annulus aperture.

        The local background is the sigma-clipped median value in the
        annulus. The background error is the standard error of the
        median, sqrt(pi / (2 * N)) * std.
        """
        bkg_aper = CircularAnnulus(
            self._xypos_finite,
            self.aperture_radii["annulus_pix"][0],
            self.aperture_radii["annulus_pix"][1],
        )
        bkg_aper_masks = bkg_aper.to_mask(method="center")
        sigclip = SigmaClip(sigma=3.0)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            warnings.simplefilter("ignore", category=AstropyUserWarning)

            nvalues = []
            bkg_median = []
            bkg_std = []
            for mask in bkg_aper_masks:
                bkg_data = mask.get_values(self.model.data)
                values = sigclip(bkg_data, masked=False)
                nvalues.append(values.size)
                med = np.median(values)
                std = np.std(values)
                if values.size == 0:
                    # handle case where source is completely masked due to
                    # forced photometry
                    med *= self.model.data.unit
                    std *= self.model.data.unit
                bkg_median.append(med)
                bkg_std.append(std)

            nvalues = np.array(nvalues)
            bkg_median = u.Quantity(bkg_median)
            bkg_std = u.Quantity(bkg_std)

            # standard error of the median
            bkg_median_err = np.sqrt(np.pi / (2.0 * nvalues)) * bkg_std

        return bkg_median.astype(np.float32), bkg_median_err.astype(np.float32)

    @lazyproperty
    def aper_bkg_flux(self):
        """
        The aperture local background flux (per pixel).
        """
        return self._aperture_background[0]

    @lazyproperty
    def aper_bkg_flux_err(self):
        """
        The aperture local background flux error (per pixel).
        """
        return self._aperture_background[1]

    @lazyproperty
    def aperture_flux_colnames(self):
        """
        The aperture flux column names.

        The column names are based on the circular aperture radii in
        tenths of arcsec. For example, the flux in a r=0.2 arcsec
        aperture is "aper02_flux".
        """
        flux_columns = []
        for radius in self.aperture_radii["circle"]:
            radius /= 0.1 * u.arcsec
            radius_str = np.round(radius.value).astype(int)
            flux_columns.append(f"aper{radius_str:02d}_flux")

        return flux_columns

    def calc_aperture_photometry(self):
        """
        Calculate the aperture photometry.

        The results are set as dynamic attributes on the class instance.
        """
        apertures = [
            CircularAperture(self._xypos_finite, radius)
            for radius in self.aperture_radii["circle_pix"]
        ]
        aper_phot = aperture_photometry(
            self.model.data, apertures, error=self.model.err
        )

        for i, aperture in enumerate(apertures):
            tmp_flux_col = f"aperture_sum_{i}"
            tmp_flux_err_col = f"aperture_sum_err_{i}"

            # subtract the local background measured in the annulus
            aper_areas = aperture.area_overlap(self.model.data)
            aper_phot[tmp_flux_col] -= self.aper_bkg_flux * aper_areas

            # set the flux and error attributes
            flux_col = self.aperture_flux_colnames[i]
            flux_err_col = f"{flux_col}_err"
            setattr(self, flux_col, aper_phot[tmp_flux_col].astype(np.float32))
            setattr(self, flux_err_col, aper_phot[tmp_flux_err_col].astype(np.float32))

    @lazyproperty
    def is_extended(self):
        """
        Boolean indicating whether the source is extended.

        Algorithm TBD.
        """
        return np.zeros(len(self), dtype=bool)

    @lazyproperty
    def image_flags(self):
        """
        Data quality bit flags (0 = good).

        Currently sets the bits:

        * 1:
          - L2: sources whose rounded centroid pixel is not finite or has
                DO_NOT_USE set in the model DQ
          - L3: sources whose rounded centroid pixel is not finite or has
                a weight of 0
        """
        xymask = np.isfinite(self._xypos[:, 0]) & np.isfinite(self._xypos[:, 1])
        xyidx = np.round(self._xypos[xymask, :]).astype(int)
        flags = np.full(self._xypos.shape[0], 0, dtype=np.int32)

        # sources whose centroid pixel is not finite
        flags[~xymask] = pixel.DO_NOT_USE

        try:
            # L2 images have a dq array
            dqflags = self.model.dq[xyidx[:, 1], xyidx[:, 0]]
            # if dqflags contains the DO_NOT_USE flag, set to DO_NOT_USE
            # (dq=1), otherwise 0
            flags[xymask] = dqflags & pixel.DO_NOT_USE

        except AttributeError:
            # L3 images
            mask = self.model.weight == 0
            flags[xymask] = mask[xyidx[:, 1], xyidx[:, 0]].astype(np.int32)

        return flags

    @lazyproperty
    def _daofind_kernel_size(self):
        """
        The DAOFind kernel size (in both x and y dimensions).
        """
        # always odd
        return 2 * int(max(2.0, 1.5 * self.kernel_sigma)) + 1

    @lazyproperty
    def _daofind_kernel_center(self):
        """
        The DAOFind kernel x/y center.
        """
        return (self._daofind_kernel_size - 1) // 2

    @lazyproperty
    def _daofind_kernel_mask(self):
        """
        The DAOFind kernel circular mask.

        NOTE: 1=good pixels, 0=masked pixels
        """
        yy, xx = np.mgrid[0 : self._daofind_kernel_size, 0 : self._daofind_kernel_size]
        radius = np.sqrt(
            (xx - self._daofind_kernel_center) ** 2
            + (yy - self._daofind_kernel_center) ** 2
        )
        return (radius <= max(2.0, 1.5 * self.kernel_sigma)).astype(int)

    @lazyproperty
    def _daofind_kernel(self):
        """
        The DAOFind kernel, a 2D circular Gaussian normalized to have
        zero sum.
        """
        size = self._daofind_kernel_size
        kernel = Gaussian2DKernel(self.kernel_sigma, x_size=size, y_size=size).array
        kernel *= self._daofind_kernel_mask
        kernel /= np.max(kernel)

        # normalize the kernel to zero sum
        npixels = self._daofind_kernel_mask.sum()
        denom = np.sum(kernel**2) - (np.sum(kernel) ** 2 / npixels)
        return ((kernel - (kernel.sum() / npixels)) / denom) * self._daofind_kernel_mask

    @lazyproperty
    def _daofind_convolved_data(self):
        """
        The DAOFind convolved data.
        """
        return ndimage.convolve(
            self.model.data.value, self._daofind_kernel, mode="constant", cval=0.0
        )

    @lazyproperty
    def _daofind_cutout(self):
        """
        3D array containing 2D cutouts centered on each source from the
        input data.

        The cutout size always matches the size of the DAOFind kernel,
        which has odd dimensions.
        """
        cutout = []
        for xcen, ycen in zip(*np.transpose(self._xypos_finite), strict=False):
            try:
                cutout_ = extract_array(
                    self.model.data,
                    self._daofind_kernel.shape,
                    (ycen, xcen),
                    fill_value=0.0,
                )
            except NoOverlapError:
                cutout_ = np.zeros(self._daofind_kernel.shape)
            cutout.append(cutout_)

        return np.array(cutout)  # all cutouts are the same size

    @lazyproperty
    def _daofind_cutout_conv(self):
        """
        3D array containing 2D cutouts centered on each source from the
        DAOFind convolved data.

        The cutout size always matches the size of the DAOFind kernel,
        which has odd dimensions.
        """
        cutout = []
        for xcen, ycen in zip(*np.transpose(self._xypos_finite), strict=False):
            try:
                cutout_ = extract_array(
                    self._daofind_convolved_data,
                    self._daofind_kernel.shape,
                    (ycen, xcen),
                    fill_value=0.0,
                )
            except NoOverlapError:
                cutout_ = np.zeros(self._daofind_kernel.shape)
            cutout.append(cutout_)

        return np.array(cutout)  # all cutouts are the same size

    @lazyproperty
    def sharpness(self):
        """
        The DAOFind source sharpness statistic.

        The sharpness statistic measures the ratio of the difference
        between the height of the central pixel and the mean of the
        surrounding non-bad pixels to the height of the best fitting
        Gaussian function at that point.

        Stars generally have a ``sharpness`` between 0.2 and 1.0.
        """
        npixels = self._daofind_kernel_mask.sum() - 1  # exclude the peak pixel
        data_masked = self._daofind_cutout * self._daofind_kernel_mask
        data_peak = self._daofind_cutout[
            :, self._daofind_kernel_center, self._daofind_kernel_center
        ]
        conv_peak = self._daofind_cutout_conv[
            :, self._daofind_kernel_center, self._daofind_kernel_center
        ]

        data_mean = (np.sum(data_masked, axis=(1, 2)) - data_peak) / npixels

        with warnings.catch_warnings():
            # ignore 0 / 0 for non-finite xypos
            warnings.simplefilter("ignore", category=RuntimeWarning)
            return (data_peak - data_mean) / conv_peak

    @lazyproperty
    def roundness(self):
        """
        The DAOFind source roundness statistic based on symmetry.

        The roundness characteristic computes the ratio of a measure of
        the bilateral symmetry of the object to a measure of the
        four-fold symmetry of the object.

        "Round" objects have a ``roundness`` close to 0, generally
        between -1 and 1.
        """
        # set the central (peak) pixel to zero
        cutout = self._daofind_cutout_conv.copy()
        cutout[:, self._daofind_kernel_center, self._daofind_kernel_center] = 0.0

        # calculate the four roundness quadrants
        quad1 = cutout[
            :, 0 : self._daofind_kernel_center + 1, self._daofind_kernel_center + 1 :
        ]
        quad2 = cutout[
            :, 0 : self._daofind_kernel_center, 0 : self._daofind_kernel_center + 1
        ]
        quad3 = cutout[
            :, self._daofind_kernel_center :, 0 : self._daofind_kernel_center
        ]
        quad4 = cutout[
            :, self._daofind_kernel_center + 1 :, self._daofind_kernel_center :
        ]

        axis = (1, 2)
        sum2 = (
            -quad1.sum(axis=axis)
            + quad2.sum(axis=axis)
            - quad3.sum(axis=axis)
            + quad4.sum(axis=axis)
        )
        sum2[sum2 == 0] = 0.0

        sum4 = np.abs(cutout).sum(axis=axis)
        sum4[sum4 == 0] = np.nan

        with warnings.catch_warnings():
            # ignore 0 / 0 for non-finite xypos
            warnings.simplefilter("ignore", category=RuntimeWarning)
            return 2.0 * sum2 / sum4

    @lazyproperty
    def _kdtree_query(self):
        """
        The distance in pixels to the nearest neighbor and its index.
        """
        if len(self) == 1:
            return [np.nan], [np.nan]

        # non-finite xypos causes memory errors on linux, but not MacOS
        tree = KDTree(self._xypos_finite)
        qdist, qidx = tree.query(self._xypos_finite, k=[2])
        return np.transpose(qdist)[0], np.transpose(qidx)[0]

    @lazyproperty
    def nn_label(self):
        """
        The label number of the nearest neighbor.

        A label value of -1 is returned if there is only one detected
        source and for sources with a non-finite centroid.
        """
        if len(self) == 1:
            return np.int32(-1)

        nn_label = self.label[self._kdtree_query[1]].astype(np.int32)
        # assign a label of -1 for non-finite xypos
        nn_label[self._xypos_nonfinite_mask] = -1

        return nn_label

    @lazyproperty
    def nn_dist(self):
        """
        The distance in arcsec to the nearest neighbor.

        NaN is returned for non-finite centroid positions or when
        the catalog contains only one source.
        """
        nn_dist = self._kdtree_query[0]
        if len(self) != 1:
            # assign a distance of np.nan for non-finite xypos
            nn_dist[self._xypos_nonfinite_mask] = np.nan
        return (nn_dist * self._pixel_scale).astype(np.float32)

    def _update_metadata(self):
        """
        Update the metadata dictionary with the package version
        information and aperture parameters.
        """
        ver_key = "versions"
        if "version" in self.meta:
            # depends on photutils version
            self.meta[ver_key] = self.meta.pop("version")

        ver_dict = self.meta.get("versions", None)
        if ver_dict is not None:
            ver_dict["romancal"] = romancal_version
            ver_dict["stpsf"] = stpsf_version
            packages = [
                "Python",
                "numpy",
                "scipy",
                "astropy",
                "photutils",
                "gwcs",
                "stpsf",
                "romancal",
            ]
            ver_dict = {key: ver_dict[key] for key in packages if key in ver_dict}
            self.meta[ver_key] = ver_dict

        # reorder the metadata dictionary
        self.meta.pop("localbkg_width")
        keys = ["date", "versions", "apermask_method", "kron_params"]
        old_meta = self.meta.copy()
        self.meta = {key: old_meta[key] for key in keys}

        # reformat the aperture radii for the metadata to remove
        # Quantity objects
        aper_radii = self.aperture_radii.copy()
        aper_radii["circle_arcsec"] = aper_radii.pop("circle").value
        aper_radii["annulus_arcsec"] = aper_radii.pop("annulus").value
        self.meta["aperture_radii"] = aper_radii

    @lazyproperty
    def psf_model(self):
        """
        A gridded PSF model based on instrument and detector
        information.

        The `~photutils.psf.GriddedPSF` model is created using the
        STPSF library.
        """
        log.info("Constructing a gridded PSF model.")
        if hasattr(self.model.meta, "instrument"):
            # ImageModel (L2 datamodel)
            filt = self.model.meta.instrument.optical_element
            detector = self.model.meta.instrument.detector.replace("WFI", "SCA")
        else:
            # MosaicModel (L3 datamodel)
            filt = self.model.meta.basic.optical_element
            detector = "SCA02"

        gridded_psf_model, _ = psf.create_gridded_psf_model(
            filt=filt,
            detector=detector,
        )

        return gridded_psf_model

    def calc_psf_photometry(self) -> None:
        """
        Perform PSF photometry by fitting PSF models to detected sources
        for refined astrometry.
        """
        log.info("Fitting a PSF model to sources for improved astrometric precision.")
        xinit, yinit = np.transpose(self._xypos)
        psf_photometry_table, _ = psf.fit_psf_to_image_model(
            image_model=self.model,
            mask=self.mask,
            psf_model=self.psf_model,
            x_init=xinit,
            y_init=yinit,
            exclude_out_of_bounds=True,
        )

        # map photutils column names to the output catalog names
        name_map = {}
        name_map["flags"] = "psf_flags"
        name_map["x_fit"] = "x_psf"
        name_map["x_err"] = "x_psf_err"
        name_map["y_fit"] = "y_psf"
        name_map["y_err"] = "y_psf_err"
        name_map["flux_fit"] = "psf_flux"
        name_map["flux_err"] = "psf_flux_err"

        # set these columns as attributes of this instance
        for old_name, new_name in name_map.items():
            setattr(self, new_name, psf_photometry_table[old_name])

    @lazyproperty
    def catalog_descriptions(self):
        """
        A dictionary of the output catalog column descriptions.

        The order is not important.
        """
        col = {}
        col["label"] = "Label of the source segment in the segmentation image"
        col["x_centroid"] = (
            "Column coordinate of the source centroid in the detection "
            "image from image moments (0 indexed)"
        )
        col["y_centroid"] = (
            "Row coordinate of the source centroid in the detection image "
            "from image moments (0 indexed)"
        )
        col["ra_centroid"] = "Right ascension (ICRS) of the image centroid"
        col["dec_centroid"] = "Declination (ICRS) of the image centroid"
        col["semimajor_sigma"] = (
            "1-sigma standard deviation along the semimajor axis of the 2D Gaussian function that has the same second-order central moments as the source"
        )
        col["semiminor_sigma"] = (
            "1-sigma standard deviation along the semiminor axis of the 2D Gaussian function that has the same second-order central moments as the source"
        )
        col["ellipticity"] = "1 - (semimajor_sigma / semiminor_sigma)"
        col["orientation_pix"] = (
            "The angle measured counter-clockwise from the positive X axis to the major axis computed from image moments"
        )
        col["orientation_sky"] = (
            "The position angle from North of the major axis computed from "
            "image moments"
        )

        col["segment_flux"] = "Isophotal flux"
        col["segment_flux_err"] = "Uncertainty in segment_flux"
        col["segment_area"] = "Area of the source segment"
        col["kron_flux"] = "Flux within the elliptical Kron aperture"
        col["kron_flux_err"] = "Uncertainty in kron_flux"
        col["aper_bkg_flux"] = "The local background estimate for aperture photometry"
        col["aper_bkg_flux_err"] = "Uncertainty in aper_bkg_flux"

        for i, colname in enumerate(self.aperture_flux_colnames):
            desc = (
                "within a circular aperture of radius="
                f"{self.aperture_radii['circle'][i]:0.1f}"
            )
            col[colname] = f"Flux {desc}"
            col[f"{colname}_err"] = f"Uncertainty in {colname}"

        col["x_psf"] = "Column position of the source from PSF fitting (0 indexed)"
        col["x_psf_err"] = "Uncertainty in x_psf"
        col["y_psf"] = "Row position of the source from PSF fitting (0 indexed)"
        col["y_psf_err"] = "Uncertainty in y_psf"
        col["psf_flux"] = "Total PSF flux"
        col["psf_flux_err"] = "Uncertainty in psf_flux"
        col["psf_flags"] = "PSF fitting bit flags"

        col["image_flags"] = "Data quality bit flags"
        col["is_extended"] = "Flag indicating whether the source is extended"
        col["sharpness"] = "The DAOFind source sharpness statistic"
        col["roundness"] = "The DAOFind source roundness statistic"
        col["nn_label"] = "The label number of the nearest neighbor in this skycell"
        col["nn_dist"] = "The distance to the nearest neighbor in this skycell"

        return col

    @lazyproperty
    def catalog_colnames(self):
        """
        An ordered list of the output catalog column names.
        """
        aper_colnames = []
        for colname in self.aperture_flux_colnames:
            aper_colnames.append(colname)
            aper_colnames.append(f"{colname}_err")

        if self.detection_cat is None:
            colnames = [
                "label",
                "x_centroid",
                "y_centroid",
                "ra_centroid",
                "dec_centroid",
                "aper_bkg_flux",
                "aper_bkg_flux_err",
            ]

            colnames.extend(aper_colnames)

            colnames2 = [
                "image_flags",
                "is_extended",
                "sharpness",
                "roundness",
                "nn_label",
                "nn_dist",
                "segment_flux",
                "segment_flux_err",
                "segment_area",
                "kron_flux",
                "kron_flux_err",
                "semimajor_sigma",
                "semiminor_sigma",
                "ellipticity",
                "orientation_pix",
                "orientation_sky",
            ]
            colnames.extend(colnames2)

        else:
            colnames = [
                "label",
                "image_flags",
                "is_extended",
                "sharpness",
                "roundness",
                "segment_flux",
                "segment_flux_err",
                "kron_flux",
                "kron_flux_err",
            ]
            colnames.extend(aper_colnames)

        if self.fit_psf:
            psf_colnames = [
                "psf_flags",
                "x_psf",
                "x_psf_err",
                "y_psf",
                "y_psf_err",
                "psf_flux",
                "psf_flux_err",
            ]
            colnames.extend(psf_colnames)

        return colnames

    @lazyproperty
    def catalog(self):
        """
        The final source catalog.
        """
        # convert data to flux units
        if isinstance(self.model, ImageModel):
            self.convert_l2_to_sb()
        self.convert_sb_to_flux_density()

        # make measurements
        self.calc_segment_properties()
        self.calc_aperture_photometry()
        if self.fit_psf:
            self.calc_psf_photometry()

        # put the measurements into a Table
        catalog = QTable()
        for column in self.catalog_colnames:
            catalog[column] = getattr(self, column)
            descrip = self.catalog_descriptions.get(column, None)
            catalog[column].info.description = descrip
        self._update_metadata()
        catalog.meta.update(self.meta)

        # convert QTable to Table to avoid having Quantity columns
        catalog = Table(catalog)

        return catalog
