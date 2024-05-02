"""
Module to calculate the source catalog.
"""

import logging
import warnings

import astropy.units as u
import numpy as np
from astropy.convolution import Gaussian2DKernel
from astropy.coordinates import SkyCoord
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

from romancal import __version__ as romancal_version

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

    convolved_data : data : 2D `~numpy.ndarray`
        The 2D array used to calculate the source centroid and shape
        measurements.

    aperture_params : `dict`
        A dictionary containing the parameters (radii, aperture
        corrections, and background annulus inner and outer radii) used
        to perform aperture photometry.

    ci_star_thresholds : array-like of 2 floats
        The concentration index (CI) thresholds for determining whether
        a source is a star. The first threshold corresponds to the
        concentration index calculated from the smallest and middle
        aperture radii (see ``aperture_params``). The second threshold
        corresponds to the concentration index calculated from the
        middle and largest aperture radii. An object is considered
        extended if both concentration indices are greater than the
        corresponding thresholds, otherwise it is considered a star.

    kernel_fwhm : float
        The full-width at half-maximum (FWHM) of the 2D Gaussian kernel.
        This is needed to calculate the DAOFind sharpness and roundness
        properties (DAOFind uses a special kernel that sums to zero).

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
        aperture_params,
        ci_star_thresholds,
        kernel_fwhm,
        flux_unit="uJy",
    ):

        if not isinstance(model, (ImageModel, MosaicModel)):
            raise ValueError("The input model must be an ImageModel or MosaicModel.")
        self.model = model  # background was previously subtracted

        self.segment_img = segment_img
        self.convolved_data = convolved_data
        self.aperture_params = aperture_params
        self.kernel_sigma = kernel_fwhm * gaussian_fwhm_to_sigma
        self.flux_unit = flux_unit

        self.sb_unit = "MJy/sr"
        self.l2_unit = "DN/s"
        self.l2_conv_factor = (
            self.model.meta.photometry.conversion_megajanskys / self.l2_unit
        )

        if len(ci_star_thresholds) != 2:
            raise ValueError("ci_star_thresholds must contain only 2 items")
        self.ci_star_thresholds = ci_star_thresholds

        self.n_sources = len(self.segment_img.labels)
        self.aperture_ee = self.aperture_params["aperture_ee"]
        self.n_aper = len(self.aperture_ee)
        self.wcs = self.model.meta.wcs
        self.column_desc = {}
        self.meta = {}

    @lazyproperty
    def _pixscale_angle(self):
        ysize, xsize = self.model.data.shape
        ycen = (ysize - 1) / 2.0
        xcen = (xsize - 1) / 2.0
        skycoord = self.wcs.pixel_to_world(xcen, ycen)
        _, pixscale, angle = pixel_scale_angle_at_skycoord(skycoord, self.wcs)
        return pixscale, angle

    @lazyproperty
    def _pixel_scale(self):
        """
        The pixel scale in arcseconds at the center of the image.
        """
        return self._pixscale_angle[0]

    @lazyproperty
    def _wcs_angle(self):
        """
        The angle (in degrees) measured counterclockwise from the
        positive x axis to the "North" axis of the celestial coordinate
        system.

        Measured at the center of the image.
        """
        return self._pixscale_angle[1]

    @lazyproperty
    def pixel_area(self):
        """
        The pixel area in steradians.

        If the value is a negative placeholder (e.g., -999999 sr), the
        value is calculated from the WCS at the center of the image.
        """
        pixel_area = self.model.meta.photometry.pixelarea_steradians
        if pixel_area < 0:
            pixel_area = (self._pixel_scale**2).to(u.sr)
        return pixel_area

    def convert_l2_to_sb(self):
        """
        Convert a level-2 image from units of DN/s to MJy/sr (surface
        brightness).
        """
        if self.model.data.unit != self.l2_unit or self.model.err.unit != self.l2_unit:
            raise ValueError(
                f"data and err are expected to be in units of {self.l2_unit}"
            )

        # the conversion in done in-place to avoid making copies of the data;
        # use a dictionary to set the value to avoid on-the-fly validation
        self.model["data"] *= self.l2_conv_factor
        self.model["err"] *= self.l2_conv_factor
        self.convolved_data *= self.l2_conv_factor

    def convert_sb_to_l2(self):
        """
        Convert the data and error Quantity arrays from MJy/sr (surface
        brightness) to DN/s (level-2 units).

        This is the inverse operation of `convert_l2_to_sb`.
        """
        if self.model.data.unit != self.sb_unit or self.model.err.unit != self.sb_unit:
            raise ValueError(
                f"data and err are expected to be in units of {self.sb_unit}"
            )

        # the conversion in done in-place to avoid making copies of the data;
        # use a dictionary to set the value to avoid on-the-fly validation
        self.model["data"] /= self.l2_conv_factor
        self.model["err"] /= self.l2_conv_factor
        self.convolved_data /= self.l2_conv_factor

    def convert_sb_to_flux_density(self):
        """
        Convert the data and error Quantity arrays from MJy/sr (surface
        brightness) to flux density units.

        The flux density unit is defined by self.flux_unit.
        """
        if self.model.data.unit != self.sb_unit or self.model.err.unit != self.sb_unit:
            raise ValueError(
                f"data and err are expected to be in units of {self.sb_unit}"
            )

        # the conversion in done in-place to avoid making copies of the data;
        # use a dictionary to set the value to avoid on-the-fly validation
        self.model["data"] *= self.pixel_area
        self.model["data"] <<= self.flux_unit
        self.model["err"] *= self.pixel_area
        self.model["err"] <<= self.flux_unit
        self.convolved_data *= self.pixel_area
        self.convolved_data <<= self.flux_unit

    def convert_flux_density_to_sb(self):
        """
        Convert the data and error Quantity arrays from flux density units to
        MJy/sr (surface brightness).

        This is the inverse operation of `convert_sb_to_flux_density`.
        """
        if (
            self.model.data.unit != self.flux_unit
            or self.model.err.unit != self.flux_unit
        ):
            raise ValueError(
                f"data and err are expected to be in units of {self.flux_unit}"
            )

        self.model["data"] /= self.pixel_area
        self.model["data"] <<= self.sb_unit
        self.model["err"] /= self.pixel_area
        self.model["err"] <<= self.sb_unit
        self.convolved_data /= self.pixel_area
        self.convolved_data <<= self.sb_unit

    def convert_flux_to_abmag(self, flux, flux_err):
        """
        Convert flux (and error) to AB magnitude (and error).

        Parameters
        ----------
        flux, flux_err : `~astropy.unit.Quantity`
            The input flux and error arrays.

        Returns
        -------
        abmag, abmag_err : `~astropy.ndarray`
            The output AB magnitude and error arrays.
        """
        # ignore RunTimeWarning if flux or flux_err contains NaNs
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)

            # exact AB mag zero point
            flux_zpt = 10 ** (-0.4 * 48.60) * u.erg / u.s / u.cm**2 / u.Hz
            flux_zpt <<= self.flux_unit

            abmag_zpt = 2.5 * np.log10(flux_zpt.value)
            abmag = -2.5 * np.log10(flux.value) + abmag_zpt
            abmag_err = 2.5 * np.log10(1.0 + (flux_err.value / flux.value))

            # handle negative fluxes
            idx = flux.value < 0
            abmag[idx] = np.nan
            abmag_err[idx] = np.nan

        return abmag, abmag_err

    @lazyproperty
    def segment_colnames(self):
        """
        A dictionary of the output table column names and descriptions
        for the segment catalog.
        """
        desc = {}
        desc["label"] = "Unique source identification label number"
        desc["xcentroid"] = "X pixel value of the source centroid (0 indexed)"
        desc["ycentroid"] = "Y pixel value of the source centroid (0 indexed)"
        desc["sky_centroid"] = " Sky coordinate (ICRS) of the source centroid"
        desc["isophotal_flux"] = "Isophotal flux"
        desc["isophotal_flux_err"] = "Isophotal flux error"
        # isophotal_flux and isophotal_flux_err must be listed before isophotal_abmag
        # TEMP: do not include ABmags
        # desc["isophotal_abmag"] = "Isophotal AB magnitude"
        # desc["isophotal_abmag_err"] = "Isophotal AB magnitude error"
        desc["isophotal_area"] = "Isophotal area"
        desc["semimajor_sigma"] = (
            "1-sigma standard deviation along the semimajor axis of the 2D Gaussian function that has the same second-order central moments as the source"
        )
        desc["semiminor_sigma"] = (
            "1-sigma standard deviation along the semiminor axis of the 2D Gaussian function that has the same second-order central moments as the source"
        )
        desc["ellipticity"] = (
            "1 minus the ratio of the 1-sigma lengths of the semimajor and semiminor axes"
        )
        desc["orientation"] = (
            "The angle (degrees) between the positive X axis and the major axis (increases counter-clockwise)"
        )
        # orientation must be listed before sky_orientation
        desc["sky_orientation"] = (
            "The position angle (degrees) from North of the major axis"
        )
        desc["sky_bbox_ll"] = (
            "Sky coordinate (ICRS) of the lower-left vertex of the minimal bounding box of the source"
        )
        desc["sky_bbox_ul"] = (
            "Sky coordinate (ICRS) of the upper-left vertex of the minimal bounding box of the source"
        )
        desc["sky_bbox_lr"] = (
            "Sky coordinate (ICRS) of the lower-right vertex of the minimal bounding box of the source"
        )
        desc["sky_bbox_ur"] = (
            "Sky coordinate (ICRS) of the upper-right vertex of the minimal bounding box of the source"
        )

        self.column_desc.update(desc)

        return list(desc.keys())

    def set_segment_properties(self):
        """
        Calculate the segment-based source photometry and morphologies.

        The results are set as dynamic attributes on the class instance.
        """
        segm_cat = SourceCatalog(
            self.model.data,
            self.segment_img,
            convolved_data=self.convolved_data,
            error=self.model.err,
            wcs=self.wcs,
        )

        self.meta.update(segm_cat.meta)

        # rename some columns in the output catalog
        prop_names = {}
        prop_names["isophotal_flux"] = "segment_flux"
        prop_names["isophotal_flux_err"] = "segment_fluxerr"
        prop_names["isophotal_area"] = "area"

        # dynamically set the attributes
        for column in self.segment_colnames:
            # use the renamed column name if it exists in prop_names
            prop_name = prop_names.get(column, column)
            try:
                value = getattr(segm_cat, prop_name)
            except AttributeError:
                # isophotal_abmag, isophotal_abmag_err, sky_orientation
                value = getattr(self, prop_name)
            setattr(self, column, value)

    @lazyproperty
    def _xypos(self):
        """
        The (x, y) source positions, defined from the segmentation
        image.
        """
        return np.transpose((self.xcentroid, self.ycentroid))

    @lazyproperty
    def _xypos_aper(self):
        """
        The (x, y) source positions for the circular apertures and
        annuli.

        Non-finite positions are set to a large negative value. At this
        position the aperture will not overlap the data, thus returning
        NaN fluxes and errors.
        """
        xypos = self._xypos.copy()
        nanmask = ~np.isfinite(xypos)
        xypos[nanmask] = -1000.0
        return xypos

    @lazyproperty
    def _xypos_nonfinite_mask(self):
        """
        A 1D boolean mask where `True` values denote sources where
        either the xcentroid or the ycentroid is not finite.
        """
        return ~np.isfinite(self._xypos).all(axis=1)

    @lazyproperty
    def _isophotal_abmag(self):
        """
        The isophotal AB magnitude and error.
        """
        return self.convert_flux_to_abmag(self.isophotal_flux, self.isophotal_flux_err)

    @lazyproperty
    def isophotal_abmag(self):
        """
        The isophotal AB magnitude.
        """
        return self._isophotal_abmag[0]

    @lazyproperty
    def isophotal_abmag_err(self):
        """
        The isophotal AB magnitude error.
        """
        return self._isophotal_abmag[1]

    @lazyproperty
    def sky_orientation(self):
        """
        The orientation of the source major axis as the position angle
        in degrees measured East of North.
        """
        return (180.0 * u.deg) - self._wcs_angle + self.orientation

    def _make_aperture_colnames(self, name):
        """
        Make the aperture column names.

        There are separate columns for the flux/magnitude for each of
        the encircled energies and the total flux/magnitude.

        Parameters
        ----------
        name : {'flux', 'abmag'}
            The name type of the column.

        Returns
        -------
        colnames : list of str
            A list of the output column names.
        """
        colnames = []
        for aper_ee in self.aperture_ee:
            basename = f"aper{aper_ee}_{name}"
            colnames.append(basename)
            colnames.append(f"{basename}_err")
        colnames.extend([f"aper_total_{name}", f"aper_total_{name}_err"])

        return colnames

    def _make_aperture_descriptions(self, name):
        """
        Make aperture column descriptions.

        Parameters
        ----------
        name : {'flux', 'abmag'}
            The name type of the column.

        Returns
        -------
        descriptions : list of str
            A list of the output column descriptions.
        """
        if name == "flux":
            ftype = "Flux"
            ftype2 = "flux"
        elif name == "abmag":
            ftype = ftype2 = "AB magnitude"

        desc = []
        for aper_ee in self.aperture_ee:
            desc.append(
                f"{ftype} within the {aper_ee}% encircled energy circular aperture"
            )
            desc.append(
                f"{ftype} error within the {aper_ee}% encircled energy circular aperture"
            )

        desc.append(
            f"Total aperture-corrected {ftype2} based on the {self.aperture_ee[-1]}% encircled energy circular aperture; should be used only for unresolved sources."
        )
        desc.append(
            f"Total aperture-corrected {ftype2} error based on the {self.aperture_ee[-1]}% encircled energy circular aperture; should be used only for unresolved sources."
        )

        return desc

    @lazyproperty
    def aperture_flux_colnames(self):
        """
        The aperture flux column names.
        """
        return self._make_aperture_colnames("flux")

    @lazyproperty
    def aperture_flux_descriptions(self):
        """
        The aperture flux column descriptions.
        """
        return self._make_aperture_descriptions("flux")

    @lazyproperty
    def aperture_abmag_colnames(self):
        """
        The aperture AB magnitude column names.
        """
        return self._make_aperture_colnames("abmag")

    @lazyproperty
    def aperture_abmag_descriptions(self):
        """
        The aperture AB magnitude column descriptions.
        """
        return self._make_aperture_descriptions("abmag")

    @lazyproperty
    def aperture_colnames(self):
        """
        A dictionary of the output table column names and descriptions
        for the aperture catalog.
        """
        desc = {}
        desc["aper_bkg_flux"] = (
            "The local background value calculated as the sigma-clipped median value in the background annulus aperture"
        )
        desc["aper_bkg_flux_err"] = (
            "The standard error of the sigma-clipped median background value"
        )

        for idx, colname in enumerate(self.aperture_flux_colnames):
            desc[colname] = self.aperture_flux_descriptions[idx]
        # TEMP: do not include ABmags
        # for idx, colname in enumerate(self.aperture_abmag_colnames):
        #     desc[colname] = self.aperture_abmag_descriptions[idx]

        self.column_desc.update(desc)

        return list(desc.keys())

    @lazyproperty
    def _aper_local_background(self):
        """
        Estimate the local background and error using a circular annulus
        aperture.

        The local background is the sigma-clipped median value in the
        annulus.  The background error is the standard error of the
        median, sqrt(pi / 2N) * std.
        """
        bkg_aper = CircularAnnulus(
            self._xypos_aper,
            self.aperture_params["bkg_aperture_inner_radius"],
            self.aperture_params["bkg_aperture_outer_radius"],
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
                bkg_data = mask.get_values(self.model.data.value)
                values = sigclip(bkg_data, masked=False)
                nvalues.append(values.size)
                bkg_median.append(np.median(values))
                bkg_std.append(np.std(values))

            nvalues = np.array(nvalues)
            bkg_median = np.array(bkg_median)
            # standard error of the median
            bkg_median_err = np.sqrt(np.pi / (2.0 * nvalues)) * np.array(bkg_std)

        bkg_median <<= self.model.data.unit
        bkg_median_err <<= self.model.data.unit

        return bkg_median, bkg_median_err

    @lazyproperty
    def aper_bkg_flux(self):
        """
        The aperture local background flux (per pixel).
        """
        return self._aper_local_background[0]

    @lazyproperty
    def aper_bkg_flux_err(self):
        """
        The aperture local background flux error (per pixel).
        """
        return self._aper_local_background[1]

    def set_aperture_properties(self):
        """
        Calculate the aperture photometry.

        The results are set as dynamic attributes on the class instance.
        """
        apertures = [
            CircularAperture(self._xypos_aper, radius)
            for radius in self.aperture_params["aperture_radii"]
        ]
        aper_phot = aperture_photometry(
            self.model.data, apertures, error=self.model.err
        )

        for i, aperture in enumerate(apertures):
            flux_col = f"aperture_sum_{i}"
            flux_err_col = f"aperture_sum_err_{i}"

            # subtract the local background measured in the annulus
            aper_phot[flux_col] -= self.aper_bkg_flux * aperture.area

            flux = aper_phot[flux_col]
            flux_err = aper_phot[flux_err_col]
            abmag, abmag_err = self.convert_flux_to_abmag(flux, flux_err)

            idx0 = 2 * i
            idx1 = (2 * i) + 1
            setattr(self, self.aperture_flux_colnames[idx0], flux)
            setattr(self, self.aperture_flux_colnames[idx1], flux_err)
            setattr(self, self.aperture_abmag_colnames[idx0], abmag)
            setattr(self, self.aperture_abmag_colnames[idx1], abmag_err)

    @lazyproperty
    def extras_colnames(self):
        """
        A dictionary of the output table column names and descriptions
        for the additional catalog values.
        """
        desc = {}
        for idx, colname in enumerate(self.ci_colnames):
            desc[colname] = self.ci_colname_descriptions[idx]

        desc["flags"] = "Data quality flags"
        desc["is_extended"] = "Flag indicating whether the source is extended"
        desc["sharpness"] = "The DAOFind source sharpness statistic"
        desc["roundness"] = "The DAOFind source roundness statistic"
        desc["nn_label"] = "The label number of the nearest neighbor"
        desc["nn_dist"] = "The distance in pixels to the nearest neighbor"
        self.column_desc.update(desc)

        return list(desc.keys())

    @lazyproperty
    def flags(self):
        """
        Data quality flags.
        """
        xyidx = np.round(self._xypos).astype(int)

        try:
            # L2 images have a dq array
            dqflags = self.model.dq[xyidx[:, 1], xyidx[:, 0]]
            # if dqflags contains the DO_NOT_USE flag, set to DO_NOT_USE
            # (dq=1), otherwise 0
            flags = dqflags & pixel.DO_NOT_USE

        except AttributeError:
            # L3 images
            mask = self.model.weight == 0
            flags = mask[xyidx[:, 1], xyidx[:, 0]].astype(int)

        return flags

    @lazyproperty
    def _ci_ee_indices(self):
        """
        The EE indices for the concentration indices.
        """
        # NOTE: the EE values are always in increasing order
        return (0, 1), (1, 2), (0, 2)

    @lazyproperty
    def ci_colnames(self):
        """
        The column names of the three concentration indices.
        """
        return [
            f"CI_{self.aperture_ee[j]}_{self.aperture_ee[i]}"
            for (i, j) in self._ci_ee_indices
        ]

    @lazyproperty
    def ci_colname_descriptions(self):
        """
        The concentration indices column descriptions.
        """
        return [
            "Concentration index calculated as "
            f"({self.aperture_flux_colnames[2 * j]} / "
            f"{self.aperture_flux_colnames[2 * i]})"
            for (i, j) in self._ci_ee_indices
        ]

    @lazyproperty
    def concentration_indices(self):
        """
        A list of concentration indices, calculated as the flux
        ratios of:

            * the (middle / smallest) aperture flux ratio
              e.g., CI_50_30 = aper50_flux / aper30_flux
            * the (largest / middle) aperture flux ratio
              e.g., CI_70_50 = aper70_flux / aper50_flux
            * the (largest / smallest) aperture flux ratio
              e.g., CI_70_30 = aper70_flux / aper30_flux
        """
        fluxes = [
            (self.aperture_flux_colnames[2 * j], self.aperture_flux_colnames[2 * i])
            for (i, j) in self._ci_ee_indices
        ]
        return [
            getattr(self, flux1).value / getattr(self, flux2).value
            for flux1, flux2 in fluxes
        ]

    def set_ci_properties(self):
        """
        Set the concentration indices as dynamic attributes on the class
        instance.
        """
        for name, value in zip(self.ci_colnames, self.concentration_indices):
            setattr(self, name, value)

    @lazyproperty
    def is_extended(self):
        """
        Boolean indicating whether the source is extended.
        """
        mask1 = self.concentration_indices[0] > self.ci_star_thresholds[0]
        mask2 = self.concentration_indices[1] > self.ci_star_thresholds[1]
        return np.logical_and(mask1, mask2)

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
        for xcen, ycen in zip(*np.transpose(self._xypos_aper)):
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
        for xcen, ycen in zip(*np.transpose(self._xypos_aper)):
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
        if self.n_sources == 1:
            return [np.nan], [np.nan]

        # non-finite xypos causes memory errors on linux, but not MacOS
        tree = KDTree(self._xypos_aper)
        qdist, qidx = tree.query(self._xypos_aper, k=[2])
        return np.transpose(qdist)[0], np.transpose(qidx)[0]

    @lazyproperty
    def nn_label(self):
        """
        The label number of the nearest neighbor.

        A label value of -1 is returned if there is only one detected
        source and for sources with a non-finite xcentroid or ycentroid.
        """
        if self.n_sources == 1:
            return -1

        nn_label = self.label[self._kdtree_query[1]]
        # assign a label of -1 for non-finite xypos
        nn_label[self._xypos_nonfinite_mask] = -1

        return nn_label

    @lazyproperty
    def nn_dist(self):
        """
        The distance in pixels to the nearest neighbor.
        """
        nn_dist = self._kdtree_query[0]
        if self.n_sources == 1:
            # NaN if only one detected source
            return nn_dist * u.pixel

        # assign a distance of np.nan for non-finite xypos
        nn_dist[self._xypos_nonfinite_mask] = np.nan
        return nn_dist * u.pixel

    @lazyproperty
    def aper_total_flux(self):
        """
        The aperture-corrected total flux for sources, based on the flux
        in largest aperture.

        The aperture-corrected total flux should be used only for
        unresolved sources.
        """
        idx = self.n_aper - 1  # use apcorr for the largest EE (largest radius)
        flux = self.aperture_params["aperture_corrections"][idx] * getattr(
            self, self.aperture_flux_colnames[idx * 2]
        )
        return flux

    @lazyproperty
    def aper_total_flux_err(self):
        """
        The aperture-corrected total flux error for sources,
        based on the flux in largest aperture.

        The aperture-corrected total flux error should be used only for
        unresolved sources.
        """
        idx = self.n_aper - 1  # use apcorr for the largest EE (largest radius)
        flux_err = self.aperture_params["aperture_corrections"][idx] * getattr(
            self, self.aperture_flux_colnames[idx * 2 + 1]
        )
        return flux_err

    @lazyproperty
    def _abmag_total(self):
        """
        The total AB magnitude and error.
        """
        return self.convert_flux_to_abmag(
            self.aper_total_flux, self.aper_total_flux_err
        )

    @lazyproperty
    def aper_total_abmag(self):
        """
        The aperture-corrected total AB magnitude.

        The aperture-corrected total magnitude should be used only for
        unresolved sources.
        """
        return self._abmag_total[0]

    @lazyproperty
    def aper_total_abmag_err(self):
        """
        The aperture-corrected total AB magnitude error.

        The aperture-corrected total magnitude error should be used only
        for unresolved sources.
        """
        return self._abmag_total[1]

    @lazyproperty
    def colnames(self):
        """
        The column name order for the final source catalog.
        """
        colnames = self.segment_colnames[0:4]
        colnames.extend(self.aperture_colnames)
        colnames.extend(self.extras_colnames)
        colnames.extend(self.segment_colnames[4:])
        return colnames

    def _update_metadata(self):
        """
        Update the metadata dictionary with the package version
        information and aperture parameters.
        """
        ver_key = "versions"
        if "version" in self.meta.keys():
            # depends on photutils version
            self.meta[ver_key] = self.meta.pop("version")

        ver_dict = self.meta.get("versions", None)
        if ver_dict is not None:
            ver_dict["romancal"] = romancal_version
            packages = [
                "Python",
                "numpy",
                "scipy",
                "astropy",
                "photutils",
                "gwcs",
                "romancal",
            ]
            ver_dict = {key: ver_dict[key] for key in packages if key in ver_dict}
            self.meta[ver_key] = ver_dict

        self.meta["aperture_params"] = self.aperture_params

    def _split_skycoord(self, table):
        """
        Split SkyCoord columns into separate RA and Dec columns.

        Parameters
        ----------
        table : `~astropy.table.Table`
            The input table.

        Returns
        -------
        table : `~astropy.table.Table`
            The output table with separate RA and Dec columns.
        """
        for colname in table.colnames:
            if isinstance(table[colname], SkyCoord):
                idx = table.colnames.index(colname)
                ra = table[colname].ra
                dec = table[colname].dec
                ra_colname = colname.replace("sky", "ra")
                dec_colname = colname.replace("sky", "dec")
                table.remove_column(colname)
                table.add_columns(
                    [ra, dec], names=[ra_colname, dec_colname], indexes=[idx, idx]
                )
                desc = self.column_desc[colname]
                ra_desc = desc.replace("Sky coordinate", "Right ascension")
                dec_desc = desc.replace("Sky coordinate", "Declination")
                table[ra_colname].info.description = ra_desc
                table[dec_colname].info.description = dec_desc

        return table

    @lazyproperty
    def catalog(self):
        """
        The final source catalog.
        """
        # convert L2 data units back to MJy/sr
        if isinstance(self.model, ImageModel):
            self.convert_l2_to_sb()

        self.convert_sb_to_flux_density()
        self.set_segment_properties()
        self.set_aperture_properties()
        self.set_ci_properties()

        catalog = QTable()
        for column in self.colnames:
            catalog[column] = getattr(self, column)
            catalog[column].info.description = self.column_desc[column]
        self._update_metadata()
        catalog.meta.update(self.meta)

        # convert QTable to Table to change Quantity columns to regular
        # columns with units
        catalog = Table(catalog)

        # split SkyCoord columns into separate RA and Dec columns
        catalog = self._split_skycoord(catalog)

        # restore units on input model back to MJy/sr
        self.convert_flux_density_to_sb()

        # restore L2 data units back to DN/s
        if isinstance(self.model, ImageModel):
            self.convert_sb_to_l2()

        return catalog
