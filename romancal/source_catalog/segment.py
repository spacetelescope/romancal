"""
Module to calculate the segmentation-based properties.
"""

import inspect
import warnings

import astropy.units as u
import numpy as np
from astropy.utils import lazyproperty
from photutils.segmentation import SourceCatalog


class SegmentCatalog:
    """
    Class to calculate the segmentation-based properties.

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

    pixel_area : `~astropy.units.Quantity`
        The pixel area in steradians. This is used to convert various
        measuments from pixels to arcseconds.

    wcs_angle : `~astropy.units.Quantity`
        The angle (as a Quantity in degrees) measured counterclockwise
        from the positive x axis to the "North" axis of the celestial
        coordinate system.

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
        pixel_area,
        wcs_angle,
        detection_cat=None,
        flux_unit="nJy",
    ):
        self.model = model
        self.segment_img = segment_img
        self.convolved_data = convolved_data
        self.pixel_area = pixel_area
        self.wcs_angle = wcs_angle
        self.detection_cat = detection_cat
        self.flux_unit = flux_unit

        self.names = []
        self.wcs = self.model.meta.wcs
        self.meta = {}

        # needed for detection_cat
        self.source_cat = None

        # calculate the segment properties
        self.calc_segment_properties()

        # lazyproperties are not set until accessed so we need to
        # manually append them
        for name in self._lazyproperties:
            self.names.append(name)

        # add the placeholder attributes
        self.add_placeholders()

    @property
    def _lazyproperties(self):
        """
        A list of all class lazyproperties.
        """

        def islazyproperty(obj):
            return isinstance(obj, lazyproperty)

        return [
            i[0] for i in inspect.getmembers(self.__class__, predicate=islazyproperty)
        ]

    @staticmethod
    def convert_flux_to_abmag(flux, flux_err):
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

            abmag = flux.to(u.ABmag)
            abmag_err = 2.5 * np.log10(1.0 + (flux_err / flux))
            abmag_err = abmag_err.value * u.ABmag

            # handle negative fluxes
            idx = flux.value < 0
            abmag[idx] = np.nan
            abmag_err[idx] = np.nan

        return abmag, abmag_err

    @lazyproperty
    def pixel_scale(self):
        """
        The pixel scale in arcseconds.
        """
        # assumes square pixels
        return np.sqrt(self.pixel_area).to(u.arcsec)

    def calc_segment_properties(self):
        """
        Calculate the segment-based properties provided by
        `~photutils.segmentation.SourceCatalog`.

        The results are set as dynamic attributes on the class instance.
        """
        # detection_cat is a RomanSourceCatalog object
        # segment_cat is a SegmentCatalog object
        # source_cat is a SourceCatalog object
        if self.detection_cat is not None:
            detection_cat = self.detection_cat.segment_cat.source_cat
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
        self.source_cat = segm_cat
        self.meta.update(segm_cat.meta)

        # Extract the properties from the segment catalog. These
        # names are the SourceCatalog property names and the order
        # is not important.
        photutils_names = (
            "label",
            "xcentroid",
            "ycentroid",
            "xcentroid_win",
            "ycentroid_win",
            "sky_centroid",
            "sky_centroid_win",
            "bbox_xmin",
            "bbox_xmax",
            "bbox_ymin",
            "bbox_ymax",
            "area",
            "semimajor_sigma",
            "semiminor_sigma",
            "fwhm",
            "orientation",
            "ellipticity",
            "cxx",
            "cxy",
            "cyy",
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
        name_map["xcentroid_win"] = "x_centroid_win"
        name_map["ycentroid_win"] = "y_centroid_win"
        name_map["area"] = "segment_area"
        name_map["semimajor_sigma"] = "semimajor"
        name_map["semiminor_sigma"] = "semiminor"
        name_map["orientation"] = "orientation_pix"
        name_map["segment_fluxerr"] = "segment_flux_err"
        name_map["kron_fluxerr"] = "kron_flux_err"

        # set the source properties as attributes of this instance
        for name in photutils_names:
            new_name = name_map.get(name, name)
            value = getattr(segm_cat, name)

            # handle any unit conversions
            if new_name in (
                "x_centroid",
                "y_centroid",
                "x_centroid_win",
                "y_centroid_win",
            ):
                value *= u.pix

            # pix**2 -> arcsec**2
            if new_name == "segment_area":
                value = value.value * self.pixel_area.to(u.arcsec**2)

            # pix -> arcsec
            if new_name in ("semimajor", "semiminor", "fwhm", "kron_radius"):
                value = value.value * self.pixel_scale

            # 1 / pix**2 -> 1 / arcsec**-2
            if new_name in ("cxx", "cxy", "cyy"):
                value = value.value / self.pixel_area.to(u.arcsec**2)

            # remove dimensionless units
            if new_name == "ellipticity":
                value = value.value

            # change the photutils dtypes
            if new_name not in ("sky_centroid", "sky_centroid_win"):
                if np.issubdtype(value.dtype, np.integer):
                    value = value.astype(np.int32)
                elif np.issubdtype(value.dtype, np.floating):
                    value = value.astype(np.float32)

            # split the sky_centroid values into separate RA and Dec
            # values
            if new_name == "sky_centroid":
                self.ra_centroid = value.ra
                self.dec_centroid = value.dec
                self.names.extend(["ra_centroid", "dec_centroid"])
            elif new_name == "sky_centroid_win":
                self.ra_centroid_win = value.ra
                self.dec_centroid_win = value.dec
                self.names.extend(["ra_centroid_win", "dec_centroid_win"])
            else:
                setattr(self, new_name, value)
                self.names.append(new_name)

    def add_placeholders(self):
        """
        Add placeholder attributes to the class instance.
        """
        # pixel columns
        pix_columns = (
            "x_centroid_err",
            "y_centroid_err",
            "x_centroid_win_err",
            "y_centroid_win_err",
        )
        pix_value = np.zeros(self.source_cat.nlabels, dtype=np.float32) * u.pix
        for name in pix_columns:
            if not hasattr(self, name):
                setattr(self, name, pix_value)

        # sky columns
        sky_columns = (
            "ra_centroid_err",
            "dec_centroid_err",
            "ra_centroid_win_err",
            "dec_centroid_win_err",
        )
        sky_value = np.zeros(self.source_cat.nlabels, dtype=np.float32) * u.arcsec
        for name in sky_columns:
            if not hasattr(self, name):
                setattr(self, name, sky_value)

        for name in [*pix_columns, *sky_columns]:
            self.names.append(name)

    @lazyproperty
    def orientation_sky(self):
        """
        The position angle of the source major axis in degrees measured
        East of North.
        """
        return ((180.0 * u.deg) - self.wcs_angle + self.orientation_pix).astype(
            np.float32
        )

    @lazyproperty
    def _kron_abmag(self):
        return self.convert_flux_to_abmag(self.kron_flux, self.kron_flux_err)

    @lazyproperty
    def kron_abmag(self):
        """
        The Kron magnitude in AB magnitudes.
        """
        return self._kron_abmag[0]

    @lazyproperty
    def kron_abmag_err(self):
        """
        The Kron magnitude error in AB magnitudes.
        """
        return self._kron_abmag[1]

    @lazyproperty
    def fluxfrac_radius_50(self):
        """
        The radius (in arcsec) at which the flux fraction is 50%.
        """
        value = self.source_cat.fluxfrac_radius(0.5)
        return (value.value * self.pixel_scale).astype(np.float32)
