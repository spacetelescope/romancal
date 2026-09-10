"""
Module to calculate the segmentation-based properties.
"""

import inspect
import warnings
from typing import ClassVar

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

    mask : 2D bool `~numpy.ndarray` or `None`
        A 2D boolean array with the same shape as the input data, where
        `True` values indicate masked pixels. If `None`, then no pixels
        are masked.

    pixel_area : `~astropy.units.Quantity`
        The pixel area in steradians. This is used to convert various
        measurements from pixels to arcseconds.

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

    requested_properties : iterable of str or `None`, optional
        The output column names actually needed by the caller. If
        provided, only properties whose output names are in this
        iterable are computed and exposed via ``self.properties``. If
        `None` (default), all properties this class can produce are
        computed and exposed.

    Notes
    -----
    ``model.err`` is assumed to be the total error array corresponding
    to the input science ``model.data`` array. It is assumed to include
    *all* sources of error, including the Poisson error of the sources,
    and have the same shape and units as the science data array.
    """

    # Mapping of photutils SourceCatalog property names to the output
    # catalog name(s) they produce. A photutils property may produce
    # more than one output column (e.g. sky_centroid -> ra/dec).
    _photutils_to_outputs: ClassVar = {
        "label": ("label",),
        "x_centroid": ("x_centroid",),
        "y_centroid": ("y_centroid",),
        "x_centroid_win": ("x_centroid_win",),
        "y_centroid_win": ("y_centroid_win",),
        "sky_centroid": ("ra_centroid", "dec_centroid"),
        "sky_centroid_win": ("ra_centroid_win", "dec_centroid_win"),
        "bbox_xmin": ("bbox_xmin",),
        "bbox_xmax": ("bbox_xmax",),
        "bbox_ymin": ("bbox_ymin",),
        "bbox_ymax": ("bbox_ymax",),
        "area": ("segment_area",),
        "semimajor_axis": ("semimajor",),
        "semiminor_axis": ("semiminor",),
        "orientation": ("orientation_pix",),
        "ellipticity": ("ellipticity",),
        "ellipse_cxx": ("cxx",),
        "ellipse_cxy": ("cxy",),
        "ellipse_cyy": ("cyy",),
        "fwhm": ("fwhm",),
        "segment_flux": ("segment_flux",),
        "segment_flux_err": ("segment_flux_err",),
        "kron_radius": ("kron_radius",),
        "kron_flux": ("kron_flux",),
        "kron_flux_err": ("kron_flux_err",),
    }

    # Lazy property output names mapped to the upstream photutils
    # property names they depend on. Used to ensure that the upstream
    # properties get computed when only the derived lazy property is
    # requested.
    _lazy_dependencies: ClassVar = {
        "orientation_sky": ("orientation",),
        "kron_abmag": ("kron_flux", "kron_flux_err"),
        "kron_abmag_err": ("kron_flux", "kron_flux_err"),
    }

    # Placeholder columns (zero-valued) to be added until the proper
    # values are computed.
    _pix_placeholder_columns = (
        "x_centroid_err",
        "y_centroid_err",
        "x_centroid_win_err",
        "y_centroid_win_err",
    )
    _sky_placeholder_columns = (
        "ra_centroid_err",
        "dec_centroid_err",
        "ra_centroid_win_err",
        "dec_centroid_win_err",
    )

    def __init__(
        self,
        model,
        segment_img,
        convolved_data,
        mask,
        pixel_area,
        wcs_angle,
        detection_cat=None,
        *,
        requested_properties=None,
    ):
        self.model = model
        self.segment_img = segment_img
        self.convolved_data = convolved_data
        self.mask = mask
        self.pixel_area = pixel_area
        self.wcs_angle = wcs_angle
        self.detection_cat = detection_cat
        self._requested_properties = (
            None if requested_properties is None else set(requested_properties)
        )

        self.properties = []
        self.wcs = self.model.meta.wcs
        self.meta = {}

        self.source_cat = None  # Needed for detection_cat

        # Calculate the segment properties
        self.calc_segment_properties()

        # Lazy properties are not set until accessed so we need to
        # manually append them
        for name in self._lazyproperties:
            if not self._is_requested(name):
                continue
            self.properties.append(name)

        # Add the placeholder attributes
        self.add_placeholders()

    def _is_requested(self, name):
        """
        Whether the given output column name was requested.
        """
        return self._requested_properties is None or name in self._requested_properties

    @property
    def available_properties(self):
        """
        The full set of source-property column names this catalog can
        produce.
        """
        names = []
        for outputs in self._photutils_to_outputs.values():
            names.extend(outputs)
        names.extend(self._lazyproperties)
        names.extend(self._pix_placeholder_columns)
        names.extend(self._sky_placeholder_columns)
        return tuple(names)

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

        For non-positive fluxes (``flux <= 0``), both ``abmag`` and
        ``abmag_err`` are set to NaN.

        Parameters
        ----------
        flux, flux_err : `~numpy.ndarray`
            The input flux and error arrays.

        Returns
        -------
        abmag, abmag_err : `~astropy.ndarray`
            The output AB magnitude and error arrays.
        """
        # Mask non-positive fluxes
        invalid = flux.value <= 0

        # Ignore RunTimeWarnings from masked/non-finite values.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)

            abmag = flux.to(u.ABmag)
            abmag_err = 2.5 * np.log10(1.0 + (flux_err / flux))
            abmag_err = abmag_err.value * u.ABmag

            abmag[invalid] = np.nan
            abmag_err[invalid] = np.nan

        return abmag, abmag_err

    @lazyproperty
    def pixel_scale(self):
        """
        The pixel scale in arcseconds (assuming square pixels).
        """
        return np.sqrt(self.pixel_area).to(u.arcsec)

    def calc_segment_properties(self):
        """
        Calculate the segment-based properties provided by
        `~photutils.segmentation.SourceCatalog`.

        The results are set as dynamic attributes on the class instance.
        """
        # detection_cat is a romancal RomanSourceCatalog object
        # segment_cat is a romancal SegmentCatalog object
        # source_cat is a photutils SourceCatalog object
        if self.detection_cat is not None:
            detection_cat = self.detection_cat.segment_cat.source_cat
        else:
            detection_cat = None

        segm_cat = SourceCatalog(
            self.model.data,
            self.segment_img,
            convolved_data=self.convolved_data,
            error=self.model.err,
            mask=self.mask,
            wcs=self.wcs,
            aperture_mask_method="mask",
            detection_catalog=detection_cat,
        )

        self.source_cat = segm_cat
        self.meta.update(segm_cat.meta)

        # Determine which photutils properties to compute. A property
        # is computed if any of its output column names was requested
        # (or if no specific columns were requested, in which case
        # we compute all properties). ``label`` is always required
        # because it is needed for downstream table joins. We also
        # pull in upstream photutils dependencies for any requested
        # lazy property (e.g., ``orientation_sky`` depends on
        # ``orientation``).
        if self._requested_properties is None:
            needed = None
        else:
            needed = set(self._requested_properties)
            for lazy_name, deps in self._lazy_dependencies.items():
                if lazy_name in needed:
                    needed.update(deps)

        photutils_names = [
            pname
            for pname, outputs in self._photutils_to_outputs.items()
            if pname == "label"
            or needed is None
            or pname in needed
            or any(out in needed for out in outputs)
        ]

        # Map photutils names to the output catalog names
        name_map = {}
        name_map["area"] = "segment_area"
        name_map["semimajor_axis"] = "semimajor"
        name_map["semiminor_axis"] = "semiminor"
        name_map["orientation"] = "orientation_pix"
        name_map["ellipse_cxx"] = "cxx"
        name_map["ellipse_cxy"] = "cxy"
        name_map["ellipse_cyy"] = "cyy"

        # Set the source properties as attributes of this instance
        for name in photutils_names:
            new_name = name_map.get(name, name)
            value = getattr(segm_cat, name)

            # Handle any unit conversions needed for specific columns
            # (photutils -> romancal).

            # unitless -> pix
            if new_name in (
                "x_centroid",
                "y_centroid",
                "x_centroid_win",
                "y_centroid_win",
            ):
                if not isinstance(value, u.Quantity):
                    value *= u.pix

            # pix**2 -> arcsec**2
            if new_name == "segment_area":
                if isinstance(value, u.Quantity) and value.unit.is_equivalent(u.pix**2):
                    value = value.value * self.pixel_area.to(u.arcsec**2)

            # pix -> arcsec
            if new_name in ("semimajor", "semiminor", "fwhm", "kron_radius"):
                if isinstance(value, u.Quantity) and value.unit.is_equivalent(u.pix):
                    value = value.value * self.pixel_scale

            # 1 / pix**2 -> 1 / arcsec**-2
            if new_name in ("cxx", "cxy", "cyy"):
                if isinstance(value, u.Quantity) and value.unit.is_equivalent(
                    1 / u.pix**2
                ):
                    value = value.value / self.pixel_area.to(u.arcsec**2)

            # Remove dimensionless units
            if new_name == "ellipticity":
                if isinstance(value, u.Quantity) and value.unit.is_equivalent(
                    u.dimensionless_unscaled
                ):
                    value = value.value

            # Change dtypes for romancal catalog
            if new_name not in ("sky_centroid", "sky_centroid_win"):
                if np.issubdtype(value.dtype, np.integer):
                    value = value.astype(np.int32)
                elif np.issubdtype(value.dtype, np.floating):
                    value = value.astype(np.float32)

            # Split the sky_centroid values into separate RA and Dec
            # values
            if new_name == "sky_centroid":
                self.ra_centroid = value.ra
                self.dec_centroid = value.dec
                for sub in ("ra_centroid", "dec_centroid"):
                    if self._is_requested(sub):
                        self.properties.append(sub)
            elif new_name == "sky_centroid_win":
                self.ra_centroid_win = value.ra
                self.dec_centroid_win = value.dec
                for sub in ("ra_centroid_win", "dec_centroid_win"):
                    if self._is_requested(sub):
                        self.properties.append(sub)
            else:
                setattr(self, new_name, value)
                if self._is_requested(new_name):
                    self.properties.append(new_name)

    def add_placeholders(self):
        """
        Add placeholder (zero-valued) attributes for columns whose
        proper values are not yet computed.

        Each placeholder gets its own array so that downstream in-place
        updates to one column do not silently affect others.
        """
        n_labels = self.source_cat.n_labels

        placeholder_groups = (
            (self._pix_placeholder_columns, u.pix),
            (self._sky_placeholder_columns, u.arcsec),
        )

        for columns, unit in placeholder_groups:
            for name in columns:
                if not hasattr(self, name):
                    setattr(self, name, np.zeros(n_labels, dtype=np.float32) << unit)

                if self._is_requested(name) and name not in self.properties:
                    self.properties.append(name)

    @lazyproperty
    def orientation_sky(self):
        """
        The position angle of the source major axis in degrees measured
        East of North.
        """
        angle = ((180.0 * u.deg) - self.wcs_angle + self.orientation_pix) % (
            360.0 * u.deg
        )
        return angle.astype(np.float32)

    @lazyproperty
    def _kron_abmag(self):
        """
        The Kron magnitude and error in AB magnitudes.

        This avoids calling the flux-to-ABmag conversion twice: once for
        the Kron magnitude and once for its error.
        """
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
        The circular radius (in arcsec) at which the total Kron flux
        fraction is 50%.
        """
        value = self.source_cat.flux_radius(0.5)
        return (value.value * self.pixel_scale).astype(np.float32)
