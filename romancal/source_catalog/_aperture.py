"""
Module to calculate aperture photometry.
"""

import logging
import warnings

import numpy as np
from astropy import units as u
from astropy.stats import SigmaClip
from astropy.utils.decorators import lazyproperty
from astropy.utils.exceptions import AstropyUserWarning
from photutils.aperture import CircularAnnulus, CircularAperture, aperture_photometry

from romancal.source_catalog._wcs_utils import pixel_area_at

log = logging.getLogger(__name__)


class ApertureCatalog:
    """
    Class to calculate circular-aperture photometry and a local
    background estimate for each source.

    The catalog produces one ``aperXX_flux`` / ``aperXX_flux_err``
    column pair per radius in `CIRCLE_APERTURE_RADII_ARCSEC`, plus the
    ``aper_bkg_flux`` / ``aper_bkg_flux_err`` columns from a circular
    annulus around each source.

    Parameters
    ----------
    model : `~roman_datamodels.datamodels.ImageModel` or \
            `~roman_datamodels.datamodels.MosaicModel`
        The data model providing ``model.data`` and ``model.err`` (both
        carrying surface-brightness units).

    xypos_finite : 2D `~numpy.ndarray`
        Source positions with shape ``(N, 2)`` of ``(x, y)`` pixel
        coordinates. Non-finite positions must already have been
        replaced with finite placeholders by the caller.

    pixel_area_map : `~astropy.units.Quantity`
        2D array of per-pixel solid angles matching the shape of
        ``model.data``. Used to convert the annulus background back to
        a surface brightness.

    ee_spline : callable, optional
        Encircled-energy spline mapping aperture radius (pixels) to
        encircled-energy fraction. Required for ``is_extended`` and
        the per-aperture ``ee_fraction_*`` attributes; if `None`, the
        extendedness criterion falls back to all-False with a warning.

    requested_properties : iterable of str, optional
        If given, restrict ``self.properties`` to this subset of
        `available_properties`.
    """

    # Define the circular aperture radii in arcsec. The output flux
    # column names are deterministic from these (e.g. 0.2 arcsec ->
    # ``aper02_flux``) and so do not depend on the pixel scale.
    CIRCLE_APERTURE_RADII_ARCSEC = (0.1, 0.2, 0.4, 0.8, 1.6)
    ANNULUS_RADII_ARCSEC = (2.4, 2.8)

    # photutils apertures take a single scalar radius, so sources are
    # grouped into bins of similar pixel scale and measured in batches.
    # These bound the resulting fractional error in the aperture radius.
    RADIUS_BIN_TOLERANCE = 1e-4
    MAX_RADIUS_BINS = 64

    @classmethod
    def aperture_flux_colnames_for_radii(cls, radii_arcsec=None):
        """
        Compute the aperture flux column names from a set of radii in
        arcsec.

        This is a class method so callers can determine the column names
        without instantiating the catalog.
        """
        if radii_arcsec is None:
            radii_arcsec = cls.CIRCLE_APERTURE_RADII_ARCSEC
        return [f"aper{round(r / 0.1):02d}_flux" for r in radii_arcsec]

    def __init__(
        self,
        model,
        xypos_finite,
        pixel_area_map,
        *,
        ee_spline=None,
        requested_properties=None,
    ):
        self.model = model
        self.xypos_finite = xypos_finite
        self.pixel_area_map = pixel_area_map
        self.ee_spline = ee_spline

        self.fractions = []
        self._requested_properties = (
            None if requested_properties is None else set(requested_properties)
        )

        # `available_properties` and `properties` depend on the
        # dynamically-generated aperture flux column names, which are
        # derived from `aperture_radii` (a lazyproperty).
        self.calc_aperture_photometry()

    @lazyproperty
    def available_properties(self):
        """
        The full set of source-property column names this catalog can
        produce.

        These names are dynamic because they depend on the aperture
        radii.
        """
        names = []
        for colname in self.aperture_flux_colnames:
            names.append(colname)
            names.append(f"{colname}_err")
        names.extend(["aper_bkg_flux", "aper_bkg_flux_err"])
        return tuple(names)

    @lazyproperty
    def properties(self):
        """
        The source-property column names this catalog will produce,
        filtered by the caller's requested properties.
        """
        if self._requested_properties is None:
            return list(self.available_properties)
        return [
            prop
            for prop in self.available_properties
            if prop in self._requested_properties
        ]

    @lazyproperty
    def _source_pixel_scale(self):
        """
        The angular size (in arcsec) of a pixel at each source position,
        defined as the square root of the local pixel solid angle.
        """
        return np.sqrt(self._source_pixel_area)

    @lazyproperty
    def aperture_radii(self):
        """
        A dictionary of the aperture radii used for aperture photometry.

        The dictionary has four keys:

        * ``'circle'`` and ``'annulus'`` contain radii as
          `~astropy.units.Quantity` arrays in arcsec
        * ``'circle_pix'`` and ``'annulus_pix'`` contain the same radii
          as plain floating-point arrays in pixels, with shape
          ``(n_radii, n_sources)``.

        The pixel radii vary from source to source because the pixel
        solid angle varies across the detector; a single pixel radius
        would subtend a different sky area in different parts of the
        image.

        These are the exact radii the apertures represent. The
        photometry itself measures sources in batches that share a
        common radius, which reproduces these values to within
        `RADIUS_BIN_TOLERANCE`; see `_radius_bins`.
        """
        params = {}
        radii = np.array(self.CIRCLE_APERTURE_RADII_ARCSEC) << u.arcsec
        annulus_radii = np.array(self.ANNULUS_RADII_ARCSEC) << u.arcsec
        params["circle"] = radii.copy()
        params["annulus"] = annulus_radii.copy()

        scale = self._source_pixel_scale[np.newaxis, :]
        params["circle_pix"] = (radii[:, np.newaxis] / scale).to_value(
            u.dimensionless_unscaled
        )
        params["annulus_pix"] = (annulus_radii[:, np.newaxis] / scale).to_value(
            u.dimensionless_unscaled
        )

        return params

    @lazyproperty
    def _source_pixel_area(self):
        """
        The solid angle (in arcsec**2) of the pixel containing each
        source.

        ``model.data`` carries the per-pixel solid angle, so recovering
        a surface brightness from it requires dividing by the local
        area rather than by a single image-wide value.
        """
        return pixel_area_at(
            self.pixel_area_map, self.xypos_finite[:, 0], self.xypos_finite[:, 1]
        ).to(u.arcsec**2)

    @lazyproperty
    def _radius_bins(self):
        """
        Group sources into bins of similar local pixel scale.

        `~photutils.aperture.CircularAperture` takes a single scalar
        radius, so sources are measured in batches that share a radius.
        Bins are spaced uniformly between the smallest and largest
        pixel scale, and each source is assigned to the nearest bin.
        The number of bins is chosen so that the aperture radius is
        accurate to roughly `RADIUS_BIN_TOLERANCE`; an image with
        uniform pixels needs only one bin.

        Returns
        -------
        bins : list of (`~numpy.ndarray`, float)
            Source indices and the pixel scale (in arcsec) to use for
            that group. A bin that no source falls in is returned empty
            rather than skipped; it simply measures nothing.
        """
        scale = self._source_pixel_scale.to_value(u.arcsec)
        if scale.size == 0:
            return []

        spread = np.ptp(scale) / scale.mean()
        n_bins = int(np.ceil(spread / self.RADIUS_BIN_TOLERANCE))
        n_bins = int(np.clip(n_bins, 1, min(self.MAX_RADIUS_BINS, scale.size)))
        if n_bins == 1:
            return [(np.arange(scale.size), scale.mean())]

        centers = np.linspace(scale.min(), scale.max(), n_bins)
        spacing = centers[1] - centers[0]
        nearest = np.round((scale - centers[0]) / spacing).astype(int)

        bins = []
        for value in range(n_bins):
            bins.append((np.flatnonzero(nearest == value), centers[value]))
        return bins

    @lazyproperty
    def _aperture_background(self):
        """
        The local background and error estimated using a circular
        annulus aperture.

        The local background is the sigma-clipped median value in the
        annulus. The background error is the standard error of the
        median, sqrt(pi / (2 * N)) * std.
        """
        r_in_arcsec, r_out_arcsec = self.ANNULUS_RADII_ARCSEC
        bkg_aper_masks = [None] * self.xypos_finite.shape[0]
        for idx, scale in self._radius_bins:
            annulus = CircularAnnulus(
                self.xypos_finite[idx], r_in_arcsec / scale, r_out_arcsec / scale
            )
            for source, mask in zip(idx, annulus.to_mask(method="center"), strict=True):
                bkg_aper_masks[source] = mask
        sigclip = SigmaClip(sigma=3.0)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            warnings.simplefilter("ignore", category=AstropyUserWarning)

            unit = self.model.data.unit
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
                    # Handle case where source is completely masked due to
                    # forced photometry
                    med <<= unit
                    std <<= unit
                bkg_median.append(med)
                bkg_std.append(std)

            nvalues = np.array(nvalues)
            pixel_area = self._source_pixel_area
            bkg_median = u.Quantity(bkg_median) / pixel_area
            bkg_std = u.Quantity(bkg_std) / pixel_area

            # Standard error of the median
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
    def aperture_radius_name(self):
        """
        Two-digit aperture radius labels in tenths of arcsec.

        For example, the aperture with radius ``0.2 arcsec`` produces
        the label ``"02"`` and contributes the column ``aper02_flux``.
        """
        return [
            f"{np.round((radius / 0.1 * u.arcsec).value).astype(int):02d}"
            for radius in self.aperture_radii["circle"]
        ]

    @lazyproperty
    def aperture_flux_colnames(self):
        """
        The aperture flux column names.

        The column names are based on the circular aperture radii in
        tenths of arcsec. For example, the flux in a r=0.2 arcsec
        aperture is "aper02_flux".
        """
        return [f"aper{name}_flux" for name in self.aperture_radius_name]

    def calc_ee_fractions(self):
        """
        Calculate the encircled energy fractions for each aperture
        radius.
        """
        if self.ee_spline is not None:
            for radius, name in zip(
                self.aperture_radii["circle"],
                self.aperture_radius_name,
                strict=True,
            ):
                attr_name = f"ee_fraction_{name}"
                self.fractions.append(attr_name)
                setattr(self, attr_name, self.ee_spline(radius))

    @lazyproperty
    def is_extended(self):
        """
        Boolean array indicating if the source is extended.

        If no encircled-energy spline is available, all sources are
        flagged as not extended (with a warning), since the
        extendedness criterion depends on the EE correction model.
        """
        if self.ee_spline is None:
            log.warning(
                "Cannot determine if source is extended without an "
                "encircled-energy spline (ee_spline). Returning all "
                "False."
            )
            return np.zeros(self.xypos_finite.shape[0], dtype=bool)

        ee_ratio = self.ee_fraction_04 / self.ee_fraction_02
        return self.aper04_flux > (self.aper02_flux * 1.2 * ee_ratio)

    def calc_aperture_photometry(self, subtract_local_bkg=False):
        """
        Calculate circular-aperture photometry for all sources.

        For each radius in `aperture_radii['circle_pix']` the flux
        and flux error are stored on this instance as attributes
        named ``aperXX_flux`` and ``aperXX_flux_err`` (see
        `aperture_radius_name` for the ``XX`` formatting). When an
        encircled-energy spline is configured, ``ee_fraction_XX``
        attributes are also populated by `calc_ee_fractions`.

        Parameters
        ----------
        subtract_local_bkg : bool, optional
            If `True`, subtract the local annulus background
            (`aper_bkg_flux`) scaled by the geometric overlap of each
            aperture with the data array before assigning the flux
            attributes.
        """
        radii_arcsec = np.array(self.CIRCLE_APERTURE_RADII_ARCSEC)
        n_radii = radii_arcsec.size
        n_sources = self.xypos_finite.shape[0]
        unit = getattr(self.model.data, "unit", None)

        flux = np.zeros((n_radii, n_sources))
        flux_err = np.zeros((n_radii, n_sources))

        # Sources are measured in batches sharing a common radius; see
        # `_radius_bins` for how the batches are chosen.
        for idx, scale in self._radius_bins:
            apertures = [
                CircularAperture(self.xypos_finite[idx], radius / scale)
                for radius in radii_arcsec
            ]
            phot = aperture_photometry(self.model.data, apertures, error=self.model.err)
            for i, aperture in enumerate(apertures):
                values = phot[f"aperture_sum_{i}"]
                if subtract_local_bkg:
                    # Subtract the local background measured in the annulus
                    areas = aperture.area_overlap(self.model.data)
                    values = values - self.aper_bkg_flux[idx] * areas
                flux[i, idx] = getattr(values, "value", values)
                errors = phot[f"aperture_sum_err_{i}"]
                flux_err[i, idx] = getattr(errors, "value", errors)

        for i in range(n_radii):
            flux_col = self.aperture_flux_colnames[i]
            values = flux[i].astype(np.float32)
            errors = flux_err[i].astype(np.float32)
            if unit is not None:
                values = values << unit
                errors = errors << unit
            setattr(self, flux_col, values)
            setattr(self, f"{flux_col}_err", errors)

        self.calc_ee_fractions()
