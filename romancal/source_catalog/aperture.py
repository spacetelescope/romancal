"""
Module to calculate aperture photometry.
"""

import warnings

import numpy as np
from astropy import units as u
from astropy.stats import SigmaClip
from astropy.utils.decorators import lazyproperty
from astropy.utils.exceptions import AstropyUserWarning
from photutils.aperture import CircularAnnulus, CircularAperture, aperture_photometry


class ApertureCatalog:
    """
    Class to calculate aperture photometry.
    """

    def __init__(self, model, pixel_scale, xypos_finite):
        self.model = model
        self.pixel_scale = pixel_scale
        self.xypos_finite = xypos_finite

        self.names = []

        self.calc_aperture_photometry()
        self.names.extend(["aper_bkg_flux", "aper_bkg_flux_err"])

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

        radii /= self.pixel_scale
        annulus_radii /= self.pixel_scale
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
            self.xypos_finite,
            self.aperture_radii["annulus_pix"][0],
            self.aperture_radii["annulus_pix"][1],
        )
        bkg_aper_masks = bkg_aper.to_mask(method="center")
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
                    # handle case where source is completely masked due to
                    # forced photometry
                    med <<= unit
                    std <<= unit
                bkg_median.append(med)
                bkg_std.append(std)

            nvalues = np.array(nvalues)
            bkg_median = u.Quantity(bkg_median) / u.pix
            bkg_std = u.Quantity(bkg_std) / u.pix

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

    @lazyproperty
    def aperture_flux_descriptions(self):
        """
        The aperture flux descriptions.

        The descriptions are based on the circular aperture radii in
        tenths of arcsec. For example, the flux in a r=0.2 arcsec
        aperture is "aper02_flux".
        """
        descriptions = {}
        for i, colname in enumerate(self.aperture_flux_colnames):
            desc = (
                "Flux within a circular aperture of radius="
                f"{self.aperture_radii['circle'][i]:0.1f}"
            )
            descriptions[colname] = desc
        return descriptions

    def calc_aperture_photometry(self, subtract_local_bkg=False):
        """
        Calculate the aperture photometry.

        The results are set as dynamic attributes on the class instance.
        """
        apertures = [
            CircularAperture(self.xypos_finite, radius)
            for radius in self.aperture_radii["circle_pix"]
        ]
        aper_phot = aperture_photometry(
            self.model.data, apertures, error=self.model.err
        )

        for i, aperture in enumerate(apertures):
            tmp_flux_col = f"aperture_sum_{i}"
            tmp_flux_err_col = f"aperture_sum_err_{i}"

            if subtract_local_bkg:
                # subtract the local background measured in the annulus
                aper_areas = aperture.area_overlap(self.model.data)
                aper_phot[tmp_flux_col] -= self.aper_bkg_flux * aper_areas

            # set the flux and error attributes
            flux_col = self.aperture_flux_colnames[i]
            flux_err_col = f"{flux_col}_err"
            setattr(self, flux_col, aper_phot[tmp_flux_col].astype(np.float32))
            setattr(self, flux_err_col, aper_phot[tmp_flux_err_col].astype(np.float32))
            self.names.append(flux_col)
            self.names.append(flux_err_col)
