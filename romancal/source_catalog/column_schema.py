"""
Column schema for the Roman source catalog.

The schema defines, for a given catalog type and configuration, the
ordered list of output columns and how the ``forced_`` prefix is applied
to a forced catalog.
"""

from astropy.utils import lazyproperty


class CatalogSchema:
    """
    Ordered column lists for the Roman source catalog.

    Parameters
    ----------
    cat_type : str
        The catalog type. One of ``"prompt"``, ``"dr_det"``,
        ``"dr_band"``, ``"psf_matched"``, ``"forced_full"``, or
        ``"forced_det"``.

    fit_psf : bool
        If `True`, include PSF-related columns in the schema.

    aperture_flux_colnames : sequence of str
        The base aperture flux column names (e.g. ``"aper02_flux"``).
        Each base name contributes a ``..._err`` column to the aperture
        column list.
    """

    def __init__(self, cat_type, fit_psf, aperture_flux_colnames):
        self.cat_type = cat_type
        self.fit_psf = bool(fit_psf)
        self._aperture_flux_colnames = list(aperture_flux_colnames)

    @lazyproperty
    def aper_colnames(self):
        """
        An ordered list of the aperture column names.
        """
        # Define the aperture background flux column names
        aper_colnames = [
            "aper_bkg_flux",
            "aper_bkg_flux_err",
        ]
        for colname in self._aperture_flux_colnames:
            aper_colnames.append(colname)
            aper_colnames.append(f"{colname}_err")
        return aper_colnames

    @lazyproperty
    def flux_colnames(self):
        """
        An ordered list of the flux column names.
        """
        other_colnames = [
            "segment_flux",
            "segment_flux_err",
            "kron_flux",
            "kron_flux_err",
            "kron_abmag",
            "kron_abmag_err",
        ]
        psf_colnames = ["psf_flux", "psf_flux_err"]
        # PSF-matched catalogs omit abmag columns (redundant, since they
        # can be derived from flux via -2.5*log10; we provide at least
        # one magnitude per source but not all variants) and PSF columns
        # (fit_psf is always False for matched images, but be explicit
        # for clarity).
        matched_other_colnames = [
            "segment_flux",
            "segment_flux_err",
            "kron_flux",
            "kron_flux_err",
        ]

        if self.cat_type in ("prompt", "forced_full", "dr_band"):
            flux_colnames = list(self.aper_colnames)
            if self.fit_psf:
                flux_colnames.extend(psf_colnames)
            flux_colnames.extend(other_colnames)

        elif self.cat_type == "psf_matched":
            flux_colnames = list(self.aper_colnames)
            flux_colnames.extend(matched_other_colnames)

        elif self.cat_type in ("dr_det", "forced_det"):
            flux_colnames = []

        else:
            raise ValueError(f"Unknown catalog type: {self.cat_type}")

        return flux_colnames

    @lazyproperty
    def band_colnames(self):
        """
        Band-specific columns for the multiband catalog.

        Also used for the forced catalog.
        """
        xypsf_colnames = ["x_psf", "x_psf_err", "y_psf", "y_psf_err"]
        skypsf_colnames = ["ra_psf", "dec_psf", "ra_psf_err", "dec_psf_err"]
        psf_flags_colnames = ["psf_flags", "psf_gof"]
        othershape_colnames = [
            "sharpness",
            "roundness1",
            "is_extended",
            "fluxfrac_radius_50",
        ]

        band_colnames = ["label"]  # needed to join the filter catalogs
        if self.fit_psf:
            band_colnames.extend(xypsf_colnames)
            band_colnames.extend(skypsf_colnames)
            band_colnames.extend(psf_flags_colnames)
        band_colnames.extend(othershape_colnames)
        band_colnames.extend(self.flux_colnames)
        return band_colnames

    @lazyproperty
    def column_names(self):
        """
        An ordered list of the output catalog column names.

        This list determines which values are calculated in the output
        catalog.
        """
        base_colnames = [
            "label",
            "flagged_spatial_id",
            "x_centroid",
            "y_centroid",
            "x_centroid_err",
            "y_centroid_err",
        ]
        xywin_colnames = [
            "x_centroid_win",
            "y_centroid_win",
            "x_centroid_win_err",
            "y_centroid_win_err",
        ]
        skybest_colnames = [
            "ra",
            "dec",
        ]
        sky_colnames = [
            "ra_centroid",
            "dec_centroid",
            "ra_centroid_err",
            "dec_centroid_err",
        ]
        skywin_colnames = [
            "ra_centroid_win",
            "dec_centroid_win",
            "ra_centroid_win_err",
            "dec_centroid_win_err",
        ]
        skypsf_colnames = [
            "ra_psf",
            "dec_psf",
            "ra_psf_err",
            "dec_psf_err",
        ]
        segm_colnames = [
            "bbox_xmin",
            "bbox_xmax",
            "bbox_ymin",
            "bbox_ymax",
            "segment_area",
        ]
        shape_colnames = [
            "semimajor",
            "semiminor",
            "fwhm",
            "ellipticity",
            "orientation_pix",
            "orientation_sky",
            "cxx",
            "cxy",
            "cyy",
            "kron_radius",
        ]
        nn_colnames = [
            "nn_label",
            "nn_distance",
        ]
        othershape_colnames = [
            "sharpness",
            "roundness1",
            "is_extended",
            "fluxfrac_radius_50",
        ]
        xypsf_colnames = [
            "x_psf",
            "x_psf_err",
            "y_psf",
            "y_psf_err",
        ]
        flag_columns = [
            "warning_flags",
            "image_flags",
        ]
        psf_flags_colnames = [
            "psf_flags",
            "psf_gof",
        ]

        det_colnames = []
        det_colnames.extend(segm_colnames)
        det_colnames.extend(shape_colnames)
        det_colnames.extend(nn_colnames)

        if self.cat_type in ("prompt", "forced_full"):
            colnames = []
            colnames.extend(base_colnames)
            colnames.extend(xywin_colnames)
            if self.fit_psf:
                colnames.extend(xypsf_colnames)
            colnames.extend(skybest_colnames)
            colnames.extend(sky_colnames)
            colnames.extend(skywin_colnames)
            if self.fit_psf:
                colnames.extend(skypsf_colnames)
            colnames.extend(det_colnames)
            colnames.extend(othershape_colnames)
            colnames.extend(self.flux_colnames)

            colnames.extend(flag_columns)
            if self.fit_psf:
                colnames.extend(psf_flags_colnames)

        elif self.cat_type == "forced_det":
            colnames = ["label"]  # needed to join the forced catalogs
            colnames.extend(base_colnames)
            colnames.extend(skybest_colnames)
            colnames.extend(sky_colnames)
            colnames.extend(shape_colnames)
            colnames.extend(nn_colnames)

        elif self.cat_type == "dr_det":
            colnames = []
            colnames.extend(base_colnames)
            colnames.extend(xywin_colnames)
            colnames.extend(skybest_colnames)
            colnames.extend(sky_colnames)
            colnames.extend(skywin_colnames)
            colnames.extend(det_colnames)
            colnames.extend(flag_columns)

        elif self.cat_type == "dr_band":
            colnames = self.band_colnames.copy()

        elif self.cat_type == "psf_matched":
            # Similar to dr_band but without the othershape_colnames
            # (sharpness, roundness1, is_extended, fluxfrac_radius_50)
            colnames = ["label"]
            colnames.extend(self.flux_colnames)

        return colnames

    def prefix_forced(self, catalog):
        """
        Prefix select columns of the catalog with ``"forced_"`` for the
        forced catalog.

        Parameters
        ----------
        catalog : `~astropy.table.Table`
            The source catalog to prefix (modified in place).

        Returns
        -------
        catalog : `~astropy.table.Table`
            The same catalog, returned for chaining.
        """
        # Prefix all columns (except "label") with "forced_"
        if self.cat_type == "forced_det":
            for colname in catalog.colnames:
                if colname != "label":
                    catalog.rename_column(colname, f"forced_{colname}")

        # Prefix the band_colnames with "forced_"
        if self.cat_type == "forced_full":
            for colname in [*self.band_colnames, "warning_flags"]:
                if colname in catalog.colnames and colname != "label":
                    catalog.rename_column(colname, f"forced_{colname}")

        return catalog
