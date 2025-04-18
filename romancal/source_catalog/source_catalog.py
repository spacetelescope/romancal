"""
Module to calculate the source catalog.
"""

import logging

import astropy.units as u
import numpy as np
from astropy.table import QTable, Table
from astropy.utils import lazyproperty
from roman_datamodels.datamodels import ImageModel, MosaicModel
from roman_datamodels.dqflags import pixel
from stpsf import __version__ as stpsf_version

from romancal import __version__ as romancal_version
from romancal.source_catalog.aperture import ApertureCatalog
from romancal.source_catalog.daofind import DAOFindCatalog
from romancal.source_catalog.neighbors import NNCatalog
from romancal.source_catalog.psf import PSFCatalog
from romancal.source_catalog.segment import SegmentCatalog

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
        The full-width at half-maximum (FWHM) of the DAOFind 2D Gaussian
        kernel. This kernel is used to calculate the DAOFind sharpness
        and roundness properties. DAOFind uses a special kernel that
        sums to zero.

    fit_psf : bool, optional
        Whether to fit a PSF model to the sources.

    mask : 2D `~numpy.ndarray` or `None`, optional
        A 2D boolean mask image with the same shape as the input data.
        This mask is used for PSF photometry. The mask should be the
        same one used to create the segmentation image.

    detection_cat : `None` or `RomanSourceCatalog`, optional
        A `RomanSourceCatalog` object for the detection image. The
        segmentation image used to create the detection catalog must
        be the same one input to ``segment_image``. If input, then the
        detection catalog source centroids and morphological/shape
        properties will be returned instead of calculating them from
        the input ``data``. The detection catalog centroids and shape
        properties will also be used to perform aperture photometry
        (i.e., circular and Kron).

    flux_unit : str, optional
        The unit of the flux density. Default is 'nJy'.

    Notes
    -----
    ``model.err`` is assumed to be the total error array corresponding
    to the input science ``model.data`` array. It is assumed to include
    all sources of error, including the Poisson error of the sources. It
    must also have the same shape and units as the science data array.
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
        self.kernel_fwhm = kernel_fwhm
        self.fit_psf = fit_psf
        self.mask = mask
        self.detection_cat = detection_cat
        self.flux_unit = flux_unit

        self.n_sources = len(segment_img.labels)
        self.wcs = self.model.meta.wcs
        self.meta = {}

        # define flux unit conversion factors
        self.l2_to_sb = self.model.meta.photometry.conversion_megajanskys
        self.sb_to_flux = (1.0 * (u.MJy / u.sr) * self._pixel_area).to(
            u.Unit(self.flux_unit)
        )

    def __len__(self):
        return self.n_sources

    def convert_l2_to_sb(self):
        """
        Convert level-2 data from units of DN/s to MJy/sr (surface
        brightness).
        """
        # the conversion in done in-place to avoid making copies of the data;
        # use dictionary syntax to set the value to avoid on-the-fly validation
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
        # use dictionary syntax to set the value to avoid on-the-fly validation
        self.model["data"] *= self.sb_to_flux.value
        self.model["data"] <<= self.sb_to_flux.unit
        self.model["err"] *= self.sb_to_flux.value
        self.model["err"] <<= self.sb_to_flux.unit
        if self.convolved_data is not None:
            self.convolved_data *= self.sb_to_flux.value
            self.convolved_data <<= self.sb_to_flux.unit

    @lazyproperty
    def _pixscale_angle(self):
        """
        The pixel scale in arcseconds and the angle in degrees measured
        counterclockwise from the positive x axis to the "North" axis of
        the celestial coordinate system.

        The pixel is returns as a Quantity in arcsec and the angle is
        returned as a Quantity in degrees.

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
    def _pixel_area(self):
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
            xycen = (self.detection_cat.x_centroid, self.detection_cat.y_centroid)
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

    def calc_segment_properties(self):
        """
        Calculate the segment-based properties provided by
        `~photutils.segmentation.SourceCatalog`.

        The results are set as dynamic attributes on the class instance.
        """
        segment_cat = SegmentCatalog(
            self.model,
            self.segment_img,
            self.convolved_data,
            self._pixel_area,
            self._wcs_angle,
            detection_cat=self.detection_cat,
            flux_unit=self.flux_unit,
        )

        self.meta.update(segment_cat.meta)
        for name in segment_cat.names:
            setattr(self, name, getattr(segment_cat, name))

        # needed for detection_cat
        self.segment_cat = segment_cat

    def calc_aperture_photometry(self):
        """
        Calculate aperture photometry.

        The results are set as dynamic attributes on the class instance.
        """
        aperture_cat = ApertureCatalog(
            self.model,
            self._pixel_scale,
            self._xypos_finite,
        )
        for name in aperture_cat.names:
            setattr(self, name, getattr(aperture_cat, name))

        # needed to get aperture flux column names and descriptions
        self.aperture_cat = aperture_cat

    def calc_psf_photometry(self):
        """
        Perform PSF photometry on the sources.

        The results are set as dynamic attributes on the class instance.
        """
        psf_cat = PSFCatalog(self.model, self._xypos, self.mask)
        for name in psf_cat.names:
            setattr(self, name, getattr(psf_cat, name))

    def calc_daofind_properties(self):
        """
        Calculate the DAOFind sharpness and roundness1 statistics.

        The results are set as dynamic attributes on the class instance.
        """
        daofind_cat = DAOFindCatalog(
            self.model.data, self._xypos_finite, self.kernel_fwhm
        )
        for name in daofind_cat.names:
            setattr(self, name, getattr(daofind_cat, name))

    def calc_nn_properties(self):
        """
        Calculate the nearest neighbor properties.

        The results are set as dynamic attributes on the class instance.
        """
        nn_cat = NNCatalog(
            self.label, self._xypos, self._xypos_finite, self._pixel_scale
        )
        for name in nn_cat.names:
            setattr(self, name, getattr(nn_cat, name))

    @lazyproperty
    def is_extended(self):
        """
        Boolean indicating whether the source is extended.
        """
        # TODO: replace with filter-dependent ee_ratio from reference file
        # ee_ratio = ee_fraction_04 / ee_fraction_02
        ee_ratio = 1.4  # F213 placeholder value
        return self.aper04_flux > (self.aper02_flux * 1.1 * ee_ratio)

    @lazyproperty
    def warning_flags(self):
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

    def update_metadata(self):
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
        aper_radii = self.aperture_cat.aperture_radii.copy()
        aper_radii["circle_arcsec"] = aper_radii.pop("circle").value
        aper_radii["annulus_arcsec"] = aper_radii.pop("annulus").value
        self.meta["aperture_radii"] = aper_radii

    @lazyproperty
    def column_descriptions(self):
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
        col["bbox_xmin"] = "Column index of the left edge of the source bounding box"
        col["bbox_xmax"] = "Column index of the right edge of the source bounding box"
        col["bbox_ymin"] = "Row index of the bottom edge of the source bounding box"
        col["bbox_ymax"] = "Row index of the top edge of the source bounding box"
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
        col["segment_area"] = "Area of the source segment"
        col["kron_flux"] = "Flux within the elliptical Kron aperture"
        col["aper_bkg_flux"] = "The local background estimate for aperture photometry"

        col["x_psf"] = "Column position of the source from PSF fitting (0 indexed)"
        col["y_psf"] = "Row position of the source from PSF fitting (0 indexed)"
        col["psf_flux"] = "Total PSF flux"
        col["psf_flags"] = "PSF fitting bit flags"

        col["warning_flags"] = "Warning bit flags"
        col["is_extended"] = "Flag indicating whether the source is extended"
        col["sharpness"] = "The DAOFind sharpness statistic"
        col["roundness"] = "The DAOFind roundness1 statistic"
        col["nn_label"] = "The label number of the nearest neighbor in this skycell"
        col["nn_dist"] = "The distance to the nearest neighbor in this skycell"

        # add the aperture flux column descriptions
        col.update(self.aperture_cat.aperture_flux_descriptions)

        # add the "*_err" column descriptions
        for column in self.column_names:
            if column.endswith("_err"):
                base_column = column.replace("_err", "")
                col[column] = f"Uncertainty in {base_column}"

        return col

    @lazyproperty
    def column_names(self):
        """
        An ordered list of the output catalog column names.
        """
        # define the aperture flux column names
        aper_colnames = []
        for colname in self.aperture_cat.aperture_flux_colnames:
            aper_colnames.append(colname)
            aper_colnames.append(f"{colname}_err")

        # define the flux column names
        flux_colnames = [
            "aper_bkg_flux",
            "aper_bkg_flux_err",
        ]
        flux_colnames.extend(aper_colnames)
        if self.fit_psf:
            flux_colnames.extend(
                [
                    "psf_flux",
                    "psf_flux_err",
                ]
            )
        flux_colnames.extend(
            [
                "segment_flux",
                "segment_flux_err",
                "kron_flux",
                "kron_flux_err",
            ]
        )

        if self.detection_cat is None:
            colnames = [
                "label",
                "x_centroid",
                "y_centroid",
            ]
            if self.fit_psf:
                colnames.extend(
                    [
                        "x_psf",
                        "x_psf_err",
                        "y_psf",
                        "y_psf_err",
                    ]
                )
            colnames.extend(
                [
                    "ra_centroid",
                    "dec_centroid",
                    "bbox_xmin",
                    "bbox_xmax",
                    "bbox_ymin",
                    "bbox_ymax",
                    "semimajor_sigma",
                    "semiminor_sigma",
                    "ellipticity",
                    "orientation_pix",
                    "orientation_sky",
                    "segment_area",
                    "is_extended",
                    "sharpness",
                    "roundness",
                    "nn_label",
                    "nn_dist",
                ]
            )
            colnames.extend(flux_colnames)
            colnames.append("warning_flags")
            if self.fit_psf:
                colnames.append("psf_flags")

        else:
            # these are band-specific columns for the multiband catalog
            colnames = ["label"]  # label is needed to join the filter catalogs
            colnames.extend(flux_colnames)

        return colnames

    @lazyproperty
    def catalog(self):
        """
        The final source catalog as an Astropy Table.
        """
        # convert data to flux units
        if isinstance(self.model, ImageModel):
            self.convert_l2_to_sb()
        self.convert_sb_to_flux_density()

        # make measurements - the order of these calculations is important
        self.calc_segment_properties()
        self.calc_aperture_photometry()
        self.calc_daofind_properties()
        self.calc_nn_properties()
        if self.fit_psf:
            self.calc_psf_photometry()

        # put the measurements into a Table
        catalog = QTable()
        for column in self.column_names:
            catalog[column] = getattr(self, column)
            descrip = self.column_descriptions.get(column, None)
            catalog[column].info.description = descrip
        self.update_metadata()
        catalog.meta.update(self.meta)

        # convert QTable to Table to avoid having Quantity columns
        catalog = Table(catalog)

        return catalog
