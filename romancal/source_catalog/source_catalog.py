"""
Module to calculate the source catalog.
"""

import logging
import re
from collections.abc import Callable
from dataclasses import dataclass

import astropy.units as u
import numpy as np
from astropy.table import QTable, Table
from astropy.utils import lazyproperty
from roman_datamodels.datamodels import ImageModel, MosaicModel
from roman_datamodels.dqflags import pixel

from romancal import __version__ as romancal_version
from romancal.skycell import skymap
from romancal.source_catalog.aperture import ApertureCatalog
from romancal.source_catalog.column_schema import CatalogSchema
from romancal.source_catalog.daofind import DAOFindCatalog
from romancal.source_catalog.neighbors import NNCatalog
from romancal.source_catalog.psf import PSFCatalog
from romancal.source_catalog.segment import SegmentCatalog
from romancal.source_catalog.unit_conversion import (
    validate_and_convert_to_flux_density,
)

from ._wcs_helpers import pixel_scale_angle_at_skycoord

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


@dataclass(frozen=True)
class MeasurementStep:
    """
    A single step in the measurement registry used to assemble the
    source catalog.

    Attributes
    ----------
    label : str
        Human-readable name used in log messages (e.g. ``"aperture
        photometry"``).

    factory : callable
        Zero-argument callable returning a sub-catalog instance for
        this step (e.g. an `ApertureCatalog`). Called only when
        ``condition()`` is True.

    store_as : str or None
        Attribute name on `RomanSourceCatalog` under which to retain the
        sub-catalog after measurement, or `None` if the sub-catalog is
        consumed during attribute copy and does not need to be kept.

    condition : callable
        Zero-argument callable returning `True` if this step should
        run for the current catalog (used to skip steps whose columns
        were not requested or which do not apply to the current catalog
        type).
    """

    label: str
    factory: Callable[[], object]
    store_as: str | None
    condition: Callable[[], bool]


class RomanSourceCatalog:
    """
    Class for the Roman source catalog.

    Parameters
    ----------
    model : `ImageModel` or `MosaicModel`
        The input data model. The image data is assumed to be background
        subtracted. For PSF-matched photometry in multiband catalogs,
        input the PSF-matched model.

    cat_model : `ImageSourceCatalogModel` or other catalog model
        The output catalog model.

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

    cat_type : str, optional
        The type of catalog to create. The default is 'prompt'. Allowed
        values are 'prompt', 'dr_det', 'dr_band', 'psf_matched',
        'forced_full', and 'forced_det'. This determines the columns
        in the output catalog. The 'dr_det' and 'dr_band' catalogs are
        band-specific catalogs for the multiband source detection. The
        'psf_matched' catalog is similar to 'dr_band' but excludes
        sharpness, roundness1, is_extended, and fluxfrac_radius_50
        properties for performance optimization.

    ee_spline : `~astropy.modeling.models.Spline1D` or `None`, optional
        The PSF aperture correction model, built from the reference
        file.

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
        cat_model,
        segment_img,
        convolved_data,
        kernel_fwhm,
        *,
        fit_psf=True,
        psf_model=None,
        mask=None,
        detection_cat=None,
        flux_unit="nJy",
        cat_type="prompt",
        ee_spline=None,
    ):
        if not isinstance(model, ImageModel | MosaicModel):
            raise ValueError("The input model must be an ImageModel or MosaicModel.")

        self.model = model  # input model is background-subtracted
        self.cat_model = cat_model
        self.segment_img = segment_img
        self.convolved_data = convolved_data
        self.kernel_fwhm = kernel_fwhm
        self.fit_psf = fit_psf
        self.psf_model = psf_model
        self.mask = mask
        self.detection_cat = detection_cat
        self.flux_unit_str = flux_unit
        self.flux_unit = u.Unit(self.flux_unit_str)
        self.cat_type = cat_type
        self.ee_spline = ee_spline

        self.n_sources = len(segment_img.labels)
        self.wcs = self.model.meta.wcs
        self.meta = {}
        self.aperture_cat = None

        # Define flux unit conversion factors
        if "photometry" in self.model.meta:
            self.l2_to_sb = self.model.meta.photometry.conversion_megajanskys
        else:
            self.l2_to_sb = 1.0
        self.sb_to_flux = (1.0 * (u.MJy / u.sr) * self._pixel_area).to(self.flux_unit)

        if self.fit_psf and self.psf_model is None:
            log.error(
                "PSF fitting is requested but no PSF reference model is provided. Skipping PSF photometry."
            )
            self.fit_psf = False

        # Output column schema (depends only on cat_type, fit_psf, and
        # the aperture flux column names — all known at this point).
        # Encapsulated in CatalogSchema for ease of maintenance.
        self._schema = CatalogSchema(
            cat_type=self.cat_type,
            fit_psf=self.fit_psf,
            aperture_flux_colnames=ApertureCatalog.aperture_flux_colnames_for_radii(),
        )

    def __len__(self):
        return self.n_sources

    @lazyproperty
    def _pixscale_angle(self):
        """
        The pixel scale in arcseconds and the angle in degrees measured
        counterclockwise from the positive x axis to the "North" axis of
        the celestial coordinate system.

        The pixel scale is returned as a Quantity in arcsec and the
        angle is returned as a Quantity in degrees.

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
        if (
            "photometry" in self.model.meta
            and self.model.meta.photometry.pixel_area > 0
        ):
            pixel_area = self.model.meta.photometry.pixel_area * u.sr
        else:
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
        replaced with a large negative value.

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
    def flagged_spatial_id(self):
        """
        The spatial index bit flag encoding the projection, skycell,
        and (x, y) pixel coordinate of the source and whether the
        object is inside the skycell core region.
        """
        bad_return = np.zeros(self.n_sources, dtype=np.int64)

        try:
            skycell_name = self.model.meta.wcsinfo.skycell_name
            pixel_scale = self.model.meta.wcsinfo.pixel_scale_ref
        except AttributeError:
            # L2 image or unrecognized schema: skycell-based spatial id
            # is not defined.
            log.warning(
                "meta.wcsinfo missing skycell_name/pixel_scale_ref; "
                "flagged_spatial_id will be zero.",
            )
            return bad_return

        try:
            sc = skymap.SkyCells.from_names([skycell_name])
        except KeyError:
            log.warning(
                f"Could not find skycell {skycell_name}; "
                "flagged_spatial_id will be zero.",
            )
            return bad_return

        core_indices = sc.cores_containing(np.array([self.ra, self.dec]).T)
        if len(core_indices) > 1:
            raise RuntimeError(
                f"Expected sources from at most one skycell, but found "
                f"{len(core_indices)} skycells in flagged_spatial_id "
                f"computation for skycell {skycell_name}."
            )
        in_core = np.zeros(len(self.ra), dtype="bool")
        if len(core_indices) == 1:
            this_cell_idx = next(iter(core_indices.keys()))
            core_indices = core_indices[this_cell_idx]
            in_core[core_indices] = True

        projection_idx = sc.projection_regions[0]
        pattern = r"x(\d+)y(\d+)"
        match = re.search(pattern, skycell_name)
        if not match:
            log.warning(f"Invalid skycell name: {skycell_name}")
            return bad_return
        skycell_x_idx = int(match.group(1))
        skycell_y_idx = int(match.group(2))

        def convert_to_pixel_idx(val):
            virtual_scale = 0.05  # virtual pixel scale used for id computation
            idx = (val.value * pixel_scale * 3600 / virtual_scale).astype("i4")
            return np.clip(idx, 0, 2**16)

        pixel_x_idx = convert_to_pixel_idx(self.x_centroid)
        pixel_y_idx = convert_to_pixel_idx(self.y_centroid)
        spatial_id = (
            ~in_core * 2**59
            + projection_idx * 2**46
            + skycell_y_idx * 2**39
            + skycell_x_idx * 2**32
            + pixel_y_idx * 2**16
            + pixel_x_idx
        )
        spatial_id = spatial_id.astype(np.int64)
        return spatial_id

    @lazyproperty
    def ra(self):
        """
        The best estimate of the right ascension (ICRS).
        """
        # FIXME: just return ra_centroid as a placeholder
        return self.ra_centroid.copy()

    @lazyproperty
    def dec(self):
        """
        The best estimate of the declination (ICRS).
        """
        # FIXME: just return dec_centroid as a placeholder
        return self.dec_centroid.copy()

    @lazyproperty
    def is_extended(self):
        """
        Boolean indicating whether the source is extended.
        """
        return self.aperture_cat.is_extended

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

        # Sources whose centroid pixel is not finite
        flags[~xymask] = pixel.DO_NOT_USE

        # Dispatch on data model type. L2 (ImageModel) has a ``dq``
        # array, while L3 mosaics have a ``weight`` array.
        if hasattr(self.model, "dq"):
            # L2 images: propagate the DO_NOT_USE flag from the DQ array.
            dqflags = self.model.dq[xyidx[:, 1], xyidx[:, 0]]
            flags[xymask] = dqflags & pixel.DO_NOT_USE
        else:
            # L3 mosaics: zero-weight pixels are unusable.
            mask = self.model.weight == 0
            flags[xymask] = mask[xyidx[:, 1], xyidx[:, 0]].astype(np.int32)

        return flags

    @lazyproperty
    def image_flags(self):
        """
        Data quality bit flag.

        Non-zero if a pixel within the segment was flagged in one of the
        input images.
        """
        return np.zeros(self.n_sources, dtype=np.int32)

    def _validate_and_convert_units(self):
        """
        Validate that model data arrays have compatible units and
        convert them to ``self.flux_unit`` if needed.
        """
        self.convolved_data = validate_and_convert_to_flux_density(
            self.model,
            self.convolved_data,
            flux_unit=self.flux_unit,
            l2_to_sb=self.l2_to_sb,
            sb_to_flux=self.sb_to_flux,
        )

    @property
    def column_names(self):
        """
        An ordered list of the output catalog column names.
        """
        return self._schema.column_names

    def _make_segment_cat(self):
        return SegmentCatalog(
            self.model,
            self.segment_img,
            self.convolved_data,
            self.mask,
            self._pixel_area,
            self._wcs_angle,
            detection_cat=self.detection_cat,
            requested_properties=self.column_names,
        )

    def _make_aperture_cat(self):
        return ApertureCatalog(
            self.model,
            self._xypos_finite,
            self._pixel_scale,
            ee_spline=self.ee_spline,
            requested_properties=self.column_names,
        )

    def _make_daofind_cat(self):
        return DAOFindCatalog(
            self.model.data,
            self._xypos_finite,
            self.kernel_fwhm,
            requested_properties=self.column_names,
        )

    def _make_nn_cat(self):
        return NNCatalog(
            self.label,
            self._xypos,
            self._xypos_finite,
            self._pixel_scale,
            requested_properties=self.column_names,
        )

    def _make_psf_cat(self):
        return PSFCatalog(
            self.model,
            self.psf_model,
            self._xypos,
            self.mask,
            requested_properties=self.column_names,
        )

    @property
    def _measurement_registry(self):
        """
        Ordered list of `MeasurementStep` entries used to build the
        catalog.

        Order matters: segment must be first because downstream steps
        depend on its outputs (``label``, ``x_centroid``, etc.) via
        ``self._xypos``. Aperture must come before steps that may
        reference its outputs.

        Returns
        -------
        result : list of `MeasurementStep`
            Each entry's ``label``, ``factory``, ``store_as`` and
            ``condition`` fields control how the corresponding
            sub-catalog is constructed and stored. See `MeasurementStep`
            for field semantics.
        """
        cols = set(self.column_names)
        return [
            MeasurementStep(
                label="segment properties",
                factory=self._make_segment_cat,
                store_as="segment_cat",
                condition=lambda: True,
            ),
            MeasurementStep(
                label="aperture photometry",
                factory=self._make_aperture_cat,
                store_as="aperture_cat",
                condition=lambda: self.cat_type != "forced_det",
            ),
            MeasurementStep(
                label="DAOFind properties",
                factory=self._make_daofind_cat,
                store_as=None,
                condition=lambda: bool({"sharpness", "roundness1"} & cols),
            ),
            MeasurementStep(
                label="nearest neighbor properties",
                factory=self._make_nn_cat,
                store_as=None,
                condition=lambda: any(col.startswith("nn_") for col in cols),
            ),
            MeasurementStep(
                label="PSF photometry",
                factory=self._make_psf_cat,
                store_as=None,
                condition=lambda: any("psf" in col for col in cols),
            ),
        ]

    def _apply_subcatalog(self, subcatalog, *, store_as=None):
        """
        Copy the named columns from a sub-catalog onto ``self`` and,
        optionally, also store the sub-catalog itself for later use.

        Parameters
        ----------
        subcatalog : object
            A sub-catalog instance with ``properties`` (iterable of
            attribute names to copy) and optionally a ``meta`` dict.

        store_as : str or `None`
            If given, also store the sub-catalog as ``self.<store_as>``.
        """
        meta = getattr(subcatalog, "meta", None)
        if meta is not None:
            self.meta.update(meta)
        for name in subcatalog.properties:
            setattr(self, name, getattr(subcatalog, name))
        if store_as is not None:
            setattr(self, store_as, subcatalog)

    def _run_measurements(self):
        """
        Run each enabled `MeasurementStep` in registry order.

        Order matters: each step may depend on attributes set by earlier
        steps (e.g. ``label`` and ``x_centroid`` from ``segment_cat``
        are needed before nearest-neighbor and PSF photometry). Steps
        whose ``condition()`` returns `False` are skipped silently.
        """
        for step in self._measurement_registry:
            if not step.condition():
                continue
            log.info(f"Calculating {step.label}")
            self._apply_subcatalog(step.factory(), store_as=step.store_as)

    def update_metadata(self):
        """
        Update the metadata dictionary with the package version
        information and aperture parameters.
        """
        ver_key = "versions"
        if "version" in self.meta:
            # Depends on photutils version
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

        # Reorder the metadata dictionary
        self.meta.pop("local_bkg_width", None)
        keys = ["date", "versions", "aperture_mask_method", "kron_params"]
        old_meta = self.meta.copy()
        self.meta = {key: old_meta.get(key, None) for key in keys}

        # Reformat the aperture radii for the metadata to remove
        # Quantity objects
        if self.aperture_cat is not None:
            aper_radii = self.aperture_cat.aperture_radii.copy()
            aper_radii["circle_arcsec"] = aper_radii.pop("circle").value
            aper_radii["annulus_arcsec"] = aper_radii.pop("annulus").value
            self.meta["aperture_radii"] = aper_radii

            if self.ee_spline:
                fractions = []
                for name in self.aperture_cat.fractions:
                    fraction = getattr(self.aperture_cat, name)
                    setattr(self, name, fraction)
                    fractions.append(fraction)

                self.meta["ee_fractions"] = np.array(fractions).astype(np.float32)

    def _assemble_catalog(self):
        """
        Build the output `~astropy.table.Table` from this catalog's
        per-source attributes.

        Returns
        -------
        catalog : `~astropy.table.Table`
            One column per name in `column_names`, with column
            descriptions sourced from ``self.cat_model``, ``forced_``
            prefixes applied where applicable, and metadata
            copied from `update_metadata`. The result is a plain
            `~astropy.table.Table` (not `QTable`) so columns are not
            `~astropy.units.Quantity`.
        """
        catalog = QTable()
        for column in self.column_names:
            catalog[column] = getattr(self, column)
            definition = self.cat_model.get_column_definition(column)
            catalog[column].info.description = definition["description"]
        self.update_metadata()
        catalog.meta.update(self.meta)

        # Prefix select columns with "forced_" for the forced catalog
        catalog = self._schema.prefix_forced(catalog)

        # Convert QTable to Table to avoid having Quantity columns
        return Table(catalog)

    @lazyproperty
    def catalog(self):
        """
        The final source catalog as an Astropy Table.
        """
        # Validate and convert units for all data arrays.
        self._validate_and_convert_units()

        # Run all enabled measurement steps in registry order.
        self._run_measurements()

        # Assemble the column values populated by the measurement steps
        # into the output Table.
        return self._assemble_catalog()
