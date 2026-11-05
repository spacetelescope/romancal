"""
Measure PSF and aperture photometry at input catalog positions.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import numpy as np
from roman_datamodels import datamodels

from romancal.datamodels.fileio import open_dataset
from romancal.source_catalog._calibrate import (
    convert_to_flux_density,
    make_catalog_model,
    make_image_mask,
    open_psf_model,
    pixel_scale_arcsec,
    subtract_background,
)
from romancal.source_catalog._utils import copy_model_arrays
from romancal.source_catalog._wcs_utils import pixel_area_map
from romancal.stpipe import RomanStep
from romancal.time_domain_catalog._time_domain_catalog import (
    add_pixel_positions,
    assemble_catalog,
    empty_catalog,
    measure_sources,
    read_forced_catalogs,
)

if TYPE_CHECKING:
    from typing import ClassVar

__all__ = ["TimeDomainCatalogStep"]

log = logging.getLogger(__name__)


class TimeDomainCatalogStep(RomanStep):
    """
    Measure PSF and aperture photometry of an image at the positions of
    sources in input catalogs.

    The image is masked, background subtracted, and converted to flux
    density exactly as in `SourceCatalogStep`. Each input source then
    gets a PSF fit with a free position, a PSF fit with the position
    fixed at the catalog position, and circular-aperture photometry at
    the catalog position.

    Parameters
    ----------
    input : str or `ImageModel`
        Path to an ASDF file, or an `ImageModel`.

    Other Parameters
    ----------------
    forced_catalogs : list of str
        Catalogs supplying the source positions. Any format readable by
        `astropy.table.Table.read` with ``ra`` and ``dec`` columns in
        degrees (ICRS).

    id_columns : list of str, optional
        The column of each catalog giving the source ID. A single value
        applies to all catalogs. Sources outside their skycell core are
        dropped when the column is ``flagged_spatial_id``.

    bkg_boxsize : float, optional
        Edge length, in arcsec, of the background mesh boxes.

    data_release_id : str, optional
        The data release identifier recorded in the catalog metadata.
        Not recorded if not given.

    suffix : str, optional
        Suffix appended to the output filename.  Default ``'fcat'``.
    """

    class_alias = "time_domain_catalog"

    reference_file_types: ClassVar = ["epsf"]

    spec = """
        forced_catalogs = force_list(default=list())  # catalogs of positions
        id_columns = force_list(default=list("flagged_spatial_id"))  # id column of each catalog
        bkg_boxsize = float(default=55.0)     # background mesh box size in arcsec
        data_release_id = string(default=None)  # data release identifier
        suffix = string(default='fcat')       # Default suffix for output files
    """

    def save_model(self, model, **kwargs):
        kwargs["ext"] = "parquet"
        return super().save_model(model, **kwargs)

    def process(self, dataset):
        input_model = open_dataset(dataset, update_version=self.update_version)
        if not isinstance(input_model, datamodels.ImageModel):
            raise TypeError("TimeDomainCatalogStep requires an ImageModel.")

        if len(self.forced_catalogs) == 0:
            raise ValueError("At least one forced catalog is required.")
        id_columns = list(self.id_columns)
        if len(id_columns) == 1:
            id_columns *= len(self.forced_catalogs)
        if len(id_columns) != len(self.forced_catalogs):
            raise ValueError(
                "id_columns must have one entry or one per forced catalog."
            )

        mask = make_image_mask(input_model)
        model = copy_model_arrays(input_model)

        cat_model = make_catalog_model(model, datamodels.TimeDomainSourceCatalogModel)
        cat_model.meta.forced_psf_catalogs = list(self.forced_catalogs)
        cat_model.meta.id_columns = id_columns
        if self.data_release_id is not None:
            cat_model.meta.data_release_id = self.data_release_id

        sources = read_forced_catalogs(self.forced_catalogs, id_columns)
        add_pixel_positions(sources, model)
        sources = sources[np.isfinite(sources["x"])]
        log.info("%d input sources land on the image", len(sources))

        if len(sources) == 0 or np.all(mask):
            log.warning("No sources to measure; returning an empty catalog.")
            cat_model.source_catalog = empty_catalog(cat_model)
            return cat_model

        bkg_boxsize_px = max(round(self.bkg_boxsize / pixel_scale_arcsec(model)), 1)
        subtract_background(model, mask, bkg_boxsize_px)
        area_map = pixel_area_map(model.meta.wcs, model.data.shape)
        convert_to_flux_density(model, area_map=area_map)

        ref_file = self.get_reference_file(input_model, "epsf")
        with open_psf_model(ref_file, input_model) as psf_ref:
            columns = measure_sources(model, mask, sources, psf_ref, area_map)

        cat_model.source_catalog = assemble_catalog(cat_model, columns)
        return cat_model
