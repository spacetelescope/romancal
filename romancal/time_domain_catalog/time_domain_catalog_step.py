"""
Measure PSF and aperture photometry at input catalog positions.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import astropy.units as u
import numpy as np
from astropy.table import QTable, Table
from roman_datamodels import datamodels

from romancal.datamodels.fileio import open_dataset
from romancal.source_catalog._aperture import ApertureCatalog
from romancal.source_catalog._calibrate import (
    convert_to_flux_density,
    make_catalog_model,
    make_image_mask,
    open_psf_model,
    subtract_background,
)
from romancal.source_catalog._utils import copy_model_arrays
from romancal.source_catalog._wcs_utils import pixel_area_map
from romancal.source_catalog.psf import _PSFCatalog
from romancal.stpipe import RomanStep

if TYPE_CHECKING:
    from typing import ClassVar

__all__ = ["TimeDomainCatalogStep"]

log = logging.getLogger(__name__)

# Bit of flagged_spatial_id marking sources outside their skycell core
SKYCELL_OVERLAP_BIT = 2**59

# Half-width, in pixels, of the box over which DQ flags are combined
# into image_flags; matches the default PSF fit shape of (5, 5)
FLAG_HALF_WIDTH = 2

FLOATING_PSF_COLUMNS = [
    "ra_psf",
    "ra_psf_err",
    "dec_psf",
    "dec_psf_err",
    "x_psf",
    "x_psf_err",
    "y_psf",
    "y_psf_err",
    "psf_flux",
    "psf_flux_err",
    "psf_flags",
    "psf_gof",
]
FIXED_PSF_COLUMNS = [
    "forced_psf_flux",
    "forced_psf_flux_err",
    "forced_psf_flags",
    "forced_psf_gof",
]


def _aperture_columns():
    """The aperture column names, without the forced_ prefix."""
    names = []
    for name in ApertureCatalog.aperture_flux_colnames_for_radii():
        names.extend([name, f"{name}_err"])
    names.extend(["aper_bkg_flux", "aper_bkg_flux_err"])
    return names


COLUMNS = [
    "id",
    "source_catalog_index",
    "ra",
    "dec",
    "x",
    "y",
    "image_flags",
    "image_flags2",
    *FLOATING_PSF_COLUMNS,
    *FIXED_PSF_COLUMNS,
    *[f"forced_{name}" for name in _aperture_columns()],
]


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

    bkg_boxsize : int, optional
        Edge length, in pixels, of the background mesh boxes.

    suffix : str, optional
        Suffix appended to the output filename.  Default ``'fcat'``.
    """

    class_alias = "time_domain_catalog"

    reference_file_types: ClassVar = ["epsf"]

    spec = """
        forced_catalogs = force_list(default=list())  # catalogs of positions
        id_columns = force_list(default=list("flagged_spatial_id"))  # id column of each catalog
        bkg_boxsize = integer(default=1000)   # background mesh box size in pixels
        suffix = string(default='fcat')       # Default suffix for output files
    """

    def save_model(self, model, **kwargs):
        kwargs["ext"] = "parquet"
        kwargs.pop("idx", None)
        return super().save_model(model, **kwargs)

    def process(self, dataset):
        input_model = open_dataset(dataset, update_version=self.update_version)
        if not isinstance(input_model, datamodels.ImageModel):
            raise TypeError("TimeDomainCatalogStep requires an ImageModel.")

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

        sources = read_forced_catalogs(self.forced_catalogs, id_columns)
        sources = add_pixel_positions(sources, model)
        log.info("%d input sources land on the image", len(sources))

        if len(sources) == 0 or np.all(mask):
            log.warning("No sources to measure; returning an empty catalog.")
            cat_model.source_catalog = empty_catalog(cat_model)
            return cat_model

        subtract_background(model, mask, self.bkg_boxsize)
        area_map = pixel_area_map(model.meta.wcs, model.data.shape)
        convert_to_flux_density(model, area_map=area_map)

        columns = {name: sources[name] for name in sources.colnames}
        columns.update(image_flags(model, sources["x"], sources["y"]))

        xypos = np.transpose([sources["x"], sources["y"]])
        ref_file = self.get_reference_file(input_model, "epsf")
        with open_psf_model(ref_file, input_model) as psf_ref:
            log.info("Fitting PSFs with free positions")
            floating = _PSFCatalog(model, psf_ref, xypos, mask)
            log.info("Fitting PSFs with fixed positions")
            fixed = _PSFCatalog(
                model, psf_ref, xypos, mask, fixed_positions=True, prefix="forced_"
            )
        for name in FLOATING_PSF_COLUMNS:
            columns[name] = getattr(floating, name)
        for name in FIXED_PSF_COLUMNS:
            columns[name] = getattr(fixed, name)

        log.info("Calculating aperture photometry")
        aperture = ApertureCatalog(model, xypos, area_map, None, None)
        for name in _aperture_columns():
            columns[f"forced_{name}"] = getattr(aperture, name)

        cat_model.source_catalog = assemble_catalog(cat_model, columns)
        return cat_model


def read_forced_catalogs(filenames, id_columns):
    """
    Read and concatenate the catalogs of forced positions.

    Parameters
    ----------
    filenames : list of str
        The catalog files.

    id_columns : list of str
        The ID column of each catalog.

    Returns
    -------
    sources : `~astropy.table.Table`
        Columns ``id``, ``source_catalog_index``, ``ra``, and ``dec``.
    """
    columns = {
        "id": [np.zeros(0, dtype=np.int64)],
        "source_catalog_index": [np.zeros(0, dtype=np.int16)],
        "ra": [np.zeros(0)],
        "dec": [np.zeros(0)],
    }
    for index, (filename, id_column) in enumerate(
        zip(filenames, id_columns, strict=True)
    ):
        catalog = Table.read(filename)
        ids = np.asarray(catalog[id_column], dtype=np.int64)
        ra = _as_degrees(catalog["ra"])
        dec = _as_degrees(catalog["dec"])
        keep = np.isfinite(ra) & np.isfinite(dec)
        if id_column == "flagged_spatial_id":
            keep &= (ids & SKYCELL_OVERLAP_BIT) == 0
        log.info("Read %d of %d sources from %s", keep.sum(), len(keep), filename)
        columns["id"].append(ids[keep])
        columns["source_catalog_index"].append(
            np.full(keep.sum(), index, dtype=np.int16)
        )
        columns["ra"].append(ra[keep])
        columns["dec"].append(dec[keep])

    sources = Table({name: np.concatenate(arrays) for name, arrays in columns.items()})
    return sources


def _as_degrees(column):
    """Return a column as float64 degrees, with masked values as NaN."""
    if column.unit is not None:
        values = u.Quantity(column).to_value(u.deg)
    else:
        values = np.asarray(column, dtype=np.float64)
    values = np.array(values, dtype=np.float64)
    values[np.ma.getmaskarray(column)] = np.nan
    return values


def add_pixel_positions(sources, model):
    """
    Add pixel positions of the sources and drop those off the image.

    Parameters
    ----------
    sources : `~astropy.table.Table`
        Table with ``ra`` and ``dec`` columns in degrees.

    model : `ImageModel`
        The image, whose WCS maps sky to pixel positions.

    Returns
    -------
    sources : `~astropy.table.Table`
        The sources on the image, with ``x`` and ``y`` columns added.
    """
    x, y = model.meta.wcs.invert(
        sources["ra"].value, sources["dec"].value, with_bounding_box=False
    )
    x = np.atleast_1d(x)
    y = np.atleast_1d(y)
    ny, nx = model.data.shape
    on_image = (x >= -0.5) & (x < nx - 0.5) & (y >= -0.5) & (y < ny - 0.5)
    sources = sources[on_image]
    sources["x"] = x[on_image].astype(np.float32)
    sources["y"] = y[on_image].astype(np.float32)
    return sources


def image_flags(model, x, y):
    """
    Combine the DQ flags near each source.

    Parameters
    ----------
    model : `ImageModel`
        The image.

    x, y : array-like
        Source pixel positions.

    Returns
    -------
    flags : dict
        ``image_flags`` and ``image_flags2``: the bitwise or of ``dq``
        and ``dq2`` over a box around each source.
    """
    ny, nx = model.data.shape
    xi = np.round(np.asarray(x)).astype(int)
    yi = np.round(np.asarray(y)).astype(int)
    flags = {}
    for name, dqname in (("image_flags", "dq"), ("image_flags2", "dq2")):
        values = np.zeros(len(xi), dtype=np.uint32)
        dq = model.get(dqname, None)
        if dq is not None:
            offsets = range(-FLAG_HALF_WIDTH, FLAG_HALF_WIDTH + 1)
            for dy in offsets:
                for dx in offsets:
                    xx = xi + dx
                    yy = yi + dy
                    good = (xx >= 0) & (xx < nx) & (yy >= 0) & (yy < ny)
                    values[good] |= dq[yy[good], xx[good]].astype(np.uint32)
        # reinterpret rather than cast so that bit 31 survives
        flags[name] = values.view(np.int32)
    return flags


def assemble_catalog(cat_model, columns):
    """
    Build the output table in the standard column order.

    Parameters
    ----------
    cat_model : `TimeDomainSourceCatalogModel`
        The catalog model, which supplies the column descriptions.

    columns : dict
        Column values by name.

    Returns
    -------
    catalog : `~astropy.table.Table`
    """
    catalog = QTable()
    for name in COLUMNS:
        catalog[name] = columns[name]
        definition = cat_model.get_column_definition(name)
        if definition is not None:
            catalog[name].info.description = definition["description"]
    catalog["ra"].unit = u.deg
    catalog["dec"].unit = u.deg
    catalog["x"].unit = u.pix
    catalog["y"].unit = u.pix
    return Table(catalog)


def empty_catalog(cat_model):
    """
    Build an empty output table.

    Parameters
    ----------
    cat_model : `TimeDomainSourceCatalogModel`
        The catalog model, which supplies the column dtypes.

    Returns
    -------
    catalog : `~astropy.table.Table`
    """
    columns = {}
    for name in COLUMNS:
        definition = cat_model.get_column_definition(name)
        dtype = np.int32 if definition is None else definition["datatype"]
        columns[name] = np.zeros(0, dtype=dtype)
    return assemble_catalog(cat_model, columns)
