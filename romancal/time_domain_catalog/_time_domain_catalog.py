"""
Photometry at input catalog positions for the time domain catalog.
"""

import logging

import astropy.units as u
import numpy as np
from astropy.table import QTable, Table

from romancal.source_catalog._aperture import ApertureCatalog
from romancal.source_catalog._source_catalog import SKYCELL_OVERLAP_BIT
from romancal.source_catalog._utils import flag_column
from romancal.source_catalog.psf import _PSFCatalog

log = logging.getLogger(__name__)

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


def measure_sources(model, mask, sources, psf_ref, area_map):
    """
    Measure PSF and aperture photometry at the source positions.

    Parameters
    ----------
    model : `ImageModel`
        The background-subtracted image in flux density units.

    mask : 2D `~numpy.ndarray` of bool
        Bad-pixel mask; True pixels are excluded from the PSF fits.

    sources : `~astropy.table.Table`
        The sources, with ``x`` and ``y`` pixel positions.

    psf_ref : ePSF reference model
        The ePSF reference model, with jitter applied.

    area_map : `~astropy.units.Quantity`
        The per-pixel solid angle map.

    Returns
    -------
    columns : dict
        Column values by name, for every name in `COLUMNS`.
    """
    columns = {name: sources[name] for name in sources.colnames}
    columns.update(image_flags(model, sources["x"], sources["y"]))

    xypos = np.transpose([sources["x"], sources["y"]])
    log.info("Fitting PSFs with free positions")
    floating = _PSFCatalog(model, psf_ref, xypos, mask)
    for name in FLOATING_PSF_COLUMNS:
        columns[name] = getattr(floating, name)

    log.info("Fitting PSFs with fixed positions")
    fixed = _PSFCatalog(
        model, psf_ref, xypos, mask, fixed_positions=True, prefix="forced_"
    )
    for name in FIXED_PSF_COLUMNS:
        columns[name] = getattr(fixed, name)

    log.info("Calculating aperture photometry")
    aperture = ApertureCatalog(model, xypos, area_map, None, None)
    for name in _aperture_columns():
        columns[f"forced_{name}"] = getattr(aperture, name)

    return columns


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
    Add pixel positions of the sources, in place.

    Sources that land off the image get NaN positions.

    Parameters
    ----------
    sources : `~astropy.table.Table`
        Table with ``ra`` and ``dec`` columns in degrees.

    model : `ImageModel`
        The image, whose WCS maps sky to pixel positions.
    """
    # The bounding box makes gwcs drop sources outside the ra/dec range of
    # the image footprint before inverting.  Without it, distant sources
    # (~1 deg away, or near the antipode) can invert onto the image.
    x, y = model.meta.wcs.invert(
        sources["ra"].value, sources["dec"].value, with_bounding_box=True
    )
    x = np.atleast_1d(x)
    y = np.atleast_1d(y)
    ny, nx = model.data.shape
    on_image = (x >= -0.5) & (x < nx - 0.5) & (y >= -0.5) & (y < ny - 0.5)
    sources["x"] = np.where(on_image, x, np.nan).astype(np.float32)
    sources["y"] = np.where(on_image, y, np.nan).astype(np.float32)


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
        ``image_flags`` and ``image_flags2``: the uint32 bitwise or of
        ``dq`` and ``dq2`` over a box around each source.
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
        flags[name] = values
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
        if name in ("image_flags", "image_flags2"):
            catalog[name] = flag_column(columns[name], cat_model, name)
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
