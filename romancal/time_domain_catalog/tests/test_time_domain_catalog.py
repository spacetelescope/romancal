"""Tests for the time domain catalog step."""

import numpy as np
import pytest
from astropy.table import Table
from numpy.testing import assert_allclose, assert_equal
from roman_datamodels.datamodels import TimeDomainSourceCatalogModel

from romancal.source_catalog._source_catalog import SKYCELL_OVERLAP_BIT
from romancal.source_catalog.source_catalog_step import SourceCatalogStep
from romancal.source_catalog.tests.test_source_catalog import (  # noqa: F401
    image_model,
)
from romancal.time_domain_catalog import TimeDomainCatalogStep
from romancal.time_domain_catalog._time_domain_catalog import (
    read_forced_catalogs,
)

BKG_BOXSIZE = 5.0


def _source_catalog(image_model):  # noqa: F811
    cat_model, _ = SourceCatalogStep.call(
        image_model,
        bkg_boxsize=BKG_BOXSIZE,
        kernel_fwhm=0.2,
        snr_threshold=5,
        npixels=10,
    )
    cat = cat_model.source_catalog
    cat["ra"] = cat["ra_psf"]
    cat["dec"] = cat["dec_psf"]
    cat[["label", "ra", "dec", "psf_flux"]].write("sources.ecsv")
    return cat


def test_time_domain_catalog(image_model, function_jail):  # noqa: F811
    """Fluxes at source catalog positions match the source catalog."""
    cat = _source_catalog(image_model)
    result = TimeDomainCatalogStep.call(
        image_model,
        forced_catalogs=["sources.ecsv"],
        id_columns=["label"],
        bkg_boxsize=BKG_BOXSIZE,
        data_release_id="r1",
        save_results=True,
        output_file="test.parquet",
    )
    assert isinstance(result, TimeDomainSourceCatalogModel)
    assert result.meta.forced_psf_catalogs == ["sources.ecsv"]
    assert result.meta.id_columns == ["label"]
    assert result.meta.data_release_id == "r1"

    # saving the catalog has validated it against the schema
    td = Table.read("test_fcat.parquet")
    assert_equal(td["id"], cat["label"])
    assert_equal(td["source_catalog_index"], 0)

    # The source positions are the PSF fit positions on the same image,
    # so both fits should reproduce the source catalog PSF fluxes
    assert_allclose(td["psf_flux"], cat["psf_flux"].value, rtol=1e-2)
    assert_allclose(td["forced_psf_flux"], cat["psf_flux"].value, rtol=1e-2)
    assert_allclose(td["x"], cat["x_psf"].value, atol=1e-3)


def test_no_sources_on_image(image_model, function_jail):  # noqa: F811
    """A catalog with no sources on the image gives an empty catalog."""
    Table({"id": [1], "ra": [0.0], "dec": [-89.0]}).write("far.ecsv")
    result = TimeDomainCatalogStep.call(
        image_model,
        forced_catalogs=["far.ecsv"],
        id_columns=["id"],
        save_results=True,
        output_file="test.parquet",
    )
    assert len(result.source_catalog) == 0
    assert "data_release_id" not in result.meta
    assert len(Table.read("test_fcat.parquet")) == 0


def test_no_forced_catalogs(image_model):  # noqa: F811
    """At least one forced catalog is required."""
    with pytest.raises(ValueError, match="forced catalog"):
        TimeDomainCatalogStep.call(image_model)


def test_read_forced_catalogs(function_jail):
    """Catalogs are concatenated and skycell overlap sources dropped."""
    Table(
        {
            "flagged_spatial_id": [1, 2 | SKYCELL_OVERLAP_BIT, 3],
            "ra": [1.0, 2.0, 3.0],
            "dec": [0.0, 0.0, 0.0],
        }
    ).write("a.ecsv")
    Table({"name": [7], "ra": [4.0], "dec": [0.0]}).write("b.ecsv")
    sources = read_forced_catalogs(["a.ecsv", "b.ecsv"], ["flagged_spatial_id", "name"])
    assert_equal(sources["id"], [1, 3, 7])
    assert_equal(sources["source_catalog_index"], [0, 0, 1])
    assert_equal(sources["ra"], [1.0, 3.0, 4.0])
    assert sources["id"].dtype == np.int64
