"""Tests for L3→L3 multiband catalog metadata blending (RCAL-1438)."""

from astropy.time import Time
from roman_datamodels.datamodels import MosaicModel

from romancal.datamodels import ModelLibrary
from romancal.multiband_catalog._metadata import (
    blend_image_metadata,
    finalize_catalog_metadata,
)
from romancal.multiband_catalog._multiband_catalog import initialize_catalog_model


def _make_mosaic(
    *,
    optical_element="F184",
    time_first="2027-01-01T00:00:00",
    time_last="2027-01-01T01:00:00",
    time_mean="2027-01-01T00:30:00",
    exposure_time=300.0,
    max_exposure_time=400.0,
    observation_pass=1,
    observation_segment=1,
    ra_ref=270.0,
    filename="mosaic.asdf",
):
    """Build a MosaicModel with controlled L3 metadata for blend tests."""
    model = MosaicModel.create_fake_data(shape=(10, 10))
    model.meta.filename = filename
    model.meta.instrument.optical_element = optical_element
    model.meta.coadd_info.time_first = Time(time_first)
    model.meta.coadd_info.time_last = Time(time_last)
    model.meta.coadd_info.time_mean = Time(time_mean)
    model.meta.coadd_info.exposure_time = exposure_time
    model.meta.coadd_info.max_exposure_time = max_exposure_time
    model.meta.observation["pass"] = observation_pass
    model.meta.observation.segment = observation_segment
    model.meta.wcsinfo.ra_ref = ra_ref
    return model


def _blend_models(models, asn_table_name=None, product_name=None):
    """Initialize a catalog and blend the given L3 models."""
    library = ModelLibrary(list(models))
    if asn_table_name is not None or product_name is not None:
        asn = dict(library.asn)
        if asn_table_name is not None:
            asn["table_name"] = asn_table_name
        products = [dict(asn["products"][0])]
        products[0]["members"] = list(products[0]["members"])
        if product_name is not None:
            products[0]["name"] = product_name
        asn["products"] = products
        library._asn = asn

    cat_model = initialize_catalog_model(library, models[0])
    time_means = []
    exposure_times = []
    for model in models:
        blend_image_metadata(model, cat_model, time_means, exposure_times)
    finalize_catalog_metadata(cat_model, time_means, exposure_times)
    return cat_model


def test_multifilter_optical_element_is_none():
    """
    Purpose: Multiband catalogs must not report a single filter at top level.
    Top-level optical_element is None; each image_metas entry keeps its filter.
    """
    model_f184 = _make_mosaic(optical_element="F184")
    model_f158 = _make_mosaic(optical_element="F158")

    cat_model = _blend_models([model_f184, model_f158])

    assert cat_model.meta.instrument.optical_element is None
    assert len(cat_model.meta.image_metas) == 2
    assert cat_model.meta.image_metas[0]["instrument"]["optical_element"] == "F184"
    assert cat_model.meta.image_metas[1]["instrument"]["optical_element"] == "F158"


def test_single_filter_keeps_optical_element():
    """
    Purpose: A single-filter catalog keeps that filter as optical_element.
    """
    model_f184 = _make_mosaic(optical_element="F184")
    model_f184_b = _make_mosaic(optical_element="F184")

    cat_model = _blend_models([model_f184, model_f184_b])

    assert cat_model.meta.instrument.optical_element == "F184"


def test_coadd_info_blend_min_mean_max():
    """
    Purpose: L3 coadd_info blends as min/mean/max across input mosaics.
    time_first=min, time_last=max, time_mean=mean, exposure_time=mean,
    max_exposure_time=max.
    """
    model_a = _make_mosaic(
        time_first="2027-01-01T00:00:00",
        time_last="2027-01-01T02:00:00",
        time_mean="2027-01-01T01:00:00",
        exposure_time=200.0,
        max_exposure_time=250.0,
    )
    model_b = _make_mosaic(
        optical_element="F158",
        time_first="2027-01-02T00:00:00",
        time_last="2027-01-02T04:00:00",
        time_mean="2027-01-02T02:00:00",
        exposure_time=400.0,
        max_exposure_time=500.0,
    )

    cat_model = _blend_models([model_a, model_b])
    coadd = cat_model.meta.coadd_info

    assert coadd.time_first == Time("2027-01-01T00:00:00")
    assert coadd.time_last == Time("2027-01-02T04:00:00")
    assert (
        coadd.time_mean == Time(["2027-01-01T01:00:00", "2027-01-02T02:00:00"]).mean()
    )
    assert coadd.exposure_time == 300.0
    assert coadd.max_exposure_time == 500.0


def test_observation_mismatch_nulls_pass_keeps_shared_segment():
    """
    Purpose: Observation fields disagree→None; shared values are retained.
    Mismatched pass becomes None while matching segment stays set.
    """
    model_a = _make_mosaic(observation_pass=1, observation_segment=1)
    model_b = _make_mosaic(
        optical_element="F158", observation_pass=2, observation_segment=1
    )

    cat_model = _blend_models([model_a, model_b])

    assert cat_model.meta.observation["pass"] is None
    assert cat_model.meta.observation.segment == 1


def test_image_filename_uses_association_table_name():
    """
    Purpose: meta.image.filename records the association name, not a single
    input mosaic filename (stand-in for contributing images).
    """
    model_a = _make_mosaic(filename="f184_mosaic.asdf")
    model_b = _make_mosaic(optical_element="F158", filename="f158_mosaic.asdf")
    asn_name = "r00342_ref-hlwas-deep-26a_full_multiband_asn.json"

    cat_model = _blend_models(
        [model_a, model_b],
        asn_table_name=asn_name,
        product_name="r00342_ref-hlwas-deep-26a_full_007m42x23y55_cat",
    )

    assert cat_model.meta.image.filename == asn_name


def test_wcsinfo_keeps_first_input():
    """
    Purpose: wcsinfo is taken from the first input and is not disagree-blended.
    """
    model_a = _make_mosaic(ra_ref=10.0)
    model_b = _make_mosaic(optical_element="F158", ra_ref=99.0)

    cat_model = _blend_models([model_a, model_b])

    assert cat_model.meta.wcsinfo.ra_ref == 10.0
    # Per-image wcsinfo is still preserved on image_metas
    assert cat_model.meta.image_metas[0]["wcsinfo"]["ra_ref"] == 10.0
    assert cat_model.meta.image_metas[1]["wcsinfo"]["ra_ref"] == 99.0
