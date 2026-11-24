"""Purpose: mosaic catalog metadata does not invent fake placeholder values."""

from astropy.time import Time
from roman_datamodels.datamodels import MosaicModel, MosaicSourceCatalogModel

from romancal.source_catalog._utils import copy_mosaic_meta


def test_copy_mosaic_meta_avoids_fake_placeholders():
    """
    Purpose: prd/sdf/survey placeholders are blank or taken from the input,
    not hardcoded junk like 8.8.8 / ?.
    """
    model = MosaicModel.create_fake_data(
        defaults={
            "meta": {
                "coadd_info": {
                    "time_first": Time("2024-01-01T00:00:00"),
                    "time_last": Time("2024-01-02T00:00:00"),
                    "time_mean": Time("2024-01-01T12:00:00"),
                    "exposure_time": 100.0,
                },
                "observation": {
                    "visit": 1,
                    "segment": 1,
                    "pass": 1,
                    "program": 1,
                },
                "instrument": {"optical_element": "F158"},
                "wcsinfo": {"skycell_name": "testcell"},
                "product_type": "l3",
            }
        }
    )
    # strip common fake defaults if present
    for key in ("prd_version", "sdf_software_version", "survey"):
        if key in model.meta:
            model.meta[key] = "?" if key == "survey" else ""

    cat = MosaicSourceCatalogModel.create_minimal({"meta": {}})
    # ensure basic/photometry containers exist for assignment
    if cat.meta.get("basic") is None:
        cat.meta["basic"] = {}
    if cat.meta.get("photometry") is None:
        cat.meta["photometry"] = {
            "conversion_megajanskys": None,
            "conversion_megajanskys_uncertainty": None,
            "pixel_area": None,
        }

    copy_mosaic_meta(model, cat)

    assert cat.meta.prd_version != "8.8.8"
    assert cat.meta.sdf_software_version != "7.7.7"
    assert cat.meta.basic.survey != "?"
    assert cat.meta.basic.survey == ""
    assert cat.meta.photometry.conversion_megajanskys is None
    assert cat.meta.photometry.pixel_area is None


def test_copy_mosaic_meta_preserves_real_versions():
    """Purpose: real input version metadata is copied through."""
    model = MosaicModel.create_fake_data(
        defaults={
            "meta": {
                "prd_version": "1.2.3",
                "sdf_software_version": "4.5.6",
                "survey": "HLWAS",
                "coadd_info": {
                    "time_first": Time("2024-01-01T00:00:00"),
                    "time_last": Time("2024-01-02T00:00:00"),
                    "time_mean": Time("2024-01-01T12:00:00"),
                    "exposure_time": 100.0,
                },
                "observation": {
                    "visit": 1,
                    "segment": 1,
                    "pass": 1,
                    "program": 1,
                },
                "instrument": {"optical_element": "F158"},
                "wcsinfo": {"skycell_name": "testcell"},
                "product_type": "l3",
            }
        }
    )
    cat = MosaicSourceCatalogModel.create_minimal({"meta": {}})
    if cat.meta.get("basic") is None:
        cat.meta["basic"] = {}
    if cat.meta.get("photometry") is None:
        cat.meta["photometry"] = {
            "conversion_megajanskys": None,
            "conversion_megajanskys_uncertainty": None,
            "pixel_area": None,
        }

    copy_mosaic_meta(model, cat)
    assert cat.meta.prd_version == "1.2.3"
    assert cat.meta.sdf_software_version == "4.5.6"
    assert cat.meta.basic.survey == "HLWAS"
