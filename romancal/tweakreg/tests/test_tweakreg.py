from roman_datamodels import datamodels as rdm
from romancal.tweakreg.tweakreg_step import TweakRegStep
from roman_datamodels import maker_utils
import os
import csv
import asdf
import pytest


def load_base_image_wcs(input_dm):
    # create a base WCS using data from a romanisim simulated image
    # Note: this is the command used to create the simulated image:
    #   romanisim-make-image \
    #     --catalog ./rsim_cat_F158.ecsv \
    #     --radec 270.0 66.0 \
    #     --bandpass F158 --sca 7 \
    #     --usecrds --webbpsf --date 2026 1 1 \
    #     --level 1 l1-270-66-gaia-2016-sca7.asdf

    # update meta.wcsinfo
    input_dm.meta.wcsinfo.v2_ref = 0.42955128282521254
    input_dm.meta.wcsinfo.v3_ref = -0.2479976768255853
    input_dm.meta.wcsinfo.vparity = -1
    input_dm.meta.wcsinfo.ra_ref = 270.0
    input_dm.meta.wcsinfo.dec_ref = 66.0
    input_dm.meta.wcsinfo.roll_ref = 60

    # read saved WCS object and add it to meta
    full_path_to_wcs_file = os.path.join(
        os.path.dirname(__file__), "data/base_image_wcs.asdf"
    )
    asdf_file = asdf.open(full_path_to_wcs_file)
    wcs = asdf_file.tree["wcs"]
    input_dm.meta["wcs"] = wcs


def create_base_image_source_catalog(tmpdir, output_filename):
    # create a temp CSV file to be used as source catalog
    header = ["x", "y"]
    data = [
        [100.00, 50.00],
        [200.00, 100.00],
        [400.00, 200.00],
    ]
    output = os.path.join(tmpdir, output_filename)
    with open(output, "w", encoding="UTF8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(data)


def add_tweakreg_catalog_attribute(tmpdir, input_dm):
    # create and add a mock source detection catalog
    tweakreg_catalog_filename = "base_image_sources.csv"
    create_base_image_source_catalog(tmpdir, tweakreg_catalog_filename)
    input_dm.meta["tweakreg_catalog"] = os.path.join(tmpdir, tweakreg_catalog_filename)
    return input_dm


@pytest.fixture
def base_image(tmpdir):
    l2 = maker_utils.mk_level2_image(shape=(4088, 4088))
    # update wcsinfo and add WCS from a simulated image
    load_base_image_wcs(l2)
    l2im = rdm.ImageModel(l2)
    return l2im


def test_tweakreg_raises_attributeerror_on_missing_tweakreg_catalog(base_image):
    """Test that TweakReg raises an AttributeError if meta.tweakreg_catalog is missing."""
    img = base_image
    with pytest.raises(Exception) as exec_info:
        TweakRegStep.call([img])

    assert type(exec_info.value) == AttributeError


def test_tweakreg_returns_modelcontainer(tmpdir, base_image):
    """Test that TweakReg returns a ModelContainer."""
    img = base_image
    add_tweakreg_catalog_attribute(tmpdir, img)
    res = TweakRegStep.call([img])

    assert type(res) == rdm.ModelContainer


def test_tweakreg_updates_cal_step(tmpdir, base_image):
    """Test that TweakReg updates meta.cal_step with tweakreg = COMPLETE."""
    img = base_image
    add_tweakreg_catalog_attribute(tmpdir, img)
    res = TweakRegStep.call([img])

    assert hasattr(res[0].meta.cal_step, "tweakreg")
    assert res[0].meta.cal_step.tweakreg == "COMPLETE"


def test_tweakreg_updates_group_id(tmpdir, base_image):
    """Test that TweakReg updates 'group_id' with a non-zero length string."""
    img = base_image
    add_tweakreg_catalog_attribute(tmpdir, img)
    res = TweakRegStep.call([img])

    assert hasattr(res[0].meta, "group_id")
    assert len(res[0].meta.group_id) > 0
