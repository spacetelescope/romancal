from roman_datamodels import datamodels as rdm
from romancal.tweakreg.tweakreg_step import TweakRegStep
from roman_datamodels.testing import utils as testutil
import os
import csv
import asdf


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


def create_base_image_and_catalog(tmpdir):
    # data from romanisim simulated image 'l1-270-66-gaia-2016-sca7_cal.asdf'
    l2 = testutil.mk_level2_image(shape=(4088, 4088))

    l2.meta.wcsinfo.v2_ref = 0.42955128282521254
    l2.meta.wcsinfo.v3_ref = -0.2479976768255853
    l2.meta.wcsinfo.vparity = -1
    l2.meta.wcsinfo.ra_ref = 270.0
    l2.meta.wcsinfo.dec_ref = 66.0
    l2.meta.wcsinfo.roll_ref = 60

    full_path_to_wcs_file = os.path.join(
        os.path.dirname(__file__), "data/base_image_wcs.asdf"
    )
    asdf_file = asdf.open(full_path_to_wcs_file)
    wcs = asdf_file.tree["wcs"]
    l2.meta["wcs"] = wcs

    # add source detection catalog name
    tweakreg_catalog_filename = "base_image_sources.csv"
    create_base_image_source_catalog(tmpdir, tweakreg_catalog_filename)
    l2.meta["tweakreg_catalog"] = os.path.join(tmpdir, tweakreg_catalog_filename)

    l2im = rdm.ImageModel(l2)
    return l2im


def test_tweakreg(tmpdir):
    img = create_base_image_and_catalog(tmpdir)

    res = TweakRegStep.call([img])

    # TweakRegStep should return a ModelContainer
    assert type(res) == rdm.ModelContainer

    # TweakRegStep should create meta.cal_step.tweakreg
    assert hasattr(res[0].meta.cal_step, "tweakreg")
    assert res[0].meta.cal_step.tweakreg == "COMPLETE"

    # TweakRegStep should create at least one group_id
    assert hasattr(res[0].meta, "group_id")
    assert len(res[0].meta.group_id)

    print("DONE")
