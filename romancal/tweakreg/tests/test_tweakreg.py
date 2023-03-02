import pytest
from roman_datamodels import datamodels as rdm
from romancal.tweakreg.tweakreg_step import TweakRegStep
from roman_datamodels.testing import utils as testutil
import os
import csv

def create_base_image_source_catalog(tmp, output_filename):
    # create a temp CSV file to be used as source catalog
    header = ['x', 'y']
    data = [
        [100.00, 50.00],
        [200.00, 100.00],
        [400.00, 200.00],
    ]
    output = os.path.join(tmp, output_filename)
    with open(output, 'w', encoding='UTF8', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(data)


def create_base_image(tmp):
    l2 = testutil.mk_level2_image(shape=(100, 100))
    # data from romanisim simulated image 'l1-270-66-gaia-2016-sca7_cal.asdf'
    l2.meta.wcsinfo = {
        "v2_ref": 0.42955128282521254,
        "v3_ref": -0.2479976768255853,
        "vparity": -999999,
        "v3yangle": -999999,
        "ra_ref": 270.0,
        "dec_ref": 66.0,
        "roll_ref": 60,
        "s_region": "POLYGON IRCS 269.3318903230621 65.56866666048172 269.32578768154605 65.69246311613287 269.02457173246125 65.69201346248587 269.0333096074621 65.56870823657276 ",
    }
    # add source detection catalog name
    tweakreg_catalog_filename = "base_image_sources.csv"
    create_base_image_source_catalog(tmp, tweakreg_catalog_filename)
    l2.meta.tweakreg_catalog = os.path.join(tmp, tweakreg_catalog_filename)

    l2im = rdm.ImageModel(l2)
    return l2im

def test_shift():
    pass
