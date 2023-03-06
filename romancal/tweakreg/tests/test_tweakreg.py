from roman_datamodels import datamodels as rdm
from romancal.tweakreg.tweakreg_step import TweakRegStep
from roman_datamodels import maker_utils
import os
import csv
import pytest
from astropy import units as u
from gwcs import coordinate_frames as cf
from astropy import coordinates as coord
from gwcs import wcs
from romancal.assign_wcs import pointing


def create_dummy_wcs(input_dm):
    # create a dummy WCS object

    # create necessary transformations
    distortion = rdm.DistortionRefModel(
        maker_utils.mk_distortion()
    ).coordinate_distortion_transform
    tel2sky = pointing.v23tosky(input_dm)

    # create required frames
    detector = cf.Frame2D(name="detector", axes_order=(0, 1), unit=(u.pix, u.pix))
    v2v3 = cf.Frame2D(
        name="v2v3",
        axes_order=(0, 1),
        axes_names=("v2", "v3"),
        unit=(u.arcsec, u.arcsec),
    )
    world = cf.CelestialFrame(reference_frame=coord.ICRS(), name="world")

    # create pipeline
    pipeline = [
        wcs.Step(detector, distortion),
        wcs.Step(v2v3, tel2sky),
        wcs.Step(world, None),
    ]

    wcs_obj = wcs.WCS(pipeline)

    input_dm.meta["wcs"] = wcs_obj


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
    l2 = maker_utils.mk_level2_image(shape=(100, 100))
    # add a mock WCS object
    create_dummy_wcs(l2)
    # update vparity
    l2.meta.wcsinfo.vparity = -1
    l2_im = rdm.ImageModel(l2)
    return l2_im


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
