import numpy as np
import pytest
from astropy.modeling import models
from gwcs.wcstools import grid_from_bounding_box
from numpy.testing import assert_allclose
from roman_datamodels import datamodels as rdm
from roman_datamodels import stnode as st
from stcal.alignment.util import wcs_bbox_from_shape

from romancal.assign_wcs.assign_wcs_step import AssignWcsStep, load_wcs

DATA_SHAPE = (100, 100)
TEST_TRANSFORM = models.Shift(1) & models.Shift(2)


def create_image():
    l2 = rdm.ImageModel.create_fake_data(shape=DATA_SHAPE)

    l2.meta.cal_step = st.L2CalStep.create_fake_data()
    l2.meta.cal_logs = st.CalLogs.create_fake_data()

    l2.meta.wcsinfo.v2_ref = -503
    l2.meta.wcsinfo.v3_ref = -318
    l2.meta.wcsinfo.ra_ref = 156
    l2.meta.wcsinfo.dec_ref = 54.2
    l2.meta.wcsinfo.vparity = -1
    l2.meta.wcsinfo.roll_ref = 0.15

    return l2


def create_distortion():
    distortions = []

    dist = rdm.DistortionRefModel.create_fake_data()
    dist.coordinate_distortion_transform = TEST_TRANSFORM.copy()
    distortions.append(dist)

    dist = rdm.DistortionRefModel.create_fake_data()
    dist.coordinate_distortion_transform = TEST_TRANSFORM.copy()
    dist.coordinate_distortion_transform.bounding_box = wcs_bbox_from_shape(DATA_SHAPE)
    distortions.append(dist)

    return distortions


def create_step():
    def load_wcs_step(image, file_name):
        return load_wcs(image, {"distortion": file_name})

    def assign_wcs_step(image, file_name):
        return AssignWcsStep.call(image, override_distortion=file_name)

    return [load_wcs_step, assign_wcs_step]


@pytest.mark.parametrize("distortion", create_distortion())
@pytest.mark.parametrize("step", create_step())
def test_wcs(tmp_path, distortion, step):
    file_name = str(tmp_path / "distortion.asdf")
    dist = rdm.DistortionRefModel(distortion)
    dist.save(file_name)

    l2im = create_image()
    l2_wcs = step(l2im, file_name)

    assert l2_wcs.meta.wcs is not None
    assert l2_wcs.meta.cal_step.assign_wcs == "COMPLETE"

    x, y = grid_from_bounding_box(l2_wcs.meta.wcs.bounding_box)

    ra, dec = l2_wcs.meta.wcs(x, y)
    assert_allclose(ra, l2im.meta.wcsinfo.ra_ref, atol=2.3)
    assert_allclose(dec, l2im.meta.wcsinfo.dec_ref, atol=1.3)

    # verify inputs outside the bounding box return NaNs
    ra, dec = l2_wcs.meta.wcs(-5, 3)
    assert np.isnan(ra)
    assert np.isnan(dec)

    # check S_REGION length and format
    s_region_list = l2_wcs.meta.wcsinfo.s_region.split()
    assert len(s_region_list) == 10
    assert s_region_list[0].lower() == "polygon"
    assert s_region_list[1].lower() == "icrs"

    # check if BBOX is not None
    assert l2_wcs.meta.wcs.bounding_box is not None

    # check if footprint solution for each detector is within 10% of (RA_REF, DEC_REF)
    s_region_alpha_list = [
        float(x) for i, x in enumerate(s_region_list[2:]) if i % 2 == 0
    ]
    s_region_delta_list = [
        float(x) for i, x in enumerate(s_region_list[2:]) if i % 2 != 0
    ]
    assert_allclose(s_region_alpha_list, l2_wcs.meta.wcsinfo.ra_ref, rtol=1e-1, atol=0)
    assert_allclose(s_region_delta_list, l2_wcs.meta.wcsinfo.dec_ref, rtol=1e-1, atol=0)


@pytest.mark.parametrize("step", create_step())
def test_crds_getbestref(step):
    l2im = create_image()
    l2_wcs = AssignWcsStep.call(l2im, override_distortion="N/A")

    assert l2_wcs.meta.cal_step.assign_wcs == "SKIPPED"
