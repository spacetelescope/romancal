import os

import pytest
from astropy import units as u
from roman_datamodels import maker_utils
from roman_datamodels.datamodels import ImageModel

from romancal.outlier_detection import OutlierDetectionStep


@pytest.mark.skipif(
    os.environ.get("CI") == "true",
    reason=(
        "Roman CRDS servers are not currently available outside the internal network"
    ),
)
@pytest.mark.parametrize(
    "instrument, exptype",
    [
        ("WFI", "WFI_IMAGE"),
    ],
)
@pytest.mark.skip(reason="Outlier detection is not yet ready for testing.")
def test_outlier_detection_step_interface(instrument, exptype):
    """Test that the basic inferface works for data requiring outlier detection"""

    shape = (20, 20)

    wfi_image = maker_utils.mk_level2_image(shape=shape)
    wfi_image.meta.instrument.name = instrument
    wfi_image.meta.instrument.detector = "WFI01"
    wfi_image.meta.instrument.optical_element = "F158"
    wfi_image.meta.exposure.type = exptype
    wfi_image.data[3:5, 7:10] = 4 * u.electron / u.s

    wfi_image_model = ImageModel(wfi_image)

    result = OutlierDetectionStep.call(wfi_image_model)

    assert (result.data == wfi_image.data).all()
    # Uncomment after adding outlier_detection added to RAD & RDM
    # assert result.meta.cal_step.outlier_detection == "COMPLETE"
