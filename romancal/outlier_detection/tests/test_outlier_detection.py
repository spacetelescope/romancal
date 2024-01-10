import os
from pathlib import Path

import pytest

from romancal.datamodels import ModelContainer
from romancal.outlier_detection import OutlierDetectionStep, outlier_detection


@pytest.mark.parametrize(
    "input_models",
    [
        list(),
        [""],
        "",
        "testing",
        {},
    ],
)
def test_outlier_raises_error_on_invalid_input_models(input_models):
    """Test that OutlierDetection raises an error if input is not a ModelContainer."""

    step = OutlierDetectionStep()

    step.process(input_models)

    assert step.skip is True


def test_outlier_skips_step_on_invalid_number_of_elements_in_input(base_image):
    img = base_image()

    step = OutlierDetectionStep()
    step.process(ModelContainer([img]))

    assert step.skip is True
    assert step.input_models[0].meta.cal_step.outlier_detection == "SKIPPED"


def test_outlier_skips_step_on_exposure_type_different_from_wfi_image(base_image):
    img_1 = base_image()
    img_1.meta.exposure.type = "WFI_PRISM"
    img_2 = base_image()
    img_2.meta.exposure.type = "WFI_PRISM"

    step = OutlierDetectionStep()
    step.process(ModelContainer([img_1, img_2]))

    assert step.input_models[0].meta.cal_step.outlier_detection == "SKIPPED"
    assert step.skip


def test_outlier_valid_input_datamodels(tmp_path, base_image):
    img_1 = base_image()
    img_1.meta.filename = "img_1.asdf"
    img_2 = base_image()
    img_2.meta.filename = "img_2.asdf"

    step = OutlierDetectionStep()
    step.process(ModelContainer([img_1, img_2]))

    assert step.skip is False
    assert img_1.meta.cal_step.outlier_detection == "COMPLETE"
    assert img_2.meta.cal_step.outlier_detection == "COMPLETE"

    # clean up
    [x.unlink() for x in Path(os.getcwd()).glob("*.asdf")]


def test_outlier_valid_input_asn(tmp_path, base_image, create_mock_asn_file):
    img_1 = base_image()
    img_1.meta.filename = "img_1.asdf"
    img_1.save(tmp_path / "img_1.asdf")
    img_2 = base_image()
    img_1.meta.filename = "img_2.asdf"
    img_2.save(tmp_path / "img_2.asdf")

    asn_filepath = create_mock_asn_file(tmp_path)

    step = OutlierDetectionStep()
    step.process(ModelContainer(asn_filepath))

    assert step.skip is False
    assert all(
        x.meta.cal_step.outlier_detection == "COMPLETE" for x in step.input_models
    )

    # clean up
    [x.unlink() for x in Path(os.getcwd()).glob("*.asdf")]


@pytest.mark.parametrize(
    "pars",
    [
        {
            "weight_type": "exptime",
            "pixfrac": 1.0,
            "kernel": "square",
            "fillval": "INDEF",
            "nlow": 0,
            "nhigh": 0,
            "maskpt": 0.7,
            "grow": 1,
            "snr": "4.0 3.0",
            "scale": "0.5 0.4",
            "backg": 0.0,
            "kernel_size": "7 7",
            "save_intermediate_results": False,
            "resample_data": True,
            "good_bits": 0,
            "allowed_memory": None,
            "in_memory": False,
            "make_output_path": None,
            "resample_suffix": "i2d",
        },
        {
            "weight_type": "exptime",
            "save_intermediate_results": True,
            "make_output_path": None,
            "resample_suffix": "some_other_suffix",
        },
    ],
)
def test_outlier_init_default_parameters(pars, base_image):
    img_1 = base_image()
    img_1.meta.filename = "img_1.asdf"
    input_models = ModelContainer([img_1])

    step = outlier_detection.OutlierDetection(input_models, **pars)

    assert step.input_models == input_models
    assert step.outlierpars == pars
    assert step.make_output_path == pars["make_output_path"]
    assert step.resample_suffix == f"_outlier_{pars['resample_suffix']}.asdf"

    # clean up
    [x.unlink() for x in Path(os.getcwd()).glob("*.asdf")]


def test_outlier_do_detection(base_image):
    img_1 = base_image()
    img_1.meta.filename = "img_1.asdf"
    img_2 = base_image()
    img_2.meta.filename = "img_2.asdf"
    input_models = ModelContainer([img_1, img_2])

    pars = {
        "weight_type": "exptime",
        "pixfrac": 1.0,
        "kernel": "square",
        "fillval": "INDEF",
        "nlow": 0,
        "nhigh": 0,
        "maskpt": 0.7,
        "grow": 1,
        "snr": "4.0 3.0",
        "scale": "0.5 0.4",
        "backg": 0.0,
        "kernel_size": "7 7",
        "save_intermediate_results": False,
        "resample_data": True,
        "good_bits": 0,
        "allowed_memory": None,
        "in_memory": False,
        "make_output_path": OutlierDetectionStep().make_output_path,
        "resample_suffix": "i2d",
    }

    blot_path_1 = Path(os.getcwd()) / img_1.meta.filename.replace(".asdf", "_blot.asdf")
    blot_path_2 = Path(os.getcwd()) / img_2.meta.filename.replace(".asdf", "_blot.asdf")
    median_path = Path(os.getcwd()) / "img_median.asdf"
    output_path = Path(os.getcwd()) / f"img_outlier_{pars['resample_suffix']}.asdf"

    outlier_files_path = [
        blot_path_1,
        blot_path_2,
        median_path,
        output_path,
    ]

    detection_step = outlier_detection.OutlierDetection
    step = detection_step(input_models, **pars)

    step.do_detection()

    assert all(x.exists() for x in outlier_files_path)

    # clean up
    [x.unlink() for x in Path(os.getcwd()).glob("*.asdf")]
