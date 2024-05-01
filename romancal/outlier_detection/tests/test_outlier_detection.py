import os
from pathlib import Path

import numpy as np
import pytest
from astropy.units import Quantity

from romancal.datamodels import ModelContainer
from romancal.outlier_detection import OutlierDetectionStep, outlier_detection


@pytest.fixture()
def clean_up_after_test():
    """
    Clean up working directory after test.
    """

    def _clean_up_working_directory(pattern="*.asdf"):
        """
        Clean up the working directory by removing files matching the specified pattern.
        Parameters
        ----------
        pattern : str, optional
            The pattern of the file pattern to match (default is "*.asdf")

        Returns
        -------
        None

        """
        for x in (Path(os.getcwd())).glob(pattern):
            if x.exists():
                x.unlink()

    return _clean_up_working_directory


@pytest.mark.parametrize(
    "input_models",
    [
        list(),
        "",
        None,
    ],
)
def test_outlier_raises_error_on_invalid_input_models(input_models, caplog):
    """Test that OutlierDetection logs out a WARNING if input is invalid."""

    OutlierDetectionStep.call(input_models)

    assert "WARNING" in [x.levelname for x in caplog.records]


def test_outlier_skips_step_on_invalid_number_of_elements_in_input(base_image):
    """Test that OutlierDetection skips processing when provided with an invalid number of elements in the input,
    and sets the appropriate metadata for the skipped step."""
    img = base_image()

    res = OutlierDetectionStep.call(ModelContainer([img]))

    assert all(x.meta.cal_step.outlier_detection == "SKIPPED" for x in res)


def test_outlier_skips_step_on_exposure_type_different_from_wfi_image(base_image):
    """
    Test if the outlier detection step is skipped when the exposure type is different from WFI image.
    """
    img_1 = base_image()
    img_1.meta.exposure.type = "WFI_PRISM"
    img_2 = base_image()
    img_2.meta.exposure.type = "WFI_PRISM"

    res = OutlierDetectionStep.call(ModelContainer([img_1, img_2]))

    assert all(x.meta.cal_step.outlier_detection == "SKIPPED" for x in res)


def test_outlier_valid_input_asn(tmp_path, base_image, create_mock_asn_file):
    """
    Test that OutlierDetection runs with valid ASN file as input.
    """
    img_1 = base_image()
    img_1.meta.filename = "img_1.asdf"
    img_1.save(tmp_path / "img_1.asdf")
    img_2 = base_image()
    img_1.meta.filename = "img_2.asdf"
    img_2.save(tmp_path / "img_2.asdf")

    asn_filepath = create_mock_asn_file(tmp_path)

    res = OutlierDetectionStep.call(
        asn_filepath,
        output_dir=tmp_path.as_posix(),
        in_memory=True,
        resample_data=False,
    )

    # assert step.skip is False
    assert all(x.meta.cal_step.outlier_detection == "COMPLETE" for x in res)


def test_outlier_valid_input_modelcontainer(tmp_path, base_image):
    """
    Test that OutlierDetection runs with valid ModelContainer as input.
    """
    img_1 = base_image()
    img_1.meta.filename = "img_1.asdf"
    img_2 = base_image()
    img_2.meta.filename = "img_2.asdf"

    mc = ModelContainer([img_1, img_2])

    res = OutlierDetectionStep.call(
        mc,
        in_memory=True,
        resample_data=False,
    )

    assert all(x.meta.cal_step.outlier_detection == "COMPLETE" for x in res)


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
            "in_memory": True,
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
    """
    Test parameter setting on initialization for OutlierDetection.
    """
    img_1 = base_image()
    img_1.meta.filename = "img_1.asdf"
    input_models = ModelContainer([img_1])

    step = outlier_detection.OutlierDetection(input_models, **pars)

    assert step.input_models == input_models
    assert step.outlierpars == pars
    assert step.make_output_path == pars["make_output_path"]
    assert step.resample_suffix == f"_outlier_{pars['resample_suffix']}.asdf"


def test_outlier_do_detection_write_files_to_custom_location(
    tmp_path, base_image, clean_up_after_test
):
    """
    Test that OutlierDetection can create files on disk in a custom location.
    """
    img_1 = base_image()
    img_1.meta.filename = "img_1.asdf"
    img_2 = base_image()
    img_2.meta.filename = "img_2.asdf"
    input_models = ModelContainer([img_1, img_2])

    outlier_step = OutlierDetectionStep()
    # set output dir for all files created by the step
    outlier_step.output_dir = tmp_path.as_posix()
    # make sure files are written out to disk
    outlier_step.in_memory = False

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
        "resample_data": False,
        "good_bits": 0,
        "allowed_memory": None,
        "in_memory": outlier_step.in_memory,
        "make_output_path": outlier_step.make_output_path,
        "resample_suffix": "i2d",
    }

    # meta.filename for the median image created by OutlierDetection.do_detection()
    median_path = tmp_path / "drizzled_median.asdf"

    outlier_files_path = [
        median_path,
    ]

    step = outlier_detection.OutlierDetection(input_models, **pars)
    step.do_detection()

    assert all(x.exists() for x in outlier_files_path)

    clean_up_after_test("*.asdf")


def test_outlier_do_detection_find_outliers(tmp_path, base_image, clean_up_after_test):
    """
    Test that OutlierDetection can find outliers.
    """
    img_1 = base_image()
    img_1.meta.filename = "img_1.asdf"
    img_2 = base_image()
    img_2.meta.filename = "img_2.asdf"

    # add outliers
    img_1_input_coords = np.array(
        [(5, 45), (25, 25), (45, 85), (65, 65), (85, 5)], dtype=[("x", int), ("y", int)]
    )
    img_2_input_coords = np.array(
        [(15, 25), (35, 5), (75, 65), (95, 45), (99, 5)], dtype=[("x", int), ("y", int)]
    )
    img_1.data[img_1_input_coords["x"], img_1_input_coords["y"]] = Quantity(
        100000, "DN / s"
    )
    img_2.data[img_2_input_coords["x"], img_2_input_coords["y"]] = Quantity(
        100000, "DN / s"
    )

    input_models = ModelContainer([img_1, img_2])

    outlier_step = OutlierDetectionStep()
    # set output dir for all files created by the step
    outlier_step.output_dir = tmp_path.as_posix()
    # make sure files are written out to disk
    outlier_step.in_memory = False

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
        "resample_data": False,
        "good_bits": 0,
        "allowed_memory": None,
        "in_memory": outlier_step.in_memory,
        "make_output_path": outlier_step.make_output_path,
        "resample_suffix": "i2d",
    }

    detection_step = outlier_detection.OutlierDetection
    step = detection_step(input_models, **pars)

    step.do_detection()

    # get flagged outliers coordinates from DQ array
    img_1_outlier_output_coords = np.where(step.input_models[0].dq > 0)

    # reformat output and input coordinates and sort by x coordinate
    outliers_output_coords = np.array(
        list(zip(*img_1_outlier_output_coords)), dtype=[("x", int), ("y", int)]
    )
    outliers_input_coords = np.concatenate((img_1_input_coords, img_2_input_coords))

    outliers_output_coords.sort(axis=0)
    outliers_input_coords.sort(axis=0)

    # assert all(outliers_input_coords == outliers_output_coords) doesn't work with python 3.9
    assert all(o == i for i, o in zip(outliers_input_coords, outliers_output_coords))

    clean_up_after_test("*.asdf")


def test_outlier_do_detection_do_not_find_outliers_in_identical_images(
    tmp_path, base_image, clean_up_after_test, caplog
):
    """
    Test that OutlierDetection does not flag any outliers in the DQ array if images are identical.
    """
    img_1 = base_image()
    img_1.meta.filename = "img_1.asdf"
    # add outliers
    img_1_input_coords = np.array(
        [(5, 45), (25, 25), (45, 85), (65, 65), (85, 5)], dtype=[("x", int), ("y", int)]
    )
    img_1.data[img_1_input_coords["x"], img_1_input_coords["y"]] = Quantity(
        100000, "DN / s"
    )

    img_2 = img_1.copy()
    img_2.meta.filename = "img_2.asdf"
    img_3 = img_1.copy()
    img_3.meta.filename = "img_3.asdf"

    input_models = ModelContainer([img_1, img_2, img_3])

    outlier_step = OutlierDetectionStep()
    # set output dir for all files created by the step
    outlier_step.output_dir = tmp_path.as_posix()
    # make sure files are written out to disk
    outlier_step.in_memory = False

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
        "resample_data": False,
        "good_bits": 0,
        "allowed_memory": None,
        "in_memory": outlier_step.in_memory,
        "make_output_path": outlier_step.make_output_path,
        "resample_suffix": "i2d",
    }

    detection_step = outlier_detection.OutlierDetection
    step = detection_step(input_models, **pars)

    step.do_detection()

    # assert that log shows no new outliers detected
    assert "New pixels flagged as outliers: 0 (0.00%)" in {
        x.message for x in caplog.records
    }
    # assert that DQ array has nothing flagged as outliers
    assert [np.count_nonzero(x.dq) for x in step.input_models] == [0, 0, 0]


@pytest.mark.parametrize(
    "input_type",
    [
        "ModelContainer",
        "ASNFile",
        "DataModelList",
        "ASDFFilenameList",
    ],
)
def test_skymatch_always_returns_modelcontainer_with_updated_datamodels(
    input_type,
    base_image,
    tmp_path,
    create_mock_asn_file,
):
    """Test that the OutlierDetectionStep always returns a ModelContainer
    with updated data models after processing different input types."""

    os.chdir(tmp_path)
    img_1 = base_image()
    img_1.meta.filename = "img_1.asdf"
    img_2 = base_image()
    img_2.meta.filename = "img_2.asdf"

    mc = ModelContainer([img_1, img_2])

    mc.save(dir_path=tmp_path)

    step_input_map = {
        "ModelContainer": mc,
        "ASNFile": create_mock_asn_file(
            tmp_path,
            members_mapping=[
                {"expname": img_1.meta.filename, "exptype": "science"},
                {"expname": img_2.meta.filename, "exptype": "science"},
            ],
        ),
        "DataModelList": [img_1, img_2],
        "ASDFFilenameList": [img_1.meta.filename, img_2.meta.filename],
    }

    step_input = step_input_map.get(input_type)

    res = OutlierDetectionStep.call(step_input)

    assert isinstance(res, ModelContainer)
    assert all(x.meta.cal_step.outlier_detection == "COMPLETE" for x in res)
