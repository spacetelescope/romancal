import numpy as np
import pytest

from romancal.datamodels import ModelLibrary
from romancal.outlier_detection import OutlierDetectionStep


def test_outlier_skips_step_on_invalid_number_of_elements_in_input(base_image):
    """Test that OutlierDetection skips processing when provided with an invalid number of elements in the input,
    and sets the appropriate metadata for the skipped step."""
    img = base_image()

    res = OutlierDetectionStep().run(ModelLibrary([img]))

    with res:
        for i, m in enumerate(res):
            assert m.meta.cal_step.outlier_detection == "SKIPPED"
            res.shelve(m, i, modify=False)


def test_outlier_raises_exception_on_exposure_type_different_from_wfi_image(base_image):
    """
    Test if the outlier detection step raises an Exception for non-image inputs.
    """
    img_1 = base_image()
    img_1.meta.exposure.type = "WFI_PRISM"
    img_2 = base_image()
    img_2.meta.exposure.type = "WFI_PRISM"

    with pytest.raises(ValueError, match="only supports WFI_IMAGE"):
        OutlierDetectionStep().run(ModelLibrary([img_1, img_2]))


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

    res = OutlierDetectionStep(
        output_dir=tmp_path.as_posix(),
        in_memory=True,
        resample_data=False,
    ).run(asn_filepath)

    # assert step.skip is False
    with res:
        for i, m in enumerate(res):
            assert m.meta.cal_step.outlier_detection == "COMPLETE"
            res.shelve(m, i, modify=False)


def test_outlier_valid_input_modelcontainer(tmp_path, base_image):
    """
    Test that OutlierDetection runs with valid ModelLibrary as input.
    """
    img_1 = base_image()
    img_1.meta.filename = "img_1.asdf"
    img_2 = base_image()
    img_2.meta.filename = "img_2.asdf"

    library = ModelLibrary([img_1, img_2])

    res = OutlierDetectionStep(
        in_memory=True,
        resample_data=False,
    ).run(library)

    with res:
        for i, m in enumerate(res):
            assert m.meta.cal_step.outlier_detection == "COMPLETE"
            res.shelve(m, i, modify=False)


def test_outlier_do_detection_write_files_to_custom_location(tmp_path, base_image):
    """
    Test that OutlierDetection can create files on disk in a custom location.
    """
    img_1 = base_image()
    img_1.meta.filename = "img1_cal.asdf"
    img_1.meta.background.level = 0
    img_2 = base_image()
    img_2.meta.filename = "img2_cal.asdf"
    img_2.meta.background.level = 0
    input_models = ModelLibrary([img_1, img_2])

    outlier_step = OutlierDetectionStep(
        in_memory=False,
        save_intermediate_results=True,
    )
    # set output dir for all files created by the step
    outlier_step.output_dir = tmp_path.as_posix()

    outlier_step.run(input_models)

    # meta.filename for the median image created by OutlierDetection.do_detection()
    median_path = tmp_path / "drizzled_median.asdf"

    outlier_files_path = [
        median_path,
    ]

    assert all(x.exists() for x in outlier_files_path)


@pytest.mark.parametrize("on_disk", (True, False))
def test_find_outliers(tmp_path, base_image, on_disk):
    """
    Test that OutlierDetection can find outliers.
    """
    cr_value = 100
    source_value = 10
    err_value = 10  # snr=1

    imgs = []
    for i in range(3):
        img = base_image()
        img.data[42, 72] = source_value
        img.err[:] = err_value
        img.meta.filename = str(tmp_path / f"img{i}_suffix.asdf")
        img.meta.observation.observation_id = str(i)
        img.meta.background.level = 0
        # populate var_rnoise with non-zero values
        img.var_rnoise = np.full(img.data.shape, 1.0, dtype=np.float32)
        imgs.append(img)

    # add outliers
    img_0_input_coords = np.array(
        [[5, 25, 45, 65, 85], [45, 25, 85, 65, 5]], dtype=np.int64
    )
    img_1_input_coords = np.array(
        [[15, 35, 75, 95, 99], [25, 5, 65, 45, 5]], dtype=np.int64
    )

    imgs[0].data[img_0_input_coords[0], img_0_input_coords[1]] = cr_value
    imgs[1].data[img_1_input_coords[0], img_1_input_coords[1]] = cr_value

    if on_disk:
        # write out models and asn to disk
        for img in imgs:
            img.save(img.meta.filename)
        asn_dict = {
            "asn_id": "a3001",
            "asn_pool": "none",
            "products": [
                {
                    "name": "test_asn",
                    "members": [
                        {"expname": m.meta.filename, "exptype": "science"} for m in imgs
                    ],
                }
            ],
        }
        input_models = ModelLibrary(asn_dict, on_disk=True)
    else:
        input_models = ModelLibrary(imgs)

    outlier_step = OutlierDetectionStep()
    # set output dir for all files created by the step
    outlier_step.output_dir = tmp_path.as_posix()
    # make sure files are written out to disk
    outlier_step.in_memory = not on_disk

    result = outlier_step.run(input_models)

    expected_crs = [img_0_input_coords, img_1_input_coords, None]
    with result:
        for cr_coords, flagged_img in zip(expected_crs, result, strict=False):
            if cr_coords is None:
                assert not np.any(flagged_img.dq > 0)
            else:
                flagged_coords = np.where(flagged_img.dq > 0)
                np.testing.assert_equal(cr_coords, flagged_coords)
            result.shelve(flagged_img, modify=False)


def test_identical_images(tmp_path, base_image, caplog):
    """
    Test that OutlierDetection does not flag any outliers in the DQ array if images are identical.
    """
    img_1 = base_image()
    img_1.meta.filename = "img1_suffix.asdf"
    img_1.meta.background.level = 0
    # add outliers
    img_1_input_coords = np.array(
        [(5, 45), (25, 25), (45, 85), (65, 65), (85, 5)], dtype=[("x", int), ("y", int)]
    )
    img_1.data[img_1_input_coords["x"], img_1_input_coords["y"]] = 100000

    img_2 = img_1.copy()
    img_2.meta.filename = "img2_suffix.asdf"
    img_3 = img_1.copy()
    img_3.meta.filename = "img3_suffix.asdf"

    input_models = ModelLibrary([img_1, img_2, img_3])

    outlier_step = OutlierDetectionStep()
    # set output dir for all files created by the step
    outlier_step.output_dir = tmp_path.as_posix()
    # make sure files are written out to disk
    outlier_step.in_memory = False

    result = outlier_step.run(input_models)

    # assert that log shows no new outliers detected
    assert "0 pixels marked as outliers" in {x.message for x in caplog.records}
    # assert that DQ array has nothing flagged as outliers
    with result:
        for i, model in enumerate(result):
            assert np.count_nonzero(model.dq) == 0
            result.shelve(model, i)


@pytest.mark.parametrize(
    "input_type",
    [
        "ModelLibrary",
        "ASNFile",
        "DataModelList",
        "ASDFFilenameList",
    ],
)
def test_outlier_detection_always_returns_modelcontainer_with_updated_datamodels(
    input_type,
    base_image,
    create_mock_asn_file,
    function_jail,
):
    """Test that the OutlierDetectionStep always returns a ModelLibrary
    with updated data models after processing different input types."""

    img_1 = base_image()
    img_1.meta.filename = "img_1.asdf"
    img_1.data *= img_1.meta.photometry.conversion_megajanskys
    img_2 = base_image()
    img_2.meta.filename = "img_2.asdf"
    img_2.data *= img_2.meta.photometry.conversion_megajanskys

    library = ModelLibrary([img_1, img_2])
    library._save(function_jail)

    step_input_map = {
        "ModelLibrary": library,
        "ASNFile": create_mock_asn_file(
            function_jail,
            members_mapping=[
                {"expname": img_1.meta.filename, "exptype": "science"},
                {"expname": img_2.meta.filename, "exptype": "science"},
            ],
        ),
        "DataModelList": [img_1, img_2],
        "ASDFFilenameList": [img_1.meta.filename, img_2.meta.filename],
    }

    step_input = step_input_map.get(input_type)

    res = OutlierDetectionStep().run(step_input)

    assert isinstance(res, ModelLibrary)
    with res:
        for i, model in enumerate(res):
            assert model.meta.cal_step.outlier_detection == "COMPLETE"
            res.shelve(model, i, modify=False)
