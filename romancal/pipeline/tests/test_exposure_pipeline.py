import numpy as np
import pytest
import roman_datamodels.datamodels as rdm
from astropy.time import Time

from romancal.associations.asn_from_list import asn_from_list
from romancal.datamodels.library import ModelLibrary
from romancal.pipeline import ExposurePipeline


@pytest.fixture
def fake_science_raw():
    model = rdm.ScienceRawModel.create_fake_data()
    model.meta.filename = "test_uncal.asdf"
    model.meta.exposure.read_pattern = [
        [1],
        [2],
        [3],
        [4],
    ]  # truncated for 4 groups below
    model.meta.exposure.ma_table_number = 5
    model.meta.exposure.start_time = Time(
        "2024-01-03T00:00:00.0", format="isot", scale="utc"
    )
    model.data = np.zeros((4, 4096, 4096), dtype=model.data.dtype)
    model.amp33 = np.zeros((4, 4096, 128), dtype=model.amp33.dtype)
    return model


@pytest.fixture(scope="function")
def input_value(request, tmp_path, fake_science_raw):
    match request.param:
        case "datamodel_fn":
            fn = tmp_path / "test_uncal.asdf"
            fake_science_raw.save(fn)
            return fn
        case "datamodel":
            return fake_science_raw
        case "asn_fn":
            fake_science_raw.meta.filename = "test_uncal.asdf"
            fake_science_raw.save(tmp_path / fake_science_raw.meta.filename)
            asn = asn_from_list(
                [fake_science_raw.meta.filename], product_name="foo_out"
            )
            base_fn, contents = asn.dump()
            asn_filename = tmp_path / base_fn
            with open(asn_filename, "w") as f:
                f.write(contents)
            return asn_filename
        case "library":
            return ModelLibrary([fake_science_raw])
        case value:
            raise Exception(f"Invalid parametrization: {value}")


@pytest.mark.parametrize(
    "input_value, expected_output_type",
    [
        ("datamodel_fn", rdm.DataModel),
        ("datamodel", rdm.DataModel),
        ("asn_fn", ModelLibrary),
        ("library", ModelLibrary),
    ],
    indirect=["input_value"],
)
def test_input_to_output(function_jail, input_value, expected_output_type):
    """
    Test that for a particular input_value (as parametrized indirectly
    through the input_value fixture) the output is the expected type.
    """
    pipeline = ExposurePipeline()
    # don't fetch references
    pipeline.prefetch_references = False
    # skip all steps except dq_init to allow the input ScienceRawModel to be converted to a RampModel
    [
        setattr(getattr(pipeline, k), "skip", True)
        for k in pipeline.step_defs
        if k != "dq_init"
    ]
    output_value, _, _ = pipeline.run(input_value)
    assert isinstance(output_value, expected_output_type)


@pytest.mark.parametrize(
    "input_value", ["datamodel", "datamodel_fn", "asn_fn", "library"], indirect=True
)
@pytest.mark.parametrize("save_results", [True, False])
def test_elp_save_results(function_jail, input_value, save_results, monkeypatch):
    """
    Test that the elp respects save_results.
    """
    output_path = function_jail / "output"
    output_path.mkdir()

    pipeline = ExposurePipeline()
    pipeline.output_dir = str(output_path)
    pipeline.save_results = save_results

    # don't try to actually run tweakreg as it will fail for an empty model
    monkeypatch.setattr(pipeline.tweakreg, "run", lambda init: init)
    allowed_files = set(p.name for p in function_jail.iterdir())

    pipeline.run(input_value)
    # check that the current directory doesn't contain extra files
    assert len(set(p.name for p in function_jail.iterdir()) - allowed_files) == 0

    output_files = set(p.name for p in output_path.iterdir())
    if save_results:
        assert output_files == {
            "test_cal.asdf",
            "test_wcs.asdf",
            "test_cat.parquet",
            "test_segm.asdf",
        }
    else:
        assert not output_files


@pytest.mark.parametrize("on_disk", [True, False])
def test_on_disk(function_jail, fake_science_raw, on_disk):
    fake_science_raw.meta.filename = "foo.asdf"
    fake_science_raw.save(fake_science_raw.meta.filename)
    asn = asn_from_list([fake_science_raw.meta.filename], product_name="foo_out")
    base_fn, contents = asn.dump()
    asn_filename = base_fn
    with open(asn_filename, "w") as f:
        f.write(contents)

    pipeline = ExposurePipeline()
    pipeline.on_disk = on_disk
    # don't fetch references
    pipeline.prefetch_references = False
    # skip all steps
    [setattr(getattr(pipeline, k), "skip", True) for k in pipeline.step_defs]
    pipeline.dq_init.skip = False  # unskip dqinit
    output_value, _, _ = pipeline.run(asn_filename)
    assert output_value._on_disk == on_disk
