import numpy as np
import pytest
import roman_datamodels.datamodels as rdm

from romancal.associations.asn_from_list import asn_from_list
from romancal.datamodels.library import ModelLibrary
from romancal.pipeline import ExposurePipeline


@pytest.fixture(scope="function")
def input_value(request, tmp_path):
    model = rdm.RampModel.create_fake_data(shape=(2, 20, 20))
    match request.param:
        case "datamodel_fn":
            fn = tmp_path / "model.asdf"
            model.save(fn)
            return fn
        case "datamodel":
            return model
        case "asn_fn":
            model.meta.filename = "foo.asdf"
            model.save(tmp_path / model.meta.filename)
            asn = asn_from_list([model.meta.filename], product_name="foo_out")
            base_fn, contents = asn.dump(format="json")
            asn_filename = tmp_path / base_fn
            with open(asn_filename, "w") as f:
                f.write(contents)
            return asn_filename
        case "library":
            return ModelLibrary([model])
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
    # skip all steps
    [setattr(getattr(pipeline, k), "skip", True) for k in pipeline.step_defs]
    output_value = pipeline.run(input_value)
    assert isinstance(output_value, expected_output_type)


@pytest.mark.parametrize("save_results", [True, False])
def test_elp_save_results(function_jail, tmp_path, save_results):
    """
    Test that the elp respects save_results.
    """
    output_path = tmp_path / "output"
    output_path.mkdir()

    model = rdm.ScienceRawModel.create_fake_data()
    model.meta.filename = "test_uncal.asdf"
    model.meta.exposure.read_pattern = [[1], [2]]  # truncated for 2 groups below
    model.meta.exposure.ma_table_number = 5
    model.data = np.zeros((2, 4096, 4096), dtype=model.data.dtype)
    model.amp33 = np.zeros((2, 4096, 128), dtype=model.amp33.dtype)

    pipeline = ExposurePipeline()
    pipeline.output_dir = str(output_path)
    pipeline.save_results = save_results

    # skip tweakreg as it fails for an empty catalog and so
    # that the wfiwcs files are not produced
    pipeline.tweakreg.skip = True

    pipeline.run(model)
    # check that the current directory doesn't contain extra files
    assert set(p.name for p in tmp_path.iterdir()) == {"output"}

    output_files = set(p.name for p in output_path.iterdir())
    if save_results:
        assert output_files == {"test_segm.asdf", "test_cal.asdf", "test_cat.parquet"}
    else:
        # TODO shouldn't this be empty?
        assert output_files == {"test_segm.asdf", "test_cat.parquet"}
