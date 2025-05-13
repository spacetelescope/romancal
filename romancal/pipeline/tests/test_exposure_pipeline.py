import pytest
import roman_datamodels.datamodels as rdm

from romancal.associations.asn_from_list import asn_from_list
from romancal.datamodels.library import ModelLibrary
from romancal.pipeline import ExposurePipeline


@pytest.fixture(scope="function")
def input_value(request, tmp_path):
    model = rdm.RampModel.fake_data(shape=(2, 20, 20))
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
def test_input_to_output(input_value, expected_output_type):
    """
    Test that for a particular input_value (as parametrized indirectly
    through the input_value fixtrue) the output is the expected type.
    """
    pipeline = ExposurePipeline()
    # don't fetch references
    pipeline.prefetch_references = False
    # skip all steps
    [setattr(getattr(pipeline, k), "skip", True) for k in pipeline.step_defs]
    output_value = pipeline.run(input_value)
    assert isinstance(output_value, expected_output_type)
