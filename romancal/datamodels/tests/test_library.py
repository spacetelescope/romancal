import pytest
import roman_datamodels.datamodels as dm

from romancal.associations.asn_from_list import asn_from_list
from romancal.datamodels.library import ModelLibrary

# for the example association, set 2 different observation numbers
# so the association will have 2 groups (since all other group_id
# determining meta is the same, see `example_asn_path`)
_OBSERVATION_NUMBERS = [1, 1, 2]
_N_MODELS = len(_OBSERVATION_NUMBERS)
_N_GROUPS = len(set(_OBSERVATION_NUMBERS))
_PRODUCT_NAME = "foo_out"


@pytest.fixture
def example_asn_path(tmp_path):
    """
    Fixture that creates a simple association, saves it (and the models)
    to disk, and returns the path of the saved association
    """
    fns = []
    for i in range(_N_MODELS):
        m = dm.ImageModel.create_fake_data(shape=(2, 2))
        m.meta.observation.program = 1
        m.meta.observation.observation = _OBSERVATION_NUMBERS[i]
        m.meta.observation.visit = 1
        m.meta.observation.visit_file_group = 1
        m.meta.observation.visit_file_sequence = 1
        m.meta.observation.visit_file_activity = "01"
        m.meta.observation.exposure = 1
        base_fn = f"{i}.asdf"
        m.meta.filename = base_fn
        m.save(str(tmp_path / base_fn))
        fns.append(base_fn)
    asn = asn_from_list(fns, product_name=_PRODUCT_NAME)
    base_fn, contents = asn.dump(format="json")
    asn_filename = tmp_path / base_fn
    with open(asn_filename, "w") as f:
        f.write(contents)
    return asn_filename


@pytest.fixture
def example_library(example_asn_path):
    """
    Fixture that builds off of `example_asn_path` and returns a
    library created from the association with default options
    """
    return ModelLibrary(example_asn_path)


@pytest.mark.parametrize("attr", ["group_names", "group_indices"])
def test_group_with_no_datamodels_open(example_asn_path, attr, monkeypatch):
    """
    Test that the "grouping" methods do not call datamodels.open
    """

    # patch datamodels.open to always raise an exception
    # this will serve as a smoke test to see if any of the attribute
    # accesses (or instance creation) attempts to open models
    def no_open(*args, **kwargs):
        raise Exception()

    monkeypatch.setattr(dm, "open", no_open)

    # use example_asn_path here to make the instance after we've patched
    # datamodels.open
    library = ModelLibrary(example_asn_path)
    getattr(library, attr)


def test_asn_data(example_library):
    """
    Test that `asn` returns the association information
    """
    assert example_library.asn["products"][0]["name"] == _PRODUCT_NAME
