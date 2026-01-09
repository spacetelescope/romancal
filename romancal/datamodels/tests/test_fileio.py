import pytest
import roman_datamodels.datamodels as rdm

from romancal.datamodels.fileio import open_dataset
from romancal.datamodels.library import ModelLibrary

OLD_TAG = "asdf://stsci.edu/datamodels/roman/tags/wfi_image-1.4.0"
OLD_TAG_VERSION = OLD_TAG.split("-")[1]


@pytest.fixture()
def model():
    model = rdm.ImageModel.create_fake_data(tag=OLD_TAG)
    model.meta.filename = "test.asdf"
    return model


@pytest.fixture()
def library(model):
    return ModelLibrary([model])


@pytest.fixture()
def model_filename(model, tmp_path):
    fn = tmp_path / "test.asdf"
    model.save(fn)
    return fn


@pytest.fixture()
def library_filename(library, tmp_path):
    library._save(tmp_path)
    return tmp_path / "asn.json"


@pytest.fixture()
def list_of_models(model):
    return [model]


@pytest.fixture()
def list_of_filenames(model_filename):
    return [model_filename]


# @pytest.fixture(params=[
#     "model",
#     "library",
#     "model_filename",
#     "library_filename",
#     "list_of_models",
#     "list_of_filenames",
# ])
# def valid_dataset(request):
#     return request.getfixturevalue(request.param)


@pytest.mark.parametrize("return_type", [True, False])
@pytest.mark.parametrize("as_library", [True, False])
@pytest.mark.parametrize(
    "dataset_fixture, dataset_type, expected_return_type",
    [
        ("model", "DataModel", rdm.DataModel),
        ("library", "ModelLibrary", ModelLibrary),
        ("model_filename", "asdf", rdm.DataModel),
        ("library_filename", "asn", ModelLibrary),
        ("list_of_models", "unknown", ModelLibrary),
        ("list_of_filenames", "unknown", ModelLibrary),
    ],
)
def test_open_dataset(
    request,
    return_type,
    as_library,
    dataset_fixture,
    dataset_type,
    expected_return_type,
):
    result = open_dataset(
        request.getfixturevalue(dataset_fixture),
        return_type=return_type,
        as_library=as_library,
    )
    if return_type:
        return_value, detected_dataset_type = result
        assert detected_dataset_type == dataset_type
    else:
        return_value = result
    if as_library:
        expected_return_type = ModelLibrary
    assert isinstance(return_value, expected_return_type)


@pytest.mark.parametrize(
    "dataset",
    [
        "model_filename",
        "library_filename",
        "list_of_models",
        "list_of_filenames",
    ],
)
def test_open_kwargs(dataset, monkeypatch):
    class TestException(Exception):
        pass

    def patched_open(self, *args, **kwargs):
        assert "test" in kwargs
        raise TestException()

    monkeypatch.setattr(rdm, "open", patched_open)
    monkeypatch.setattr(ModelLibrary, "__init__", patched_open)

    with pytest.raises(TestException):
        open_dataset(dataset, open_kwargs={"test": 1})


@pytest.mark.parametrize("update_version", [True, False])
@pytest.mark.parametrize("as_library", [True, False])
@pytest.mark.parametrize(
    "dataset",
    [
        "model",
        "model_filename",
        "library_filename",
        "list_of_filenames",
    ],
)
def test_update_version(request, update_version, as_library, dataset):
    opened_dataset = open_dataset(
        request.getfixturevalue(dataset),
        update_version=update_version,
        as_library=as_library,
    )
    if isinstance(opened_dataset, rdm.DataModel):
        model = opened_dataset
    else:
        with opened_dataset:
            model = opened_dataset.borrow(0)
            opened_dataset.shelve(model, modify=False)

    tag_version = model.tag.split("-")[1]
    if update_version:
        assert tag_version > OLD_TAG_VERSION
    else:
        assert tag_version == OLD_TAG_VERSION
