import json
from pathlib import Path

import pytest
import roman_datamodels as rdm

from romancal.datamodels import ModelLibrary, filetype

DATA_DIRECTORY = Path(__file__).parent / "data"


@pytest.fixture
def model():
    model = rdm.datamodels.ImageModel.create_fake_data(shape=(20, 20))
    model.meta.filename = "test.asdf"
    return model


@pytest.fixture
def model_filename(model, tmp_path):
    fn = tmp_path / "test.asdf"
    model.save(fn)
    return fn


@pytest.fixture
def model_file_handle(model_filename):
    with open(model_filename, "rb") as fh:
        yield fh


@pytest.fixture
def library(model):
    return ModelLibrary([model])


@pytest.fixture
def association_filename(library, function_jail):
    return library._save(function_jail)


@pytest.fixture
def association_file_handle(association_filename):
    with open(association_filename) as fh:
        yield fh


@pytest.fixture
def association_dict(association_filename):
    with open(association_filename) as f:
        asn = json.load(f)
    return asn


@pytest.mark.parametrize(
    "init_fixture, expected",
    [
        ("model", "DataModel"),
        ("model_filename", "asdf"),
        ("model_file_handle", "asdf"),
        ("library", "ModelLibrary"),
        ("association_filename", "asn"),
        ("association_file_handle", "asn"),
        ("association_dict", "asn"),
    ],
)
def test_filetype(init_fixture, expected, request):
    assert filetype.check(request.getfixturevalue(init_fixture)) == expected


@pytest.mark.parametrize("init_value", [2, "missing.txt", "test"])
def test_failures(init_value, function_jail):
    with pytest.raises(ValueError):
        filetype.check(init_value)
