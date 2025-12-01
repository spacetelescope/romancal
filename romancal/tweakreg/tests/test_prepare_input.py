import pytest

from romancal.datamodels import ModelLibrary
from romancal.tweakreg.tweakreg_step import TweakRegStep


@pytest.fixture()
def datamodel(tweakreg_image):
    model = tweakreg_image()
    model.meta.filename = "test.asdf"
    return model


@pytest.fixture()
def datamodel_filename(datamodel, tmp_path):
    fn = tmp_path / "test.asdf"
    datamodel.save(fn)
    return fn


@pytest.fixture()
def list_of_models(datamodel):
    return [datamodel]


@pytest.fixture()
def list_of_filenames(datamodel_filename):
    return [datamodel_filename]


@pytest.fixture()
def library(list_of_models):
    return ModelLibrary(list_of_models)


@pytest.fixture()
def association(library, tmp_path):
    fn = tmp_path / "asn.json"
    library._save(fn.parent)
    return fn


@pytest.mark.parametrize(
    "init",
    [
        "datamodel",
        "datamodel_filename",
        "list_of_models",
        "list_of_filenames",
        "library",
        "association",
    ],
)
def test_tweakreg_handle_input(init, request):
    assert isinstance(
        TweakRegStep()._prepare_input(request.getfixturevalue(init)), ModelLibrary
    )
