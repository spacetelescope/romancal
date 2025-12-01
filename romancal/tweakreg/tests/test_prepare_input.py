import pytest

from romancal.datamodels import ModelLibrary
from romancal.tweakreg.tweakreg_step import TweakRegStep


@pytest.mark.parametrize(
    "init_type",
    [
        "datamodel",
        "datamodel_filename",
        "list_of_models",
        "list_of_filenames",
        "library",
        "asn",
    ],
)
def test_tweakreg_handle_input(base_image, tmp_path, init_type):
    match init_type:
        case "datamodel":
            init = base_image()
        case "datamodel_filename":
            init = tmp_path / "test.asdf"
            base_image().save(init)
        case "list_of_models":
            init = [base_image()]
        case "list_of_filenames":
            init = [tmp_path / "test.asdf"]
            base_image().save(init[0])
        case "library":
            init = ModelLibrary([base_image()])
        case "asn":
            init = tmp_path / "asn.json"
            model = base_image()
            model.meta.filename = "test.asdf"
            lib = ModelLibrary([model])
            lib._save(init.parent)
    assert isinstance(TweakRegStep()._prepare_input(init), ModelLibrary)
