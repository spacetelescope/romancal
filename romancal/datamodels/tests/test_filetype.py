from pathlib import Path

import pytest
import roman_datamodels as rdm

from romancal.datamodels import ModelLibrary, filetype

DATA_DIRECTORY = Path(__file__).parent / "data"


def test_filetype():
    file_1 = filetype.check(DATA_DIRECTORY / "empty.json")
    file_2 = filetype.check(DATA_DIRECTORY / "example_schema.json")
    with open(DATA_DIRECTORY / "fake.json") as file_h:
        file_3 = filetype.check(file_h)
    file_4 = filetype.check(DATA_DIRECTORY / "empty.asdf")
    file_5 = filetype.check(DATA_DIRECTORY / "pluto.asdf")
    with open(DATA_DIRECTORY / "pluto.asdf", "rb") as file_h:
        file_6 = filetype.check(file_h)
    file_7 = filetype.check(DATA_DIRECTORY / "fake.asdf")
    with open(DATA_DIRECTORY / "fake.json") as file_h:
        file_8 = filetype.check(file_h)
    file_9 = filetype.check(str(DATA_DIRECTORY / "pluto.asdf"))
    image_node = rdm.maker_utils.mk_level2_image(shape=(20, 20))
    im1 = rdm.datamodels.ImageModel(image_node)
    file_11 = filetype.check(im1)
    model_library = ModelLibrary([im1])
    file_10 = filetype.check(model_library)

    assert file_1 == "asn"
    assert file_2 == "asn"
    assert file_3 == "asn"
    assert file_4 == "asdf"
    assert file_5 == "asdf"
    assert file_6 == "asdf"
    assert file_7 == "asdf"
    assert file_8 == "asn"
    assert file_9 == "asdf"
    assert file_10 == "ModelLibrary"
    assert file_11 == "DataModel"

    with pytest.raises(ValueError):
        filetype.check(DATA_DIRECTORY / "empty.txt")

    with pytest.raises(ValueError):
        filetype.check(2)

    with pytest.raises(ValueError):
        filetype.check("test")
