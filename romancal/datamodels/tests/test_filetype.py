from pathlib import Path

import pytest

from romancal.datamodels import filetype

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

    assert file_1 == "asn"
    assert file_2 == "asn"
    assert file_3 == "asn"
    assert file_4 == "asdf"
    assert file_5 == "asdf"
    assert file_6 == "asdf"
    assert file_7 == "asdf"
    assert file_8 == "asn"
    assert file_9 == "asdf"

    with pytest.raises(ValueError):
        filetype.check(DATA_DIRECTORY / "empty.txt")

    with pytest.raises(ValueError):
        filetype.check(2)

    with pytest.raises(ValueError):
        filetype.check("test")
