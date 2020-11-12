import pytest

from astropy.io import fits
import asdf
import numpy as np

from romancal import datamodels


def test_asdf_file_input():
    with asdf.AsdfFile() as af:
        af["meta"] = {"telescope": "binoculars"}

        model = datamodels.open(af)
        assert model.meta.telescope == "binoculars"

        # When an already open file is passed in, the model
        # is not responsible for closing it.
        model.close()
        assert not af._closed


def test_path_input(tmp_path):
    file_path = tmp_path/"test.asdf"
    with asdf.AsdfFile() as af:
        af["meta"] = {"telescope": "magnifying glass"}
        af.write_to(file_path)

    # Test with PurePath input:
    with datamodels.open(file_path) as model:
        assert model.meta.telescope == "magnifying glass"
        af = model._asdf

    # When open creates the file pointer, it should be
    # closed when the model is closed:
    assert af._closed

    # Test with string input:
    with datamodels.open(str(file_path)) as model:
        assert model.meta.telescope == "magnifying glass"
        af = model._asdf

    assert af._closed

    # Appropriate error when file is missing:
    with pytest.raises(FileNotFoundError):
        with datamodels.open(tmp_path/"missing.asdf"):
            pass


def test_invalid_input():
    with pytest.raises(TypeError):
        datamodels.open(fits.HDUList())


def test_memmap(tmp_path):
    file_path = tmp_path/"test.asdf"
    with asdf.AsdfFile() as af:
        af["data"] = np.zeros((1024,))
        af["meta"] = {}
        af.write_to(file_path)

    with datamodels.open(file_path, memmap=True) as model:
        assert isinstance(model.data.base, np.memmap)

    with datamodels.open(file_path, memmap=False) as model:
        assert not isinstance(model.data.base, np.memmap)

    # Default should be false:
    with datamodels.open(file_path) as model:
        assert not isinstance(model.data.base, np.memmap)
