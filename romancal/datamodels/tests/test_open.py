import pytest

from astropy.io import fits
import asdf
import numpy as np
from numpy.testing import assert_array_equal

from roman_datamodels import datamodels
from roman_datamodels.tests import util

def test_asdf_file_input():
    tree = util.mk_level2_image()
    with asdf.AsdfFile() as af:
        af.tree = {'roman': tree}
        model = datamodels.open(af)
        assert model.meta.telescope == 'ROMAN'
        model.close()
        # should the asdf file be closed by model.close()?

def test_path_input(tmp_path):
    file_path = tmp_path/"test.asdf"
    with asdf.AsdfFile() as af:
        tree = util.mk_level2_image()
        af.tree = {'roman': tree}
        af.write_to(file_path)

    # Test with PurePath input:
    with datamodels.open(file_path) as model:
        assert model.meta.telescope == "ROMAN"
        af = model._asdf

    # When open creates the file pointer, it should be
    # closed when the model is closed:
    assert af._closed

    # Test with string input:
    with datamodels.open(str(file_path)) as model:
        assert model.meta.telescope == "ROMAN"
        af = model._asdf

    assert af._closed

    # Appropriate error when file is missing:
    with pytest.raises(FileNotFoundError):
        with datamodels.open(tmp_path/"missing.asdf"):
            pass


def test_model_input(tmp_path):
    file_path = tmp_path/"test.asdf"
    data = np.random.uniform(size=(1024, 1024)).astype(np.float32)
    with asdf.AsdfFile() as af:
        af.tree = {'roman': util.mk_level2_image()}
        af.tree['roman'].meta['bozo'] = 'clown'
        af.tree['roman'].data = data
        af.write_to(file_path)

    original_model = datamodels.open(file_path)
    reopened_model = datamodels.open(original_model)

    # It's essential that we get a new instance so that the original
    # model can be closed without impacting the new model.
    assert reopened_model is not original_model

    assert_array_equal(original_model.data, data)
    original_model.close()
    assert_array_equal(reopened_model.data, data)
    reopened_model.close()


def test_invalid_input():
    with pytest.raises(TypeError):
        datamodels.open(fits.HDUList())


def test_memmap(tmp_path):
    file_path = tmp_path/"test.asdf"
    with asdf.AsdfFile() as af:
        af.tree = {'roman': util.mk_level2_image()}
        af.tree['roman'].data = np.zeros((400, 400,), dtype=np.float32)
        af.write_to(file_path)

    with datamodels.open(file_path, memmap=True) as model:
        assert isinstance(model.data.base, np.memmap)

    with datamodels.open(file_path, memmap=False) as model:
        assert not isinstance(model.data.base, np.memmap)

    # Default should be false:
    with datamodels.open(file_path) as model:
        assert not isinstance(model.data.base, np.memmap)
