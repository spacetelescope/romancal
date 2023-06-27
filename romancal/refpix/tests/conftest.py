from enum import IntEnum

import numpy as np
import pytest
from astropy import units as u
from roman_datamodels.datamodels import RampModel, RefpixRefModel
from roman_datamodels.maker_utils import mk_ramp, mk_refpix

from romancal.refpix.data import Coefficients, Const, StandardView

RNG = np.random.default_rng(42)


class Dims(IntEnum):
    """
    Assumed data dimensions for dat (aid clarity of testing).
        - N_FRAMES: the number of frames used in the test data (arbitrary) choice
        - N_ROWS: the number of rows of the test data (the correction is
                  independent of this)
        - N_COLS: the number of columns of the test data (the correction is very
                  dependent on this)
    """

    N_FRAMES = 8
    N_ROWS = 12
    N_COLS = Const.N_COLUMNS + Const.CHAN_WIDTH
    N_CHAN = N_COLS // Const.CHAN_WIDTH
    PADDED_WIDTH = Const.CHAN_WIDTH + Const.PAD

    N_FFT_COLS = (N_ROWS * PADDED_WIDTH) // 2 + 1


@pytest.fixture(scope="module")
def data():
    return RNG.uniform(1, 100, size=(Dims.N_FRAMES, Dims.N_ROWS, Dims.N_COLS)).astype(
        np.float32
    )


@pytest.fixture(scope="module")
def datamodel(data):
    datamodel = mk_ramp(shape=(2, 2, 2))

    detector = data[:, :, : Const.N_COLUMNS]
    amp33 = data[:, :, Const.N_COLUMNS :]

    # Sanity check, that this is in line with the datamodel maker_utils
    assert detector.shape == (Dims.N_FRAMES, Dims.N_ROWS, Const.N_COLUMNS)
    assert amp33.shape == (Dims.N_FRAMES, Dims.N_ROWS, Const.CHAN_WIDTH)

    datamodel.data = u.Quantity(
        detector, unit=datamodel.data.unit, dtype=datamodel.data.dtype
    )
    datamodel.amp33 = u.Quantity(
        amp33, unit=datamodel.amp33.unit, dtype=datamodel.amp33.dtype
    )

    datamodel.border_ref_pix_left = datamodel.data[:, :, : Const.REF]
    datamodel.border_ref_pix_right = datamodel.data[:, :, -Const.REF :]

    return RampModel(datamodel)


@pytest.fixture(scope="module")
def standard(data):
    return StandardView(data)


@pytest.fixture(scope="module")
def channels(standard):
    return standard.channels


@pytest.fixture(scope="module")
def coeffs() -> Coefficients:
    coeffs = RNG.uniform(1, 100, size=(6, Const.N_DETECT_CHAN, Dims.N_FFT_COLS))

    coeffs = coeffs[:3, :, :] + 1.0j * coeffs[3:, :, :]

    return Coefficients(coeffs[0], coeffs[1], coeffs[2])


@pytest.fixture(scope="module")
def offset() -> np.ndarray:
    return RNG.uniform(0, 100, size=(Dims.N_ROWS, Dims.N_COLS)).astype(np.float32)


@pytest.fixture(scope="module")
def ref_pix_ref(coeffs):
    refpix_ref = mk_refpix(gamma=coeffs.gamma, zeta=coeffs.zeta, alpha=coeffs.alpha)

    return RefpixRefModel(refpix_ref)
