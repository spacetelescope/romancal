from enum import IntEnum

import numpy as np
import pytest
from astropy import units as u
from roman_datamodels.maker_utils import mk_ramp

from romancal.refpix.data import Aligned, Standard, Width


class Dims(IntEnum):
    N_FRAMES = 8
    N_ROWS = 4096
    N_COLS = N_ROWS + Width.CHANNEL
    N_CHAN = N_COLS // Width.CHANNEL
    N_DETECT_CHAN = N_CHAN - 1
    N_ALIGN_COLS = Width.CHANNEL + Width.PAD


RNG = np.random.default_rng(42)


@pytest.fixture(scope="module")
def ramp_model():

    datamodel = mk_ramp()
    assert datamodel.data.shape == (Dims.N_FRAMES, Dims.N_ROWS, Dims.N_ROWS)
    assert datamodel.amp33.shape == (Dims.N_FRAMES, Dims.N_ROWS, Width.CHANNEL)

    data = RNG.uniform(1, 100, size=(Dims.N_FRAMES, Dims.N_ROWS, Dims.N_ROWS)).astype(
        datamodel.data.dtype
    )
    amp33 = RNG.uniform(
        1, 100, size=(Dims.N_FRAMES, Dims.N_ROWS, Width.CHANNEL)
    ).astype(datamodel.amp33.dtype)

    datamodel.data = u.Quantity(
        data, unit=datamodel.data.unit, dtype=datamodel.data.dtype
    )
    datamodel.border_ref_pix_left = datamodel.data[:, :, : Width.REF]
    datamodel.border_ref_pix_right = datamodel.data[:, :, -Width.REF :]

    datamodel.amp33 = u.Quantity(
        amp33, unit=datamodel.amp33.unit, dtype=datamodel.amp33.dtype
    )

    return datamodel


@pytest.fixture(scope="module")
def standard_data(ramp_model) -> Standard:
    return Standard.from_datamodel(ramp_model)


@pytest.fixture(scope="module")
def aligned_data(standard_data) -> Aligned:
    return Aligned.from_standard(standard_data)
