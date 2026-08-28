"""Test features of transforms"""

from pathlib import Path

import numpy as np
import pytest

from romancal.orientation import _transforms as tlib


@pytest.mark.parametrize(
    "reffile, expected",
    [
        (
            None,
            np.array(
                [
                    [2.77555756e-17, -5.00000000e-01, -8.66025404e-01],
                    [-8.65672497e-03, 8.65992954e-01, -4.99981265e-01],
                    [9.99962530e-01, 7.49694373e-03, -4.32836248e-03],
                ]
            ),
        ),
        (
            "bam_ref_good.asdf",
            np.array(
                [
                    [2.77555756e-17, -5.00000000e-01, -8.66025404e-01],
                    [-8.65672497e-03, 8.65992954e-01, -4.99981265e-01],
                    [9.99962530e-01, 7.49694373e-03, -4.32836248e-03],
                ]
            ),
        ),
        (
            "bam_ref_alternate.asdf",
            np.array(
                [
                    [0.0, -0.5, 0.8660254],
                    [0.00865672, -0.86599295, -0.49998126],
                    [0.99996253, 0.00749694, 0.00432836],
                ]
            ),
        ),
    ],
)
def test_bam_override(reffile, expected, caplog):
    """Test local reference file reading

    Parameters
    ----------
    reffile : str
        Name of the local reference file.

    expected : numpy.array(3, 3)
        B-to-FGS matrix
    """
    header = {
        "roman.meta.exposure.start_time": "2026-10-01",
        "roman.meta.instrument.name": "wfi",
    }
    if reffile is None:
        path = None
    else:
        path = Path(__file__).parent / "data" / reffile
    m = tlib.calc_m_b2fgs(path, header)

    assert np.allclose(m, expected)

    assert "Unable to retrieve BAM reference" not in caplog.text
