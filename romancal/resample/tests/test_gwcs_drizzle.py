from unittest.mock import MagicMock

import numpy as np
import pytest

from ..gwcs_drizzle import GWCSDrizzle


@pytest.mark.parametrize(
    "product, outwcs, wt_scl, pixfrac, kernel, fillval, expected_exception",
    [
        (
            MagicMock(
                data=np.array([[1, 2], [3, 4]]),
                weight=np.array([[0.5, 0.5], [0.5, 0.5]]),
                meta=MagicMock(
                    wcs="WCS", resample=MagicMock(product_exposure_time=1.0)
                ),
                context=np.zeros((1, 2, 2), dtype=np.int32),
            ),
            None,
            None,
            1.0,
            "square",
            "INDEF",
            None,
        ),
        (
            MagicMock(
                data=np.array([[1, 2], [3, 4]]),
                weight=np.array([[1, 1], [1, 1]]),
                meta=MagicMock(
                    wcs="WCS", resample=MagicMock(product_exposure_time=2.0)
                ),
                context=np.zeros((1, 2, 2), dtype=np.int32),
            ),
            "CustomWCS",
            "exptime",
            0.5,
            "gaussian",
            "0",
            None,
        ),
        (
            MagicMock(
                data=np.array([[0, 1], [2, 3]]),
                weight=np.array([[1, 0], [0, 1]]),
                meta=MagicMock(
                    wcs="WCS", resample=MagicMock(product_exposure_time=0.5)
                ),
                context=np.zeros((1, 2, 2), dtype=np.int32),
            ),
            None,
            "expsq",
            0.8,
            "lanczos3",
            "NaN",
            None,
        ),
    ],
)
def test_gwcs_drizzle_init_happy_path(
    product, outwcs, wt_scl, pixfrac, kernel, fillval, expected_exception
):
    if expected_exception:
        with pytest.raises(expected_exception):
            GWCSDrizzle(product, outwcs, wt_scl, pixfrac, kernel, fillval)
    else:
        drizzle = GWCSDrizzle(product, outwcs, wt_scl, pixfrac, kernel, fillval)

    if not expected_exception:
        assert drizzle.outwcs == (outwcs or product.meta.wcs)
        assert drizzle.wt_scl == ("" if wt_scl is None else wt_scl)
        assert drizzle.pixfrac == pixfrac
        assert drizzle.kernel == kernel
        assert drizzle.fillval == fillval


@pytest.mark.parametrize(
    "pixfrac, kernel, expected_pixfrac, expected_kernel",
    [
        (0, "square", 0, "square"),
        (
            2.0,
            "nonexistent",
            2.0,
            "nonexistent",
        ),  # Invalid kernel, but no validation in __init__
    ],
)
def test_gwcs_drizzle_init_edge_cases(
    pixfrac, kernel, expected_pixfrac, expected_kernel
):
    # Arrange
    product = MagicMock(
        data=np.array([[1, 2], [3, 4]]),
        weight=np.array([[1, 1], [1, 1]]),
        meta=MagicMock(wcs="WCS", resample=MagicMock(product_exposure_time=1.0)),
        context=np.zeros((1, 2, 2), dtype=np.int32),
    )

    # Act
    drizzle = GWCSDrizzle(product, pixfrac=pixfrac, kernel=kernel)

    # Assert
    assert drizzle.pixfrac == expected_pixfrac
    assert drizzle.kernel == expected_kernel
