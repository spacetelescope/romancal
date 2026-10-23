"""Unit tests for the epsf module"""

import numpy as np
import pytest
import roman_datamodels as rdm
from astropy.table import Table

from romancal.source_catalog import _epsf


@pytest.mark.parametrize(
    "attr_test",
    [
        lambda inputs, expected: inputs["result"].converged,
        lambda inputs, expected: len(inputs["result"].excluded_star_indices) == 0,
    ],
    ids=["converted", "excluded stars"],
)
def test_build_epsf_attributes(build_epsf, attr_test):
    """Test various attributes of the results from pick_stars"""
    result, expected = build_epsf

    inputs = {"result": result}
    assert attr_test(inputs, expected)


@pytest.mark.parametrize(
    "attr_test",
    [
        lambda inputs, expected: inputs["model"].validate() is None,
    ],
    ids=["validation"],
)
def test_epsf_from_model_attributes(epsf_from_model, attr_test):
    """Test building an ePSF from a model."""
    result, expected = epsf_from_model

    inputs = {"model": result}
    assert attr_test(inputs, expected)


def test_make_epsf_reference(build_epsf):
    """Validate the epsf reference creation."""
    result, _ = build_epsf

    model = _epsf.make_epsf_reference(result, "wfi01", "f158")

    model.validate()


@pytest.mark.parametrize(
    "attr_test",
    [
        lambda inputs, expected: len(inputs["stars"]) == expected["len"],
        lambda inputs, expected: np.isclose(
            inputs["max_psferr"], expected["max_psferr"]
        ),
    ],
    ids=["n_stars", "max_psferr"],
)
def test_pick_stars_attributes(pick_stars, attr_test):
    """Test various attributes of the results from pick_stars"""
    stars, max_psferr, expected = pick_stars

    inputs = {"stars": stars, "max_psferr": max_psferr}
    assert attr_test(inputs, expected)


# ######################
# Fixtures and utilities
# ######################
@pytest.fixture(
    params=[
        ({}, {}),  # use all defaults
    ],
    ids=["default kwargs"],
)
def build_epsf(pick_stars, request):
    """Build the ePSF using various inputs."""
    stars, max_psferr, _ = pick_stars
    kwargs, expected = request.param

    result = _epsf.build_epsf(stars, max_psferr=max_psferr, **kwargs)

    return result, expected


@pytest.fixture(
    params=[
        ({}, {}),  # use all defaults
    ],
    ids=["default kwargs"],
)
def epsf_from_model(l2_model, request):
    """Create an ePSF reference from an L2 model"""
    l2_model, l2_cat = l2_model
    kwargs, expected = request.param

    result = _epsf.epsf_from_model(l2_model, l2_cat, **kwargs)

    return result, expected


@pytest.fixture(scope="module")
def l2_model(rtdata_module):
    """Retrieve an L2 model

    Returns
    -------
    model, catalog : roman_datamodels.datamodels.ImageModel, astropy.table.Table
        L2 data and associated source catalog.
    """
    rtdata = rtdata_module

    model_fname = "r0000101001001001001_0001_wfi01_f158_cal.asdf"
    cat_fname = "r0000101001001001001_0001_wfi01_f158_cat.parquet"
    rtdata.get_data(f"WFI/image/{model_fname}")
    rtdata.get_data(f"WFI/image/{cat_fname}")

    cat = Table.read(cat_fname)
    with rdm.open(model_fname) as model:
        yield model, cat


@pytest.fixture(
    params=[
        ({}, {"len": 21, "max_psferr": 0.052558005}),  # use all defaults
    ],
    ids=["default kwargs"],
)
def pick_stars(l2_model, request):
    """Generate list of stars based on various inputs."""
    model, catalog = l2_model
    pick_stars_args, expected = request.param

    stars, max_psferr = _epsf.pick_stars(model, catalog, **pick_stars_args)

    return stars, max_psferr, expected
