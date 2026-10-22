"""Unit tests for the epsf module"""

import pytest

from astropy.table import Table
import numpy as np
import roman_datamodels as rdm

from romancal.source_catalog import _epsf as epsf


@pytest.mark.parametrize('pick_star_args, expected',
                         [
                             ({}, {'len': 21, 'max_psferr': 0.052558005}),  # use all defaults
                         ])
def test_pick_stars(l2_model, pick_star_args, expected):
    """Test functionlity of pick_stars"""
    model, catalog = l2_model

    stars, max_psferr = epsf.pick_stars(model, catalog, **pick_star_args)

    assert len(stars) == expected['len']
    assert np.isclose(expected['max_psferr'], max_psferr)


@pytest.fixture(scope='module')
def l2_model(rtdata_module):
    """Retrieve an L2 model

    Returns
    -------
    model, catalog : roman_datamodels.datamodels.ImageModel, astropy.table.Table
        L2 data and associated source catalog.
    """
    rtdata = rtdata_module

    model_fname = "r0000101001001001001_0001_wfi01_f158_cal.asdf"
    cat_fname = 'r0000101001001001001_0001_wfi01_f158_cat.parquet'
    rtdata.get_data(f"WFI/image/{model_fname}")
    rtdata.get_data(f"WFI/image/{cat_fname}")

    cat = Table.read(cat_fname)
    with rdm.open(model_fname) as model:
        yield model, cat
