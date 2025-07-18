"""Test module v1_calculate"""

from pathlib import Path

import numpy as np
import pytest
from astropy.table import Table
from astropy.time import Time
from roman_datamodels.datamodels import ScienceRawModel

import romancal.orientation.v1_calculate as v1c

DATA_PATH = Path(__file__).parent / "data"

# Engineering parameters
MAST_GOOD_STARTTIME = Time("2027-03-23T19:20:40", format="isot")
MAST_GOOD_ENDTIME = Time("2027-03-23T19:21:36", format="isot")


def test_from_models_mast(tmp_path):
    """Test v1_calculate_from_models for basic running"""
    model = ScienceRawModel.create_fake_data(
        {
            "meta": {
                "exposure": {
                    "start_time": MAST_GOOD_STARTTIME,
                    "end_time": MAST_GOOD_ENDTIME,
                },
            }
        }
    )

    try:
        v1_table = v1c.v1_calculate_from_models([model])
    except ValueError as exception:
        pytest.skip(
            f"MAST engineering database not available, possibly no token specified: {exception}"
        )
    v1_formatted = v1c.simplify_table(v1_table)

    # Save for post-test examination and read in to remove object infomation
    del v1_formatted["source"]
    v1_formatted.write(tmp_path / "test_from_models_mast.ecsv", format="ascii.ecsv")

    truth = Table.read(DATA_PATH / "test_from_models_mast.ecsv")
    errors = v1_compare_simplified_tables(v1_formatted, truth)
    errors_str = "\n".join(errors)
    assert len(errors) == 0, f"V1 tables are different: {errors_str}"


def test_over_time_mast(tmp_path):
    """Test v1_calculate_over_time for basic running"""
    try:
        v1_table = v1c.v1_calculate_over_time(MAST_GOOD_STARTTIME, MAST_GOOD_ENDTIME)
    except ValueError as exception:
        pytest.skip(
            f"MAST engineering database not available, possibly no token specified: {exception}"
        )
    v1_formatted = v1c.simplify_table(v1_table)

    # Save for post-test examination
    v1_formatted.write(tmp_path / "test_over_time_mast.ecsv", format="ascii.ecsv")

    truth = Table.read(DATA_PATH / "test_over_time_mast.ecsv")
    errors = v1_compare_simplified_tables(v1_formatted, truth)
    errors_str = "\n".join(errors)
    assert len(errors) == 0, f"V1 tables are different: {errors_str}"


# ######################
# Fixtures and utilities
# ######################
def v1_compare_simplified_tables(a, b, rtol=1e-05):
    """Calculate diff between tables generated by v1_calculate"""
    errors = []
    if len(a) != len(b):
        errors.append(f"len(a)={len(a)} not equal to len(b)={len(b)}")

    if not all(a["obstime"] == b["obstime"]):
        errors.append("obstimes different")

    if not np.allclose(a["ra"], b["ra"], rtol=rtol):
        errors.append("RAs are different")

    if not np.allclose(a["dec"], b["dec"], rtol=rtol):
        errors.append("DECs are different")

    if not np.allclose(a["pa_v3"], b["pa_v3"], rtol=rtol):
        errors.append("PA_V3s are different")

    return errors
