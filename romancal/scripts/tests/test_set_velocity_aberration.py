"""Test set_velocity_aberration"""

import subprocess

import roman_datamodels as rdm
from numpy import isclose

# Testing constants
GOOD_VELOCITY = (100.0, 100.0, 100.0)
GOOD_POS = (359.0, -2.0)
GOOD_SCALE_FACTOR = 1.000316017905845
GOOD_APPARENT_RA = 359.01945099823
GOOD_APPARENT_DEC = -1.980247580394956


def test_velocity_aberration_script(tmp_path):
    """Test the whole script on a FITS file"""
    path = tmp_path / "velocity_aberration_tmpfile.asdf"

    model = rdm.datamodels.ScienceRawModel.create_fake_data(
        {
            "meta": {
                "ephemeris": {
                    "velocity_x": GOOD_VELOCITY[0],
                    "velocity_y": GOOD_VELOCITY[1],
                    "velocity_z": GOOD_VELOCITY[2],
                },
                "wcsinfo": {"ra_ref": GOOD_POS[0], "dec_ref": GOOD_POS[1]},
            }
        }
    )
    model.save(path)

    subprocess.check_call(["roman_set_velocity_aberration", path])  # noqa: S603 S607

    with rdm.open(path) as model:
        va = model.meta.velocity_aberration
        assert isclose(va.ra_reference, GOOD_APPARENT_RA, rtol=0, atol=1e-7)
        assert isclose(va.dec_reference, GOOD_APPARENT_DEC, rtol=0, atol=1e-7)
        assert isclose(va.scale_factor, GOOD_SCALE_FACTOR, rtol=0, atol=1e-7)
