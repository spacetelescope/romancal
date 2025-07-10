"""Test module pointing_summary"""

import sys
from pathlib import Path

import pytest
from astropy.table import Table
from astropy.time import Time
from astropy.utils.diff import report_diff_values
from roman_datamodels.datamodels import ScienceRawModel

import romancal.orientation.pointing_summary as ps
from romancal.lib.engdb import engdb_tools

DATA_PATH = Path(__file__).parent / "data"

# Engineering parameters
GOOD_STARTTIME = Time("2027-02-23T01:00:00", format="isot")
GOOD_ENDTIME = Time("2027-02-23T01:00:05", format="isot")

# Header defaults
TARG_RA = 270.0
TARG_DEC = 66.01


@pytest.fixture
def engdb():
    """Setup the service to operate through the mock service"""
    try:
        engdb = engdb_tools.engdb_service()
    except RuntimeError as exception:
        pytest.skip(f"Engineering database unvailable: {exception}.")
    yield engdb


@pytest.fixture(scope="module")
def model_path(module_jail):
    """Create data file with needed header parameters"""
    model = ScienceRawModel.create_fake_data(
        {
            "meta": {
                "exposure": {"start_time": GOOD_STARTTIME, "end_time": GOOD_ENDTIME},
                "pointing": {"target_ra": TARG_RA, "target_dec": TARG_DEC},
            }
        }
    )

    model.meta.pointing.target_ra = TARG_RA
    model.meta.pointing.target_dec = TARG_DEC
    model.meta.pointing.ra_v1 = 91.08142005
    model.meta.pointing.dec_v1 = -66.60547869
    model.meta.wcsinfo.ra_ref = 90.70377653
    model.meta.wcsinfo.dec_ref = -66.59540224

    model_path = module_jail / "model.asdf"
    model.save(model_path)
    return model_path


def test_calc_pointing_deltas(engdb, model_path):
    """Test `calc_pointing_deltas` basic running"""
    truth = (
        "Delta(target=<SkyCoord (ICRS): (ra, dec) in deg"
        "\n    (270., 66.01)>, v1=<SkyCoord (ICRS): (ra, dec) in deg"
        "\n    (91.08142005, -66.60547869)>, refpoint=<SkyCoord (ICRS): (ra, dec) in deg"
        "\n    (90.70377653, -66.59540224)>, delta_v1=<Angle 179.26285175 deg>, delta_refpoint=<Angle 179.34985532 deg>)"
    )
    with ScienceRawModel(str(model_path)) as model:
        deltas = ps.calc_pointing_deltas(model)

    assert truth == str(deltas)


def test_calc_deltas(engdb, model_path, tmp_path):
    """Test `calc_deltas` basic running"""
    with ScienceRawModel(model_path) as model:
        deltas = ps.calc_deltas([model])

    # Save the calculated deltas for test debugging
    deltas["exposure"] = [m.meta.filename for m in deltas["exposure"]]
    deltas.write(tmp_path / "calc_deltas_truth.ecsv", format="ascii.ecsv")

    truth = Table.read(DATA_PATH / "calc_deltas_truth.ecsv")

    # round the delta values to a reasonable level
    deltas[0][4] = round(deltas[0][4], 8)
    deltas[0][5] = round(deltas[0][5], 8)
    truth[0][4] = round(truth[0][4], 8)
    truth[0][5] = round(truth[0][5], 8)

    assert report_diff_values(truth, deltas, fileobj=sys.stderr)
