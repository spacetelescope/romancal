"""Roman tests for time domain catalog creation"""

import pytest
from astropy.table import Table

from romancal.stpipe import RomanStep

# mark all tests in this module
pytestmark = [pytest.mark.bigdata, pytest.mark.soctests]


def test_time_domain_catalog(rtdata_module, resource_tracker, request):
    """Forced photometry of an L2 image at multiband catalog positions."""
    rtdata = rtdata_module
    inputfn = "r0000101001001001001_0001_wfi01_f158_cal.asdf"
    forcedfn = "r00001_p_v01001001001001_270p65x70y49_cat.parquet"
    outputfn = "r0000101001001001001_0001_wfi01_f158_fcat.parquet"
    rtdata.get_data(f"WFI/image/{forcedfn}")
    rtdata.get_data(f"WFI/image/{inputfn}")
    rtdata.input = inputfn
    rtdata.output = outputfn
    rtdata.get_truth(f"truth/WFI/image/{outputfn}")

    args = [
        "romancal.step.TimeDomainCatalogStep",
        inputfn,
        f"--forced_catalogs={forcedfn}",
    ]
    with resource_tracker.track(log=request):
        RomanStep.from_cmdline(args)

    cat = Table.read(outputfn)
    truth = Table.read(rtdata.truth)
    assert len(cat) > 0

    # Same columns with the same types and units as the truth file
    assert cat.colnames == truth.colnames
    for name in cat.colnames:
        assert cat[name].dtype == truth[name].dtype
        assert cat[name].unit == truth[name].unit
