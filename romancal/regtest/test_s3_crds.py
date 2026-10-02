import crds.config
import crds.utils
import pytest

from romancal.pipeline.exposure_pipeline import ExposurePipeline
from romancal.pipeline.mosaic_pipeline import MosaicPipeline

pytestmark = pytest.mark.bigdata


@pytest.fixture(scope="module")
def tmp_cache(tmp_path_factory):
    yield tmp_path_factory.mktemp("crds")


@pytest.fixture
def s3_crds(tmp_cache, monkeypatch):
    old_state = crds.config.get_crds_state()
    crds.utils.clear_function_caches()
    # this reproduces configuration from the crds_s3_set script and from crds tests
    monkeypatch.setenv("CRDS_MODE", "s3")
    monkeypatch.setenv("CRDS_S3_BUCKET", "stpubdata")
    monkeypatch.setenv("CRDS_S3_PREFIX", "/roman/crds")
    monkeypatch.setenv("CRDS_S3_ENABLED", "1")
    monkeypatch.setenv("CRDS_S3_RETURN_URI", "0")
    monkeypatch.setenv(
        "CRDS_DOWNLOAD_PLUGIN",
        "crds_s3_get ${FILENAME} -d ${OUTPUT_PATH} -s ${FILE_SIZE} -c ${FILE_SHA1SUM}",
    )
    monkeypatch.setenv("CRDS_DOWNLOAD_MODE", "plugin")
    monkeypatch.setenv("CRDS_MAPPING_URI", "s3://stpubdata/roman/crds/mappings/roman")
    monkeypatch.setenv(
        "CRDS_REFERENCE_URI", "s3://stpubdata/roman/crds/references/roman"
    )
    monkeypatch.setenv("CRDS_CONFIG_URI", "s3://stpubdata/roman/crds/config/roman")
    monkeypatch.setenv("CRDS_REF_SUBDIR_MODE", "flat")
    monkeypatch.setenv("CRDS_SERVER_URL", "https://roman-crds-serverless.stsci.edu")
    monkeypatch.setenv("CRDS_OBSERVATORY", "roman")
    monkeypatch.setenv("CRDS_PATH", str(tmp_cache))

    yield
    crds.config.set_crds_state(old_state)
    crds.utils.clear_function_caches()


def test_s3_elp(s3_crds, rtdata):
    input_data = "r0000101001001001001_0001_wfi01_f158_uncal.asdf"
    rtdata.get_data(f"WFI/image/{input_data}")
    rtdata.input = input_data

    # Test Pipeline
    ExposurePipeline.call(rtdata.input)

    # no truth comparison here since the context may differ
    # TODO do some basic checks


def test_s3_mos(s3_crds, rtdata):
    rtdata.get_asn("WFI/image/L3_regtest_asn.json")
    MosaicPipeline.call(rtdata.input)

    # no truth comparison here since the context may differ
    # TODO do some basic checks
