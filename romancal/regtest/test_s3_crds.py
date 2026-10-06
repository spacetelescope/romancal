import crds.config
import crds.log
import crds.utils
import pytest
from roman_datamodels import datamodels as rdm

from romancal.pipeline.exposure_pipeline import ExposurePipeline

pytestmark = [pytest.mark.bigdata]
pytest.importorskip("awscli", reason="crds[aws] dependencies are needed for s3 tests")


@pytest.fixture(scope="module")
def tmp_cache(tmp_path_factory):
    yield tmp_path_factory.mktemp("crds")


@pytest.fixture
def s3_crds(tmp_cache, monkeypatch):
    # patch the download plugin?
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
    old_level = crds.log.set_verbose()
    yield
    crds.log.set_verbose(old_level)
    crds.config.set_crds_state(old_state)
    crds.utils.clear_function_caches()


def test_get_reference(s3_crds, caplog):
    pipeline = ExposurePipeline()
    model = rdm.ImageModel.create_fake_data()
    try:
        pipeline.get_reference_file(model, "flat")
    except crds.CrdsLookupError:
        # ok to have a lookup error as fetching config/mapping is enough
        pass
    assert [
        m
        for m in caplog.messages
        if "Loading config from URI 's3://stpubdata/roman/crds/config/roman/server_config'."
        in m
    ]
