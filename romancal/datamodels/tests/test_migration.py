import pytest
import roman_datamodels.datamodels as rdm

from romancal.datamodels.migration import MigrationWarning, update_model_version


@pytest.fixture
def latest_model():
    return rdm.ImageModel.create_fake_data()


@pytest.fixture
def latest_L3_model():
    return rdm.MosaicModel.create_fake_data()


@pytest.fixture
def old_model():
    yield rdm.ImageModel.create_fake_data(
        tag="asdf://stsci.edu/datamodels/roman/tags/wfi_image-1.4.0"
    )


@pytest.fixture
def old_L3_model():
    yield rdm.MosaicModel.create_fake_data(
        tag="asdf://stsci.edu/datamodels/roman/tags/wfi_mosaic-1.4.0"
    )


@pytest.mark.parametrize("close_on_update", [True, False])
def test_old_open_model(old_model, close_on_update, monkeypatch):
    close_called = False

    def close_watcher():
        nonlocal close_called
        close_called = True

    # check to see if model.close is called
    # patch the backing _asdf since we can't patch the DataModel
    monkeypatch.setattr(old_model._asdf, "close", close_watcher)
    with pytest.warns(MigrationWarning, match="hga_move"):
        update_model_version(old_model, close_on_update=close_on_update)
    assert close_on_update == close_called


@pytest.mark.parametrize(
    "vfs_value, wp_bool",
    [
        (1, False),
        (2, True),
    ],
)
def test_update(old_model, latest_model, vfs_value, wp_bool):
    old_model.meta.observation.visit_file_sequence = vfs_value
    with pytest.warns(MigrationWarning, match="hga_move"):
        new_model = update_model_version(old_model)
    assert new_model is not old_model
    assert new_model.tag != old_model.tag
    assert new_model.tag == latest_model.tag
    assert new_model.meta.observation.wfi_parallel == wp_bool
    assert not new_model.meta.exposure.hga_move


def test_update_adds_dq2(old_model):
    old_model.pop("dq2", None)
    with pytest.warns(MigrationWarning, match="empty dq2 array"):
        new_model = update_model_version(old_model)
    assert new_model.dq2.shape == old_model.dq.shape
    assert new_model.dq2.dtype == old_model.dq.dtype
    assert not new_model.dq2.any()
    assert "dq2" not in old_model


def test_update_adds_pixeldq2():
    old_ramp = rdm.RampModel.create_fake_data(
        tag="asdf://stsci.edu/datamodels/roman/tags/ramp-1.4.0"
    )
    old_ramp.pop("pixeldq2", None)
    with pytest.warns(MigrationWarning, match="empty pixeldq2 array"):
        new_ramp = update_model_version(old_ramp)
    assert new_ramp.pixeldq2.shape == old_ramp.pixeldq.shape
    assert new_ramp.pixeldq2.dtype == old_ramp.pixeldq.dtype
    assert not new_ramp.pixeldq2.any()


def test_update_keeps_existing_dq2(old_model):
    old_model["dq2"] = old_model.dq + 4
    with pytest.warns(MigrationWarning) as record:
        new_model = update_model_version(old_model)
    assert not any("dq2" in str(w.message) for w in record)
    assert (new_model.dq2 == 4).all()


def test_L3_update(old_L3_model, latest_L3_model):
    old_L3_model.meta.psf_match_reference_filter = "f158"
    new_L3_model = update_model_version(old_L3_model)
    assert new_L3_model is not old_L3_model
    assert new_L3_model.tag != old_L3_model.tag
    assert new_L3_model.tag == latest_L3_model.tag
    assert new_L3_model.meta.psf_match_reference_filter == "F158"


def test_no_update(latest_model):
    new_model = update_model_version(latest_model)
    assert new_model is latest_model
