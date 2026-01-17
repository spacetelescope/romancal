import pytest
import roman_datamodels.datamodels as rdm

from romancal.datamodels.migration import update_model_version


@pytest.fixture
def latest_model():
    return rdm.ImageModel.create_fake_data()


@pytest.fixture
def old_model():
    return rdm.ImageModel.create_fake_data(
        tag="asdf://stsci.edu/datamodels/roman/tags/wfi_image-1.4.0"
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
    new_model = update_model_version(old_model)
    assert new_model is not old_model
    assert new_model.tag != old_model.tag
    assert new_model.tag == latest_model.tag
    assert new_model.meta.observation.wfi_parallel == wp_bool


def test_no_update(latest_model):
    new_model = update_model_version(latest_model)
    assert new_model is latest_model
