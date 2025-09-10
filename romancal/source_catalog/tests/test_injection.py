"""
Unit tests for the Roman source injection step code
"""

import pytest

pytest.importorskip("romanisim")

import numpy as np
from astropy import table
from astropy.time import Time
from roman_datamodels import stnode
from roman_datamodels.datamodels import ImageModel, MosaicModel
from romanisim import bandpass, parameters

from romancal.skycell.tests.test_skycell_match import mk_gwcs
from romancal.source_catalog.injection import inject_sources

# Set parameters
RA = 270.0
DEC = 66.0
ROLL = 0
SHAPE = (200, 200)
XPOS_IDX = [50, 50, 150, 150]
YPOS_IDX = [50, 150, 50, 150]
MAG_FLUX = 1e-9
DETECTOR = "WFI04"
SCA = 4
FILTER = "F158"
RNG_SEED = 42
MATABLE = 4

# Create gaussian noise generators
# sky should generate ~0.2 electron / s / pix.
# MJy / sr has similar magnitude to electron / s (i.e., within a factor
# of several), so just use 0.2 here.
MEANFLUX = 0.2


def make_test_data():
    # Create Four-quadrant pattern of gaussian noise, centered around one
    # Each quadrant's gaussian noise scales like total exposure time
    # (total files contributed to each quadrant)
    quarter_shape = (int(SHAPE[0] / 2), int(SHAPE[1] / 2))
    data = err = np.zeros(shape=SHAPE)
    noise_rng = np.random.default_rng(RNG_SEED)

    # Populate the data array with gaussian noise
    data[0 : quarter_shape[0], 0 : quarter_shape[1]] = noise_rng.normal(
        scale=(0.01 * MEANFLUX), size=quarter_shape
    )
    data[0 : quarter_shape[0], quarter_shape[1] : SHAPE[1]] = noise_rng.normal(
        scale=(0.02 * MEANFLUX), size=quarter_shape
    )
    data[quarter_shape[0] : SHAPE[0], 0 : quarter_shape[1]] = noise_rng.normal(
        scale=(0.05 * MEANFLUX), size=quarter_shape
    )
    data[quarter_shape[0] : SHAPE[0], quarter_shape[1] : SHAPE[1]] = noise_rng.normal(
        scale=(0.10 * MEANFLUX), size=quarter_shape
    )

    # Define Poisson Noise
    err[0 : quarter_shape[0], 0 : quarter_shape[1]] = 0.01 * MEANFLUX
    err[0 : quarter_shape[0], quarter_shape[1] : SHAPE[1]] = 0.02 * MEANFLUX
    err[quarter_shape[0] : SHAPE[0], 0 : quarter_shape[1]] = 0.05 * MEANFLUX
    err[quarter_shape[0] : SHAPE[0], quarter_shape[1] : SHAPE[1]] = 0.1 * MEANFLUX

    return data, err


@pytest.fixture
def image_model():
    defaults = {"meta": parameters.default_parameters_dictionary}
    defaults["meta"]["instrument"]["detector"] = DETECTOR
    defaults["meta"]["instrument"]["optical_element"] = FILTER
    defaults["meta"]["wcsinfo"]["ra_ref"] = RA
    defaults["meta"]["wcsinfo"]["dec_ref"] = DEC
    defaults["meta"]["wcsinfo"]["roll_ref"] = ROLL

    model = ImageModel.create_fake_data(defaults=defaults, shape=SHAPE)
    model.meta.filename = "none"
    model.meta.exposure.read_pattern = parameters.read_pattern[MATABLE]
    model.meta.cal_step = stnode.L2CalStep.create_fake_data()
    model.meta.cal_logs = stnode.CalLogs.create_fake_data()

    data, err = make_test_data()
    model.data = data
    model.err = err

    # Create WCS
    model.meta.wcs = mk_gwcs(
        model.meta.wcsinfo.ra_ref,
        model.meta.wcsinfo.dec_ref,
        model.meta.wcsinfo.roll_ref,
        bounding_box=((-0.5, SHAPE[0] - 0.5), (-0.5, SHAPE[1] - 0.5)),
        shape=SHAPE,
    )

    return model


@pytest.fixture
def mosaic_model():
    defaults = {
        "meta": {
            "coadd_info": {
                "time_first": Time("2024-01-01T12:00:00.000", format="isot"),
            },
            "instrument": {
                "optical_element": FILTER,
            },
            "resample": {"pixfrac": 1.0},
            "wcsinfo": {
                "pixel_scale": 1.5277777769528157e-05,
                "ra_ref": RA,
                "dec_ref": DEC,
                "roll_ref": ROLL,
            },  # Taken from regtest test L3 mosaic.
        }
    }
    model = MosaicModel.create_fake_data(defaults=defaults, shape=SHAPE)
    model.meta.filename = "none"
    model.meta.ref_file = stnode.RefFile.create_fake_data()
    model.meta.cal_step = stnode.L3CalStep.create_fake_data()
    model.cal_logs = stnode.CalLogs.create_fake_data()
    data, err = make_test_data()
    model.data = data
    model.err = err
    model.var_poisson = err**2
    model.weight = 1.0 / err

    # Create WCS
    model.meta.wcs = mk_gwcs(
        model.meta.wcsinfo.ra_ref,
        model.meta.wcsinfo.dec_ref,
        model.meta.wcsinfo.roll_ref,
        bounding_box=((-0.5, SHAPE[0] - 0.5), (-0.5, SHAPE[1] - 0.5)),
        shape=SHAPE,
    )

    return model


def make_catalog(metadata):
    # Create WCS
    wcsobj = mk_gwcs(
        metadata.wcsinfo.ra_ref,
        metadata.wcsinfo.dec_ref,
        metadata.wcsinfo.roll_ref,
        bounding_box=((-0.5, SHAPE[0] - 0.5), (-0.5, SHAPE[1] - 0.5)),
        shape=SHAPE,
    )

    # Create normalized psf source catalog (same source in each quadrant)
    ra, dec = wcsobj.pixel_to_world_values(np.array(XPOS_IDX), np.array(YPOS_IDX))

    tabcat = table.Table()
    tabcat["ra"] = ra
    tabcat["dec"] = dec
    tabcat[FILTER] = len(XPOS_IDX) * [MAG_FLUX]
    tabcat["type"] = len(XPOS_IDX) * ["PSF"]
    tabcat["n"] = len(XPOS_IDX) * [-1]
    tabcat["half_light_radius"] = len(XPOS_IDX) * [-1]
    tabcat["pa"] = len(XPOS_IDX) * [-1]
    tabcat["ba"] = len(XPOS_IDX) * [-1]

    return tabcat

# Ignore this warning - it comes from stpsf,
#   and they are aware of the need for an upgrade.
@pytest.mark.filterwarnings("ignore:Python 3.14 will, by default,"
" filter extracted tar archives and"
" reject files or modify their metadata."
" Use the filter argument to control this behavior.:DeprecationWarning")
def test_inject_sources(image_model, mosaic_model):
    for si_model in (image_model, mosaic_model):
        """Test simple source injection"""
        cat = make_catalog(si_model.meta)

        data_orig = si_model.copy()

        si_model = inject_sources(si_model, cat)

        # Ensure that sources were actually injected
        for x_val, y_val in zip(XPOS_IDX, YPOS_IDX, strict=False):
            assert np.all(
                si_model.data[y_val - 1 : y_val + 2, x_val - 1 : x_val + 2]
                > data_orig.data[y_val - 1 : y_val + 2, x_val - 1 : x_val + 2]
            )

        # Test that pixels far from the injected source are close to the original image
        # Numpy isclose is needed to determine equality, due to float precision issues
        assert np.all(
            np.isclose(
                si_model.data[90:110, 90:110],
                data_orig.data[90:110, 90:110],
                rtol=1e-06,
            )
        )

        # Ensure that every pixel's poisson variance has increased or
        # remained the same with the new sources injected
        # Numpy isclose is needed to determine equality,
        # due to float precision issues
        close_mask = np.isclose(si_model.var_poisson, data_orig.var_poisson, rtol=1e-06)
        assert False in close_mask
        assert np.all(
            si_model.var_poisson[~close_mask] > data_orig.var_poisson[~close_mask]
        )

        # Ensure that every data pixel value has increased or
        # remained the same with the new sources injected
        assert np.all(si_model.data[~close_mask] >= data_orig.data[~close_mask])

        # Ensure total added flux matches expected added flux
        # maggies to counts (large number)
        cps_conv = bandpass.get_abflux(FILTER, 2)
        # electrons to mjysr (roughly order unity in scale)
        unit_factor = bandpass.etomjysr(FILTER, 2)
        total_rec_flux = np.sum(si_model.data - data_orig.data)  # MJy / sr
        total_theo_flux = len(cat) * MAG_FLUX * cps_conv * unit_factor  # u.MJy / u.sr
        if isinstance(si_model, ImageModel):
            total_theo_flux /= parameters.reference_data["gain"].value
        assert np.isclose(total_rec_flux, total_theo_flux, rtol=4e-02)
