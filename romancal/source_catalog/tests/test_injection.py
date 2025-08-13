"""
Unit tests for the Roman source injection step code
"""

import numpy as np
import pytest
import galsim
from astropy.time import Time
from astropy import table
from roman_datamodels import stnode
from roman_datamodels.datamodels import ImageModel, MosaicModel
from romanisim import parameters, wcs, bandpass
from romancal.source_catalog.injection import inject_sources

# Set parameters
RA = 270.0
DEC = 66.0
SHAPE = (200, 200)
XPOS_IDX = [50, 50, 150, 150]
YPOS_IDX = [50, 150, 50, 150]
MAG_FLUX = 1e-9
DETECTOR = "WFI04"
SCA = 4
FILTER = "F158"
RNG_SEED = 42
MATABLE = 4


def make_catalog(metadata):
    # Create WCS
    twcs = wcs.GWCS(wcs.get_mosaic_wcs(metadata, shape=SHAPE))

    # Create normalized psf source catalog (same source in each quadrant)
    ra, dec = twcs._radec(np.array(XPOS_IDX), np.array(YPOS_IDX))

    tabcat = table.Table()
    tabcat['ra'] = np.degrees(ra)
    tabcat['dec'] = np.degrees(dec)
    tabcat[FILTER] = 4 * [MAG_FLUX]
    tabcat['type'] = 4 * ['PSF']
    tabcat['n'] = 4 * [-1]
    tabcat['half_light_radius'] = 4 * [-1]
    tabcat['pa'] = 4 * [-1]
    tabcat['ba'] = 4 * [-1]

    return tabcat


def make_test_image():
    rng = np.random.default_rng(seed=123)
    noise_scale = 2.5
    noise = rng.normal(0, noise_scale, size=SHAPE)
    data = noise
    err = np.zeros_like(data) + noise_scale

    return data, err


def make_test_mosaic():
    # Create Four-quadrant pattern of gaussian noise, centered around one
    # Each quadrant's gaussian noise scales like total exposure time
    # (total files contributed to each quadrant)
    data = err = np.zeros(shape=SHAPE)

    # Create gaussian noise generators
    # sky should generate ~0.2 electron / s / pix.
    # MJy / sr has similar magnitude to electron / s (i.e., within a factor
    # of several), so just use 0.2 here.
    meanflux = 0.2
    g1 = galsim.GaussianDeviate(RNG_SEED, mean=meanflux, sigma=0.01 * meanflux)
    g2 = galsim.GaussianDeviate(RNG_SEED, mean=meanflux, sigma=0.02 * meanflux)
    g3 = galsim.GaussianDeviate(RNG_SEED, mean=meanflux, sigma=0.05 * meanflux)
    g4 = galsim.GaussianDeviate(RNG_SEED, mean=meanflux, sigma=0.10 * meanflux)

    # Populate the mosaic data array with gaussian noise from generators
    g1.generate(data[0:100, 0:100])
    g2.generate(data[0:100, 100:200])
    g3.generate(data[100:200, 0:100])
    g4.generate(data[100:200, 100:200])

    # Define Poisson Noise of mosaic
    err[0:100, 0:100] = 0.01 * meanflux
    err[0:100, 100:200] = 0.02 * meanflux
    err[100:200, 0:100] = 0.05 * meanflux
    err[100:200, 100:200] = 0.1 * meanflux

    return data, err


@pytest.fixture
def image_model():
    defaults = {
        "meta" : parameters.default_parameters_dictionary
    }
    defaults["meta"]['instrument']['detector'] = DETECTOR
    defaults["meta"]['instrument']['optical_element'] = FILTER
    defaults["meta"]['wcsinfo']['ra_ref'] = RA
    defaults["meta"]['wcsinfo']['dec_ref'] = DEC

    model = ImageModel.create_fake_data(defaults=defaults, shape=SHAPE)
    model.meta.filename = "none"
    model.meta.exposure.read_pattern = parameters.read_pattern[MATABLE]
    model.meta.cal_step = stnode.L2CalStep.create_fake_data()
    model.meta.cal_logs = stnode.CalLogs.create_fake_data()
    data, err = make_test_image()
    model.data = data
    model.err = err

    # Create WCS
    twcs = wcs.GWCS(wcs.get_mosaic_wcs(model.meta, shape=SHAPE))
    model.meta.wcs = twcs.wcs

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
                'ra_ref': RA,
                'dec_ref': DEC,
            },  # Taken from regtest test L3 mosaic.
        }
    }
    model = MosaicModel.create_fake_data(defaults=defaults, shape=SHAPE)
    model.meta.filename = "none"
    model.meta.ref_file = stnode.RefFile.create_fake_data()
    model.meta.cal_step = stnode.L3CalStep.create_fake_data()
    model.cal_logs = stnode.CalLogs.create_fake_data()
    data, err = make_test_mosaic()
    model.data = data
    model.err = err
    model.var_poisson = err**2
    model.weight = 1.0 / err

    # Create WCS
    twcs = wcs.GWCS(wcs.get_mosaic_wcs(model.meta, shape=SHAPE))

    model.meta.wcs = twcs.wcs

    return model


def test_inject_sources(image_model, mosaic_model):
    for si_model in (image_model, mosaic_model):
        """Test simple source injection"""
        cat = make_catalog(si_model.meta)

        data_orig = si_model.copy()

        si_model = inject_sources(si_model, cat)

        # Ensure that sources were actually injected
        for x_val, y_val in zip(XPOS_IDX, YPOS_IDX):
            assert np.all(si_model.data[y_val - 1:y_val + 2, x_val - 1:x_val + 2] !=
                        data_orig.data[y_val - 1:y_val + 2, x_val - 1: x_val + 2])

        # Test that pixels far from the injected source are close to the original image
        # Numpy isclose is needed to determine equality, due to float precision issues
        assert np.all(np.isclose(si_model.data[90:110, 90:110],
            data_orig.data[90:110, 90:110], rtol=1e-06))

        # IMAGE L2 TEST
        if isinstance(si_model, ImageModel):
            # Test that the amount of added flux makes sense
            fluxeps = MAG_FLUX * bandpass.get_abflux(FILTER,
                int(si_model.meta['instrument']['detector'][3:]))  # u.electron / u.s
            assert np.abs((np.sum(si_model.data - data_orig.data) *
                parameters.reference_data['gain'].value / (4 * fluxeps)) - 1) < 0.1

        # MOSAIC L3 TESTS
        elif isinstance(si_model, MosaicModel):
            # Ensure that every pixel's poisson variance has increased or
            # remained the same with the new sources injected
            # Numpy isclose is needed to determine equality,
            # due to float precision issues
            close_mask = np.isclose(si_model.var_poisson, data_orig.var_poisson,
                                    rtol=1e-06)
            assert False in close_mask
            assert np.all(si_model.var_poisson[~close_mask] >
                          data_orig.var_poisson[~close_mask])

            # Ensure that every data pixel value has increased or
            # remained the same with the new sources injected
            assert np.all(si_model.data[~close_mask] >= data_orig.data[~close_mask])

            # Ensure total added flux matches expected added flux
            # maggies to counts (large number)
            cps_conv = bandpass.get_abflux(FILTER, 2)
            # electrons to mjysr (roughly order unity in scale)
            unit_factor = bandpass.etomjysr(FILTER, 2)
            total_rec_flux = np.sum(si_model.data - data_orig.data)  # MJy / sr
            total_theo_flux = 4 * MAG_FLUX * cps_conv * unit_factor  # u.MJy / u.sr
            assert np.isclose(total_rec_flux, total_theo_flux, rtol=4e-02)
