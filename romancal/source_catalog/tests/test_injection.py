"""
Unit tests for the Roman source injection step code
"""

import pytest

pytest.importorskip("romanisim")

import numpy as np
from astropy import table
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
from roman_datamodels.datamodels import ImageModel, MosaicModel
from romanisim import bandpass, parameters

from romancal.skycell.tests.test_skycell_match import mk_gwcs
from romancal.source_catalog import injection
from romancal.source_catalog.injection import inject_sources, make_cosmoslike_catalog

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
FILTERS = ["F062", "F158", "F213", "F146"]
RNG_SEED = 42
MATABLE = 4
BANDPASSES = set(bandpass.galsim2roman_bandpass.values())

# Create gaussian noise generators
# sky should generate ~0.2 electron / s / pix.
# MJy / sr has similar magnitude to electron / s (i.e., within a factor
# of several), so just use 0.2 here.
MEANFLUX = 0.2


def make_test_data():
    noise_rng = np.random.default_rng(RNG_SEED)

    # Populate the data array with gaussian noise
    data = noise_rng.normal(loc=MEANFLUX, scale=0.01 * MEANFLUX, size=SHAPE)
    err = np.ones_like(data) * 0.01 * MEANFLUX

    return data, err


@pytest.fixture
def image_model(filter=FILTERS[0]):
    defaults = {"meta": parameters.default_parameters_dictionary}
    defaults["meta"]["instrument"]["detector"] = DETECTOR
    defaults["meta"]["instrument"]["optical_element"] = filter
    defaults["meta"]["wcsinfo"]["ra_ref"] = RA
    defaults["meta"]["wcsinfo"]["dec_ref"] = DEC
    defaults["meta"]["wcsinfo"]["roll_ref"] = ROLL

    model = ImageModel.create_fake_data(defaults=defaults, shape=SHAPE)
    model.meta.exposure.read_pattern = parameters.read_pattern[MATABLE]

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
def mosaic_model(filter=FILTERS[0]):
    defaults = {
        "meta": {
            "coadd_info": {
                "time_first": Time("2024-01-01T12:00:00.000", format="isot"),
            },
            "instrument": {
                "optical_element": filter,
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
        metadata["wcsinfo"]["ra_ref"],
        metadata["wcsinfo"]["dec_ref"],
        metadata["wcsinfo"]["roll_ref"],
        bounding_box=((-0.5, SHAPE[0] - 0.5), (-0.5, SHAPE[1] - 0.5)),
        shape=SHAPE,
    )

    # Create normalized psf source catalog (same source in each quadrant)
    ra, dec = wcsobj.pixel_to_world_values(np.array(XPOS_IDX), np.array(YPOS_IDX))

    tabcat = table.Table()
    tabcat["ra"] = ra
    tabcat["dec"] = dec
    tabcat[metadata["instrument"]["optical_element"]] = len(XPOS_IDX) * [MAG_FLUX]
    tabcat["type"] = len(XPOS_IDX) * ["PSF"]
    tabcat["n"] = len(XPOS_IDX) * [-1]
    tabcat["half_light_radius"] = len(XPOS_IDX) * [-1] * u.arcsec
    tabcat["pa"] = len(XPOS_IDX) * [-1]
    tabcat["ba"] = len(XPOS_IDX) * [-1]

    return tabcat


# Ignore this warning - it comes from stpsf,
#   and they are aware of the need for an upgrade.
@pytest.mark.filterwarnings(
    "ignore:Python 3.14 will, by default,"
    " filter extracted tar archives and"
    " reject files or modify their metadata."
    " Use the filter argument to control this behavior.:DeprecationWarning"
)
def test_inject_sources(image_model, mosaic_model):
    for si_model in (image_model, mosaic_model):
        """Test simple source injection"""
        # Set filter
        test_filter = FILTERS[0]

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
        cps_conv = bandpass.get_abflux(test_filter, 2)
        # electrons to mjysr (roughly order unity in scale)
        if isinstance(si_model, ImageModel):
            unit_factor = 1 / parameters.reference_data["gain"].value
        else:
            unit_factor = bandpass.etomjysr(test_filter, 2)
        total_rec_flux = np.sum(si_model.data - data_orig.data)  # MJy / sr
        total_theo_flux = len(cat) * MAG_FLUX * cps_conv * unit_factor  # u.MJy / u.sr
        assert np.isclose(total_rec_flux, total_theo_flux, rtol=0.1)


def test_create_cosmoscat():
    # Set a test filter
    test_filter = FILTERS[0]

    # Pointing
    cen = SkyCoord(ra=RA * u.deg, dec=DEC * u.deg)

    # WCS object for ra & dec conversion
    wcsobj = mk_gwcs(
        RA,
        DEC,
        ROLL,
        bounding_box=((-0.5, SHAPE[0] - 0.5), (-0.5, SHAPE[1] - 0.5)),
        shape=SHAPE,
    )

    # Convert x,y to ra, dec
    ra, dec = wcsobj.pixel_to_world_values(np.array(XPOS_IDX), np.array(YPOS_IDX))

    # Exposure times (s)
    exptimes = {}
    for bp in FILTERS:
        exptimes[bp] = 300

    # Generate cosmos-like catalog
    cat = make_cosmoslike_catalog(
        cen, ra, dec, exptimes, filters=FILTERS, seed=RNG_SEED
    )

    # Set simple wcs metadata for mcat
    meta = {
        "wcsinfo": {
            "ra_ref": RA,
            "dec_ref": DEC,
            "roll_ref": ROLL,
        },
        "instrument": {
            "optical_element": test_filter,
        },
    }
    mcat = make_catalog(meta)

    # Ensure that locations are as expected
    assert np.allclose(np.sort(cat["ra"]), np.sort(mcat["ra"]), rtol=1.0e-6)
    assert np.allclose(np.sort(cat["dec"]), np.sort(mcat["dec"]), rtol=1.0e-6)

    # Ensure correct number of point sources
    assert np.sum(cat["type"] == "PSF") == int(len(XPOS_IDX) / 4)
    assert np.sum(cat["n"] == -1) == int(len(XPOS_IDX) / 4)

    # Set the point magnitude limit
    point_band_mag_limit = []
    for bp in FILTERS:
        # Normalize the mag limit to exptimes
        if bp in exptimes:
            point_band_mag_limit.append(
                injection.HRPOINTMAGLIMIT[bp]
                + (1.25 * np.log10((exptimes[bp] * u.s).to(u.hour).value))
            )
    point_mag_limit = max(point_band_mag_limit)

    for bp in FILTERS:
        # Ensure point fluxes in range
        assert np.all(
            cat[cat["type"] == "PSF"][bp] < 10.0 ** (-(point_mag_limit - 4) / 2.5)
        )
        assert np.all(
            cat[cat["type"] == "PSF"][bp] > 10.0 ** (-(point_mag_limit + 1) / 2.5)
        )

        # Ensure points lack color
        assert np.all(cat[cat["type"] == "PSF"][bp] == cat[cat["type"] == "PSF"][bp])

    # Ensure galaxy sizes are reasonable
    assert np.all(cat[cat["type"] == "SER"]["half_light_radius"].value < 1)
    assert np.all(cat[cat["type"] == "SER"]["half_light_radius"].value >= 0.036)

    # Set the galaxy magnitude limit
    gal_band_mag_limit = []
    for bp in FILTERS:
        # Normalize the mag limit to exptimes
        if bp in exptimes:
            gal_band_mag_limit.append(
                injection.HRGALMAGLIMIT[bp]
                + (1.25 * np.log10((exptimes[bp] * u.s).to(u.hour).value))
            )
    gal_mag_limit = max(gal_band_mag_limit)

    # Ensure galaxy fluxes in range
    for bp in FILTERS:
        assert np.all(
            cat[cat["type"] == "SER"][bp] < 10.0 ** (-(gal_mag_limit - 4) / 2.5)
        )
        assert np.all(cat[cat["type"] == "SER"][bp] >= 0)
