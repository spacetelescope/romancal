import math
import warnings

import pytest
import numpy as np

from astropy import units as u

from romancal.photom import photom, PhotomStep
from roman_datamodels.datamodels import ImageModel, WfiImgPhotomRefModel
from roman_datamodels.testing import utils as testutil

MJSR_TO_UJA2 = (u.megajansky / u.steradian).to(u.microjansky / (u.arcsecond**2))

# Multiply by this to convert from square arcseconds to steradians
A2_TO_SR = (np.pi / (180. * 3600.))**2


def create_photom_wfi_image(min_r=3.1, delta=0.1):
    """Create a photom table for WFI.

    Parameters
    ----------
    min_r : float
        Minimum value to assign when populating the photometry array.

    max_r : float
        Maximum value to assign when populating the photometry array.

    Returns
    -------
    photom_model : Roman Photom reference datamodel
        photom roman datamodel

    """

    # Keyword and eleents list for phot_table construction
    optical_element = ["F062", "F087", "F106", "F129", "W146", "F158", "F184", "F213",
                       "GRISM", "PRISM", "DARK"]
    none_type_elements = ["GRISM", "PRISM", "DARK"]
    keyword = ["photmjsr", "uncertainty", "pixelareasr"]
    nrows = len(optical_element)

    # Create sample photometry keyword values
    photmjsr = np.linspace(min_r, min_r + (nrows - 1.) * delta, nrows) * u.megajansky / u.steradian
    uncertainty = np.linspace(min_r/20.0, min_r/20.0 + (nrows - 1.) * delta/20.0, nrows) * u.megajansky / u.steradian

    # Create sample area keyword values
    area_ster = 2.31307642258977E-14 * u.steradian
    pixelareasr = np.ones(nrows, dtype=np.float32) * area_ster

    # Bundle values into a list
    values = list(zip(photmjsr, uncertainty, pixelareasr))

    #Create dictionary containing all values
    reftab = {}
    for element_idx, element in enumerate(optical_element):
        key_dict = {}
        for key_idx, key in enumerate(keyword):
            key_dict[key] = values[element_idx][key_idx]
        # GRISM, PRISM, and DARK optical elements shuld have their photomeetry keywords set to None
        if element in none_type_elements:
            key_dict["photmjsr"] = None
            key_dict["uncertainty"] = None
        reftab[element] = key_dict

    # Create default datamodel
    photom_model = testutil.mk_wfi_img_photom()

    # Copy values above into defautl datamodel
    photom_model.phot_table=reftab

    return photom_model


def test_find_photom_parameters():
    """ Test that find_photom_parameters correctly finds proper values
    (and reports mismatch errors)
    """

    # Create sample WFI Level 2 science datamodel
    input_model = testutil.mk_level2_image()

    # Create photom reference datamodel and set some new values
    photom_model = create_photom_wfi_image(min_r=3.1, delta=0.1)
    photom_model.phot_table['W146']['photmjsr'] = 2.71828 * u.megajansky / u.steradian
    photom_model.phot_table['W146']['uncertainty'] = 0.0031415 * u.megajansky / u.steradian
    photom_model.phot_table['W146']['pixelareasr'] = 1.764 * u.steradian

    # Obtain matched parameters for W146
    photom_parameters = photom.find_photom_parameters(photom_model.phot_table, {"optical_element": 'W146'})

    # Test that the expected values were obtained
    assert photom_parameters['photmjsr'] == 2.71828 * u.megajansky / u.steradian
    assert photom_parameters['uncertainty'] == 0.0031415 * u.megajansky / u.steradian
    assert photom_parameters['pixelareasr'] == 1.764 * u.steradian

    # Same tests utilizing lower case
    photom_parameters = photom.find_photom_parameters(photom_model.phot_table, {"optical_element": 'w146'})
    assert photom_parameters['photmjsr'] == 2.71828 * u.megajansky / u.steradian
    assert photom_parameters['uncertainty'] == 0.0031415 * u.megajansky / u.steradian
    assert photom_parameters['pixelareasr'] == 1.764 * u.steradian

    # Test for expected None values
    photom_model.phot_table['GRISM']['pixelareasr'] = 1.764 * u.steradian
    photom_parameters = photom.find_photom_parameters(photom_model.phot_table, {"optical_element": "GRISM"})
    assert photom_parameters['photmjsr'] is None
    assert photom_parameters['uncertainty'] is None
    assert photom_parameters['pixelareasr'] == 1.764 * u.steradian

    # Test for "no match" warnings
    with warnings.catch_warnings(record=True) as caught:
        # Look for non existent p314 optical element
        photom_parameters = photom.find_photom_parameters(photom_model.phot_table, {"optical_element": 'p314'})

        assert photom_parameters is None
        assert len(caught) == 1


def test_apply_photom1():
    """Test apply_photom applies correct metadata"""

    # Create sample WFI Level 2 science datamodel
    input_model = testutil.mk_level2_image()

    # Create photom reference datamodel
    photom_model = create_photom_wfi_image(min_r=3.1, delta=0.1)

    shape = input_model.data.shape

    input_model.meta.instrument.optical_element = "W146"

    # Apply photom correction for optical element F062
    output_model = photom.apply_photom(input_model, photom_model)

    # Set reference photometry
    area_ster = 2.31307642258977E-14 * u.steradian
    area_a2 = 0.000984102303070964 * u.arcsecond * u.arcsecond

    # Tests for pixel areas
    assert(np.isclose(output_model.meta.photometry.pixelarea_steradians.value,
                        area_ster.value, atol=1.e-7))
    assert output_model.meta.photometry.pixelarea_steradians.unit == area_ster.unit
    assert(np.isclose(output_model.meta.photometry.pixelarea_arcsecsq.value,
                        area_a2.value, atol=1.e-7))
    assert output_model.meta.photometry.pixelarea_arcsecsq.unit == area_a2.unit

    # Set reference photometry
    phot_ster = 3.5 * u.megajansky / u.steradian
    phot_a2 = phot_ster.value * MJSR_TO_UJA2 * u.microjansky / (u.arcsecond**2)

    # Tests for photometry
    assert (np.isclose(output_model.meta.photometry.conversion_megajanskys.value,
                       phot_ster.value, atol=1.e-7))
    assert output_model.meta.photometry.conversion_megajanskys.unit == phot_ster.unit
    assert (np.isclose(output_model.meta.photometry.conversion_microjanskys.value,
                       phot_a2.value, atol=1.e-7))
    assert output_model.meta.photometry.conversion_microjanskys.unit == phot_a2.unit

    # Set reference photometric uncertainty
    muphot_ster = 0.175 * u.megajansky / u.steradian
    muphot_a2 = muphot_ster.value * MJSR_TO_UJA2 * u.microjansky / (u.arcsecond ** 2)

    # Tests for photometric uncertainty
    assert (np.isclose(output_model.meta.photometry.conversion_megajanskys_uncertainty.value,
                       muphot_ster.value, atol=1.e-7))
    assert output_model.meta.photometry.conversion_megajanskys_uncertainty.unit == muphot_ster.unit
    assert (np.isclose(output_model.meta.photometry.conversion_microjanskys_uncertainty.value,
                       muphot_a2.value, atol=1.e-7))
    assert output_model.meta.photometry.conversion_microjanskys_uncertainty.unit == muphot_a2.unit


def test_apply_photom2():
    """Test apply_photom does not change data values"""

    # Create sample WFI Level 2 science datamodel
    input_model = testutil.mk_level2_image()

    # Create photom reference datamodel
    photom_model = create_photom_wfi_image(min_r=3.1, delta=0.1)

    # Apply photom correction
    output_model = photom.apply_photom(input_model, photom_model)

    # Select pixel for comparison
    shape = input_model.data.shape
    ix = shape[1] // 2
    iy = shape[0] // 2

    # Test that the data has not changed
    assert (np.allclose(output_model.data[iy, ix], input_model.data[iy, ix], rtol=1.e-7))


@pytest.mark.parametrize(
    "instrument, exptype",
    [
        ("WFI", "WFI_IMAGE"),
    ]
)
def test_photom_step_interface(instrument, exptype):
    """Test that the basic inferface works for data requiring a photom reffile"""

    # Create a small area for the file
    shape = (20, 20)

    # Create input model
    wfi_image = testutil.mk_level2_image(shape=shape)
    wfi_image_model = ImageModel(wfi_image)

    # Create photom model
    photom =  testutil.mk_wfi_img_photom()
    photom_model = WfiImgPhotomRefModel(photom)

    # Run photom correction step
    result = PhotomStep.call(wfi_image_model, override_photom=photom_model)

    assert (result.data == wfi_image.data).all()
    assert result.data.shape == shape
    assert result.meta.cal_step.photom == 'COMPLETE'
