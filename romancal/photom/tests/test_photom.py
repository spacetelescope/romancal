import os
import pytest
import numpy as np
import warnings
from astropy import units as u

from romancal.photom import photom, PhotomStep
from roman_datamodels.datamodels import ImageModel, WfiImgPhotomRefModel
from roman_datamodels.testing import utils as testutil


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
    optical_element = ["F062", "F087", "F106", "F129", "F146", "F158", "F184", "F213",
                       "GRISM", "PRISM", "DARK"]
    none_type_elements = ["GRISM", "PRISM", "DARK"]
    keyword = ["photmjsr", "uncertainty", "pixelareasr"]
    nrows = len(optical_element)

    # Create sample photometry keyword values
    photmjsr = np.linspace(min_r, min_r + (nrows - 1.) * delta, nrows) * u.megajansky / u.steradian
    uncertainty = np.linspace(min_r/20.0, min_r/20.0 + (nrows - 1.) * delta/20.0, nrows) * \
                  u.megajansky / u.steradian

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
        # GRISM, PRISM, and DARK optical elements shuld have their
        # photomeetry keywords set to None
        if element in none_type_elements:
            key_dict["photmjsr"] = None
            key_dict["uncertainty"] = None
        reftab[element] = key_dict

    # Create default datamodel
    photom_model = testutil.mk_wfi_img_photom()

    # Copy values above into defautl datamodel
    photom_model.phot_table=reftab

    return photom_model


def test_no_photom_match():
    """Test apply_photom warning for no match"""

    # Create sample WFI Level 2 science datamodel
    input_model = testutil.mk_level2_image()

    # Create photom reference datamodel
    photom_model = create_photom_wfi_image(min_r=3.1, delta=0.1)

    # Remove key for failed test (that won't fail validation)
    photom_model.phot_table.pop("F146")

    # Select optical element
    input_model.meta.instrument.optical_element = "F146"

    # Set bad values which would be overwritten by apply_photom
    input_model.meta.photometry.pixelarea_steradians = -1.0 * u.sr
    input_model.meta.photometry.conversion_megajanskys = -1.0 * u.megajansky / u.steradian
    input_model.meta.photometry.conversion_microjanskys_uncertainty = \
        -1.0 * u.microjansky / u.arcsecond ** 2

    with warnings.catch_warnings(record=True) as caught:
        # Look for now non existent F146 optical element
        output_model = photom.apply_photom(input_model, photom_model)

        # Assert warning key matches that of the input file
        assert str(caught[0].message).split()[-1] == input_model.meta.instrument.optical_element

        # Assert that photom elements are not updated
        assert output_model.meta.photometry.pixelarea_steradians == -1.0 * u.sr
        assert output_model.meta.photometry.conversion_megajanskys == \
               -1.0 * u.megajansky / u.steradian
        assert output_model.meta.photometry.conversion_microjanskys_uncertainty == \
               -1.0 * u.microjansky / u.arcsecond ** 2


def test_apply_photom1():
    """Test apply_photom applies correct metadata"""

    # Create sample WFI Level 2 science datamodel
    input_model = testutil.mk_level2_image()

    # Create photom reference datamodel
    photom_model = create_photom_wfi_image(min_r=3.1, delta=0.1)

    # Select optical element
    input_model.meta.instrument.optical_element = "F146"

    # Apply photom correction for optical element F146
    output_model = photom.apply_photom(input_model, photom_model)

    # Set reference photometry
    area_ster = 2.31307642258977E-14 * u.steradian
    area_a2 = 0.000984102303070964 * u.arcsecond * u.arcsecond

    # Tests for pixel areas
    assert (np.isclose(output_model.meta.photometry.pixelarea_steradians.value,
                        area_ster.value, atol=1.e-7))
    assert output_model.meta.photometry.pixelarea_steradians.unit == area_ster.unit
    assert (np.isclose(output_model.meta.photometry.pixelarea_arcsecsq.value,
                        area_a2.value, atol=1.e-7))
    assert output_model.meta.photometry.pixelarea_arcsecsq.unit == area_a2.unit

    # Set reference photometry
    phot_ster = 3.5 * u.megajansky / u.steradian
    phot_a2 = phot_ster.to(u.microjansky / u.arcsecond ** 2)

    # Tests for photometry
    assert (np.isclose(output_model.meta.photometry.conversion_megajanskys.value,
                       phot_ster.value, atol=1.e-7))
    assert output_model.meta.photometry.conversion_megajanskys.unit == phot_ster.unit
    assert (np.isclose(output_model.meta.photometry.conversion_microjanskys.value,
                       phot_a2.value, atol=1.e-7))
    assert output_model.meta.photometry.conversion_microjanskys.unit == phot_a2.unit

    # Set reference photometric uncertainty
    muphot_ster = 0.175 * u.megajansky / u.steradian
    muphot_a2 = muphot_ster.to(u.microjansky / u.arcsecond ** 2)

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
@pytest.mark.skipif(
    os.environ.get("CI") == "true",
    reason="Roman CRDS servers are not currently available outside the internal network"
)
def test_photom_step_interface(instrument, exptype):
    """Test that the basic inferface works for data requiring a photom reffile"""

    # Create a small area for the file
    shape = (20, 20)

    # Create input model
    wfi_image = testutil.mk_level2_image(shape=shape)
    wfi_image_model = ImageModel(wfi_image)

    # Create photom model
    photom = testutil.mk_wfi_img_photom()
    photom_model = WfiImgPhotomRefModel(photom)

    # Run photom correction step
    result = PhotomStep.call(wfi_image_model, override_photom=photom_model)

    assert (result.data == wfi_image.data).all()
    assert result.data.shape == shape
    if exptype == "WFI_IMAGE":
        assert result.meta.cal_step.photom == 'COMPLETE'
    else:
        assert result.meta.cal_step.photom == 'SKIPPED'


@pytest.mark.parametrize(
    "instrument, exptype",
    [
        ("WFI", "WFI_PRISM"),
    ]
)
@pytest.mark.skipif(
    os.environ.get("CI") == "true",
    reason="Roman CRDS servers are not currently available outside the internal network"
)
def test_photom_step_interface_spectroscopic(instrument, exptype):
    """Test apply_photom properly populates photometric keywords for spectroscopic data"""

    # Create a small area for the file
    shape = (20, 20)

    # Create input node
    wfi_image = testutil.mk_level2_image(shape=shape)

    # Select exposure type and optical element
    wfi_image.meta.exposure.type = "WFI_PRISM"
    wfi_image.meta.instrument.optical_element = "PRISM"

    # Set photometric values for spectroscopic data
    wfi_image.meta.photometry.pixelarea_steradians = 2.31307642258977E-14 * u.steradian
    wfi_image.meta.photometry.pixelarea_arcsecsq = 0.000984102303070964 * u.arcsecond * u.arcsecond
    wfi_image.meta.photometry.conversion_megajanskys = None
    wfi_image.meta.photometry.conversion_megajanskys_uncertainty = None
    wfi_image.meta.photometry.conversion_microjanskys = None
    wfi_image.meta.photometry.conversion_microjanskys_uncertainty = None

    # Create input model
    wfi_image_model = ImageModel(wfi_image)

    # Create photom model
    photom = testutil.mk_wfi_img_photom()
    photom_model = WfiImgPhotomRefModel(photom)

    # Run photom correction step
    result = PhotomStep.call(wfi_image_model, override_photom=photom_model)

    # Test that the data has not changed
    assert (np.allclose(result.data, wfi_image_model.data, rtol=1.e-7))

    # Test that keywords are properly preserved
    assert result.meta.photometry.conversion_megajanskys is None
    assert result.meta.photometry.conversion_microjanskys is None
    assert result.meta.photometry.conversion_megajanskys_uncertainty is None
    assert result.meta.photometry.conversion_microjanskys_uncertainty is None

    # Set reference pixel areas
    area_ster = 2.31307642258977E-14 * u.steradian
    area_a2 = 0.000984102303070964 * u.arcsecond * u.arcsecond

    # Tests for pixel areas
    assert (np.isclose(result.meta.photometry.pixelarea_steradians.value,
                        area_ster.value, atol=1.e-7))
    assert result.meta.photometry.pixelarea_steradians.unit == area_ster.unit
    assert (np.isclose(result.meta.photometry.pixelarea_arcsecsq.value,
                        area_a2.value, atol=1.e-7))
    assert result.meta.photometry.pixelarea_arcsecsq.unit == area_a2.unit
