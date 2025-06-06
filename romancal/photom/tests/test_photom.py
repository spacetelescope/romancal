import numpy as np
import pytest
from astropy import units as u
from roman_datamodels import stnode
from roman_datamodels.datamodels import ImageModel, WfiImgPhotomRefModel

from romancal.photom import PhotomStep, photom


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
    optical_element = [
        "F062",
        "F087",
        "F106",
        "F129",
        "F146",
        "F158",
        "F184",
        "F213",
        "GRISM",
        "PRISM",
        "DARK",
    ]
    none_type_elements = ["GRISM", "PRISM", "DARK"]
    keyword = ["photmjsr", "uncertainty", "pixelareasr"]
    nrows = len(optical_element)

    # Create sample photometry keyword values
    photmjsr = np.linspace(min_r, min_r + (nrows - 1.0) * delta, nrows)
    uncertainty = np.linspace(
        min_r / 20.0, min_r / 20.0 + (nrows - 1.0) * delta / 20.0, nrows
    )

    # Create sample area keyword values
    area_ster = 2.31307642258977e-14
    pixelareasr = np.ones(nrows, dtype=np.float64) * area_ster

    # Bundle values into a list
    values = list(zip(photmjsr, uncertainty, pixelareasr, strict=False))

    # Create dictionary containing all values
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
    photom_model = stnode.WfiImgPhotomRef.create_fake_data()

    # Copy values above into defautl datamodel
    photom_model.phot_table = reftab

    return photom_model


def test_no_photom_match():
    """Test apply_photom warning for no match"""

    # Create sample WFI Level 2 science datamodel
    input_model = stnode.WfiImage.create_fake_data(shape=(20, 20))

    # Create photom reference datamodel
    photom_model = create_photom_wfi_image(min_r=3.1, delta=0.1)

    # Remove key for failed test (that won't fail validation)
    photom_model.phot_table.pop("F146")

    # Select optical element
    input_model.meta.instrument.optical_element = "F146"

    # Set bad values which would be overwritten by apply_photom
    input_model.meta.photometry.pixel_area = -1.0
    input_model.meta.photometry.conversion_megajanskys = -1.0

    with pytest.warns(
        UserWarning,
        match=f"No matching photom parameters for {input_model.meta.instrument.optical_element}",
    ):
        # Look for now non existent F146 optical element
        output_model = photom.apply_photom(input_model, photom_model)

    # Assert that photom elements are not updated
    assert output_model.meta.photometry.pixel_area == -1.0
    assert output_model.meta.photometry.conversion_megajanskys == -1.0


def test_apply_photom1():
    """Test apply_photom applies correct metadata"""

    # Create sample WFI Level 2 science datamodel
    input_model = stnode.WfiImage.create_fake_data(shape=(20, 20))

    # Create photom reference datamodel
    photom_model = create_photom_wfi_image(min_r=3.1, delta=0.1)

    # Select optical element
    input_model.meta.instrument.optical_element = "F146"

    # Apply photom correction for optical element F146
    output_model = photom.apply_photom(input_model, photom_model)

    # Set reference photometry
    area_ster = 2.31307642258977e-14

    # Tests for pixel areas
    assert np.isclose(
        output_model.meta.photometry.pixel_area,
        area_ster,
        atol=1.0e-7,
    )

    # Set reference photometry
    phot_ster = 3.5

    # Tests for photometry
    assert np.isclose(
        output_model.meta.photometry.conversion_megajanskys,
        phot_ster,
        atol=1.0e-7,
    )

    # Set reference photometric uncertainty
    muphot_ster = 0.175

    # Tests for photometric uncertainty
    assert np.isclose(
        output_model.meta.photometry.conversion_megajanskys_uncertainty,
        muphot_ster,
        atol=1.0e-7,
    )


def test_apply_photom2():
    """Test apply_photom does not change data values"""

    # Create sample WFI Level 2 science datamodel
    input_model = stnode.WfiImage.create_fake_data(shape=(20, 20))

    # Create photom reference datamodel
    photom_model = create_photom_wfi_image(min_r=3.1, delta=0.1)

    # Apply photom correction
    output_model = photom.apply_photom(input_model, photom_model)

    # Select pixel for comparison
    shape = input_model.data.shape
    ix = shape[1] // 2
    iy = shape[0] // 2

    # Test that the data has not changed
    assert np.allclose(output_model.data[iy, ix], input_model.data[iy, ix], rtol=1.0e-7)


@pytest.mark.parametrize(
    "instrument, exptype",
    [
        ("WFI", "WFI_IMAGE"),
    ],
)
def test_photom_step_interface(instrument, exptype):
    """Test that the basic inferface works for data requiring a photom reffile"""

    # Create a small area for the file
    shape = (20, 20)

    # Create input model
    wfi_image = stnode.WfiImage.create_fake_data(shape=shape)
    wfi_image_model = ImageModel(wfi_image)
    wfi_image_model.meta.cal_step = stnode.L2CalStep.create_fake_data()
    wfi_image_model.meta.cal_logs = stnode.CalLogs.create_fake_data()

    # Create photom model
    photom = stnode.WfiImgPhotomRef.create_fake_data()
    photom_model = WfiImgPhotomRefModel(photom)

    # Run photom correction step
    result = PhotomStep.call(wfi_image_model, override_photom=photom_model)

    assert (result.data == wfi_image.data).all()
    assert result.data.shape == shape
    if exptype == "WFI_IMAGE":
        assert result.meta.cal_step.photom == "COMPLETE"
    else:
        assert result.meta.cal_step.photom == "SKIPPED"

    # Run photom correction step with reffile as N/A
    result = PhotomStep.call(wfi_image_model, override_photom="N/A")

    assert result.meta.cal_step.photom == "SKIPPED"


@pytest.mark.parametrize(
    "instrument, exptype",
    [
        ("WFI", "WFI_PRISM"),
    ],
)
def test_photom_step_interface_spectroscopic(instrument, exptype):
    """
    Test apply_photom properly populates photometric keywords for spectroscopic data
    """

    # Create a small area for the file
    shape = (20, 20)

    # Create input node
    wfi_image = stnode.WfiImage.create_fake_data(shape=shape)

    # Select exposure type and optical element
    wfi_image.meta.exposure.type = "WFI_PRISM"
    wfi_image.meta.instrument.optical_element = "PRISM"

    # Set photometric values for spectroscopic data
    wfi_image.meta.photometry.pixel_area = (2.31307642258977e-14 * u.steradian).value
    wfi_image.meta.photometry.conversion_megajanskys = (
        -99999 * u.megajansky / u.steradian
    ).value
    wfi_image.meta.photometry.conversion_megajanskys_uncertainty = (
        -99999 * u.megajansky / u.steradian
    ).value

    # Create input model
    wfi_image_model = ImageModel(wfi_image)
    wfi_image_model.meta.cal_step = stnode.L2CalStep.create_fake_data()
    wfi_image_model.meta.cal_logs = stnode.CalLogs.create_fake_data()

    # Create photom model
    photom = stnode.WfiImgPhotomRef.create_fake_data()
    photom_model = WfiImgPhotomRefModel(photom)

    # Run photom correction step
    result = PhotomStep.call(wfi_image_model, override_photom=photom_model)

    # Test that the data has not changed
    assert np.allclose(result.data, wfi_image_model.data, rtol=1.0e-7)

    # Test that keywords are properly overwritten
    assert result.meta.photometry.conversion_megajanskys is None
    assert result.meta.photometry.conversion_megajanskys_uncertainty is None
    assert result.meta.photometry.pixel_area is None
