import logging
import functools
import warnings

import numpy as np
from astropy import units as u
from romancal.lib import dqflags

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

# Conversion factor from MJy/sr to uJy/arcsec^2
MJSR_TO_UJA2 = (u.megajansky / u.steradian).to(u.microjansky / u.arcsecond / u.arcsecond)

# Conversion factor from square arcseconds to steradians
A2_TO_SR = (np.pi / (180. * 3600.))**2


def photom_io(input_model, photom_metadata):
    """
    Short Summary
    -------------
    Combine photometric scalar conversion factors and add to the science metadata.

    Parameters
    ----------
    input_model : Roman level 2 image datamodel
        input roman datamodel

    photom_metadata : list
        Set of matched photom meta data keys

    Returns
    -------

    """
    # Get the scalar conversion factor.
    conversion = photom_metadata['photmjsr']  # unit is MJy / sr

    # Store the conversion factor in the meta data
    log.info(f'photmjsr value: {conversion:.6g}')
    input_model.meta.photometry.conversion_megajanskys = conversion
    input_model.meta.photometry.conversion_microjanskys = conversion.value * MJSR_TO_UJA2 * u.microjansky / (u.arcsecond**2)

    # Get the scalar conversion uncertainty factor
    uncertainty_conv = photom_metadata['uncertainty']

    # Store the uncertainty conversion factor in the meta data
    log.info(f'uncertainty value: {conversion:.6g}')
    input_model.meta.photometry.conversion_megajanskys_uncertainty = uncertainty_conv
    input_model.meta.photometry.conversion_microjanskys_uncertainty = uncertainty_conv.value * MJSR_TO_UJA2 * u.microjansky / (u.arcsecond**2)

    # Return updated input model
    return input_model


def find_photom_parameters(photom_table, match_fields):
    """
    Find metadata keyword object matching fields.

    Parameters
    ----------
    photom_table : dict
        "Optical Element":
            {"photmjsr": value * u.MJ / u.sr,
             "uncertainty": value* u.MJ / u.sr,
             "pixelareasr": value * u.sr}
    match_fields : dict
        {field_name: value} pair to use as a matching criteria.

    Raises
    ------
    Warning
        When a field name is not in the table.
    MatchFitsTableRowError
        When more than one rows match.

    Returns
    -------
    row : dict
        Matched metadata keyword object
    """

    # Find photom entries that match the optical element in match_fields
    results = [
        item[1].upper() in dict(photom_table).keys() for item in
        match_fields.items()
    ]

    # Trim to matched rows and test that we have exactly one match
    photom_parameters_match = list(filter(bool, results))
    if len(photom_parameters_match) > 1:
        raise MatchFitsTableRowError(f"Expected to find one matching row in table, found {len(photom_parameters_match)}.")
    if len(photom_parameters_match) == 0:
        warnings.warn("Expected to find one matching row in table, found 0.")
        return None

    # Return photom parameters for the selected optical element
    return photom_table[match_fields['optical_element'].upper()]



def save_area_info(input_model, photom_parameters):
    """
    Short Summary
    -------------
    Read the pixel area value in the photom parameters, then convert and
    copy them to the metadata of the input datamodel.

    Parameters
    ----------
    input_model : Roman level 2 image datamodel
        input roman datamodel

    photom_parameters : dict
        Photom parameter object
    """

    # Load the average pixel area values from the photom reference file header
    area_ster = photom_parameters["pixelareasr"]
    area_a2 = photom_parameters["pixelareasr"] * ((u.arcsecond * u.arcsecond) / u.steradian) / A2_TO_SR

    # Copy the pixel area values to the input model
    log.debug('pixelarea_steradians = %s, pixelarea_arcsecsq = %s', str(area_ster), str(area_a2))
    input_model.meta.photometry.pixelarea_arcsecsq = area_a2
    input_model.meta.photometry.pixelarea_steradians = area_ster

    # Return updated input model
    return input_model


def apply_photom(input_model, photom):
    """
    Short Summary
    -------------
    Open the reference file, retrieve the conversion factor from the reference
    file that is appropriate to the instrument mode. The scalar factor is written
    to the photmjsr and uncertainty keywords in the model.

    For WFI, matching is based on optical_element.

    Parameters
    ----------
    input_model : Roman level 2 image datamodel
        input roman datamodel

    photom : Roman Photom reference datamodel
        photom roman datamodel

    Returns
    -------
    output_model : Roman level 2 image datamodel
        output roman datamodel with the flux calibration keywords updated
    """

    # Obtain photom parameters for the selected optical element
    fields_to_match = {'optical_element': input_model.meta.instrument.optical_element}
    photom_parameters = find_photom_parameters(photom.phot_table, fields_to_match)
    if photom_parameters is None:
        return

    # Copy pixel area information to input datamodel
    input_model = save_area_info(input_model, photom_parameters)

    # Copy conversions to input model
    input_model = photom_io(input_model, photom_parameters)

    # Return updated input model
    return input_model
