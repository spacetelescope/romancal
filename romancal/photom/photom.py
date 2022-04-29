import logging
import warnings
from astropy import units as u

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def photom_io(input_model, photom_metadata):
    """
    Combine photometric scalar conversion factors and add to the science metadata.

    Parameters
    ----------
    input_model : Roman level 2 image datamodel
        Input Roman datamodel

    photom_metadata : dict
        Set of matched photom meta data keys

    Returns
    -------

    """
    # Get the scalar conversion factor.
    conversion = photom_metadata['photmjsr']  # unit is MJy / sr

    # Store the conversion factor in the meta data
    log.info(f'photmjsr value: {conversion:.6g}')
    input_model.meta.photometry.conversion_megajanskys = conversion
    input_model.meta.photometry.conversion_microjanskys = conversion.to(
        u.microjansky / u.arcsecond ** 2)

    # Get the scalar conversion uncertainty factor
    uncertainty_conv = photom_metadata['uncertainty']

    # Store the uncertainty conversion factor in the meta data
    log.info(f'uncertainty value: {uncertainty_conv:.6g}')
    input_model.meta.photometry.conversion_megajanskys_uncertainty = uncertainty_conv
    input_model.meta.photometry.conversion_microjanskys_uncertainty = uncertainty_conv.to(
        u.microjansky / u.arcsecond ** 2)

    # Return updated input model
    return input_model


def save_area_info(input_model, photom_parameters):
    """
    Read the pixel area value in the photom parameters, then convert and
    copy them to the metadata of the input datamodel.

    Parameters
    ----------
    input_model : Roman level 2 image datamodel
        Input Roman datamodel

    photom_parameters : dict
        Photom parameter object
    """

    # Load the average pixel area values from the photom reference file header
    area_ster = photom_parameters["pixelareasr"]
    area_a2 = photom_parameters["pixelareasr"].to(u.arcsecond**2)

    # Copy the pixel area values to the input model
    log.debug(f'pixelarea_steradians = {area_ster}, pixelarea_arcsecsq = {area_a2}')
    input_model.meta.photometry.pixelarea_arcsecsq = area_a2
    input_model.meta.photometry.pixelarea_steradians = area_ster

    # Return updated input model
    return input_model


def apply_photom(input_model, photom):
    """
    Retrieve the conversion factor from the photom reference datamodel
    that is appropriate to the instrument mode. The scalar factor is written
    to the photmjsr and uncertainty keywords in the model.

    For WFI, matching is based on optical_element.

    Parameters
    ----------
    input_model : Roman level 2 image datamodel
        Input Roman datamodel

    photom : Roman Photom reference datamodel
        Photom Roman datamodel

    Returns
    -------
    output_model : Roman level 2 image datamodel
        Output Roman datamodel with the flux calibration keywords updated
    """
    # Obtain photom parameters for the selected optical element
    try:
        photom_parameters = photom.phot_table[input_model.meta.instrument.optical_element.upper()]
    except KeyError:
        warnings.warn(f'No matching photom parameters for '
                      f'{input_model.meta.instrument.optical_element}')
        return input_model

    # Copy pixel area information to output datamodel
    output_model = save_area_info(input_model, photom_parameters)

    # Copy conversions to output model
    output_model = photom_io(output_model, photom_parameters)

    # Return updated output model
    return output_model
