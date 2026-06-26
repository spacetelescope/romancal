"""Set Telescope Pointing from Observatory Engineering Telemetry.

Calculate and update the pointing-related and world coordinate system-related
keywords. Given a time period, usually defined by an exposure, the engineering
mnemonic database is queried for observatory orientation. The orientation
defines the sky coordinates that a particular point on the observatory is
pointed to. Then, using a set of matrix transformations, the sky coordinates of
the reference pixel of a desired aperture is calculated.

The transformations are defined by the STScI Innerspace (non-public) document
titles "Quaternion Transforms for Coarse Pointing WCS". The code itself follows
a demonstrative jupyter notebook.

**Interface**

The primary usage is through the command line interface
``roman_set_telescope_pointing``. Operating on a list of Roman Level 1 exposures,
this command updates the world coordinate system keywords with the values
necessary to translate from aperture pixel to sky coordinates.

Access to the Roman Engineering Mnemonic database is required. See the
:ref:`Engineering Database Interface <engdb>` for more information.

Programmatically, the command line is implemented by the function
`~roman.orientation.set_telescope_pointing.add_wcs`, which calls the basic function
`~roman.orientation.set_telescope_pointing.calc_wcs`.

There are two data structures used to maintain the state of the transformation.
`~roman.orientation.set_telescope_pointing.TransformParameters` contains the parameters
needed to perform the transformations.
`~roman.orientation.set_telescope_pointing.Transforms` contains the calculated
transformation matrices.

**Transformation Matrices**

All the transformation matrices, as defined by
`~roman.orientation.set_telescope_pointing.Transforms`, are Direction Cosine Matrices
(DCM). A DCM contains the Euler rotation angles that represent the sky
coordinates for a particular frame-of-reference. The initial DCM is provided
through the engineering telemetry and represents how the observatory is orientated.

**META Affected**

The following meta values are populated:

    - meta.guide_star.hv_position
    - meta.pointing.dec_v1
    - meta.pointing.pa_aperture
    - meta.pointing.pa_v3
    - meta.pointing.quaternion
    - meta.pointing.ra_v1
    - meta.pointing.target_dec
    - meta.pointing.target_ra
    - meta.wcsinfo.dec_ref
    - meta.wcsinfo.ra_ref
    - meta.wcsinfo.roll_ref
    - meta.wcsinfo.s_region

"""

import logging

import roman_datamodels as rdm

from . import _transforms as tlib
from . import _wcs as wlib

__all__ = [
    "add_wcs",
    "update_wcs",
]

# Setup logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())
DEBUG_FULL = logging.DEBUG - 1
LOGLEVELS = [logging.INFO, logging.DEBUG, DEBUG_FULL]

# Datamodels that can be updated, normally
EXPECTED_MODELS = rdm.datamodels.ScienceRawModel


def add_wcs(
    filename,
    dry_run=False,
    save_transforms=None,
    **transform_kwargs,
):
    """Add WCS information to a Roman DataModel.

    Telescope orientation is attempted to be obtained from
    the engineering database. Failing that, a default pointing
    is used based on input default information.

    The file is updated in-place.

    Parameters
    ----------
    filename : Path-like
        The path to a data file.

    dry_run : bool
        Run through the calculations but do not modify the file.

    save_transforms : Path-like or None
        File to save the calculated transforms to.

    transform_kwargs : dict
        dict to use to initialize the `TransformParameters` object.
        See `TransformParameters` for more information.`

    Notes
    -----
    This function adds absolute pointing information to the Roman datamodels
    provided. By default, only Level 1 exposures are allowed to be updated.
    These have the suffixes of "uncal" representing datamodel ScienceRawModel.
    Any higher level product, from Level 2 and beyond, that has had the
    `assign_wcs` step applied, have improved WCS information. Running this task
    on such files will potentially corrupt the WCS.
    """
    logger.info("Updating WCS info for file %s", filename)

    with rdm.open(filename) as model:
        if not isinstance(model, EXPECTED_MODELS):
            logger.warning("Input %s is not of an expected type (uncal)", model)
            logger.warning(
                "    Updating pointing may have no effect or detrimental effects on the "
                "WCS information,"
            )
            logger.warning(
                "    especially if the input is the result of Level2b or higher calibration."
            )
            raise TypeError(
                f"Input model {model} is not one of {EXPECTED_MODELS}."
                "\n\tFailing WCS processing."
            )

        t_pars, transforms = update_wcs(
            model,
            **transform_kwargs,
        )

        if dry_run:
            logger.info("Dry run requested; results are not saved.")
        else:
            logger.info("Saving updated model %s", filename)
            model.save(filename)
            if transforms and save_transforms:
                logger.info("Saving transform matrices to %s", save_transforms)
                transforms.write_to_asdf(save_transforms)

    logger.info("...update completed")


def update_wcs(
    model,
    **transform_kwargs,
):
    """
    Update WCS pointing information.

    Given a `roman.datamodels.DataModel`, determine the simple WCS parameters
    from the SIAF keywords in the model and the engineering parameters
    that contain information about the telescope orientation.

    It presumes all the accessed keywords are present.

    Parameters
    ----------
    model : `~roman.datamodels.DataModel`
        The model to update.

    **transform_kwargs : dict
        dict to use to initialize the `TransformParameters` object.
        See `TransformParameters` for more information.`

    Returns
    -------
    t_pars, transforms : TransformParameters, Transforms
        The parameters and transforms calculated. May be
        None for either if telemetry calculations were not
        performed.
    """
    # Configure transformation parameters.
    t_pars = tlib.TransformParameters(**transform_kwargs)
    tlib.t_pars_from_model(model, t_pars)
    logger.log(DEBUG_FULL, 'TransformParameters: %s', t_pars)

    # Calculate WCS.
    transforms = wlib.update_wcs_from_telem(model, t_pars)

    return t_pars, transforms
