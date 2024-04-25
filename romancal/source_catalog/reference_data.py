"""
Module for parsing APCORR and ABVEGAOFFSET reference file data.
"""

import logging

import numpy as np
from astropy.utils import lazyproperty
from roman_datamodels.datamodels import MosaicModel

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


class ReferenceData:
    """
    Class for reference file data needed by `SourceCatalogStep`.

    .. note::
        This class currently uses placeholder values as reference file
        data is not yet available.

    Parameters
    ----------
    model : `ImageModel`
        An `ImageModel` of drizzled image.

    aperture_ee : tuple of 3 int
        The aperture encircled energies to be used for aperture
        photometry. The values must be 3 strictly-increasing integers.

        Valid values will be defined in the APCORR reference files
        (e.g., 20, 30, 40, 50, 60, 70, or 80).

    Attributes
    ----------
    aperture_params : `dict`
        A dictionary containing the aperture parameters (radii, aperture
        corrections, and background annulus inner and outer radii).
    """

    def __init__(self, model, aperture_ee):
        if not isinstance(model, MosaicModel):
            raise ValueError("The input model must be a MosaicModel.")
        self.model = model

        self.aperture_ee = self._validate_aperture_ee(aperture_ee)

        self.instrument = self.model.meta.basic.instrument
        self.optical_element = self.model.meta.basic.optical_element
        log.info(f"Instrument: {self.instrument}")
        log.info(f"Optical Element: {self.optical_element}")

    @staticmethod
    def _validate_aperture_ee(aperture_ee):
        """
        Validate the input ``aperture_ee``.
        """
        aperture_ee = np.array(aperture_ee).astype(int)
        if not np.all(aperture_ee[1:] > aperture_ee[:-1]):
            raise ValueError("aperture_ee values must be strictly increasing")
        if len(aperture_ee) != 3:
            raise ValueError("aperture_ee must contain only 3 values")
        if np.any(np.logical_or(aperture_ee <= 0, aperture_ee >= 100)):
            raise ValueError("aperture_ee values must be between 0 and 100")
        return aperture_ee

    @lazyproperty
    def aperture_params(self):
        """
        A dictionary containing the aperture parameters (radii, aperture
        corrections, and background annulus inner and outer radii) for
        performing aperture photometry.
        """
        log.warning(
            "Reference file data mapping aperture encircled energies to "
            "radii (in pixels), background annuli sizes, and aperture "
            "corrections are not yet available. Using placeholder values."
        )

        return {
            "aperture_radii": np.array((1.0, 2.0, 3.0)),
            "aperture_corrections": np.array((1.0, 1.0, 1.0)),
            "aperture_ee": self.aperture_ee,
            "bkg_aperture_inner_radius": 5.0,
            "bkg_aperture_outer_radius": 10.0,
        }
