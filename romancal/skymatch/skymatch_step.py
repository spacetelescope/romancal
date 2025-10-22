"""
Roman step for sky matching.
"""

from __future__ import annotations

import logging
from copy import deepcopy
from typing import TYPE_CHECKING

import numpy as np
from astropy.nddata.bitmask import bitfield_to_boolean_mask, interpret_bit_flags
from roman_datamodels.dqflags import pixel
from stcal.skymatch import SkyImage, SkyStats, skymatch

from romancal.datamodels import ModelLibrary
from romancal.stpipe import RomanStep

if TYPE_CHECKING:
    from typing import ClassVar

__all__ = ["SkyMatchStep"]

log = logging.getLogger(__name__)


class SkyMatchStep(RomanStep):
    """
    SkyMatchStep: Subtraction or equalization of sky background in science images.
    """

    class_alias = "skymatch"

    spec = """
        # General sky matching parameters:
        skymethod = option('local', 'global', 'match', 'global+match', default='match') # sky computation method
        match_down = boolean(default=True) # adjust sky to lowest measured value?
        subtract = boolean(default=False) # subtract computed sky from image data?

        # Image's bounding polygon parameters:
        stepsize = integer(default=None) # Max vertex separation

        # Sky statistics parameters:
        skystat = option('median', 'midpt', 'mean', 'mode', default='mode') # sky statistics
        dqbits = string(default='~DO_NOT_USE+NON_SCIENCE') # "good" DQ bits
        lower = float(default=None) # Lower limit of "good" pixel values
        upper = float(default=None) # Upper limit of "good" pixel values
        nclip = integer(min=0, default=5) # number of sky clipping iterations
        lsigma = float(min=0.0, default=4.0) # Lower clipping limit, in sigma
        usigma = float(min=0.0, default=4.0) # Upper clipping limit, in sigma
        binwidth = float(min=0.0, default=0.1) # Bin width for 'mode' and 'midpt' `skystat`, in sigma
    """

    reference_file_types: ClassVar = []

    def process(self, input):
        log.setLevel(logging.DEBUG)

        if isinstance(input, ModelLibrary):
            library = input
        else:
            library = ModelLibrary(input)

        self._dqbits = interpret_bit_flags(self.dqbits, flag_name_map=pixel)

        # set sky statistics:
        self._skystat = SkyStats(
            skystat=self.skystat,
            lower=self.lower,
            upper=self.upper,
            nclip=self.nclip,
            lsig=self.lsigma,
            usig=self.usigma,
            binwidth=self.binwidth,
        )

        # create a list of "Sky" Images and/or Groups:
        images = []
        with library:
            for model in library:
                images.append(self._imodel2skyim(model))

            # match/compute sky values:
            skymatch(
                images,
                skymethod=self.skymethod,
                match_down=self.match_down,
                subtract=self.subtract,
            )

            # set sky background value in each image's meta:
            for im in images:
                if isinstance(im, SkyImage):
                    self._set_sky_background(
                        im, "COMPLETE" if im.is_sky_valid else "SKIPPED"
                    )
                else:
                    for gim in im:
                        self._set_sky_background(
                            gim, "COMPLETE" if gim.is_sky_valid else "SKIPPED"
                        )
            for index, image in enumerate(images):
                library.shelve(image.meta["image_model"], index)

        return library

    def _imodel2skyim(self, image_model):
        input_image_model = image_model

        if "background" not in image_model.meta:
            image_model.meta["background"] = {
                "level": None,
                "subtracted": None,
                "method": None,
            }

        if self._dqbits is None:
            dqmask = np.isfinite(image_model.data).astype(dtype=np.uint8)
        else:
            dqmask = bitfield_to_boolean_mask(
                image_model.dq, self._dqbits, good_mask_value=1, dtype=np.uint8
            ) * np.isfinite(image_model.data)

        # see if 'skymatch' was previously run and raise an exception
        # if 'subtract' mode has changed compared to the previous pass:
        level = image_model.meta.background.level

        if image_model.meta.background.subtracted and level is None:
            # NOTE: In principle we could assume that level is 0 and
            # possibly add a log entry documenting this, however,
            # at this moment I think it is safer to quit and...
            #
            # report inconsistency:
            raise ValueError(
                "Background level was subtracted but the "
                "'level' property is undefined (None)."
            )

        if image_model.meta.background.subtracted and self.subtract:
            # cannot run 'skymatch' step on already "skymatched" images
            # when 'subtract' spec is inconsistent with
            # meta.background.subtracted:
            raise ValueError(
                "'subtract' step's specification is "
                "inconsistent with background info already "
                f"present in image '{image_model.meta.filename:s}' meta."
            )

        wcs = deepcopy(image_model.meta.wcs)

        sky_im = SkyImage(
            image=image_model.data,
            wcs_fwd=wcs.forward_transform,
            wcs_inv=wcs.backward_transform,
            pix_area=1.0,  # TODO: pixel area
            convf=1.0,  # TODO: conv. factor to brightness
            mask=dqmask,
            sky_id=image_model.meta.filename,  # file name?
            skystat=self._skystat,
            stepsize=self.stepsize,
            reduce_memory_usage=False,
            meta={"image_model": input_image_model},
        )

        if self.subtract and level is not None:
            sky_im.sky = level

        return sky_im

    def _set_sky_background(self, sky_image, step_status):
        image = sky_image.meta["image_model"]
        sky = sky_image.sky if sky_image.sky is not None else 0

        image.meta.background.method = str(self.skymethod)
        image.meta.background.subtracted = self.subtract
        # In numpy 2, the dtypes are more carefully controlled, so to match the
        # schema the data type needs to be re-cast to float64.
        image.meta.background.level = sky

        if step_status == "COMPLETE" and self.subtract:
            image.data[...] = sky_image.image[...]

        image.meta.cal_step.skymatch = step_status
