"""
Roman step for sky matching.
"""

import logging
from copy import deepcopy

import numpy as np
from astropy.nddata.bitmask import bitfield_to_boolean_mask, interpret_bit_flags
from roman_datamodels import maker_utils
from roman_datamodels.dqflags import pixel

from romancal.datamodels import ModelLibrary
from romancal.stpipe import RomanStep

from .skyimage import SkyImage
from .skymatch import match
from .skystatistics import SkyStats

__all__ = ["SkyMatchStep"]


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
    """  # noqa: E501

    reference_file_types = []

    def process(self, input):
        self.log.setLevel(logging.DEBUG)

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
            for index, model in enumerate(library):
                images.append(self._imodel2skyim(model))

            # match/compute sky values:
            match(
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
            image_model.meta["background"] = maker_utils.mk_sky_background(
                level=None, subtracted=None, method="None"
            )

        if self._dqbits is None:
            dqmask = np.isfinite(image_model.data).astype(dtype=np.uint8)
        else:
            dqmask = bitfield_to_boolean_mask(
                image_model.dq, self._dqbits, good_mask_value=1, dtype=np.uint8
            ) * np.isfinite(image_model.data)

        # see if 'skymatch' was previously run and raise an exception
        # if 'subtract' mode has changed compared to the previous pass:
        level = image_model.meta.background.level
        if image_model.meta.background.subtracted is None:
            if level is not None:
                # report inconsistency:
                raise ValueError(
                    "Background level was set but the "
                    "'subtracted' property is undefined (None)."
                )
            # Assume level is zero:
            level = 0.0

        else:
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
                    "present in image '{:s}' meta.".format(image_model.meta.filename)
                )

        wcs = deepcopy(image_model.meta.wcs)

        sky_im = SkyImage(
            image=image_model.data,
            wcs_fwd=wcs.forward_transform,
            wcs_inv=wcs.backward_transform,
            pix_area=1.0,  # TODO: pixel area
            convf=1.0,  # TODO: conv. factor to brightness
            mask=dqmask,
            id=image_model.meta.filename,  # file name?
            skystat=self._skystat,
            stepsize=self.stepsize,
            reduce_memory_usage=False,
            meta={"image_model": input_image_model},
        )

        if self.subtract:
            sky_im.sky = level

        return sky_im

    def _set_sky_background(self, sky_image, step_status):
        image = sky_image.meta["image_model"]
        sky = sky_image.sky
        if sky == 0 or sky is None:
            sky = 0 * image.data.unit

        image.meta.background.method = str(self.skymethod)
        image.meta.background.subtracted = self.subtract
        image.meta.background.level = sky

        if step_status == "COMPLETE" and self.subtract:
            image.data[...] = sky_image.image[...]

        image.meta.cal_step.skymatch = step_status
