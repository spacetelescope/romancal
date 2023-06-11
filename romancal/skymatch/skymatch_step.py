"""
Roman step for sky matching.
"""

import logging
from copy import deepcopy
from itertools import chain

import numpy as np
from astropy.nddata.bitmask import bitfield_to_boolean_mask, interpret_bit_flags
from roman_datamodels import datamodels as rdd

from romancal.datamodels import ModelContainer
from romancal.lib.dqflags import pixel
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
        self._is_asn = False

        img = ModelContainer(
            input, save_open=not self._is_asn, return_open=not self._is_asn
        )

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

        # group images by their "group id":
        grp_img = chain.from_iterable(img.models_grouped)

        # create a list of "Sky" Images and/or Groups:
        images = [self._imodel2skyim(g) for grp_id, g in enumerate(grp_img, start=1)]

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

        return input if self._is_asn else img

    def _imodel2skyim(self, image_model):
        input_image_model = image_model

        if "background" not in image_model.meta:
            # TODO: remove this block when ``rad`` has a background schema.
            # This is a temporary workaround to insert a 'background'
            # entry into the metadata, which we'll later put into ``rad``:
            image_model.meta["background"] = dict(
                level=None, subtracted=None, method=None
            )
        if self._is_asn:
            image_model = rdd.open(image_model)

        if self._dqbits is None:
            dqmask = np.isfinite(image_model.data).astype(dtype=np.uint8)
        else:
            dqmask = bitfield_to_boolean_mask(
                image_model.dq, self._dqbits, good_mask_value=1, dtype=np.uint8
            ) * np.isfinite(image_model.data)

        # see if 'skymatch' was previously run and raise an exception
        # if 'subtract' mode has changed compared to the previous pass:
        level = image_model.meta["background"]["level"]
        if image_model.meta["background"]["subtracted"] is None:
            if level is not None:
                if self._is_asn:
                    image_model.close()

                # report inconsistency:
                raise ValueError(
                    "Background level was set but the "
                    "'subtracted' property is undefined (None)."
                )
            # Assume level is zero:
            level = 0.0

        else:
            if image_model.meta["background"]["subtracted"] and level is None:
                # NOTE: In principle we could assume that level is 0 and
                # possibly add a log entry documenting this, however,
                # at this moment I think it is safer to quit and...
                #
                # report inconsistency:
                if self._is_asn:
                    image_model.close()

                raise ValueError(
                    "Background level was subtracted but the "
                    "'level' property is undefined (None)."
                )

            if image_model.meta["background"]["subtracted"] and self.subtract:
                # cannot run 'skymatch' step on already "skymatched" images
                # when 'subtract' spec is inconsistent with
                # meta.background.subtracted:
                if self._is_asn:
                    image_model.close()

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
            reduce_memory_usage=self._is_asn,
            meta={"image_model": input_image_model},
        )

        if self._is_asn:
            image_model.close()

        if self.subtract:
            sky_im.sky = level

        return sky_im

    def _set_sky_background(self, sky_image, step_status):
        image = sky_image.meta["image_model"]
        sky = sky_image.sky

        if self._is_asn:
            dm = rdd.open(image)
        else:
            dm = image

        if step_status == "COMPLETE":
            # TODO: remove this block when ``rad`` has a background schema:
            # https://github.com/spacetelescope/rad/issues/247
            # This is a temporary workaround to access a 'background'
            # entry into metadata as a Python dict, which we'll later define with
            # a schema in ``rad``:
            dm.meta["background"]["method"] = str(self.skymethod)
            dm.meta["background"]["level"] = sky
            dm.meta["background"]["subtracted"] = (
                self.subtract or dm.meta["background"]["subtracted"]
            )

            if self.subtract:
                dm.data[...] = sky_image.image[...]

        dm.meta.cal_step.skymatch = step_status

        if self._is_asn:
            dm.save(image)
            dm.close()
