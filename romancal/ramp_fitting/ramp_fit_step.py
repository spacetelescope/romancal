#! /usr/bin/env python
#
from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import numpy as np
from roman_datamodels import datamodels as rdd
from roman_datamodels import stnode as rds
from roman_datamodels.dqflags import group, pixel
from stcal.ramp_fitting import ols_cas22_fit
from stcal.ramp_fitting.ols_cas22 import Parameter, Variance

from romancal.stpipe import RomanStep

if TYPE_CHECKING:
    from typing import ClassVar

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

__all__ = ["RampFitStep"]


class RampFitStep(RomanStep):
    """
    This step fits a straight line to the value of counts vs. time to
    determine the mean count rate for each pixel.
    """

    class_alias = "ramp_fit"

    spec = """
        algorithm = option('ols_cas22', default='ols_cas22')  # Algorithm to use to fit.
        save_opt = boolean(default=False) # Save optional output
        suffix = string(default='rampfit')  # Default suffix of results
        use_ramp_jump_detection = boolean(default=True) # Use jump detection during ramp fitting
        threshold_intercept = float(default=None) # Override the intercept parameter for the threshold function in the jump detection algorithm.
        threshold_constant = float(default=None) # Override the constant parameter for the threshold function in the jump detection algorithm.
    """

    weighting = "optimal"  # Only weighting allowed for OLS

    reference_file_types: ClassVar = ["readnoise", "gain", "dark"]

    def process(self, input):
        with rdd.open(input, mode="rw") as input_model:
            # Retrieve reference info
            readnoise_filename = self.get_reference_file(input_model, "readnoise")
            gain_filename = self.get_reference_file(input_model, "gain")
            log.info("Using READNOISE reference file: %s", readnoise_filename)
            readnoise_model = rdd.open(readnoise_filename, mode="r")
            log.info("Using GAIN reference file: %s", gain_filename)
            gain_model = rdd.open(gain_filename, mode="r")

            # Do the fitting.
            algorithm = self.algorithm.lower()
            if algorithm == "ols_cas22":
                dark_filename = self.get_reference_file(input_model, "dark")
                dark_model = rdd.open(dark_filename, mode="r")
                out_model = self.ols_cas22(
                    input_model, readnoise_model, gain_model, dark_model
                )
                out_model.meta.cal_step.ramp_fit = "COMPLETE"
            else:
                log.error("Algorithm %s is not supported. Skipping step.")
                out_model = input
                out_model.meta.cal_step.ramp_fit = "SKIPPED"

        return out_model

    def ols_cas22(self, input_model, readnoise_model, gain_model, dark_model):
        """Peform Optimal Linear Fitting on arbitrarily space resulants

        Parameters
        ----------
        input_model : RampModel
            Model containing ramps.

        readnoise_model : ReadnoiseRefModel
            Model with the read noise reference information.

        gain_model : GainRefModel
            Model with the gain reference information.

        dark_model : DarkRefModel
            Model with the dark reference information

        Returns
        -------
        out_model : ImageModel
            Model containing a count-rate image.
        """
        use_jump = self.use_ramp_jump_detection

        if use_jump:
            log.info("Jump detection as part of ramp fitting is enabled.")
        else:
            log.info("Jump detection as part of ramp fitting is disabled.")

        kwargs = {}
        if self.threshold_intercept is not None:
            kwargs["threshold_intercept"] = self.threshold_intercept
        if self.threshold_constant is not None:
            kwargs["threshold_constant"] = self.threshold_constant

        resultants = input_model.data
        dq = input_model.groupdq
        read_noise = readnoise_model.data
        gain = gain_model.data
        read_time = input_model.meta.exposure.frame_time

        # Force read pattern to be pure lists not LNodes
        read_pattern = [list(reads) for reads in input_model.meta.exposure.read_pattern]
        if len(read_pattern) != resultants.shape[0]:
            raise RuntimeError("mismatch between resultants shape and read_pattern.")

        # add dark current back into resultants so that Poisson noise is
        # properly accounted for
        tbar = np.array([np.mean(reads) * read_time for reads in read_pattern])
        resultants += dark_model.dark_slope[None, ...] * tbar[:, None, None]

        # account for the gain
        resultants *= gain
        read_noise *= gain

        # Fit the ramps
        output = ols_cas22_fit.fit_ramps_casertano(
            resultants,
            dq,
            read_noise,
            read_time,
            read_pattern=read_pattern,
            use_jump=use_jump,
            **kwargs,
        )

        # Break out the information and fix units
        slopes = output.parameters[..., Parameter.slope]
        var_rnoise = output.variances[..., Variance.read_var]
        var_poisson = output.variances[..., Variance.poisson_var]
        err = np.sqrt(var_poisson + var_rnoise)
        dq = output.dq.astype(np.uint32)

        # remove dark current contribution to slopes
        slopes -= dark_model.dark_slope * gain

        # Propagate DQ flags forward.
        ramp_dq = get_pixeldq_flags(dq, input_model.pixeldq, slopes, err, gain)

        # Create the image model
        image_info = (slopes, ramp_dq, var_poisson, var_rnoise, err)
        image_model = create_image_model(input_model, image_info)

        # Rescale by the gain back to DN/s
        image_model.data /= gain[4:-4, 4:-4]
        image_model.err /= gain[4:-4, 4:-4]
        image_model.var_poisson /= gain[4:-4, 4:-4] ** 2
        image_model.var_rnoise /= gain[4:-4, 4:-4] ** 2

        # That's all folks
        return image_model


# #########
# Utilities
# #########
def create_image_model(input_model, image_info):
    """
    Creates an ImageModel from the computed arrays from ramp_fit.

    Parameters
    ----------
    input_model : `~roman_datamodels.datamodels.RampModel`
        Input ``RampModel`` for which the output ``ImageModel`` is created.

    image_info : tuple
        The ramp fitting arrays needed for the ``ImageModel``.

    Returns
    -------
    out_model : `~roman_datamodels.datamodels.ImageModel`
        The output ``ImageModel`` to be returned from the ramp fit step.
    """
    data, dq, var_poisson, var_rnoise, err = image_info

    if dq is None:
        dq = np.zeros(data.shape, dtype="u4")

    # Create output datamodel
    # ... and add all keys from input
    meta = dict(wcs=None)  # default empty WCS
    meta.update(input_model.meta)
    meta["cal_step"]["ramp_fit"] = "INCOMPLETE"
    meta["cal_logs"] = rds.CalLogs()
    meta["photometry"] = rds.Photometry(
        {
            "conversion_megajanskys": None,
            "conversion_megajanskys_uncertainty": None,
            "pixel_area": None,
        }
    )
    inst = {
        "meta": meta,
        "data": data,
        "dq": dq,
        "var_poisson": var_poisson,
        "var_rnoise": var_rnoise,
        "err": err,
        "amp33": input_model.amp33.copy(),
        "border_ref_pix_left": input_model.border_ref_pix_left.copy(),
        "border_ref_pix_right": input_model.border_ref_pix_right.copy(),
        "border_ref_pix_top": input_model.border_ref_pix_top.copy(),
        "border_ref_pix_bottom": input_model.border_ref_pix_bottom.copy(),
        "dq_border_ref_pix_left": input_model.dq_border_ref_pix_left.copy(),
        "dq_border_ref_pix_right": input_model.dq_border_ref_pix_right.copy(),
        "dq_border_ref_pix_top": input_model.dq_border_ref_pix_top.copy(),
        "dq_border_ref_pix_bottom": input_model.dq_border_ref_pix_bottom.copy(),
    }
    out_node = rds.WfiImage(inst)
    im = rdd.ImageModel(out_node)

    # trim off border reference pixels from science data, dq, err
    # and var_poisson/var_rnoise
    im.data = im.data[4:-4, 4:-4]
    im.dq = im.dq[4:-4, 4:-4]
    im.err = im.err[4:-4, 4:-4]
    im.var_poisson = im.var_poisson[4:-4, 4:-4]
    im.var_rnoise = im.var_rnoise[4:-4, 4:-4]

    return im


def get_pixeldq_flags(groupdq, pixeldq, slopes, err, gain):
    """
    Construct pixeldq for ramp fit output from input dqs and ramp slopes.

    The algorithm is:
    - pass forward existing pixeldq flags
    - if we flagged a jump, flag the pixel as containing a jump
    - if everything is saturated, flag the pixel as saturated
    - if everything is saturated or do not use, flag the pixel as do not use
    - add NO_GAIN_VALUE if gain is not finite or less than zero

    Parameters
    ----------
    groupdq : np.ndarray
        dq flags for each resultant
    pixeldq : np.ndarray
        dq flags for each pixel
    slopes : np.ndarray
        derived slopes for each pixel
    err : np.ndarray
        derived total uncertainty for each pixel
    gain : np.ndarray
        gains for each pixel

    Returns
    -------
    pixeldq : np.ndarray
        Updated pixeldq array combining information from input dq and slopes.
    """
    outpixeldq = pixeldq.copy()
    # jump flagging
    m = np.any(groupdq & group.JUMP_DET, axis=0)
    outpixeldq |= (m * pixel.JUMP_DET).astype(np.uint32)
    # all saturated flagging
    m = np.all(groupdq & group.SATURATED, axis=0)
    outpixeldq |= (m * pixel.SATURATED).astype(np.uint32)
    # all either saturated or do not use or NaN slope flagging
    satordnu = group.SATURATED | group.DO_NOT_USE
    m = np.all(groupdq & satordnu, axis=0)
    m |= ~np.isfinite(slopes) | (err <= 0)
    outpixeldq |= (m * pixel.DO_NOT_USE).astype(np.uint32)
    m = (gain < 0) | ~np.isfinite(gain)
    outpixeldq |= (m * pixel.NO_GAIN_VALUE).astype(np.uint32)

    return outpixeldq
