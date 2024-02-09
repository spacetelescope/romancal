#! /usr/bin/env python
#
import logging

import numpy as np
from astropy import units as u
from roman_datamodels import datamodels as rdd
from roman_datamodels import maker_utils
from roman_datamodels import stnode as rds
from stcal.ramp_fitting import ols_cas22_fit, ramp_fit
from stcal.ramp_fitting.ols_cas22 import Parameter, Variance

from romancal.lib import dqflags
from romancal.stpipe import RomanStep

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

__all__ = ["RampFitStep"]


class RampFitStep(RomanStep):
    """
    This step fits a straight line to the value of counts vs. time to
    determine the mean count rate for each pixel.
    """

    spec = """
        algorithm = option('ols','ols_cas22', default='ols_cas22')  # Algorithm to use to fit.
        save_opt = boolean(default=False) # Save optional output
        opt_name = string(default='')
        maximum_cores = option('none','quarter','half','all',default='none') # max number of processes to create
        suffix = string(default='rampfit')  # Default suffix of results
        use_ramp_jump_detection = boolean(default=True) # Use jump detection during ramp fitting
        threshold_intercept = float(default=None) # Override the intercept parameter for the threshold function in the jump detection algorithm.
        threshold_constant = float(default=None) # Override the constant parameter for the threshold function in the jump detection algorithm.
    """  # noqa: E501

    weighting = "optimal"  # Only weighting allowed for OLS

    reference_file_types = ["readnoise", "gain"]

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
            if algorithm == "ols":
                out_model = self.ols(input_model, readnoise_model, gain_model)
                out_model.meta.cal_step.ramp_fit = "COMPLETE"
            elif algorithm == "ols_cas22":
                out_model = self.ols_cas22(input_model, readnoise_model, gain_model)
                out_model.meta.cal_step.ramp_fit = "COMPLETE"
                if self.use_ramp_jump_detection:
                    out_model.meta.cal_step.jump = "COMPLETE"
            else:
                log.error("Algorithm %s is not supported. Skipping step.")
                out_model = input
                out_model.meta.cal_step.ramp_fit = "SKIPPED"

        return out_model

    def ols(self, input_model, readnoise_model, gain_model):
        """Perform Optimal Linear Fitting on evenly-spaced resultants

        The OLS algorithm used is the same used by JWST for it's ramp fitting.

        Parameters
        ----------
        input_model : RampModel
            Model containing ramps.

        readnoise_model : ReadnoiseRefModel
            Model with the read noise reference information.

        gain_model : GainRefModel
            Model with the gain reference information.

        Returns
        -------
        out_model : ImageModel
            Model containing a count-rate image.
        """
        max_cores = self.maximum_cores
        input_model.data = input_model.data[np.newaxis, :]
        input_model.groupdq = input_model.groupdq[np.newaxis, :]
        input_model.err = input_model.err[np.newaxis, :]

        log.info(f"Using algorithm = {self.algorithm}")
        log.info(f"Using weighting = {self.weighting}")

        buffsize = ramp_fit.BUFSIZE
        image_info, integ_info, opt_info, gls_opt_model = ramp_fit.ramp_fit(
            input_model,
            buffsize,
            self.save_opt,
            readnoise_model.data.value,
            gain_model.data.value,
            self.algorithm,
            self.weighting,
            max_cores,
            dqflags.pixel,
        )

        if image_info is not None:
            # JWST has a separate GainScale step.  This isn't in Roman.
            # We can adjust the gains here instead.
            # data, dq, var_poisson, var_rnoise, err = image_info
            # This step isn't needed for ols_cas22, which works directly in
            # electrons.
            image_info[0][...] *= gain_model.data.value  # data
            # no dq multiplication!
            image_info[2][...] *= gain_model.data.value**2  # var_poisson
            image_info[3][...] *= gain_model.data.value**2  # var_rnoise
            image_info[4][...] *= gain_model.data.value  # err

        readnoise_model.close()

        # Save the OLS optional fit product, if it exists
        if opt_info is not None:
            # We should scale these by the gain, too.
            # opt_info = (slope, sigslope, var_poisson, var_readnoise, yint,
            # sig_yint,
            # ped_int, weights, cr_mag_seg
            opt_info[0][...] *= gain_model.data.value  # slope
            opt_info[1][...] *= gain_model.data.value  # sigslope
            opt_info[2][...] *= gain_model.data.value**2  # var_poisson
            opt_info[3][...] *= gain_model.data.value**2  # var_readnoise
            opt_info[4][...] *= gain_model.data.value  # yint
            opt_info[5][...] *= gain_model.data.value  # sigyint
            opt_info[6][...] *= gain_model.data.value  # pedestal
            # no change to weights due to gain
            opt_info[8][...] *= gain_model.data.value  # crmag
            opt_model = create_optional_results_model(input_model, opt_info)
            self.save_model(opt_model, "fitopt", output_file=self.opt_name)

        gain_model.close()

        # All pixels saturated, therefore returning an image file with zero data
        if image_info is None:
            log.info("All pixels are saturated. Returning a zeroed-out image.")

            # Image info order is: data, dq, var_poisson, var_rnoise, err
            image_info = (
                np.zeros(input_model.data.shape[2:], dtype=input_model.data.dtype),
                input_model.pixeldq
                | input_model.groupdq[0][0]
                | dqflags.group["SATURATED"],
                np.zeros(input_model.err.shape[2:], dtype=input_model.err.dtype),
                np.zeros(input_model.err.shape[2:], dtype=input_model.err.dtype),
                np.zeros(input_model.err.shape[2:], dtype=input_model.err.dtype),
            )

        out_model = create_image_model(input_model, image_info)
        return out_model

    def ols_cas22(self, input_model, readnoise_model, gain_model):
        """Peform Optimal Linear Fitting on arbitrarily space resulants

        Parameters
        ----------
        input_model : RampModel
            Model containing ramps.

        readnoise_model : ReadnoiseRefModel
            Model with the read noise reference information.

        gain_model : GainRefModel
            Model with the gain reference information.

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

        resultants = input_model.data.value
        dq = input_model.groupdq
        read_noise = readnoise_model.data.value
        gain = gain_model.data.value
        read_time = input_model.meta.exposure.frame_time

        # account for the gain
        resultants *= gain
        read_noise *= gain

        # Force read pattern to be pure lists not LNodes
        read_pattern = [list(reads) for reads in input_model.meta.exposure.read_pattern]

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

        # Propagate DQ flags forward.
        ramp_dq = get_pixeldq_flags(dq, input_model.pixeldq, slopes, err, gain)

        # Create the image model
        image_info = (slopes, ramp_dq, var_poisson, var_rnoise, err)
        image_model = create_image_model(input_model, image_info)

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

    data = u.Quantity(data, u.electron / u.s, dtype=data.dtype)
    var_poisson = u.Quantity(
        var_poisson, u.electron**2 / u.s**2, dtype=var_poisson.dtype
    )
    var_rnoise = u.Quantity(var_rnoise, u.electron**2 / u.s**2, dtype=var_rnoise.dtype)
    err = u.Quantity(err, u.electron / u.s, dtype=err.dtype)
    if dq is None:
        dq = np.zeros(data.shape, dtype="u4")

    # Create output datamodel
    # ... and add all keys from input
    meta = dict(wcs=None)  # default empty WCS
    meta.update(input_model.meta)
    meta["cal_step"]["ramp_fit"] = "INCOMPLETE"
    meta["photometry"] = maker_utils.mk_photometry()
    inst = {
        "meta": meta,
        "data": u.Quantity(data, u.electron / u.s, dtype=data.dtype),
        "dq": dq,
        "var_poisson": u.Quantity(
            var_poisson, u.electron**2 / u.s**2, dtype=var_poisson.dtype
        ),
        "var_rnoise": u.Quantity(
            var_rnoise, u.electron**2 / u.s**2, dtype=var_rnoise.dtype
        ),
        "err": u.Quantity(err, u.electron / u.s, dtype=err.dtype),
        "amp33": input_model.amp33.copy(),
        "border_ref_pix_left": input_model.border_ref_pix_left.copy(),
        "border_ref_pix_right": input_model.border_ref_pix_right.copy(),
        "border_ref_pix_top": input_model.border_ref_pix_top.copy(),
        "border_ref_pix_bottom": input_model.border_ref_pix_bottom.copy(),
        "dq_border_ref_pix_left": input_model.dq_border_ref_pix_left.copy(),
        "dq_border_ref_pix_right": input_model.dq_border_ref_pix_right.copy(),
        "dq_border_ref_pix_top": input_model.dq_border_ref_pix_top.copy(),
        "dq_border_ref_pix_bottom": input_model.dq_border_ref_pix_bottom.copy(),
        "cal_logs": rds.CalLogs(),
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


def create_optional_results_model(input_model, opt_info):
    """
    Creates the optional output from the computed arrays from ramp_fit.

    Parameters
    ----------
    input_model : `~roman_datamodels.datamodels.RampModel`
        The input data model.
    opt_info : tuple
        The ramp fitting arrays needed for the ``RampFitOutputModel``.

    Returns
    -------
    opt_model : `~roman_datamodels.datamodels.RampFitOutputModel`
        The optional ``RampFitOutputModel`` to be returned from the ramp fit step.
    """
    (
        slope,
        sigslope,
        var_poisson,
        var_rnoise,
        yint,
        sigyint,
        pedestal,
        weights,
        crmag,
    ) = opt_info
    meta = {}
    meta.update(input_model.meta)
    crmag.shape = crmag.shape[1:]
    crmag.dtype = np.float32

    inst = {
        "meta": meta,
        "slope": u.Quantity(np.squeeze(slope), u.electron / u.s, dtype=slope.dtype),
        "sigslope": u.Quantity(
            np.squeeze(sigslope), u.electron / u.s, dtype=sigslope.dtype
        ),
        "var_poisson": u.Quantity(
            np.squeeze(var_poisson), u.electron**2 / u.s**2, dtype=var_poisson.dtype
        ),
        "var_rnoise": u.Quantity(
            np.squeeze(var_rnoise), u.electron**2 / u.s**2, dtype=var_rnoise.dtype
        ),
        "yint": u.Quantity(np.squeeze(yint), u.electron, dtype=yint.dtype),
        "sigyint": u.Quantity(np.squeeze(sigyint), u.electron, dtype=sigyint.dtype),
        "pedestal": u.Quantity(np.squeeze(pedestal), u.electron, dtype=pedestal.dtype),
        "weights": np.squeeze(weights),
        "crmag": u.Quantity(crmag, u.electron, dtype=pedestal.dtype),
    }

    out_node = rds.RampFitOutput(inst)
    opt_model = rdd.RampFitOutputModel(out_node)
    opt_model.meta.filename = input_model.meta.filename

    return opt_model


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
    m = np.any(groupdq & dqflags.group["JUMP_DET"], axis=0)
    outpixeldq |= (m * dqflags.pixel["JUMP_DET"]).astype(np.uint32)
    # all saturated flagging
    m = np.all(groupdq & dqflags.group["SATURATED"], axis=0)
    outpixeldq |= (m * dqflags.pixel["SATURATED"]).astype(np.uint32)
    # all either saturated or do not use or NaN slope flagging
    satordnu = dqflags.group["SATURATED"] | dqflags.group["DO_NOT_USE"]
    m = np.all(groupdq & satordnu, axis=0)
    m |= ~np.isfinite(slopes) | (err <= 0)
    outpixeldq |= (m * dqflags.pixel["DO_NOT_USE"]).astype(np.uint32)
    m = (gain < 0) | ~np.isfinite(gain)
    outpixeldq |= (m * dqflags.pixel["NO_GAIN_VALUE"]).astype(np.uint32)

    return outpixeldq
