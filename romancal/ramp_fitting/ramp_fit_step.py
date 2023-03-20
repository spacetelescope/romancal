#! /usr/bin/env python
#
import logging

import numpy as np
from astropy import units as u
from roman_datamodels import datamodels as rdd
from roman_datamodels import maker_utils
from roman_datamodels import stnode as rds
from stcal.ramp_fitting import ramp_fit

from romancal.lib import dqflags
from romancal.stpipe import RomanStep

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


__all__ = ["RampFitStep"]


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


def create_image_model(input_model, image_info):
    """
    Creates an ImageModel from the computed arrays from ramp_fit.

    Parameters
    ----------
    input_model : `~roman_datamodels.datamodels.RampModel`
        Input ``RampModel`` for which the output ``ImageModel`` is created.
    image_info : tuple
        The ramp fitting arrays needed for the ``ImageModel``.
    refpix_info : tuple
        The reference pixel arrays.

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
    var_rnoise = u.Quantity(
        var_rnoise, u.electron**2 / u.s**2, dtype=var_rnoise.dtype
    )
    err = u.Quantity(err, u.electron / u.s, dtype=err.dtype)

    # Create output datamodel
    # ... and add all keys from input
    meta = {}
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
        "amp33": input_model.amp33,
        "border_ref_pix_left": input_model.border_ref_pix_left,
        "border_ref_pix_right": input_model.border_ref_pix_right,
        "border_ref_pix_top": input_model.border_ref_pix_top,
        "border_ref_pix_bottom": input_model.border_ref_pix_bottom,
        "dq_border_ref_pix_left": input_model.dq_border_ref_pix_left,
        "dq_border_ref_pix_right": input_model.dq_border_ref_pix_right,
        "dq_border_ref_pix_top": input_model.dq_border_ref_pix_top,
        "dq_border_ref_pix_bottom": input_model.dq_border_ref_pix_bottom,
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


class RampFitStep(RomanStep):

    """
    This step fits a straight line to the value of counts vs. time to
    determine the mean count rate for each pixel.
    """

    spec = """
        opt_name = string(default='')
        maximum_cores = option('none','quarter','half','all',default='none') # max number of processes to create
        save_opt = boolean(default=False) # Save optional output
    """  # noqa: E501
    algorithm = "ols"  # Only algorithm allowed

    weighting = "optimal"  # Only weighting allowed

    reference_file_types = ["readnoise", "gain"]

    def process(self, input):
        with rdd.open(input, mode="rw") as input_model:
            max_cores = self.maximum_cores
            readnoise_filename = self.get_reference_file(input_model, "readnoise")
            gain_filename = self.get_reference_file(input_model, "gain")
            input_model.data = input_model.data[np.newaxis, :]
            input_model.groupdq = input_model.groupdq[np.newaxis, :]
            input_model.err = input_model.err[np.newaxis, :]

            log.info("Using READNOISE reference file: %s", readnoise_filename)
            readnoise_model = rdd.open(readnoise_filename, mode="rw")
            log.info("Using GAIN reference file: %s", gain_filename)
            gain_model = rdd.open(gain_filename, mode="rw")

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
            readnoise_model.close()
            gain_model.close()

        # Save the OLS optional fit product, if it exists
        if opt_info is not None:
            opt_model = create_optional_results_model(input_model, opt_info)
            self.save_model(opt_model, "fitopt", output_file=self.opt_name)

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
        out_model.meta.cal_step.ramp_fit = "COMPLETE"

        if self.save_results:
            try:
                self.suffix = "rampfit"
            except AttributeError:
                self["suffix"] = "rampfit"

        return out_model
