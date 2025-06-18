#! /usr/bin/env python
#
from __future__ import annotations
import pdb

import copy
import logging
from typing import TYPE_CHECKING

import asdf
import numpy as np
from roman_datamodels import datamodels as rdm
from roman_datamodels import stnode
import scipy
from roman_datamodels.dqflags import group, pixel
from stcal.ramp_fitting import ols_cas22_fit, utils
from stcal.ramp_fitting.ols_cas22 import Parameter, Variance
from stcal.ramp_fitting.likely_fit import (get_readtimes, compute_image_info,
                                           mask_jumps, fit_ramps)
from stcal.ramp_fitting.likely_algo_classes import IntegInfo, RampResult, Covar

from romancal.stpipe import RomanStep

SQRT2 = np.sqrt(2)

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
        algorithm = option('ols_cas22', 'likely', default='ols_cas22')  # Algorithm to use to fit.
        save_opt = boolean(default=False) # Save optional output
        suffix = string(default='rampfit')  # Default suffix of results
        use_ramp_jump_detection = boolean(default=True) # Use jump detection during ramp fitting
        threshold_intercept = float(default=None) # Override the intercept parameter for the threshold function in the jump detection algorithm.
        threshold_constant = float(default=None) # Override the constant parameter for the threshold function in the jump detection algorithm.
    """

    weighting = "optimal"  # Only weighting allowed for OLS

    reference_file_types: ClassVar = ["readnoise", "gain", "dark"]

    def process(self, input):
        with rdm.open(input, mode="rw") as input_model:
            # Retrieve reference info
            readnoise_filename = self.get_reference_file(input_model, "readnoise")
            gain_filename = self.get_reference_file(input_model, "gain")
            log.info("Using READNOISE reference file: %s", readnoise_filename)
            readnoise_model = rdm.open(readnoise_filename, mode="r")
            log.info("Using GAIN reference file: %s", gain_filename)
            gain_model = rdm.open(gain_filename, mode="r")

            # Do the fitting.
            algorithm = self.algorithm.lower()
            if algorithm == "ols_cas22":
                dark_filename = self.get_reference_file(input_model, "dark")
                dark_model = rdm.open(dark_filename, mode="r")
                out_model = self.ols_cas22(
                    input_model, readnoise_model, gain_model, dark_model
                )
                out_model.meta.cal_step.ramp_fit = "COMPLETE"
            elif  algorithm == "likely":
                input_model['flags_do_not_use'] = pixel.DO_NOT_USE
                input_model['flags_saturated'] = pixel.SATURATED
                input_model['rejection_threshold'] = None
                input_model['flags_jump_det'] = pixel.JUMP_DET
                input_model.data = input_model.data[np.newaxis, :, :, :]
                out_model = self.likely_fit(
                    input_model, readnoise_model.data, gain_model.data)
                out_model.meta.cal_step.ramp_fit = "COMPLETE"               
            else:
                log.error("Algorithm %s is not supported. Skipping step.")
                out_model = input
                out_model.meta.cal_step.ramp_fit = "SKIPPED"

        return out_model

    def likely_fit(self,  ramp_model, readnoise_2d, gain_2d):
        """Ramp fiting using the likelyhood algorithm from stcal
        Parameters
        ----------
        ramp_data : RampData
            Input data necessary for computing ramp fitting.

        readnoise_2d : ndarray
            readnoise for all pixels

        gain_2d : ndarray
            gain for all pixels

        Returns
        -------
        image_info : tuple
            The tuple of computed ramp fitting arrays.

        """
        
        image_info = None

        nints, nresultants, nrows, ncols = ramp_model.data.shape


        readtimes = get_readtimes(ramp_model)

        covar = Covar(readtimes)
        integ_class = IntegInfo(nints, nrows, ncols)

        readnoise_2d = readnoise_2d / SQRT2

        for integ in range(nints):
            data = ramp_model.data[integ, :, :, :]
            gdq = ramp_model.groupdq.copy()
            pdq = ramp_model.pixeldq[:, :].copy()

            # Eqn (5)
            diff = (data[1:] - data[:-1]) / covar.delta_t[:, np.newaxis, np.newaxis]
            alldiffs2use = np.ones(diff.shape, np.uint8)

            for row in range(nrows):
                d2use = determine_diffs2use(ramp_model, integ, row, diff)
                d2use_copy = d2use.copy()  # Use to flag jumps
                #if row == 1024: pdb.set_trace()
                if ramp_model.rejection_threshold is not None:
                    threshold_one_omit = ramp_model.rejection_threshold**2
                    pval = scipy.special.erfc(ramp_model.rejection_threshold/SQRT2)
                    threshold_two_omit = scipy.stats.chi2.isf(pval, 2)
                    if np.isinf(threshold_two_omit):
                        threshold_two_omit = threshold_one_omit + 10
                    d2use, countrates = mask_jumps(
                        diff[:, row], covar, readnoise_2d[row], gain_2d[row],
                        threshold_one_omit=threshold_one_omit,
                        threshold_two_omit=threshold_two_omit,
                        diffs2use=d2use
                    )
                else:
                    d2use, countrates = mask_jumps(
                        diff[:, row], covar, readnoise_2d[row], gain_2d[row],
                        diffs2use=d2use
                    )

                # Set jump detection flags
                jump_locs = d2use_copy ^ d2use
                jump_locs[jump_locs > 0] = ramp_model.flags_jump_det
                gdq[1:, row] |= jump_locs

                alldiffs2use[:, row] = d2use

                #rateguess = countrates * (countrates > 0) + ramp_data.average_dark_current[row, :]
                rateguess = countrates * (countrates > 0)
                result = fit_ramps(
                    diff[:, row],
                    covar,
                    gain_2d[row],
                    readnoise_2d[row],
                    diffs2use=d2use,
                    count_rate_guess=rateguess,
                )
                integ_class.get_results(result, integ, row)

            pdq = utils.dq_compress_sect(ramp_model, integ, gdq, pdq)
            integ_class.dq[integ, :, :] = pdq

            del gdq

        image_info = compute_image_info(integ_class, ramp_model)

        image_model = create_image_model(ramp_model, image_info)
        # Rescale by the gain back to DN/s
        image_model.data *= gain_2d[4:-4, 4:-4]
        image_model.err /= gain_2d[4:-4, 4:-4]
        image_model.var_poisson /= gain_2d[4:-4, 4:-4] ** 2
        image_model.var_rnoise /= gain_2d[4:-4, 4:-4] ** 2
        
        return image_model


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
    im = rdm.ImageModel()
    # use getitem here to avoid copying the DNode
    im.meta = copy.deepcopy(input_model["meta"])
    im.meta.model_type = "ImageModel"
    # since we've copied nodes let's remove any "read" tags
    for node in asdf.treeutil.iter_tree(im):
        if hasattr(node, "_read_tag"):
            del node._read_tag
    im.meta.product_type = "l2"
    im.meta.cal_step = stnode.L2CalStep.create_minimal(input_model.meta.cal_step)
    im.meta.cal_logs = stnode.CalLogs.create_minimal()
    im.meta.photometry = stnode.Photometry.create_minimal(
        {
            "conversion_megajanskys": -999999,
            "conversion_megajanskys_uncertainty": -999999,
            "pixel_area": -999999,
        }
    )
    im.amp33 = input_model.amp33.copy()
    im.border_ref_pix_left = input_model.border_ref_pix_left.copy()
    im.border_ref_pix_right = input_model.border_ref_pix_right.copy()
    im.border_ref_pix_top = input_model.border_ref_pix_top.copy()
    im.border_ref_pix_bottom = input_model.border_ref_pix_bottom.copy()
    im.dq_border_ref_pix_left = input_model.dq_border_ref_pix_left.copy()
    im.dq_border_ref_pix_right = input_model.dq_border_ref_pix_right.copy()
    im.dq_border_ref_pix_top = input_model.dq_border_ref_pix_top.copy()
    im.dq_border_ref_pix_bottom = input_model.dq_border_ref_pix_bottom.copy()

    # trim off border reference pixels from science data, dq, err
    # and var_poisson/var_rnoise
    im.data = data[4:-4, 4:-4].copy()
    if dq is not None:
        im.dq = dq[4:-4, 4:-4].copy()
    else:
        im.dq = np.zeros(im.data.shape, dtype="u4")
    im.err = err[4:-4, 4:-4].copy()
    im.var_poisson = var_poisson[4:-4, 4:-4].copy()
    im.var_rnoise = var_rnoise[4:-4, 4:-4].copy()

    return im


def determine_diffs2use(ramp_data, integ, row, diffs):
    """
    Compute the diffs2use mask based on DQ flags of a row.

    Parameters
    ----------
    ramp_data : RampData
        Input data necessary for computig ramp fitting.

    integ : int
        The current integration being processed.

    row : int
        The current row being processed.

    diffs : ndarray
        The group differences of the data array for a given integration and row
        (ngroups-1, ncols).

    Returns
    -------
    d2use : ndarray
        A boolean array definined the segmented ramps for each pixel in a row.
        (ngroups-1, ncols)
    """
    # import ipdb; ipdb.set_trace()
    _, nresultants, _, ncols = ramp_data.data.shape
    dq = np.zeros(shape=(nresultants, ncols), dtype=np.uint8)
    dq[:, :] = ramp_data.groupdq[integ, row, :]
    d2use_tmp = np.ones(shape=diffs.shape, dtype=np.uint8)
    d2use = d2use_tmp[:, row]

    d2use[dq[1:, :] != 0] = 0
    d2use[dq[:-1, :] != 0] = 0

    return d2use


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

def get_readtimes(ramp_data):
    """
    Get the read times needed to compute the covariance matrices.

    If there is already a read_pattern in the ramp_data class, then just get it.
    If not, then one needs to be constructed.  If one needs to be constructed it
    is assumed the groups are evenly spaced in time, as are the frames that make
    up the group.  If each group has only one frame and no group gap, then a list
    of the group times is returned.  If nframes > 0, then a list of lists of each
    frame time in each group is returned with the assumption:
        group_time = (nframes + groupgap) * frame_time

    Parameters
    ----------
    ramp_data : RampData
        Input data necessary for computing ramp fitting.

    Returns
    -------
    readtimes : list
        A list of frame times for each frame used in the computation of the ramp.
    """
    nresultants = ramp_data.meta.exposure.nresultants
    log.info("Number of resultants: %d ", nresultants)

    ngroups = ramp_data.data.shape[1]
    #tot_frames = ramp_data.meta.exposure.nresultants + 0 #ramp_data.groupgap
    read_numbers = [x for xs in ramp_data.meta.exposure.read_pattern for x in xs]
    nskips = np.max(read_numbers) - len(read_numbers)
    tot_frames = len(read_numbers) + nskips
    tot_nreads = np.arange(1, ramp_data.meta.exposure.nresultants + 1)
    rtimes = [
        (tot_nreads + k * tot_frames) * ramp_data.meta.exposure.frame_time for k in range(ngroups)
    ]

    return rtimes
