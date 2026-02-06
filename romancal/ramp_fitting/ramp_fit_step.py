#! /usr/bin/env python
#
from __future__ import annotations

import copy
import logging
from typing import TYPE_CHECKING

import asdf
import numpy as np
from roman_datamodels import datamodels as rdm
from roman_datamodels.dqflags import group, pixel
from stcal.jump.jump_class import JumpData
from stcal.ramp_fitting import ols_cas22_fit
from stcal.ramp_fitting.likely_fit import likely_ramp_fit
from stcal.ramp_fitting.ols_cas22 import Parameter, Variance

from romancal.datamodels.fileio import open_dataset
from romancal.stpipe import RomanStep

SQRT2 = np.sqrt(2)

if TYPE_CHECKING:
    from typing import ClassVar

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

__all__ = ["RampFitStep"]


class RampFitStep(RomanStep):
    """This step fits a straight line to the value of counts vs. time to
    determine the mean count rate for each pixel.
    """

    class_alias = "ramp_fit"

    spec = """
        algorithm = option('ols_cas22', 'likely', default='likely')  # Algorithm to use to fit.
        suffix = string(default='rampfit')  # Default suffix of results
        use_ramp_jump_detection = boolean(default=True) # Use jump detection during ramp fitting
        threshold_intercept = float(default=None) # Override the intercept parameter for the threshold function in the jump detection algorithm.
        threshold_constant = float(default=None) # Override the constant parameter for the threshold function in the jump detection algorithm.
        include_var_rnoise = boolean(default=False) # include var_rnoise in output (can be reconstructed from err and other variances)
        maximum_cores = string(default='1') # cores for multiprocessing. Can be an integer, 'half', 'quarter', or 'all'
        expand_large_events = boolean(default=False) # Turns on Snowball detection
        min_sat_area = float(default=1.0) # minimum required area for the central saturation of snowballs
        min_jump_area = float(default=6.0) # minimum area to trigger large events processing
        expand_factor = float(default=1.9) # The expansion factor for the enclosing circles or ellipses
        sat_required_snowball = boolean(default=True) # Require the center of snowballs to be saturated
        min_sat_radius_extend = float(default=0.5) # The min radius of the sat core to trigger the extension of the core
        sat_expand = integer(default=2) # Number of pixels to add to the radius of the saturated core of snowballs
        edge_size = integer(default=0) # Distance from detector edge where a saturated core is not required for snowball detection
     """

    weighting = "optimal"  # Only weighting allowed for OLS

    reference_file_types: ClassVar = ["readnoise", "gain"]

    def process(self, dataset):
        input_model = open_dataset(dataset, update_version=self.update_version)

        # Retrieve reference info
        readnoise_filename = self.get_reference_file(input_model, "readnoise")
        gain_filename = self.get_reference_file(input_model, "gain")
        log.info("Using READNOISE reference file: %s", readnoise_filename)
        readnoise_model = rdm.open(readnoise_filename, mode="r")
        log.info("Using GAIN reference file: %s", gain_filename)
        gain_model = rdm.open(gain_filename, mode="r")

        # Do the fitting based on the algorithm selected.
        algorithm = self.algorithm.lower()
        if algorithm == "ols_cas22":
            out_model = self.ols_cas22(
                input_model,
                readnoise_model,
                gain_model,
                include_var_rnoise=self.include_var_rnoise,
            )
            out_model.meta.cal_step.ramp_fit = "COMPLETE"
        elif algorithm == "likely":
            out_model = self.likely(
                input_model,
                readnoise_model,
                gain_model,
            )
        else:
            log.error("Algorithm %s is not supported. Skipping step.")
            out_model = input
            out_model.meta.cal_step.ramp_fit = "SKIPPED"

        return out_model

    def ols_cas22(
        self, input_model, readnoise_model, gain_model, include_var_rnoise=False
    ):
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

        gain = gain_model.data
        read_time = input_model.meta.exposure.frame_time

        # Account for the gain.  Do not modify input_model.
        resultants = input_model.data * gain
        read_noise = readnoise_model.data * gain
        dq = input_model.groupdq

        # Force read pattern to be pure lists not LNodes
        read_pattern = [list(reads) for reads in input_model.meta.exposure.read_pattern]
        if len(read_pattern) != resultants.shape[0]:
            raise RuntimeError("mismatch between resultants shape and read_pattern.")

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

        # Break out the information and fix units back to DN/s
        slopes = output.parameters[..., Parameter.slope] / gain
        var_rnoise = output.variances[..., Variance.read_var] / gain**2
        var_poisson = output.variances[..., Variance.poisson_var] / gain**2
        err = np.sqrt(var_poisson + var_rnoise)
        dq = output.dq.astype(np.uint32)

        # Propagate DQ flags forward.
        ramp_dq = get_pixeldq_flags(dq, input_model.pixeldq, slopes, err, gain)

        # Create the image model.  Rescale by the gain back to DN/s
        image_info = {
            "slope": slopes,
            "dq": ramp_dq,
            "var_poisson": var_poisson,
            "var_rnoise": var_rnoise,
            "err": err,
        }
        image_model = create_image_model(
            input_model, image_info, include_var_rnoise=include_var_rnoise
        )

        # That's all folks
        return image_model

    def likely(self, input_model, readnoise_model, gain_model):
        """Perform Maximum Likelihood Algorithm

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
        # Add the needed components to the input model.
        input_model["flags_do_not_use"] = pixel.DO_NOT_USE
        input_model["flags_saturated"] = pixel.SATURATED
        input_model["rejection_threshold"] = None
        input_model["flags_jump_det"] = pixel.JUMP_DET
        # Add an axis to match the JWST data cube
        input_model.data = input_model.data[np.newaxis, :, :, :]
        input_model.groupdq = input_model.groupdq[np.newaxis, :, :, :]
        # add ancillary information needed by likelihood fitting
        input_model.read_pattern = get_readtimes(input_model)
        input_model.zeroframe = None
        input_model.average_dark_current = np.zeros(
            [input_model.data.shape[2], input_model.data.shape[3]]
        )

        # Flag pixels that have only a single resultant.
        input_model.pixeldq |= np.sum(
            input_model.groupdq[0] & (group.SATURATED | group.DO_NOT_USE), axis=0
        )

        # Setup jump data to handle snowballs and other special situations handled by
        # the likelihood algorithm
        jump_data = self._setup_jump_data(input_model, readnoise_model, gain_model)

        image_info, _, _ = likely_ramp_fit(
            input_model, readnoise_model.data, gain_model.data, jump_data=jump_data
        )

        out_model = create_image_model(
            input_model, image_info, include_var_rnoise=self.include_var_rnoise
        )
        out_model.meta.cal_step.ramp_fit = "COMPLETE"

        return out_model

    def _setup_jump_data(self, result, rnoise_m, gain_m):
        """
        Create a JumpData instance to be used by STCAL jump.

        Parameters
        ----------
        result : RampModel
            The ramp model input from the previous step.

        rnoise_m : ReadNoise model
            Readnoise reference model

        gain_m : GainModel
            Gain reference model

        Returns
        -------
        jump_data : JumpData
            The data container to be used to run the STCAL detect_jumps_data.
        """

        # Instantiate a JumpData class and populate it based on the input RampModel.
        jump_data = JumpData(result, gain_m.data, rnoise_m.data, pixel)

        # Setup snowball detection
        sat_expand = self.sat_expand * 2
        jump_data.set_snowball_info(
            self.expand_large_events,
            self.min_jump_area,
            self.min_sat_area,
            self.expand_factor,
            self.sat_required_snowball,
            self.min_sat_radius_extend,
            sat_expand,
            self.edge_size,
        )

        # Performance setup
        jump_data.max_cores = self.maximum_cores

        return jump_data


# #########
# Utilities
# #########
def create_image_model(input_model, image_info, include_var_rnoise=False):
    """Creates an ImageModel from the computed arrays from ramp_fit.

    Parameters
    ----------
    input_model : `~roman_datamodels.datamodels.RampModel`
        Input ``RampModel`` for which the output ``ImageModel`` is created.

    image_info : dict
        The ramp fitting arrays needed for the ``ImageModel``.

    include_var_rnoise : bool
        If True, include var_rnoise in the output model.

    Returns
    -------
    out_model : `~roman_datamodels.datamodels.ImageModel`
        The output ``ImageModel`` to be returned from the ramp fit step.

    """
    im = rdm.ImageModel()
    # use getitem here to avoid copying the DNode
    im.meta = copy.deepcopy(input_model["meta"])
    im.meta.model_type = "ImageModel"
    # since we've copied nodes let's remove any "read" tags
    for node in asdf.treeutil.iter_tree(im):
        if node._read_tag is not None:
            node._read_tag = None
    im.meta.product_type = "l2"
    im.meta.cal_step = {}
    for step_name in im.schema_info("required")["roman"]["meta"]["cal_step"][
        "required"
    ].info:
        im.meta.cal_step[step_name] = input_model.meta.cal_step.get(
            step_name, "INCOMPLETE"
        )
    im.meta.cal_logs = []
    im.meta.photometry = {
        "conversion_megajanskys": -999999,
        "conversion_megajanskys_uncertainty": -999999,
        "pixel_area": -999999,
    }
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
    im.data = image_info["slope"][4:-4, 4:-4].copy()
    if image_info["dq"] is not None:
        im.dq = image_info["dq"][4:-4, 4:-4].copy()
    else:
        im.dq = np.zeros(im.data.shape, dtype="u4")

    im.err = image_info["err"][4:-4, 4:-4].copy().astype("float16")
    im.var_poisson = image_info["var_poisson"][4:-4, 4:-4].copy().astype("float16")
    if include_var_rnoise:
        im.var_rnoise = image_info["var_rnoise"][4:-4, 4:-4].copy().astype("float16")

    # Add required chisq and dumo fields.  chisq will be populated with zeroes
    # if the likelihood algorithm was not run, so that "chisq" is not a key in
    # image_info.

    if "chisq" in image_info.keys() and image_info["chisq"] is not None:
        im.chisq = image_info["chisq"][4:-4, 4:-4].copy().astype("float16")
    else:
        im.chisq = np.zeros(im.data.shape, dtype=np.float16)

    slopes_alt = slopes_uniform_weights(input_model)

    # Add this to the optimal-weighted slopes to get the uniform-weighted slopes
    im.dumo = (slopes_alt[4:-4, 4:-4] - im.data).astype("float16")

    return im


def slopes_uniform_weights(input_model):
    """
    Compute ramp slopes using uniform (read-noise-limited) weights.

    Parameters
    ----------
    input_model : RampModel
        Model containing ramps.

    Returns
    -------
    slopes : ndarray
        The slope for each pixel under uniform weighting, which is optimal
        in the read noise limit.  All flags, including saturation and
        jump, will be ignored.
    """

    # The lines below compute the weight for each resultant in the case
    # of uniform weighting (a diagonal covariance matrix consisting only
    # of read noise).

    readtimes = get_readtimes(input_model)

    ni = np.array([len(t) for t in readtimes])
    ti = np.array([np.mean(t) for t in readtimes])
    N = np.sum(ni)
    Nt = np.sum(ni * ti)
    Ntt = np.sum(ni * ti**2)
    weights = (N * ni * ti - Nt * ni) / (N * Ntt - Nt**2)

    # We want the weighted sum over reads.
    if len(input_model.data.shape) == 3:
        return np.sum(weights[:, None, None] * input_model.data, axis=0)
    elif len(input_model.data.shape) == 4:
        return np.sum(weights[:, None, None] * input_model.data[0], axis=0)
    else:
        raise ValueError("Unexpected shape for input_model.data")


def get_pixeldq_flags(groupdq, pixeldq, slopes, err, gain):
    """Construct pixeldq for ramp fit output from input dqs and ramp slopes.

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
    """Get the read times needed to compute the covariance matrices.

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

    rtimes = [
        list(np.array(r) * ramp_data.meta.exposure.frame_time)
        for r in ramp_data.meta.exposure.read_pattern
    ]

    return rtimes
