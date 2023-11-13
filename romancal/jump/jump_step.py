"""
Detect jumps in a science image
"""

import logging
import time

import numpy as np
from roman_datamodels import datamodels as rdd
from stcal.jump.jump import detect_jumps

from romancal.lib import dqflags
from romancal.stpipe import RomanStep

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


__all__ = ["JumpStep"]


class JumpStep(RomanStep):
    """
    JumpStep: Performs CR/jump detection. The 2-point difference method is
    applied.
    """

    spec = """
        rejection_threshold = float(default=180.0,min=0) # CR sigma rej thresh
        three_group_rejection_threshold = float(default=185.0,min=0) # CR sigma rej thresh
        four_group_rejection_threshold = float(default=190.0,min=0) # CR sigma rej thresh
        maximum_cores = option('none', 'quarter', 'half', 'all', default='none') # max number of processes to create
        flag_4_neighbors = boolean(default=True) # flag the four perpendicular neighbors of each CR
        max_jump_to_flag_neighbors = float(default=1000) # maximum jump sigma that will trigger neighbor flagging
        min_jump_to_flag_neighbors = float(default=10) # minimum jump sigma that will trigger neighbor flagging
        min_sat_area = float(default=1.0) # minimum required area for the central saturation of snowballs
        min_jump_area = float(default=5.0) # minimum area to trigger large events processing
        expand_factor = float(default=2.0) # The expansion factor for the enclosing circles or ellipses
        use_ellipses = boolean(default=False) # Use an enclosing ellipse rather than a circle for MIRI showers
        sat_required_snowball = boolean(default=True) # Require the center of snowballs to be saturated
        expand_large_events = boolean(default=False) # must be True to trigger snowball and shower flagging
        use_ramp_jump_detection = boolean(default=True) # Use jump detection during ramp fitting
    """  # noqa: E501

    reference_file_types = ["gain", "readnoise"]

    def process(self, input):
        # Open input as a Roman DataModel (single integration; 3D arrays)

        with rdd.open(input, lazy_load=False) as input_model:
            # Extract the needed info from the Roman Data Model
            meta = input_model.meta
            r_data = input_model.data.value
            r_gdq = input_model.groupdq
            r_pdq = input_model.pixeldq
            r_err = input_model.err.value
            result = input_model

            # If the ramp fitting jump detection is enabled, then skip this step
            if self.use_ramp_jump_detection:
                result.meta.cal_step.jump = "SKIPPED"
                return result

            frames_per_group = meta.exposure.nframes

            # Modify the arrays for input into the 'common' jump (4D)
            data = np.copy(r_data[np.newaxis, :])
            gdq = r_gdq[np.newaxis, :]
            pdq = r_pdq[np.newaxis, :]
            err = np.copy(r_err[np.newaxis, :])

            tstart = time.time()

            # Check for an input model with NGROUPS<=2
            ngroups = data.shape[1]

            if ngroups <= 2:
                self.log.warning("Cannot apply jump detection as NGROUPS<=2;")
                self.log.warning("Jump step will be skipped")

                result = input_model

                result.meta.cal_step.jump = "SKIPPED"
                return result

            # Retrieve the parameter values
            rej_thresh = self.rejection_threshold
            three_grp_rej_thresh = self.three_group_rejection_threshold
            four_grp_rej_thresh = self.four_group_rejection_threshold
            max_cores = self.maximum_cores
            max_jump_to_flag_neighbors = self.max_jump_to_flag_neighbors
            min_jump_to_flag_neighbors = self.min_jump_to_flag_neighbors
            flag_4_neighbors = self.flag_4_neighbors
            min_sat_area = self.min_sat_area
            min_jump_area = self.min_jump_area
            expand_factor = self.expand_factor
            use_ellipses = self.use_ellipses
            sat_required_snowball = self.sat_required_snowball
            expand_large_events = self.expand_large_events

            self.log.info("CR rejection threshold = %g sigma", rej_thresh)
            if self.maximum_cores != "none":
                self.log.info("Maximum cores to use = %s", max_cores)

            # Get the gain and readnoise reference files
            gain_filename = self.get_reference_file(input_model, "gain")
            self.log.info("Using GAIN reference file: %s", gain_filename)
            gain_model = rdd.GainRefModel(gain_filename)
            gain_2d = gain_model.data.value

            readnoise_filename = self.get_reference_file(input_model, "readnoise")
            self.log.info("Using READNOISE reference file: %s", readnoise_filename)
            readnoise_model = rdd.ReadnoiseRefModel(readnoise_filename)
            # This is to clear the WRITEABLE=False flag?
            readnoise_2d = np.copy(readnoise_model.data.value)

            # DG 0810/21:  leave for now; make dqflags changes in a later,
            #              separate PR
            dqflags_d = {}  # Dict of DQ flags
            dqflags_d = {
                "GOOD": dqflags.group["GOOD"],
                "DO_NOT_USE": dqflags.group["DO_NOT_USE"],
                "SATURATED": dqflags.group["SATURATED"],
                "JUMP_DET": dqflags.group["JUMP_DET"],
                "NO_GAIN_VALUE": dqflags.pixel["NO_GAIN_VALUE"],
            }

            gdq, pdq, *_ = detect_jumps(
                frames_per_group,
                data,
                gdq,
                pdq,
                err,
                gain_2d,
                readnoise_2d,
                rej_thresh,
                three_grp_rej_thresh,
                four_grp_rej_thresh,
                max_cores,
                max_jump_to_flag_neighbors,
                min_jump_to_flag_neighbors,
                flag_4_neighbors,
                dqflags_d,
                min_sat_area=min_sat_area,
                min_jump_area=min_jump_area,
                expand_factor=expand_factor,
                use_ellipses=use_ellipses,
                sat_required_snowball=sat_required_snowball,
                expand_large_events=expand_large_events,
            )

            gdq = gdq[0, :, :, :]
            pdq = pdq[0, :, :]
            result.groupdq = gdq
            result.pixeldq = pdq
            gain_model.close()
            readnoise_model.close()
            tstop = time.time()
            self.log.info("The execution time in seconds: %f", tstop - tstart)

        result.meta.cal_step.jump = "COMPLETE"

        if self.save_results:
            try:
                self.suffix = "jump"
            except AttributeError:
                self["suffix"] = "jump"

        return result
