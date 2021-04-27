"""
Apply Jump-detect to a science image
"""
import logging
import time
import numpy as np

from ..stpipe import RomanStep
from .. import datamodels

from . jump import detect_jumps

import jwst

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

__all__ = ["JumpStep"]


class JumpStep(RomanStep):
    """
    JumpStep: Performs CR/jump detection on each ramp integration within an
    exposure. The 2-point difference method is applied.
    """

    spec = """
        rejection_threshold = float(default=4.0,min=0) # CR sigma rejection threshold
        maximum_cores = option('none', 'quarter', 'half', 'all', default='none') # max number of processes to create
        flag_4_neighbors = boolean(default=True) # flag the four perpendicular neighbors of each CR
        max_jump_to_flag_neighbors = float(default=200) # maximum jump sigma that will trigger neighbor flagging
        min_jump_to_flag_neighbors = float(default=10) # minimum jump sigma that will trigger neighbor flagging
    """

    reference_file_types = ['gain', 'readnoise']

    def process(self, input):

       # Open input as a Roman DataModel (single integration; 3D arrays)
        with datamodels.RampModel(input) as R_input_model:
            # Create a JWST DataModel (4D arrays) to populate
            input_model = jwst.datamodels.RampModel()

            # Populate JWST model with correct sized arrays from Roman model
            input_model.meta = R_input_model.meta
            input_model.data = np.broadcast_to(R_input_model.data, (1,) +
                                               R_input_model.data.shape)

            input_model.pixeldq = R_input_model.pixeldq

            input_model.groupdq = np.broadcast_to(R_input_model.groupdq, (1,) +
                                                  R_input_model.groupdq.shape)

            input_model.err = np.broadcast_to(R_input_model.err, (1,) +
                                              R_input_model.err.shape)

            input_model.refout = np.broadcast_to(R_input_model.refout, (1,) +
                                                 R_input_model.refout.shape)

            tstart = time.time()
            # Check for an input model with NGROUPS<=2
            ngroups = input_model.data.shape[1]
            if ngroups <= 4:
                self.log.warning('Can not apply jump detection when NGROUPS<=4;')
                self.log.warning('Jump step will be skipped')
                result = input_model.copy()
                result.meta.cal_step.jump = 'SKIPPED'
                return result

            # Retrieve the parameter values
            rej_thresh = self.rejection_threshold
            max_cores = self.maximum_cores
            max_jump_to_flag_neighbors = self.max_jump_to_flag_neighbors
            min_jump_to_flag_neighbors = self.min_jump_to_flag_neighbors
            flag_4_neighbors = self.flag_4_neighbors

            self.log.info('CR rejection threshold = %g sigma', rej_thresh)
            if self.maximum_cores != 'none':
                self.log.info('Maximum cores to use = %s', max_cores)

            # Get the gain and readnoise reference files
            gain_filename = self.get_reference_file(input_model, 'gain')
            self.log.info('Using GAIN reference file: %s', gain_filename)
            gain_model = datamodels.GainModel(gain_filename)

            readnoise_filename = self.get_reference_file(input_model,
                                                         'readnoise')
            self.log.info('Using READNOISE reference file: %s',
                          readnoise_filename)
            readnoise_model = datamodels.ReadnoiseModel(readnoise_filename)


            result = detect_jumps(input_model, gain_model, readnoise_model, rej_thresh,
                                  max_cores, max_jump_to_flag_neighbors,
                                  min_jump_to_flag_neighbors, flag_4_neighbors)

            gain_model.close()
            readnoise_model.close()
            tstop = time.time()
            self.log.info('The execution time in seconds: %f', tstop - tstart)

        result.meta.cal_step.jump = 'COMPLETE'

        return result
