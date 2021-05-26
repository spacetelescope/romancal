#! /usr/bin/env python

from romancal.stpipe import RomanStep
#from stcal.dq_init.dq_init_step import DQInitStep
from stcal.dq_init import dq_init_step

__all__ = ["DQInitStep"]

class DQInitStep(RomanStep,dq_init_step.DQInitStep):
    """Initialize the Data Quality extension from the
    mask reference file.

    The dq_init step initializes the pixeldq attribute of the
    input datamodel using the MASK reference file.  For some
    Guiding exp_types, initalize the dq attribute of the input model
    instead.  The dq attribute of the MASK model is bitwise OR'd
    with the pixeldq (or dq) attribute of the input model.
    """
    pass




# class DQInitStep(dq_init_step.DQInitStep):
#     """Initialize the Data Quality extension from the
#     mask reference file.
#
#     The dq_init step initializes the pixeldq attribute of the
#     input datamodel using the MASK reference file.  For some
#     Guiding exp_types, initalize the dq attribute of the input model
#     instead.  The dq attribute of the MASK model is bitwise OR'd
#     with the pixeldq (or dq) attribute of the input model.
#     """
#
#     def __init__(self, RomanStep):
#         self.__init_subclass__()class__ = RomanStep
#
#     pass

# class DQInitStep(RomanStep):
#     """Initialize the Data Quality extension from the
#     mask reference file.
#
#     The dq_init step initializes the pixeldq attribute of the
#     input datamodel using the MASK reference file.  For some
#     Guiding exp_types, initalize the dq attribute of the input model
#     instead.  The dq attribute of the MASK model is bitwise OR'd
#     with the pixeldq (or dq) attribute of the input model.
#     """
#
#     def process(self):
#         """Perform the dq_init calibration step
#
#         Parameters
#         ----------
#         input : science datamodel
#             input science datamodel
#
#         Returns
#         -------
#         output_model : science datamodel
#             result science datamodel
#         """
#         return RomanStep(dq_init_step.DQInitStep())
