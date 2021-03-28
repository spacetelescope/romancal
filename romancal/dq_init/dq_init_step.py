#! /usr/bin/env python

from romancal.stpipe import RomanStep
from stcal.dq_init.dq_init_step import DQInitStep

__all__ = ["DQInitStep"]


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
#     return dq_init.DQInitStep(RomanStep)
