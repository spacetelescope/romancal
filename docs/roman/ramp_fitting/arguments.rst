.. _rampfit-arguments:

Arguments
=========
The ramp fitting step has the following optional argument that can be set by the user:

* ``--algorithm``: Algorithm to use. At the moment the only supported option
  is  `ols_cas22` (the default).

* ``--use_ramp_jump_detection``: A True/False value that specifies whether to use
  the unevenly-spaced jump detection integrated into the ramp fitting algorithm.
  If True, then the jump detection step will be skipped and then revisited alongside
  the ramp fitting step. If False, then the jump detection step will be run. The
  default is True.
