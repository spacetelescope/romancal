.. _rampfit-arguments:

Arguments
=========
The ramp fitting step has the following optional argument that can be set by the user:

* ``--algorithm``: Algorithm to use. Possible values are `ols` and `ols_cas22`.
  `ols` is the same algorithm used by JWST and can only be used with even ramps.
  `ols_cas22` needs to be used for uneven ramps. `ols_cas22` is the default.

The following optional arguments are valid only if using the `ols` algorithm.

* ``--save_opt``: A True/False value that specifies whether to write
  the optional output product. Default if False.

* ``--opt_name``: A string that can be used to override the default name
  for the optional output product.

* ``--maximum_cores``: The fraction of available cores that will be
  used for multi-processing in this step. The default value is 'none' which does not use
  multi-processing. The other options are 'quarter', 'half', and 'all'. Note that these
  fractions refer to the total available cores and on most CPUs these include physical
  and virtual cores. The clock time for the step is reduced
  almost linearly by the number of physical cores used on all machines. For example, on an Intel CPU with
  six real cores and 6 virtual cores setting maximum_cores to 'half' results in a
  decrease of a factor of six in the clock time for the step to run. Depending on the system
  the clock time can also decrease even more with maximum_cores is set to 'all'.

* ``--use_ramp_jump_detection``: A True/False value that specifies whether to use
  the unevenly-spaced jump detection integrated into the ramp fitting algorithm.
  If True, then the jump detection step will be skipped and then revisited alongside
  the ramp fitting step. If False, then the jump detection step will be run. The
  default is True.
