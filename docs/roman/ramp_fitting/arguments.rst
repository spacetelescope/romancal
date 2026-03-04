.. _rampfit-arguments:

Arguments
=========
The ramp fitting step has the following optional argument that can be set by the user:

* ``--algorithm``: Algorithm to use, `ols_cas22` (default) or `likely`.

* ``--use_ramp_jump_detection``: A True/False value that specifies whether to use
  the unevenly-spaced jump detection integrated into the ramp fitting algorithm.
  If True, then the jump detection step will be skipped and then revisited alongside
  the ramp fitting step. If False, then the jump detection step will be run. The
  default is True.

* ``--include_var_rnoise``: boolean indicating whether to include var_rnoise in output (can be reconstructed from err and other variances)

* ``--rejection_threshold``: float specifying the CR sigma rejection threshold.

* ``--threshold_intercept``: float overriding the intercept parameter for the threshold function in the jump detection algorithm.

* ``--threshold_constant``: float overriding the constant parameter for the threshold function in the jump detection algorithm.

* ``--maximum_cores``: string/integer specifying the cores for multiprocessing. Can be an integer, 'half', 'quarter', or 'all'

* ``--expand_large_events``: boolean to turn on Snowball detection.

* ``--min_sat_area``: float specifying minimum required area for the central saturation of snowballs.

* ``--min_jump_area``: float specifying minimum area to trigger large events processing.

* ``--expand_factor``: float specifying the expansion factor for the enclosing circles or ellipses.

* ``--sat_required_snowball``: boolean specifying whether to  require the center of snowballs to be saturated.

* ``--min_sat_radius_extend``: float specifying the min radius of the sat core to trigger the extension of the core.

* ``--sat_expand``: integer specifying the number of pixels to add to the radius of the saturated core of snowballs.

* ``--edge_size``: integer specifying the distance from detector edge where a saturated core is not required for snowball detection.
