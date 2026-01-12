.. _wfi18_transient_arguments:

Step Arguments
==============
The ``wfi18_transient`` step uses the following optional arguments:

``--mask_rows`` (boolean, default=False)
  If True, the rows in the first resultant that are most affected by
  the transient anomaly will be set to DO_NOT_USE in the GROUPDQ
  array.  If False, a fit to the transient anomaly is attempted, and
  the fit correction is subtracted from the first resultant.
