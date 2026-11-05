Step Arguments
==============

The ``time_domain_catalog`` step has the following arguments:

* ``--forced_catalogs``: A list of the catalog files providing positions and
  IDs for forced photometry.

* ``--id_columns``: A list of the ID column names of the catalogs, one
  per catalog, or a single name used for all of them (default is
  ``flagged_spatial_id``, appropriate for romancal multiband catalogs).

* ``--bkg_boxsize``: A float value giving the background mesh box
  size in arcsec (default is 55.0, as in the source catalog step).

* ``--data_release_id``: The data release identifier to record in the
  catalog metadata (``meta.data_release_id``).  Not recorded if not given.

* ``--suffix``: A string value giving the file name suffix to use for
  the output catalog file (default is ``'fcat'``).
