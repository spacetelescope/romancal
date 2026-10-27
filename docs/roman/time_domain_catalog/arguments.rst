Step Arguments
==============

The ``time_domain_catalog`` step has the following arguments:

* ``--forced_catalogs``: A list of the catalog files supplying the
  source positions.

* ``--id_columns``: A list of the ID column names of the catalogs, one
  per catalog, or a single name used for all of them (default is
  ``flagged_spatial_id``).

* ``--bkg_boxsize``: An integer value giving the background mesh box
  size in pixels.

* ``--suffix``: A string value giving the file name suffix to use for
  the output catalog file (default is ``'fcat'``).
