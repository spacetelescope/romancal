Description
===========

:Class: `romancal.time_domain_catalog.TimeDomainCatalogStep`
:Alias: time_domain_catalog

This step measures the fluxes of sources in a Level 2 image at fixed
sky positions taken from input catalogs, such as catalogs of sources
detected in deep coadds.  Running it on each exposure of a field gives
light curves of those sources.

The image is masked, background subtracted, and converted to flux
density in the same way as in the :ref:`source catalog step
<source_catalog_step>`.  The ``ra`` and ``dec`` of each input source
(ICRS, degrees; proper motion is ignored) are converted to pixel
positions with the image WCS, and sources that do not land on the
image are dropped.  Each remaining source then gets three
measurements:

* a PSF fit in which the position and flux are both free
  (``x_psf``, ``y_psf``, ``ra_psf``, ``dec_psf``, ``psf_flux``, and
  their uncertainties, flags, and goodness of fit),
* a PSF fit in which the position is fixed at the input position and
  only the flux is free (``forced_psf_flux`` and its uncertainty, flags,
  and goodness of fit), and
* circular-aperture photometry and a local background at the input
  position (``forced_aperXX_flux`` and ``forced_aper_bkg_flux``).  No
  pixels of neighboring sources are masked and no aperture corrections
  are applied.

Nearby sources are fit simultaneously in both PSF fits, as in the
source catalog step.  The ``image_flags`` and ``image_flags2`` columns
are the bitwise or of the ``dq`` and ``dq2`` arrays over the 5x5 pixel
box around each input position.

Input catalogs
--------------

Any number of input catalogs may be given, in any format that
`astropy.table.Table.read` can read.  Each must have ``ra`` and ``dec``
columns and an ID column, named by ``id_columns``.  The output ``id``
column holds the ID of each source, and ``source_catalog_index`` is
the index of the catalog it came from in
``meta.forced_psf_catalogs``.  When the ID column is
``flagged_spatial_id``, sources outside the core of their skycell are
dropped so that sources in skycell overlap regions are measured only
once.

Output
------

The output is a ``TimeDomainSourceCatalogModel`` saved as a parquet
file with the suffix ``fcat``.
