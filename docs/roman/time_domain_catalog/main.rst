Description
===========

:Class: `romancal.time_domain_catalog.TimeDomainCatalogStep`
:Alias: time_domain_catalog

This step measures the fluxes of sources in a Level 2 image at arbitrary
locations provided in an input catalog.  For example, it can be used to
find the fluxes of sources at the location of sources identified in deep
coadds, to produce light curves.

The step begins by calibrating the image in the same way as the
:ref:`source catalog step <source_catalog_step>`, performing masking,
background subtraction, and flux conversion.  The ``ra`` and ``dec`` of
each input source (ICRS, degrees) are converted to pixel
positions with the image WCS.  Sources that land off the image are dropped.
The remaining sources get three measurements:

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

Each source is fit independently in both PSF fits; unlike in the
source catalog step, neighboring sources are not fit simultaneously
or masked.  The
``image_flags`` and ``image_flags2`` columns
are the bitwise or of the ``dq`` and ``dq2`` arrays over the 5x5 pixel
box around each input position.

Input catalogs
--------------

Any number of input catalogs may be given, in any format that
`astropy.table.Table.read` can read.  Each must have ``ra`` and ``dec``
columns and an ID column, named by ``id_columns``.  The output ``id``
column records the ID of each source.  The list of input PSF catalogs
used for forced photometry is stored in ``meta.forced_psf_catalogs``,
and ``source_catalog_index`` records which input catalog each source
came from.  When the ID column is
``flagged_spatial_id``, sources outside the core of their skycell are
dropped so that sources in skycell overlap regions are measured only
once.

Output
------

The output is a ``TimeDomainSourceCatalogModel`` saved as a parquet
file with the suffix ``fcat``.
