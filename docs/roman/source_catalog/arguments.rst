Step Arguments
==============

The ``source_catalog`` step has the following arguments:

* ``--bkg_boxsize``: An integer value giving the background mesh box
  size in pixels

* ``--kernel_fwhm``: A floating-point value giving the FWHM in arcsec
  of the point-source detection template.  The larger templates are
  angular sizes as well, so the bank searches the same physical scales
  whatever the pixel scale of the image.

* ``--snr_threshold``: A floating-point value that sets the
  signal-to-noise ratio threshold above the background for source
  detection.

* ``--npixels``: An integer value that sets the minimum number of
  pixels a source segment must retain to be kept, counting only the
  unmasked pixels that are positive in the convolved image the moments
  are measured on

* ``--deblend``: A boolean indicating whether to deblend sources (default
  is ``True``)

* ``--max_sources``: An integer value giving the maximum number of
  sources to keep, retaining the most significant; zero means no limit

* ``--suffix``: A string value giving the file name suffix to use for
  the output catalog file (default is ``'cat'``).

* ``--fit_psf``: A boolean value indicating whether to perform PSF
  photometry (default is ``True``)

* ``--forced_segmentation``: A string value indicating the filename of
  the segmentation map to use for forced segmentation

* ``--compute_skyvals``: A boolean value indicating whether to compute
  and attach HEALPix sky summary arrays (``skyvals`` and
  ``healpix11_cov``) to L2 segmentation-map outputs (default is
  ``True``)
