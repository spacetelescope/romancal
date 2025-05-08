Step Arguments
==============

The ``source_catalog`` step has the following arguments:

* ``--bkg_boxsize``: An integer value giving the background mesh box
  size in pixels

* ``--kernel_fwhm``: A floating-point value giving the Gaussian kernel
  FWHM in pixels

* ``--snr_threshold``: A floating-point value that sets the
  signal-to-noise ratio (SNR) threshold above the background for source
  detection.

* ``--npixels``: An integer value that sets the minimum number of
  pixels in a source

* ``--deblend``: A boolean indicating whether to deblend sources (default
  is ``False``)

* ``--suffix``: A string value giving the file name suffix to use for
  the output catalog file (default is ``'cat'``).

* ``--fit_psf``: A boolean value indicating whether to perform PSF
  photometry (default is ``True``)

* ``--forced_segmentation``: A string value indicating the filename of
  the segmentation map to use for forced segmentation
