PSF Fitting
===========

A few PSF fitting utilities are included to interface between observations
within Roman datamodels and methods within dependencies that generate and
fit PSF models to observations.

Create PSF models
-----------------

`~romancal.source_catalog.psf.create_gridded_psf_model`
computes a gridded PSF model for a given detector using
the reference files in CRDS.

These PSF models account for jitter by deconvolving the amount of jitter
indicated to be present in the PSF reference file, and reconvolving with
the amount of jitter present in the individual image.  Because the amount of
jitter present in the reference file is small (~8 mas), this deconvolution
does not introduce significant noise.  The convolution process works in
the Fourier domain using the approach of Lang (2020) because the jitter
kernel would otherwise be badly undersampled.


Fit model PSFs to an ImageModel
-------------------------------

Once PSF models are generated, you can fit those
models to observations in an ImageModel with
`~romancal.source_catalog.psf.fit_psf_to_image_model` using `photutils
<https://photutils.readthedocs.io/en/stable/psf.html>`_. By default
the fits are done with `~photutils.psf.PSFPhotometry`, and crowded
fields may benefit from using `~photutils.psf.IterativePSFPhotometry`.
For neighboring sources that are near one another on the detector,
grouping the sources and fitting their PSFs simultaneously
improves the fit quality. Initial guesses for target centroids
can be given or source detection can be performed with, e.g.,
`~photutils.detection.DAOStarFinder`.
