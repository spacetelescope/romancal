PSF Fitting
===========

A few PSF fitting utilities are included to interface between observations
within Roman datamodels and methods within dependencies that generate and
fit PSF models to observations.

Create PSF models
-----------------

`~romancal.lib.psf.create_gridded_psf_model` computes a gridded PSF model for
a given detector using `~webbpsf.gridded_library.CreatePSFLibrary` from
`WebbPSF <https://webbpsf.readthedocs.io/>`_. The defaults are chosen to
balance more accurate PSF models with the cost of increased runtime. For
further reading on the WebbPSF approach to ePSFs, see the WebbPSF docs on
`Using PSF Grids <https://webbpsf.readthedocs.io/en/latest/psf_grids.html>`_.

Fit model PSFs to an ImageModel
-------------------------------

Once PSF models are generated, you can fit those models to observations
in an ImageModel with `~romancal.lib.psf.fit_psf_to_image_model` using
`photutils <https://photutils.readthedocs.io/en/stable/psf.html>`_.
By default the fits are done with `~photutils.psf.PSFPhotometry`, and
crowded fields may benefit from using `~photutils.psf.IterativePSFPhotometry`.
For neighboring sources that are near one another on the detector, grouping
the sources and fitting their PSFs simultaneously improves the fit quality.
Initial guesses for target centroids can be given or source
detection can be performed with, e.g., `~photutils.detection.DAOStarFinder`.
