Arguments
=========

The source detection fitting step has several arguments. These can be specified
by the user by passing them to the step in a Python session, or setting them
in a parameter file.

* ``--kernel_fwhm``: A parameter for DAOStarFinder: size of Gaussian kernel in
                     pixels. By default the FWHM is assumed to be the diffraction
                     limited PSF, given the filter used for this observation.
* ``--sharplo``: A parameter for DAOStarFinder: lower bound for sharpness.
                 Default is 0.0.
* ``--sharphi``: A parameter for DAOStarFinder: upper bound for sharpness.
                 Default is 1.0.
* ``--roundlo``: A parameter for DAOStarFinder: lower bound for roundness.
                 Default is -1.0. A circular source will have a zero roundness.
                 A source extended in x or y will have a negative or positive
                 roundness, respectively.
* ``--roundhi``: A parameter for DAOStarFinder: upper bound for roundness.
                 Default is 1.0. A circular source will have a zero roundness.
                 A source extended in x or y will have a negative or positive
                 roundness, respectively.
* ``--peakmax``: A parameter for DAOStarFinder: upper limit on brightest pixel
                 in sources. Default is 1000.0.
* ``--max_sources``: The maximum number of sources in the output catalog,
                     choosing brightest. Default is None, which will return all
                     detected sources.
* ``--scalar_threshold``: If specified, the absolute detection threshold to be
                          used for the entire image. Units are assumed to be the
                          same as input data. One of `scalar_threshold`,
                          `calc_scalar_threshold` must be chosen. Default is
                          None.
* ``--calc_scalar_threshold``: If specified, a single scalar threshold will be
                               determined for the entire image. This is done by
                               calculating a 2D background image, and using that
                               to determine a single threshold value for the
                               entire image. One of `scalar_threshold` or
                               `calc_scalar_threshold` must be chosen.
                               must be chosen. Default is True.
* ``--snr_threshold``: If using `calc_threshold_img`, the SNR for the threshold
                       image. Default is 3.0.
* ``--bkg_estimator``: If using `calc_threshold_img`, choice of mean, median, or
                        mode. Default is median.
* ``--bkg_boxsize``: If using `calc_threshold_img` size of box in pixels for
                     2D background / threshold images and if using
                     calc_threshold_2d, the size of the box used when detecting
                     sources. Default is 9.
* ``--bkg_sigma``: If using `calc_threshold_img`, n sigma for sigma clipping
                   for background calculation. Default is 2.0.
* ``--bkg_filter_size``: If using `calc_threshold_img` or `calc_threshold_2d`,
                         size of square gaussian kernel for background
                         calculation. Default is 3.
* ``--save_catalogs``: A True/False value that specifies whether to write
                      the optional output catalog. Default is False.
* ``--output_cat_filetype``: If `save_catalogs` is True, file type of output
                             catalog from choice of asdf and escv. Catalog
                             will be saved as a numpy array with four dimensions.
                             In order, these represent source ID, x centroid
                             position, y centroid position, and flux.
* ``--fit_psf``: If True, fit a PSF model to each detected source for more precise
                 source centroids and fluxes.
