Description
============

The source detection step produces catalogs of point-like sources for use by the
Tweakreg step for image alignment. It uses DAOStarFinder to detect point sources
in the image, with an option to subsequently fit PSF models to the detected
sources for more precise centroids and fluxes.

Detecting Sources
-----------------

Sources are detected using `~photutils.detection.DAOStarFinder` from
`photutils <https://photutils.readthedocs.io/en/stable/>`_, which is an
implementation of the method `DAOFIND` from
`Stetson (1987) <https://ui.adsabs.harvard.edu/abs/1987PASP...99..191S/abstract>`_.
The algorithm can be provided limits on the source flux, radius, roundness,
sharpness, and background.

PSF Fitting
-----------

Star finding algorithms like `~photutils.detection.DAOStarFinder` provide
approximate stellar centroids. More precise centroids may be inferred by
fitting model PSFs to the observations. Setting the SourceDetectionStep's
option `fit_psf` to True will generate model Roman PSFs with
`WebbPSF <https://webbpsf.readthedocs.io/en/latest/roman.html>`_, and fit
those models to each of the sources detected by
`~photutils.detection.DAOStarFinder`. More details are in :doc:`psf`.

Outputs / Returns
-----------------

By default, the resulting source catalog will be temporarily attached to the
output ImageModel in the `meta.source_catalog.tweakreg_catalog` attribute as
numpy array representing, in order, source ID, x centroid position, y centroid
position, and flux. This catalog will then be deleted from the model in the
Tweakreg step.

Optionally, the catalog can be saved to disk in which case a
`meta.source_catalog.tweakreg_catalog_name` attribute will be added to the file
to point Tweakreg to the catalog on disk. To do this, set `save_catalogs` to
True. Output catalogs will be saved in the same directory as input files, and
are also 4D numpy arrays representing, in order, source ID, x centroid position,
y centroid position, and flux. Output catalogs can be in ASDF or ECSV format.

NOTE: The intermediate resulting ImageModel from SourceDetectionStep can
only be saved if it does not contain an attached catalog - to do this, use the
`save_catalogs` option to separate the catalog from the file and save them
separately.

Options for Thresholding
------------------------

The DAOStarFinder routine detects point-like sources in an image that are above
a certain, specified floating point threshold. This step provides several options
for calculating this threshold, either using one value for the entire image,
or by detecting sources in segments of the image and using a different appropriate
threshold for each (useful if background varies across the image).

The first option is to set `scalar_threshold` - this will use the specified
threshold as the detection threshold for the entire image.

The second option is to use `calc_threshold` - this will calculate a single
threshold value for the entire image based on the sigma-clipped average
(mean, median, or mode) background level of the whole image.

Other Options
-------------

Limiting maximum number of sources
++++++++++++++++++++++++++++++++++

By default, all detected sources will be returned in the final output catalog.
If you wish to limit the number of sources, this can be done with the
`max_sources` argument, which will sort the output catalog by flux and return
only the N brightest.
