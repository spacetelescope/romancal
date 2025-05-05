Description
===========

:Class: `romancal.source_catalog.SourceCatalogStep`
:Alias: source_catalog

This step detects sources in an image and generates a catalog of their
properties, including photometric and shape measurements.


Background Subtraction
----------------------

A two-dimensional background is estimated and subtracted from the
data. The background and background noise are estimated using the
:py:class:`photutils.background.Background2D` class from `Photutils
<https://photutils.readthedocs.io/en/stable/index.html>`_. This class
calculates the background by measuring the sigma-clipped median within
user-defined boxes of a specified size (``bkg_boxsize``). The background
RMS noise is then estimated using the sigma-clipped standard deviation
within the same boxes.


Source Detection
----------------

Sources are detected using `image segmentation
<https://en.wikipedia.org/wiki/Image_segmentation>`_, a process that
assigns an integer label to each pixel in an image, such that pixels
with the same label correspond to the same source. Source extraction is
then performed using the `Photutils segmentation <https://photutils.readthedocs.io/en/latest/user_guide/segmentation.html>`_ tools.

The background-subtracted image is convolved with a Gaussian kernel
to smooth the noise and enhance the detectability of objects with a
shape similar to the kernel. The Gaussian kernel is defined by
the ``kernel_fwhm`` parameter, which specifies the full-width at half-maximum
(FWHM) of the kernel. The kernel is normalized to have a total
flux of 1.0.

Detected sources must consist of a minimum number of connected pixels
(``npixels``), each exceeding a specified threshold value in
the convolved image. The threshold level is defined as a per-pixel
multiple (``snr_threshold``) of the background RMS image.

Overlapping sources are deblended using the `Photutils deblender
<https://photutils.readthedocs.io/en/latest/user_guide/segmentation.html
#source-deblending>`_. The deblending algorithm first applies
a multi-thresholding approach to identify potentially
overlapping sources, then uses `watershed segmentation
<https://en.wikipedia.org/wiki/Watershed_(image_processing)>`_
to separate them. For successful deblending, the sources must be
sufficiently separated so that a saddle exists between them. Currently,
the ``deblend`` keyword must be set to deblend sources.


Source Photometry and Properties
--------------------------------

After detecting sources using image segmentation, we can measure their
photometry, centroids, and shape/morphological properties.

The source centroids and shape properties are derived from 2D image
moments of the pixel values within the source segments. These properties
include the semimajor and semiminor axes, ellipticity, and orientation
of the major axis.

Circular aperture photometry is performed at several aperture sizes
(:math:`r` = 0.1, 0.2, 0.4, 0.8, 1.6 arcsec) for each source. Elliptical
Kron aperture photometry is also performed, where the aperture size is
determined by the source shape.

Isophotal photometry is measured using the total flux within the source
segment.

Optionally, Point Spread Function (PSF) photometry can be
performed by setting the ``fit_psf`` keyword. Enabling
this option fits a model PSF to each source to measure its
position and flux. The PSF model is generated using the
`STPSF <https://stpsf.readthedocs.io/en/latest/roman.html>`_
package. PSF photometry is performed using the
:py:class:`photutils.psf.PSFPhotometry` class from Photutils.

A local background is estimated using a circular annulus around the
source. The annulus is defined by the inner and outer radii of 2.4 and
2.8 arcsec, respectively. The local background flux is calculated as
the sigma-clipped median value within the annulus. Although this local
background value is included in the source catalog, it is not subtracted
from any of the measured fluxes.

All fluxes are reported in nJy. To calculate AB magnitudes, use the
following formula:

.. math::

    m_{\rm AB} = -2.5 \log_{10}(f_{\rm nJy}) + 31.4

Photometric errors are calculated from the resampled total-error array
contained in the ``model.err`` array. Note that this array includes
source Poisson noise.


Output Products
---------------

Source Catalog Table
^^^^^^^^^^^^^^^^^^^^

The output source catalog table is saved to a file in the `Parquet
<https://parquet.apache.org/>`_ format.

The table contains one row for each source, with the columns listed
below (assuming PSF-photometry is requested).

All pixel coordinates are 0-indexed, following Python's 0-based
indexing. This means pixel coordinate 0 corresponds to the center of the
first pixel.


+------------------------+----------------------------------------------------+
| Column                 | Description                                        |
+========================+====================================================+
| label                  | Unique source identification label number          |
+------------------------+----------------------------------------------------+
| xcentroid              | X pixel value of the source centroid (0 indexed)   |
+------------------------+----------------------------------------------------+
| ycentroid              | Y pixel value of the source centroid (0 indexed)   |
+------------------------+----------------------------------------------------+
| ra/dec_centroid        | ra/dec coordinate of the source centroid           |
+------------------------+----------------------------------------------------+
| aper_bkg_flux          | The local background value calculated as the       |
|                        | sigma-clipped median value in the background       |
|                        | annulus aperture                                   |
+------------------------+----------------------------------------------------+
| aper_bkg_flux_err      | The standard error of the sigma-clipped median     |
|                        | background value                                   |
+------------------------+----------------------------------------------------+
| aper30_flux            | Flux within the 30% encircled energy circular      |
|                        | aperture                                           |
+------------------------+----------------------------------------------------+
| aper30_flux_err        | Flux error within the 30% encircled energy         |
|                        | circular aperture                                  |
+------------------------+----------------------------------------------------+
| aper50_flux            | Flux within the 50% encircled energy circular      |
|                        | aperture                                           |
+------------------------+----------------------------------------------------+
| aper50_flux_err        | Flux error within the 50% encircled energy         |
|                        | circular aperture                                  |
+------------------------+----------------------------------------------------+
| aper70_flux            | Flux within the 70% encircled energy circular      |
|                        | aperture                                           |
+------------------------+----------------------------------------------------+
| aper70_flux_err        | Flux error within the 70% encircled energy         |
|                        | circular aperture                                  |
+------------------------+----------------------------------------------------+
| aper_total_flux        | Total aperture-corrected flux based on the 70%     |
|                        | encircled energy circular aperture; should be used |
|                        | only for unresolved sources                        |
+------------------------+----------------------------------------------------+
| aper_total_flux_err    | Total aperture-corrected flux error based on the   |
|                        | 70% encircled energy circular aperture; should be  |
|                        | used only for unresolved sources                   |
+------------------------+----------------------------------------------------+
| flags                  | Flag recording DQ value of source central pixel    |
+------------------------+----------------------------------------------------+
| is_extended            | Flag indicating whether the source is extended     |
+------------------------+----------------------------------------------------+
| sharpness              | The DAOFind source sharpness statistic             |
+------------------------+----------------------------------------------------+
| roundness              | The DAOFind source roundness statistic             |
+------------------------+----------------------------------------------------+
| nn_label               | The label number of the nearest neighbor           |
+------------------------+----------------------------------------------------+
| nn_dist                | The distance in pixels to the nearest neighbor     |
+------------------------+----------------------------------------------------+
| isophotal_flux         | Isophotal flux                                     |
+------------------------+----------------------------------------------------+
| isophotal_flux_err     | Isophotal flux error                               |
+------------------------+----------------------------------------------------+
| isophotal_area         | Isophotal area                                     |
+------------------------+----------------------------------------------------+
| kron_flux              | Kron flux                                          |
+------------------------+----------------------------------------------------+
| kron_flux_err          | Kron flux error                                    |
+------------------------+----------------------------------------------------+
| semimajor_sigma        | 1-sigma standard deviation along the semimajor     |
|                        | axis of the 2D Gaussian function that has the same |
|                        | second-order central moments as the source         |
+------------------------+----------------------------------------------------+
| semiminor_sigma        | 1-sigma standard deviation along the semiminor     |
|                        | axis of the 2D Gaussian function that has the same |
|                        | second-order central moments as the source         |
+------------------------+----------------------------------------------------+
| ellipticity            | 1 minus the ratio of the 1-sigma lengths of the    |
|                        | semimajor and semiminor axes                       |
+------------------------+----------------------------------------------------+
| orientation            | The angle (degrees) between the positive X axis    |
|                        | and the major axis (increases counter-clockwise)   |
+------------------------+----------------------------------------------------+
| sky_orientation        | The position angle (degrees) from North of the     |
|                        | major axis                                         |
+------------------------+----------------------------------------------------+
| ra_bbox_ll, dec_bbox_ll| Sky coordinate of the lower-left vertex of the     |
|                        | minimal bounding box of the source                 |
+------------------------+----------------------------------------------------+
| ra_bbox_ul, dec_bbox_ul| Sky coordinate of the upper-left vertex of the     |
|                        | minimal bounding box of the source                 |
+------------------------+----------------------------------------------------+
| ra_bbox_lr, dec_bbox_lr| Sky coordinate of the lower-right vertex of the    |
|                        | minimal bounding box of the source                 |
+------------------------+----------------------------------------------------+
| ra_bbox_ur, dec_bbox_ur| Sky coordinate of the upper-right vertex of the    |
|                        | minimal bounding box of the source                 |
+------------------------+----------------------------------------------------+


If ``fit_psf=True``, the following columns will also be available:

+------------------------+----------------------------------------------------+
| Column                 | Description                                        |
+========================+====================================================+
| x_psf, x_psf_err       | X pixel value of the source and its associated     |
|                        | error as determined by PSF fitting                 |
+------------------------+----------------------------------------------------+
| y_psf, y_psf_err       | Y pixel value of the source and its associated     |
|                        | error as determined by PSF fitting                 |
+------------------------+----------------------------------------------------+
| flux_psf, flux_psf_err | Flux of the source and its associated error as     |
|                        | determined by PSF fitting                          |
+------------------------+----------------------------------------------------+
| flag_psf               | DQ flag of the resulting PSF fitting.              |
|                        | Possible values are [1]_:                          |
|                        |                                                    |
|                        | - 1 : one or more pixels in the fitting region     |
|                        |   were masked                                      |
|                        | - 2 : the fit x and/or y position lies outside of  |
|                        |   the input data                                   |
|                        | - 4 : the fit flux is less than or equal to zero   |
|                        | - 8 : the fitter may not have converged            |
|                        | - 16 : the fitter parameter covariance matrix was  |
|                        |   not returned                                     |
+------------------------+----------------------------------------------------+


Segmentation Map
^^^^^^^^^^^^^^^^

The segmentation map generated during the
source-finding process is saved as an `ASDF
<https://en.wikipedia.org/wiki/Advanced_Scientific_Data_Format>`_ file.
Each pixel in the image contains an integer value corresponding to a
source label in the source catalog. Pixels that do not belong to any
source are assigned a value of zero.


Multiband Catalogs
------------------

Multiband catalogs combine multiple images to create a deep detection
image, which is used to detect sources and identify segments. The
measured positions and shapes of the sources in these deep images are
then used to perform aperture, Kron, isophotal, and PSF photometry for
each filter.

The catalog fields are similar to those in the source catalog schema,
but with the following differences:

* Fields derived from individual filter images include the
  filter name from which they were derived. For example, fields
  like ``aper_flux_<filter>``, ``segment_flux_<filter>``,
  ``kron_flux_<filter>``, and ``psf_flux_<filter>`` provide the aperture
  and PSF flux for each filter, respectively.

* Fields derived from the detection image and segmentation map do not
  include the filter name.

Multiband catalogs are generated by the
:py:class:`~romancal.multiband_catalog.MultibandCatalogStep`, which
takes an association file as input. This file lists the images that need
to be photometered simultaneously.


Forced Source Catalogs
----------------------

Source catalogs can optionally be generated by using the segmentation
image from one image (the "forcing" image) and computing shapes and
fluxes for those same segments in another image (the "forced" image).
For this to work, the two images must be perfectly aligned in pixel
space.

Forced source catalogs can be generated by specifying a segmentation
image with the ``forced_segmentation`` keyword when running the source
catalog step.

In this mode, the source catalog contains fields with the ``forced``
prefix, in addition to the fields described above. Fields without the
"forced" prefix contain position and shape information derived from
the forcing image, indicating where measurements were taken on the
forced image. Fields with the forced prefix represent values computed
on the forced image, using the information from the forcing image.
For example, the field ``forced_kron_flux`` represents the Kron flux
measured on the forced image, using the centroid and shape information
from the ``x_centroid``, ``y_centroid``, ``semimajor``, ``semiminor``,
and ``orientation_pix`` fields.
