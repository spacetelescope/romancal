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


Source Catalog Table
--------------------

The source catalog table contains one row for each source, with the
columns listed below (assuming PSF-photometry is requested).

All pixel coordinates are 0-indexed, following Python's 0-based
indexing. This means pixel coordinate 0 corresponds to the center of the
first pixel.

All sky coordinates are in decimal degrees in the International
Celestial Reference System (ICRS) reference frame.

Uncertainties are reported as the 1-sigma (68.27% confidence) errors.


+-----------------------+-----------------------------------------------------+
| Column                | Description                                         |
+=======================+=====================================================+
| label                 | Label of the source segment in the segmentation     |
|                       | image                                               |
+-----------------------+-----------------------------------------------------+
| flagged_spatial_index | Bit flag encoding the overlap flag, projection,     |
|                       | skycell, and pixel coordinates of the source        |
+-----------------------+-----------------------------------------------------+
| x_centroid            | Column coordinate of the source centroid in the     |
|                       | detection image from image moments (0 indexed)      |
+-----------------------+-----------------------------------------------------+
| y_centroid            | Row coordinate of the source centroid in the        |
|                       | detection image from image moments (0 indexed)      |
+-----------------------+-----------------------------------------------------+
| x_centroid_err        | Uncertainty in x_centroid                           |
+-----------------------+-----------------------------------------------------+
| y_centroid_err        | Uncertainty in y_centroid                           |
+-----------------------+-----------------------------------------------------+
| x_centroid_win        | Column coordinate of the windowed source centroid   |
|                       | in the detection image from image moments (0        |
|                       | indexed)                                            |
+-----------------------+-----------------------------------------------------+
| y_centroid_win        | Row coordinate of the windowed source centroid in   |
|                       | the detection image from image moments (0 indexed)  |
+-----------------------+-----------------------------------------------------+
| x_centroid_win_err    | Uncertainty in x_centroid_win                       |
+-----------------------+-----------------------------------------------------+
| y_centroid_win_err    | Uncertainty in y_centroid_win                       |
+-----------------------+-----------------------------------------------------+
| x_psf                 | Column coordinate of the source from PSF fitting (0 |
|                       | indexed)                                            |
+-----------------------+-----------------------------------------------------+
| x_psf_err             | Uncertainty in x_psf                                |
+-----------------------+-----------------------------------------------------+
| y_psf                 | Row coordinate of the source from PSF fitting (0    |
|                       | indexed)                                            |
+-----------------------+-----------------------------------------------------+
| y_psf_err             | Uncertainty in y_psf                                |
+-----------------------+-----------------------------------------------------+
| ra                    | Best estimate of the right ascension (ICRS)         |
+-----------------------+-----------------------------------------------------+
| dec                   | Best estimate of the declination (ICRS)             |
+-----------------------+-----------------------------------------------------+
| ra_centroid           | Right ascension (ICRS) of the source centroid       |
+-----------------------+-----------------------------------------------------+
| dec_centroid          | Declination (ICRS) of the source centroid           |
+-----------------------+-----------------------------------------------------+
| ra_centroid_err       | Uncertainty in ra_centroid                          |
+-----------------------+-----------------------------------------------------+
| dec_centroid_err      | Uncertainty in dec_centroid                         |
+-----------------------+-----------------------------------------------------+
| ra_centroid_win       | Right ascension (ICRS) of the windowed source       |
|                       | centroid                                            |
+-----------------------+-----------------------------------------------------+
| dec_centroid_win      | Declination (ICRS) of the windowed source centroid  |
+-----------------------+-----------------------------------------------------+
| ra_centroid_win_err   | Uncertainty in ra_centroid_win                      |
+-----------------------+-----------------------------------------------------+
| dec_centroid_win_err  | Uncertainty in dec_centroid_win                     |
+-----------------------+-----------------------------------------------------+
| ra_psf                | Right ascension (ICRS) of the PSF-fitted position   |
+-----------------------+-----------------------------------------------------+
| dec_psf               | Declination (ICRS) of the PSF-fitted position       |
+-----------------------+-----------------------------------------------------+
| ra_psf_err            | Uncertainty in ra_psf                               |
+-----------------------+-----------------------------------------------------+
| dec_psf_err           | Uncertainty in dec_psf                              |
+-----------------------+-----------------------------------------------------+
| bbox_xmin             | Column index of the left edge of the source         |
|                       | bounding box (0 indexed)                            |
+-----------------------+-----------------------------------------------------+
| bbox_xmax             | Column index of the right edge of the source        |
|                       | bounding box (0 indexed)                            |
+-----------------------+-----------------------------------------------------+
| bbox_ymin             | Row index of the bottom edge of the source bounding |
|                       | box (0 indexed)                                     |
+-----------------------+-----------------------------------------------------+
| bbox_ymax             | Row index of the top edge of the source bounding    |
|                       | box (0 indexed)                                     |
+-----------------------+-----------------------------------------------------+
| segment_area          | Area of the source segment                          |
+-----------------------+-----------------------------------------------------+
| semimajor             | Length of the source semimajor axis computed from   |
|                       | image moments                                       |
+-----------------------+-----------------------------------------------------+
| semiminor             | Length of the source semiminor axis computed from   |
|                       | image moments                                       |
+-----------------------+-----------------------------------------------------+
| fwhm                  | Circularized full width at half maximum (FWHM)      |
|                       | calculated from the semimajor and semiminor axes as |
|                       | 2*sqrt(ln(2) * (semimajor**2 + semiminor**2))       |
+-----------------------+-----------------------------------------------------+
| ellipticity           | Source ellipticity as 1 - (semimajor / semiminor)   |
+-----------------------+-----------------------------------------------------+
| orientation_pix       | Angle measured counter-clockwise from the positive  |
|                       | X axis to the source major axis computed from image |
|                       | moments                                             |
+-----------------------+-----------------------------------------------------+
| orientation_sky       | Position angle from North of the source major axis  |
|                       | computed from image moments                         |
+-----------------------+-----------------------------------------------------+
| cxx                   | Coefficient for the x**2 term in the generalized    |
|                       | quadratic ellipse equation                          |
+-----------------------+-----------------------------------------------------+
| cxy                   | Coefficient for the x*y term in the generalized     |
|                       | quadratic ellipse equation                          |
+-----------------------+-----------------------------------------------------+
| cyy                   | Coefficient for the y**2 term in the generalized    |
|                       | quadratic ellipse equation                          |
+-----------------------+-----------------------------------------------------+
| kron_radius           | Unscaled first-moment Kron radius                   |
+-----------------------+-----------------------------------------------------+
| nn_label              | Segment label of the nearest neighbor in this       |
|                       | skycell                                             |
+-----------------------+-----------------------------------------------------+
| nn_distance           | Distance to the nearest neighbor in this skycell    |
+-----------------------+-----------------------------------------------------+
| sharpness             | Photutils DAOStarFinder sharpness statistic         |
+-----------------------+-----------------------------------------------------+
| roundness1            | Photutils DAOStarFinder roundness1 statistic        |
+-----------------------+-----------------------------------------------------+
| is_extended           | Flag indicating that the source appears to be more  |
|                       | extended than a point source                        |
+-----------------------+-----------------------------------------------------+
| fluxfrac_radius_50    | Radius of a circle centered on the source centroid  |
|                       | that encloses 50% of the Kron flux                  |
+-----------------------+-----------------------------------------------------+
| aper_bkg_flux         | Local background estimated within a circular        |
|                       | annulus                                             |
+-----------------------+-----------------------------------------------------+
| aper_bkg_flux_err     | Uncertainty in aper_bkg_flux                        |
+-----------------------+-----------------------------------------------------+
| aper01_flux           | Flux within a circular aperture of radius=0.1       |
|                       | arcsec                                              |
+-----------------------+-----------------------------------------------------+
| aper01_flux_err       | Uncertainty in aper01_flux                          |
+-----------------------+-----------------------------------------------------+
| aper02_flux           | Flux within a circular aperture of radius=0.2       |
|                       | arcsec                                              |
+-----------------------+-----------------------------------------------------+
| aper02_flux_err       | Uncertainty in aper02_flux                          |
+-----------------------+-----------------------------------------------------+
| aper04_flux           | Flux within a circular aperture of radius=0.4       |
|                       | arcsec                                              |
+-----------------------+-----------------------------------------------------+
| aper04_flux_err       | Uncertainty in aper04_flux                          |
+-----------------------+-----------------------------------------------------+
| aper08_flux           | Flux within a circular aperture of radius=0.8       |
|                       | arcsec                                              |
+-----------------------+-----------------------------------------------------+
| aper08_flux_err       | Uncertainty in aper08_flux                          |
+-----------------------+-----------------------------------------------------+
| aper16_flux           | Flux within a circular aperture of radius=1.6       |
|                       | arcsec                                              |
+-----------------------+-----------------------------------------------------+
| aper16_flux_err       | Uncertainty in aper16_flux                          |
+-----------------------+-----------------------------------------------------+
| psf_flux              | Total PSF flux                                      |
+-----------------------+-----------------------------------------------------+
| psf_flux_err          | Uncertainty in psf_flux                             |
+-----------------------+-----------------------------------------------------+
| segment_flux          | Isophotal flux                                      |
+-----------------------+-----------------------------------------------------+
| segment_flux_err      | Uncertainty in segment_flux                         |
+-----------------------+-----------------------------------------------------+
| kron_flux             | Flux within the elliptical Kron aperture            |
+-----------------------+-----------------------------------------------------+
| kron_flux_err         | Uncertainty in kron_flux                            |
+-----------------------+-----------------------------------------------------+
| kron_abmag            | AB magnitude within the elliptical Kron aperture    |
+-----------------------+-----------------------------------------------------+
| kron_abmag_err        | Uncertainty in kron_abmag                           |
+-----------------------+-----------------------------------------------------+
| warning_flags         | Warning bit flags (0 = good)                        |
+-----------------------+-----------------------------------------------------+
| image_flags           | Image quality bit flags (0 = good)                  |
+-----------------------+-----------------------------------------------------+
| psf_flags             | PSF fitting bit flags (0 = good)                    |
+-----------------------+-----------------------------------------------------+
| psf_gof               | PSF goodness of fit metric                          |
+-----------------------+-----------------------------------------------------+


Detailed descriptions of many of the columns can be found in the
`Photutils documentation <https://photutils.readthedocs.io/en/latest/>`_:

* `SourceCatalog
  <https://photutils.readthedocs.io/en/latest/api/photutils.segmentation.SourceCatalog.html>`_

* `PSFPhotometry
  <https://photutils.readthedocs.io/en/latest/api/photutils.psf.PSFPhotometry.html>`_

* `DAOStarFinder
  <https://photutils.readthedocs.io/en/latest/api/photutils.detection.DAOStarFinder.html>`_

Further details for some of the columns are provided below.

``flagged_spatial_index`` is a bit flag encoding the overlap flag,
projection, skycell, and pixel coordinates of the source. From high to
low, bit 64 is 1 if the object was outside of the core region of this
skycell or projection region. There is likely to be a better measurement
of the object in a different skycell with this bit set to 0. This bit
is the same as bit **TBD** of ``warning_flags``. Bits 49-63 encode the
primary projection region for this object. Bits 33-40 and 41-48 encode
the (x, y) skycell indices within this projection region, starting from
(0, 0) at the lower left. Bits 1-15 and 16-31 encode the x & y pixel
coordinate of the object within this skycell in virtual 0.05" pixels
(regardless of the pixel scale of the skycell).


Flag Columns
^^^^^^^^^^^^

The ``warning_flags`` column contains the following bit flags:

- 0 : good
- 1 :

  * Level 2: sources whose rounded centroid pixel is not finite or has
    DO_NOT_USE set in the model DQ

  * Level 3: sources whose rounded centroid pixel is not finite or has a
    weight of 0

The ``image_flags`` column contains the following bit flags:

- 0 : good
- 1 : one or more pixels in the source segment was flagged

The ``psf_flags`` column contains the following bit flags defined by the
:py:class:`photutils.psf.PSFPhotometry` class:

- 0 : good
- 1 : one or more pixels in the ``fit_shape`` region were masked
- 2 : the fit x and/or y position lies outside of the input data
- 4 : the fit flux is less than or equal to zero
- 8 : the fitter may not have converged
- 16 : the fitter parameter covariance matrix was not returned
- 32 : the fit x or y position is at the bounded value


Output Products
---------------

Source Catalog Table
^^^^^^^^^^^^^^^^^^^^

The output source catalog table is saved to a file in the `Parquet
<https://parquet.apache.org/>`_ format.


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
