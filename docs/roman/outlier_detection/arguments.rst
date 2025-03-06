.. _outlier_detection_step_args:

For more details about step arguments (including datatypes, possible values
and defaults) see :py:obj:`romancal.outlier_detection.OutlierDetectionStep.spec`.

Step Arguments
==============
The ``outlier_detection`` step has the following optional arguments that control the
behavior of the processing:

``--weight_type``
  The type of data weighting to use during resampling the images for creating the
  median image used for detecting outliers; options are `'ivm'`, `'exptime'`,
  and `None` (see :ref:`weight_type_options_details_section` for details).

``--pixfrac``
  Fraction by which input pixels are “shrunk” before being drizzled onto the output
  image grid, given as a real number between 0 and 1. This specifies the size of the
  footprint, or “dropsize”, of a pixel in units of the input pixel size. If `pixfrac`
  is set to less than 0.001, the `kernel` parameter will be reset to `'point'`` for more
  efficient processing. In the step of drizzling each input image onto a separate
  output image, the default value of 1.0 is best in order to ensure that each
  output drizzled image is fully populated with pixels from the input image.
  Valid values range from 0.0 to 1.0.

``--kernel``
  This parameter specifies the form of the kernel function used to distribute
  flux onto the separate output images, for the initial separate drizzling
  operation only. The value options for this parameter include:

      * `'square'`: original classic drizzling kernel

      * `'tophat'`: this kernel is a circular "top hat" shape of width
        `pixfrac`. It effects only output pixels within a radius of
        `pixfrac/2` from the output position.

      * `'lanczos3'`: a Lanczos style kernel, extending a radius of
        3 pixels from the center of the detection. The Lanczos kernel is
        a damped and bounded form of the "sinc" interpolator, and is very
        effective for resampling single images when ``scale=pixfrac=1``.
        It leads to less resolution loss than other kernels, and typically
        results in reduced correlated noise in outputs.

      .. warning:: The `'lanczos3'` kernel tends to result in much slower
         processing as compared to other kernel options. This option
         should never be used for ``pixfrac != 1.0``, and is not recommended
         for ``scale!=1.0``.

``--fillval``
    The value for this parameter is to be assigned to the output pixels that
    have zero weight or which do not receive flux from any input pixels during
    drizzling. This parameter corresponds to the ``fillval`` parameter of the
    `drizzle` task and will be converted to a float.

``--maskpt``
  Percentage of weight image values below which they are flagged as bad and rejected
  from the median image. Valid values range from 0.0 to 1.0.

``--snr``
  The signal-to-noise values to use for bad pixel identification. Since cosmic rays
  often extend across several pixels the user must specify two cut-off values for
  determining whether a pixel should be masked: the first for detecting the primary
  cosmic ray, and the second (typically lower threshold) for masking lower-level bad
  pixels adjacent to those found in the first pass.  Valid values are a pair of
  floating-point values in a single string (for example "5.0 4.0").

``--scale``
  The scaling factor applied to derivative used to identify bad pixels. Since cosmic
  rays often extend across several pixels the user must specify two cut-off values for
  determining whether a pixel should be masked: the first for detecting the primary
  cosmic ray, and the second (typically lower threshold) for masking lower-level bad
  pixels adjacent to those found in the first pass.  Valid values are a pair of
  floating-point values in a single string (for example "1.2 0.7").

``--backg``
  User-specified background value (scalar) to subtract during final identification
  step of outliers in `driz_cr` computation.

``--save_intermediate_results``
  Boolean specifying whether or not to write out intermediate products such as median image or
  resampled individual input exposures to disk. Typically, only used to track down
  problems with final results when too many or too few pixels are flagged as outliers.

``--resample_data``
  Boolean specifying whether or not to resample the input images when performing outlier
  detection.

``--resample_on_skycell``
  If input association contains skycell information use it for the resampling
  wcs. If `False` (or if the association contains no skycell information) the
  resampled wcs will be the combined wcs of all input models.

``--good_bits``
  The DQ bit values from the input image DQ arrays that should be considered 'good'
  when creating masks of bad pixels during outlier detection when resampling the data.
  See `Roman's Data Quality Flags
  <https://github.com/spacetelescope/romancal/blob/main/romancal/lib/dqflags.py>`_
  for details.

``--in_memory``
  Boolean specifying whether or not to keep all intermediate products and datamodels in
  memory at the same time during the processing of this step.  If set to `False`,
  any `ModelLibrary` opened by this step will use ``on_disk=True`` and use temporary
  files to store model modifications. Additionally any resampled images will
  be kept in memory (as long as needed). This can result in much lower memory
  usage (at the expense of file I/O) to process large associations.

.. _weight_type_options_details_section:

Weighting types
===============
``weight_type`` specifies the type of weighting image to apply with the bad pixel
mask for the final drizzle step.  The options for this parameter include:

* `ivm`: allows the user to either supply their own inverse-variance weighting map,
  or allow ``drizzle`` to generate one automatically on-the-fly during the final
  drizzle step. This parameter option may be necessary for specific purposes.
  For example, to create a drizzled weight file for software such as ``SExtractor``,
  it is expected that a weight image containing all of the background noise sources
  (sky level, read-noise, dark current, etc), but not the Poisson noise from the
  objects themselves will be available. The user can create the inverse variance
  images and then specify their names using the ``input`` parameter for ``drizzle``
  to specify an '\@file'. This would be a single ``ASCII`` file containing the list
  of input calibrated exposure filenames (one per line), with a second column
  containing the name of the ``IVM`` file corresponding to each calibrated exposure.
  Each ``IVM`` file must have the same file format as the input file.

* `exptime`: the images will be weighted according to their exposure time, which is the
  standard behavior for drizzle. This weighting is a good approximation in the regime
  where the noise is dominated by photon counts from the sources, while contributions
  from sky background, read-noise and dark current are negligible. This option is
  provided as the default since it produces reliable weighting for all types of data.

* ``None``: In this case, a bit mask will be generated based on the DQ array and a
  bit flag set to 0 (i.e. `GOOD`; see `Roman's Data Quality Flags
  <https://github.com/spacetelescope/romancal/blob/main/romancal/lib/dqflags.py>`_
  for details).
