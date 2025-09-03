.. _resample_step_args:

For more details about step arguments (including datatypes, possible values
and defaults) see :py:obj:`romancal.resample.ResampleStep.spec`.

Step Arguments
==============
The ``resample`` step has the following optional arguments that control
the behavior of the processing and the characteristics of the resampled
image.

``--pixfrac``
    The fraction by which input pixels are "shrunk" before being drizzled
    onto the output image grid, given as a real number between 0 and 1.

``--kernel``
    The form of the kernel function used to distribute flux onto the output
    image.  Available kernels are `square`, `gaussian`, `point`, `tophat`,
    `turbo`, `lanczos2`, and `lanczos3`.

``--pixel_scale_ratio``
    Ratio of input to output pixel scale.  A value of 0.5 means the output
    image would have 4 pixels sampling each input pixel.
    Ignored when ``pixel_scale`` or ``output_wcs`` are provided.

``--pixel_scale``
    Absolute pixel scale in ``degrees``. When provided, overrides
    ``pixel_scale_ratio``. Ignored when ``output_wcs`` is provided.

``--rotation``
    Position angle (in degrees) of output image’s Y-axis relative to North.
    A value of 0.0 would orient the final output image to be North up.
    The value of `None` specifies that the images will not be rotated,
    but will instead be resampled in the default orientation for the camera
    with the x and y axes of the resampled image corresponding
    approximately to the detector axes. Ignored when ``pixel_scale``
    or ``output_wcs`` are provided.

``--crpix``
    Position of the reference pixel in the image array in the ``x, y`` order.
    If ``crpix`` is not specified, it will be set to the center of the bounding
    box of the returned WCS object. When supplied from command line, it should
    be a comma-separated list of floats. Ignored when ``output_wcs``
    is provided.

``--crval``
    Right ascension and declination of the reference pixel. Automatically
    computed if not provided. When supplied from command line, it should be a
    comma-separated list of floats. Ignored when ``output_wcs`` is provided.

``--output_shape``
    Shape of the image (data array) using "standard" ``nx`` first and ``ny``
    second (as opposite to the ``numpy.ndarray`` convention - ``ny`` first and
    ``nx`` second). This value will be assigned to
    ``pixel_shape`` and ``array_shape`` properties of the returned
    WCS object. When supplied from command line, it should be a comma-separated
    list of integers ``nx, ny``.

    .. note::
        Specifying ``output_shape`` *is required* when the WCS in
        ``output_wcs`` does not have ``bounding_box`` property set.

``--output_wcs``
    File name of a ``ASDF`` file with a GWCS stored under the ``"wcs"`` key
    under the root of the file. The output image size is determined from the
    bounding box of the WCS (if any). Argument ``output_shape`` overrides
    computed image size and it is required when output WCS does not have
    ``bounding_box`` property set.

    .. note::
        When ``output_wcs`` is specified, WCS-related arguments such as
        ``pixel_scale_ratio``, ``pixel_scale``, ``rotation``, ``crpix``,
        and ``crval`` will be ignored.

``--resample_on_skycell``
    If input association contains skycell information use it to compute
    the output wcs. If ``output_wcs`` is defined it will be used instead.
    If ``resample_on_skycell`` is `False` the output wcs will be the combined
    wcs of all input models.

``--fillval``
    The value to assign to output pixels that have zero weight or do not
    receive any flux from any input pixels during drizzling.

``--weight_type``
    The weighting type for each input image.
    If `weight_type=ivm`, the scaling value
    will be determined per-pixel using the inverse of the read noise
    (VAR_RNOISE) array stored in each input image. If the VAR_RNOISE array does
    not exist, the variance is set to 1 for all pixels (equal weighting).
    If `weight_type=ivsky`, the scaling value
    will be determined per-pixel using the inverse of the sky variance
    (VAR_SKY) array stored in each input image. If the VAR_SKY array does
    not exist, the variance is set to 1 for all pixels (equal weighting).
    If `weight_type=exptime`, the scaling value will be set equal to the
    exposure time found in the image header.

    .. note::
        VAR_SKY is calculated as follows:

        VAR_SKY = (VAR_RNOISE + VAR_POISSON) / (DATA * MEDIAN(DATA))

        where VAR_RNOISE is the read noise variance, VAR_POISSON is the
        Poisson noise variance, DATA is the science data, and MEDIAN(DATA)
        is the median of the science data.

``--in_memory``
    If set to `False`, write output datamodel to disk.

``--good_bits``
    Specifies the bits to use when creating the resampling mask.
    Either a single bit value or a combination of them can be provided.
    If the string starts with a tilde (`~`), then the provided bit(s)
    will be excluded when creating the resampling mask.
    A value of ``~DO_NOT_USE+NON_SCIENCE`` will exclude pixels
    flagged with ``DO_NOT_USE`` and ``NON_SCIENCE``.

    The bit value can be provided in a few different ways, but always as
    a string type. For example, if the user deems OK to use pixels with
    low QE and highly nonlinear, then any of the ways listed below will
    work to set ``good_bits``:

    - ``good_bits = 'LOW_QE+NON_LINEAR'`` (concatenated DQ flag labels);
    - ``good_bits = '8192+65536'`` (concatenated DQ flag bit values);
    - ``good_bits = '8192,65536'`` (comma-separated DQ flag bit values);
    - ``good_bits = '73728'`` (sum of DQ flag bit values).

    .. note::
        Adding a tilde (`~`) to the beginning of the string will flip the
        bit(s) and actually `exclude` the provided bit(s). This is the same
        as providing the bad bits instead of the good bits.
