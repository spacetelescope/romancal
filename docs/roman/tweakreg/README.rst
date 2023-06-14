Description
===========

:Class: `roman.tweakreg.TweakRegStep`
:Alias: tweakreg

Overview
--------
This step uses the coordinates of point-like sources from an input catalog
(i.e. the result from `SourceDetectionStep` saved in the
`meta.tweakreg_catalog` attribute) and compares them with the
coordinates from a Gaia catalog to compute corrections to
the WCS of the input images such that sky catalogs obtained from the image catalogs
using the corrected WCS will align on the sky.

Custom Source Catalogs
----------------------
The default catalog used by ``tweakreg`` step can be disabled by
providing a file name to a custom source catalog in the
``meta.tweakreg_catalog`` attribute of input data models.
The catalog must be in a format automatically recognized by
:py:meth:`~astropy.table.Table.read`. The catalog must contain
either ``'x'`` and ``'y'`` or ``'xcentroid'`` and ``'ycentroid'`` columns which
indicate source *image* coordinates (in pixels). Pixel coordinates are
0-indexed.

For the ``tweakreg`` step to use user-provided input source catalogs,
``use_custom_catalogs`` parameter of the ``tweakreg`` step must be set to
`True`.

In addition to setting the ``meta.tweakreg_catalog`` attribute of input data
models to the custom catalog file name, the ``tweakreg_step`` also supports two
other ways of supplying custom source catalogs to the step:

1. Adding ``tweakreg_catalog`` attribute to the ``members`` of the input ASN
   table - see `~roman.datamodels.ModelContainer` for more details.
   Catalog file names are relative to ASN file path.

2. Providing a simple two-column text file, specified via step's parameter
   ``catfile``, that contains input data models' file names in the first column
   and the file names of the corresponding catalogs in the second column.
   Catalog file names are relative to ``catfile`` file path.

Specifying custom source catalogs via either the input ASN table or
``catfile``, will update input data models' ``meta.tweakreg_catalog``
attributes to the catalog file names provided in either in the ASN table or
``catfile``.

.. note::
    When custom source catalogs are provided via both ``catfile`` and
    ASN table members' attributes, the ``catfile`` takes precedence and
    catalogs specified via ASN table are ignored altogether.

.. note::
    1. Providing a data model file name in the ``catfile`` and leaving
       the corresponding source catalog file name empty -- same as setting
       ``'tweakreg_catalog'`` in the ASN table to an empty string ``""`` --
       would set corresponding input data model's ``meta.tweakreg_catalog``
       attribute to `None`. In this case, ``tweakreg_step`` will automatically
       generate a source catalog for that data model.

    2. If an input data model is not listed in the ``catfile`` or does not
       have ``'tweakreg_catalog'`` attribute provided in the ASN table,
       then the catalog file name in that model's ``meta.tweakreg_catalog``
       attribute will be used. If ``model.meta.tweakreg_catalog`` is `None`,
       ``tweakreg_step`` will automatically generate a source catalog for
       that data model.

Alignment
---------
The source catalog (either created by `SourceDetectionStep` or provided by the user)
gets cross-matched and fit to an astrometric reference catalog
(set by ``TweakRegStep.abs_refcat``) and the results are stored in
``model.meta.wcs_fit_results``. The pipeline initially supports fitting to any
Gaia Data Release (defaults to `GAIADR3`).

An example of the content of ``model.meta.wcs_fit_results`` is as follows:

  .. code-block:: python

          model.meta.wcs_fit_results = {
            "status": "SUCCESS",
            "fitgeom": "rshift",
            "matrix": array([[ 1.00000000e+00,  1.04301609e-13],
                  [-1.04301609e-13,  1.00000000e+00]]),
            "shift": array([ 7.45523163e-11, -1.42718944e-10]),
            "center": array([-183.87997841, -119.38467775]),
            "proper_rot": 5.9760419875149846e-12,
            "proper": True,
            "rot": (5.9760419875149846e-12, 5.9760419875149846e-12),
            "<rot>": 5.9760419875149846e-12,
            "scale": (1.0, 1.0),
            "<scale>": 1.0,
            "skew": 0.0,
            "rmse": 2.854152848489525e-10,
            "mae": 2.3250544963289652e-10,
            "nmatches": 22
          }

Details about most of the parameters available in ``model.meta.wcs_fit_results`` can be
found on the TweakWCS_'s webpage, under its linearfit_ module.



WCS Correction
--------------
The linear coordinate transformation computed in the previous step
is used to define tangent-plane corrections that need to be applied
to the GWCS pipeline in order to correct input image WCS.
This correction is implemented by inserting a ``v2v3corr`` frame with
tangent plane corrections into the GWCS pipeline of the image's WCS.

Step Arguments
--------------
``TweakRegStep`` has the following arguments:

**Catalog parameters:**

* ``use_custom_catalogs``: A boolean that indicates whether
  to ignore source catalog in the input data model's ``meta.tweakreg_catalog``
  attribute (Default=`False`).

  .. note::
    If `True`, the user must provide a valid custom catalog that will be assigned to
    `meta.tweakreg_catalog` and used throughout the step.

* ``catalog_format``: A `str` indicating one of the catalog output file format
  supported by :py:class:`astropy.table.Table` (Default='ascii.ecsv').

  .. note::
    - This option must be provided whenever `use_custom_catalogs = True`.

    - The full list of supported formats can be found on
      the `astropy`'s `Built-In Table Readers/Writers`_ webpage.

.. _`Built-In Table Readers/Writers`: https://docs.astropy.org/en/stable/io/unified.html#built-in-table-readers-writers

* ``catfile``: Name of the file with a list of custom user-provided catalogs
  (Default='').

  .. note::
    - This option must be provided whenever `use_custom_catalogs = True`.

* ``catalog_path``: A `str` indicating the catalogs output file path (Default='').

  .. note::
      All catalogs will be saved to this path.
      The default value is the current working directory.

**Reference Catalog parameters:**

* ``expand_refcat``: A boolean indicating whether or not to expand reference
  catalog with new sources from other input images that have been already
  aligned to the reference image (Default=False).

**Object matching parameters:**

* ``minobj``: A positive `int` indicating minimum number of objects acceptable
  for matching (Default=15).

* ``searchrad``: A `float` indicating the search radius in arcsec for a match
  (Default=2.0).

* ``use2dhist``: A boolean indicating whether to use 2D histogram to find
  initial offset (Default=True).

* ``separation``: Minimum object separation in arcsec (Default=1.0).

* ``tolerance``: Matching tolerance for ``xyxymatch`` in arcsec (Default=0.7).

**Catalog fitting parameters:**

* ``fitgeometry``: A `str` value indicating the type of affine transformation
  to be considered when fitting catalogs. Allowed values:

  - ``'shift'``: x/y shifts only
  - ``'rshift'``: rotation and shifts
  - ``'rscale'``: rotation and scale
  - ``'general'``: shift, rotation, and scale

  The default value is "rshift".

  .. note::
      Mathematically, alignment of images observed in different tangent planes
      requires ``fitgeometry='general'`` in order to fit source catalogs
      in the different images even if mis-alignment is caused only by a shift
      or rotation in the tangent plane of one of the images.

      However, under certain circumstances, such as small alignment errors or
      minimal dithering during observations that keep tangent planes of the
      images to be aligned almost parallel, then it may be more robust to
      use a ``fitgeometry`` setting with fewer degrees of freedom such as
      ``'rshift'``, especially for "ill-conditioned" source catalogs such as
      catalogs with very few sources, or large errors in source positions, or
      sources placed along a line or bunched in a corner of the image (not
      spread across/covering the entire image).

* ``nclip``: A non-negative integer number of clipping iterations
  to use in the fit (Default = 3).

* ``sigma``: A positive `float` indicating the clipping limit, in sigma units,
  used when performing fit (Default=3.0).

**Absolute Astrometric fitting parameters:**

Parameters used for absolute astrometry to a reference catalog.

* ``abs_refcat``: String indicating what astrometric catalog should be used.
  Currently supported options are (Default='GAIADR3'): ``'GAIADR1'``, ``'GAIADR2'``,
  or ``'GAIADR3'``.

  .. note::
    If `None` or an empty string is passed in, `TweakRegStep` will
    use the default catalog as set by `tweakreg_step.DEFAULT_ABS_REFCAT`.

* ``abs_minobj``: A positive `int` indicating minimum number of objects
  acceptable for matching. (Default=15).

* ``abs_searchrad``: A `float` indicating the search radius in arcsec for
  a match. It is recommended that a value larger than ``searchrad`` be used for
  this parameter (e.g. 3 times larger) (Default=6.0).

* ``abs_use2dhist``: A boolean indicating whether to use 2D histogram to find
  initial offset. It is strongly recommended setting this parameter to `True`.
  Otherwise the initial guess for the offsets will be set to zero
  (Default=True).

* ``abs_separation``: Minimum object separation in arcsec. It is recommended
  that a value smaller than ``separation`` be used for this parameter
  (e.g. 10 times smaller) (Default=0.1).

* ``abs_tolerance``: Matching tolerance for ``xyxymatch`` in arcsec (Default=0.7).

* ``abs_fitgeometry``: A `str` value indicating the type of affine
  transformation to be considered when fitting catalogs. Allowed values:

  - ``'shift'``: x/y shifts only
  - ``'rshift'``: rotation and shifts
  - ``'rscale'``: rotation and scale
  - ``'general'``: shift, rotation, and scale

  The default value is "rshift". Note that the same conditions/restrictions
  that apply to ``fitgeometry`` also apply to ``abs_fitgeometry``.

* ``abs_nclip``: A non-negative integer number of clipping iterations
  to use in the fit (Default = 3).

* ``abs_sigma``: A positive `float` indicating the clipping limit, in sigma
  units, used when performing fit (Default=3.0).

* ``save_abs_catalog``: A boolean specifying whether or not to write out the
  astrometric catalog used for the fit as a separate product (Default=False).


Further Documentation
---------------------
The underlying algorithms as well as formats of source catalogs are described
in more detail on the TweakWCS_ webpage.

.. _TweakWCS: https://tweakwcs.readthedocs.io/en/latest/
.. _linearfit: https://tweakwcs.readthedocs.io/en/latest/source/linearfit.html
    #tweakwcs.linearfit.iter_linear_fit
