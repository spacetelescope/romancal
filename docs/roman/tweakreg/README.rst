Description
===========

:Class: `roman.tweakreg.TweakRegStep`
:Alias: tweakreg

Overview
--------
This step uses the coordinates of point-like sources from an input catalog
(i.e. the result from `SourceCatalogStep` saved in the
`meta.source_catalog.tweakreg_catalog_name` attribute) and compares them with the
coordinates from a Gaia catalog to compute corrections to
the WCS of the input images such that sky catalogs obtained from the image catalogs
using the corrected WCS will align on the sky.

Custom Source Catalogs
----------------------
Custom source catalogs can be supplied through
``meta.source_catalog.tweakreg_catalog`` (in-memory table) or
``meta.source_catalog.tweakreg_catalog_name`` (file path). The catalog must
be in a format supported by :py:meth:`~astropy.table.Table.read` (or parquet
for ``.parquet`` files). The catalog must contain
either ``'x'`` and ``'y'`` or ``'x_psf'`` and ``'y_psf'`` columns which
indicate source *image* coordinates (in pixels). Pixel coordinates are
0-indexed.

Association files can also be used as ``tweakreg`` input for custom catalogs.
When an association is provided, ``tweakreg`` reads the custom catalog
information for each member from that member's ``tweakreg_catalog`` attribute
and sets it as the value for that member's
``meta.source_catalog.tweakreg_catalog_name`` metadata.
For example, the following association contains two members
(``image1.asdf`` and ``image2.asdf``) with custom catalogs and one member
(``image3.asdf``) without a custom catalog:

  .. code-block:: json

    {
      "asn_type": "tweakreg",
      "asn_id": "tweakreg_12345678",
      "members": [
        {
          "expname": "image1.asdf",
          "tweakreg_catalog": "/path/to/image1_catalog.parquet"
        },
        {
          "expname": "image2.asdf",
          "tweakreg_catalog": "/path/to/image2_catalog.parquet"
        },
        {
          "expname": "image3.asdf",
        }
      ]
    }

In this case, ``tweakreg`` will read the custom catalogs for ``image1.asdf`` and
``image2.asdf`` from the specified file paths and use them for alignment, while
it will attempt to read the source catalog for ``image3.asdf`` from the file path
specified in its ``meta.source_catalog.tweakreg_catalog_name`` metadata
(which is expected to be set by a previous step such as `SourceCatalogStep`).

.. note::
    ``tweakreg`` requires ``meta.source_catalog`` to be present.

Alignment
---------
The source catalog (either created by `SourceCatalogStep` or provided by the user)
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


* ``catalog_format``: A `str` indicating one of the catalog output file format
  supported by :py:class:`astropy.table.Table` (Default='ascii.ecsv').

  .. note::
    - The full list of supported formats can be found on
      the `astropy`'s `Built-In Table Readers/Writers`_ webpage.

.. _`Built-In Table Readers/Writers`: https://docs.astropy.org/en/stable/io/unified.html#built-in-table-readers-writers


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
  - ``'rscale'``: rotation, shifts, and scale
  - ``'general'``: rotation, shifts, scale, and skew

  The default value is "rshift".

  .. note::
      Mathematically, alignment of images observed in different tangent planes
      requires ``fitgeometry='general'`` in order to fit source catalogs
      in the different images even if misalignment is caused only by a shift
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
  (e.g. 10 times smaller) (Default=1.0).

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

* ``update_source_catalog_coordinates``: A boolean indicating whether to update
  the source catalog coordinates after applying the WCS corrections (Default=False).

  .. note::
    If `True`, the source catalog coordinates will be updated to reflect
    the new positions based on the corrected WCS of the input image.

Further Documentation
---------------------
The underlying algorithms as well as formats of source catalogs are described
in more detail on the TweakWCS_ webpage.

.. _TweakWCS: https://tweakwcs.readthedocs.io/en/latest/
.. _linearfit: https://tweakwcs.readthedocs.io/en/latest/source/linearfit.html
    #tweakwcs.linearfit.iter_linear_fit
