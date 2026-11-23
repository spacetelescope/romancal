.. _asn-roman-conventions:

=================
Roman Conventions
=================

.. _asn-roman-naming:

Association Naming
==================

Association files produced by ``romancal`` are JSON files whose names end in
``_asn.json``. Two association generators live in the package and define the
naming conventions below:

* :ref:`skycell_asn` — single-band, skycell-based Level 3 associations
* :ref:`multiband_asn` — multi-filter associations per skycell

Associations that group the 18 WFI detectors for exposure-level processing are
produced by SDP and are not controlled by ``romancal``; their naming is outside
the scope of this page.

Skycell associations
--------------------

``skycell_asn`` builds one association per unique combination of program root,
data release ID, product type, skycell, and filter. The association file name
follows::

  rPPPPP_<data_release_id>_<product_type>_<skycell>_<filter>_asn.json

where:

  * ``r``: All Roman-related products begin with ``r``
  * ``PPPPP``: 5-digit program number
  * ``data_release_id``: Data release identifier (default ``p`` for prompt;
    e.g. ``dr01`` for a numbered release). Override with
    ``--data-release-id``.
  * ``product_type``: Product granularity requested via ``--product-type``.
    Values used in the file name are:

    * ``visit`` → ``v`` followed by execution plan, pass, segment, observation,
      and visit IDs (for example ``v01001001003001``)
    * ``pass`` → ``p`` followed by execution plan and pass IDs
    * ``full`` → ``full``
    * ``daily`` → ``d`` followed by execution plan, pass, and segment IDs
    * ``user`` → ``user``

  * ``skycell``: Skycell name on the celestial sphere (for example
    ``270p65x48y69``)
  * ``filter``: Optical element in lower case (for example ``f158``, ``f146``)

Example (visit-level, prompt release)::

  r00991_p_v01001001003001_270p65x48y69_f158_asn.json

Another visit-level example::

  r01055_p_v02001001008001_249p38x31y66_f146_asn.json

See :ref:`skycell_asn` for command-line usage and input file selection.

Multiband associations
----------------------

``multiband_asn`` groups coadd products that share a skycell (and program /
release / product-type prefix) across filters. Multiband association file names
**omit the filter**, which distinguishes them from single-band skycell
associations::

  rPPPPP_<data_release_id>_<product_type>_<skycell>_asn.json

Example inputs::

  r00001_p_full_270p65x48y69_f184_coadd.asdf
  r00001_p_full_270p65x48y69_f213_coadd.asdf

produce::

  r00001_p_full_270p65x48y69_asn.json

See :ref:`multiband_asn` for further examples and the optional
``psf_match_reference_filter`` association key.

.. _asn-roman-association-types:

Association Types
=================

Each association carries an ``asn_type`` meta keyword that describes the kind of
science product the association is intended to produce. Pipeline steps and
tools use this value together with the association member list to decide how to
process the set.

Common values set by ``romancal`` tools include:

  * ``image``: Set by :ref:`skycell_asn` for single-band mosaic / Level 3 image
    associations.
  * ``user_built``: Default for associations created with
    :ref:`asn_from_list <asn-from-list>` when no type is supplied.

Other ``asn_type`` values may appear in associations constructed outside
``romancal`` (for example by SDP).

.. _asn-roman-field-guide:

Field Guide to File Names
=========================

The patterns below identify the main classes of Roman pipeline products.
They are not a full catalog of every suffix; individual step and pipeline pages
document additional product types. Exposure-level naming is also described in
:ref:`file_naming_schemes`.

Exposure-level (Level 1 / Level 2) products
-------------------------------------------

Calibrated and uncalibrated detector products use ASDF and follow::

  rPPPPPCCAAASSSOOOVVV_eeee_wfiDD[_<filter>]_<prodType>.asdf

where ``PPPPPCCAAASSSOOOVVV`` is the 19-digit visit identifier, ``eeee`` is the
exposure number, ``wfiDD`` is the detector (``wfi01`` … ``wfi18``), the optional
``filter`` token is the optical element in lower case (for example ``f146``),
and ``prodType`` is the product suffix (for example ``uncal`` or ``cal``).

A compact regular expression for these names is::

  r[0-9]{19}_[0-9]{4}_wfi[0-9]{2}(?:_f[0-9]+)?_[A-Za-z0-9]+\.asdf

Examples::

  r0000101001001001001_0001_wfi01_uncal.asdf
  r0105502001001008001_0069_wfi11_f146_cal.asdf

Intermediate optional step suffixes are listed under
:ref:`pipeline_step_suffix_definitions`.

Association files
-----------------

Association JSON files produced by ``romancal`` end in ``_asn.json`` and match
one of:

* Skycell (single-band)::

    r[0-9]{5}_.+_.+_.+_f[0-9]+_asn\.json

* Multiband (no filter token)::

    r[0-9]{5}_.+_.+_.+_asn\.json

  where the final underscore-delimited field before ``_asn.json`` is a skycell
  id rather than a filter.

Examples::

  r00991_p_v01001001003001_270p65x48y69_f158_asn.json
  r00001_p_full_270p65x48y69_asn.json

Level 3 / mosaic coadd products
-------------------------------

Resampled skycell products from high-level processing typically look like::

  rPPPPP_<data_release_id>_<product_type>_<skycell>_<filter>_coadd.asdf

Example::

  r00001_p_full_270p65x48y69_f184_coadd.asdf

A compact regular expression is::

  r[0-9]{5}_.+_.+_.+_f[0-9]+_coadd\.asdf

Related documentation
=====================

* :ref:`file_naming_schemes` — exposure file name fields
* :ref:`pipeline_step_suffix_definitions` — calibration step output suffixes
* :ref:`skycell_asn` — creating skycell associations
* :ref:`multiband_asn` — creating multiband associations
* :ref:`association-commands` — association command reference
