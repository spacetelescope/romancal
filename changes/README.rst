Writing news fragments for the change log
#########################################

This ``changes/`` directory contains "news fragments": small ReStructured Text files describing a change in a few sentences.
When making a release, run ``towncrier build --version <VERSION>`` to consume existing fragments in ``changes/`` and insert them as a full change log entry at the top of ``CHANGES.rst`` for the released version.

News fragment filenames consist of the pull request number and the change log category (see below). A single change can have more than one news fragment, if it spans multiple categories:

.. code-block::

  2071.skycell.rst
  2077.associations.rst
  2089.skycell.rst
  2128.other.rst
  2034.docs.rst

Change log categories
*********************
- ``<PR#>.breaking.rst``: Also add this fragment if your change breaks existing functionality

General Pipeline Changes
========================

- ``changes/<PR#>.stpipe.rst``
- ``changes/<PR#>.associations.rst``
- ``changes/<PR#>.scripts.rst``
- ``changes/<PR#>.mosaic_pipeline.rst``
- ``changes/<PR#>.skycell.rst``

Step Changes
============

- ``changes/<PR#>.dq_init.rst``
- ``changes/<PR#>.saturation.rst``
- ``changes/<PR#>.refpix.rst``
- ``changes/<PR#>.wfi18_transient.rst``
- ``changes/<PR#>.linearity.rst``
- ``changes/<PR#>.dark_current.rst``
- ``changes/<PR#>.jump_detection.rst``
- ``changes/<PR#>.ramp_fitting.rst``
- ``changes/<PR#>.assign_wcs.rst``
- ``changes/<PR#>.flatfield.rst``
- ``changes/<PR#>.photom.rst``
- ``changes/<PR#>.flux.rst``
- ``changes/<PR#>.source_detection.rst``
- ``changes/<PR#>.tweakreg.rst``
- ``changes/<PR#>.skymatch.rst``
- ``changes/<PR#>.outlier_detection.rst``
- ``changes/<PR#>.resample.rst``
- ``changes/<PR#>.source_catalog.rst``

Other Changes
=============

- ``changes/<PR#>.other.rst``: infrastructure or miscellaneous change
- ``changes/<PR#>.docs.rst``

.. note:: This README was adapted from the Astropy changelog readme under the terms of BSD license, which in turn adapted from the Numpy changelog readme under the terms of the MIT licence.
