Changelog
=========

This directory contains "news fragments" which are short files that contain a
small **ReST**-formatted text that will be added to the full changelog.

Make sure to use full sentences with correct case and punctuation.

Consuming news fragments in `changes/` into a new change log entry
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Running `towncrier build` will read all existing fragment files in `changes/`
and create a new entry at the top of `CHANGES.rst` with the specified version number.

```shell
pip install towncrier
towncrier build --version <VERSION>
```

News fragment change types
--------------------------

- ``<PR#>.breaking.rst``: Also add a fragment of this type if your change breaks existing functionality

General Pipeline Changes
""""""""""""""""""""""""

- ``changes/<PR#>.stpipe.rst``
- ``changes/<PR#>.associations.rst``
- ``changes/<PR#>.scripts.rst``
- ``changes/<PR#>.mosaic_pipeline.rst``
- ``changes/<PR#>.skycell.rst``

Step Changes
""""""""""""

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
"""""""""""""

- ``changes/<PR#>.other.rst``: infrastructure or miscellaneous change
- ``changes/<PR#>.docs.rst``

Note
----

This README was adapted from the Astropy changelog readme under the terms
of BSD license, which in turn adapted from the Numpy changelog readme under
the terms of the MIT licence.
