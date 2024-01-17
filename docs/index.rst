.. romancal documentation master file
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

:ref:`genindex`  |  :ref:`modindex`

==============================================
The Roman Space Telescope Calibration Pipeline
==============================================

.. image:: _static/roman_logo_black_w200px.png
   :align: center
   :alt: Nancy Roman Space Telescope

Welcome to the documentation for the Roman calibration software,
`romancal <https://github.com/spacetelescope/romancal>`__.
This package contains the Python software suite for the
Roman Space Telescope (RST) calibration pipeline, which processes data
from the Roman Wide-Field Instrument (WFI) by applying various corrections
to produce science-ready, calibrated output products including fully calibrated
individual exposures as well as high-level data products (mosaics,
catalogs, etc.). The tools in this package allow users to run and
configure the pipeline to custom process their Roman data.
Additionally, the romancal package contains the interface to
Roman datamodels, the recommended method of reading and writing
Roman data files in Python.


If you have questions or concerns regarding the software, please contact the Roman Help
desk at `Roman Help Desk <https://stsci.service-now.com/roman>`_.

--------------------------------

.. include:: roman/introduction.rst

.. toctree::
   :caption: RomanCal Pipeline
   :maxdepth: 2

   roman/pipeline_installation.rst
   roman/pipeline_levels.rst
   roman/pipeline_run.rst
   roman/pipeline_steps.rst
   roman/pipeline_ref_files.rst
   roman/pipeline_parameters.rst
   roman/pipeline_naming_conventions.rst
   roman/error_propagation/index.rst

.. toctree::
   :maxdepth: 2
   :caption: Data Products Documentation

   roman/data_products/index.rst

.. toctree::
   :caption: RomanCal Package Index
   :maxdepth: 3

   roman/package_index.rst

.. toctree::
   :caption: Datamodels
   :maxdepth: 3

   roman/datamodels/index.rst

.. toctree::
   :caption: Additional Information
   :maxdepth: 1

   roman/pipeline_static_preview.rst
