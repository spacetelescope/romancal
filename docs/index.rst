:ref:`genindex`  |  :ref:`modindex`

.. image:: _static/stsci_logo.png
   :width: 15%
   :alt: STScI Logo
   :target: https://stsci.edu

.. image:: _static/stsci_name.png
   :width: 68%
   :alt: STScI Logo
   :target: https://stsci.edu

.. image:: _static/A_Rst_1020_Logo_Box_Color_Seafoamondarkblue_Logo_Full_Color_RGB_850px@72ppi.png
   :width: 15%
   :alt: Nancy Grace Roman Space Telescope
   :target: https://science.nasa.gov/mission/roman-space-telescope/

.. _roman-pipeline-doc-index:

==============================================
The Roman Space Telescope Calibration Pipeline
==============================================

**Version**: |release|

This package (``romancal``) processes uncalibrated data from the `Nancy Grace Roman Space Telescope (Roman) <https://science.nasa.gov/mission/roman-space-telescope/>`_.
The pipeline performs a series of calibration steps that result in standard data products,
applying various corrections to produce science-ready, calibrated output products including
individual exposures and high-level data products (mosaics, catalogs etc.).

`See README.md for installation and usage instructions <https://github.com/spacetelescope/romancal?tab=readme-ov-file#installation>`_.

This package allows users to run and configure the calibration pipeline themselves for custom processing of Roman Telescope data,
either :ref:`from the command line <running-the-pipeline-from-command-line>` with ``strun``
or from Python with :ref:`pipeline and step functions and classes <running-the-pipeline-from-python>` in the ``romancal`` package.
Additionally, the ``romancal`` package provides :ref:`Roman Telescope datamodel classes <datamodels>`,
the recommended method for reading and writing Roman Telescope data files in Python.

.. note::

   If you have trouble installing this package, have encountered a bug while running the pipeline, or wish to request a new feature,
   please `open an issue on GitHub <https://github.com/spacetelescope/romancal/issues>`_ or `contact the Roman Telescope Help Desk <https://romanhelp.stsci.edu>`_.

Detailed explanations of specific calibration stages, reference files, and pipeline builds can be found on `RDox <https://roman-docs.stsci.edu/data-handbook-home/roman-data-pipelines>`_.

============
Contributing
============

``romancal`` is an open source package written in Python.
The source code is `available on GitHub <https://github.com/spacetelescope/romancal>`_.
New contributions and contributors are very welcome!

Please read `CONTRIBUTING.md <https://github.com/spacetelescope/romancal/blob/main/CONTRIBUTING.md>`_.

We strive to provide a welcoming community by abiding with our `CODE_OF_CONDUCT.md <https://github.com/spacetelescope/romancal/blob/main/CODE_OF_CONDUCT.md>`_.

--------------------------------

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
   roman/changes.rst
