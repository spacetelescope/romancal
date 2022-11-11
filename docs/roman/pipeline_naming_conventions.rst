I/O File Conventions
=================================

.. _intro_input_file_discussion:

Input Files
-----------

There are two general types of input to any step or pipeline: references files
and data files.  The references files, unless explicitly
overridden, are provided through CRDS.

Data files are the science input, such as exposure `ASDF files <https://asdf-standard.readthedocs.io/>`_. All files are
assumed to be co-resident in the directory where the primary input file is
located.

.. _intro_output_file_discussion:

Output Files
------------

Output files will be created either in the current working directory, or where
specified by the :ref:`output_dir <intro_output_directory>` parameter.

File names for the outputs from pipelines and steps come from
two different sources:

- The name of the input file
- As specified by the :ref:`output_file <intro_output_file>` parameter

Regardless of the source, each pipeline/step uses the name as a base
name, onto which several different suffixes are appended, which
indicate the type of data in that particular file. A list of the main suffixes
can be :ref:`found below <pipeline_step_suffix_definitions>`.

The pipelines do not manage versions. When re-running a pipeline, previous files
will be overwritten.

Individual Step Outputs
^^^^^^^^^^^^^^^^^^^^^^^

If individual steps are executed without an output file name specified via
the ``output_file`` parameter, the ``stpipe`` infrastructure
automatically uses the input file name as the root of the output file name
and appends the name of the step as an additional suffix to the input file
name. If the input file name already has a known suffix, that suffix
will be replaced. For example:
::

   $ strun romancal.dq_init.DQInitStep r0008308002010007027_06311_0019_WFI01_uncal.asdf

produces an output file named
``r0008308002010007027_06311_0019_WFI01_dq_init.asdf``.


.. _pipeline_step_suffix_definitions:

Suffix Definitions
^^^^^^^^^^^^^^^^^^

However the output file name is determined (:ref:`see above
<intro_output_file_discussion>`), the various pipeline modules
will use that file name, along with a set of predetermined suffixes, to compose
output file names. The output file name suffix will always replace any known
suffix of the input file name. Each pipeline module uses the appropriate suffix
for the product(s) it is creating. The list of suffixes is shown in the
following table. Replacement occurs only if the suffix is one known to the
calibration code. Otherwise, the new suffix will simply be appended to the
basename of the file.

=============================================  ============
Product                                        Suffix
=============================================  ============
Uncalibrated raw input                         `_uncal`
DQ initialization                              `_dq_init`
Saturation detection                           `_saturation`
Linearity correction                           `_linearity`
Dark current                                   `_dark_current`
Jump detection                                 `_jump`
Corrected ramp data                            `_rampfit`
Optional fitting results from ramp_fit step    `_fitopt`
Assign WCS                                     `_assignwcs`
Flat field                                     `_flat`
Photometric calibration			               `_phot`
Calibrated image                               `_cal`
=============================================  ============
