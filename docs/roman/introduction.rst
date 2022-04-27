Introduction
============

This document provides instructions on running the Roman Science Calibration
Pipeline (referred to as "the pipeline") and individual pipeline steps.


Reference Files
===============

Many pipeline steps rely on the use of reference files that contain different types of
calibration data or information necessary for processing the data. The reference files are
instrument-specific and are periodically updated as the data processing evolves and the
understanding of the instruments improves. They are created, tested, and validated by the
Roman Instrument Teams. They ensure all the files are in the correct format and have all
required attributes. The files are then delivered to the Reference Data for Calibration
and Tools (ReDCaT) Management Team. The result of this process is the files being ingested
into the Roman Calibration Reference Data System (CRDS), and made available to the pipeline
team and any other ground subsystem that needs access to them.

Information about all the reference files used by the Calibration Pipeline can be found at
:ref:`reference_file_information`,
as well as in the documentation for each Calibration Step that uses a reference file.

CRDS
====

CRDS reference file mappings are usually set by default to always give access
to the most recent reference file deliveries and selection rules. On
occasion it might be necessary or desirable to use one of the non-default
mappings in order to, for example, run different versions of the pipeline
software or use older versions of the reference files. This can be
accomplished by setting the environment variable ``CRDS_CONTEXT`` to the
desired project mapping version, e.g.
::

$ export CRDS_CONTEXT='roman_0017.pmap'

Within STScI, the current storage location for all Roman CRDS reference files is:
::

/grp/crds/roman/references/roman/

Running From the Command Line
=============================

Individual steps and pipelines (consisting of a series of steps) can be run
from the command line using the ``strun`` command:
::

    $ strun <pipeline_name, class_name, or input_file>

The first argument to ``strun`` must be one of either a pipeline name, python
class of the step or pipeline to be run. The second argument to
``strun`` is the name of the input data file to be processed.

For example, the Stage 1 pipeline is implemented by the class
:ref:`romancal.pipeline.ExposurePipeline <exposure_pipeline>`. The command to
run this pipeline is:
::

  $ strun romancal.pipeline.ExposurePipeline r0008308002010007027_06311_0019_WFI01_uncal.asdf


Pipeline classes also have a **pipeline name**, or **alias**, that can be used
instead of thefull class specification. For example,
``calroman.pipeline.ExposurePipeline`` has the alias ``roman_elp`` and
can be run as
::

 $ strun roman_elp r0008308002010007027_06311_0019_WFI01_uncal.asdf

Exit Status
-----------
 ``strun`` produces the following exit status codes:

 - 0: Successful completion of the step/pipeline
 - 1: General error occurred
 - 64: No science data found

 .. _intro_file_conventions:

Input and Output File Conventions
=================================

.. _intro_input_file_discussion:

Input Files
-----------

There are two general types of input to any step or pipeline: references files
and data files.  The references files, unless explicitly
overridden, are provided through CRDS.

Data files are the science input, such as exposure ASDF files. All files are
assumed to be co-resident in the directory where the primary input file is
located.

.. _intro_output_file_discussion:

Output Files
------------

Output files will be created either in the current working directory, or where
specified by the :ref:`output_dir <intro_output_directory>` parameter.

File names for the outputs from pipelines and steps come from
three different sources:

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

See :ref:`pipeline_step_suffix_definitions` for a list of the more common
suffixes used.

Parameters
==========

All pipelines and steps have **parameters** that can be set to change various
aspects of how they execute. To see what parameters are available for any given
pipeline or step, use the ``-h`` option on ``strun``. Some examples are:
::

   $ strun roman_elp -h
   $ strun calroman.dq_init.DQInitStep -h

To set a parameter, simply specify it on the command line. For example, to have
:ref:`roman_elp <exposure_pipeline>` save the calibrated ramp files, the
``strun`` command would be as follows:
::

   $ strun roman_elp r0008308002010007027_06311_0019_WFI01_uncal.asdf --save_calibrated_ramp=true

To specify parameter values for an individual step when running a pipeline
use the syntax ``--steps.<step_name>.<parameter>=value``.
For example, to override the default selection of a dark current reference
file from CRDS when running a pipeline:
::

    $ strun roman_elp r0008308002010007027_06311_0019_WFI01_uncal.asdf
          --steps.dark_current.override_dark='my_dark.asdf'

Universal Parameters
--------------------

The set of parameters that are common to all pipelines and steps are referred to
as **universal parameters** and are described below.

.. _intro_output_directory:

Output Directory
^^^^^^^^^^^^^^^^

By default, all pipeline and step outputs will drop into the current
working directory, i.e., the directory in which the process is
running. To change this, use the ``output_dir`` parameter. For example, to
have all output from ``roman_elp``, including any saved
intermediate steps, appear in the sub-directory ``calibrated``, use
::

    $ strun roman_elp r0008308002010007027_06311_0019_WFI01_uncal.asdf
        --output_dir=calibrated

``output_dir`` can be specified at the step level, overriding what was
specified for the pipeline. From the example above, to change the name
and location of the ``dark_current`` step, use the following
::

    $ strun roman_elp r0008308002010007027_06311_0019_WFI01_uncal.asdf
        --output_dir=calibrated
        --steps.dark_current.output_file='dark_sub.asdf'
        --steps.dark_current.output_dir='dark_calibrated'

.. _intro_output_file:

Output File
^^^^^^^^^^^

When running a pipeline, the ``stpipe`` infrastructure automatically passes the
output data model from one step to the input of the next step, without
saving any intermediate results to disk.

.. _pipeline_step_suffix_definitions:

Pipeline/Step Suffix Definitions
--------------------------------

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
Uncalibrated raw input                         uncal
DQ initialization                              dq_init
Saturation detection                           saturation
Linearity correction                           linearity
Dark current                                   dark_current
Jump detection                                 jump
Corrected ramp data                            rampfit
Optional fitting results from ramp_fit step    fitopt
Assign WCS                                     assign_wcs
Flat field                                     flat
Photometric calibration			       phot
Calibrated image                               cal
=============================================  ============
