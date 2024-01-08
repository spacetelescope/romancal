Parameters
==========

All pipelines and steps have **parameters** that can be set to change various
aspects of how they execute. To see what parameters are available for any given
pipeline or step, use the ``-h`` option on ``strun``. Some examples are:
::

   $ strun roman_elp -h
   $ strun romancal.dq_init.DQInitStep -h

   $strun roman_hlp -h
   $strun romancal.skymatch.SkyMatchStep -h

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
