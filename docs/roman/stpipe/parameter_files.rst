.. _parameter_files:

Parameter Files
===============

Parameter files can be used to specify parameter values when running a
pipeline or individual steps. These values can be
overridden either by the command line options and/or a
local parameter file. See :ref:`Parameter Precedence` for a full description of
how a parameter gets its final value.

A parameter file should be used when there are parameters a user wishes to
change from the default version for a custom run of the step. To create a
parameter file add ``--save-parameters <filename.asdf>`` to the command:
::

$ strun <step.class> <required-input-files> --save-parameters <filename.asdf>

For example, to save the parameters used for a run of the ``ExposurePipeline``
pipeline, use:
::

$ strun roman_elp r0000101001001001001_01101_0001_WFI01_uncal.asdf --save-parameters my_exppars.asdf

Once saved, the file can be edited, removing parameters that should be left
at their default values, and setting the remaining parameters to the
desired values. Once modified, the new parameter file can be used:
::

$ strun my_exppars2.asdf r0000101001001001001_01101_0001_WFI01_uncal.asdf

Note that the parameter values will reflect whatever was set on the
command-line, or through a specified local parameter file. In short, the
values will be those actually used in the running of the step.

For more information about and editing of parameter files, see
:ref:`config_asdf_files`.


More information on parameter files can be found in the ``stpipe`` User's
Guide at :ref:`stpipe-user-steps`.
