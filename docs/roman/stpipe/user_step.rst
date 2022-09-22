=====
Steps
=====

.. _configuring-a-step:

Configuring a Step
==================

This section describes how to instantiate a Step and set configuration
parameters on it.

Steps can be configured by:

    - Instantiating the Step directly from Python

    - Reading the input from a parameter file

.. _running_a_step_from_a_configuration_file:

Running a Step from a parameter file
====================================

A parameter file contains one or more of a ``Step``'s parameters. Any parameter
not specified in the file will take its value from the defaults coded
directly into the ``Step``. Note that any parameter specified on the command
line overrides all other values.

The format of parameter files is the :ref:`config_asdf_files` format.
Refer to the :ref:`minimal example<asdf_minimal_file>` for a complete
description of the contents. The rest of this document will focus on the step
parameters themselves.

Every parameter file must contain the key ``class``, followed by
the optional ``name`` followed by any parameters that are specific to the step
being run.

``class`` specifies the Python class to run.  It should be a
fully-qualified Python path to the class.  Step classes can ship with
``stpipe`` itself, they may be part of other Python packages, or they
exist in freestanding modules alongside the configuration file.

``name`` defines the name of the step.  This is distinct from the
class of the step, since the same class of Step may be configured in
different ways, and it is useful to be able to have a way of
distinguishing between them.  For example, when Steps are combined
into :ref:`stpipe-user-pipelines`, a Pipeline may use the same Step class
multiple times, each with different configuration parameters.

The parameters specific to the Step all reside under the key ``parameters``. The
set of accepted parameters is defined in the Step’s spec member.  The easiest
way to get started on a parameter file is to call ``Step.export_config`` and
then edit the file that is created.  This will generate an ASDF config file
that includes every available parameter, which can then be trimmed to the
parameters that require customization.

Here is an example parameter file (``do_cleanup.asdf``) that runs the (imaginary)
step ``stpipe.cleanup`` to clean up an image.

.. code-block::

    #ASDF 1.0.0
    #ASDF_STANDARD 1.3.0
    %YAML 1.1
    %TAG ! tag:stsci.edu:asdf/
    --- !core/asdf-1.1.0
    class: stpipe.cleanup
    name: MyCleanup
    parameters:
      threshold: 42.0
      scale: 0.01
    ...

.. _strun:

Running a Step from the commandline
-----------------------------------
The ``strun`` command can be used to run Steps from the commandline.

The first argument may be either:

    - The a parameter file

   - A Python class

Additional parameters may be passed on the commandline. These parameters
override any defaults. Any extra positional
parameters on the commandline are passed to the step's process method. This will
often be input filenames.

To display a list of the parameters that are accepted for a given Step
class, pass the ``-h`` parameter, and the name of a Step class or
parameter file::

    $ strun -h romancal.dq_init.DQInitStep
    usage: strun [-h] [--logcfg LOGCFG] [--verbose] [--debug] [--save-parameters SAVE_PARAMETERS] [--disable-crds-steppars]
                 [--pre_hooks] [--post_hooks] [--output_file] [--output_dir] [--output_ext] [--output_use_model] [--output_use_index]
                 [--save_results] [--skip] [--suffix] [--search_output_file] [--input_dir] [--override_mask]
                 cfg_file_or_class [args ...]

    (selected) optional arguments:
      -h, --help       show this help message and exit
      --logcfg LOGCFG  The logging configuration file to load
      --verbose, -v    Turn on all logging messages
      --debug            When an exception occurs, invoke the Python debugger, pdb
      --save-parameters SAVE_PARAMETERS Save step parameters to specified file.
      --disable-crds-steppars Disable retrieval of step parameter references files from CRDS
      --output_file   File to save the output to

Every step has an `--output_file` parameter.  If one is not provided,
the output filename is determined based on the input file by appending
the name of the step.  For example, in this case, `foo.asdf` is output
to `foo_cleanup.asdf`.

Finally, the parameters a ``Step`` actually ran with can be saved to a new
parameter file using the `--save-parameters` option. This file will have all
the parameters, specific to the step, and the final values used.

.. _`Parameter Precedence`:

Parameter Precedence
````````````````````

There are a number of places where the value of a parameter can be specified.
The order of precedence, from most to least significant, for parameter value
assignment is as follows:

    1. Value specified on the command-line: ``strun step.asdf --par=value_that_will_be_used``
    2. Value found in the user-specified parameter file
    3. ``Step``-coded default, determined by the parameter definition ``Step.spec``

For pipelines, if a pipeline parameter file specifies a value for a step in the
pipeline, that takes precedence over any step-specific value found from
a step-specific parameter file.
The full order of precedence for a pipeline and its sub steps is as follows:

    1. Value specified on the command-line: ``strun pipeline.asdf --steps.step.par=value_that_will_be_used``
    2. Value found in the user-specified pipeline parameter file: ``strun pipeline.asdf``
    3. Value found in the parameter file specified in a pipeline parameter file
    4. ``Pipeline``-coded default for itself and all sub-steps
    5. ``Step``-coded default for each sub-step


Debugging
`````````

To output all logging output from the step, add the `--verbose` option
to the commandline.  (If more fine-grained control over logging is
required, see :ref:`user-logging`).

To start the Python debugger if the step itself raises an exception,
pass the `--debug` option to the commandline.


.. _run_step_from_python:

Running a Step in Python
------------------------

There are a number of methods to run a step within a Python interpreter,
depending on how much control one needs.

Step.from_cmdline()
```````````````````

For individuals who are used to using the ``strun`` command, `Step.from_cmdline`
is the most direct method of executing a step or pipeline. The only argument is
a list of strings, representing the command line arguments one would have used
for ``strun``. The call signature is::

    Step.from_cmdline([string,...])

For example, given the following command-line::

    $ strun romancal.pipeline.ExposurePipeline r0000101001001001001_01101_0001_WFI01_uncal.asdf \
            --steps.jump.override_gain=roman_wfi_gain_0033.asdf

the equivalent `from_cmdline` call would be::

    from romancal.pipeline import ExposurePipeline
    ExposurePipeline.from_cmdline([' r0000101001001001001_01101_0001_WFI01_uncal.asdf',
                                   'steps.jump.override_gain', 'roman_wfi_gain_0033.asdf'])


call()
``````

Class method `Step.call` is the slightly more programmatic, and preferred,
method of executing a step or pipeline. When using ``call``, one gets the full
configuration initialization that
one gets with the ``strun`` command or ``Step.from_cmdline`` method. The call
signature is::

    Step.call(input, logcfg=None, **parameters)

The positional argument ``input`` is the data to be operated on, usually a
string representing a file path or a :ref:`DataModel<datamodels>`. The optional
keyword argument ``config_file`` is used to specify a local parameter file. The
optional keyword argument ``logcfg`` is used to specify a logging configuration file.
Finally, the remaining optional keyword arguments are the parameters that the
particular step accepts. The method returns the result of the step. A basic
example is::

    from romancal.jump import JumpStep
    output = JumpStep.call('r0000101001001001001_01101_0001_WFI01_uncal.asdf')

makes a new instance of `JumpStep` and executes using the specified exposure
file. `JumpStep` has a parameter ``rejection_threshold``. To use a different
value than the default, the statement would be::

    output = JumpStep.call('r0000101001001001001_01101_0001_WFI01_uncal.asdf',
                           rejection_threshold=42.0)

If one wishes to use a :ref:`parameter file<parameter_files>`, specify the path
to it using the ``config_file`` argument::

    output = JumpStep.call('r0000101001001001001_01101_0001_WFI01_uncal.asdf',
                           config_file='my_jumpstep_config.asdf')

run()
`````

The instance method `Step.run()` is the lowest-level method to executing a step
or pipeline. Initialization and parameter settings are left up to the user. An
example is::

    from romancal.flatfield import FlatFieldStep

    mystep = FlatFieldStep()
    mystep.override_sflat = 'sflat.asdf'
    output = mystep.run(input)

`input` in this case can be a asdf file containing the appropriate data, or the output
of a previously run step/pipeline, which is an instance of a particular :ref:`datamodel<datamodels>`.

Unlike the ``call`` class method, there is no parameter initialization that
occurs, either by a local parameter file or from a CRDS-retrieved parameter
reference file. Parameters can be set individually on the instance, as is shown
above. Parameters can also be specified as keyword arguments when instantiating
the step. The previous example could be re-written as::

    from romancal.flatfield import FlatFieldStep

    mystep = FlatFieldStep(override_sflat='sflat.asdf')
    output = mystep.run(input)

Using the ``.run()`` method is the same as calling the instance directly.
They are equivalent::

    output = mystep(input)
