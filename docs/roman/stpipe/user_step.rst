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

.. _strun:

Running a Step from the commandline
-----------------------------------
The ``strun`` command can be used to run Steps from the commandline.

The first argument should be:

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
    3. CRDS-retrieved parameter reference
    4. ``Step``-coded default, determined by the parameter definition ``Step.spec``

For pipelines, if a pipeline parameter file specifies a value for a step in the
pipeline, that takes precedence over any step-specific value found, either from
a step-specific parameter file or CRDS-retrieved step-specific parameter file.
The full order of precedence for a pipeline and its sub steps is as follows:

    1. Value specified on the command-line: ``strun pipeline.asdf --steps.step.par=value_that_will_be_used``
    2. Value found in the user-specified pipeline parameter file: ``strun pipeline.asdf``
    3. Value found in the parameter file specified in a pipeline parameter file
    4. CRDS-retrieved parameter reference for the pipeline
    5. CRDS-retrieved parameter reference for each sub-step
    6. ``Pipeline``-coded default for itself and all sub-steps
    7. ``Step``-coded default for each sub-step


Debugging
`````````

To output all logging output from the step, add the `--verbose` option
to the commandline.  (If more fine-grained control over logging is
required, see :ref:`user-logging`).

To start the Python debugger if the step itself raises an exception,
pass the `--debug` option to the commandline.


CRDS Retrieval of Step Parameters
`````````````````````````````````

In general, CRDS uses the input to a ``Step`` to determine which reference files
to use. Nearly all Roman-related steps take only a single input file. However,
the input may be an association file. Since step parameters are
configured only once per execution of a step or pipeline, only the first
qualifying member, usually of type ``science`` is used.

Retrieval of ``Step`` parameters from CRDS can be completely disabled by
using the ``--disable-crds-steppars`` command-line switch, or setting the
environment variable ``STPIPE_DISABLE_CRDS_STEPPARS`` to ``true``.

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

One can implement parameter reference file retrieval and use of a local
parameter file as follows::

    from stpipe import config_parser
    from romancal.flatfield import FlatFieldStep

    config = FlatFieldStep.get_config_from_reference(input)
    local_config = config_parser.load_config_file('my_flatfield_config.asdf')
    config_parser.merge_config(config, local_config)

    flat_field_step = FlatFieldStep.from_config_section(config)
    output = flat_field_step.run(input)

Using the ``.run()`` method is the same as calling the instance directly.
They are equivalent::

    output = mystep(input)
