.. _config_asdf_files:

ASDF Parameter Files
====================

ASDF is the format of choice for parameter files. `ASDF
<https://asdf-standard.readthedocs.io/>`_ stands for "Advanced Scientific Data
Format", a general purpose, non-proprietary, and system-agnostic format for the
dissemination of data. Built on `YAML <https://yaml.org/>`_, the most basic file
is text-based requiring minimal formatting.

.. _asdf_minimal_file:

To create a parameter file, the most direct way is to choose the Pipeline
class, Step class, or already existing .asdf or .cfg file, and run that step
using the ``--save-parameters`` option. For example, to get the parameters for
the ``ExposurePipeline`` pipeline, do the following: ::

   $ strun --save-parameters=exp_pars.asdf roman_elp r0000101001001001001_01101_0001_WFI01_uncal.asdf

Once created and modified as necessary, the file can now be used by ``strun``
to run the step/pipeline with the desired parameters:
::

   $ strun exp_pars.asdf r0000101001001001001_01101_0001_WFI01_uncal.asdf

The remaining sections will describe the file format and contents.

File Contents
-------------

To describe the contents of an ASDF file, the configuration for the step
``roman_elp`` will be used as the example:

.. code-block::

    #ASDF 1.0.0
    #ASDF_STANDARD 1.5.0
    %YAML 1.1
    %TAG ! tag:stsci.edu:asdf/
    --- !core/asdf-1.1.0
    asdf_library: !core/software-1.0.0 {author: The ASDF Developers, homepage: 'http://github.com/asdf-format/asdf',
      name: asdf, version: 2.13.0}
    history:
      extensions:
      - !core/extension_metadata-1.0.0
        extension_class: asdf.extension.BuiltinExtension
        software: !core/software-1.0.0 {name: asdf, version: 2.13.0}
    class: romancal.pipeline.exposure_pipeline.ExposurePipeline
    meta:
      author: <SPECIFY>
      date: '2022-09-15T13:59:54'
      description: Parameters for calibration step romancal.pipeline.exposure_pipeline.ExposurePipeline
      instrument: {name: <SPECIFY>}
      origin: <SPECIFY>
      pedigree: <SPECIFY>
      reftype: <SPECIFY>
      telescope: <SPECIFY>
      useafter: <SPECIFY>
    name: ExposurePipeline
    parameters:
      input_dir: ''
      output_dir: null
      output_ext: .asdf
      output_file: null
      output_use_index: true
      output_use_model: false
      post_hooks: []
      pre_hooks: []
      save_calibrated_ramp: false
      save_results: true
      search_output_file: true
      skip: false
      suffix: null
    steps:
    - class: romancal.jump.jump_step.JumpStep
      name: jump
      parameters:
        flag_4_neighbors: true
        four_group_rejection_threshold: 190.0
        input_dir: ''
        max_jump_to_flag_neighbors: 1000.0
        maximum_cores: none
        min_jump_to_flag_neighbors: 10.0
        output_dir: null
        output_ext: .asdf
        output_file: null
        output_use_index: true
        output_use_model: false
        post_hooks: []
        pre_hooks: []
        rejection_threshold: 180.0
        save_results: false
        search_output_file: true
        skip: false
        suffix: null
        three_group_rejection_threshold: 185.0
    ...

Required Components
~~~~~~~~~~~~~~~~~~~

Preamble
++++++++

The first 5 lines, up to and including the "---" line, define the file as an
ASDF file. The rest of the file is formatted as one would format YAML data.
Being YAML, the last line, containing the three ``...`` is essential.

class and name
++++++++++++++

There are two required keys at the top level: ``class`` and ``parameters``.
``parameters`` is discussed below.

``class`` specifies the Python class to run.  It should be a
fully-qualified Python path to the class.  Step classes can ship with
``stpipe`` itself, they may be part of other Python packages, or they
exist in freestanding modules alongside the configuration file.  For
example, to use the ``SystemCall`` step included with ``stpipe``, set
``class`` to ``stpipe.subprocess.SystemCall``.  To use a class called
``Custom`` defined in a file ``mysteps.py`` in the same directory as
the configuration file, set ``class`` to ``mysteps.Custom``.

``name`` defines the name of the step.  This is distinct from the
class of the step, since the same class of Step may be configured in
different ways, and it is useful to be able to have a way of
distinguishing between them.  For example, when Steps are combined
into :ref:`stpipe-user-pipelines`, a Pipeline may use the same Step class
multiple times, each with different configuration parameters.

Parameters
++++++++++

``parameters`` contains all the parameters to pass onto the step. The order of
the parameters does not matter. It is not necessary to specify all parameters
either. If not defined, the default, as defined in the code or values from CRDS
parameter references, will be used.

Formatting
**********

YAML has two ways of formatting a list of key/value pairs. In the above example,
each key/value pair is on separate line. The other way is using a form that is similar to a Python ``dict``.
For example, the ``parameters`` block above could also have been formatted as:

.. code-block::

    parameters: {flag_4_neighbors: true, four_group_rejection_threshold: 190.0,
      max_jump_to_flag_neighbors: 1000.0, maximum_cores: none,
      min_jump_to_flag_neighbors: 10.0, output_dir: null, output_ext: .asdf,
      output_file: null, output_use_index: true, output_use_model: false,
      rejection_threshold: 180.0, three_group_rejection_threshold: 185.0}

Optional Components
~~~~~~~~~~~~~~~~~~~

The ``asdf_library`` and ``history`` blocks are necessary only when a parameter
file is to be used as a parameter reference file in CRDS which is not currently
implemented in the Roman pipeline.

.. _`Completeness`:

Completeness
~~~~~~~~~~~~

For any parameter file, it is not necessary to specify all step/pipeline
parameters. Any parameter left unspecified will get, at least, the default value
define in the step's code. If a parameter is defined without a default value,
and the parameter is never assigned a value, an error will be produced when the
step is executed.

Remember that parameter values can come from numerous sources. Refer to
:ref:`Parameter Precedence` for a full listing of how parameters can be set.

From the ``JumpStep`` example, if all that needed to change is the
``rejection_threshold`` parameter with a setting of ``80.0``,
the ``parameters`` block need only contain the following:

.. code-block::

    parameters:
      rejection_threshold: 80.0


Pipeline Configuration
~~~~~~~~~~~~~~~~~~~~~~

Pipelines are essentially steps that refer to sub-steps. As in the original cfg
format, parameters for sub-steps can also be specified. All sub-step parameters
appear in a key called `steps`. Sub-step parameters are specified by using the
sub-step name as the key, then underneath and indented, the parameters to change
for that sub-step. For example, to define the ``rejection_threshold`` of the
``JumpStep`` step in a ``ExposurePipeline`` parameter file, the parameter
block would look as follows:

.. code-block::

   class: romancal.pipeline.exposure_pipeline.ExposurePipeline
   parameters: {}
   steps:
   - class: romancal.jump.jump_step.JumpStep
     parameters:
       rejection_threshold: 80.0

As with step parameter files, not all sub-steps need to be specified. If left
unspecified, the sub-steps will be run with their default parameter sets. For
the example above, the other steps of ``ExposurePipeline``, such as ``assign_wcs``
and ``photom`` would still be executed.

Similarly, to skip a particular step, one would specify ``skip: true`` for that
substep. Continuing from the above example, to skip the ``flatfield`` step,
the parameter file would look like:

.. code-block::

   class: romancal.pipeline.exposure_pipeline.ExposurePipeline
   parameters: {}
   steps:
   - class: romancal.flatfield.flat_field_step.FlatFieldStep
     name: flatfield
     parameters:
       skip: true

.. note::

    In the previous examples, one may have noted the line parameters: {}. Often
    when configuring a pipeline, one needs not set any parameters for the pipeline
    itself. However, the keyword ``parameters`` is required. As such,
    the value for ``parameters`` is defined as an empty dictionary, ``{}``.

Python API
----------

There are a number of ways to create an ASDF parameter file. From the
command line utility ``strun``, the option ``--save-parameters`` can be used.

Within a Python script, the method ``Step.export_config(filename: str)`` can be
used. For example, to create a parameter file for ``JumpStep``, use the
following:

.. doctest-skip::

   >>> from romancal.jump import JumpStep
   >>> step = JumpStep()
   >>> step.export_config('jump_step.asdf')



History
~~~~~~~

Parameter reference files also require at least one history entry. This can be found in the ``history`` block under ``entries``:

.. code-block::

    history:
      extensions:
      - !core/extension_metadata-1.0.0
        extension_class: asdf.extension.BuiltinExtension
        software: !core/software-1.0.0 {name: asdf, version: 2.13.0}
    history:
      entries:
      - !core/history_entry-1.0.0 {description: Base values, time: !!timestamp '2019-10-29
          21:20:50'}

It is highly suggested to use the ASDF API to add history entries:

.. doctest-skip::

   >>> import asdf
   >>> cfg = asdf.open('config.asdf')
       #
       # Modify `parameters` and `meta` as necessary.
       #
   >>> cfg.add_history_entry('Parameters modified for some reason')
   >>> cfg.write_to('config_modified.asdf')

Roman, Parameters and Parameter References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In general, the default parameters for any pipeline or step are valid for nearly
all instruments and observing modes. This means that when a pipeline or step is
run without any explicit parameter setting, that pipeline or step will usually
do the desired operation. Hence, most of the time there is no need for a
parameter reference to be provided by the user. Only for a
small set of observing mode/step combinations, will there be need to create a
parameter reference. Even then, nearly all cases will involve changing a subset
of a pipeline or step parameters.

Keeping this sparse-population philosophy in mind, for most parameter
references, only those parameters that are explicitly changed should be
specified in the reference. If adhered to, when a pipeline/step default value
for a particular parameter needs to change, the change will be immediately
available. Otherwise, all references that mistakenly set said parameter will
need to be updated. See :ref:`Completeness` for more information.

Furthermore, every pipeline/step have a common set of parameters, listed
below. These parameters generally affect the infrastructure operation of
pipelines/steps, and should not be included in a parameter reference.

- input_dir
- output_ext
- output_use_index
- output_use_model
- post_hooks
- pre_hooks
- save_results
- search_output_file
