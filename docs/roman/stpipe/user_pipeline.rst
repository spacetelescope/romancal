.. _stpipe-user-pipelines:

=========
Pipelines
=========

It is important to note that a Pipeline is also a Step, so everything
that applies to a Step in the :ref:`stpipe-user-steps` chapter also
applies to Pipelines.

Configuring a Pipeline
======================

This section describes how to set parameters on the individual steps
in a pipeline.  To change the order of steps in a pipeline, one must
write a Pipeline subclass in Python.  That is described in the
:ref:`devel-pipelines` section of the developer documentation.

Just as with Steps, Pipelines can by configured either by a
parameter file or directly from Python.

From a parameter file
---------------------

A Pipeline parameter file follows the same format as a Step parameter file:
:ref:`config_asdf_files`

Here is an example pipeline parameter file for the `ExposurePipeline`
class:

.. code-block:: yaml

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
    - class: romancal.dq_init.dq_init_step.DQInitStep
    name: dq_init
    parameters:
     input_dir: ''
     output_dir: null
     output_ext: .asdf
     output_file: null
     output_use_index: true
     output_use_model: false
     post_hooks: []
     pre_hooks: []
    save_results: false
    search_output_file: true
    skip: false
    suffix: null
    - class: romancal.saturation.saturation_step.SaturationStep
    ...


Just like a ``Step``, it must have ``name`` and ``class`` values.
Here the ``class`` must refer to a subclass of `stpipe.Pipeline`.

Following ``name`` and ``class`` is the ``steps`` section.  Under
this section is a subsection for each step in the pipeline.  The easiest
way to get started on a parameter file is to call ``Step.export_config`` and
then edit the file that is created.  This will generate an ASDF config file
that includes every available parameter, which can then be trimmed to the
parameters that require customization.

For each Step’s section, the parameters for that step may either be
specified inline, or specified by referencing an external
parameter file just for that step.  For example, a pipeline
parameter file that contains:

.. code-block:: yaml

    - class: romancal.jump.jump_step.JumpStep
      name: jump
      parameters:
        flag_4_neighbors: true
        four_group_rejection_threshold: 190.0
        input_dir: ''
        max_jump_to_flag_neighbors: 1000.0
        maximum_cores: none
        min_jump_to_flag_neighbors: 10.0

is equivalent to:

.. code-block:: yaml

   steps:
   - class: romancal.jump.jump_step.JumpStep
     name: jump
     parameters:
        config_file = myjump.asdf

with the file ``myjump.asdf.`` in the same directory:

.. code-block:: yaml

   class: romancal.jump.jump_step.JumpStep
   name: jump
   parameters:
     flag_4_neighbors: true
     four_group_rejection_threshold: 190.0

If both a ``config_file`` and additional parameters are specified, the
``config_file`` is loaded, and then the local parameters override
them.

Any optional parameters for each Step may be omitted, in which case
defaults will be used.


From Python
-----------

A pipeline may be configured from Python by passing a nested
dictionary of parameters to the Pipeline’s constructor.  Each key is
the name of a step, and the value is another dictionary containing
parameters for that step.  For example, the following is the
equivalent of the parameter file above:

.. code-block:: python

    from stpipe.pipeline import Image2Pipeline

    steps = {
        'jump':{'rejection_threshold': 180.,
                'three_group_rejection_threshold': 190.,
                'four_group_rejection_threshold':195.
    }

    pipe = ExposurePipeline(steps=steps)

Running a Pipeline
==================

From the commandline
--------------------

The same ``strun`` script used to run Steps from the commandline can
also run Pipelines.

The only wrinkle is that any parameters overridden from the
commandline use dot notation to specify the parameter name.  For
example, to override the ``rejection_threshold`` value on the ``jump``
step in the example above, one can do::

    > strun romancal.pipeline.ExposurePipeline --steps.jump.rejection_threshold=180.

From Python
-----------

Once the pipeline has been configured (as above), just call the
instance to run it.

    pipe(input_data)

Caching details
---------------

The results of a Step are cached using Python pickles.  This allows
virtually most of the standard Python data types to be cached.  In
addition, any ASDF models that are the result of a step are saved as
standalone ASDF files to make them more easily used by external tools.
The filenames are based on the name of the substep within the
pipeline.

Hooks
=====

Each Step in a pipeline can also have pre- and post-hooks associated.
Hooks themselves are Step instances, but there are some conveniences
provided to make them easier to specify in a parameter file.

Pre-hooks are run right before the Step.  The inputs to the pre-hook
are the same as the inputs to their parent Step.
Post-hooks are run right after the Step.  The inputs to the post-hook
are the return value(s) from the parent Step. The return values are
always passed as a list. If the return value from the parent Step is a
single item, a list of this single item is passed to the post hooks.
This allows the post hooks to modify the return results, if necessary.

Hooks are specified using the ``pre_hooks`` and ``post_hooks`` parameters
associated with each step. More than one pre- or post-hook may be assigned, and
they are run in the order they are given. There can also be ``pre_hooks`` and
``post_hooks`` on the Pipeline as a whole (since a Pipeline is also a Step).
Each of these parameters is a list of strings, where each entry is one of:

   - An external commandline application.  The arguments can be
     accessed using {0}, {1} etc.  (See
     `stpipe.subproc.SystemCall`).

   - A dot-separated path to a Python Step class.

   - A dot-separated path to a Python function.
