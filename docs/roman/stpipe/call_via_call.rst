.. _call_examples:

Executing a pipeline or pipeline step via call()
================================================

The ``call`` method will create an instance and run a pipeline or pipeline step
in a single call.

::

 from romancal.pipeline import ExposurePipeline
 result = ExposurePipeline.call('r0000101001001001001_01101_0001_WFI01_uncal.asdf')

 from romancal.linearity import LinearityStep
 result = LinearityStep.call('r0000101001001001001_01101_0001_WFI01_dqinit.asdf')


To set custom parameter values when using the ``call`` method, set the
parameters in the pipeline or parameter file and then supply the file using the
``config_file`` keyword: ::

 # Calling a pipeline
 result = ExposurePipeline.call('r0000101001001001001_01101_0001_WFI01_uncal.asdf', config_file='exp_pars.asdf'))

 # Calling a step
 result = LinearityStep.call('r0000101001001001001_01101_0001_WFI01_dqinit.asdf', config_file = 'linearity_pars.asdf')


When running a pipeline, parameter values can also be supplied in the call to ``call`` itself by using a nested dictionary of step and
parameter names:

::

 result = ExposurePipeline.call('r0000101001001001001_01101_0001_WFI01_uncal.asdf', config_file='exp_pars.asdf', steps={"jump":{"rejection_threshold": 200}})

When running a single step with ``call``, parameter values can be supplied more simply:

::

 result = JumpStep.call("r0000101001001001001_01101_0001_WFI01_dqinit.asdf", rejection_threshold=200)

Running steps and pipelines with ``call`` also allows for the specification of a logging
configuration file using the keyword ``logcfg``:

::

 result = ExposurePipeline.call("r0000101001001001001_01101_0001_WFI01_dqinit.asdf",
                                 config_file="exp_pars.asdf",
                                 logcfg="my-logging-config.cfg")


Where are the results?
----------------------

A fundamental difference between running steps and pipelines in Python as
opposed to from the command line using ``strun`` is whether files are created or
not. When using ``strun``, results are automatically saved to files because that
is the only way to access the results.

However, when running within a Python interpreter or script, the presumption is
that results will be used within the script. As such, results are not
automatically saved to files. It is left to the user to decide when to save.

If one wishes for results to be saved by a particular ``call``, use the
parameter ``save_results=True``::

 result = JumpStep.call("r0000101001001001001_01101_0001_WFI01_dqinit.asdf",
                        rejection_threshold=200, save_results=True)

If one wishes to specify a different file name, rather than a system-generated
one, set :ref:`output_file<intro_output_file>` and/or
:ref:`output_dir<intro_output_directory>`.
