.. _run_examples:

Executing a pipeline or pipeline step directly, or via run()
============================================================

When calling a pipeline or step instance directly, or using the ``run`` method,
you can specify individual parameter values manually. In this case, parameter
files are not used. If you use ``run`` after instantiating with a parameter
file, the parameter file will be ignored.

::

 # Instantiate the class. Do not provide a parameter file.
 pipe = ExposurePipeline()

 # Manually set any desired non-default parameter values
 pipe.assign_wcs.skip = True
 pipe.jump.rejection_threshold = 5
 pipe.ramp_fit.override_gain = 'my_gain_file.asdf'
 pipe.save_result = True
 pipe.output_dir = '/my/data/pipeline_outputs'

 # Run the pipeline
 result = pipe.run('r0000101001001001001_01101_0001_WFI01_uncal.asdf')

To run a single step:

::

 from romancal.jump import JumpStep

 # Instantiate the step
 step = JumpStep()

 # Set parameter values
 step.rejection_threshold = 5
 step.save_results = True
 step.output_dir = '/my/data/jump_data'

 # Execute using the run method
 result = step.run('r0000101001001001001_01101_0001_WFI01_linearity.asdf')
