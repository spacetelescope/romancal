Running the Pipeline
====================

From the Command Line
-----------------------------

Individual steps and pipelines (consisting of a series of steps) can be run
from the command line using the ``strun`` command:
::

    $ strun <pipeline_name or class_name> input_file

The first argument to ``strun`` must be one of either a pipeline name, python
class of the step or pipeline to be run. The second argument to
``strun`` is the name of the input data file to be processed.
For a list of all the options available for ``strun``, please read the
`STPIPE Documentation <https://roman-pipeline.readthedocs.io/en/latest/roman/stpipe/index.html>`_.

For example, the exposure level  pipeline is implemented by the class
:ref:`romancal.pipeline.ExposurePipeline <exposure_pipeline>`. The command to
run this pipeline is:
::

  $ strun romancal.pipeline.ExposurePipeline r0008308002010007027_06311_0019_WFI01_uncal.asdf


Pipeline classes also have a **pipeline name**, or **alias**, that can be used
instead of the full class specification. For example,
``romancal.pipeline.ExposurePipeline`` has the alias ``roman_elp`` and
can be run as
::

 $ strun roman_elp r0008308002010007027_06311_0019_WFI01_uncal.asdf

The mosaic level pipeline can be run in a similar manner and is implemented using the class
:ref:`romancal.pipeline.MosaicPipeline <mosaic_pipeline>`.
The command to run this pipeline is:
::

  $ strun romancal.pipeline.MosaicPipeline r0008308002010007027_asn.json

An important point is that the mosaic level pipeline needs multiple exposures to run correctly. The
most convenient method to supply the input is to use an association. Instructions on how to create
an input association an be found at :ref:`asn-from-list`.

The mosaic level pipeline also has an alias, ``roman_mos``, and can be run as
::

 $ strun roman_mos r0008308002010007027_asn.json


Exit Status
```````````
 ``strun`` produces the following exit status codes:

 - 0: Successful completion of the step/pipeline
 - 1: General error occurred
 - 64: No science data found

 .. _intro_file_conventions:


From the Python Prompt
------------------------------

You can execute a pipeline or a step from within python by using the
``call`` method of the class.

The ``call`` method creates a new instance of the class and runs the pipeline or
step. Optional parameter settings can be specified by via keyword arguments or
supplying a parameter file. Some examples are shown below.

For the exposure pipeline and steps,

::

 from romancal.pipeline import ExposurePipeline
 result = ExposurePipeline.call('r0000101001001001001_01101_0001_WFI01_uncal.asdf')

 from romancal.linearity import LinearityStep
 result = LinearityStep.call('r0000101001001001001_01101_0001_WFI01_uncal.asdf')

One difference between the mosaic level pipeline and the exposure level pipeline is that the
mosaic level pipeline is generally designed to run on multiple overlapping exposures. To achieve
that the input to the pipeline is a list of images, usually an association.
For the mosaic level pipeline and steps,

::

 from romancal.pipeline import MosaicPipeline
 result = ExposurePipeline.call('r0000101001001001001_asn.json')

 from romancal.skymatch import SkyMatchStep
 result = SkyMatchStep.call('r0000101001001001001_asn.json')


For more information, see :ref:`Execute via call()<call_examples>`

For details on the different ways to run a pipeline step, see
the :ref:`Configuring a Step<configuring-a-step>` page.
