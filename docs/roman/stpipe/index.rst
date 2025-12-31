======
Stpipe
======

:doc:`Stpipe <stpipe:index>` is used for constructing and running
data pipelines and steps.

.. _stpipe-user-steps:

For Users
=========

See :ref:`running-the-pipeline` for an overview of using the romancal
pipeline.

TODO remove user_step (stpipe describes the same thing)
TODO remove user_pipeline (stpipe describes the same thing)
TODO remove user_logging (stpipe describes the same thing)
TODO remove config_asdf (stpipe describes ALMOST the same thing) except Roman parameters (which might not be correct)
TODO remove call_via_call (stpipe describes the same thing)
TODO remove call_via_run (stpipe describes the same thing)
TODO remove parameter_files (stpipe describes the same thing)

Maybe completely remove docs/stpipe and instead expand `running-the-pipeline`

.. toctree::
   :maxdepth: 2

   user_step.rst
   user_pipeline.rst
   user_logging.rst
   config_asdf.rst
   call_via_call.rst
   call_via_run.rst
   parameter_files.rst


.. _stpipe-devel-steps:

For Developers
==============

Developers should consult the :doc:`stpipe documentation <stpipe:index>`
for details about constructing new steps, pipelines and using logging.

.. automodapi:: romancal.stpipe
