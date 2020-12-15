Introduction
============

This document provides instructions on running the Roman Science Calibration
Pipeline (referred to as "the pipeline") and individual pipeline steps.


Reference Files
===============

Many pipeline steps rely on the use of reference files that contain different types of
calibration data or information necessary for processing the data. The reference files are
instrument-specific and are periodically updated as the data processing evolves and the
understanding of the instruments improves. They are created, tested, and validated by the
Roman Instrument Teams. They ensure all the files are in the correct format and have all
required attributes. The files are then delivered to the Reference Data for Calibration
and Tools (ReDCaT) Management Team. The result of this process is the files being ingested
into the Roman Calibration Reference Data System (CRDS), and made available to the pipeline
team and any other ground subsystem that needs access to them.

Information about all the reference files used by the Calibration Pipeline can be found at
:ref:`reference_file_information`,
as well as in the documentation for each Calibration Step that uses a reference file.

CRDS
====

CRDS reference file mappings are usually set by default to always give access
to the most recent reference file deliveries and selection rules. On
occasion it might be necessary or desirable to use one of the non-default
mappings in order to, for example, run different versions of the pipeline
software or use older versions of the reference files. This can be
accomplished by setting the environment variable ``CRDS_CONTEXT`` to the
desired project mapping version, e.g.
::

$ export CRDS_CONTEXT='roman_0421.pmap'

Within STScI, the current storage location for all Roman CRDS reference files is:
::

/grp/crds/roman/references/roman/
