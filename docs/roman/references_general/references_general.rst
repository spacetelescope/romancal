Introduction
============

This document is intended to be a core reference guide to the formats, naming convention and
data quality flags used by the reference files for pipeline steps requiring them, and is not
intended to be a detailed description of each of those pipeline steps. It also does not give
details on pipeline steps that do not use reference files.
The present manual is the living document for the reference file specifications.

Reference File Naming Convention
================================

Before reference files are ingested into CRDS, they are renamed following a
convention used by the pipeline. As with any other changes undergone by the reference files,
the previous names are kept in header keywords, so the Instrument Teams
can easily track which delivered file is being used by the pipeline in each step.

The naming of reference files uses the following syntax::

 roman_<instrument>_<reftype>_<version>.<extension>

where

- ``instrument`` is currently "WFI"
- ``reftype`` is one of the type names listed in the table below
- ``version`` is a 4-digit version number (e.g. 0042)
- ``extension`` gives the file format, "asdf"

An example WFI FLAT reference file name would be "roman_wfi_flat_0042.asdf".


Reference File Types
====================

Most reference files have a one-to-one relationship with calibration steps, e.g.
there is one step that uses one type of reference file. Some steps, however, use
several types of reference files and some reference file types are used by more
than one step. The tables below show the correspondence between pipeline steps and
refernece file types. The first table is ordered by pipeline step, while the second
is ordered by reference file type. Links to the reference file types provide detailed
documentation on each reference file.

+---------------------------------------------+--------------------------------------------------+
| Pipeline Step                               | Reference File Type (REFTYPE)                    |
+=============================================+==================================================+
+---------------------------------------------+--------------------------------------------------+
| :ref:`flatfield <flatfield_step>`           | :ref:`FLAT <flat_reffile>`                       |
+---------------------------------------------+--------------------------------------------------+


+--------------------------------------------------+---------------------------------------------+
| Reference File Type (REFTYPE)                    | Pipeline Step                               |
+==================================================+=============================================+
| :ref:`FLAT <flat_reffile>`                       | :ref:`flatfield <flatfield_step>`           |
+--------------------------------------------------+---------------------------------------------+

.. _`Standard ASDF metadata`:

Standard ASDF metadata
======================

Al Roman science and reference files are ASDF files.

The required Keywords Documenting Contents of Reference Files are:

=========== ==================================================================================
Keyword     Comment
=========== ==================================================================================
reftype     `FLAT    Required values are listed in the discussion of each pipeline step.`
description `Summary of file content and/or reason for delivery`
author      `Fred Jones     Person(s) who created the file`
useafter    `YYYY-MM-DDThh:mm:ss Date and time after the reference files will
            be used. The T is required. Time string may NOT be omitted;
            use T00:00:00 if no meaningful value is available.
            Astropy Time objects are allowed.`
pedigree    `Options are
            'SIMULATION'
            'GROUND'
            'DUMMY'
            'INFLIGHT YYYY-MM-DD YYYY-MM-DD'`
history     `Description of Reference File Creation`
telescope   `ROMAN   Name of the telescope/project.`
instrument  `WFI   Instrument name.`
=========== ==================================================================================
