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
the previous names are kept in header, so the Instrument Teams
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
reference file types. The first table is ordered by pipeline step, while the second
is ordered by reference file type. Links to the reference file types provide detailed
documentation on each reference file.

+---------------------------------------------+--------------------------------------------------+
| Pipeline Step                               | Reference File Type (reftype)                    |
+=============================================+==================================================+
| :ref:`assign_wcs <assign_wcs_step>`         | :ref:`DISTORTION <distortion_reffile>`           |
+---------------------------------------------+--------------------------------------------------+
| :ref:`dark_current <dark_current_step>`     | :ref:`DARK <dark_reffile>`                       |
+---------------------------------------------+--------------------------------------------------+
| :ref:`dq_init <dq_init_step>`               | :ref:`MASK <mask_reffile>`                       |
+---------------------------------------------+--------------------------------------------------+
| :ref:`flatfield <flatfield_step>`           | :ref:`FLAT <flat_reffile>`                       |
+---------------------------------------------+--------------------------------------------------+
| :ref:`jump_detection <jump_step>`           | :ref:`GAIN <gain_reffile>`                       |
+                                             +--------------------------------------------------+
|                                             | :ref:`READNOISE <readnoise_reffile>`             |
+---------------------------------------------+--------------------------------------------------+
| :ref:`linearity <linearity_step>`           | :ref:`LINEARITY <linearity_reffile>`             |
+---------------------------------------------+--------------------------------------------------+
| :ref:`photom <photom_step>`                 | :ref:`PHOTOM <photom_reffile>`                   |
+---------------------------------------------+--------------------------------------------------+
| :ref:`ramp_fitting <ramp_fitting_step>`     | :ref:`GAIN <gain_reffile>`                       |
+                                             +--------------------------------------------------+
|                                             | :ref:`READNOISE <readnoise_reffile>`             |
+---------------------------------------------+--------------------------------------------------+
| :ref:`saturation <saturation_step>`         | :ref:`SATURATION <saturation_reffile>`           |
+---------------------------------------------+--------------------------------------------------+


+--------------------------------------------------+---------------------------------------------+
| Reference File Type (reftype)                    | Pipeline Step                               |
+==================================================+=============================================+
| :ref:`DARK <dark_reffile>`                       | :ref:`dark_current <dark_current_step>`     |
+--------------------------------------------------+---------------------------------------------+
| :ref:`DISTORTION <distortion_reffile>`           | :ref:`assign_wcs <assign_wcs_step>`         |
+--------------------------------------------------+---------------------------------------------+
| :ref:`FLAT <flat_reffile>`                       | :ref:`flatfield <flatfield_step>`           |
+--------------------------------------------------+---------------------------------------------+
| :ref:`GAIN <gain_reffile>`                       | :ref:`jump_detection <jump_step>`           |
+                                                  +---------------------------------------------+
|                                                  | :ref:`ramp_fitting <ramp_fitting_step>`     |
+--------------------------------------------------+---------------------------------------------+
| :ref:`LINEARITY <linearity_reffile>`             | :ref:`linearity <linearity_step>`           |
+--------------------------------------------------+---------------------------------------------+
| :ref:`MASK <mask_reffile>`                       | :ref:`dq_init <dq_init_step>`               |
+--------------------------------------------------+---------------------------------------------+
| :ref:`PHOTOM <photom_reffile>`                   | :ref:`photom <photom_step>`                 |
+--------------------------------------------------+---------------------------------------------+
| :ref:`READNOISE <readnoise_reffile>`             | :ref:`jump_detection <jump_step>`           |
+                                                  +---------------------------------------------+
|                                                  | :ref:`ramp_fitting <ramp_fitting_step>`     |
+--------------------------------------------------+---------------------------------------------+
| :ref:`SATURATION <saturation_reffile>`           | :ref:`saturation <saturation_step>`         |
+--------------------------------------------------+---------------------------------------------+

.. _`Standard ASDF metadata`:

Standard ASDF metadata
======================

Al Roman science and reference files are ASDF files.

The required attributes Documenting Contents of Reference Files are:

=========== ==================================================================================
Attribute     Comment
=========== ==================================================================================
reftype     `FLAT    Required values are listed in the discussion of each pipeline step.`
description `Summary of file content and/or reason for delivery.`
author      `Fred Jones     Person(s) who created the file.`
useafter    `YYYY-MM-DDThh:mm:ss Date and time after the reference files will
            be used. The T is required. Time string may NOT be omitted;
            use T00:00:00 if no meaningful value is available.
            Astropy Time objects are allowed.`
pedigree    `Options are
            'SIMULATION'
            'GROUND'
            'DUMMY'
            'INFLIGHT YYYY-MM-DD YYYY-MM-DD'`
history     `Description of Reference File Creation`.
telescope   `ROMAN   Name of the telescope/project.`
instrument  `WFI   Instrument name.`
=========== ==================================================================================

Observing Mode Attributes
=========================

A pipeline module may require separate reference files for each instrument, detector,
optical element, observation date, etc.  The values of these parameters must be included in the
reference file attributes.  The observing-mode attributes are vital to the process of
ingesting reference files into CRDS, as they are used to establish the mapping between
observing modes and specific reference files. Some observing-mode attributes are also
used in the pipeline processing steps.

The Keywords Documenting the Observing Mode are:

===============  ==================  ==============================================================================
Keyword          Sample Value        Comment
===============  ==================  ==============================================================================
detector         WFI01               Allowed values WFI01, WFI02, ... WFI18

optical element  F158                Name of the filter element and includes PRISM and GRISM

exposure type    WFI_IMAGE           Allowed values WFI_IMAGE, WFI_GRATING, WFI_PRISM, WFI_DARK, WFI_FLAT, WFI_WFSC
===============  ==================  ==============================================================================

Tracking Pipeline Progress
++++++++++++++++++++++++++

As each pipeline step is applied to a sciece data product, it will record a status
indicator in a cal_step attribute of the science data product. These statuses
may be included in the primary header of reference files, in order to maintain
a history of the data that went into creating the reference file.
Allowed values for the status Attribute are  'INCOMPLETE', 'COMPLETE'
and 'SKIPPED'. The default value is set to 'INCOMPLETE'. The pipeline modules
will set the value to 'COMPLETE' or 'SKIPPED'. If the pipeline steps are run
manually and you skip a step the cal_step will remain 'INCOMPLETE'.

Data Quality Flags
==================

Within science data files, the PIXELDQ flags are stored as 32-bit integers;
the GROUPDQ flags are 8-bit integers. All calibrated data from a particular
instrument and observing mode have the same set of DQ flags in the same (bit)
order. The table below lists the allowed DQ flags. Only the first eight entries
in the table below are relevant to the GROUPDQ array.

Flags for the DQ, PIXELDQ, and GROUPDQ Arrays.

===  ==========    ================  ===========================================
Bit  Value         Name              Description
===  ==========    ================  ===========================================
0    1             DO_NOT_USE        Bad pixel. Do not use.
1    2             SATURATED         Pixel saturated during exposure
2    4             JUMP_DET          Jump detected during exposure
3    8             DROPOUT           Data lost in transmission
4    16            RESERVED_1
5    32            PERSISTENCE       High persistence (was RESERVED_2)
6    64            AD_FLOOR          Below A/D floor (0 DN, was RESERVED_3)
7    128           RESERVED_4
8    256           UNRELIABLE_ERROR  Uncertainty exceeds quoted error
9    512           NON_SCIENCE       Pixel not on science portion of detector
10   1024          DEAD              Dead pixel
11   2048          HOT               Hot pixel
12   4096          WARM              Warm pixel
13   8192          LOW_QE            Low quantum efficiency
15   32768         TELEGRAPH         Telegraph pixel
16   65536         NONLINEAR         Pixel highly nonlinear
17   131072        BAD_REF_PIXEL     Reference pixel cannot be used
18   262144        NO_FLAT_FIELD     Flat field cannot be measured
19   524288        NO_GAIN_VALUE     Gain cannot be measured
20   1048576       NO_LIN_CORR       Linearity correction not available
21   2097152       NO_SAT_CHECK      Saturation check not available
22   4194304       UNRELIABLE_BIAS   Bias variance large
23   8388608       UNRELIABLE_DARK   Dark variance large
24   16777216      UNRELIABLE_SLOPE  Slope variance large (i.e., noisy pixel)
25   33554432      UNRELIABLE_FLAT   Flat variance large
26   67108864      RESERVED_5
27   134217728     RESERVED_6
28   268435456     UNRELIABLE_RESET  Sensitive to reset anomaly
29   536870912     RESERVED_7
30   1073741824    OTHER_BAD_PIXEL   A catch-all flag
31   2147483648    REFERENCE_PIXEL   Pixel is a reference pixel
===  ==========    ================  ===========================================
