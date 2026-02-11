DARKDECAYSIGNAL reference file
------------------------------

:REFTYPE: DARKDECAYSIGNAL
:Data model: `~roman_datamodels.datamodels.DarkdecaysignalRefModel`

The DARKDECAYSIGNAL reference file contains amplitudes and time
constants for the exponential model for dark decay for each of
the 18 WFI detectors.  Time zero is referenced to the middle of
the first science read, so the typical amplitude of the signal
in the first read should correspond to the value of the amplitude
given in the reference file.  Decay times
are in seconds and amplitudes are in DN.

.. include:: ../references_general/darkdecaysignal_selection.inc

.. include:: ../includes/standard_keywords.inc

Reference File Format
+++++++++++++++++++++

DARKDECAYSIGNAL reference files are asdf files with a
decay_table dictionary.  This dictionary contains entries for
each WFI detector, WFI01--WFI18.  Each entry in the dictionary
in turn contains a dictionary with two keywords,
``amplitude`` and ``time_constant``, both floats.  These are
the amplitudes (in DN) and time constants (in seconds) of the
exponentially decaying dark decay model that needs to be
subtracted from uncalibrated Roman exposures.
