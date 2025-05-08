PSF Reference File
------------------

:REFTYPE: PSF
:Data model: `~roman_datamodels.datamodels.EpsfRefModel`

The PSF reference file contains a grid of point spread
images for a given detector position, filter, and 
spectral type.

.. include:: ../references_general/psf_selection.inc

.. include:: ../includes/standard_keywords.inc

Type Specific Keywords for the PSF
+++++++++++++++++++++++++++++++++++++
In addition to the standard reference file keywords listed above,
the following keywords are *required* in PSF reference files,
because they are used as CRDS selectors.
(see :ref:`psf_selectors`):

===============  ===================================== ==========
Keyword          Data Model Name                       Instrument
===============  ===================================== ==========
detector         model.meta.instrument.detector        WFI
optical_element  model.meta.instrument.optical_element WFI
===============  ===================================== ==========

Reference File Format
+++++++++++++++++++++
PSF reference files are ASDF format, with four data objects.
The format and content of the file is as follows:

===================  ============ ======================  =============
Data                 Object Type  Dimensions              Data type
===================  ============ ======================  =============
psf                  NDArray      3 x 3 x 9 x 391 x 391   float32
extended_psf         NDArray      3641 x 3641             float32  
extended_psf_noipc   NDArray      3641 x 3641             float32
psf_noipc            NDArray      3 x 3 x 9 x 391 x 391   float32
===================  ============ ======================  =============

The values in the ``psf`` and ``psf_noipc`` arrays contain a grid of
images with three focus positions, three spectral types and nine
positions on the detector. The focus is the number of defocus waves in
the reference files (0, 1, or 2) and 0 corresponds to in-focus position. The
spectral type is in order of [A0V, G2V, and M5V], and grid position is
an index related to the meta.pixel_x and meta.pixel_y values in the
file metadata. The field of view of the ePSF is 361 oversampled pixels
(0.11 native arcseconds / pixel oversampled by a factor of 4 yielding
approximately 10 arcseconds on a side), while the extended PSF is 3641
pixels (or approximately 100 arcseconds on a side). The detector
positions where the PSFs are simulated are: x = [4.0, 2047.5, 4091.0,
4.0, 2047.5, 4091.0, 4.0, 2047.5, 4091.0], y = [4.0, 4.0, 4.0, 2047.5,
2047.5, 2047.5, 4091.0, 4091.0, 4091.0]. The psf and extended_psf
arrays have been convolved with an interpixel capacitance (IPC) kernel
to make them more like empirical PSFs (ePSFs) that will be observed
on-orbit.

The ASDF file contains four data arrays.
