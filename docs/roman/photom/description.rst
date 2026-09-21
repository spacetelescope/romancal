Description
============

Algorithm
---------
The ``photom`` step adds flux (photometric) calibrations to the metadata of
a data product. The calibration information is read from a photometric reference
file, and the exact nature of the calibration information loaded from the reference file
is described below. This step does not affect the pixel values of the data product.

Upon successful completion of the photometric correction, the "photom" keyword in
"cal_step" in the metadata is set to "COMPLETE".


Photom and Pixel Area Data
--------------------------------

The photom reference file contains a table of
exposure parameters that define the flux
conversion and pixel area data for each optical element. The table contains one row
for each optical_element, and the photom step searches the
table for the row that matches the parameters of the science exposure and
then loads the calibration information from that row of the table.



For these table-based PHOTOM reference files, the calibration information in each
row includes a scalar flux conversion constant, the conversion uncertainty, and the nominal pixel area.

The scalar conversion constant is copied to
``meta.photometry.conversion_megajanskys``, which gives the conversion
from DN/s to megaJy/steradian, and the uncertainty is copied to
``meta.photometry.conversion_megajanskys_uncertainty``.

The step also populates ``meta.photometry.pixel_area`` in the science
data product, which gives the average pixel area over the detector in
steradians. This is a single nominal value and does not describe how the
pixel area varies across the detector; code that needs the solid angle
of a particular pixel should compute it from the WCS instead.
