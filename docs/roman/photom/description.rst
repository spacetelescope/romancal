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

The scalar conversion constant is copied to the header keyword "conversion_megajanskys", which
gives the conversion from DN/s to megaJy/steradian, and converted to microJy/square arcseconds and saved to
the header keyword "conversion_microjanskys". The same process is performed for the uncertainty, with
the values saved in "conversion_megajanskys_uncertainty" and "conversion_microjanskys_uncertainty",
respectively.

The step also populates the metadata keywords "pixelarea_steradians" and "pixelarea_arcsecsq" in the
science data product, which give the average pixel area in units of
steradians and square arcseconds, respectively.
