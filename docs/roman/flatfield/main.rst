Description
===========
At its basic level this step flat-fields an input science data set by dividing
by a flat-field reference image. In particular, the SCI array from the
flat-field reference file is divided into both the SCI and ERR arrays of the
science data set, and the flat-field DQ array is combined with the science DQ
array using a bitwise OR operation.

Upon completion of the step, the cal_step attribute "flat_field" gets set
to "COMPLETE" in the output science data.

Imaging Data
------------
Imaging data use a straight-forward approach that involves applying a single
flat-field reference file to the science image. The processing steps are:

- Find pixels that have a value of NaN or zero in the FLAT reference file
  SCI array and set their DQ values to "NO_FLAT_FIELD."

- Reset the values of pixels in the flat that have DQ="NO_FLAT_FIELD" to
  1.0, so that they have no effect when applied to the science data.

- Apply the flat by dividing it into the science exposure SCI and ERR arrays.

- Propagate the FLAT reference file DQ values into the science exposure
  DQ array using a bitwise OR operation.

Error Propagation
-----------------
The VAR_POISSON and VAR_RNOISE variance arrays of the science exposure
are divided by the square of the flat-field value for each pixel.
A flat-field variance array, VAR_FLAT, is created from the science exposure
and flat-field reference file data using the following formula:

 The flat-field is applied to the science data, in-place, according to:

.. math::
   SCI^{\prime}_{science} = SCI_{science} / SCI_{flat}

The flat-field data is also applied to the VAR_POISSON and VAR_RNOISE
variance arrays,

.. math::
   VAR\_POISSON_{science} = VAR\_POISSON_{science} / SCI_{flat}^2

.. math::
   VAR\_RNOISE_{science} = VAR\_RNOISE_{science} / SCI_{flat}^2

The variance for the flat-fielded data associated with the science
data is determined using,

.. math::
   VAR\_FLAT_{science} = ( (SCI^{\prime}_{science})^{2} / SCI_{flat}^{2} ) * ERR_{flat}^{2}

and finally the error that is associated with the science data is given by,

.. math::
   ERR_{science} = \sqrt{VAR\_POISSON + VAR\_RNOISE + VAR\_FLAT}

The total ERR array in the science exposure is updated as the square root
of the quadratic sum of VAR_POISSON, VAR_RNOISE, and VAR_FLAT.
