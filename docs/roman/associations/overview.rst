.. _asn-overview:

====================
Association Overview
====================

.. _asn-what-are-associations:

What are Associations?
======================

Associations are basically just lists of things, mostly exposures,
that are somehow related. With respect to Roman Mission and the Data Management
System (DMS), associations have the following characteristics:

  * Relationships between multiple exposures are captured in an association.

  * An association is a means of identifying a set of exposures that belong
    together and may be dependent upon one another.

  * The association concept permits exposures to be calibrated, archived,
    retrieved, and reprocessed as a set rather than as individual objects.

  * For each association, DMS will generate the most combined and least combined
    data products.

.. _asn-associations-and-roman:

Associations and Roman
======================

The basic chunk in which science data arrives from the observatory is
termed an *exposure*. An exposure contains the data from detector reads that
for the Roman mission are set by the MA table (Multiaccum Table). These
resultants are the product transmitted to the ground and a set of these
constitutes an exposure for the detector. In general, it takes many
exposures to make up a single observation, and a whole program is made
up of a large number of observations.

On first arrival, an exposure is termed to be at *Level 0*: The only
transformation that has occurred is the extraction of the science data
from the observatory telemetry into a ASDF file. At this point, the
science exposures enter the calibration pipeline.

The pipeline consists of the ELP (Exposure Level Processing) and
the HLP (High Level Processing) which together comprise three levels of data generation and processing:
Level 1, Level 2, and Level 3. Level 1 data consist of uncalibrated individual
exposures, containing raw pixel information, formatted into the shape of
the detectors. Level 2 data have been  processed to correct for instrument artifacts and
have appropriate astrometric and geometric distortion information attached,
and with the exception of grism data, are in units that have known scaling
with flux. The resulting files contain flux
and spatially calibrated data, called *Level 2* data. The information
contained in these files are  still related to  individual exposures.

In order to combine or compare exposures, the data are resampled to a
regularized grid, removing the geometric distortion of the original pixels.
This process creates  *Level 3* data.

Utilities
=========

There are a number of utilities to create user-specific associations that are
documented under :ref:`Association Commands<association-commands>`.
