Description
===========

Algorithm
---------

Processing multiple datasets together allows for the identification of bad pixels
or cosmic-rays that remain in each of the input images, many times at levels which
were not detectable by the :ref:`jump <jump_step>` step. The ``outlier_detection`` step
implements the following algorithm to identify and flag any remaining cosmic-rays or
other artifacts left over from previous calibrations:

  - TBD

The outlier detection step serves as a single interface to apply this general
process to any Roman data, with specific variations of this algorithm for each
type of data.  Sub-classes of the outlier detection algorithm have been developed
specifically for

  - TBD
