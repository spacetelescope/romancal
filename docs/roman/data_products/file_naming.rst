.. _file_naming_schemes:

File Naming Schemes
-------------------

.. _exp_file_names:

Exposure file names
^^^^^^^^^^^^^^^^^^^
The names of the exposure level data are constructed with information from the
science header of the exposure, allowing users to map it to the observation in their corresponding
APT files. The ASDF file naming scheme for the Exposure products is::

 r<ppppp><cc><aaa><sss><ooo><vvv>_<gg><s><aa>_<eeeee>_<detector>_<prodType>.fits

where

 - ppppp: program ID number
 - cc: Execution Plan number
 - aaa: Pass Number (within execution plan)
 - sss: Segment Number (within pass)
 - ooo: observation number
 - vvv: visit number
 - gg: visit group
 - s: sequence ID (1=prime, >1 parallel)
 - aa: activity number (base 36)
 - eeeee: exposure number
 - detector: detector name (e.g. WFI01, WFI02, ...)
 - prodType: product type identifier (e.g. 'uncal', 'cal')

The standard <prodType> for the pipeline are uncal and cal, for the input products and resulting
calibrated product. There are optional suffixes for intermediate products that are not routinely
produced by the pipeline and are based of the processing level and can include dqinit, saturation,
linearity, jump, darkcurrent, rampfit, assignwcs, flat, and photom.

An example Exposure Level product ASDF file name is::

 r0000101001001001001_01101_0001_WFI01_cal.asdf
