The dq_def extension contains the bit assignments used in the DQ array.
It contains the following 4 columns:

===========  ===========  =================================================
Flags         Data Type    Description
===========  ===========  =================================================
BIT          int32         The bit number, starting at zero
VALUE        uint32        The equivalent base-10 value of BIT
NAME         string        The mnemonic name of the data quality condition
DESCRIPTION  string        A description of the data quality condition
===========  ===========  =================================================

**NOTE:** For more information on standard bit definitions see: :ref:`Data Quality Flags`.


