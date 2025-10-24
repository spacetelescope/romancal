.. _skycell_asn:

skycell_asn
===========

Create an association using either the command line tool
``skycell_asn`` or through the Python API using
:func:`romancal.associations.skycellasn.skycell_asn`

Associations
^^^^^^^^^^^^

Refer to TBD for a full description of associations.

To create an association, use the following command:

.. code-block:: python

		skycell_asn file_list -o my_program

where the individual calibrated detector files are in the current directory and the -o is the root
for the output associations. To get a complete list of options you can run the command with the
\-h option

.. code-block:: python

		skycell_asn -h


By knowing the structure of the conventional file name you can generate custom associations based
on the information included in the file names.
The current naming convention for the input files consists of a Visit_ID + Exposure_ID + detector + suffix
Where the Visit_ID is a 19 digit string
PPPPPCCAAASSSOOOVVV, with

PPPPP = Program number
CC = Execution plan number
AAA = Pass number (within execution plan)
SSS = Segment Number (within pass)
OOO = Observation number (within the leg, not to be confused with the Observation ID)
VVV = Visit Number (within the observation)

and the Exposure_ID is a four digits designating the exposures

eeee = Exposure number (within the visit)

The detector is wfi01, wfi02, ... wfi18

The suffix indicates the processing level of the file and to create products based on the
skycells the suffix should generally be 'cal'.

and the file name is constructed as,
rPPPPPCCAAASSSOOOVVV_eeee_detector_suffix.asdf

To give a more concrete example we'll use the APT example program for the
High Latitude Wide Angle Survey (hlwas). The program is 00991, the execution plan number is 01,
there are 12 passes, 169 segments, 5 observations and a various number of visits.
For this example we want to select a single filter, say F158, and that is observation 003.
So to generate the visit level associations for observation 001 we would select the files using the bash
command line,

.. code-block:: text

		skycell_asn r0099101001001003001_*_f158_cal.asdf  -o  r0099101001001003001 --product-type visit

Where the wildcard selects all the exposures for visit 001 and generates associations based on the skycells
the exposures will populate. This will generate associations based the skycells that the exposures can
contribute data to. The association files will be json files with names based

.. code-block:: text

	r00991_<data release ID>_v<visit_id>_<skycell name>_<filter>_asn.json

or for the selections above

.. code-block:: text

	r00991_p_v101001001003001_270p65x48y69_f158_asn.json

where the skycell name can vary based on the location on the celestial sphere.
The data release ID name can be changed from the default ('p')
by adding the optional argument --data-release-ID <new name> to the command line.

An analogous command to generate the pass level products, again setting observation to 003 to only select
the F158 filter.

.. code-block:: text

		skycell_asn r0099101001003001_*_f158_cal.asdf  -o  r00991 --product-type pass
