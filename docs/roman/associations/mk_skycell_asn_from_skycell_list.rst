.. _mk_skycell_asn_from_skycell_list:

mk_skycell_asn_from_skycell_list
=================================

The function reads a list of ascii match files that were created using the ``mk_skycell_list``
program and generates a list of associations for all the skycells that touch the given exposure.
using either the command line tool
``mk_skycell_asn_from_skycell_list`` or through the Python API using
:func:`romancal.associations.mk_skycell_asn_from_skycell_list.mk_skycell_asn_from_skycell_list`

Association Files
^^^^^^^^^^^^^^^^^^^

To create the list of association files, use the following command:

.. code-block:: python

		mk_skycell_asn_from_skycell_list file_list

where the individual calibrated detector files are in the current directory and the file_list is a
list of the match files that will be used to generate the list of associations.
To get a complete list of options you can run the command with the
\-h option

.. code-block:: python

		mk_skycell_asn_from_skycell_list -h
                usage: mk_skycell_asn_from_skycell_list *.match

                Create level 3 associations from a list of match files

                positional arguments:
                filelist              A list of match files to generate level 3 asn's

                options:
                -h, --help            show this help message and exit
                --release-product RELEASE_PRODUCT
                                    The release product when creating the association
                --optical_element OPTICAL_ELEMENT
                                    The optical element used for the visit

The match file names are based on the input calibrated file names with the .asdf extention replaced
with .match and the match file contains the name of the calibrated exposure file and a list of index
values of the skycells that touch the input exposure.

.. code-block:: text

		mk_skycell_asn_from_skycell_list r0099101001001003001_*_cal.match

Where the wildcard selects all the exposures for visit 001 and generates associations based on the
list of skycell matches for each of the exposures. There will be one association for each unique
skycell index in the match files.
