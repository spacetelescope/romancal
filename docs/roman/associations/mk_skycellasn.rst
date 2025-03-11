.. _mk_skycellasn:

mk_skycellasn
=============

THe function reads a list ascii match files that that were created using the mk_patchlist
program and generates a list of associations for all the skycells that touch the given exposure.
using either the command line tool
``mk_skycellasn`` or through the Python API using
:func:`romancal.associations.mk_skycellasn.mk_skycellasn`

Running this command requires that you have the patch to the
file containing the division of the sky into patches. This is done
by setting the environment variable PATCH_TABLE_PATH.
As an example:
::
   
   export PATCH_TABLE_PATH=<location of my patch table>

.. Note::

   The patch table will be available in the CRDS system soon.

   
.. Note::

   **For STScI Users Only:**
    Users at STScI may access the required
    data files from the Central Storage network. Set the following
    environment variables in your ``bash`` shell. (You will probably
    want to add this to your bash setup file.) ::
      
      export PATCH_TABLE_PATH="/grp/roman/scsb/tesselation/patches.asdf"


Association Files
^^^^^^^^^^^^^^^^^^^

To create the list of association files, use the following command:

.. code-block:: python

		mk_skycellasn file_list 

where the individual calibrated detector files are in the current directory and the file_list is a
list of the match files that will be used to generate the list of associations. 
To get a complete list of options you can run the command with the
\-h option

.. code-block:: python

		mk_skycellasn -h
                usage: mk_skycellasn *.match

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

		mk_skycellasn r0099101001001003001_*_cal.match

Where the wildcard selects all the exposures for visit 001 and generates associations based on the
list of skycell matches for each of the exposures. There will be one association for each unique
skycell index in the match files. 

