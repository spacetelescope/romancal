.. _mk_patchfile:

mk_patchlist
============

Creates a list of files that contain the calibrated file name with the
skycell index for all the skycells that touch the given exposure in
ascii format. These are created using either the command line tool
``mk_patchfile`` or through the Python API using
:func:`romancal.associations.mk_patchlist.mk_patchlist`

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

      
Match Files
^^^^^^^^^^^

To create the list of match files, use the following command:

.. code-block:: python

		mk_patchlist file_list 

where the individual calibrated detector files are in the current directory and the -o is the root
for the output associations. To get a complete list of options you can run the command with the
\-h option

.. code-block:: python

		mk_patchlist -h
                usage: mk_patchlist *_cal.asdf

                Create list of patches from a level 2 file or files

                positional arguments:
                filelist    Input file list to generate a list of patchs

                options:
                -h, --help  show this help message and exit

The match file names are based on the input calibrated file names with the .asdf extention replaced
with .match and the match file contains the name of the calibrated exposure file and a list of index
values of the skycells that touch the input exposure. 

As an example to create a list of match files for all the exposures in a visit the command is

.. code-block:: text

		mk_patchlist r0099101001001003001_*_cal.asdf

Where the wildcard selects all the exposures for visit 001 and generates a list of skycell matches for each 
of the exposures. There will be one match file for each input exposure file. For the example above this
will create match files named, with the `*` replaced by the exposure number of the file, 

.. code-block:: text

		r0099101001001003001_*_cal.match
