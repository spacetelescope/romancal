.. _mk_skycell_list:

mk_skycell_list
===============

Creates a list of files that contain the calibrated file name with the
skycell index for all the skycells that touch the given exposure in
ascii format. These are created using either the command line tool
``mk_skycell_list`` or through the Python API using
:func:`romancal.associations.mk_skycell_list.mk_skycell_list`


Match Files
^^^^^^^^^^^

To create the list of match files, use the following command:

.. code-block:: python

		mk_skycell_list file_list

where the individual calibrated detector files are in a directory.
To get a complete list of options you can run the command with the
\-h option

.. code-block:: python

		mk_skycell_list -h
                usage: mk_skycell_list *_cal.asdf

                Create list of skycells from a level 2 file or files

                positional arguments:
                filelist    Input file list to generate a list of skycells

                options:
                 -h, --help            show this help message and exit
                 --output_dir OUTPUT_DIR
                           The optional directory to write the list of skycells

The match file names are based on the input Level 2 calibrated file names with the .asdf extention replaced
with .match and the match file contains the name of the calibrated exposure file and a list of index
values of the skycells that touch the input exposure. The output directory for the match files is
determined from the input directory of the calibrated exposures. You can override this by using
the flag to set the output directory.

As an example to create a list of match files for all the exposures in a visit the command is

.. code-block:: text

		mk_skycell_list r0099101001001003001_*_cal.asdf

Where the wildcard selects all the exposures for visit 001 and generates a list of skycell matches for each
of the exposures. There will be one match file for each input exposure file. For the example above this
will create match files named, with the `*` replaced by the exposure number of the file,

.. code-block:: text

		r0099101001001003001_*_cal.match
