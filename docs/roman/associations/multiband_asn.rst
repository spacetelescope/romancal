
.. _multiband_asn:

multiband_asn
=============

Create multiband associations using either the command line tool
``multiband_asn`` or through the Python API using
:class:`romancal.associations.multiband_asn.MultibandAssociation`.

Multiband Associations
^^^^^^^^^^^^^^^^^^^^^^

This module groups input files by their skycell identifier, creating an
association file for each unique skycell. Each association file contains all
filters (observations) that were used for that specific skycell, allowing for
multiband data products to be generated per skycell. This enables efficient
organization and processing of data across multiple filters for the same region
of the sky.

To create a multiband association, use the following command:

.. code-block:: bash

    multiband_asn r00001_*full*_coadd.asdf

where the input files are in the current directory and the wildcard expands to
all relevant files. The tool will group files by their skycell identifier and
generate association files for each group. To get a complete list of options
you can run the command with the ``-h`` option:

.. code-block:: bash

    multiband_asn -h

The input filenames should follow the convention:

.. code-block:: text

    rPPPPPCCAAASSSOOOVVV_<data_release_id>_<product_type>_<skycell_id>_coadd.asdf

Where:
    PPPPP = Program number
    CC = Execution plan number
    AAA = Pass number (within execution plan)
    SSS = Segment Number (within pass)
    OOO = Observation number (within the leg)
    VVV = Visit Number (within the observation)
    data_release_id = Data release identifier (e.g., 'p' for prompt)
    product_type = Product type (e.g., 'full')
    skycell_id = Skycell identifier (e.g., '270p65x48y69')

The association files will be JSON files named:

.. code-block:: text

    rPPPPPCCAAASSSOOOVVV_<data_release_id>_<product_type>_<skycell_id>_asn.json

For example, to generate associations for all skycells in a program:

.. code-block:: bash

    multiband_asn r00001_*full*_coadd.asdf

This will create association files for each unique skycell, grouping the
corresponding input files.

You can also use the Python API:

.. code-block:: python

    from romancal.associations.multiband_asn import MultibandAssociation
    files = ["r00001_p_full_270p65x48y69_coadd.asdf", ...]
    multiband = MultibandAssociation(files)
    multiband.create_multiband_asn()

The data release ID and product type are extracted from the filenames.
Association files are generated for each skycell group.
