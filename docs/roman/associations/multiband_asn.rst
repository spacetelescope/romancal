
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

Multiband catalog products always omit the optical element (filter) from their
output filenames, regardless of whether they contain one or multiple filters.
This naming convention distinguishes multiband products from single-band products.

Usage
^^^^^

Command Line Interface
""""""""""""""""""""""

To create a multiband association, use the following command:

.. code-block:: bash

    multiband_asn r00001_*full*_coadd.asdf

where the input files are in the current directory and the wildcard expands to
all relevant files. The tool will group files by their skycell identifier and
generate association files for each group. To get a complete list of options
you can run the command with the ``-h`` option:

.. code-block:: bash

    multiband_asn -h

Python API
""""""""""

You can also use the Python API:

.. code-block:: python

    from romancal.associations.multiband_asn import MultibandAssociation

    # Example with multiple filters for the same skycell
    files = [
        "r00001_p_full_270p65x48y69_f184_coadd.asdf",
        "r00001_p_full_270p65x48y69_f213_coadd.asdf"
    ]
    multiband = MultibandAssociation(files)

    # Or use a wildcard pattern
    # multiband = MultibandAssociation(["r00001_p_full_270p65x48y69_f*_coadd.asdf"])
    
    multiband.create_multiband_asn()
    # This will create: r00001_p_full_270p65x48y69_asn.json

The data release ID and product type are extracted from the filenames.
Association files are generated for each skycell group.

File Naming Conventions
^^^^^^^^^^^^^^^^^^^^^^^^

Input Filenames
"""""""""""""""

The input filenames should follow the convention:

.. code-block:: text

    rPPPPP_<data_release_id>_<product_type>_<skycell_id>_<filter>_coadd.asdf

Where:
    - PPPPP = Program number (5 digits, zero-padded)
    - data_release_id = Data release identifier (e.g., 'p' for prompt, 'dr01' for data release 1)
    - product_type = Product type (e.g., 'full', 'visit', 'pass')
    - skycell_id = Skycell identifier (e.g., '270p65x48y69')
    - filter = Filter identifier (e.g., 'f184', 'f213')

Output Filenames
""""""""""""""""

The association files will be JSON files named:

.. code-block:: text

    rPPPPP_<data_release_id>_<product_type>_<skycell_id>_asn.json

Note that the filter is **always omitted** from multiband association filenames,
regardless of the number of filters included in the association.

Examples
^^^^^^^^

Example 1: Prompt Multi-Filter Product
"""""""""""""""""""""""""""""""""""""""

**Input files:**

.. code-block:: text

    r00001_p_visit_x001y001_f184_coadd.asdf
    r00001_p_visit_x001y001_f213_coadd.asdf

**Command:**

.. code-block:: bash

    multiband_asn r00001_p_visit_x001y001_f*_coadd.asdf

**Output association:**

.. code-block:: text

    r00001_p_visit_x001y001_asn.json

Example 2: Prompt Single-Filter Product
""""""""""""""""""""""""""""""""""""""""

**Input files:**

.. code-block:: text

    r00001_p_full_x002y003_f184_coadd.asdf

**Command:**

.. code-block:: bash

    multiband_asn r00001_p_full_x002y003_f184_coadd.asdf

**Output association:**

.. code-block:: text

    r00001_p_full_x002y003_asn.json

Note: Even with a single filter, the filter identifier is omitted from the
multiband association filename to distinguish it from single-band products.

Example 3: Data Release Product
""""""""""""""""""""""""""""""""

**Input files:**

.. code-block:: text

    r00001_dr01_full_270p65x48y69_f184_coadd.asdf
    r00001_dr01_full_270p65x48y69_f213_coadd.asdf
    r00001_dr01_full_270p65x48y69_f158_coadd.asdf

**Command:**

.. code-block:: bash

    multiband_asn r00001_dr01_full_270p65x48y69_f*_coadd.asdf

**Output association:**

.. code-block:: text

    r00001_dr01_full_270p65x48y69_asn.json

Example 4: Multiple Skycells and Product Types
"""""""""""""""""""""""""""""""""""""""""""""""

**Input files:**

.. code-block:: text

    r00001_p_visit_x001y001_f184_coadd.asdf
    r00001_p_visit_x001y001_f213_coadd.asdf
    r00001_p_visit_x002y003_f184_coadd.asdf
    r00001_p_full_x001y001_f184_coadd.asdf
    r00001_p_full_x001y001_f213_coadd.asdf

**Command:**

.. code-block:: bash

    multiband_asn r00001_p_*_x*_f*_coadd.asdf

**Output associations:**

.. code-block:: text

    r00001_p_visit_x001y001_asn.json
    r00001_p_visit_x002y003_asn.json
    r00001_p_full_x001y001_asn.json

The tool automatically groups files by skycell and product type, creating
separate associations for each unique combination.
