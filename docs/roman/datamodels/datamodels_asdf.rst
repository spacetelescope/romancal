.. _datamodels_asdf:


Working with Roman datamodels and ASDF files
============================================

If you've installed the roman calibration pipeline you should also have access
to the standalone tool asdfinfo which allows access to asdf (and roman) files
from the terminal prompt,::

    asdftool info r0000101001001001001_01101_0001_WFI16_cal.asdf
    root (AsdfObject)
    ├─asdf_library (Software)
    │ ├─author (str): The ASDF Developers
    ...

also useful is::

    asdftool help
    usage: asdftool [-h] [--verbose] {help,explode,implode,extract,defragment,diff,edit,remove-hdu,info,extensions,tags,to_yaml} ...

    Commandline utilities for managing ASDF files.

    optional arguments:
      -h, --help            show this help message and exit
      --verbose, -v         Increase verbosity

    subcommands:

      {help,explode,implode,extract,defragment,diff,edit,remove-hdu,info,extensions,tags,to_yaml}
        help                Display usage information
        explode             Explode a ASDF file.
        implode             Implode a ASDF file.
        extract             Extract ASDF extensions in ASDF-in-FITS files into pure ASDF files
        defragment          Defragment an ASDF file..
        diff                Report differences between two ASDF files
        remove-hdu          Remove ASDF extension from ASDF-in-FITS file
        info                Print a rendering of an ASDF tree.
        extensions          Show information about installed extensions
        tags                List currently available tags
        to_yaml             Convert as ASDF file to pure YAML.


which gives a list of possible actions one of the more useful can be::

    asdftool edit file.asdf

Which will open the file in an editor you have set via the EDITOR environment variable.
A more complete description of the `asdftool` can be found at _.

.. _a link: https://asdf.readthedocs.io/en/stable/asdf/asdf_tool.html

To access the files via a python session,

.. code-block:: python

    import roman_datamodels as rdm
    import asdf
    rdm_a = rdm.open('r0000101001001001001_01101_0001_WFI16_cal.asdf')
    asdf_a = asdf.open('r0000101001001001001_01101_0001_WFI16_cal.asdf', mode='rw')

Once the files are loaded you can access various attributes. Below is a table
showing how to access various properties using the roman_datamodels and the
asdf.open methods,

+--------------------------------------+---------------------------------------------------------------+
| Roman Datamodels                     | ASDF                                                          |
+--------------------------------------+---------------------------------------------------------------+
| .. code-block:: python               | .. code-block:: python                                        |
|                                      |                                                               |
|   rdm_a.meta                         |    asdf_a.tree['roman']['meta']                               |
|   rdm_a.meta.aperture                |    asdf_a.tree['roman']['meta']['aperture']                   |
|   rdm_a.meta.aperture.position_angle |    asdf_a.tree['roman']['meta']['aperture']['position_angle'] |
|   120                                |    120                                                        |
+--------------------------------------+---------------------------------------------------------------+

You can also update or modify the metadata in Roman datamodels

.. code-block:: python

    rdm_a.meta.aperture.position_angle = 120.21
    rdm_a.meta.aperture.position_angle
    120.21

The ASDF equivalent is

.. code-block:: python

    asdf_a.tree['roman']['meta']['aperture']['position_angle'] = 120.05
    asdf_a.tree['roman']['meta']['aperture']['position_angle']
    120.05

.. HINT::
    If you trigger an error,
    "ValueError: assignment destination is read-only"
    make sure the asdf file was opened with mode='rw'

You can also access and modify the data arrays

.. code-block:: python
    :caption: Roman Datamodels

    rdm_a.data
    <array (unloaded) shape: [4096, 4096] dtype: float32>

    rdm_a.data[10,11]
    0.0

    rdm_a.data[10,11] = 122.1
    rdm_a.data[10,11]
    122.1

.. code-block:: python
    :caption: ASDF

    asdf_a.tree['roman']['data']
    <array (unloaded) shape: [4096, 4096] dtype: float32>

    asdf_a.tree['roman']['data'][10,11]
    0.0

    asdf_a.tree['roman']['data'][10,11] = 3.14159
    asdf_a.tree['roman']['data'][10,11]
    3.14159
