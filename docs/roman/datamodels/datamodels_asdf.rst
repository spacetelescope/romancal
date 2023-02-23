.. _datamodels_asdf:


Working with Roman datamodels and ASDF files
============================================

Please refer to `Roman Documentation <https://roman-datamodels.readthedocs.io/en/latest/>`_
for more details about `roman_datamodels`.

This section assumes that you are familiar with the ASDF standard format.
If that's not the case, a good starting point would be to go over the `ASDF Documentation <https://asdf-standard.readthedocs.io/>`_.

If you have installed the roman calibration pipeline you should also have access
to the standalone tool asdfinfo which allows access to `ASDF <https://asdf-standard.readthedocs.io/>`_ (and roman) files
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
A more complete description of the options can be found in the
`asdftool <https://asdf.readthedocs.io/en/stable/asdf/asdf_tool.html>`_
documentation.

To access the files via a python session,

.. code-block:: python

    import roman_datamodels as rdm
    import asdf
    with rdm.open('r0000101001001001001_01101_0001_WFI16_cal.asdf') as model:
        <Manipulate the files>

    with asdf.open('r0000101001001001001_01101_0001_WFI16_cal.asdf', copy_arrays=True) as model:
        <Manipulate the files>

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
    make sure the asdf file was opened with copy_arrays=True, or
    with mode='rw'

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

or by modifying the ASDF tree,

.. code-block:: python
    :caption: ASDF

    asdf_a.tree['roman']['data']
    <array (unloaded) shape: [4096, 4096] dtype: float32>

    asdf_a.tree['roman']['data'][10,11]
    0.0

    asdf_a.tree['roman']['data'][10,11] = 3.14159
    asdf_a.tree['roman']['data'][10,11]
    3.14159

Using the info method
---------------------

You can examine a roman data model using the info and search methods provided
from the asdf package. The info fuction will print a representation of the
asdf tree.

.. code:: python

    >>> from roman_datamodels import datamodels as rdm  # doctest: +SKIP
    >>> d_uncal = rdm.open('r0000101001001001001_01101_0001_WFI01_uncal.asdf')  # doctest: +SKIP
    >>> d_uncal.info()  # doctest: +SKIP
    root (AsdfObject)
    ├─asdf_library (Software)
    │ ├─author (str): The ASDF Developers
    │ ├─homepage (str): http://github.com/asdf-format/asdf
    │ ├─name (str): asdf
    │ └─version (str): 2.8.1
    ├─history (dict)
    │ └─extensions (list)
    │   ├─[0] (ExtensionMetadata) ...
    │   ├─[1] (ExtensionMetadata) ...
    │   └─[2] (ExtensionMetadata) ...
    └─roman (WfiScienceRaw)
      ├─meta (dict)
      │ ├─aperture (Aperture) ...
      │ ├─cal_step (CalStep) ...
      │ ├─calibration_software_version (str): 0.4.3.dev89+gca5771d
      │ ├─coordinates (Coordinates) ...
      │ ├─crds_context_used (str): roman_0020.pmap
      │ ├─crds_software_version (str): 11.5.0
      │ ├─ephemeris (Ephemeris) ...
      │ ├─exposure (Exposure) ...
      │ └─17 not shown
      └─data (NDArrayType): shape=(8, 4096, 4096), dtype=uint16
    Some nodes not shown.

The info command also gives you control over the number of lines displayed
by passing the argument ``max_rows``. As an integer, ``max_rows``
will be interpreted as an overall limit on the number of displayed lines.
If ``max_rows`` is a tuple, then each member limits lines per node at the
depth corresponding to its tuple index.
For example, to show all top-level nodes and 5 of each's children:

.. code:: python

    >>> d_uncal.info(max_rows=(None,5))  # doctest: +SKIP
    root (AsdfObject)
    ├─asdf_library (Software)
    │ ├─author (str): The ASDF Developers
    │ ├─homepage (str): http://github.com/asdf-format/asdf
    │ ├─name (str): asdf
    │ └─version (str): 2.8.1
    ├─history (dict)
    │ └─extensions (list) ...
    └─roman (WfiScienceRaw)
      ├─meta (dict) ...
      └─data (NDArrayType): shape=(8, 4096, 4096), dtype=uint16
    Some nodes not shown.

Or you can use the asdf.info method to view the contents of the tree

.. code:: python

    import asdf
    asdf.info(d_uncal)

Will print the same information as the above `d_uncal.info` command but also
gives you enhanced capabilities. For instance you can display the first three
lines for each of the meta entries,

.. code:: python

    >>> asdf.info(d_uncal.meta,max_rows=(None, 3))  # doctest: +SKIP
    root (DNode)
    ├─aperture (Aperture)
    │ ├─name (str): WFI_CEN
    │ └─position_angle (int): 120
    ├─cal_step (CalStep)
    │ ├─assign_wcs (str): INCOMPLETE
    │ ├─flat_field (str): INCOMPLETE
    │ └─6 not shown
    ├─calibration_software_version (str): 0.4.3.dev89+gca5771d
    ├─coordinates (Coordinates)
    │ └─reference_frame (str): ICRS
    ├─crds_context_used (str): roman_0020.pmap
    ├─crds_software_version (str): 11.5.0
    ├─ephemeris (Ephemeris)
    │ ├─earth_angle (float): 3.3161255787892263
    │ ├─moon_angle (float): 3.3196162372932148
    │ └─10 not shown
    ...

or you can concentrate on a given attribute. To list all the attributes
in `cal_step` without listing the values,

.. code:: python

    >>> asdf.info(d_uncal.meta.cal_step,max_rows=(None, 3),show_values=False)  # doctest: +SKIP
    root (CalStep)
    ├─assign_wcs (str)
    ├─flat_field (str)
    ├─dark (str)
    ├─dq_init (str)
    ├─jump (str)
    ├─linearity (str)
    ├─ramp_fit (str)
    └─saturation (str)

More information on the info method can be found in the ASDF documentation at
`rendering the ASDF trees. <https://asdf.readthedocs.io/en/stable/asdf/features.html#endering-asdf-trees>`_

Using the search method
-----------------------

You can also use the search method to find attributes,

.. code:: python

    >>> d_uncal.search('cal_step')  # doctest: +SKIP
    root (AsdfObject)
    └─roman (WfiScienceRaw)
      └─meta (dict)
        └─cal_step (CalStep)

or a a general search for all attributes with cal in the name

.. code:: python

    >>> d_uncal.search('cal')  # doctest: +SKIP
    root (AsdfObject)
    └─roman (WfiScienceRaw)
     └─meta (dict)
       ├─cal_step (CalStep)
       ├─calibration_software_version (str): 0.4.3.dev89+gca5771d
       ├─instrument (WfiMode)
       │ └─optical_element (str): F158
       └─velocity_aberration (VelocityAberration)
         └─scale_factor (float): 0.9999723133902021

This will do a regular expression search for `cal` in the attribute name. More
information on using regular expressions in the search method can be found
in the ASDF documentation linked below.

To search only within the meta tree,

.. code:: python

    >>> d_uncal.search('cal_')['roman']['meta']  # doctest: +SKIP
    meta (dict)
    ├─cal_step (CalStep)
    └─instrument (WfiMode)
      └─optical_element (str): F158

You can also use the search method to find attributes by type in the asdf tree.
For instance, you can find all integers, floats, or booleans by using the type
keyword,

.. code:: python

    >>> d_uncal.search(type=bool)  # doctest: +SKIP
    root (AsdfObject)
    └─roman (WfiScienceRaw)
      └─meta (dict)
        ├─exposure (Exposure)
        │ └─data_problem (bool): False
        └─visit (Visit)
          ├─internal_target (bool): False
          └─target_of_opportunity (bool): False

    >>> d_uncal.search(type=bool, value=True)  # doctest: +SKIP
    No results found.

More information and options for the search method can be found in the
ASDF documentation
`here. <https://asdf.readthedocs.io/en/stable/asdf/features.html#searching-the-asdf-tree>`_
