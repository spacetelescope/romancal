.. _asn-from-list:

asn_from_list
=============
Create an association using either the command line tool
``asn_from_list`` or through the Python API using
:func:`romancal.associations.asn_from_list.asn_from_list`

Associations
^^^^^^^^^^^^
Refer to TBD for a full description of associations.

To create an association programmatically, use:

.. code-block:: python

    import romancal.associations.asn_from_list as asn_from_list
    product_name = 'test_product'
    items = {
        'r0000101001001001001_0001_wfi01_uncal.asdf': 'science',
        'r0000101001001001001_3_0001_wfi01_uncal.asdf': 'guide_window',
        'c': 'somethingelse'
    }
    # Items should be a list of (filename, exptype) tuples if using with_exptype=True
    asn = asn_from_list.asn_from_list(
        [(item, type_) for item, type_ in items.items()],
        product_name=product_name,
        with_exptype=True
    )
    asn['asn_rule']
    # 'DMS_ELPP_Base'

An example product that has both science and guide window exposures would look like:

.. code-block:: python

    asn['products']
    [   {   'members': [   {   'expname': 'r0000101001001001001_0001_wfi01_uncal.asdf',
                               'exptype': 'science'},
                           {   'expname': 'r0000101001001001001_3_0001_wfi01_uncal.asdf',
                               'exptype': 'guide_window'},
                           {'expname': 'c', 'exptype': 'somethingelse'}],
            'name': 'test_product'}]


To create an association with all the detectors for a given exposure from the command line:

.. code-block:: shell

    asn_from_list -o detector_asn.json --product-name r0000101001001001001_0001_wfi data/*_cal.asdf

where the individual calibrated detector files are in a data subdirectory.

You can also specify a target for the association file:

.. code-block:: shell

    asn_from_list -o level3_mosaic_asn.json --product-name level3_mosaic --target 270p65x48y69 data/*_cal.asdf

To specify a data release ID in the metadata and filename:

.. code-block:: shell

    asn_from_list -o my_asn.json --product-name my_product --data-release-id r0 data/*_cal.asdf

Notes
^^^^^
- The CLI expects a list of filenames as positional arguments.
- The Python API allows more advanced usage, such as passing (filename, exptype) tuples.
- The `"data_release_id"` field will always appear immediately after `"program"` in the output ASN metadata for consistency.
- The `--rule` and `--ruledefs` options allow you to specify custom association rules if needed.
