.. _asn-from-list:

asn_from_list
=============

Create an association using either the command line tool
``asn_from_list`` or through the Python API using
:func:`romancal.associations.asn_from_list.asn_from_list`


Associations
^^^^^^^^^^^^

Refer to TBD for a full description of associations.

To create an association, use the following command:

.. code-block:: python

 import romancal.associations.asn_from_list as asn_from_list
 product_name = 'test_product'
 items = {'r0000101001001001001_0001_wfi01_uncal.asdf': 'science', 'r0000101001001001001_3_0001_wfi01_uncal.asdf': 'guide_window', 'c': 'somethingelse'}
 asn = asn_from_list.asn_from_list([(item, type_) for item, type_ in items.items()], product_name=product_name, with_exptype=True)
 asn['asn_rule']
 'DMS_ELPP_Base'


an example product that has both a science and guide window exposures
would look like the following::

    asn['products']
    [   {   'members': [   {   'expname': 'r0000101001001001001_0001_wfi01_uncal.asdf',
                               'exptype': 'science'},
                           {   'expname': 'r0000101001001001001_3_0001_wfi01_uncal.asdf',
                               'exptype': 'guide_window'},
                           {'expname': 'c', 'exptype': 'somethingelse'}],
            'name': 'test_product'}]


To create a association with all the detectors for a given exposure from the command line,

.. code-block:: python

		asn_from_list -o detector_asn.json --product-name r0000101001001001001_0001_wfi data/*_cal.asdf

where the individual calibrated detector files are in a data subdirectory.

In addition you can also specify a target for the association file,

.. code-block:: python

		asn_from_list -o level3_mosaic_asn.json --product-name level3_mosaic --target r274dp63x32y80  data/*_cal.asdf
