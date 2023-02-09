.. _asn-from-list:

asn_from_list
=============

Create an association using either the command line tool
``asn_from_list`` or through the Python API using either
:class:`romancal.associations.asn_from_list.Main` or
:func:`romancal.associations.asn_from_list.asn_from_list`


Level2 Associations
^^^^^^^^^^^^^^^^^^^

Refer to TBD for a full description of Level2
associations.

To create a Level2 association, use the following command:

.. code-block:: python

 import romancal.associations.asn_from_list as asn_from_list
 product_name = 'test_product'
 items = {'r0000101001001001001_01101_0001_WFI01_uncal.asdf': 'science', 'r0000101001001001001_3_01101_0001_WFI01_uncal.asdf': 'guide_window', 'c': 'somethingelse'}
 asn = asn_from_list.asn_from_list([(item, type_) for item, type_ in items.items()], product_name=product_name, with_exptype=True)
 asn['asn_rule']
 'DMS_ELPP_Base'


an example product that has both a science and guide window exposures
would look like the following::

    asn['products']
    [   {   'members': [   {   'expname': 'r0000101001001001001_01101_0001_WFI01_uncal.asdf',
                               'exptype': 'science'},
                           {   'expname': 'r0000101001001001001_3_01101_0001_WFI01_uncal.asdf',
                               'exptype': 'guide_window'},
                           {'expname': 'c', 'exptype': 'somethingelse'}],
            'name': 'test_product'}]


To create a association with all the detectors for a given exposure from the command line,

.. code-block:: python

		asn_from_list -o detector_asn.json --product-name r0000101001001001001_01101_0001_WFI_cal.asdf data/*_cal.asdf

where the individual calibrated detector files are in a data subdirectory.
