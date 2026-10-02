.. _data_products_versioning:

Versioning
----------

Roman ASDF data products are highly versioned, allowing for careful tracking of datamodel changes over time. With each new version of romancal it is expected that the versions of datamodels will advance. The changes for each new version may be minor, or may have significant impact on scientific analysis. Processing files from several versions requires understanding how the file contents differ and what options romancal provides for handling different versions of files.

.. _data_products_backwards_compatibility:

Backward Compatibility
^^^^^^^^^^^^^^^^^^^^^^

By default, romancal expects to process the latest datamodel version(s) at the time of release. Attempting to process an older file may require using the `intro_update_version` parameter which will attempt to update old files to the newest version (providing warnings about which migrations were performed).

.. _data_products_forward_compatibility:

Forward Compatibility
^^^^^^^^^^^^^^^^^^^^^

Processing a newer-than-known datamodel version may be possible. However, unlike with `data_products_backwards_compatibility` the old version of romancal knows nothing about the newer version and is unable to perform migrations that are informed by known datamodel changes. The preferred option in this situation is to upgrade the environment to use a newer romancal that is aware of the datamodel version. This is likely to pull in other bugfixes, improvements and features and is highly recommended.

The `intro_downgrade_version` parameter can be used to perform a basic conversion of the newer-than-known datamodel to the latest known version. This can fail (if the versions are incompatible) but can also succeed and produce incorrect results. It is recommended that users of this feature become familiar with the differences between the versions to understand how this will impact interpretation of these files.
