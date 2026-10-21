Populate the ``image_flags`` column by combining the ``dq`` bits of every pixel
in a source segment, and add the corresponding ``image_flags2`` column for
``dq2``. Columns without a schema definition no longer raise when the catalog
is built.
