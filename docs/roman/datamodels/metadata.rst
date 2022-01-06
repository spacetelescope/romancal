.. _metadata:


Metadata
========

Metadata information associated with a data model is accessed through
its `meta` member.  For example, to access the date that an
observation was made::

    print(model.meta.observation.start_time)

Metadata values are automatically type-checked against the schema when
they are set. Therefore, setting a attribute which expects a number to a
string will raise an exception.

.. code-block:: python

        from roman_datamodels.testing.factories import create_wfi_image
        from roman_datamodels import datamodels as rdmfrom romancal.datamodels import ImageModel
        model = rdm.ImageModel(create_wfi_image())
        model.meta.target.ra = "foo"
        Traceback (most recent call last):
          File "<stdin>", line 1, in <module>
          File "/Users/ddavis/miniconda3/envs/rcal_dev/lib/python3.9/site-packages/roman_datamodels/stnode.py", line 183, in __setattr__
            if schema is None or _validate(key, value, schema, self.ctx):
          File "/Users/ddavis/miniconda3/envs/rcal_dev/lib/python3.9/site-packages/roman_datamodels/stnode.py", line 97, in _validate
            return _value_change(attr, tagged_tree, schema, False, strict_validation, ctx)
          File "/Users/ddavis/miniconda3/envs/rcal_dev/lib/python3.9/site-packages/roman_datamodels/stnode.py", line 68, in _value_change
            raise jsonschema.ValidationError(errmsg)
        jsonschema.exceptions.ValidationError: While validating ra the following error occurred:
        'foo' is not of type 'number'

        Failed validating 'type' in schema:
            {'$schema': 'http://stsci.edu/schemas/asdf-schema/0.1.0/asdf-schema',
             'archive_catalog': {'datatype': 'float',
                                 'destination': ['ScienceCommon.ra']},
             'sdf': {'source': {'origin': 'PSS:fixed_target.ra_computed'},
                     'special_processing': 'VALUE_REQUIRED'},
             'title': 'Target RA at mid time of exposure',
             'type': 'number'}

        On instance:
            'foo'

The set of available metadata elements is defined in a YAML Schema
that is installed with `roman_datamodels <https://github.com/spacetelescope/roman_datamodels>`_
from the
`RAD <https://github.com/spacetelescope/RAD>`_ (Roman Attribute Dictionary).

There is also a utility method for finding the schema associated with a given
model.

.. code-block:: python

    from roman_datamodels import datamodels as rdm
    from roman_datamodels.testing.factories import create_wfi_science_raw
    # Create a model of the desired type
    raw = create_wfi_science_raw()
    raw_science = rdm.ScienceRawModel(raw)
    # find the associated Schema
    raw_science.schema_uri
    'asdf://stsci.edu/datamodels/roman/schemas/wfi_science_raw-1.0.0'


An alternative method to get and set metadata values is to use a
dot-separated name as a dictionary lookup.  This is useful for
databases, such as CRDS, where the path to the metadata element is
most conveniently stored as a string.  The following two lines are
equivalent::

    print(raw_science.meta['observation']['start_time'])
    print(raw_science.meta.observation.start_time)

In addtion the times are stored as Astropy time objects and so the date can be
displayed using various formats::

    print(raw_science.meta.observation.start_time.iso)
    2028-12-22 05:17:56.203
    print(raw_science.meta.observation.start_time.mjd)
    62127.22078938165
    print(raw_science.meta.observation.start_time.yday)
    2028:357:05:17:56.203
