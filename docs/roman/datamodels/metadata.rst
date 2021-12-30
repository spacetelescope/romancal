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

        from romancal.datamodels import ImageModel
        model = ImageModel()
        model.meta.target.ra = "foo"
        Traceback (most recent call last):
          File "<stdin>", line 1, in <module>
          File "site-packages/romancal.datamodels/schema.py", line 672, in __setattr__
            object.__setattr__(self, attr, val)
          File "site-packages/romancal.datamodels/schema.py", line 490, in __set__
            val = self.to_basic_type(val)
          File "site-packages/romancal.datamodels/schema.py", line 422, in to_basic_type
            raise ValueError(e.message)
        ValueError: 'foo' is not of type u'number'

The set of available metadata elements is defined in a YAML Schema
that is installed with `roman_datamodels` from the `RAD` (Roman Attribute
Dictionary).

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

In addtion the times are stored as python time objects and so the date can be
displayed using various formats::

    print(raw_science.meta.observation.start_time.iso)
    2028-12-22 05:17:56.203
    print(raw_science.meta.observation.start_time.mjd)
    62127.22078938165
    print(raw_science.meta.observation.start_time.yday)
    2028:357:05:17:56.203
