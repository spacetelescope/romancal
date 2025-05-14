.. _metadata:


Metadata
========

Metadata information associated with a data model is accessed through
its `meta` member.  For example, to access the date that an
observation was made::

    print(model.meta.observation.start_time)

Metadata values will be validated against the schema when ``validate``
is called, when the data model is saved or when a file is read.

.. code-block:: python

        >>> from roman_datamodels import datamodels as rdm

        >>> model = rdm.ImageModel.fake_data()
        >>> model.meta.pointing.target_ra = "foo"
        >>> model.validate()  # doctest: +IGNORE_EXCEPTION_DETAIL
        Traceback (most recent call last):
          File "<stdin>", line 1, in <module>
          File "/path/to/python/lib/python3.12/site-packages/roman_datamodels/stnode/_node.py", line 251, in __setattr__
            _validate(key, value, schema, self.ctx)
          File "/path/to/python/lib/python3.12/site-packages/roman_datamodels/stnode/_node.py", line 78, in _validate
            return _value_change(attr, tagged_tree, schema, False, will_strict_validate(), ctx)
                  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
          File "/path/to/python/lib/python3.12/site-packages/roman_datamodels/stnode/_node.py", line 47, in _value_change
            raise ValidationError(errmsg)
        asdf._jsonschema.exceptions.ValidationError: While validating target_ra the following error occurred:
        'foo' is not of type 'number'

        Failed validating 'type' in schema:
            {'$schema': 'http://stsci.edu/schemas/asdf-schema/0.1.0/asdf-schema',
            'archive_catalog': {'datatype': 'float',
                                'destination': ['WFIExposure.target_ra',
                                                'GuideWindow.target_ra']},
            'description': 'Right ascension in units of degrees at the location '
                            'of\n'
                            'the target aperture.\n',
            'sdf': {'source': {'origin': 'TBD'},
                    'special_processing': 'VALUE_REQUIRED'},
            'title': 'Right Ascension of the Target Aperture (deg)',
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

    >>> raw_science = mk_datamodel(rdm.ScienceRawModel)  # Create a model of the desired type
    >>> print(raw_science.schema_uri)  # find the associated Schema # doctest: +SKIP
    asdf://stsci.edu/datamodels/roman/schemas/wfi_science_raw-1.0.0


An alternative method to get and set metadata values is to use a
dot-separated name as a dictionary lookup.  This is useful for
databases, such as CRDS, where the path to the metadata element is
most conveniently stored as a string.  The following two lines are
equivalent

.. code-block:: python

    >>> print(raw_science.meta['visit']['start_time'])
    2020-01-01T00:00:00.000
    >>> print(raw_science.meta.visit.start_time)
    2020-01-01T00:00:00.000

In addition the times are stored as Astropy time objects and so the date can be
displayed using various formats

.. code-block:: python

    >>> print(raw_science.meta.visit.start_time.iso)
    2020-01-01 00:00:00.000
    >>> print(raw_science.meta.visit.start_time.mjd)
    58849.0
    >>> print(raw_science.meta.visit.start_time.yday)
    2020:001:00:00:00.000
