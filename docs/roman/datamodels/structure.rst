The Structure of DataModels
===========================

Datamodels allows for the creation of a new model though the usual
method of calling the __init__ method. Each type of model has its own
class and schema. The schema is specified through the class variable
schema_url. The schema gives the binding between the metadata
and/or extension and the datamodels attribute name. The chief
distinction between the two is that ASDF has a flat model, datamodels
supports a hierarchical data model. The typical structure of a
datamodels class is that it first calls the __init__ method of the
base class and then initializes the required arrays of the models with
lines that look like they shouldn't do anything, for example

    self.dq = self.dq

The reason why a line like the above initializes the array is that the
access to the array on the right side of the assignment will
initialize the array to a default value if it is not already defined
and one is found in the schema. The sequence of calls is that the dot
notation invokes __getattr__, which calls _make_default if the
attribute is not defined, which in turn calls _make_default_array if
the schema says the attribute is an array. All these methods can be
found in properties.py.

The base class for Datamodels is RomanDataModel, in model_base.py. It takes
several arguments, the most important of which is init, which as the
name suggests, specifies how to initialize the primary data array of
the model. Init is most usually the name of a file, but can be an
already opened asdf file, a numpy array, a shape tuple, or
None. If init is a shape tuple the primary data array is initialized
to its default value.

Optional arguments to __init__  can give a schema which overrides the
class schema, extensions to the schema, two flages pass_invalid_values
and strict_validation, which control the data validation, and numpy arrays
which are used to initialized the model arrays by using parameters of the
same name.

As an alternative to creating a model by initializing an object of the
specific class, you can call the open function, which is in
util.py. This function takes the same arguments as the __init_
method. If it is called with the name of aa ASDF file, it looks in the
metadata for a attribute named datamodel that contains the name of
the class to use to open the model. If that attribute is not found,
checks the dimensionality of the image and uses a generic model type
to open the image.

The base class for Datamodels loads the schema from the a file in the
schemas subdirectory. If the base class is passed a descriptor of an
already open model, it returns a shallow copy of the already open
image. This is done to speed the code, as re-opening already open
models is a common operation in the pipeline. If it is passed the
name of a file, it peeks at the first several bytes of the file to
determine the file type. This test is in filetype.py.

To write a model back to a file, call the save method on the file. It
first calls validate_required to check the schema to see if all the
required fields are present in the model. (More information needed)

Items within a model are accessed as attributes, that is, with dot
notation. The code which handles getting and setting attributes is
found in properties.py. Datamodels distinguishes between items at the
endpoints of the asdf tree and subtrees within the asdf tree. The
former are returned as scalars or numpy arrays, depending on whether
the endpoint represents a attribute or data array. Subtrees are
returned as nodes. A node is an object containing the subtree as well
as the subschema which describes the subtree.  If one tries to get an
attribute that does not exist in the asdf tree, one of several things
may happen. If the attribute is not mentioned in the schema, the value
of the attribute is set to None. If it is in the schema and the schema
has a default value, the code creates the item with the default value
and then returns it. The functions that do this are _make_default and
_make_default_array, which it calls. If not only the item, but the
subtree containing the item is missing, the code throws an
AttributeError. When an attribute representing an array is accessed,
the type of the array is compared to the type in the schema and if
they are different, the array is cast to the type in the schema.

When setting or deleting an attribute, the code validates the
change. The code which does the validation can be found in
validate.py. The validator checks the values of pass_invalid_values,
which allows values inconsistent with the schema to be set, and
strict_validation, which throws an exception if the value does not
match the schema.
