=====
Steps
=====

.. _writing-a-step:

Writing a step
==============

Writing a new step involves writing a class that has a `process`
method to perform work and a `spec` member to define its configuration
parameters.  (Optionally, the `spec` member may be defined in a
separate `spec` file).

Inputs and outputs
------------------

A `Step` provides a full framework for handling I/O.

Steps get their inputs from two sources:

    - Configuration parameters come from the parameter file or the
      command line and are set as member variables on the Step object
      by the stpipe framework.

    - Arguments are passed to the Step’s `process` function as regular
      function arguments.

Parameters should be used to specify things that must be determined outside of
the code by a user using the class. Arguments should be used to pass data that
needs to go from one step to another as part of a larger pipeline. Another way
to think about this is: if the user would want to examine or change the value,
use a parameter.

The parameters are defined by the :ref:`Step.spec <the-spec-member>` member.

Input Files, Associations, and Directories
``````````````````````````````````````````

It is presumed that all input files are co-resident in the same
directory. This directory is whichever directory the first input file
is found in. This is particularly important for associations. It is
assumed that all files referenced by an association are in the same
directory as the association file itself.

Output Files and Directories
````````````````````````````

The step will generally return its output as a data model. Every step has
implicitly created parameters `output_dir` and `output_file` which the user can
use to specify the directory and file to save this model to. Since the `stpipe`
architecture generally creates output file names, in general, it is expected
that `output_file` be rarely specified, and that different sets of outputs be
separated using `output_dir`.

Output Suffix
-------------

There are three ways a step's results can be written to a file:

1. Implicitly when a step is run from the command line or with
   `Step.from_cmdline`

2. Explicitly by specifying the parameter `save_results`

3. Explicitly by specifying a file name with the parameter
   `output_file`

In all cases, the file, or files, is/are created with an added suffix
at the end of the base file name. By default this suffix is the class
name of the step that produced the results. Use the `suffix` parameter
to explicitly change the suffix.


The Python class
----------------

At a minimum, the Python Step class should inherit from `stpipe.Step`, implement
a ``process`` method to do the actual work of the step and have a `spec` member
to describe its parameters.

1. Objects from other Steps in a pipeline are passed as arguments to
   the ``process`` method.

2. The parameters described in :ref:`configuring-a-step`
   are available as member variables on ``self``.

3. To support the caching suspend/resume feature of pipelines, images
   must be passed between steps as model objects.  To ensure you’re
   always getting a model object, call the model constructor on the
   parameter passed in.  It is good idea to use a `with` statement
   here to ensure that if the input is a file path that the file will
   be appropriately closed.

4. Objects to pass to other Steps in the pipeline are simply returned
   from the function.  To return multiple objects, return a tuple.

5. The parameters for the step are described in the `spec` member in the
   `configspec` format.

6. Declare any CRDS reference files used by the Step.  (See
   :ref:`interfacing-with-crds`).

::

    from romancal.stpipe import RomanStep

    from roman_datamodels.datamodels import ImageModel
    from my_awesome_astronomy_library import combine

    class ExampleStep(RomanStep):
        """
        Every step should include a docstring.  This docstring will be
        displayed by the `strun --help`.
        """

        # 1.
        def process(self, image1, image2):
            self.log.info("Inside ExampleStep")

            # 2.
            threshold = self.threshold

            # 3.
            with ImageModel(image1) as image1, ImageModel(image2) as image2:
                # 4.
                with self.get_reference_file_model(image1, "flat_field") as flat:
                    new_image = combine(image1, image2, flat, threshold)

            # 5.
            return new_image

       # 6.
       spec = """
       # This is the configspec file for ExampleStep

       threshold = float(default=1.0)  # maximum flux
       """

       # 7.
       reference_file_types = ['flat_field']

The Python Step subclass may be installed anywhere that your Python
installation can find it.  It does not need to be installed in the
`stpipe` package.

.. _the-spec-member:

The spec member
---------------

The `spec` member variable is a string containing information about
the parameters.  It is in the `configspec` format
defined in the `ConfigObj` library that stpipe uses.

The `configspec` format defines the types of the parameters, as well as allowing
an optional tree structure.

The types of parameters are declared like this::

    n_iterations = integer(1, 100)  # The number of iterations to run
    factor = float()                # A multiplication factor
    author = string()               # The author of the file

Note that each parameter may have a comment.  This comment is
extracted and displayed in help messages and docstrings etc.

Parameters can be grouped into categories using
ini-file-like syntax::

    [red]
    offset = float()
    scale = float()

    [green]
    offset = float()
    scale = float()

    [blue]
    offset = float()
    scale = float()

Default values may be specified on any parameter using the `default`
keyword argument::

    name = string(default="John Doe")

While the most commonly useful parts of the configspec format are
discussed here, greater detail can be found in the `configspec
documentation
<https://configobj.readthedocs.io/en/latest/>`_.

Configspec types
````````````````

The following is a list of the commonly useful configspec types.

    `integer`: matches integer values. Takes optional `min` and `max`
    arguments::

        integer()
        integer(3, 9)  # any value from 3 to 9
        integer(min=0) # any positive value
        integer(max=9)

    `float`: matches float values Has the same parameters as the
    integer check.

    `boolean`: matches boolean values: True or False.

    `string`: matches any string. Takes optional keyword args `min`
    and `max` to specify min and max length of string.

    `list`: matches any list. Takes optional keyword args `min`, and
    `max` to specify min and max sizes of the list. The list checks
    always return a list.

    `force_list`: matches any list, but if a single value is passed in
    will coerce it into a list containing that value.

    `int_list`: Matches a list of integers. Takes the same arguments
    as list.

    `float_list`: Matches a list of floats. Takes the same arguments
    as list.

    `bool_list`: Matches a list of boolean values. Takes the same
    arguments as list.

    `string_list`: Matches a list of strings. Takes the same arguments
    as list.

    `option`: matches any from a list of options. You specify this
    test with::

        option('option 1', 'option 2', 'option 3')

    Normally, steps will receive input files as parameters and return
    output files from their process methods.  However, in cases where
    paths to files should be specified in the parameter file,
    there are some extra parameter types that stpipe provides that
    aren’t part of the core configobj library.

    `input_file`: Specifies an input file.  Relative paths are
    resolved against the location of the parameter file.  The file
    must also exist.

    `output_file`: Specifies an output file.  Identical to
    `input_file`, except the file doesn’t have to already exist.

.. _interfacing-with-crds:

Interfacing with CRDS
---------------------

If a Step uses CRDS to retrieve reference files, there are two
things to do:

1. Within the `process` method, call `self.get_reference_file` or
   `self.get_reference_file_model` to get a reference file from CRDS.
   These methods take as input a) a model for the input file, whose
   metadata is used to do a CRDS bestref lookup, and b) a reference
   file type, which is just a string to identify the kind of reference
   file.

2. Declare the reference file types used by the Step in the
   `reference_file_types` member. This information is used by the stpipe
   framework for two purposes: a) to pre-cache the reference files needed by a
   Pipeline before any of the pipeline processing actually runs, and b) to add
   override parameters to the Step's configspec.

For each reference file type that the Step declares, an `override_*` parameter
is added to the Step's configspec. For example, if a step declares the
following::

   reference_file_types = ['flat_field']

then the user can override the flat field reference file using the
parameter file::

   override_flat_field = /path/to/my_reference_file.asdf

or at the command line::

   --override_flat_field=/path/to/my_reference_file.asdf
