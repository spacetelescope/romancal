Examples
========
In the examples below, `img` is either a string with the filename of a Roman `ASDF` file
or a Roman datamodel `ImageModel`.

#. To run TweakReg in a Python session with the default parameters:

        .. code-block:: python

                from romancal.tweakreg.tweakreg_step import TweakRegStep
                step = TweakRegStep()
                step.run([img])

        .. note::
            If the input is a single Roman ``DataModel``,
            either ``step.run([img])`` or ``step.run(img)`` will work. For multiple elements as input,
            they must be passed in as either a list or a ModelLibrary.

#. To run TweakReg in a Python session on an association file with the default parameters:

        .. code-block:: python

                from romancal.tweakreg.tweakreg_step import TweakRegStep
                step = TweakRegStep()
                step.run("asn_file.json")

#. To run TweakReg on a Roman's exposure with default astrometric parameters and save
   the absolute catalog data:

        .. code-block:: python

                from romancal.tweakreg.tweakreg_step import TweakRegStep
                step = TweakRegStep()
                step.save_abs_catalog = True # save the catalog data used for absolute astrometry
                step.abs_refcat = 'GAIADR3' # use Gaia DR3 for absolute astrometry
                step.catalog_path = '/path/for/the/abs/catalog' # save the Gaia catalog to this path
                step.run([img])


