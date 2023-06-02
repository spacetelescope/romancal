Examples
========
In the examples below, `img` is either a string with the filename of a Roman `ASDF` file
or a Roman datamodel `ImageModel`.

1. To run TweakReg on a python session one image file
   (i.e. one Roman's SCA) with the default parameters:

        .. code-block:: python

                from romancal.tweakreg.tweakreg_step import TweakRegStep
                step = TweakRegStep()
                step.process([img])

2. To run TweakReg on a Roman's exposure with default astrometric parameters and save
   the absolute catalog data:

        .. code-block:: python

                from romancal.tweakreg.tweakreg_step import TweakRegStep
                step = TweakRegStep()
                step.save_abs_catalog = True # save the catalog data used for absolute astrometry
                step.abs_refcat = 'GAIADR3' # use Gaia DR3 for absolute astrometry
                step.catalog_path = '/path/for/the/abs/catalog' # save the Gaia catalog to this path
                step.process([img])

3. To run TweakReg using a custom source catalog with the default parameters:

   - make sure the source catalog is in one of the supported formats. The file content
     below is an example of a valid catalog (`ascii.csv` format):

        .. code-block:: bash

                $ cat ref_catalog_1

                x,y
                846.1321662446178,945.839358133909
                664.7073537074112,1028.4613139252003
                1036.160742774408,642.3379043578552
                935.8827367579428,594.1745467413945
                1491.9672737821606,1037.4723609624757
                735.1256651803337,1410.2791591559157
                1358.2876707625007,651.7112260833995
                526.4715950130742,751.745104066621
                1545.082698426152,703.601696337681
                374.9609365496525,972.6561578187437
                1110.3498547121228,1644.2214966576498
                341.18333252240654,891.4733849441861
                820.0520846885105,312.0088351823117
                567.7054174813052,386.8883078361564
                1447.356249085851,1620.3390168916592
                1400.4271386280673,1674.3765672924937
                1681.9744852889235,571.6748779060324
                959.7317254404431,197.8757865066898
                1806.3360866990297,769.0603031839573
                487.1560001146406,257.30706691141086
                1048.7910126076483,85.36675265982751
                1075.508595999755,29.085099663125334

   - create `catfile` containing the filename of the input Roman datamodel and
     its corresponding catalog, one per line, as shown below

        .. code-block:: bash

                $ cat /path/to/catfile/catfilename

                img1 ref_catalog_1
                img2 ref_catalog_2
                img3 ref_catalog_3

   The content of `catfile` will allow TweakReg to assign the custom catalog to the
   correct input Roman datamodel. In the example above, source catalog
   `ref_catalog_1` will be assign to `img1`, and so on.

   Now we can execute the following:

        .. code-block:: python

                from romancal.tweakreg.tweakreg_step import TweakRegStep
                step = TweakRegStep()
                step.use_custom_catalogs = True # use custom catalogs
                step.catalog_format = "ascii.ecsv" # custom catalogs format
                step.catfile = '/path/to/catfile/catfilename' # path to datamodel:catalog mapping
                step.process([img])
