"""
Roman pipeline step for image alignment.
"""
from os import path

from astropy.table import Table
from astropy import units as u
from astropy.coordinates import SkyCoord
from tweakwcs.imalign import align_wcs
from tweakwcs.correctors import JWSTWCSCorrector
from tweakwcs.matchutils import XYXYMatch

# LOCAL
from ..stpipe import RomanStep
from roman_datamodels import datamodels
from . import astrometric_utils as amutils
from roman_datamodels.util import is_association


def _oxford_or_str_join(str_list):
    nelem = len(str_list)
    if not nelem:
        return 'N/A'
    str_list = list(map(repr, str_list))
    if nelem == 1:
        return str_list
    elif nelem == 2:
        return str_list[0] + ' or ' + str_list[1]
    else:
        return ', '.join(map(repr, str_list[:-1])) + ', or ' + repr(str_list[-1])


SINGLE_GROUP_REFCAT = ['GAIADR3', 'GAIADR2', 'GAIADR1']
_SINGLE_GROUP_REFCAT_STR = _oxford_or_str_join(SINGLE_GROUP_REFCAT)
DEFAULT_ABS_REFCAT = SINGLE_GROUP_REFCAT[0]
ALIGN_TO_ABS_REFCAT = True

__all__ = ['TweakRegStep']


class TweakRegStep(RomanStep):
    """
    TweakRegStep: Image alignment based on catalogs of sources detected in
    input images.
    """

    class_alias = "tweakreg"

    spec = f"""
        save_catalogs = boolean(default=False) # Write out catalogs?
        use_custom_catalogs = boolean(default=False) # Use custom user-provided catalogs?
        catalog_format = string(default='ecsv') # Catalog output file format
        catfile = string(default='') # Name of the file with a list of custom user-provided catalogs
        enforce_user_order = boolean(default=False) # Align images in user specified order?
        expand_refcat = boolean(default=False) # Expand reference catalog with new sources?
        minobj = integer(default=15) # Minimum number of objects acceptable for matching
        searchrad = float(default=2.0) # The search radius in arcsec for a match
        use2dhist = boolean(default=True) # Use 2d histogram to find initial offset?
        separation = float(default=1.0) # Minimum object separation in arcsec
        tolerance = float(default=0.7) # Matching tolerance for xyxymatch in arcsec
        fitgeometry = option('shift', 'rshift', 'rscale', 'general', default='rshift') # Fitting geometry
        nclip = integer(min=0, default=3) # Number of clipping iterations in fit
        sigma = float(min=0.0, default=3.0) # Clipping limit in sigma units
        abs_refcat = string(default='{DEFAULT_ABS_REFCAT}')  # Catalog file name or one of:
        # {_SINGLE_GROUP_REFCAT_STR}, or None, or ''
        save_abs_catalog = boolean(default=False)  # Write out used absolute astrometric reference catalog as a separate product
        abs_minobj = integer(default=15) # Minimum number of objects acceptable for matching when performing absolute astrometry
        abs_searchrad = float(default=6.0) # The search radius in arcsec for a match when performing absolute astrometry
        # We encourage setting this parameter to True. Otherwise, xoffset and yoffset will be set to zero.
        abs_use2dhist = boolean(default=True) # Use 2D histogram to find initial offset when performing absolute astrometry?
        abs_separation = float(default=0.1) # Minimum object separation in arcsec when performing absolute astrometry
        abs_tolerance = float(default=0.7) # Matching tolerance for xyxymatch in arcsec when performing absolute astrometry
        # Fitting geometry when performing absolute astrometry
        abs_fitgeometry = option('shift', 'rshift', 'rscale', 'general', default='rshift')
        abs_nclip = integer(min=0, default=3) # Number of clipping iterations in fit when performing absolute astrometry
        abs_sigma = float(min=0.0, default=3.0) # Clipping limit in sigma units when performing absolute astrometry
        output_use_model = boolean(default=True)  # When saving use `DataModel.meta.filename`
    """

    reference_file_types = []

    def process(self, input):
        use_custom_catalogs = self.use_custom_catalogs

        if use_custom_catalogs:
            catdict = _parse_catfile(self.catfile)
            # if user requested the use of custom catalogs and provided a
            # valid 'catfile' file name that has no custom catalogs,
            # turn off the use of custom catalogs:
            if catdict is not None and not catdict:
                self.log.warning(
                    "'use_custom_catalogs' is set to True but 'catfile' "
                    "contains no user catalogs. Turning on built-in catalog "
                    "creation."
                )
                use_custom_catalogs = False

        try:
            if use_custom_catalogs and catdict:
                images = datamodels.ModelContainer()
                if isinstance(input, str):
                    asn_dir = path.dirname(input)
                    asn_data = images.read_asn(input)
                    for member in asn_data['products'][0]['members']:
                        filename = member['expname']
                        member['expname'] = path.join(asn_dir, filename)
                        if filename in catdict:
                            member['tweakreg_catalog'] = catdict[filename]
                        elif 'tweakreg_catalog' in member:
                            del member['tweakreg_catalog']

                    images.from_asn(input)

                elif is_association(input):
                    images.from_asn(input)

                else:
                    images = datamodels.ModelContainer(input)
                    for im in images:
                        filename = im.meta.filename
                        if filename in catdict:
                            print(f"setting {filename}.tweakreg_catalog = {repr(catdict[filename])}")
                            im.meta.tweakreg_catalog = catdict[filename]

            else:
                images = datamodels.ModelContainer(input)

        except TypeError as e:
            e.args = ("Input to tweakreg must be a list of DataModels, an "
                      "association, or an already open ModelContainer "
                      "containing one or more DataModels.", ) + e.args[1:]
            raise e

        if self.abs_refcat is not None and self.abs_refcat.strip() and self.abs_refcat != DEFAULT_ABS_REFCAT:
            # Set expand_refcat to True to eliminate possibility of duplicate
            # entries when aligning to absolute astrometric reference catalog
            self.expand_refcat = True

        if len(images) == 0:
            raise ValueError("Input must contain at least one image model.")

        # Build the catalogs for input images
        for i, image_model in enumerate(images):
            image_model.meta.wcsinfo["WCSAXES"] = image_model.meta.wcs.world_n_dim
            try:
                # use user-supplied catalog:
                self.log.info("Using user-provided input catalog "
                              f"'{image_model.meta.tweakreg_catalog}'")
                catalog = Table.read(image_model.meta.tweakreg_catalog)
                new_cat = False

            except OSError:
                self.log.error("Failed to read 'meta.tweakreg_catalog' from source detection step ")

            for axis in ['x', 'y']:
                if axis not in catalog.colnames:
                    long_axis = axis + 'centroid'
                    if long_axis in catalog.colnames:
                        catalog.rename_column(long_axis, axis)
                    else:
                        raise ValueError(
                            "'tweakreg' source catalogs must contain either "
                            "columns 'x' and 'y' or 'xcentroid' and "
                            "'ycentroid'."
                        )

            # filter out sources outside the WCS bounding box
            bb = image_model.meta.wcs.bounding_box
            if bb is not None:
                ((xmin, xmax), (ymin, ymax)) = bb
                x = catalog['x']
                y = catalog['y']
                mask = (x > xmin) & (x < xmax) & (y > ymin) & (y < ymax)
                catalog = catalog[mask]

            filename = image_model.meta.filename
            nsources = len(catalog)
            if nsources == 0:
                self.log.warning(f'No sources found in {filename}.')
            else:
                self.log.info('Detected {} sources in {}.'
                              .format(len(catalog), filename))

            if new_cat and self.save_catalogs:
                catalog_filename = filename.replace(
                    '.fits', f'_cat.{self.catalog_format}'
                )
                if self.catalog_format == 'ecsv':
                    fmt = 'ascii.ecsv'
                elif self.catalog_format == 'fits':
                    # NOTE: The catalog must not contain any 'None' values.
                    #       FITS will also not clobber existing files.
                    fmt = 'fits'
                else:
                    raise ValueError(
                        '\'catalog_format\' must be "ecsv" or "fits".'
                    )
                catalog.write(catalog_filename, format=fmt, overwrite=True)
                self.log.info('Wrote source catalog: {}'
                              .format(catalog_filename))
                image_model.meta.tweakreg_catalog = catalog_filename

            # Temporarily attach catalog to the image model so that it follows
            # the grouping by exposure, to be removed after use below
            image_model["catalog"] = catalog
            images[i] = image_model

        # group images by their "group id":
        grp_img = list(images.models_grouped)

        self.log.info('')
        self.log.info("Number of image groups to be aligned: {:d}."
                      .format(len(grp_img)))
        self.log.info("Image groups:")

        if len(grp_img) == 1 and not ALIGN_TO_ABS_REFCAT:
            self.log.info("* Images in GROUP 1:")
            for im in grp_img[0]:
                self.log.info(f"     {im.meta.filename}")
            self.log.info('')

            # we need at least two exposures to perform image alignment
            self.log.warning("At least two exposures are required for image "
                             "alignment.")
            self.log.warning("Nothing to do. Skipping 'TweakRegStep'...")
            self.skip = True
            for model in images:
                model.meta.cal_step["tweakreg"] = "SKIPPED"
                # Remove the attached catalogs
                del model.catalog
            return input

        elif len(grp_img) == 1 and ALIGN_TO_ABS_REFCAT:
            # create a list of WCS-Catalog-Images Info and/or their Groups:
            g = grp_img[0]
            if len(g) == 0:
                raise AssertionError("Logical error in the pipeline code.")
            group_name = _common_name(g)
            imcats = list(map(self._imodel2wcsim, g))
            # Remove the attached catalogs
            for model in g:
                model = model if isinstance(
                    model, datamodels.DataModel
                ) else datamodels.open(path.basename(model))
                if hasattr(model, "catalog"):
                    delattr(model, "catalog")
            self.log.info(f"* Images in GROUP '{group_name}':")
            for im in imcats:
                im.meta['group_id'] = group_name
                self.log.info("     {}".format(im.meta['name']))

            self.log.info('')

        elif len(grp_img) > 1:
            # create a list of WCS-Catalog-Images Info and/or their Groups:
            imcats = []
            for g in grp_img:
                if len(g) == 0:
                    raise AssertionError("Logical error in the pipeline code.")
                else:
                    group_name = _common_name(g)
                    wcsimlist = list(map(self._imodel2wcsim, g))
                    # Remove the attached catalogs
                    for model in g:
                        del model.catalog
                    self.log.info(f"* Images in GROUP '{group_name}':")
                    for im in wcsimlist:
                        im.meta['group_id'] = group_name
                        self.log.info("     {}".format(im.meta['name']))
                    imcats.extend(wcsimlist)

            self.log.info('')

            # align images:
            xyxymatch = XYXYMatch(
                searchrad=self.searchrad,
                separation=self.separation,
                use2dhist=self.use2dhist,
                tolerance=self.tolerance,
                xoffset=0,
                yoffset=0
            )

            try:
                align_wcs(
                    imcats,
                    refcat=None,
                    enforce_user_order=self.enforce_user_order,
                    expand_refcat=self.expand_refcat,
                    minobj=self.minobj,
                    match=xyxymatch,
                    fitgeom=self.fitgeometry,
                    nclip=self.nclip,
                    sigma=(self.sigma, 'rmse')
                )

            except ValueError as e:
                msg = e.args[0]
                if (msg == "Too few input images (or groups of images) with "
                        "non-empty catalogs."):
                    # we need at least two exposures to perform image alignment
                    self.log.warning(msg)
                    self.log.warning("At least two exposures are required for "
                                     "image alignment.")
                    self.log.warning("Nothing to do. Skipping 'TweakRegStep'...")
                    for model in images:
                        model.meta.cal_step.tweakreg = "SKIPPED"
                    if not ALIGN_TO_ABS_REFCAT:
                        self.skip = True
                        return images
                else:
                    raise e

            except RuntimeError as e:
                msg = e.args[0]
                if msg.startswith("Number of output coordinates exceeded allocation"):
                    # we need at least two exposures to perform image alignment
                    self.log.error(msg)
                    self.log.error("Multiple sources within specified tolerance "
                                   "matched to a single reference source. Try to "
                                   "adjust 'tolerance' and/or 'separation' parameters.")
                    self.log.warning("Skipping 'TweakRegStep'...")
                    self.skip = True
                    for model in images:
                        model.meta.cal_step.tweakreg = "SKIPPED"
                    return images
                else:
                    raise e

            for imcat in imcats:
                model = imcat.meta['image_model']
                if model.meta.cal_step.tweakreg == "SKIPPED":
                    continue
                wcs = model.meta.wcs
                twcs = imcat.wcs
                if not self._is_wcs_correction_small(wcs, twcs):
                    # Large corrections are typically a result of source
                    # mis-matching or poorly-conditioned fit. Skip such models.
                    self.log.warning(f"WCS has been tweaked by more than {10 * self.tolerance} arcsec")

                    for model in images:
                        model.meta.cal_step.tweakreg = "SKIPPED"
                    if ALIGN_TO_ABS_REFCAT:
                        self.log.warning("Skipping relative alignment (stage 1)...")
                    else:
                        self.log.warning("Skipping 'TweakRegStep'...")
                        self.skip = True
                        return images

        if ALIGN_TO_ABS_REFCAT:
            # Get catalog of GAIA sources for the field
            #
            # NOTE:  If desired, the pipeline can write out the reference
            #        catalog as a separate product with a name based on
            #        whatever convention is determined by the JWST Cal Working
            #        Group.
            if self.save_abs_catalog:
                output_name = f'fit_{self.abs_refcat.lower()}_ref.ecsv'
            else:
                output_name = None

            # initial shift to be used with absolute astrometry
            self.abs_xoffset = 0
            self.abs_yoffset = 0

            self.abs_refcat = self.abs_refcat.strip()
            gaia_cat_name = self.abs_refcat.upper()

            if gaia_cat_name in SINGLE_GROUP_REFCAT:
                ref_cat = amutils.create_astrometric_catalog(
                    images,
                    gaia_cat_name,
                    output=output_name
                )

            elif path.isfile(self.abs_refcat):
                ref_cat = Table.read(self.abs_refcat)

            else:
                raise ValueError("'abs_refcat' must be a path to an "
                                 "existing file name or one of the supported "
                                 f"reference catalogs: {_SINGLE_GROUP_REFCAT_STR}.")

            # Check that there are enough GAIA sources for a reliable/valid fit
            num_ref = len(ref_cat)
            if num_ref < self.abs_minobj:
                # Raise Exception here to avoid rest of code in this try block
                self.log.warning(
                    f"Not enough sources ({num_ref}) in the reference catalog "
                    "for the single-group alignment step to perform a fit. "
                    f"Skipping alignment to the {self.abs_refcat} reference "
                    "catalog!"
                )
            else:
                # align images:
                # Update to separation needed to prevent confusion of sources
                # from overlapping images where centering is not consistent or
                # for the possibility that errors still exist in relative overlap.
                xyxymatch_gaia = XYXYMatch(
                    searchrad=self.abs_searchrad,
                    separation=self.abs_separation,
                    use2dhist=self.abs_use2dhist,
                    tolerance=self.abs_tolerance,
                    xoffset=self.abs_xoffset,
                    yoffset=self.abs_yoffset
                )

                # Set group_id to same value so all get fit as one observation
                # The assigned value, 987654, has been hard-coded to make it
                # easy to recognize when alignment to GAIA was being performed
                # as opposed to the group_id values used for relative alignment
                # earlier in this step.
                for imcat in imcats:
                    imcat.meta['group_id'] = 987654
                    if ('fit_info' in imcat.meta and
                            'REFERENCE' in imcat.meta['fit_info']['status']):
                        del imcat.meta['fit_info']

                # Perform fit
                align_wcs(
                    imcats,
                    refcat=ref_cat,
                    enforce_user_order=True,
                    expand_refcat=False,
                    minobj=self.abs_minobj,
                    match=xyxymatch_gaia,
                    fitgeom=self.abs_fitgeometry,
                    nclip=self.abs_nclip,
                    sigma=(self.abs_sigma, 'rmse'),
                    ref_tpwcs=imcats[0]
                )

        for imcat in imcats:
            image_model = imcat.meta['image_model']
            image_model.meta.cal_step["tweakreg"] = 'COMPLETE'

            # retrieve fit status and update wcs if fit is successful:
            if 'SUCCESS' in imcat.meta.get('fit_info')['status']:

                # Update/create the WCS .name attribute with information
                # on this astrometric fit as the only record that it was
                # successful:
                if ALIGN_TO_ABS_REFCAT:
                    # NOTE: This .name attrib agreed upon by the JWST Cal
                    #       Working Group.
                    #       Current value is merely a place-holder based
                    #       on HST conventions. This value should also be
                    #       translated to the FITS WCSNAME keyword
                    #       IF that is what gets recorded in the archive
                    #       for end-user searches.
                    imcat.wcs.name = f"FIT-LVL2-{self.abs_refcat}"

                image_model.meta.wcs = imcat.wcs

        return images

    def _is_wcs_correction_small(self, wcs, twcs):
        """Check that the newly tweaked wcs hasn't gone off the rails"""
        tolerance = 10.0 * self.tolerance * u.arcsec

        ra, dec = wcs.footprint(axis_type="spatial").T
        tra, tdec = twcs.footprint(axis_type="spatial").T
        skycoord = SkyCoord(ra=ra, dec=dec, unit="deg")
        tskycoord = SkyCoord(ra=tra, dec=tdec, unit="deg")

        separation = skycoord.separation(tskycoord)

        return (separation < tolerance).all()

    def _imodel2wcsim(self, image_model):
        image_model = image_model if isinstance(
            image_model, datamodels.DataModel
            ) else datamodels.open(path.basename(image_model))
        # make sure that we have a catalog:
        if hasattr(image_model, 'catalog'):
            catalog = image_model.catalog
        else:
            catalog = image_model.meta.tweakreg_catalog

        model_name = path.splitext(image_model.meta.filename)[0].strip('_- ')

        if isinstance(catalog, Table):
            if not catalog.meta.get('name', None):
                catalog.meta['name'] = model_name
        else:
            try:
                cat_name = str(catalog)
                catalog = Table.read(catalog)
                catalog.meta['name'] = cat_name
            except OSError:
                self.log.error(f"Cannot read catalog {catalog}")

        # make sure catalog has 'x' and 'y' columns
        for axis in ['x', 'y']:
            if axis not in catalog.colnames:
                long_axis = axis + 'centroid'
                if long_axis in catalog.colnames:
                    catalog.rename_column(long_axis, axis)
                else:
                    raise ValueError(
                        "'tweakreg' source catalogs must contain either "
                        "columns 'x' and 'y' or 'xcentroid' and "
                        "'ycentroid'."
                    )

        # create WCSImageCatalog object:
        refang = image_model.meta.wcsinfo
        # TODO: create RSTWCSCorrector in tweakwcs
        im = JWSTWCSCorrector(
            wcs = image_model.meta.wcs,
            wcsinfo = {
                'roll_ref': refang['roll_ref'],
                'v2_ref': refang['v2_ref'],
                'v3_ref': refang['v3_ref'],
            },
            meta = {
                'image_model': image_model,
                'catalog': catalog,
                'name': model_name,
            }
        )

        return im


def _common_name(group):
    file_names = []
    for im in group:
        im = im if isinstance(im, datamodels.DataModel) else datamodels.open(im)
        file_names.append(path.splitext(im.meta.filename)[0].strip('_- '))
        im.close()

    cn = path.commonprefix(file_names)
    return cn


def _parse_catfile(catfile):
    if catfile is None or not catfile.strip():
        return None

    catdict = {}

    with open(catfile) as f:
        catfile_dir = path.dirname(catfile)

        for line in f.readlines():
            sline = line.strip()
            if not sline or sline[0] == '#':
                continue

            data_model, *catalog = sline.split()
            catalog = list(map(str.strip, catalog))
            if len(catalog) == 1:
                catdict[data_model] = path.join(catfile_dir, catalog[0])
            elif len(catalog) == 0:
                catdict[data_model] = None
            else:
                raise ValueError("'catfile' can contain at most two columns.")

    return catdict
