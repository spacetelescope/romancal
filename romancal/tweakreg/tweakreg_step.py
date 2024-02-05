"""
Roman pipeline step for image alignment.
"""

import os
from pathlib import Path

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table
from roman_datamodels import datamodels as rdm
from tweakwcs.correctors import JWSTWCSCorrector
from tweakwcs.imalign import align_wcs
from tweakwcs.matchutils import XYXYMatch

from romancal.lib.basic_utils import is_association

# LOCAL
from ..datamodels import ModelContainer
from ..stpipe import RomanStep
from . import astrometric_utils as amutils


def _oxford_or_str_join(str_list):
    nelem = len(str_list)
    if not nelem:
        return "N/A"
    str_list = list(map(repr, str_list))
    if nelem == 1:
        return str_list
    elif nelem == 2:
        return f"{str_list[0]} or {str_list[1]}"
    else:
        return ", ".join(map(repr, str_list[:-1])) + ", or " + repr(str_list[-1])


SINGLE_GROUP_REFCAT = ["GAIADR3", "GAIADR2", "GAIADR1"]
_SINGLE_GROUP_REFCAT_STR = _oxford_or_str_join(SINGLE_GROUP_REFCAT)
DEFAULT_ABS_REFCAT = SINGLE_GROUP_REFCAT[0]
ALIGN_TO_ABS_REFCAT = True

__all__ = ["TweakRegStep"]


class TweakRegStep(RomanStep):
    """
    TweakRegStep: Image alignment based on catalogs of sources detected in
    input images.
    """

    class_alias = "tweakreg"

    spec = f"""
        use_custom_catalogs = boolean(default=False) # Use custom user-provided catalogs?
        catalog_format = string(default='ascii.ecsv') # Catalog output file format
        catfile = string(default='') # Name of the file with a list of custom user-provided catalogs
        catalog_path = string(default='') # Catalog output file path
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
        abs_refcat = string(default='{DEFAULT_ABS_REFCAT}')  # Absolute reference
        # catalog. Options: {_SINGLE_GROUP_REFCAT_STR}
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
    """  # noqa: E501

    reference_file_types = []
    refcat = None

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
                    "contains no user catalogs."
                )
                use_custom_catalogs = False

        try:
            if use_custom_catalogs and catdict:
                images = ModelContainer()
                if isinstance(input, str):
                    asn_dir = os.path.dirname(input)
                    asn_data = images.read_asn(input)
                    for member in asn_data["products"][0]["members"]:
                        filename = member["expname"]
                        member["expname"] = os.path.join(asn_dir, filename)
                        if filename in catdict:
                            member["tweakreg_catalog"] = catdict[filename]
                        elif "tweakreg_catalog" in member:
                            del member["tweakreg_catalog"]

                    images.from_asn(asn_data)
                elif is_association(input):
                    images.from_asn(input)
                else:
                    images = ModelContainer(input)
                    for im in images:
                        filename = im.meta.filename
                        if filename in catdict:
                            self.log.info(
                                f"setting "
                                f"{filename}.source_detection.tweakreg_catalog_name ="
                                f" {repr(catdict[filename])}"
                            )
                            # set catalog name only (no catalog data at this point)
                            im.meta["source_detection"] = {
                                "tweakreg_catalog_name": catdict[filename],
                            }
            else:
                images = ModelContainer(input)
        except TypeError as e:
            e.args = (
                "Input to tweakreg must be a list of DataModels, an "
                "association, or an already open ModelContainer "
                "containing one or more DataModels.",
            ) + e.args[1:]
            raise e

        if len(self.catalog_path) == 0:
            self.catalog_path = os.getcwd()

        self.catalog_path = Path(self.catalog_path).as_posix()
        self.log.info(f"All source catalogs will be saved to: {self.catalog_path}")

        if self.abs_refcat is None or len(self.abs_refcat.strip()) == 0:
            self.abs_refcat = DEFAULT_ABS_REFCAT

        if self.abs_refcat != DEFAULT_ABS_REFCAT:
            # Set expand_refcat to True to eliminate possibility of duplicate
            # entries when aligning to absolute astrometric reference catalog
            self.expand_refcat = True

        if len(images) == 0:
            raise ValueError("Input must contain at least one image model.")

        # Build the catalogs for input images
        for i, image_model in enumerate(images):
            if image_model.meta.exposure.type != "WFI_IMAGE":
                # Check to see if attempt to run tweakreg on non-Image data
                self.log.info("Skipping TweakReg for spectral exposure.")
                # Uncomment below once rad & input data have the cal_step tweakreg
                # image_model.meta.cal_step.tweakreg = "SKIPPED"
                return image_model

            if hasattr(image_model.meta, "source_detection"):
                is_tweakreg_catalog_present = hasattr(
                    image_model.meta.source_detection, "tweakreg_catalog"
                )
                is_tweakreg_catalog_name_present = hasattr(
                    image_model.meta.source_detection, "tweakreg_catalog_name"
                )
                if is_tweakreg_catalog_present:
                    # read catalog from structured array
                    catalog = Table(
                        np.asarray(image_model.meta.source_detection.tweakreg_catalog)
                    )
                elif is_tweakreg_catalog_name_present:
                    catalog = Table.read(
                        image_model.meta.source_detection.tweakreg_catalog_name,
                        format=self.catalog_format,
                    )
                else:
                    raise AttributeError(
                        "Attribute 'meta.source_detection.tweakreg_catalog' is missing."
                        "Please either run SourceDetectionStep or provide a"
                        "custom source catalog."
                    )
                # remove 4D numpy array from meta.source_detection
                if is_tweakreg_catalog_present:
                    del image_model.meta.source_detection["tweakreg_catalog"]
            else:
                raise AttributeError(
                    "Attribute 'meta.source_detection' is missing."
                    "Please either run SourceDetectionStep or provide a"
                    "custom source catalog."
                )

            for axis in ["x", "y"]:
                if axis not in catalog.colnames:
                    long_axis = axis + "centroid"
                    if long_axis in catalog.colnames:
                        catalog.rename_column(long_axis, axis)
                    else:
                        raise ValueError(
                            "'tweakreg' source catalogs must contain a header with "
                            "columns named either 'x' and 'y' or "
                            "'xcentroid' and 'ycentroid'."
                        )

            filename = image_model.meta.filename

            # filter out sources outside the WCS bounding box
            bb = image_model.meta.wcs.bounding_box
            x = catalog["x"]
            y = catalog["y"]
            if bb is None:
                r, d = image_model.meta.wcs(x, y)
                mask = np.isfinite(r) & np.isfinite(d)
                catalog = catalog[mask]

                n_removed_src = np.sum(np.logical_not(mask))
                if n_removed_src:
                    self.log.info(
                        f"Removed {n_removed_src} sources from {filename}'s "
                        "catalog whose image coordinates could not be "
                        "converted to world coordinates."
                    )
            else:
                # assume image coordinates of all sources within a bounding box
                # can be converted to world coordinates.
                ((xmin, xmax), (ymin, ymax)) = bb
                mask = (x > xmin) & (x < xmax) & (y > ymin) & (y < ymax)
                catalog = catalog[mask]

                n_removed_src = np.sum(np.logical_not(mask))
                if n_removed_src:
                    self.log.info(
                        f"Removed {n_removed_src} sources from {filename}'s "
                        "catalog that were outside of the bounding box."
                    )

            # set meta.tweakreg_catalog
            image_model.meta["tweakreg_catalog"] = catalog.as_array()

            nsources = len(catalog)
            if nsources == 0:
                self.log.warning(f"No sources found in {filename}.")
            else:
                self.log.info(f"Detected {len(catalog)} sources in {filename}.")

            images[i] = image_model

        # group images by their "group id":
        grp_img = list(images.models_grouped)

        self.log.info("")
        self.log.info(f"Number of image groups to be aligned: {len(grp_img):d}.")
        self.log.info("Image groups:")

        if len(grp_img) == 1 and not ALIGN_TO_ABS_REFCAT:
            self.log.info("* Images in GROUP 1:")
            for im in grp_img[0]:
                self.log.info(f"     {im.meta.filename}")
            self.log.info("")

            # we need at least two exposures to perform image alignment
            self.log.warning("At least two exposures are required for image alignment.")
            self.log.warning("Nothing to do. Skipping 'TweakRegStep'...")
            self.skip = True
            for model in images:
                model.meta.cal_step["tweakreg"] = "SKIPPED"
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
                model = (
                    model
                    if isinstance(model, rdm.DataModel)
                    else rdm.open(os.path.basename(model))
                )
            self.log.info(f"* Images in GROUP '{group_name}':")
            for im in imcats:
                im.meta["group_id"] = group_name
                self.log.info(f"     {im.meta['name']}")

            self.log.info("")

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
                    # for model in g:
                    #     del model.catalog
                    self.log.info(f"* Images in GROUP '{group_name}':")
                    for im in wcsimlist:
                        im.meta["group_id"] = group_name
                        self.log.info(f"     {im.meta['name']}")
                    imcats.extend(wcsimlist)

            self.log.info("")

            # align images:
            xyxymatch = XYXYMatch(
                searchrad=self.searchrad,
                separation=self.separation,
                use2dhist=self.use2dhist,
                tolerance=self.tolerance,
                xoffset=0,
                yoffset=0,
            )

            try:
                align_wcs(
                    imcats,
                    refcat=None or self.refcat,
                    enforce_user_order=self.enforce_user_order,
                    expand_refcat=self.expand_refcat,
                    minobj=self.minobj,
                    match=xyxymatch,
                    fitgeom=self.fitgeometry,
                    nclip=self.nclip,
                    sigma=(self.sigma, "rmse"),
                )

            except ValueError as e:
                msg = e.args[0]
                if (
                    msg == "Too few input images (or groups of images) with non-empty"
                    " catalogs."
                ):
                    # we need at least two exposures to perform image alignment
                    self.log.warning(msg)
                    self.log.warning(
                        "At least two exposures are required for image alignment."
                    )
                    self.log.warning("Nothing to do. Skipping 'TweakRegStep'...")
                    for model in images:
                        model.meta.cal_step["tweakreg"] = "SKIPPED"
                    if not ALIGN_TO_ABS_REFCAT:
                        self.skip = True
                        return images
                else:
                    raise e

            for imcat in imcats:
                model = imcat.meta["image_model"]
                if model.meta.cal_step.get("tweakreg") == "SKIPPED":
                    continue
                wcs = model.meta.wcs
                twcs = imcat.wcs
                if not self._is_wcs_correction_small(wcs, twcs):
                    # Large corrections are typically a result of source
                    # mis-matching or poorly-conditioned fit. Skip such models.
                    self.log.warning(
                        "WCS has been tweaked by more than"
                        f" {10 * self.tolerance} arcsec"
                    )

                    for model in images:
                        model.meta.cal_step["tweakreg"] = "SKIPPED"
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
                output_name = os.path.join(
                    self.catalog_path, f"fit_{self.abs_refcat.lower()}_ref.ecsv"
                )
            else:
                output_name = None

            # initial shift to be used with absolute astrometry
            self.abs_xoffset = 0
            self.abs_yoffset = 0

            self.abs_refcat = self.abs_refcat.strip()
            gaia_cat_name = self.abs_refcat.upper()

            if gaia_cat_name in SINGLE_GROUP_REFCAT:
                try:
                    ref_cat = amutils.create_astrometric_catalog(
                        images, gaia_cat_name, output=output_name
                    )
                except Exception as e:
                    self.log.warning(
                        "TweakRegStep cannot proceed because of an error that "
                        "occurred while fetching data from the VO server. "
                        f"Returned error message: '{e}'"
                    )
                    self.log.warning("Skipping 'TweakRegStep'...")
                    self.skip = True
                    for model in images:
                        model.meta.cal_step["tweakreg"] = "SKIPPED"
                    return images

            elif os.path.isfile(self.abs_refcat):
                ref_cat = Table.read(self.abs_refcat)

            else:
                raise ValueError(
                    "'abs_refcat' must be a path to an "
                    "existing file name or one of the supported "
                    f"reference catalogs: {_SINGLE_GROUP_REFCAT_STR}."
                )

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
                    yoffset=self.abs_yoffset,
                )

                # Set group_id to same value so all get fit as one observation
                # The assigned value, 987654, has been hard-coded to make it
                # easy to recognize when alignment to GAIA was being performed
                # as opposed to the group_id values used for relative alignment
                # earlier in this step.
                for imcat in imcats:
                    imcat.meta["group_id"] = 987654
                    if (
                        "fit_info" in imcat.meta
                        and "REFERENCE" in imcat.meta["fit_info"]["status"]
                    ):
                        del imcat.meta["fit_info"]

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
                    sigma=(self.abs_sigma, "rmse"),
                    ref_tpwcs=imcats[0],
                )

        for imcat in imcats:
            image_model = imcat.meta["image_model"]
            image_model.meta.cal_step["tweakreg"] = "COMPLETE"

            # retrieve fit status and update wcs if fit is successful:
            if "SUCCESS" in imcat.meta.get("fit_info")["status"]:
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

                # serialize object from tweakwcs
                # (typecasting numpy objects to python types so that it doesn't cause an
                # issue when saving datamodel to ASDF)
                wcs_fit_results = {
                    k: v.tolist() if isinstance(v, (np.ndarray, np.bool_)) else v
                    for k, v in imcat.meta["fit_info"].items()
                }
                # add fit results and new WCS to datamodel
                image_model.meta["wcs_fit_results"] = wcs_fit_results
                # remove unwanted keys from WCS fit results
                for k in [
                    "eff_minobj",
                    "matched_ref_idx",
                    "matched_input_idx",
                    "fit_RA",
                    "fit_DEC",
                    "fitmask",
                ]:
                    del image_model.meta["wcs_fit_results"][k]

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
        image_model = (
            image_model
            if isinstance(image_model, rdm.DataModel)
            else rdm.open(os.path.basename(image_model))
        )
        catalog = image_model.meta.tweakreg_catalog
        model_name = os.path.splitext(image_model.meta.filename)[0].strip("_- ")

        try:
            if self.use_custom_catalogs:
                catalog_format = self.catalog_format
            else:
                catalog_format = "ascii.ecsv"

            if isinstance(catalog, str):
                # a string with the name of the catalog was provided
                catalog = Table.read(catalog, format=catalog_format)
            else:
                # catalog is a structured array, convert to astropy table:
                catalog = Table(catalog)

            catalog.meta["name"] = (
                str(catalog) if isinstance(catalog, str) else model_name
            )
        except OSError:
            self.log.error(f"Cannot read catalog {catalog}")

        # make sure catalog has 'x' and 'y' columns
        for axis in ["x", "y"]:
            if axis not in catalog.colnames:
                long_axis = axis + "centroid"
                if long_axis in catalog.colnames:
                    catalog.rename_column(long_axis, axis)
                else:
                    raise ValueError(
                        "'tweakreg' source catalogs must contain either columns 'x' and"
                        " 'y' or 'xcentroid' and 'ycentroid'."
                    )

        # create WCSImageCatalog object:
        refang = image_model.meta.wcsinfo
        # TODO: create RSTWCSCorrector in tweakwcs
        im = JWSTWCSCorrector(
            wcs=image_model.meta.wcs,
            wcsinfo={
                "roll_ref": refang["roll_ref"],
                "v2_ref": refang["v2_ref"],
                "v3_ref": refang["v3_ref"],
            },
            meta={
                "image_model": image_model,
                "catalog": catalog,
                "name": model_name,
            },
        )

        return im


def _common_name(group):
    file_names = []
    for im in group:
        if isinstance(im, rdm.DataModel):
            file_names.append(os.path.splitext(im.meta.filename)[0].strip("_- "))
        else:
            raise TypeError("Input must be a list of datamodels list.")

    cn = os.path.commonprefix(file_names)
    return cn


def _parse_catfile(catfile):
    if catfile is None or not catfile.strip():
        return None

    catdict = {}

    with open(catfile) as f:
        catfile_dir = os.path.dirname(catfile)

        for line in f.readlines():
            sline = line.strip()
            if not sline or sline[0] == "#":
                continue

            data_model, *catalog = sline.split()
            catalog = list(map(str.strip, catalog))
            if len(catalog) == 1:
                catdict[data_model] = os.path.join(catfile_dir, catalog[0])
            elif len(catalog) == 0:
                catdict[data_model] = None
            else:
                raise ValueError("'catfile' can contain at most two columns.")

    return catdict
