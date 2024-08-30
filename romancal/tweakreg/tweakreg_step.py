"""
Roman pipeline step for image alignment.
"""

import os
from pathlib import Path

import numpy as np
from astropy.table import Table
from roman_datamodels import datamodels as rdm
from stcal.tweakreg import tweakreg
from typing import List

# LOCAL
from ..datamodels import ModelLibrary
from ..stpipe import RomanStep


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
        abs_separation = float(default=1.0) # Minimum object separation in arcsec when performing absolute astrometry
        abs_tolerance = float(default=0.7) # Matching tolerance for xyxymatch in arcsec when performing absolute astrometry
        # Fitting geometry when performing absolute astrometry
        abs_fitgeometry = option('shift', 'rshift', 'rscale', 'general', default='rshift')
        abs_nclip = integer(min=0, default=3) # Number of clipping iterations in fit when performing absolute astrometry
        abs_sigma = float(min=0.0, default=3.0) # Clipping limit in sigma units when performing absolute astrometry
        output_use_model = boolean(default=True)  # When saving use `DataModel.meta.filename`
    """  # noqa: E501

    reference_file_types = []

    def process(self, input):

        # properly handle input
        images = self.parse_input(step_input=input)

        catdict = _parse_catfile(self.catfile)

        if self.use_custom_catalogs:
            self.validate_custom_catalogs(catdict, images)

        # set path where the source catalog will be saved to
        self.set_catalog_path()

        self.set_reference_catalog()

        # Build the catalogs for input images
        self.set_tweakreg_catalog_attribute(images)

        imcats = _build_image_catalogs(images)

        self.do_relative_alignment(images, imcats)

        self.do_absolute_alignment(images, imcats)

        self.finalize_step(images, imcats)

        return images

    def parse_input(self, step_input) -> ModelLibrary:
        """
        Parse the input for the class.

        This method handles various types of input, including single DataModels,
        file paths to ASDF files, and collections of DataModels. It ensures that
        the input is converted into a ModelLibrary for further processing.

        Parameters
        ----------
        step_input : Union[rdm.DataModel, str, ModelLibrary, list]
            The input to be parsed, which can be a DataModel, a file path,
            a ModelLibrary, or a list of DataModels.

        Returns
        -------
        ModelLibrary
            A ModelLibrary containing the parsed images.

        Raises
        ------
        TypeError
            If the input is not a valid type for processing.
        """
        try:
            if isinstance(step_input, rdm.DataModel):
                images = ModelLibrary([step_input])
            elif str(step_input).endswith(".asdf"):
                images = ModelLibrary([rdm.open(step_input)])
            elif isinstance(step_input, ModelLibrary):
                images = step_input
            else:
                images = ModelLibrary(step_input)
        except TypeError as e:
            e.args = (
                "Input to tweakreg must be a list of DataModels, an "
                "association, or an already open ModelLibrary "
                "containing one or more DataModels.",
            ) + e.args[1:]
            raise e

        if len(images) == 0:
            raise ValueError("Input must contain at least one image model.")

        self.log.info("")
        self.log.info(
            f"Number of image groups to be aligned: {len(images.group_indices):d}."
        )
        self.log.info("Image groups:")

        return images

    def validate_custom_catalogs(self, catdict, images):
        """
        Validate and apply custom catalogs for the tweak registration step.

        This method checks if the user has requested the use of custom catalogs
        and whether the provided catalog file contains valid entries. If valid
        catalogs are found, it updates the image models with the corresponding
        catalog names.

        Parameters
        ----------
        catdict : dict
            A dictionary mapping image filenames to custom catalog file paths.
        images : ModelLibrary
            A collection of image models to be updated with custom catalog information.

        Returns
        -------
        None
        """
        use_custom_catalogs = self.use_custom_catalogs
        # if user requested the use of custom catalogs and provided a
        # valid 'catfile' file name that has no custom catalogs,
        # turn off the use of custom catalogs:
        if catdict is not None and not catdict:
            self.log.warning(
                "'use_custom_catalogs' is set to True but 'catfile' "
                "contains no user catalogs."
            )
            use_custom_catalogs = False

        if use_custom_catalogs and catdict:
            with images:
                for i, member in enumerate(images.asn["products"][0]["members"]):
                    filename = member["expname"]
                    if filename in catdict:
                        # FIXME: I'm not sure if this captures all the possible combinations
                        # for example, meta.tweakreg_catalog is set by the container (when
                        # it's present in the association). However the code in this step
                        # checks meta.source_catalog.tweakreg_catalog. I think this means
                        # that setting a catalog via an association does not work. Is this
                        # intended? If so, the container can be updated to not support that.
                        model = images.borrow(i)
                        model.meta["source_detection"] = {
                            "tweakreg_catalog_name": catdict[filename],
                        }
                        images.shelve(model, i)
                    else:
                        images.shelve(model, i, modify=False)

    def set_catalog_path(self):
        if len(self.catalog_path) == 0:
            self.catalog_path = os.getcwd()
        self.catalog_path = Path(self.catalog_path).as_posix()
        self.log.info(f"All source catalogs will be saved to: {self.catalog_path}")

    def set_reference_catalog(self):
        if self.abs_refcat is None or len(self.abs_refcat) == 0:
            self.abs_refcat = DEFAULT_ABS_REFCAT.strip().upper()
        if self.abs_refcat != DEFAULT_ABS_REFCAT:
            # Set expand_refcat to True to eliminate possibility of duplicate
            # entries when aligning to absolute astrometric reference catalog
            self.expand_refcat = True

    def set_tweakreg_catalog_attribute(self, images):
        """
        Set the tweakreg catalog attribute for each image model.

        This method iterates through the provided image models, checking the
        exposure type and ensuring that the necessary source detection metadata
        is present. It retrieves the tweak registration catalog, validates its
        columns, filters it based on WCS, and updates the image model's metadata.

        Parameters
        ----------
        images : ModelLibrary
            A collection of image models to be updated with tweak registration catalogs.

        Returns
        -------
        None

        Raises
        ------
        AttributeError
            If the required source detection metadata is missing from an image model.

        Logs
        -----
        Information about the number of detected sources is logged for each image model.
        """

        with images:
            for i, image_model in enumerate(images):
                exposure_type = image_model.meta.exposure.type
                if exposure_type != "WFI_IMAGE":
                    self.log.info("Skipping TweakReg for spectral exposure.")
                    image_model.meta.cal_step.tweakreg = "SKIPPED"
                    images.shelve(image_model)
                    return image_model

                source_detection = getattr(image_model.meta, "source_detection", None)
                if source_detection is None:
                    images.shelve(image_model, i, modify=False)
                    raise AttributeError(
                        "Attribute 'meta.source_detection' is missing. "
                        "Please either run SourceDetectionStep or provide a custom source catalog."
                    )

                catalog = self.get_tweakreg_catalog(source_detection, image_model, i)

                for axis in ["x", "y"]:
                    _validate_catalog_columns(catalog, axis, image_model, i)

                filename = image_model.meta.filename
                catalog = self.filter_catalog_by_wcs(catalog, image_model, filename)

                if self.save_abs_catalog:
                    self.save_abs_ref_catalog(catalog)

                image_model.meta["tweakreg_catalog"] = catalog.as_array()
                nsources = len(catalog)
                self.log.info(
                    f"Detected {nsources} sources in {filename}."
                    if nsources
                    else f"No sources found in {filename}."
                )
                images.shelve(image_model, i)

    def do_relative_alignment(self, images, imcats):
        if len(images.group_indices) > 1:
            # local align images:
            tweakreg.relative_align(
                imcats,
                searchrad=self.searchrad,
                separation=self.separation,
                use2dhist=self.use2dhist,
                tolerance=self.tolerance,
                xoffset=0,
                yoffset=0,
                enforce_user_order=self.enforce_user_order,
                expand_refcat=self.expand_refcat,
                minobj=self.minobj,
                fitgeometry=self.fitgeometry,
                nclip=self.nclip,
                sigma=self.sigma,
            )

    def do_absolute_alignment(self, images, imcats):
        if self.abs_refcat in SINGLE_GROUP_REFCAT:
            ref_image = _get_reference_image(images)

            tweakreg.absolute_align(
                imcats,
                self.abs_refcat,
                ref_wcs=ref_image.meta.wcs,
                ref_wcsinfo=ref_image.meta.wcsinfo,
                epoch=ref_image.meta.exposure.mid_time.decimalyear,
                abs_minobj=self.abs_minobj,
                abs_fitgeometry=self.abs_fitgeometry,
                abs_nclip=self.abs_nclip,
                abs_sigma=self.abs_sigma,
                abs_searchrad=self.abs_searchrad,
                abs_use2dhist=self.abs_use2dhist,
                abs_separation=self.abs_separation,
                abs_tolerance=self.abs_tolerance,
                save_abs_catalog=self.save_abs_catalog,
                abs_catalog_output_dir=self.output_dir,
            )

    def finalize_step(self, images, imcats):
        """
        Finalize the tweak registration step by updating image metadata and WCS information.

        This method iterates through the provided image catalogs, marking TweakRegStep as complete,
        removing the source catalog, and updating the WCS if the fit was successful.
        It also serializes fit results for storage in the image model's metadata.

        Parameters
        ----------
        images : ModelLibrary
            A collection of image models to be updated.
        imcats : list
            A collection of image catalogs containing fit information.

        Returns
        -------
        None
        """
        with images:
            for i, imcat in enumerate(imcats):
                image_model = images.borrow(i)
                image_model.meta.cal_step["tweakreg"] = "COMPLETE"
                # remove source catalog
                del image_model.meta["tweakreg_catalog"]

                # retrieve fit status and update wcs if fit is successful:
                if "SUCCESS" in imcat.meta.get("fit_info")["status"]:
                    # Update/create the WCS .name attribute with information
                    # on this astrometric fit as the only record that it was
                    # successful:

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
                images.shelve(image_model, i)

    def read_catalog(self, catalog_name):
        """
        Reads a source catalog from a specified file.

        This function determines the format of the catalog based on the file extension.
        If the file ends with "asdf", it uses a specific method to open and read the catalog;
        otherwise, it reads the catalog using a standard table format.

        Parameters
        ----------
        catalog_name : str
            The name of the catalog file to read.

        Returns
        -------
        Table
            The read catalog as a Table object.

        Raises
        ------
        ValueError
            If the catalog format is unsupported.
        """
        if catalog_name.endswith("asdf"):
            with rdm.open(catalog_name) as source_catalog_model:
                catalog = source_catalog_model.source_catalog
        else:
            catalog = Table.read(catalog_name, format=self.catalog_format)
        return catalog

    def save_abs_ref_catalog(self, catalog_table: Table):
        """
        Save the absolute reference catalog to a specified file.

        This method writes the provided catalog table to a file in the specified
        format and location, using a naming convention based on the absolute
        reference catalog.

        Parameters
        ----------
        catalog_table : Table
            The catalog table to be saved as an output file.

        Returns
        -------
        None
        """
        output_name = os.path.join(
            self.catalog_path, f"fit_{self.abs_refcat.lower()}_ref.ecsv"
        )
        catalog_table.write(output_name, format=self.catalog_format, overwrite=True)

    def get_tweakreg_catalog(self, source_detection, image_model, index):
        """
        Retrieve the tweakreg catalog from source detection.

        This method checks the source detection metadata for the presence of a
        tweakreg catalog data or a string with its name. It returns the catalog
        as a Table object if either is found, or raises an error if neither is available.

        Parameters
        ----------
        source_detection : object
            The source detection metadata containing catalog information.
        image_model : DataModel
            The image model associated with the source detection.
        index : int
            The index of the image model in the collection.

        Returns
        -------
        Table
            The retrieved tweakreg catalog as a Table object.

        Raises
        ------
        AttributeError
            If the required catalog information is missing from the source detection.
        """
        if hasattr(source_detection, "tweakreg_catalog"):
            tweakreg_catalog = Table(np.asarray(source_detection.tweakreg_catalog))
            del image_model.meta.source_detection["tweakreg_catalog"]
            return tweakreg_catalog
        elif hasattr(source_detection, "tweakreg_catalog_name"):
            return self.read_catalog(source_detection.tweakreg_catalog_name)
        else:
            images.shelve(image_model, index, modify=False)
            raise AttributeError(
                "Attribute 'meta.source_detection.tweakreg_catalog' is missing. "
                "Please either run SourceDetectionStep or provide a custom source catalog."
            )

    def filter_catalog_by_wcs(self, catalog, image_model, filename):
        """
        Filter sources in the catalog based on WCS bounding box.

        This method removes sources from the catalog that fall outside the
        specified WCS bounding box. If no bounding box is defined, it checks
        the validity of the sources' world coordinates and filters accordingly.

        Parameters
        ----------
        catalog : Table
            The catalog containing source information to be filtered.
        image_model : DataModel
            The image model associated with the catalog, used to access WCS information.
        filename : str
            The name of the file associated with the catalog, used for logging.

        Returns
        -------
        Table
            The filtered catalog containing only sources within the bounding box.

        Logs
        -----
        Information about the number of sources removed from the catalog is logged.
        """
        bb = image_model.meta.wcs.bounding_box
        x, y = catalog["x"], catalog["y"]

        if bb is None:
            r, d = image_model.meta.wcs(x, y)
            mask = np.isfinite(r) & np.isfinite(d)
        else:
            ((xmin, xmax), (ymin, ymax)) = bb
            mask = (x > xmin) & (x < xmax) & (y > ymin) & (y < ymax)

        catalog = catalog[mask]
        n_removed_src = np.sum(np.logical_not(mask))
        if n_removed_src:
            self.log.info(
                f"Removed {n_removed_src} sources from {filename}'s "
                "catalog that were outside of the bounding box."
                if bb
                else f"catalog whose image coordinates could not be converted to world coordinates."
            )

        return catalog


def _parse_catfile(catfile):
    """
    Parse a catalog file and return a dictionary mapping data models to catalog paths.

    This function reads a specified catalog file, extracting data model names and
    their associated catalog paths. It supports a format where each line contains
    a data model followed by an optional catalog path, and it ensures that the
    file adheres to the expected structure.

    Parameters
    ----------
    catfile : str
        The path to the catalog file to be parsed.

    Returns
    -------
    dict or None
        A dictionary mapping data model names to catalog paths, or None if the
        input file is empty or invalid.

    Raises
    ------
    ValueError
        If the catalog file contains more than two columns per line.
    """

    if catfile is None or not catfile.strip():
        return None

    catdict = {}

    with open(catfile) as f:
        catfile_dir = os.path.dirname(catfile)

        for line in f:
            sline = line.strip()
            if not sline or sline[0] == "#":
                continue

            data_model, *catalog = sline.split()
            catalog = list(map(str.strip, catalog))
            if len(catalog) == 1:
                catdict[data_model] = os.path.join(catfile_dir, catalog[0])
            elif not catalog:
                catdict[data_model] = None
            else:
                raise ValueError("'catfile' can contain at most two columns.")

    return catdict


def _build_image_catalogs(images) -> List:
    """
    Build image catalogs from the provided images.

    This method constructs a list of image catalogs by extracting the necessary
    metadata from each image model. It creates a WCS corrector for each image
    based on its associated catalog and metadata.

    Parameters
    ----------
    images : ModelLibrary
        A collection of image models from which to build catalogs.

    Returns
    -------
    imcats : list
        A list of image catalogs constructed from the input images.
    """
    imcats = []
    with images:
        for i, m in enumerate(images):
            # catalog name
            catalog_name = os.path.splitext(m.meta.filename)[0].strip("_- ")
            # catalog data
            catalog_table = Table(m.meta.tweakreg_catalog)
            catalog_table.meta["name"] = catalog_name

            imcats.append(
                tweakreg.construct_wcs_corrector(
                    wcs=m.meta.wcs,
                    refang=m.meta.wcsinfo,
                    catalog=catalog_table,
                    group_id=m.meta.group_id,
                )
            )
            images.shelve(m, i, modify=False)
    return imcats


def _get_reference_image(images):
    with images:
        ref_image = images.borrow(0)
        images.shelve(ref_image, 0, modify=False)
    return ref_image


def _validate_catalog_columns(catalog, axis, image_model, index):
    """
    Validate the presence of required columns in the catalog.

    This method checks if the specified axis column exists in the catalog.
    If the axis is not found, it looks for a corresponding centroid column
    and renames it if present. If neither is found, it raises an error.

    Parameters
    ----------
    catalog : Table
        The catalog to validate, which should contain source information.
    axis : str
        The axis to check for in the catalog (e.g., 'x' or 'y').
    image_model : DataModel
        The image model associated with the catalog.
    index : int
        The index of the image model in the collection.

    Returns
    -------
    None

    Raises
    ------
    ValueError
        If the required columns are missing from the catalog.
    """
    if axis not in catalog.colnames:
        long_axis = f"{axis}centroid"
        if long_axis in catalog.colnames:
            catalog.rename_column(long_axis, axis)
        else:
            images.shelve(image_model, index, modify=False)
            raise ValueError(
                "'tweakreg' source catalogs must contain a header with "
                "columns named either 'x' and 'y' or 'xcentroid' and 'ycentroid'."
            )
