"""
Roman pipeline step for image alignment.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
from astropy.table import Table
from roman_datamodels import datamodels as rdm
from stcal.tweakreg import tweakreg
from stcal.tweakreg.tweakreg import _SINGLE_GROUP_REFCAT_STR, SINGLE_GROUP_REFCAT

from romancal.assign_wcs.utils import add_s_region
from romancal.lib.save_wcs import save_wfiwcs

# LOCAL
from ..datamodels import ModelLibrary
from ..stpipe import RomanStep

if TYPE_CHECKING:
    from typing import ClassVar

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
        update_source_catalog_coordinates = boolean(default=False) # Update source catalog file with tweaked coordinates?
        save_l1_wcs = boolean(default=True)
        vo_timeout = integer(min=0, default=120) # VO catalog service timeout.
    """

    reference_file_types: ClassVar = []

    def process(self, input):
        # properly handle input
        try:
            if isinstance(input, rdm.DataModel):
                images = ModelLibrary([input])
            elif str(input).endswith(".asdf"):
                images = ModelLibrary([rdm.open(input)])
            elif isinstance(input, ModelLibrary):
                images = input
            else:
                images = ModelLibrary(input)
        except TypeError as e:
            e.args = (
                "Input to tweakreg must be a list of DataModels, an "
                "association, or an already open ModelLibrary "
                "containing one or more DataModels.",
                *e.args[1:],
            )
            raise e

        if not images:
            raise ValueError("Input must contain at least one image model.")

        self.log.info(
            f"Number of image groups to be aligned: {len(images.group_indices):d}."
        )
        self.log.info("Image groups:")
        for name in images.group_names:
            self.log.info(f"  {name}")
        # set the first image as reference
        with images:
            ref_image = images.borrow(0)
            images.shelve(ref_image, 0, modify=False)

        catdict = _parse_catfile(self.catfile)

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
                        model.meta["source_catalog"] = {
                            "tweakreg_catalog_name": catdict[filename],
                        }
                        images.shelve(model, i)
                    else:
                        images.shelve(model, i, modify=False)

        # set path where the source catalog will be saved to
        if len(self.catalog_path) == 0:
            self.catalog_path = os.getcwd()
        self.catalog_path = Path(self.catalog_path).as_posix()
        self.log.info(f"All source catalogs will be saved to: {self.catalog_path}")

        # set reference catalog name
        if not self.abs_refcat:
            self.abs_refcat = DEFAULT_ABS_REFCAT.strip().upper()
        if self.abs_refcat != DEFAULT_ABS_REFCAT:
            self.expand_refcat = True

        # build the catalogs for input images
        imcats = []
        with images:
            for i, image_model in enumerate(images):
                exposure_type = image_model.meta.exposure.type
                if exposure_type != "WFI_IMAGE":
                    self.log.info("Skipping TweakReg for spectral exposure.")
                    image_model.meta.cal_step.tweakreg = "SKIPPED"
                else:
                    source_catalog = getattr(image_model.meta, "source_catalog", None)
                    if source_catalog is None:
                        images.shelve(image_model, i, modify=False)
                        raise AttributeError(
                            "Attribute 'meta.source_catalog' is missing. "
                            "Please either run SourceCatalogStep or provide a custom source catalog."
                        )

                    try:
                        catalog = self.get_tweakreg_catalog(source_catalog, image_model)
                    except AttributeError as e:
                        self.log.error(f"Failed to retrieve tweakreg_catalog: {e}")
                        images.shelve(image_model, i, modify=False)
                        raise e

                    try:
                        # validate catalog columns
                        _validate_catalog_columns(catalog)
                    except ValueError as e:
                        self.log.error(f"Failed to validate catalog columns: {e}")
                        images.shelve(image_model, i, modify=False)
                        raise e

                    catalog = tweakreg.filter_catalog_by_bounding_box(
                        catalog, image_model.meta.wcs.bounding_box
                    )

                    if self.save_abs_catalog:
                        output_name = os.path.join(
                            self.catalog_path, f"fit_{self.abs_refcat.lower()}_ref.ecsv"
                        )
                        catalog.write(
                            output_name, format=self.catalog_format, overwrite=True
                        )

                    image_model.meta["tweakreg_catalog"] = catalog.as_array()
                    nsources = len(catalog)
                    self.log.info(
                        f"Detected {nsources} sources in {image_model.meta.filename}."
                        if nsources
                        else f"No sources found in {image_model.meta.filename}."
                    )
                    # build image catalog
                    # catalog name
                    catalog_name = os.path.splitext(image_model.meta.filename)[0].strip(
                        "_- "
                    )
                    # catalog data
                    catalog_table = Table(image_model.meta.tweakreg_catalog)
                    catalog_table.meta["name"] = catalog_name

                    imcat = tweakreg.construct_wcs_corrector(
                        wcs=image_model.meta.wcs,
                        refang=image_model.meta.wcsinfo,
                        catalog=catalog_table,
                        group_id=images._model_to_group_id(image_model),
                    )
                    imcat.meta["model_index"] = i
                    imcats.append(imcat)
                images.shelve(image_model, i)

        # run alignment only if it was possible to build image catalogs
        if len(imcats):
            # extract WCS correctors to use for image alignment
            if len(images.group_indices) > 1:
                self.do_relative_alignment(imcats)

            if self.abs_refcat in SINGLE_GROUP_REFCAT:
                self.do_absolute_alignment(ref_image, imcats)

            # finalize step
            with images:
                for imcat in imcats:
                    image_model = images.borrow(imcat.meta["model_index"])
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
                            k: (
                                v.tolist()
                                if isinstance(v, np.ndarray | np.bool_)
                                else v
                            )
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

                        # update WCS
                        image_model.meta.wcs = imcat.wcs
                        # update S_REGION
                        add_s_region(image_model)

                    images.shelve(image_model, imcat.meta["model_index"])

        # Write out the WfiWcs products
        if self.save_l1_wcs:
            save_wfiwcs(self, images, force=True)

        return images

    def update_catalog_coordinates(self, tweakreg_catalog_name, tweaked_wcs):
        """
        Update the source catalog coordinates using the tweaked WCS.

        Parameters
        ----------
        tweakreg_catalog_name : str
            The name of the TweakReg catalog file produced by `SourceCatalog`.
        tweaked_wcs : `gwcs.wcs.WCS`
            The tweaked World Coordinate System (WCS) object.

        Returns
        -------
        None
        """
        # read in cat file
        with rdm.open(tweakreg_catalog_name) as source_catalog_model:
            # get catalog
            catalog = source_catalog_model.source_catalog

            # define mapping between pixel and world coordinates
            colname_mapping = {
                ("xcentroid", "ycentroid"): ("ra_centroid", "dec_centroid"),
                ("x_psf", "y_psf"): ("ra_psf", "dec_psf"),
            }

            for k, v in colname_mapping.items():
                # get column names
                x_colname, y_colname = k
                ra_colname, dec_colname = v

                # calculate new coordinates using tweaked WCS and update catalog coordinates
                catalog[ra_colname], catalog[dec_colname] = tweaked_wcs(
                    catalog[x_colname], catalog[y_colname]
                )

            # save updated catalog (overwrite cat file)
            self.save_model(
                source_catalog_model,
                output_file=source_catalog_model.meta.filename,
                suffix="cat",
                force=True,
            )

    def read_catalog(self, catalog_name):
        """
        Reads a source catalog from a specified file.

        This function determines the format of the catalog based on the
        file extension:

        * "asdf":  uses roman datamodels
        * "parquet":  uses pyarrow
        * otherwise:  uses astropy Table.

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
        filetype = (
            "parquet" if catalog_name.endswith("parquet") else self.catalog_format
        )
        if catalog_name.endswith("asdf"):
            # leave this for now
            with rdm.open(catalog_name) as source_catalog_model:
                catalog = source_catalog_model.source_catalog
        else:
            catalog = Table.read(catalog_name, format=filetype)
        return catalog

    def get_tweakreg_catalog(self, source_catalog, image_model):
        """
        Retrieve the tweakreg catalog from source detection.

        This method checks the source detection metadata for the presence of a
        tweakreg catalog data or a string with its name. It returns the catalog
        as a Table object if either is found, or raises an error if neither is available.

        Parameters
        ----------
        source_catalog : object
            The source catalog metadata containing catalog information.
        image_model : DataModel
            The image model associated with the source detection.

        Returns
        -------
        Table
            The retrieved tweakreg catalog as a Table object.

        Raises
        ------
        AttributeError
            If the required catalog information is missing from the source detection.
        """
        if getattr(source_catalog, "tweakreg_catalog", None):
            tweakreg_catalog = Table(np.asarray(source_catalog.tweakreg_catalog))
            del image_model.meta.source_catalog["tweakreg_catalog"]
            return tweakreg_catalog

        if getattr(source_catalog, "tweakreg_catalog_name", None):
            return self.read_catalog(source_catalog.tweakreg_catalog_name)

        raise AttributeError(
            "Attribute 'meta.source_catalog.tweakreg_catalog' is missing. "
            "Please either run SourceCatalogStep or provide a custom source catalog."
        )

    def do_relative_alignment(self, imcats):
        """
        Perform relative alignment of images.

        This method performs relative alignment with the specified parameters,
        including search radius, separation, and fitting geometry.

        Parameters
        ----------
        imcats : list
            A list of image catalogs containing source information for alignment.

        Returns
        -------
        None
        """
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
            clip_accum=True,
        )

    def do_absolute_alignment(self, ref_image, imcats):
        """
        Perform absolute alignment of images.

        This method retrieves a reference image and performs absolute alignment
        using the specified parameters, including reference WCS information and
        catalog details. It aligns the provided image catalogs to the absolute
        reference catalog.

        Parameters
        ----------
        ref_image : DataModel
            The reference image used for alignment, which contains WCS information.
        imcats : list
            A list of image catalogs containing source information for alignment.

        Returns
        -------
        None
        """
        tweakreg.absolute_align(
            imcats,
            self.abs_refcat,
            ref_wcs=ref_image.meta.wcs,
            ref_wcsinfo=ref_image.meta.wcsinfo,
            epoch=ref_image.meta.exposure.start_time.decimalyear,
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
            clip_accum=True,
            timeout=self.vo_timeout,
        )


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


def _validate_catalog_columns(catalog):
    """
    Validate the presence of required columns in the catalog.

    This method checks if the specified axis column exists in the catalog.
    If the axis is not found, it looks for a corresponding psf column
    and renames it if present. If neither is found, it raises an error.

    Parameters
    ----------
    catalog : Table
        The catalog to validate, which should contain source information.
    axis : str
        The axis to check for in the catalog (e.g., 'x' or 'y').

    Returns
    -------
    None

    Raises
    ------
    ValueError
        If the required columns are missing from the catalog.
    """
    for axis in ["x", "y"]:
        if axis not in catalog.colnames:
            long_axis = f"{axis}_psf"
            if long_axis in catalog.colnames:
                catalog.rename_column(long_axis, axis)
            else:
                raise ValueError(
                    "'tweakreg' source catalogs must contain a header with "
                    "columns named either 'x' and 'y' or 'x_psf' and 'y_psf'."
                )
    return catalog
