from astropy.table import Table


def get_direct_image_columns(table):
    """
    Get the column names associated with measurements on a direct image.

    Parameters
    ----------
    table : astropy.table.Table
        The table to search for columns computed from the direct image.

    Returns
    -------
    result: list
        The list of column names computed from the direct image.
    """
    if not isinstance(table, Table):
        raise ValueError("table must be an astropy.table.Table object")

    phot_cols = []
    for col in table.colnames:
        if (
            col.startswith("aper")
            or col.startswith("segment_flux")
            or col.startswith("kron_flux")
            or col == "sharpness"
            or col == "roundness"
        ):
            phot_cols.append(col)

    return phot_cols


def get_detection_image_columns(table):
    """
    Get the column names associated with measurements on a detection image.

    Parameters
    ----------
    table : astropy.table.Table
        The table to search for columns computed from the detection image.

    Returns
    -------
    result: list
        The list of column names computed from the detection image
    """
    if not isinstance(table, Table):
        raise ValueError("table must be an astropy.table.Table object")

    phot_cols = []
    for col in table.colnames:
        if (
            col.endswith("centroid")
            or col.startswith("nn_")
            or col in ["semimajor_sigma", "semiminor_sigma", "ellipticity"]
            or col.endswith("orientation")
            or col == "segment_area"
        ):
            phot_cols.append(col)

    return phot_cols


def prefix_colnames(table, prefix, colnames=None):
    """
    Update column names in an astropy table in-place by adding a prefix
    to specified colnames, defaulting to direct-image related columns.

    Parameters
    ----------
    table : `~astropy.table.Table`
        The table to update.

    prefix : str
        The prefix to add to the photometry-related column names. If
        "det_", the columns will be removed.

    colnames : str
        The column names to prefix

    Returns
    -------
    result : `~astropy.table.Table`
        The updated table.
    """
    if not isinstance(prefix, str):
        raise ValueError("prefix must be a string")

    if colnames is None:
        colnames = get_direct_image_columns(table)
    for col in table.colnames:
        if col in colnames:
            table.rename_column(col, prefix + col)

    return table


def remove_columns(table):
    """
    Remove photometry-related columns from an astropy table in-place.

    Parameters
    ----------
    table : `~astropy.table.Table`
        The table to update.

    Returns
    -------
    result : `~astropy.table.Table`
        The updated table.
    """
    phot_cols = get_direct_image_columns(table)
    for col in table.colnames:
        if col in phot_cols:
            table.remove_column(col)

    return table


def insert_filter(colname, filter_name):
    """
    Insert the filter name into the column name.

    Parameters
    ----------
    colname : str
        The original column name.

    filter_name : str
        The filter name to insert.

    Returns
    -------
    result : str
        The updated column name.
    """
    if colname.endswith("_err"):
        base = colname[:-4]
        return f"{base}_{filter_name}_err"
    else:
        return f"{colname}_{filter_name}"


def add_filter_to_colnames(table, filter_name):
    """
    Add a filter name to the column names in an astropy table.

    The filter name is inserted before the "_flux" part of the column
    name.

    Parameters
    ----------
    table : `~astropy.table.Table`
        The table to update.

    filter_name : str
        The filter name to add to the column names.

    Returns
    -------
    result : `~astropy.table.Table`
        The updated table.
    """
    if not isinstance(filter_name, str):
        raise ValueError("filter_name must be a string")

    filter_name = filter_name.lower()
    append_cols = ("is_extended", "sharpness", "roundness", "psf_flags")

    for colname in table.colnames:
        if "_flux" in colname or "_psf" in colname:
            new_colname = insert_filter(colname, filter_name)
            table.rename_column(colname, new_colname)
        elif colname in append_cols:
            table.rename_column(colname, f"{colname}_{filter_name}")

    return table
