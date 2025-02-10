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
            or col.startswith("CI_")
            or col.startswith("isophotal_flux")
            or col.startswith("kron_flux")
            or "_psf" in col
            or col == "is_extended"
            or col == "sharpness"
            or col == "roundness"
            or col == "flags"
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
            or col in ['semimajor_sigma', 'semiminor_sigma', 'ellipticity']
            or col == 'ellipticity'
            or col.endswith('orientation')
            or col == 'isophotal_area'
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
