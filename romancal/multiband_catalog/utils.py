from astropy.table import Table


def get_photometry_columns(table):
    """
    Get the photometry-related column names in an astropy table.

    Parameters
    ----------
    table : astropy.table.Table
        The table to search for photometry-related columns.

    Returns
    -------
    result: list
        The list of photometry-related column names.
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
        ):
            phot_cols.append(col)

    return phot_cols


def prefix_colnames(table, prefix):
    """
    Update column names in an astropy table in-place by adding a prefix
    to photometry-related columns.

    Parameters
    ----------
    table : `~astropy.table.Table`
        The table to update.

    prefix : str
        The prefix to add to the photometry-related column names. If
        "det_", the columns will be removed.

    Returns
    -------
    result : `~astropy.table.Table`
        The updated table.
    """
    if not isinstance(prefix, str):
        raise ValueError("prefix must be a string")

    phot_cols = get_photometry_columns(table)
    for col in table.colnames:
        if col in phot_cols:
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
    phot_cols = get_photometry_columns(table)
    for col in table.colnames:
        if col in phot_cols:
            table.remove_column(col)

    return table
