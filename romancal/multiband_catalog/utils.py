from astropy.table import Table


def update_colnames(table, prefix):
    """
    Update column names in an astropy table by adding a prefix to
    photometry-related columns.

    If the prefix is "det_", the function will remove these columns.

    Parameters
    ----------
    table : astropy.table.Table
        The table to update.

    prefix : str
        The prefix to add to the photometry-related column names. If
        "det_", the columns will be removed.

    Returns
    -------
    astropy.table.Table
        The updated table.
    """
    if not isinstance(table, Table):
        raise ValueError("table must be an astropy.table.Table object")
    if not isinstance(prefix, str):
        raise ValueError("prefix must be a string")

    for col in table.colnames:
        if (
            col.startswith("aper")
            or col.startswith("CI_")
            or col.startswith("isophotal_flux")
            or col.startswith("kron_flux")
            or "_psf" in col
        ):
            if prefix == "det_":
                table.remove_column(col)
            else:
                table.rename_column(col, prefix + col)

    return table
