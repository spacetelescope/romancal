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
    Add a filter name to the column names in an astropy table
    for the multiband catalog.

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
    insert_col_exts = ("_flux", "_psf", "_abmag")
    append_cols = (
        "sharpness",
        "roundness1",
        "is_extended",
        "fluxfrac_radius_50",
        "psf_gof",
        "psf_flags",
    )

    for colname in table.colnames:
        if any(ext in colname for ext in insert_col_exts):
            new_colname = insert_filter(colname, filter_name)
            table.rename_column(colname, new_colname)
        elif colname in append_cols:
            table.rename_column(colname, f"{colname}_{filter_name}")

    return table
