def insert_substring(original, insert_str, substring, before=True):
    """
    Insert ``insert_str`` into ``original`` before or after ``substring``.

    If ``substring`` is not found, then ``insert_str`` is appended to
    ``original``.

    Parameters
    ----------
    original : str
        Original string to modify.
    insert_str : str
        Substring to insert.
    substrings : str
        Substring to match.
    before : bool, optional
        If True, insert before the substring. If False, insert
        after the substring. Default is True.

    Returns
    -------
    result : str
        Modified string.
    """
    if (idx := original.find(substring)) != -1:
        if before:
            pos = idx
        else:
            pos = idx + len(substring)
        return original[:pos] + insert_str + original[pos:]

    return original + insert_str


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
        if colname in append_cols:
            table.rename_column(colname, f"{colname}_{filter_name}")
        else:
            for ext in insert_col_exts:
                if ext in colname:
                    before = False if ext == "_psf" else True
                    new_colname = insert_substring(
                        colname, "_" + filter_name, ext, before=before
                    )
                    table.rename_column(colname, new_colname)
                    break  # no need to check other ext

    return table
