from astropy.table import Table


def update_colnames(table, prefix):
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
            table.rename_column(col, prefix + col)
    return table
