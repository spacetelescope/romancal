# SIAF utilities


def open_siaf(basepath=None, filename=None):
    """Open the pysiaf object

    Parameters
    ----------
    basepath : str or None
        The folder in which the pysiaf XML data exists.
        None to use the package default

    filename : str or None
        The filename of the XML file to use in in `basepath`
        None to use the package default of `roman_siaf.xml`

    Returns
    -------
    siaf : `pysiaf.siaf.Siaf`
        The Siaf object

    """
    from pysiaf import Siaf

    siaf = Siaf("roman", basepath=basepath, filename=filename)

    return siaf
