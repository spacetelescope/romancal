# SIAF utilities

all = ['SIAF']

# Roman siaf definition
SIAF = None

def open_siaf(basepath=None, filename=None):
    """Open the pysiaf object"""
    global SIAF

    from pysiaf import Siaf
    SIAF = Siaf('roman', basepath=basepath, filename=filename)
