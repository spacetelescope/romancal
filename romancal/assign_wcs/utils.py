def wcs_bbox_from_shape(shape):
    """Create a bounding box from the shape of the data.
    This is appropriate to attach to a wcs object.

    Parameters
    ----------
    shape : tuple
        The shape attribute from a `numpy.ndarray` array.

    Returns
    -------
    bbox : tuple
        Bounding box in x, y order.
    """
    bbox = ((-0.5, shape[-1] - 0.5),
            (-0.5, shape[-2] - 0.5))
    return bbox
