__all__ = ["update_model_version"]


def update_model_version(model, *, close_on_update=False):
    """
    Update a DataModel instance to the newest version.

    If the DataModel is already the newest (tag) version return
    the same DataModel instance.

    Parameters
    ----------
    close_on_update : bool, optional
        If enabled, call model.close if an updated (and copied)
        DataModel is returned.

    Returns
    -------
    DataModel
        Either the provided DataModel or an updated (copy).
    """
    if model.tag == model._node_type._default_tag:
        return model

    updated_model = model.__class__.create_from_model(model)

    if close_on_update:
        model.close()

    return updated_model
