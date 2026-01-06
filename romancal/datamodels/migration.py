def update_model_version(model, *, close_on_update=False):
    if model.tag == model._node_type._default_tag:
        return model
    updated_model = model.__class__.create_from_model(model)
    if close_on_update:
        model.close()
    return updated_model
