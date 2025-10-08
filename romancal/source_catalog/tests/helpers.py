import pyarrow


def compare_model_and_parquet_metadata(model, parquet_filename, ignore_paths):
    """
    Compare a selection of DataModel and parquet table metadata.
    """
    tbl = pyarrow.parquet.read_table(parquet_filename)
    metadata = tbl.schema.metadata
    assert metadata[b"roman.meta.filename"].decode("ascii") == parquet_filename
    assert metadata[b"roman.meta.image.filename"].decode("ascii") == model.meta.filename
    assert (
        metadata[b"roman.meta.image.file_date"].decode("ascii") == model.meta.file_date
    )
    for key, value in metadata.items():
        str_key = key.decode("ascii")
        if str_key in ignore_paths or any(str_key.startswith(p) for p in ignore_paths):
            continue
        path = str_key.split(".")
        assert path[0] == "roman"
        target = model
        for sub_path in path[1:]:
            if sub_path.isnumeric():
                target = target[int(sub_path)]
            else:
                target = getattr(target, sub_path)
        str_target = str(target)
        str_value = value.decode("ascii")
        assert str_value == str_target, (
            f"values for {key} do not match: {str_target}, {str_value}"
        )
