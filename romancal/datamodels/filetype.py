import io
import os
from pathlib import Path
from typing import Union


def check(init: Union[os.PathLike, Path, io.FileIO]) -> str:
    """
    Determine the type of a file and return it as a string

    Parameters
    ----------

    init : str
        file path or file object

    Returns
    -------
    file_type: str
        a string with the file type ("asdf" or "asn")

    """

    supported = ("asdf", "json")

    if isinstance(init, (str, os.PathLike, Path)):
        path, ext = os.path.splitext(init)
        ext = ext.strip(".")

        if not ext:
            raise ValueError(f"Input file path does not have an extension: {init}")

        if ext not in supported:  # Could be the file is zipped; try splitting again
            path, ext = os.path.splitext(path)
            ext = ext.strip(".")

            if ext not in supported:
                raise ValueError(f"Unrecognized file type for: {init}")

        if ext == "json":  # Assume json input is an association
            return "asn"

        return ext
    elif hasattr(init, "read") and hasattr(init, "seek"):
        magic = init.read(5)
        init.seek(0, 0)

        if not magic or len(magic) < 5:
            raise ValueError(f"Cannot get file type of {str(init)}")

        if magic == b"#ASDF":
            return "asdf"

        return "asn"
    else:
        raise ValueError(f"Cannot get file type of {str(init)}")
