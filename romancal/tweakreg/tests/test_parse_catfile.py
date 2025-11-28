import pytest

from romancal.tweakreg.tweakreg_step import _parse_catfile


def test_parse_catfile_raises_error_on_invalid_content(tmp_path):
    """
    Test that _parse_catfile raises an error if catfile contains more than
    two columns (i.e. image name and corresponding catalog path).
    """
    # create custom catalog file and input datamodels
    catfile = str(tmp_path / "catfile.txt")
    with open(catfile, "w") as f:
        f.write("img1.asdf column1 column2 column3")

    with pytest.raises(ValueError):
        _parse_catfile(catfile)


def test_parse_catfile_valid_catalog(tmp_path):
    """
    Test that _parse_catfile can parse a custom catalog with valid format.
    """
    catfile = str(tmp_path / "catfile.txt")
    with open(catfile, mode="w") as f:
        f.write("img1.asdf cat1\nimg2.asdf cat2")
    catdict = _parse_catfile(catfile)

    assert catdict == {
        "img1.asdf": str(tmp_path / "cat1"),
        "img2.asdf": str(tmp_path / "cat2"),
    }


@pytest.mark.parametrize("catfile", (None, ""))
def test_parse_catfile_returns_none(catfile):
    """
    Test that _parse_catfile returns None when catfile = None or catfile = "".
    """
    catdict = _parse_catfile(catfile=catfile)

    assert catdict is None


def test_parse_catfile_returns_none_on_invalid_content(tmp_path):
    """
    Test that _parse_catfile returns a dict where all the values are None
    if only filename is present in catfile (i.e. no associated catalog).
    """
    # FIXME the tested behavior here causes the step to fail
    # create custom catalog file and input datamodels
    catfile = str(tmp_path / "catfile.txt")
    with open(catfile, mode="w") as f:
        f.write("img1.asdf\nimg2.asdf\nimg3.asdf")

    catdict = _parse_catfile(catfile)

    assert not all(catdict.values())
