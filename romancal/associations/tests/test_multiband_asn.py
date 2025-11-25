import os
import sys

from romancal.associations import multiband_asn


def test_parse_file_list_wildcard(tmp_path):
    """Test that _parse_file_list expands wildcards and returns correct file list."""
    # Create dummy files
    filenames = ["test1.asdf", "test2.asdf"]
    for fname in filenames:
        (tmp_path / fname).touch()
    files = [str(tmp_path / "*.asdf")]
    assoc = multiband_asn.MultibandAssociation(files)
    parsed_files = assoc._parse_file_list(files)
    assert set(os.path.basename(f) for f in parsed_files) == set(filenames)


def test_get_skycell_groups():
    """Test that _get_skycell_groups correctly groups files by skycell id."""
    files = [
        "r00001_p_full_270p65x48y69_f123_coadd.asdf",
        "r00001_p_full_270p65x48y69_f456_coadd.asdf",
        "r00001_p_full_271p66x49y70_f123_coadd.asdf",
    ]
    assoc = multiband_asn.MultibandAssociation(files)
    groups = assoc._get_skycell_groups(files)
    assert set(groups.keys()) == {"270p65x48y69", "271p66x49y70"}
    assert set(groups["270p65x48y69"]) == {
        "r00001_p_full_270p65x48y69_f123_coadd.asdf",
        "r00001_p_full_270p65x48y69_f456_coadd.asdf",
    }


def test_create_multiband_asn_runs(monkeypatch):
    """Test that create_multiband_asn runs without error and calls asn_from_list._cli."""
    files = [
        "r00001_p_full_270p65x48y69_f123_coadd.asdf",
        "r00001_p_full_270p65x48y69_f456_coadd.asdf",
    ]
    assoc = multiband_asn.MultibandAssociation(files)
    # Patch asn_from_list._cli to a dummy function to avoid side effects
    import romancal.associations.multiband_asn as multiband_asn_mod

    monkeypatch.setattr(multiband_asn_mod.asn_from_list, "_cli", lambda args: None)
    assoc.create_multiband_asn()


def test_parse_file_list_no_wildcard():
    """Test that _parse_file_list returns the input list if no wildcard is present."""
    files = ["file1.asdf", "file2.asdf"]
    assoc = multiband_asn.MultibandAssociation(files)
    assert assoc._parse_file_list(files) == files


def test_get_skycell_groups_empty():
    """Test that _get_skycell_groups returns empty dict when no skycell patterns match."""
    files = ["not_a_skycell_file.txt", "another_file.fits"]
    assoc = multiband_asn.MultibandAssociation(files)
    groups = assoc._get_skycell_groups(files)
    assert groups == {}


def test_multiband_association_empty_list():
    """Test MultibandAssociation with an empty file list."""
    assoc = multiband_asn.MultibandAssociation([])
    assert assoc.files == []
    assert assoc.skycell_groups == {}


def test_cli_entry_point(monkeypatch, tmp_path):
    """Test that the CLI entry point runs without error with minimal arguments."""
    # Create dummy files
    filenames = ["r00001_p_full_270p65x48y69_f123_coadd.asdf"]
    for fname in filenames:
        (tmp_path / fname).touch()
    test_args = ["multiband_asn", str(tmp_path / "*.asdf")]
    monkeypatch.setattr(sys, "argv", test_args)
    # Patch asn_from_list._cli to avoid side effects
    monkeypatch.setattr(multiband_asn.asn_from_list, "_cli", lambda args: None)
    # Should run without error
    multiband_asn._cli()
