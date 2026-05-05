"""Test suffix replacement"""

import pytest

from romancal.lib import suffix as s


def test_suffix_existence():
    """Generate current suffix list and compare"""

    calculated_suffixes = s._find_suffixes()
    found_suffixes = s._combine_suffixes(
        to_add=(calculated_suffixes, s._SUFFIXES_TO_ADD),
        to_remove=(s._SUFFIXES_TO_DISCARD,),
    )
    diff = set(found_suffixes).symmetric_difference(s._KNOW_SUFFIXES)
    assert not diff, f"Suffixes unaccounted for:\n{diff}"


@pytest.mark.parametrize("suffix", s._KNOW_SUFFIXES)
def test_suffix_removal(suffix):
    """Test suffix removal"""
    basename = "file"
    full_fpath = basename + "_" + suffix
    removed_path, separator = s.remove_suffix(full_fpath)
    assert removed_path == basename
    assert separator == "_"


@pytest.mark.parametrize("suffix", s._KNOW_SUFFIXES)
def test_suffix_replacement(suffix, base="file", new="junk", sep="_"):
    """Test suffix replacement"""
    full_path = base + sep + suffix
    replaced = s.replace_suffix(full_path, new)
    assert replaced == base + sep + new
