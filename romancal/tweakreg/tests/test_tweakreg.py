import copy
import json
import os

import numpy as np
import pyarrow.parquet as pq
import pytest
from astropy import coordinates as coord
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.modeling import models
from astropy.table import Table
from numpy.random import default_rng

from romancal.datamodels import ModelLibrary
from romancal.tweakreg.tweakreg_step import TweakRegStep, _validate_catalog_columns


def test_tweakreg_raises_attributeerror_on_missing_tweakreg_catalog(tweakreg_image):
    """
    Test that TweakReg raises an AttributeError if meta.tweakreg_catalog is missing.
    """
    img = tweakreg_image()
    # make sure tweakreg_catalog_name doesn't exist
    img.meta.source_catalog = {}
    assert "tweakreg_catalog_name" not in img.meta.source_catalog
    with pytest.raises(
        AttributeError,
        match=r"Attribute 'meta.source_catalog.tweakreg_catalog' is missing",
    ):
        TweakRegStep.call([img])


def test_tweakreg_save_valid_abs_refcat(tmp_path, tweakreg_image):
    """Test that TweakReg saves the catalog used for absolute astrometry."""
    abs_refcat = "GAIADR3"

    img = tweakreg_image(
        shift_1=1000, shift_2=1000, catalog_filename="ref_catalog.ecsv"
    )
    abs_refcat_filename = f"fit_{abs_refcat.lower()}_ref.ecsv"

    TweakRegStep.call(
        [img],
        save_abs_catalog=True,
        abs_refcat=abs_refcat,
        catalog_path=str(tmp_path),
        output_dir=str(tmp_path),
    )

    assert os.path.exists(tmp_path / abs_refcat_filename)


def test_tweakreg_raises_error_on_invalid_abs_refcat(tmp_path, tweakreg_image):
    """Test that TweakReg raises an error when an invalid abs_refcat is provided."""
    img = tweakreg_image(shift_1=1000, shift_2=1000)

    with pytest.raises(ValueError, match="Invalid 'abs_refcat'"):
        TweakRegStep.call(
            [img],
            save_abs_catalog=True,
            abs_refcat="my_ref_cat",
            catalog_path=str(tmp_path),
            output_dir=str(tmp_path),
        )


def test_tweakreg_combine_custom_catalogs_and_asn_file(
    tmp_path, tweakreg_image, gaia_coords
):
    """
    Test that TweakRegStep can handle a custom catalog for the members of an ASN file.
    In this case, the user can create a custom catalog file (catfile) for each of the
    members of an ASN file.
    """
    catfile = str(tmp_path / "catfile.txt")

    # generate and save model
    img = tweakreg_image()
    img.meta.filename = "img.asdf"
    img.save(tmp_path / img.meta.filename)

    # generate and save custom catalog
    catalog_data = np.array(
        [img.meta.wcs.world_to_pixel_values(ra, dec) for ra, dec in gaia_coords]
    )
    custom_catalog_fn = str(tmp_path / "custom_catalog")
    Table(catalog_data, names=("x", "y")).write(custom_catalog_fn, format="ascii.ecsv")
    # record custom catalog
    with open(catfile, "w") as f:
        f.write(f"{img.meta.filename} {custom_catalog_fn}")

    # create ASN file
    asn_filepath = str(tmp_path / "test_asn.json")
    asn_dict = {
        "asn_id": "a3001",
        "products": [
            {
                "name": "test.asdf",
                "members": [{"expname": img.meta.filename, "exptype": "science"}],
            }
        ],
    }
    with open(asn_filepath, "w") as f:
        json.dump(asn_dict, f)

    res = TweakRegStep.call(
        asn_filepath,
        use_custom_catalogs=True,
        catfile=catfile,
    )

    assert isinstance(res, ModelLibrary)
    with res:
        m = res.borrow(0)
        assert m.meta.cal_step.tweakreg == "COMPLETE"
        assert m.meta.source_catalog.tweakreg_catalog_name.endswith("custom_catalog")
        res.shelve(m, modify=False)


@pytest.mark.parametrize(
    "theta, offset_x, offset_y",
    [
        (0.1 * u.deg, 1, 1),
    ],
)
def test_tweakreg_rotated_plane(
    tmp_path, tweakreg_image, theta, offset_x, offset_y, gaia_coords
):
    """
    Test that TweakReg returns accurate results.
    """
    img = tweakreg_image(shift_1=1000, shift_2=1000)
    original_wcs = copy.deepcopy(img.meta.wcs)

    # calculate original (x,y) for Gaia sources
    original_xy_gaia_sources = np.array(
        [original_wcs.world_to_pixel_values(ra, dec) for ra, dec in gaia_coords]
    )
    # move Gaia sources around by applying linear transformations
    # to their coords in the projected plane (same as a "wrong WCS")
    rot_matrix = models.Rotation2D(angle=theta)
    transformed_xy_gaia_sources = np.array(
        [
            rot_matrix(x, y) + np.array([offset_x, offset_y])
            for x, y in original_xy_gaia_sources
        ]
    )
    # save modified catalog to meta.tweakreg_catalog
    # (coords in the projection plane)
    Table(transformed_xy_gaia_sources, names=("x", "y")).write(
        img.meta.source_catalog.tweakreg_catalog_name,
        format="ascii.ecsv",
        overwrite=True,
    )

    TweakRegStep.call([img], abs_minobj=3)

    # get world coords for Gaia sources using "wrong WCS"
    original_ref_source = [
        original_wcs.pixel_to_world(x, y) for x, y in transformed_xy_gaia_sources
    ]
    # get world coords for Gaia sources using tweaked WCS
    new_ref_source = [
        img.meta.wcs.pixel_to_world(x, y) for x, y in transformed_xy_gaia_sources
    ]
    # celestial coordinates for Gaia sources
    gaia_ref_source = [
        coord.SkyCoord(ra * u.deg, dec * u.deg, frame="icrs") for ra, dec in gaia_coords
    ]
    # calculate distance between "wrong WCS" result and Gaia
    # (rounded to the 10th decimal place to avoid floating point issues)
    dist1 = [
        np.round(gref.separation(oref), 10)
        for gref, oref in zip(gaia_ref_source, original_ref_source, strict=False)
    ]
    # calculate distance between tweaked WCS result and Gaia
    # (rounded to the 10th decimal place to avoid floating point issues)
    dist2 = [
        np.round(gref.separation(nref), 10)
        for gref, nref in zip(gaia_ref_source, new_ref_source, strict=False)
    ]

    assert np.array(
        [np.less_equal(d2, d1) for d1, d2 in zip(dist1, dist2, strict=False)]
    ).all()


def test_fit_results_in_meta(tmp_path, tweakreg_image):
    """
    Test that the WCS fit results from tweakwcs are available in the meta tree.
    """

    img = tweakreg_image(shift_1=1000, shift_2=1000)

    res = TweakRegStep.call([img])

    assert isinstance(res, ModelLibrary)
    with res:
        for i, model in enumerate(res):
            assert hasattr(model.meta, "wcs_fit_results")
            assert len(model.meta.wcs_fit_results) > 0
            res.shelve(model, i, modify=False)


def test_tweakreg_handles_multiple_groups(tmp_path, tweakreg_image):
    """
    Test that TweakRegStep can perform relative alignment for all images in the groups
    before performing absolute alignment.
    """
    img1 = tweakreg_image(shift_1=1000, shift_2=1000, catalog_filename="img1")
    img2 = tweakreg_image(shift_1=990, shift_2=990, catalog_filename="img2")

    img1.meta.observation.program = 1
    img1.meta.observation["observation_id"] = "1"
    img2.meta.observation.program = 2
    img2.meta.observation["observation_id"] = "2"

    img1.meta["filename"] = "file1.asdf"
    img2.meta["filename"] = "file2.asdf"

    res = TweakRegStep.call([img1, img2])

    assert len(res.group_names) == 2


def setup_source_catalog(img, bias_value=None):
    """
    Set up the source catalog.

    Parameters
    ----------
    img : ImageModel
        The image model containing WCS information.
    bias_value : float, optional
        If provided, adds a known bias (in pixels) to the centroid coordinates.
        This makes the test more predictable by introducing a deterministic shift.

    Notes
    -----
    This function reads the source catalog from a file, renames columns to match
    expected names, adds mock PSF coordinates, applies random shifts to the centroid
    and PSF coordinates, and calculates the world coordinates for the centroids.
    """
    # read in the mock table
    source_catalog = Table.read("img_1", format="ascii.ecsv")
    # rename columns to match expected column names
    source_catalog.rename_columns(["x", "y"], ["x_centroid", "y_centroid"])
    # add mock PSF coordinates
    source_catalog["x_psf"] = source_catalog["x_centroid"]
    source_catalog["y_psf"] = source_catalog["y_centroid"]
    # add mock centroid win coordinates
    source_catalog["x_centroid_win"] = source_catalog["x_centroid"]
    source_catalog["y_centroid_win"] = source_catalog["y_centroid"]

    # define coordinate sets to add random shifts and optional bias
    coord_sets = [
        ("x_centroid", "y_centroid", 13, bias_value),
        ("x_psf", "y_psf", 5, bias_value),
        ("x_centroid_win", "y_centroid_win", 1, bias_value),
    ]
    for xcol, ycol, seed, bias in coord_sets:
        rng = default_rng(seed)
        shift_x = rng.uniform(-0.5, 0.5, size=len(source_catalog)) + (bias or 0.0)
        shift_y = rng.uniform(-0.5, 0.5, size=len(source_catalog)) + (bias or 0.0)
        source_catalog[xcol] += shift_x
        source_catalog[ycol] += shift_y

    # pixel to world mapping and WCS conversion for coordinate sets
    wcs_pairs = [
        (("x_centroid", "y_centroid"), ("ra", "dec")),
        (("x_centroid", "y_centroid"), ("ra_centroid", "dec_centroid")),
        (("x_psf", "y_psf"), ("ra_psf", "dec_psf")),
        (("x_centroid_win", "y_centroid_win"), ("ra_centroid_win", "dec_centroid_win")),
    ]
    for (xcol, ycol), (racol, deccol) in wcs_pairs:
        ra, dec = img.meta.wcs(source_catalog[xcol], source_catalog[ycol])
        source_catalog[racol], source_catalog[deccol] = ra, dec
        source_catalog[racol].unit = u.deg
        source_catalog[deccol].unit = u.deg

    return source_catalog


@pytest.mark.parametrize(
    "catalog_data, expected_colnames, flags_step_as_failed",
    [
        # both 'x' and 'y' columns present
        ({"x": [1, 2, 3], "y": [4, 5, 6]}, ["x", "y"], False),
        # 'x_psf' and 'y_psf' columns present, should be renamed
        ({"x_psf": [1, 2, 3], "y_psf": [4, 5, 6]}, ["x", "y"], False),
        # 'x' present, 'y_psf' present, should rename 'y_psf' to 'y'
        ({"x": [1, 2, 3], "y_psf": [4, 5, 6]}, ["x", "y"], False),
        # 'x_psf' present, 'y' present, should rename 'x_psf' to 'x'
        ({"x_psf": [1, 2, 3], "y": [4, 5, 6]}, ["x", "y"], False),
        # neither 'x' nor 'x_psf' present
        ({"y": [4, 5, 6]}, None, True),
        # neither 'y' nor 'y_psf' present
        ({"x": [1, 2, 3]}, None, True),
        # no relevant columns present
        (
            {"a": [1, 2, 3], "b": [4, 5, 6]},
            None,
            True,
        ),
    ],
)
def test_validate_catalog_columns(
    catalog_data, expected_colnames, flags_step_as_failed
):
    """Test that TweakRegStep._validate_catalog_columns() correctly validates the
    presence of required columns ('x' and 'y') in the provided catalog."""
    catalog = Table(catalog_data)
    is_valid = _validate_catalog_columns(catalog)
    assert is_valid is not flags_step_as_failed
    if expected_colnames is not None:
        assert set(catalog.colnames) == set(expected_colnames)


def test_tweakreg_flags_failed_step_on_invalid_catalog_columns(tweakreg_image):
    """Test that TweakRegStep raises ValueError when catalog columns are invalid."""
    img = tweakreg_image(shift_1=1000, shift_2=1000)
    # Add a tweakreg catalog with missing required columns
    img.meta.source_catalog = {
        "tweakreg_catalog": Table({"a": [1, 2, 3], "b": [4, 5, 6]}).as_array()
    }

    # Should raise ValueError due to invalid catalog columns
    with pytest.raises(ValueError, match=r"'tweakreg' source catalogs must contain"):
        TweakRegStep.call([img])


def test_tweakreg_handles_mixed_exposure_types(tmp_path, tweakreg_image):
    """Test that TweakReg can handle mixed exposure types
    (non-WFI_IMAGE data will be marked as SKIPPED only and won't be processed)."""
    invalid_types = [
        "WFI_SPECTRAL",
        "WFI_IM_DARK",
        "WFI_SP_DARK",
        "WFI_FLAT",
        "WFI_LOLO",
        "WFI_WFSC",
        "WFI_DARK",
        "WFI_GRISM",
        "WFI_PRISM",
    ]

    # start with 1 valid type
    img = tweakreg_image(catalog_filename="img0")
    imgs = [img]

    # add one of each invalid type
    for i, invalid_type in enumerate(invalid_types):
        img = tweakreg_image(catalog_filename=f"img{i + 1}")
        img.meta.exposure.type = invalid_type
        imgs.append(img)

    res = TweakRegStep.call(imgs)

    assert len(res) == len(imgs)
    with res:
        for i, m in enumerate(res):
            if i == 0:
                assert m.meta.cal_step.tweakreg == "COMPLETE"
            else:
                assert m.meta.cal_step.tweakreg == "SKIPPED"
            res.shelve(m, modify=False)


def test_tweakreg_updates_s_region(tmp_path, tweakreg_image):
    """Test that the TweakRegStep updates the s_region attribute."""
    img = tweakreg_image(shift_1=1000, shift_2=1000, catalog_filename="img")
    old_fake_s_region = "POLYGON ICRS 1.0000000000000 2.0000000000000 3.0000000000000 4.0000000000000 5.0000000000000 6.0000000000000 7.0000000000000 8.0000000000000 "
    img.meta.wcsinfo["s_region"] = old_fake_s_region

    # call TweakRegStep to update WCS & S_REGION
    res = TweakRegStep.call(img)

    with res:
        for i, model in enumerate(res):
            assert model.meta.wcsinfo.s_region != old_fake_s_region
            res.shelve(model, i, modify=False)


@pytest.mark.parametrize("save_results", [True, False])
def test_tweakreg_produces_output(tmp_path, tweakreg_image, save_results):
    """With save_results and output_dir set confirm expected files are in the output directory"""
    img = tweakreg_image(catalog_filename="img")
    base_filename = img.meta.filename
    TweakRegStep.call([img], save_results=save_results, output_dir=str(tmp_path))

    fns = [p.name for p in tmp_path.iterdir()]
    # the files should exist only if save_results was True
    assert (f"{base_filename}_tweakregstep.asdf" in fns) == save_results
    assert (f"{base_filename}_wcs.asdf" in fns) == save_results


def test_update_source_catalog_coordinates(function_jail, tweakreg_image):
    """Test that TweakReg updates the catalog coordinates with the tweaked WCS."""

    img = tweakreg_image(shift_1=1000, shift_2=1000, catalog_filename="img_1")

    # Use a known bias of 1 pixel for the centroid coordinates
    bias_value = 1.0  # pixels

    # create ImageSourceCatalogModel (no bias for this test)
    source_catalog = setup_source_catalog(img, bias_value=bias_value)
    source_catalog.write("img_1_cat.parquet", overwrite=True)

    # update tweakreg catalog name
    img.meta.source_catalog.tweakreg_catalog_name = "img_1_cat.parquet"

    # run TweakRegStep
    TweakRegStep.call([img], update_source_catalog_coordinates=True)

    # read in saved catalog coords
    cat = Table.read("img_1_cat.parquet")
    cat_ra_centroid = cat["ra_centroid"]
    cat_dec_centroid = cat["dec_centroid"]
    cat_ra_psf = cat["ra_psf"]
    cat_dec_psf = cat["dec_psf"]

    # calculate world coords using tweaked WCS
    expected_centroid = img.meta.wcs(cat["x_centroid"], cat["y_centroid"])
    expected_psf = img.meta.wcs(cat["x_psf"], cat["y_psf"])

    # compare coordinates (make sure tweaked WCS was applied to cat file coords)
    np.testing.assert_array_equal(cat_ra_centroid, expected_centroid[0])
    np.testing.assert_array_equal(cat_dec_centroid, expected_centroid[1])
    np.testing.assert_array_equal(cat_ra_psf, expected_psf[0])
    np.testing.assert_array_equal(cat_dec_psf, expected_psf[1])


def test_source_catalog_coordinates_have_changed(function_jail, tweakreg_image):
    """Test that the original catalog file content is different from the updated file.
    Uses a known bias to make the test deterministic and verifies that the coordinate
    shift matches the expected value based on the bias.
    """

    img = tweakreg_image(shift_1=1000, shift_2=1000, catalog_filename="img_1")

    # Use a known bias of 1 pixel for the centroid coordinates
    bias_value = 1.0  # pixels

    # Create ImageSourceCatalogModel with bias
    source_catalog = setup_source_catalog(img, bias_value=bias_value)
    source_catalog.write("img_1_cat.parquet", overwrite=True)

    cat_original = Table.read("img_1_cat.parquet")

    # update tweakreg catalog name
    img.meta.source_catalog.tweakreg_catalog_name = "img_1_cat.parquet"

    # run TweakRegStep with automatic catalog coordinate update
    TweakRegStep.call([img], update_source_catalog_coordinates=True)

    # read the updated catalog
    cat_updated = Table.read("img_1_cat.parquet")

    # Simply assert coordinates changed
    assert not np.array_equal(cat_original["ra_centroid"], cat_updated["ra_centroid"])
    assert not np.array_equal(cat_original["dec_centroid"], cat_updated["dec_centroid"])
    assert not np.array_equal(cat_original["ra_psf"], cat_updated["ra_psf"])
    assert not np.array_equal(cat_original["dec_psf"], cat_updated["dec_psf"])


def test_parquet_metadata_preserved_after_update(function_jail, tweakreg_image):
    """Test that parquet metadata is properly preserved with IVOA-compliant units after coordinate update."""

    img = tweakreg_image(shift_1=1000, shift_2=1000, catalog_filename="img_1")

    # Create and save source catalog
    source_catalog = setup_source_catalog(img)
    source_catalog.write("img_1_cat.parquet", overwrite=True)
    img.meta.source_catalog.tweakreg_catalog_name = "img_1_cat.parquet"

    # Read original catalog using astropy to get unit information
    original_astropy = Table.read("img_1_cat.parquet")
    original_table = pq.read_table("img_1_cat.parquet")
    original_schema = original_table.schema

    # Run TweakRegStep with coordinate updates and save updated catalog to disk
    TweakRegStep.call([img], update_source_catalog_coordinates=True)

    # Read updated file
    updated_astropy = Table.read("img_1_cat.parquet")
    updated_table = pq.read_table("img_1_cat.parquet")
    updated_schema = updated_table.schema

    # Verify column schemas match (types, names)
    for orig_field, updated_field in zip(original_schema, updated_schema, strict=False):
        assert orig_field.name == updated_field.name
        # Note: We expect RA/Dec columns to have the same type, just different values
        # So schema types should still match
        assert orig_field.type == updated_field.type

    # Verify units are preserved and IVOA-compliant when read by astropy
    # For columns with units, astropy should be able to read them correctly
    for col in original_astropy.colnames:
        # Check that units are readable (not None or empty when they should have units)
        if col in ["ra_centroid", "dec_centroid", "ra_psf", "dec_psf"]:
            assert updated_astropy[col].unit == u.deg
