import copy
import json
import os

import numpy as np
import pytest
from astropy import coordinates as coord
from astropy import units as u
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


def setup_source_catalog(img):
    """
    Set up the source catalog.

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

    # generate a set of random shifts to be added to the original coordinates
    seed = 13
    rng = default_rng(seed)
    shift_x = rng.uniform(-0.5, 0.5, size=len(source_catalog))
    shift_y = rng.uniform(-0.5, 0.5, size=len(source_catalog))
    # add random fraction of a pixel shifts to the centroid coordinates
    source_catalog["x_centroid"] += shift_x
    source_catalog["y_centroid"] += shift_y

    # generate another set of random shifts to be added to the original coordinates
    seed = 5
    rng = default_rng(seed)
    shift_x = rng.uniform(-0.5, 0.5, size=len(source_catalog))
    shift_y = rng.uniform(-0.5, 0.5, size=len(source_catalog))
    # add random fraction of a pixel shifts to the centroid coordinates
    source_catalog["x_psf"] += shift_x
    source_catalog["y_psf"] += shift_y

    # calculate centroid world coordinates
    centroid = img.meta.wcs(
        source_catalog["x_centroid"],
        source_catalog["y_centroid"],
    )
    # calculate PSF world coordinates
    psf = img.meta.wcs(
        source_catalog["x_psf"],
        source_catalog["y_psf"],
    )
    # add world coordinates to catalog
    source_catalog["ra_centroid"], source_catalog["dec_centroid"] = centroid
    source_catalog["ra_psf"], source_catalog["dec_psf"] = psf
    # add units
    source_catalog["ra_centroid"].unit = u.deg
    source_catalog["dec_centroid"].unit = u.deg
    source_catalog["ra_psf"].unit = u.deg
    source_catalog["dec_psf"].unit = u.deg

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


@pytest.mark.parametrize(
    "update_coordinates, should_update",
    [
        (True, True),  # coordinates should be updated
        (False, False),  # coordinates should NOT be updated
    ],
    ids=["update_enabled", "update_disabled"],
)
def test_tweakreg_catalog_coordinate_update_behavior(
    function_jail, tweakreg_image, update_coordinates, should_update
):
    """Test that TweakRegStep updates (or doesn't update) catalog coordinates based on parameter."""

    img = tweakreg_image(shift_1=1000, shift_2=1000, catalog_filename="img_1")

    # create and save source catalog
    source_catalog = setup_source_catalog(img)
    source_catalog.write("img_1_cat.parquet", overwrite=True)
    img.meta.source_catalog.tweakreg_catalog_name = "img_1_cat.parquet"

    # save original catalog for comparison
    cat_original = Table.read("img_1_cat.parquet")

    # run TweakRegStep with specified parameter
    res = TweakRegStep.call([img], update_source_catalog_coordinates=update_coordinates)

    # read catalog after tweakreg
    cat_after = Table.read("img_1_cat.parquet")

    with res:
        dm = res.borrow(0)

        # verify step completed successfully
        assert dm.meta.cal_step.tweakreg == "COMPLETE"

        if should_update:
            # verify coordinates were updated with tweaked WCS
            expected_centroid = dm.meta.wcs(
                cat_after["x_centroid"], cat_after["y_centroid"]
            )
            expected_psf = dm.meta.wcs(cat_after["x_psf"], cat_after["y_psf"])

            # The updated catalog should match the tweaked WCS exactly
            # (ignoring floating point precision issues)
            np.testing.assert_array_almost_equal(
                cat_after["ra_centroid"], expected_centroid[0]
            )
            np.testing.assert_array_almost_equal(
                cat_after["dec_centroid"], expected_centroid[1]
            )
            np.testing.assert_array_almost_equal(cat_after["ra_psf"], expected_psf[0])
            np.testing.assert_array_almost_equal(cat_after["dec_psf"], expected_psf[1])

            # Calculate actual differences between original and updated coordinates
            diff_ra_centroid = np.abs(
                cat_original["ra_centroid"] - cat_after["ra_centroid"]
            )
            diff_dec_centroid = np.abs(
                cat_original["dec_centroid"] - cat_after["dec_centroid"]
            )
            diff_ra_psf = np.abs(cat_original["ra_psf"] - cat_after["ra_psf"])
            diff_dec_psf = np.abs(cat_original["dec_psf"] - cat_after["dec_psf"])

            # New coordinates must change by at least 1 mas
            min_change_threshold_mas = 1.0  # mas
            min_change_threshold_deg = (
                (min_change_threshold_mas * u.mas).to(u.deg).value
            )

            # Verify coordinates changed by at least 1 mas
            assert np.any(diff_ra_centroid > min_change_threshold_deg)
            assert np.any(diff_dec_centroid > min_change_threshold_deg)
            assert np.any(diff_ra_psf > min_change_threshold_deg)
            assert np.any(diff_dec_psf > min_change_threshold_deg)

            # Verify changes are reasonable (within 1/2 a pixel)
            max_expected_change_mas = 55.0  # mas
            max_expected_change_deg = (max_expected_change_mas * u.mas).to(u.deg).value

            assert np.all(diff_ra_centroid < max_expected_change_deg)
            assert np.all(diff_dec_centroid < max_expected_change_deg)
            assert np.all(diff_ra_psf < max_expected_change_deg)
            assert np.all(diff_dec_psf < max_expected_change_deg)
        else:
            # verify coordinates were NOT updated (should be identical to original)
            np.testing.assert_array_equal(
                cat_original["ra_centroid"], cat_after["ra_centroid"]
            )
            np.testing.assert_array_equal(
                cat_original["dec_centroid"], cat_after["dec_centroid"]
            )
            np.testing.assert_array_equal(cat_original["ra_psf"], cat_after["ra_psf"])
            np.testing.assert_array_equal(cat_original["dec_psf"], cat_after["dec_psf"])

        res.shelve(dm, 0, modify=False)


def test_parquet_metadata_preserved_after_update(function_jail, tweakreg_image):
    """Test that parquet metadata survives the coordinate update operation."""
    import pyarrow.parquet as pq

    img = tweakreg_image(shift_1=1000, shift_2=1000, catalog_filename="img_1")

    # Create and save source catalog
    source_catalog = setup_source_catalog(img)
    source_catalog.write("img_1_cat.parquet", overwrite=True)
    img.meta.source_catalog.tweakreg_catalog_name = "img_1_cat.parquet"

    # Read original metadata using pyarrow
    original_table = pq.read_table("img_1_cat.parquet")
    original_schema_metadata = original_table.schema.metadata
    original_schema = original_table.schema

    # Run TweakRegStep with coordinate updates
    TweakRegStep.call([img], update_source_catalog_coordinates=True)

    # Read updated file's metadata
    updated_table = pq.read_table("img_1_cat.parquet")
    updated_schema_metadata = updated_table.schema.metadata
    updated_schema = updated_table.schema

    # Verify schema metadata is preserved
    assert original_schema_metadata == updated_schema_metadata, (
        "Schema metadata was not preserved during coordinate update"
    )

    # Verify column schemas match (types, names, field metadata)
    for orig_field, updated_field in zip(original_schema, updated_schema, strict=False):
        assert orig_field.name == updated_field.name
        # Note: We expect RA/Dec columns to have the same type, just different values
        # So schema types should still match
        assert orig_field.type == updated_field.type, (
            f"Column {orig_field.name} type changed from {orig_field.type} to {updated_field.type}"
        )
