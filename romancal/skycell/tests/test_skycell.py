"""Unit tests for skycell functions"""

from pathlib import Path

import numpy as np
import pytest
import spherical_geometry.polygon as sgp
from numpy.testing import assert_allclose

from romancal.skycell import skymap

DATA_DIRECTORY = Path(__file__).parent / "data"

SAMPLE_SKYCELL_NAMES = [
    "000p86x30y34",
    "000p86x50y65",
    "000p86x59y38",
    # north pole
    "135p90x25y49",
    "135p90x30y51",
    "135p90x33y62",
    "135p90x39y33",
    "135p90x43y65",
    "135p90x48y41",
    "135p90x52y59",
    "135p90x57y35",
    "135p90x61y67",
    "135p90x67y38",
]


@pytest.fixture(scope="module")
def skymap_subset() -> skymap.SkyMap:
    """
    smaller subset to allow these tests
    to run without access to the full skymap from CRDS.
    """
    return skymap.SkyMap(DATA_DIRECTORY / "skymap_subset.asdf")


@pytest.fixture()
def sample_skycells(skymap_subset) -> skymap.SkyCells:
    return skymap.SkyCells.from_names(SAMPLE_SKYCELL_NAMES, skymap=skymap_subset)


@pytest.fixture()
def all_skycells(skymap_subset) -> skymap.SkyCells:
    """every skycell in the projection regions of the subset"""
    return skymap.SkyCells(
        np.arange(skymap_subset.model.projection_regions[-1]["skycell_end"]),
        skymap=skymap_subset,
    )


def test_skycell_from_name(skymap_subset):
    skycell = skymap.SkyCells.from_names(["135p90x50y57"], skymap=skymap_subset)

    assert skycell == skymap.SkyCells([999], skymap=skymap_subset)
    assert skycell.names == ["135p90x50y57"]
    assert_allclose(skycell.radec_centers, [[225, 89.48668040110688]])

    for name in ["r274dp63x63y81", "notaskycellname"]:
        with pytest.raises(KeyError):
            skymap.SkyCells.from_names([name], skymap=skymap_subset)


def test_skycell_from_asn(skymap_subset):
    skycell = skymap.SkyCells.from_asns(
        [DATA_DIRECTORY / "L3_mosaic_asn.json"], skymap=skymap_subset
    )
    assert skycell.names == ["000p86x69y62"]

    for asns in [
        [DATA_DIRECTORY / "L3_regtest_asn.json"],
        [DATA_DIRECTORY / "L3_skycell_mbcat_asn.json"],
        DATA_DIRECTORY.glob("*_asn.json"),
    ]:
        with pytest.raises(ValueError):
            skymap.SkyCells.from_asns(asns, skymap=skymap_subset)


def test_skycell_from_projregion(skymap_subset):
    projregion = skymap.ProjectionRegion(0, skymap=skymap_subset)

    assert skymap.SkyCells(
        projregion.skycell_indices[100], skymap=skymap_subset
    ) == skymap.SkyCells.from_names(["135p90x30y44"], skymap=skymap_subset)

    assert (
        projregion.skycell_indices[-1]
        != skymap.ProjectionRegion(1, skymap=skymap_subset).skycell_indices[0]
    )


def test_projregion_from_skycell(skymap_subset):
    skycell = skymap.SkyCells.from_names(["135p90x50y57"], skymap=skymap_subset)
    assert skycell.projection_regions.tolist() == [0]

    first_of_region_1 = skymap.ProjectionRegion(0, skymap=skymap_subset).data[
        "skycell_end"
    ]
    for skycell_index, projregion_index in [(107, 0), (0, 0), (first_of_region_1, 1)]:
        assert skymap.ProjectionRegion.from_skycell_index(
            skycell_index, skymap=skymap_subset
        ) == skymap.ProjectionRegion(projregion_index, skymap=skymap_subset)

    for skycell_index in [-1, 10000]:
        with pytest.raises(KeyError):
            skymap.ProjectionRegion.from_skycell_index(
                skycell_index, skymap=skymap_subset
            )


@pytest.mark.parametrize("name", SAMPLE_SKYCELL_NAMES)
def test_skycell_wcs(name, skymap_subset):
    skycell = skymap.SkyCells.from_names([name], skymap=skymap_subset)
    wcsobj = skycell.wcs[0]
    wcs_info = skycell.wcs_infos[0]
    nx, ny = skycell.pixel_shape
    pixel_corners = np.array(
        [(-0.5, -0.5), (nx - 0.5, -0.5), (nx - 0.5, ny - 0.5), (-0.5, ny - 0.5)]
    )

    # the pixel grid's corners are the reference file's corners, both ways
    assert_allclose(
        np.array(wcsobj(*pixel_corners.T, with_bounding_box=False)).T,
        skycell.radec_corners[0],
        rtol=1e-7,
    )
    assert_allclose(
        np.array(wcsobj.invert(*skycell.radec_corners[0].T, with_bounding_box=False)).T,
        pixel_corners,
        rtol=1e-5,
    )
    # and its center is the reference file's center
    assert_allclose(
        wcsobj(wcs_info["nx"] / 2 - 0.5, wcs_info["ny"] / 2 - 0.5),
        (wcs_info["ra_center"], wcs_info["dec_center"]),
        rtol=1e-7,
    )


def test_skycells(sample_skycells):
    count = len(SAMPLE_SKYCELL_NAMES)
    assert sorted(sample_skycells.names) == sorted(SAMPLE_SKYCELL_NAMES)
    assert sample_skycells.radec_corners.shape == (count, 4, 2)
    assert sample_skycells.vectorpoint_corners.shape == (count, 4, 3)
    assert sample_skycells.radec_centers.shape == (count, 2)
    assert sample_skycells.vectorpoint_centers.shape == (count, 3)
    assert len(sample_skycells.polygons) == count


def test_skycells_containing_centers(sample_skycells):
    own_centers = {
        int(index): [position] for position, index in enumerate(sample_skycells.indices)
    }

    # each skycell contains its own center, and perhaps other centers too
    containing = sample_skycells.containing(sample_skycells.radec_centers)
    assert containing.keys() == own_centers.keys()
    assert all(own_centers[index][0] in points for index, points in containing.items())

    # each center is in its own core, except that of 135p90x25y49, which lies
    # outside the bounds of its projection region; a skycell of the
    # neighboring region, not among these, owns it
    cores_containing = sample_skycells.cores_containing(sample_skycells.radec_centers)
    del own_centers[
        sample_skycells.indices[sample_skycells.names.index("135p90x25y49")]
    ]
    assert cores_containing == own_centers


def test_skycells_containing(all_skycells):
    rng = np.random.default_rng(3)
    # around projection regions 0 and 1, including their edges
    radec = np.stack(
        [rng.uniform(-25, 25, 200) % 360, rng.uniform(84.4, 90, 200)], axis=1
    )

    containing = all_skycells.containing(radec)

    for point, point_radec in enumerate(radec):
        # a skycell that contains the point has its center within 0.055 degrees
        nearby = np.flatnonzero(
            skymap._separation(
                all_skycells.vectorpoint_centers, skymap._vectorpoints(point_radec)
            )
            < 0.06
        )
        assert {index for index, points in containing.items() if point in points} == {
            int(index)
            for index in nearby
            if sgp.SingleSphericalPolygon(
                all_skycells.vectorpoint_corners[index],
                all_skycells.vectorpoint_centers[index],
            ).contains_lonlat(*point_radec)
        }


@pytest.mark.parametrize(
    "radec,expected",
    [
        ([68.5, 3.0], {}),
        ([17.00495323, 86.23671728], {"000p86x50y65": [0]}),
        # the pole belongs to the skycell centered on it
        ([123.0, 90.0], {"135p90x50y50": [0]}),
        ([0.0, 90.0], {"135p90x50y50": [0]}),
    ],
)
def test_skycells_cores_containing(radec, expected, all_skycells, skymap_subset):
    assert {
        skymap.SkyCells([index], skymap=skymap_subset).names[0]: points
        for index, points in all_skycells.cores_containing(radec).items()
    } == expected
