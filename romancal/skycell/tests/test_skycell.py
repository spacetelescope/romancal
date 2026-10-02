"""Unit tests for skycell functions"""

from pathlib import Path

import numpy as np
import pytest
import roman_datamodels
from astropy.coordinates import SkyCoord
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


def assert_allclose_lonlat(actual: np.ndarray, desired: np.ndarray, rtol=1e-7, atol=0):
    assert_allclose(actual, desired)


def assert_corners_on_pixel_corners(
    wcsobj, radec_corners: np.ndarray, pixel_shape: tuple[int, int]
):
    """the stored corners of a skycell lie on the corners of its pixel grid

    Which corner goes with which is not checked; that depends on the handedness
    of the skymap. See `skymap.SkyMap.vparity`.
    """
    pixels = np.array(wcsobj.invert(*radec_corners.T, with_bounding_box=False)).T
    expected = np.array(
        [
            (-0.5, -0.5),
            (pixel_shape[0] - 0.5, -0.5),
            (pixel_shape[0] - 0.5, pixel_shape[1] - 0.5),
            (-0.5, pixel_shape[1] - 0.5),
        ]
    )

    def sorted_rows(xy: np.ndarray) -> np.ndarray:
        # round the sort keys; the inverse transform is not exact, so nominally
        # equal coordinates would otherwise order arbitrarily
        # by construction passing test means we are getting
        # numbers like (-0.5, -0.5) or (N - 0.5, N - 0.5)
        # so the rounding just has to be good enough to get us to that
        # level
        keys = np.round(xy, 3)
        return xy[np.lexsort((keys[:, 1], keys[:, 0]))]

    assert_allclose(sorted_rows(pixels), sorted_rows(expected), atol=1e-4)


def sky_handedness(wcsobj, x: float, y: float, delta: float = 1.0) -> float:
    """determinant of d(east, north)/d(x, y), positive for a mirror image

    Local, so it is meaningful for skycells near a pole, where right ascension
    along a pixel row is not monotonic.
    """
    ra, dec = wcsobj(x, y, with_bounding_box=False)
    cos_dec = np.cos(np.deg2rad(dec))

    def offset(dx: float, dy: float) -> tuple[float, float]:
        ra_offset, dec_offset = wcsobj(x + dx, y + dy, with_bounding_box=False)
        east = (((ra_offset - ra + 180) % 360) - 180) * cos_dec
        return east, dec_offset - dec

    (east_x, north_x), (east_y, north_y) = offset(delta, 0), offset(0, delta)
    return east_x * north_y - east_y * north_x


@pytest.fixture(scope="module")
def skymap_subset() -> skymap.SkyMap:
    """
    smaller subset to allow these tests
    to run without access to the full skymap from CRDS.
    """
    return skymap.SkyMap(DATA_DIRECTORY / "skymap_subset.asdf")


@pytest.fixture(scope="module")
def mirrored_skymap_subset(tmp_path_factory) -> skymap.SkyMap:
    """
    the same subset with `x_tangent` mirrored within each skycell, standing in
    for a future skymap delivery in the standard handedness
    """
    model = roman_datamodels.open(DATA_DIRECTORY / "skymap_subset.asdf")
    skycells = np.array(model.skycells)
    skycells["x_tangent"] = model.meta.nxy_skycell - 1 - skycells["x_tangent"]
    model.skycells = skycells

    path = tmp_path_factory.mktemp("skymap") / "skymap_subset_mirrored.asdf"
    model.save(path)
    return skymap.SkyMap(path)


@pytest.fixture(params=["delivered", "mirrored"])
def either_skymap_subset(
    request, skymap_subset, mirrored_skymap_subset
) -> skymap.SkyMap:
    """both handedness conventions, which the WCS tests below are blind to"""
    return {
        "delivered": skymap_subset,
        "mirrored": mirrored_skymap_subset,
    }[request.param]


@pytest.fixture()
def sample_skycells(skymap_subset) -> skymap.SkyCells:
    return skymap.SkyCells.from_names(SAMPLE_SKYCELL_NAMES, skymap=skymap_subset)


def test_skycell_from_name(skymap_subset):
    skycell = skymap.SkyCells.from_names(["135p90x50y57"], skymap=skymap_subset)

    assert skycell == skymap.SkyCells([999], skymap=skymap_subset)

    assert skycell.data == np.void(
        (
            "135p90x50y57",
            224.99999999999997,
            89.48668040110688,
            -90.0,
            -31100.5,
            2499.5,
            220.4041123081352,
            89.52333943360831,
            221.0384731174898,
            89.44716843824965,
            228.96152688251019,
            89.44716843824965,
            229.59588769186476,
            89.52333943360831,
        ),
        dtype=[
            ("name", "<U16"),
            ("ra_center", "<f8"),
            ("dec_center", "<f8"),
            ("orientat", "<f4"),
            ("x_tangent", "<f8"),
            ("y_tangent", "<f8"),
            ("ra_corn1", "<f8"),
            ("dec_corn1", "<f8"),
            ("ra_corn2", "<f8"),
            ("dec_corn2", "<f8"),
            ("ra_corn3", "<f8"),
            ("dec_corn3", "<f8"),
            ("ra_corn4", "<f8"),
            ("dec_corn4", "<f8"),
        ],
    )

    with pytest.raises(KeyError):
        skymap.SkyCells.from_names(["r274dp63x63y81"], skymap=skymap_subset)

    with pytest.raises(KeyError):
        skymap.SkyCells.from_names(["notaskycellname"], skymap=skymap_subset)


def test_skycell_from_asn(skymap_subset):
    skycell = skymap.SkyCells.from_asns(
        [DATA_DIRECTORY / "L3_mosaic_asn.json"], skymap=skymap_subset
    )
    assert skycell.names == ["000p86x69y62"]

    with pytest.raises(ValueError):
        skymap.SkyCells.from_asns(
            [DATA_DIRECTORY / "L3_regtest_asn.json"], skymap=skymap_subset
        )

    with pytest.raises(ValueError):
        skymap.SkyCells.from_asns(
            [DATA_DIRECTORY / "L3_skycell_mbcat_asn.json"], skymap=skymap_subset
        )
    with pytest.raises(ValueError):
        skymap.SkyCells.from_asns(
            DATA_DIRECTORY.glob("*_asn.json"), skymap=skymap_subset
        )


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

    projregion0 = skymap.ProjectionRegion(0, skymap=skymap_subset)
    projregion1 = skymap.ProjectionRegion(1, skymap=skymap_subset)

    assert len(skycell.projection_regions) == 1
    assert skycell.projection_regions[0] == projregion0.index  # this calls CRDS!

    assert (
        skymap.ProjectionRegion.from_skycell_index(107, skymap=skymap_subset)
        == projregion0
    )

    assert (
        skymap.ProjectionRegion.from_skycell_index(0, skymap=skymap_subset)
        == projregion0
    )

    assert (
        skymap.ProjectionRegion.from_skycell_index(
            projregion0.data["skycell_end"], skymap=skymap_subset
        )
        == projregion1
    )

    with pytest.raises(KeyError):
        skymap.ProjectionRegion.from_skycell_index(-1, skymap=skymap_subset)

    with pytest.raises(KeyError):
        skymap.ProjectionRegion.from_skycell_index(10000, skymap=skymap_subset)


@pytest.mark.parametrize("name", SAMPLE_SKYCELL_NAMES)
def test_skycell_wcs_pixel_to_world(name, either_skymap_subset):
    skycell = skymap.SkyCells.from_names([name], skymap=either_skymap_subset)

    wcsobj = skycell.wcs[0]

    # forward transform of the pixel corners covers the stored corners
    corners = SkyCoord(
        *wcsobj(
            *np.array(
                [
                    (-0.5, -0.5),
                    (skycell.pixel_shape[0] - 0.5, -0.5),
                    (skycell.pixel_shape[0] - 0.5, skycell.pixel_shape[1] - 0.5),
                    (-0.5, skycell.pixel_shape[1] - 0.5),
                ]
            ).T,
            with_bounding_box=False,
        ),
        unit="deg",
    )
    stored = SkyCoord(*skycell.radec_corners[0].T, unit="deg")

    # every stored corner has a computed corner on top of it; which one depends
    # on the handedness of the skymap
    separations = stored[:, None].separation(corners[None, :])
    assert_allclose(separations.min(axis=1).to("mas").value, 0, atol=1)


@pytest.mark.parametrize("name", SAMPLE_SKYCELL_NAMES)
def test_skycell_wcs_world_to_pixel(name, either_skymap_subset):
    skycell = skymap.SkyCells.from_names([name], skymap=either_skymap_subset)

    assert_corners_on_pixel_corners(
        skycell.wcs[0], skycell.radec_corners[0], skycell.pixel_shape
    )


@pytest.mark.parametrize("name", SAMPLE_SKYCELL_NAMES)
def test_skycell_wcsinfo(name, either_skymap_subset):
    skycell = skymap.SkyCells.from_names([name], skymap=either_skymap_subset)

    wcsobj = skycell.wcs[0]
    wcs_info = skycell.wcs_infos[0]

    assert_allclose(
        wcsobj(
            (wcs_info["nx"] / 2.0) - 0.5,
            (wcs_info["ny"] / 2.0) - 0.5,
        ),
        (wcs_info["ra_center"], wcs_info["dec_center"]),
        rtol=1e-7,
    )

    assert_corners_on_pixel_corners(
        wcsobj, skycell.radec_corners[0], (wcs_info["nx"], wcs_info["ny"])
    )


def test_skymap_vparity(skymap_subset, mirrored_skymap_subset):
    """the handedness of a skymap follows from its `x_tangent` values"""

    assert skymap_subset.vparity == 1
    assert mirrored_skymap_subset.vparity == -1


@pytest.mark.parametrize("name", SAMPLE_SKYCELL_NAMES)
def test_skycell_wcs_mirrored_skymap(name, skymap_subset, mirrored_skymap_subset):
    """a mirrored skymap gives the same skycell, flipped in x"""

    delivered = skymap.SkyCells.from_names([name], skymap=skymap_subset)
    mirrored = skymap.SkyCells.from_names([name], skymap=mirrored_skymap_subset)

    nx, ny = delivered.pixel_shape
    x, y = np.meshgrid(np.linspace(0, nx - 1, 4), np.linspace(0, ny - 1, 4))
    x, y = x.ravel(), y.ravel()

    assert_allclose_lonlat(
        np.array(mirrored.wcs[0](nx - 1 - x, y, with_bounding_box=False)),
        np.array(delivered.wcs[0](x, y, with_bounding_box=False)),
    )

    # the footprint is unchanged, so the stored corners still apply
    assert_allclose(mirrored.radec_corners, delivered.radec_corners)
    assert_corners_on_pixel_corners(
        mirrored.wcs[0], mirrored.radec_corners[0], mirrored.pixel_shape
    )

    # right ascension increases to the left, unlike in the delivered skymap
    center = ((nx - 1) / 2, (ny - 1) / 2)
    assert sky_handedness(delivered.wcs[0], *center) > 0
    assert sky_handedness(mirrored.wcs[0], *center) < 0


def test_skycells(skymap_subset):
    skycells = skymap.SkyCells.from_names(SAMPLE_SKYCELL_NAMES, skymap=skymap_subset)

    assert sorted(skycells.names) == sorted(SAMPLE_SKYCELL_NAMES)

    assert skycells.radec_corners.shape == (len(SAMPLE_SKYCELL_NAMES), 4, 2)
    assert skycells.vectorpoint_corners.shape == (len(SAMPLE_SKYCELL_NAMES), 4, 3)

    assert skycells.radec_centers.shape == (len(SAMPLE_SKYCELL_NAMES), 2)
    assert skycells.vectorpoint_centers.shape == (len(SAMPLE_SKYCELL_NAMES), 3)

    assert len(skycells.polygons) == len(SAMPLE_SKYCELL_NAMES)


def test_skycells_cores_containing_center(sample_skycells):
    assert np.all(sample_skycells.containing(sample_skycells.radec_centers))
    assert sample_skycells.cores_containing(sample_skycells.radec_centers) != {}


@pytest.mark.parametrize(
    "radec,expected",
    [
        (
            [68.5, 3.0],
            {},
        ),
        (
            [17.00495323, 86.23671728],
            {3190: [0]},
        ),
    ],
)
def test_skycells_cores_containing(radec, expected, sample_skycells):
    assert sample_skycells.cores_containing(radec) == expected
