import numpy as np
import pytest
from astropy import coordinates as coord
from astropy import units as u
from astropy.modeling.models import RotationSequence3D, Scale, Shift
from astropy.table import Table
from astropy.time import Time
from gwcs import coordinate_frames as cf
from gwcs import wcs
from gwcs.geometry import CartesianToSpherical, SphericalToCartesian
from roman_datamodels import datamodels as rdm
from stcal.tweakreg.astrometric_utils import get_catalog

from romancal.tweakreg.tweakreg_step import DEFAULT_ABS_REFCAT


@pytest.fixture(scope="session")
def gaia_coords():
    gaia_cat = get_catalog(
        right_ascension=270,
        declination=66,
        search_radius=100 / 3600.0,
        catalog=DEFAULT_ABS_REFCAT,
        timeout=120,
    )
    return [(ra, dec) for ra, dec in zip(gaia_cat["ra"], gaia_cat["dec"], strict=False)]


@pytest.fixture
def tweakreg_image(tmp_path, gaia_coords):
    """
    Create a base image with a realistic WCSInfo and a WCS.

    Notes
    -----
    The size of the image needs to be relatively large in order for
    the source catalog step to find a reasonable number of sources in the image.

    shift_1 and shift_2 (units in pixel) are used to shift the WCS projection plane.
    """

    def _tweakreg_image(shift_1=0, shift_2=0, catalog_filename="catalog.ecsv"):
        l2 = rdm.ImageModel.create_fake_data(shape=(2000, 2000))
        l2.meta.filename = "none"
        l2.meta.cal_step = {}
        for step_name in l2.schema_info("required")["roman"]["meta"]["cal_step"][
            "required"
        ].info:
            l2.meta.cal_step[step_name] = "INCOMPLETE"
        l2.meta.cal_logs = []

        # to match GAIADR3 epoch
        l2.meta.exposure.start_time = Time("2016-01-01T00:00:00")

        # update wcsinfo
        l2.meta.wcsinfo.v2_ref = 0.42955128282521254
        l2.meta.wcsinfo.v3_ref = -0.2479976768255853
        l2.meta.wcsinfo.vparity = -1
        l2.meta.wcsinfo.v3yangle = -999999
        l2.meta.wcsinfo.ra_ref = 270.0
        l2.meta.wcsinfo.dec_ref = 66.0
        l2.meta.wcsinfo.roll_ref = 60
        l2.meta.wcsinfo.s_region = (
            "POLYGON ICRS "
            "269.3318903230621 65.56866666048172 "
            "269.32578768154605 65.69246311613287 "
            "269.02457173246125 65.69201346248587 "
            "269.0333096074621 65.56870823657276 "
        )

        # create necessary transformations
        distortion = Shift(-shift_1) & Shift(-shift_2)
        v2_ref = l2.meta.wcsinfo.v2_ref / 3600
        v3_ref = l2.meta.wcsinfo.v3_ref / 3600
        roll_ref = l2.meta.wcsinfo.roll_ref
        ra_ref = l2.meta.wcsinfo.ra_ref
        dec_ref = l2.meta.wcsinfo.dec_ref

        angles = np.array([v2_ref, -v3_ref, roll_ref, dec_ref, -ra_ref])
        axes = "zyxyz"
        rot = RotationSequence3D(angles, axes_order=axes)

        # The sky rotation expects values in deg.
        tel2sky = (
            (Scale(0.1 / 3600) & Scale(0.1 / 3600))
            | SphericalToCartesian(wrap_lon_at=180)
            | rot
            | CartesianToSpherical(wrap_lon_at=360)
        )
        tel2sky.name = "v23tosky"

        # create required frames
        detector = cf.Frame2D(name="detector", axes_order=(0, 1), unit=(u.pix, u.pix))
        v2v3 = cf.Frame2D(
            name="v2v3",
            axes_order=(0, 1),
            axes_names=("v2", "v3"),
            unit=(u.arcsec, u.arcsec),
        )
        world = cf.CelestialFrame(reference_frame=coord.ICRS(), name="world")

        # create pipeline
        pipeline = [
            wcs.Step(detector, distortion),
            wcs.Step(v2v3, tel2sky),
            wcs.Step(world, None),
        ]

        wcs_obj = wcs.WCS(pipeline)
        wcs_obj.bounding_box = (
            (-0.5, l2.data.shape[-2] + 0.5),
            (-0.5, l2.data.shape[-1] + 0.5),
        )

        l2.meta.wcs = wcs_obj

        catalog_path = tmp_path / catalog_filename
        catalog_data = np.array(
            [wcs_obj.world_to_pixel_values(ra, dec) for ra, dec in gaia_coords]
        )
        Table(catalog_data, names=("x", "y")).write(catalog_path, format="ascii.ecsv")
        l2.meta["source_catalog"] = {
            "tweakreg_catalog_name": str(catalog_path),
        }

        return l2

    return _tweakreg_image
