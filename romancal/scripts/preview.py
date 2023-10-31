from pathlib import Path
from typing import Optional

import asdf
import astropy.coordinates
import astropy.units
import gwcs
import numpy
import typer
from stpreview.__main__ import percentile_normalization, write_image
from stpreview.downsample import downsample_asdf_to

# try:
#     import typer
#     from stpreview.__main__ import percentile_normalization, write_image
#     from stpreview.downsample import downsample_asdf_to
# except (ImportError, ModuleNotFoundError):
#     raise ImportError(
#         'SDP requirements not installed; do `pip install "romancal[sdp]"`'
#     )

app = typer.Typer()


def north_pole_angle(wcs: gwcs.WCS, pixel=(0, 0), ddec=0.1 * astropy.units.arcsec):
    """
    Computes counterclockwise angle between positive x-axis and sky North.


    Parameters
    ----------
    wcs : `~gwcs.WCS`
    pixel : tuple of float, optional
        Reference pixel (x,y) in image from which to find angle.
        Default is (0,0).
    ddec : `~astropy.coordinates.Angle`, optional.
        Small angular offset for computing direction of North.
        No need to change this value unless pixel is very close to North pole.
        Default is 0.1 arcseconds.


    Returns
    -------
    angle : `~astropy.coordinates.Angle`
        Angle between sky North and pixel coordinate x-axis,
        along tangent line of great circle running through
        pixel and sky North.
    """

    coord = wcs.pixel_to_world(*pixel)
    coord = coord.transform_to("icrs")
    north_coord = coord.directional_offset_by(0.0 * astropy.units.deg, ddec)

    north_pixel = numpy.asarray(wcs.world_to_pixel(north_coord))
    diff = north_pixel - numpy.asarray(pixel)

    return astropy.coordinates.Angle(
        numpy.rad2deg(numpy.arctan2(diff[1], diff[0])) * astropy.units.deg
    )


@app.command()
def preview(input: Path, output: Optional[Path] = None):
    if output is None:
        output = Path.cwd()
    if output.is_dir():
        output = output / f"{input.stem}.png"

    with asdf.open(input) as file:
        wcs = file["roman"]["meta"]["wcs"]

    shape = (1000, 1000)

    data = downsample_asdf_to(input=input, shape=shape)

    write_image(
        data,
        output,
        shape=shape,
        normalization=percentile_normalization(data, percentile=90),
        colormap="afmhot",
        north_arrow_angle=north_pole_angle(wcs).degree - 90,
    )


@app.command()
def thumbnail(input: Path, output: Optional[Path] = None):
    if output is None:
        output = Path.cwd()
    if output.is_dir():
        output = output / f"{input.stem}_thumb.png"

    shape = (300, 300)

    data = downsample_asdf_to(input=input, shape=shape)

    write_image(
        data,
        output,
        shape=shape,
        normalization=percentile_normalization(data, percentile=90),
        colormap="afmhot",
    )


def command():
    app()


if __name__ == "__main__":
    command()
