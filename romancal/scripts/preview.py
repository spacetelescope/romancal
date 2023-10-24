from pathlib import Path
from typing import Union

import numpy
import typer
from astropy.visualization import (
    BaseStretch,
    ImageNormalize,
    LinearStretch,
    PercentileInterval,
)
from matplotlib import pyplot
from matplotlib.colors import Colormap
from stpreview.downsample import downsample_asdf_to

app = typer.Typer()


def percentile_normalization(
    data: numpy.ndarray, percentile: float = 90, stretch: BaseStretch = None
) -> ImageNormalize:
    """
    stretch the given data based on the given percentile
    """

    if stretch is None:
        stretch = LinearStretch()

    interval = PercentileInterval(percentile)
    vmin, vmax = interval.get_limits(data)

    normalization = ImageNormalize(vmin=vmin, vmax=vmax, stretch=stretch)

    return normalization


def write_image(
    data: numpy.ndarray,
    output: Path,
    normalization: ImageNormalize = None,
    colormap: Union[str, Colormap] = None,
):
    """
    write data as an image to the given path
    """

    if normalization is None:
        normalization = percentile_normalization(data)

    if colormap is None:
        colormap = "afmhot"

    figure = pyplot.figure()
    axis = figure.add_subplot(1, 1, 1)
    axis.imshow(data, norm=normalization, cmap=colormap)
    axis.savefig(output, bbox_inches="tight")


@app.command()
def preview(input: Path, output: Path):
    data = downsample_asdf_to(input=input, to=(1000, 1000))
    write_image(data, output)


@app.command()
def thumbnail(input: Path, output: Path):
    data = downsample_asdf_to(input=input, to=(300, 300))
    write_image(data, output)


def command():
    app()


if __name__ == "__main__":
    command()
