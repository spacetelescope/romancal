from pathlib import Path
from typing import Optional

import asdf
import numpy


def command():
    try:
        import typer
        from stpreview.downsample import downsample_asdf_to
        from stpreview.image import (
            north_pole_angle,
            percentile_normalization,
            write_image,
        )
        from typing_extensions import Annotated
    except (ImportError, ModuleNotFoundError):
        raise ImportError(
            'SDP requirements not installed; do `pip install "romancal[sdp]"`'
        )

    app = typer.Typer()

    @app.command()
    def preview(
        input: Annotated[
            Path, typer.Argument(help="path to ASDF file with 2D image data")
        ],
        output: Annotated[
            Optional[Path], typer.Argument(help="path to output image file")
        ] = None,
        shape: Annotated[
            Optional[tuple[int, int]],
            typer.Argument(help="desired pixel resolution of output image"),
        ] = (1080, 1080),
        compass: Annotated[
            Optional[bool],
            typer.Option(help="whether to draw a north arrow on the image"),
        ] = True,
    ):
        """
        create a preview image with a north arrow overlay indicating orientation
        """

        if output is None:
            output = Path.cwd()
        if output.is_dir():
            output = output / f"{input.stem}.png"

        with asdf.open(input) as file:
            model = file["roman"]["meta"]["model_type"]
            if "image" not in model.lower():
                raise NotImplementedError(f'"{model}" model not supported')
            wcs = file["roman"]["meta"]["wcs"]

        data = downsample_asdf_to(input=input, shape=shape, func=numpy.nanmean)

        write_image(
            data,
            output,
            shape=shape,
            normalization=percentile_normalization(data, percentile=90),
            colormap="afmhot",
            north_arrow_angle=north_pole_angle(wcs).degree - 90,
        )

    @app.command()
    def thumbnail(
        input: Annotated[
            Path, typer.Argument(help="path to ASDF file with 2D image data")
        ],
        output: Annotated[
            Optional[Path], typer.Argument(help="path to output image file")
        ] = None,
        shape: Annotated[
            Optional[tuple[int, int]],
            typer.Argument(help="desired pixel resolution of output image"),
        ] = (300, 300),
        compass: Annotated[
            Optional[bool],
            typer.Option(help="whether to draw a north arrow on the image"),
        ] = False,
    ):
        if output is None:
            output = Path.cwd()
        if output.is_dir():
            output = output / f"{input.stem}_thumb.png"

        with asdf.open(input) as file:
            model = file["roman"]["meta"]["model_type"]
            if "image" not in model.lower():
                raise NotImplementedError(f'"{model}" model not supported')

        data = downsample_asdf_to(input=input, shape=shape, func=numpy.nanmean)

        write_image(
            data,
            output,
            shape=shape,
            normalization=percentile_normalization(data, percentile=90),
            colormap="afmhot",
        )

    app()


if __name__ == "__main__":
    command()
