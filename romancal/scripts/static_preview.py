import argparse
from pathlib import Path

import asdf
import numpy


def command():
    try:
        from stpreview.downsample import downsample_asdf_to
        from stpreview.image import (
            north_pole_angle,
            percentile_normalization,
            write_image,
        )
    except (ImportError, ModuleNotFoundError) as err:
        raise ImportError(
            'SDP requirements not installed; do `pip install "romancal[sdp]"`'
        ) from err

    parser = argparse.ArgumentParser()
    parser.add_argument("INPUT", type=Path, help="path to ASDF file with 2D image data")
    parser.add_argument("OUTPUT", type=Path, help="path to output image file")

    subparsers = parser.add_subparsers(dest="subcommand")

    preview_parser = subparsers.add_parser(
        "preview", help="downsample the given ASDF image by the given integer factor"
    )
    preview_parser.add_argument(
        "SHAPE",
        type=int,
        nargs="+",
        default=(1080, 1080),
        help="desired pixel shape of output image",
    )
    preview_parser.add_argument(
        "--no-compass",
        action="store_true",
        help="do not draw a north arrow on the image",
    )

    thumbnail_parser = subparsers.add_parser(
        "thumbnail",
        help="downsample the given ASDF image to the desired shape (the output image may be smaller than the desired shape, if no even factor exists)",
    )
    thumbnail_parser.add_argument(
        "SHAPE",
        type=int,
        nargs="+",
        default=(300, 300),
        help="desired pixel shape of output image",
    )
    thumbnail_parser.add_argument(
        "--compass",
        action="store_true",
        help="draw a north arrow on the image",
    )

    arguments = parser.parse_args()

    output = arguments.OUTPUT
    if output is None:
        output = Path.cwd()
    if output.is_dir():
        output = (
            output / f"{input.stem}.png"
            if arguments.subcommand == "preview"
            else f"{input.stem}_thumb.png"
        )

    with asdf.open(arguments.INPUT) as file:
        model = file["roman"]["meta"]["model_type"]
        if "image" not in model.lower() and "mosaic" not in model.lower():
            raise NotImplementedError(f'"{model}" model not supported')
        wcs = file["roman"]["meta"]["wcs"]

    data = downsample_asdf_to(input=input, shape=arguments.SHAPE, func=numpy.nanmean)

    write_image(
        data,
        output,
        shape=arguments.SHAPE,
        normalization=percentile_normalization(data, percentile=90),
        colormap="afmhot",
        north_arrow_angle=north_pole_angle(wcs).degree - 90
        if arguments.subcommand == "preview"
        else None,
    )


if __name__ == "__main__":
    command()
