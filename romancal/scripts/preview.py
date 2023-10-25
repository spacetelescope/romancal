from pathlib import Path

import typer
from stpreview.__main__ import percentile_normalization, write_image
from stpreview.downsample import downsample_asdf_to

app = typer.Typer()


@app.command()
def preview(input: Path, output: Path):
    data = downsample_asdf_to(input=input, shape=(1000, 1000))
    write_image(
        data,
        output,
        shape=(1080, 1080),
        normalization=percentile_normalization(),
        colormap="afmhot",
    )


@app.command()
def thumbnail(input: Path, output: Path):
    data = downsample_asdf_to(input=input, shape=(300, 300))
    write_image(
        data,
        output,
        shape=(300, 300),
        normalization=percentile_normalization(),
        colormap="afmhot",
    )


def command():
    app()


if __name__ == "__main__":
    command()
