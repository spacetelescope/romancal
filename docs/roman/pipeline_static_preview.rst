Generating Static Previews
==========================

Roman archiving requires static preview images for viewing and selecting images, with the
following requirements for each ``ImageModel``:

- 1080p x 1080p preview image
- 300p x 300p thumbnail image
- output as PNG files
- 90th percentile linear histogram stretch
- using ``afmhot`` colormap
- overlay indicating orientation

The ``roman_static_preview`` script creates downsampled images from ASDF files containing
an ``ImageModel``, with an optional compass rose overlayed onto the image indicating orientation.

Installation
------------

The requirements for this script are not installed by default as part of ``romancal``; install with
the ``sdp`` extra to include them.

.. code-block:: shell

	pip install "romancal[sdp]"

Usage
-----

``roman_static_preview`` includes two convenience commands, ``preview`` and ``thumbnail``, that set
default options to the static preview requirements.

.. code-block:: shell

	❯ roman_static_preview preview --help
	Usage: roman_static_preview preview [OPTIONS] INPUT [OUTPUT] [SHAPE]...

	  create a preview image with a north arrow overlay indicating orientation

	Arguments:
	  INPUT       path to ASDF file with 2D image data  [required]
	  [OUTPUT]    path to output image file
	  [SHAPE]...  desired pixel resolution of output image  [default: 1080, 1080]

	Options:
	  --compass / --no-compass  whether to draw a north arrow on the image
	                            [default: compass]
	  --help                    Show this message and exit.

.. code-block:: shell

	❯ roman_static_preview thumbnail --help
	Usage: roman_static_preview thumbnail [OPTIONS] INPUT [OUTPUT] [SHAPE]...

	Arguments:
	  INPUT       path to ASDF file with 2D image data  [required]
	  [OUTPUT]    path to output image file
	  [SHAPE]...  desired pixel resolution of output image  [default: 300, 300]

	Options:
	  --compass / --no-compass  whether to draw a north arrow on the image
	                            [default: no-compass]
	  --help                    Show this message and exit.

Examples
--------

.. code-block:: shell

	roman_static_preview preview r0000501001001001001_01101_0001_WFI01_cal.asdf r0000501001001001001_01101_0001_WFI01_cal.png 400 400

.. image:: ../images/r0000501001001001001_01101_0001_WFI01_cal.png
   :alt: preview of Roman imagery, with compass rose showing orientation

.. code-block:: shell

	roman_static_preview thumbnail r0000501001001001001_01101_0001_WFI01_cal.asdf r0000501001001001001_01101_0001_WFI01_cal_thumb.png

.. image:: ../images/r0000501001001001001_01101_0001_WFI01_cal_thumb.png
   :alt: thumbnail of Roman imagery

using ``stpreview`` directly
----------------------------

The ``roman_static_preview`` script is merely a wrapper over ``stpreview to``, which
offers more options for fine-grained control of the output image. ``stpreview`` offers
the ``to`` and ``by`` commands (for resampling ``to`` a desired image shape, or ``by``
a desired factor, respectively). Refer to `the documentation <https://github.com/spacetelescope/stpreview#usage>`_
for usage instructions.
