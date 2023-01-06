
Description
===========

``romancal.assign_wcs`` is the first step run on an image, after ``romancal.ramp_fitting``.
It associates a World Coordinate System (WCS) object with each science exposure.
Note that no fitting is performed in this step; it only creates a WCS object that
transforms positions in the ``detector`` frame to positions in the ``world``
coordinate frame (ICRS) based on the telescope pointing and reference files provided by CRDS.
The constructed WCS object can be accessed as an attribute of the ``meta`` object
when the file is opened as a data model. The forward direction of the transforms is
from detector to world coordinates and the input positions are 0-based.

``romancal.assign_wcs`` uses `GWCS <https://github.com/spacetelescope/gwcs>`__ -
a package for managing the World Coordinate System of astronomical data.
It expects to find the basic WCS keywords in the
``model.meta.wcsinfo`` structure. Distortions are stored in reference files in the
`ASDF <http://asdf-standard.readthedocs.org/en/latest/>`__  format.

``assign_wcs`` retrieves reference files from CRDS and creates a pipeline of transforms from
input frame ``detector`` to the telescope frame ``v2v3`` [1]_. This part of the WCS pipeline may include
intermediate coordinate frames. The basic WCS keywords are used to create
the transform from frame ``v2v3`` to frame ``world``.

Note: in earlier builds of the pipeline the distortions are not available.

Basic WCS keywords and the transform from ``v2v3`` to ``world``
---------------------------------------------------------------

The following attributes in ``meta.wcsinfo`` are used to
define the transform from ``v2v3`` to ``world``:

``RA_REF``, ``DEC_REF`` - a fiducial point on the sky, ICRS, [deg]

``V2_REF``, ``V3_REF`` - a point in the V2V3 system which maps to ``RA_REF``, ``DEC_REF``, [arcsec]

``ROLL_REF`` - local roll angle associated with each aperture, [deg]

``RADESYS`` - standard coordinate system [ICRS]

These quantities are used to create a 3D Euler angle rotation between the V2V3 spherical system,
associated with the telescope, and a standard celestial system.

Using the WCS interactively
---------------------------

Once a science file is opened as a `DataModel` the WCS can be accessed as an attribute
of the meta object. Calling it as a function with detector positions as inputs returns the
corresponding world coordinates:

.. doctest-skip::

  >>> from roman_datamodels import datamodels as rdm
  >>> image = rdm.open('roman_assign_wcs.asdf')
  >>> ra, dec = image.meta.wcs(x, y)

The WCS provides access to intermediate coordinate frames
and transforms between any two frames in the WCS pipeline in forward or
backward direction:

.. doctest-skip::

  >>> image.meta.wcs.available_frames
      ['detector', 'v2v3', 'world']
  >>> v2world = image.meta.wcs.get_transform('v2v3', 'world')
  >>> ra, dec = v2world(v2, v3)
  >>> x1, y1 = image.meta.wcs.invert(ra, dec)

There are methods which allow the result of evaluating the WCS object
to be an ``astropy.SkyCoord`` objext (as opposed to numbers) which allows
further transformation of coordinates to different coordinate frames.

.. rubric:: Footnotes

.. [1] V2V3 is a frame defined by the two perpendicular axes that lay along the primary's mirror plane.
        For completeness, V1 is also part of the telescope frame system, being the axis perpendicular
        to the primary mirror (i.e. along the telecope's optical axis).
