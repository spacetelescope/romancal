import logging

import numpy as np
from drizzle import cdrizzle, util

from . import resample_utils

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


class Resampler:
    """
    Combine images using the drizzle algorithm
    """

    def __init__(
        self,
        outsci,
        outwht,
        outcon,
        outwcs,
        pixfrac=1.0,
        kernel="square",
        fillval="INDEF",
        rollover_context=False,
    ):
        """
        Create a new Drizzle output object and set the drizzle parameters.

        Parameters
        ----------

        outsci : TODO

        outwht : TODO

        outcon : TODO

        outwcs : `gwcs.WCS`
            The world coordinate system (WCS) of the resampled image.

        pixfrac : float, optional
            The fraction of a pixel that the pixel flux is confined to. The
            default value of 1 has the pixel flux evenly spread across the image.
            A value of 0.5 confines it to half a pixel in the linear dimension,
            so the flux is confined to a quarter of the pixel area when the square
            kernel is used.

        kernel : str, optional
            The name of the kernel used to combine the inputs. The choice of
            kernel controls the distribution of flux over the kernel. The kernel
            names are: "square", "gaussian", "point", "tophat", "turbo", "lanczos2",
            and "lanczos3". The square kernel is the default.

        fillval : str, optional
            The value a pixel is set to in the output if the input image does
            not overlap it. The default value of INDEF does not set a value.

        rollover_context : TODO
        """

        # Initialize the object fields
        self.uniqid = 0
        self.rollover_context = rollover_context

        self.kernel = kernel

        # Insure that the fillval parameter gets properly interpreted for use with tdriz
        if util.is_blank(str(fillval)):
            self.fillval = "INDEF"
        else:
            self.fillval = str(fillval)

        self.pixfrac = pixfrac

        self.sciext = "SCI"
        self.whtext = "WHT"
        self.conext = "CON"

        self.outwcs = outwcs

        self.outsci = outsci
        self.outwht = outwht
        # FIXME this will need to be re-assigned back to the model
        self.outcon = outcon  # FIXME cast to int32?

        if self.outcon.ndim == 2:
            self.outcon = np.reshape(
                self.outcon, (1, self.outcon.shape[0], self.outcon.shape[1])
            )
        elif self.outcon.ndim != 3:
            raise ValueError(
                f"Drizzle context image has wrong dimensions: {self.outcon.ndim}"
            )

        # Check field values
        if not self.outwcs:
            raise ValueError("Either an existing file or wcs must be supplied")

    def add_image(self, insci, inwcs, inwht, pixmap=None):
        """
        Combine an input image with the output drizzled image.

        Instead of reading the parameters from a fits file, you can set
        them by calling this lower level method. `Add_fits_file` calls
        this method after doing its setup.

        Parameters
        ----------

        insci : array
            A 2d numpy array containing the input image to be drizzled.
            it is an error to not supply an image.

        inwcs : wcs
            The world coordinate system of the input image. This is
            used to convert the pixels to the output coordinate system.

        inwht : array
            A 2d numpy array containing the pixel by pixel weighting.
            Must have the same dimensions as insci.

        pixmap : TODO
        """

        xmin, xmax, ymin, ymax = resample_utils.resample_range(
            insci.shape,
            inwcs.bounding_box,
        )

        self.increment_id()

        if xmax is None or xmax == xmin:
            xmax = insci.shape[1]
        if ymax is None or ymax == ymin:
            ymax = insci.shape[0]

        # Compute what plane of the context image this input would
        # correspond to:
        planeid = int((self.uniqid - 1) / 32)

        # Check if the context image has this many planes
        if self.outcon.ndim == 3:
            nplanes = self.outcon.shape[0]
        elif self.outcon.ndim == 2:
            nplanes = 1
        else:
            nplanes = 0

        if nplanes <= planeid:
            raise IndexError("Not enough planes in drizzle context image")

        # Alias context image to the requested plane if 3d
        if self.outcon.ndim == 3:
            outcon = self.outcon[planeid]

        # Compute the mapping between the input and output pixel coordinates
        # for use in drizzle.cdrizzle.tdriz
        if pixmap is None:
            pixmap = resample_utils.calc_gwcs_pixmap(inwcs, self.outwcs, insci.shape)
        # inwht[np.isnan(pixmap[:,:,0])] = 0.  # FIXME why is this commented?

        log.debug(f"Pixmap shape: {pixmap[:,:,0].shape}")
        log.debug(f"Input Sci shape: {insci.shape}")
        log.debug(f"Output Sci shape: {self.outsci.shape}")

        # Call 'drizzle' to perform image combination
        log.info(f"Drizzling {insci.shape} --> {self.outsci.shape}")

        _vers, nmiss, nskip = cdrizzle.tdriz(
            insci,
            inwht,
            pixmap,
            self.outsci,
            self.outwht,
            outcon,
            uniqid=self.uniqid,
            xmin=xmin,
            xmax=xmax,
            ymin=ymin,
            ymax=ymax,
            pixfrac=self.pixfrac,
            kernel=self.kernel,
            in_units="cps",  # was always cps
            expscale=1.0,  # was expscale set from in_units, always cps
            wtscale=1.0,  # was wt_scl and hard-coded to 1.0
            fillstr=self.fillval,
        )
        log.info(
            f"Results from cdrizzle.tdriz(): \
                '_vers'={_vers}, 'nmiss'={nmiss}, 'nskip'={nskip}"
        )
        return _vers, nmiss, nskip

    def increment_id(self):
        """
        Increment the id count and add a plane to the context image if needed

        Drizzle tracks which input images contribute to the output image
        by setting a bit in the corresponding pixel in the context image.
        The uniqid indicates which bit. So it must be incremented each time
        a new image is added. Each plane in the context image can hold 32 bits,
        so after each 32 images, a new plane is added to the context.
        """
        if self.rollover_context:
            self.uniqid = (self.uniqid + 1) % 32
            return

        # Compute what plane of the context image this input would
        # correspond to:
        planeid = int(self.uniqid / 32)

        # Add a new plane to the context image if planeid overflows

        if self.outcon.shape[0] == planeid:
            plane = np.zeros_like(self.outcon[0])
            plane = plane.reshape((1, plane.shape[0], plane.shape[1]))
            self.outcon = np.concatenate((self.outcon, plane))

        # Increment the id
        self.uniqid += 1
