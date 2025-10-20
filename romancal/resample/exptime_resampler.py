import numpy as np
from drizzle.resample import Drizzle
from roman_datamodels import dqflags
from stcal.resample.utils import (
    build_driz_weight,
)


class ExptimeResampler:
    def __init__(self, out_wcs, out_shape, good_bits, kernel):
        self.out_wcs = out_wcs
        self.out_shape = out_shape
        self.good_bits = good_bits

        self.out_img = np.zeros(out_shape, dtype="f4")
        self.out_wht = np.zeros(out_shape, dtype="f4")
        self.exptime_total = np.zeros(out_shape, dtype="f4")

        self._driz = Drizzle(
            kernel=kernel,
            fillval=0,
            out_img=self.out_img,
            out_wht=self.out_wht,
            disable_ctx=True,
        )

    def add_image(self, model, pixmap, xmin, xmax, ymin, ymax):
        data = np.full(model["data"].shape, model["effective_exposure_time"])

        # create a unit weight map for all the input pixels with science data
        inwht = build_driz_weight(
            {
                "data": model["data"],
                "dq": model["dq"],
            },
            weight_type=None,
            good_bits=self.good_bits,
            flag_name_map=dqflags.pixel,
        )

        self._driz.add_image(
            data,
            pixmap=pixmap,
            exptime=1.0,
            iscale=1.0,
            pixel_scale_ratio=1.0,
            weight_map=inwht,
            wht_scale=1.0,
            pixfrac=1.0,
            in_units="cps",
            xmin=xmin,
            xmax=xmax,
            ymin=ymin,
            ymax=ymax,
        )

        self.exptime_total += self.out_img

        # zero out arrays for next image
        self.out_img[:] = 0
        self.out_wht[:] = 0

    def finalize(self):
        return self.exptime_total
