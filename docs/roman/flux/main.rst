Description
===========

:Classes: `romancal.flux.FluxStep`
:Alias: flux

This step will apply the flux scale factor, as defined in the meta of the input
ImageModel. The step will check units and not re-apply the correction, simply
returning the original data.

The ``flux`` step can take:

  * a single 2D input image (in the format of either a string with the full
    path and filename of an ASDF file or a Roman
    Datamodel/:py:class:`~romancal.datamodels.container.ModelContainer`);
  * an association table (in JSON format).


This step takes no other parameters than the standard ``stpipe`` parameters.

The output product are ImageModels, but with the flux scale applied to all the data and variance arrays.


Flux Application
----------------
If the input data is in :math:`DN / second`, the flux scale factor, as
found in the meta ``meta.photometry.conversion_megajanskys``, is simply
multiplied to the ``data`` and ``err`` arrays. The square of the scale factor is
multiplied to all the variance arrays. The resultant units is in :math:`MJy/sr`.

.. math::
   DATA &= DATA * SCALE

   ERR &= ERR * SCALE

   VAR_{rnoise} &= VAR_{rnoise} * SCALE ^ 2

   VAR_{poisson} &= VAR_{poisson} * SCALE ^2

   VAR_{flat} &= VAR_{flat} * SCALE ^2
