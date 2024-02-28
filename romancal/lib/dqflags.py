import warnings

from roman_datamodels.dqflags import group, pixel

__all__ = ["pixel", "group"]

warnings.warn(
    "romancal.dqflags is deprecated. Please use roman_datamodels.dqflags instead.",
    DeprecationWarning,
    stacklevel=2,
)
