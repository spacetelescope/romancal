import warnings

from roman_datamodels.dqflags import group, pixel

__all__ = ["group", "pixel"]

warnings.warn(
    "romancal.dqflags is deprecated. Please use roman_datamodels.dqflags instead.",
    DeprecationWarning,
    stacklevel=2,
)
