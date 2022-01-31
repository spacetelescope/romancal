""" Roman Data Quality Flags

The definitions are documented in the Roman RTD:

[NOTE: Documentation not yet implemented. Fix this URL when completed.]

https://roman-cal-pipeline.readthedocs.io/en/latest/roman/references_general/references_general.html#data-quality-flags


Implementation
-------------

The flags are implemented as "bit flags": Each flag is assigned a bit position
in a byte, or multi-byte word, of memory. If that bit is set, the flag assigned
to that bit is interpreted as being set or active.

The data structure that stores bit flags is just the standard Python `int`,
which provides 32 bits. Bits of an integer are most easily referred to using
the formula `2**bit_number` where `bit_number` is the 0-index bit of interest.
"""
# for this file drop E241 (multiple spaces after :, for readability)
# flake8: noqa: E241
# Pixel-specific flags
pixel = {'GOOD':             0,      # No bits set, all is good
         'DO_NOT_USE':       2**0,   # Bad pixel. Do not use.
         'SATURATED':        2**1,   # Pixel saturated during exposure
         'JUMP_DET':         2**2,   # Jump detected during exposure
         'DROPOUT':          2**3,   # Data lost in transmission
         'RESERVED_1':       2**4,   #
         'PERSISTENCE':      2**5,   # High persistence (was RESERVED_2)
         'AD_FLOOR':         2**6,   # Below A/D floor (0 DN, was RESERVED_3)
         'RESERVED_4':       2**7,   #
         'UNRELIABLE_ERROR': 2**8,   # Uncertainty exceeds quoted error
         'NON_SCIENCE':      2**9,   # Pixel not on science portion of detector
         'DEAD':             2**10,  # Dead pixel
         'HOT':              2**11,  # Hot pixel
         'WARM':             2**12,  # Warm pixel
         'LOW_QE':           2**13,  # Low quantum efficiency
         'TELEGRAPH':        2**15,  # Telegraph pixel
         'NONLINEAR':        2**16,  # Pixel highly nonlinear
         'BAD_REF_PIXEL':    2**17,  # Reference pixel cannot be used
         'NO_FLAT_FIELD':    2**18,  # Flat field cannot be measured
         'NO_GAIN_VALUE':    2**19,  # Gain cannot be measured
         'NO_LIN_CORR':      2**20,  # Linearity correction not available
         'NO_SAT_CHECK':     2**21,  # Saturation check not available
         'UNRELIABLE_BIAS':  2**22,  # Bias variance large
         'UNRELIABLE_DARK':  2**23,  # Dark variance large
         'UNRELIABLE_SLOPE': 2**24,  # Slope variance large (i.e., noisy pixel)
         'UNRELIABLE_FLAT':  2**25,  # Flat variance large
         'RESERVED_5':       2**26,  #
         'RESERVED_6':       2**27,  #
         'UNRELIABLE_RESET': 2**28,  # Sensitive to reset anomaly
         'RESERVED_7':       2**29,  #
         'OTHER_BAD_PIXEL':  2**30,  # A catch-all flag
         'REFERENCE_PIXEL':  2**31,  # Pixel is a reference pixel
         }

# Group-specific flags. Once groups are combined, these flags
# are equivalent to the pixel-specific flags.
group = {'GOOD':       pixel['GOOD'],
         'DO_NOT_USE': pixel['DO_NOT_USE'],
         'SATURATED':  pixel['SATURATED'],
         'JUMP_DET':   pixel['JUMP_DET'],
         'DROPOUT':    pixel['DROPOUT'],
         'AD_FLOOR':   pixel['AD_FLOOR'],
         }
