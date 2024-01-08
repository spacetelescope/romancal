Flags for the DQ, PIXELDQ, and GROUPDQ Arrays.

===  ==========    ================  ===========================================
Bit  Value         Name              Description
===  ==========    ================  ===========================================
0    1             DO_NOT_USE        Bad pixel. Do not use.
1    2             SATURATED         Pixel saturated during exposure
2    4             JUMP_DET          Jump detected during exposure
3    8             DROPOUT           Data lost in transmission
4    16            GW_AFFECTED_DATA  Data affected by the GW read window
5    32            PERSISTENCE       High persistence (was RESERVED_2)
6    64            AD_FLOOR          Below A/D floor (0 DN, was RESERVED_3)
7    128           OUTLIER           Detected as outlier in coadded image
8    256           UNRELIABLE_ERROR  Uncertainty exceeds quoted error
9    512           NON_SCIENCE       Pixel not on science portion of detector
10   1024          DEAD              Dead pixel
11   2048          HOT               Hot pixel
12   4096          WARM              Warm pixel
13   8192          LOW_QE            Low quantum efficiency
15   32768         TELEGRAPH         Telegraph pixel
16   65536         NONLINEAR         Pixel highly nonlinear
17   131072        BAD_REF_PIXEL     Reference pixel cannot be used
18   262144        NO_FLAT_FIELD     Flat field cannot be measured
19   524288        NO_GAIN_VALUE     Gain cannot be measured
20   1048576       NO_LIN_CORR       Linearity correction not available
21   2097152       NO_SAT_CHECK      Saturation check not available
22   4194304       UNRELIABLE_BIAS   Bias variance large
23   8388608       UNRELIABLE_DARK   Dark variance large
24   16777216      UNRELIABLE_SLOPE  Slope variance large (i.e., noisy pixel)
25   33554432      UNRELIABLE_FLAT   Flat variance large
26   67108864      RESERVED_5
27   134217728     RESERVED_6
28   268435456     UNRELIABLE_RESET  Sensitive to reset anomaly
29   536870912     RESERVED_7
30   1073741824    OTHER_BAD_PIXEL   A catch-all flag
31   2147483648    REFERENCE_PIXEL   Pixel is a reference pixel
===  ==========    ================  ===========================================
