from pypret.material import BaseMaterial
import numpy as np


class SellmeierF1(BaseMaterial):
    """ Defines a dispersive material via a specific Sellmeier equation.

        This subclass supports materials with a Sellmeier equation of the
        form::

            n^2(l) - 1 = c1 + c2 * l^2 / (l2 - c3^2) + ...

        This is formula 1 from refractiveindex.info [DispersionFormulas]_.
    """

    def _func(self, x):
        c = self._coefficients
        x2 = x * x
        n2 = np.full_like(x, 1.0 + c[0])
        for i in range(1, len(c) - 1, 2):
            n2 += c[i] * x2 / (x2 - c[i + 1] * c[i + 1])
        n2[n2 < 0] = 0
        return np.sqrt(n2)


class SellmeierF2(BaseMaterial):
    """ Defines a dispersive material via a specific Sellmeier equation.

        This subclass supports materials with a Sellmeier equation of the
        form::

            n^2(l) - 1 = c1 + c2 * l^2 / (l2 - c3) + ...

        This is formula 2 from refractiveindex.info [DispersionFormulas]_.
    """

    def _func(self, x):
        c = self._coefficients
        x2 = x * x
        n2 = np.full_like(x, 1.0 + c[0])
        for i in range(1, c.size - 1, 2):
            n2 += c[i] * x2 / (x2 - c[i + 1])
        n2[n2 < 0] = 0
        return np.sqrt(n2)


class RefractiveIndexDotInfo(BaseMaterial):
    """ Defines a dispersive material via a specific Sellmeier equation.

        This subclass supports materials with a Sellmeier equation of the
        form::

            n^2(l) = c1 + c2 * l^(c3) / (l^2 - c4^(c5)) + c6 * l^(c7) / (l^2 - c8^(c9)) + c10 * l^(c11) + ...

        This is formula 4 from refractiveindex.info [DispersionFormulas]_.
    """

    def _func(self, x):
        c = self._coefficients
        x2 = x * x
        n2 = np.full_like(x, c[0])

        if len(c) > 1:
            n2 += c[1] * x ** c[2] / (x2 - c[3] ** c[4])
        if len(c) > 5:
            n2 += c[5] * x ** c[6] / (x2 - c[7] ** c[8])
        for i in range(9, len(c) - 1, 2):
            n2 += c[i] * x ** c[i + 1]
        n2[n2 < 0] = 0
        return np.sqrt(n2)


# Fused Silica dispersion with extended spectral range
FS = SellmeierF1(
    coefficients=[
        0.0000000,
        0.6961663,
        0.0684043,
        0.4079426,
        0.1162414,
        0.8974794,
        9.8961610,
    ],
    freq_range=[1e-7, 6.7e-6],
    name="FS",
    long_name="Fused silica (fused quartz) extended range",
)

# Air dispersion
Air = SellmeierF1(
    coefficients=[
        0.0000000,
        14926.44e-8,
        19.36e-6,
        41807.57e-8,
        7.434e-3,
        0.0000000,
        0.0000000,
    ],
    freq_range=[1e-7, 1e-4],
    name="Air",
    long_name="Air at 0 degrees C",
)

# BK7
BK7 = SellmeierF2(
    coefficients=[
        0.00000000000,
        1.039612120,
        0.00600069867,
        0.231792344,
        0.02001791440,
        1.010469450,
        103.560653,
    ],
    freq_range=[0.3e-6, 2.5e-6],
    name="BK7",
    long_name="N-BK7 (SCHOTT)",
)

# KDP
KDP = RefractiveIndexDotInfo(
    coefficients=[2.259276, 13.00522, 2, 400, 1, 0.01008956, 0, 0.0129426, 1],
    freq_range=[0.2138e-6, 1.529e-6],
    name="KDP",
    long_name="Potassium dihydrogen phosphate",
)

# ADP
ADP = RefractiveIndexDotInfo(
    coefficients=[2.302842, 15.102464, 2, 400, 1, 0.011125165, 0, 0.01325366, 1],
    freq_range=[0.2138e-6, 1.529e-6],
    name="ADP",
    long_name="Ammonium dihydrogen phosphate",
)
