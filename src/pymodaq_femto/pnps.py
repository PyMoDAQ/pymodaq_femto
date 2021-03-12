from pypret.pnps import CollinearPNPS
from pypret.frequencies import wl2om
import numpy as np




class DSCAN(CollinearPNPS):
    """ Implements the dispersion scan method with Fresnel reflection.

    Subclassing d-scan to support fresnel reflection.
    """

    parameter_name = "insertion"
    parameter_unit = "m"
    method = "dscan"

    def __init__(self, pulse, process, material, **kwargs):
        super().__init__(
            pulse, process, material=material, **kwargs
        )

    def _post_init(self):
        super()._post_init()
        w = self.ft.w + self.w0
        # use only the valid range of the Sellmeier equations
        w1, w2 = sorted(wl2om(np.array(self.material._range)))
        valid = (w >= w1) & (w <= w2)
        w = w[valid]
        # remove first and second Taylor order
        k = self.material.k(w, unit="om")
        k0 = self.material.k(self.w0, unit="om")
        k1 = self.material.k(self.w0 + self.ft.dw, unit="om")
        dk = (k1 - k0) / self.ft.dw
        self._k = k - k0 - dk * self.ft.w[valid]
        self._n = self.material.n(w, unit="om")
        self._mask_valid = valid

    def mask(self, insertion):
        if insertion == 0.0:
            return np.ones(self.ft.N, dtype=np.complex128)
        H = np.zeros(self.ft.N, dtype=np.complex128)
        H[self._mask_valid] = np.exp(1.0j * self._k * insertion)
        return H
