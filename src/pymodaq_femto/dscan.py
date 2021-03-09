from pypret.pnps import CollinearPNPS
import numpy as np



class DSCAN(CollinearPNPS):
    """ Implements the dispersion scan method [Miranda2012a]_ [Miranda2012b]_.
    """
    method = "dscan"
    parameter_name = "insertion"
    parameter_unit = "m"

    def __init__(self, pulse, process, material):
        """ Creates the instance.

        Parameters
        ----------
        pulse : Pulse instance
            The pulse object that defines the simulation grid.
        process : str
            The nonlinear process used in the PNPS method.
        material : str
            The wedge material.
        """
        super().__init__(pulse, process, material=material)

    def mask(self, insertion):
        w = self.ft.w + self.w0
        k = self.material.k(w, unit="om")
        return np.exp(1.0j * k * insertion)
