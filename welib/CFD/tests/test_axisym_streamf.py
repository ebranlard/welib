import unittest
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Local 
from welib.CFD.axisym_streamf import *

# --------------------------------------------------------------------------------}
# --- TEST  
# --------------------------------------------------------------------------------{
class Test(unittest.TestCase):

    def test_velocity_from_psi(self):
        # Uniform flow, psi=1/2 U0 r^2
        U0   = 10
        nr   = 16  # number of nodes
        nz   = 11
        zmin = -1.5
        zmax = 1.5
        rmin = 0
        rmax = 3
        r   = np.linspace(rmin,rmax,nr)
        z   = np.linspace(zmin,zmax,nz)
        dz  = (zmax-zmin)/(nz-1)
        dr  = (rmax-rmin)/(nr-1)
        ur_ref  = np.zeros((nz,nr))
        uz_ref  = np.zeros((nz,nr))
        Rcp,Zcp = np.meshgrid(r,z)
        uz_ref[:,:]=U0

        psi = 0.5*U0*Rcp**2
        ur,uz = velocity_from_streamfunction(psi, r, dr, dz)

        np.testing.assert_almost_equal(ur,ur_ref)
        np.testing.assert_almost_equal(uz,uz_ref)

    #def test_streamfunction(self):
    # NOTE: test moved to example


if __name__=='__main__':
    unittest.main()


