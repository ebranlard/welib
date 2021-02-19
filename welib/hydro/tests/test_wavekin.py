import unittest
import numpy as np    
import os
from welib.hydro.wavekin import *

MyDir=os.path.dirname(__file__)

class Test(unittest.TestCase):
    def test_wavenumber(self):
        # Test for one period
        T     = 12.                # period [s]
        g     = 9.81               # gravity [m/s^2]
        h     = 30.                # water depth [m]
        f=1/T


        k= wavenumber(f, h, g)
        wavelength = 2*np.pi/k 
        np.testing.assert_almost_equal(wavelength, [177.04211]*2, 5)

        # Test for an array of periods
        T = np.array([12,12]) # period [s]
        f = 1./T

        k= wavenumber(f, h, g)
        wavelength = 2*np.pi/k 
        np.testing.assert_almost_equal(wavelength, [177.04211]*2, 5)

if __name__ == '__main__':
    unittest.main()
