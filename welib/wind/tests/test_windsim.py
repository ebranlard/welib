# --- Common libraries 
import os
import unittest
import numpy as np
# --- Local
from welib.wind.windsim import *

scriptDir=os.path.dirname(__file__)

class Test(unittest.TestCase):
    """ See examples/ for more examples """

    def test_windsim(self):
        U0    = 8      # Wind speed [m/s], for Kaimal spectrum
        I     = 0.14   # Turbulence intensity [-], for Kaimal spectrum
        L     = 340.2  # Length scale [m], for Kaimal spectrum
        dt    = 0.01   # Time step [s]
        tMax  = 100    # Maximum time for time series [s]

        # NOTE: not sure if this might depend on the random seed
        for method in ['sumcos-manual', 'sumcos-irfft', 'sumcos-idft-ifft']:
            t1, u1, f1, S1 = pointTSKaimal(tMax, dt, U0, U0*I, L, method=method, seed=11)

            np.testing.assert_almost_equal(t1[:3], [0.  , 0.01, 0.02], 6)
            np.testing.assert_almost_equal(t1[-1], tMax, 6)

            df = 1/(len(t1)*dt)
            np.testing.assert_almost_equal(f1[:3], np.array([0., 1, 2])*df, 6)
            np.testing.assert_almost_equal(f1[-1], 49.99500049995 , 6)

            np.testing.assert_almost_equal(S1[:3], [6400.64   ,  25.81305 ,  10.470368], 6)

            np.testing.assert_almost_equal(u1[2],  9.05015524, 6)
            np.testing.assert_almost_equal(u1[35], 9.10152196, 6)




if __name__ == '__main__':
    unittest.main()
