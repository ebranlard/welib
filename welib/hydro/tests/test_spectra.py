import unittest
import numpy as np    
import os
from welib.hydro.spectra import *

MyDir=os.path.dirname(__file__)

class Test(unittest.TestCase):
    def test_jonswap(self):
        Hs = 8.1                      # Significant wave height [m]
        Tp = 12.7                     # Peak period [s]
        df   = 1./3600 # Frequency Delta [Hz]
        fMax = 0.5    # Maximum frequency [Hz]
        freq = np.arange(df, fMax, df)
        S = jonswap(freq, Hs, Tp=Tp)

        fmax_ref = 0.0788888
        Smax_ref = 113.86073

        iMax = np.argmax(S)
        np.testing.assert_almost_equal(freq[iMax], fmax_ref, 5 )
        np.testing.assert_almost_equal(S[iMax]   , Smax_ref, 5 )


if __name__ == '__main__':
    unittest.main()
