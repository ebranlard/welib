import unittest
import numpy as np    
from welib.wt_theory.wakeexpansion import *

class TestExpansion(unittest.TestCase):
    def test_downstreamD(self):
        # Check that analytical solution and numerical solution match when getting downstream distance
        CT=0.8
        fraction = 0.5
        rw0 = wake_expansion_momentum(CT=CT)
        
        expansion = 1 + fraction * (rw0-1)
        xa = downstreamDistanceForGivenExpansion(CT, expansion, model='cylinder', method='analytical')
        xn = downstreamDistanceForGivenExpansion(CT, expansion, model='cylinder', method='interp')
        np.testing.assert_almost_equal(xa, xn, 5)

    def test_methods(self):
        CT = 0.8
        fraction = 0.5
        rw0 = wake_expansion_momentum(CT=CT)
        np.testing.assert_almost_equal(rw0, 1.27201965, 7)

        xb = [0, 1, 20] # xb =r/R
        r = wake_expansion(xb, CT=CT, model='cylinder')
        np.testing.assert_almost_equal(r, [1, 1.17048, 1.27153], 5)


if __name__ == '__main__':
    unittest.main()
