import unittest
import numpy as np
from welib.beams.theory import *
import welib.beams.cantilever2d as c2d

class Test(unittest.TestCase):
    def test_sequential_general(self):
        # Compare general algorithm with sequential one
        from scipy.interpolate import CubicSpline

        # Data Initialization
        r0= np.array([0.30, 0.50, 1.00, 1.50, 2.00, 2.50, 2.60, 2.70, 2.80, 2.90, 2.95])
        EI1 = np.array([500.00, 468.75, 390.62, 312.50, 234.38, 167.19, 151.56, 137.50, 121.88, 107.81, 100.00]) * 100
        EI2 = np.array([12000.00, 11250.00, 9375.00, 7500.00, 5625.00, 4012.50, 3637.50, 3300.00, 2925.00, 2587.50, 2400.00]) * 100
        m = np.array([1.70, 1.59, 1.33, 1.06, 0.80, 0.57, 0.52, 0.47, 0.41, 0.37, 0.34])
        beta = np.array([21.32, 14.66, 7.69, 5.17, 3.31, 1.38, 0.90, 0.42, -0.07, 0.03, 0.74])*np.pi/180 # [rad]

        # --- Reinterpolate
        r    = np.linspace(0.3, 3, 40)
        m    = CubicSpline(r0, m)(r)
        beta = CubicSpline(r0, beta)(r)
        EI1  = CubicSpline(r0, EI1)(r)
        EI2  = CubicSpline(r0, EI2)(r)

        n_modes = 5
        freqs_gen, modes_gen  = c2d.compute_modes(n_modes, r, EI1, EI2, m, beta)
        freqs_seq, modes_seq = c2d.compute_modes_seq(r, EI1, EI2, m, beta)

        np.testing.assert_almost_equal(freqs_gen, freqs_seq, 5)
        np.testing.assert_almost_equal(modes_gen, modes_seq, 5)

if __name__=='__main__':
    unittest.main()
