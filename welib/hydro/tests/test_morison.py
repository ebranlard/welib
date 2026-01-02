import unittest
import numpy as np    
import os
from welib.hydro.wavekin import *
from welib.hydro.morison import *

class Test(unittest.TestCase):

    def test_morison(self):
        g   = 9.81 # gravity [m/s^2]
        h   = 30.  # water depth [m]
        rho = 1025 # water density
        D   = 6    # monopile diameter [m]
        Cd  = 1    # drag coefficient
        CM  = 2    # inertia coefficient
        a   = 3    # wave peak amplitude [m]
        eps = 0    # phase shift
        T   = 12.  # period [s]
        f   = 1./T
        t   = 3*T/4.

        # Get kinematics
        k    = wavenumber(f, h, g)
        eta  = elevation2d(a, f, k, eps, t)
        z = np.linspace(-h,eta,30)
        u_wave, a_wave = kinematics2d(a, f, k, eps, h, t, z, Wheeler=True, eta=eta)

        # Relative motion
        u_rel = u_wave# - u_struct
        a_rel = a_wave# - a_struct

        # Wave loads with wheeler

        p_tot = inline_load(u_rel, a_wave, a_rel, D, rho, Cd, CM=CM)[0]
#         p_tot = inline_load(u, du, D, CD  , CM  , rho)

        np.testing.assert_almost_equal(p_tot[-1], 60539.6349674)



if __name__ == '__main__':
    MyDir=os.path.dirname(__file__)
    unittest.main()
