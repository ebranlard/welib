# --- Common libraries 
import os
import unittest
import numpy as np
from welib.BEM.steadyBEM import *

MyDir=os.path.dirname(__file__)

class Test(unittest.TestCase):
    """ See examples/ for more examples """

    def test_BEM(self):

        old_settings = np.seterr()
        np.seterr(all='raise')

        # Performs simple BEM simulations of the NREL 5MW turbine for one operating condition.
        # --- Read a FAST model to get the main parameters needed
        nB,cone,r,chord,twist,polars,rho,KinVisc = FASTFile2SteadyBEM(os.path.join(MyDir,'../../../data/NREL5MW/Main_Onshore.fst'))

        # --- Run BEM on a set of operating points
        V0        = 5
        Omega     = 7
        pitch=2     #[deg]
        xdot=0      #[m/s]
        u_turb=0    #[m/s]
        BEM=SteadyBEM(Omega,pitch,V0,xdot,u_turb,
                    nB,cone,r,chord,twist,polars,
                    rho=rho,KinVisc=KinVisc,bTIDrag=False,bAIDrag=True)

        np.testing.assert_almost_equal(BEM.Power ,445183.13,1)
        np.testing.assert_almost_equal(BEM.Thrust,140978.66,1)


        np.seterr(**old_settings)

if __name__ == '__main__':
    unittest.main()
