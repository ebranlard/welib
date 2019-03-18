import unittest
import sys
import os
import numpy as np
MyDir=os.path.dirname(__file__)
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))

from Polar import * 

# --------------------------------------------------------------------------------}
# ---  
# --------------------------------------------------------------------------------{
class TestDynamicStall(unittest.TestCase):
    def assertNaN(self,x):
        self.assertTrue(np.isnan(x))

    def test_prescribed(self):
        #FFA-W3-241 airfoil Dyna Stall
        P=Polar.fromfile(os.path.join(MyDir,'../data/FFA-W3-241-Re12M.dat'),compute_params=True)

        omega       = 12.57
        T           = 2*np.pi/omega
        tau         = 0.08
        alpham      = 20
        dt          = 0.01                   # time step
        # 
        fs_prev = P.f_st_interp(alpham) # init with steady value
        Cl0 = P.cl_interp(alpham) # init with steady value
        Cl_new,fs_prev_new = P.dynaStallOye_DiscreteStep(alpham,tau,fs_prev,dt)

        # Testing that value at t=0 is equal to the steady state cl
        self.assertEqual(Cl_new,Cl0)
        self.assertEqual(fs_prev_new,fs_prev)

        # An increase of alpha from the steady value should have dCl/dt>0
        Cl_new,fs_prev_new = P.dynaStallOye_DiscreteStep(alpham+1,tau,fs_prev,dt)
        self.assertEqual( (Cl_new-Cl0)>0 ,True)
        self.assertEqual( (fs_prev_new-fs_prev)<0 ,True)

        # A decrease of alpha from the steady value should have dCl/dt<0
        Cl_new,fs_prev_new = P.dynaStallOye_DiscreteStep(alpham-1,tau,fs_prev,dt)
        self.assertEqual( (Cl_new-Cl0)<0 ,True)
        self.assertEqual( (fs_prev_new-fs_prev)>0 ,True)

if __name__ == '__main__':
    unittest.main()
