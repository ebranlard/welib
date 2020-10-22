import unittest
import numpy as np
from welib.beams.theory import *

# --------------------------------------------------------------------------------}
# --- TEST  
# --------------------------------------------------------------------------------{
class Test(unittest.TestCase):
    def test_beam_theory_bending(self):
        # --- Uniform beam, testing shape functions and frequencies
        freq_ref = [0.81404393,5.10152622 ,14.28442115,27.99176434,46.27239213,69.12294146]
        Umid_ref = [0.33952311,-0.71366583,0.01968759 ,0.70711864 ,0.00085144 ,-0.70710680]
        Vmid_ref = [0.01163054,0.00453142 ,-0.05551999,0.00045037 ,0.09996480 ,0.00003058]
        Kmid_ref = [0.00011938,0.00157253 ,0.00012147 ,-0.00854920,0.00001702 ,0.02111106]
        L  = 100                  ;
        EI = 1.868211939147334e+12;
        m  = 8.828201296825122e+03;
        freq,x,U,V,K = UniformBeamBendingModes('unloaded-clamped-free',EI,m,A=1,L=L,nModes=6)
        #import matplotlib.pyplot as plt
        #plt.figure
        #plt.plot(x,U[0,:])
        #plt.plot(x,U[1,:])
        #plt.show()
        np.testing.assert_almost_equal(freq,freq_ref)
        np.testing.assert_almost_equal(U[:,50],Umid_ref)
        np.testing.assert_almost_equal(V[:,50],Vmid_ref)
        np.testing.assert_almost_equal(K[:,50],Kmid_ref)

    def test_beam_theory_longi(self):
        # --- Longitudinal beam
        freq_ref=[ 12.93048538, 38.79145615, 64.65242691, 90.51339768]
        Umid_ref=[ 0.70710678, -0.70710678,  -0.70710678,  0.70710678]
        L   = 100
        D   = 8
        t   = 0.045
        A   = np.pi*((D/2)**2-(D/2-t)**2)
        E   = 210e9                       # Young modulus [Pa] [N/m^2]
        rho = 7850
        freq,x,U =  UniformBeamLongiModes('unloaded-clamped-free',E,rho,A,L,nModes=4,norm='tip_norm')

        np.testing.assert_almost_equal(freq,freq_ref)
        np.testing.assert_almost_equal(U[:,50],Umid_ref)
        #import matplotlib.pyplot as plt
        #plt.figure()
        #for i in np.arange(4):
        #    plt.plot(x,U[i,:])
        #print(freq)
        #plt.show()
    
    def test_beam_theory_torsion(self):
        # --- Torsion of a uniform beam
        freq_ref=([5.61858268, 16.85574804, 28.09291340, 39.33007876])
        Vmid_ref=[ 0.70710678,-0.70710678, -0.70710678,   0.70710678]

        L   = 100
        G   = 79.3e9                   #% Shear modulus. Steel: 79.3  [Pa] [N/m^2]
        D   = 8
        t   = 0.045
        A   = np.pi * ((D / 2) ** 2 - (D/ 2 - t) ** 2)
        rho = 7850
        Ip  = np.pi / 32 * (D ** 4 - (D - 2 * t) ** 4)
        Kt  = np.pi / 64 * (D ** 4 - (D - 2 * t) ** 4)
        freq,x,V,_ = UniformBeamTorsionModes('unloaded-clamped-free',G,Kt,Ip,rho,A,L)
        #import matplotlib.pyplot as plt
        #plt.figure()
        #for i in np.arange(4):
        #    plt.plot(x,V[i,:])
        #print(freq)
        #plt.show()
        np.testing.assert_almost_equal(freq,freq_ref)
        np.testing.assert_almost_equal(V[:,50],Vmid_ref)

if __name__=='__main__':
    unittest.main()
