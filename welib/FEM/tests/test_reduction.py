import unittest
import os
import numpy as np
from welib.FEM.reduction import *

MyDir=os.path.dirname(__file__)

class Test(unittest.TestCase):
    def test_MeKe(self):
        # --- Test reduction on a beam2d element matrix with interface at the top (DOF 2 and 3)
        # Without boundary conditions, the Guyan modes should be rigid body modes
        np.set_printoptions(linewidth=500)
        L = 100
        EI = 1868211939147.334
        Maff = L * 8828.201296825122
        KK = EI / (L ** 3) * np.array([[12,6 * L,- 12,6 * L],[6 * L,4 * L ** 2,- 6 * L,2 * L ** 2],[- 12,- 6 * L,12,- 6 * L],[6 * L,2 * L ** 2,- 6 * L,4 * L ** 2]])
        MM = Maff / 420 * np.array([[156,22 * L,54,- 13 * L],[22 * L,4 * L ** 2,13 * L,- 3 * L ** 2],[54,13 * L,156,- 22 * L],[- 13 * L,- 3 * L ** 2,- 22 * L,4 * L ** 2]])
        Mr,Kr,Phi_G,Phi_CB,f_G,f_CB,_,_ = CraigBampton(MM,KK,[2,3], nModesCB=2, fullModesOut=True)
        #print(MM)
        #print(Mr)
        #print(Kr)
        #print(Phi_G)
        #print(Phi_CB)
        #print('f_CB',f_CB)
        #print('f_G',f_G)

        U1 = Phi_G[0::2,0]
        V1 = Phi_G[1::2,0]
        U2 = Phi_G[0::2,1]
        V2 = Phi_G[1::2,1]

        # --- Solve EVA (has rigid body modes...)
        #__,Lambda = eig(Kr,Mr)
        #f= np.sqrt(np.sort(np.diag(Lambda)))/(2*np.pi)
        #print(f)

        np.testing.assert_almost_equal(f_G, [0, 0])
        np.testing.assert_almost_equal(f_CB, [0.817914, 8.058651],5)

        np.testing.assert_almost_equal(U1, [1.0, 1.0],5)
        np.testing.assert_almost_equal(V1, [0.0, 0.0],5)

        np.testing.assert_almost_equal(V2, [1.0, 1.0],5)
        np.testing.assert_almost_equal(U2, [-L, 0.0],5)

        np.testing.assert_almost_equal(Mr[0,0], 882820.12968,5)

        # --- Test reduction on a beam2d element matrix with interface at the top (DOF 2 and 3)
        # Pin at bottom (one rigid body mode)
        KK = EI / (L ** 3) * np.array([[12,6 * L,- 12,6 * L],[6 * L,4 * L ** 2,- 6 * L,2 * L ** 2],[- 12,- 6 * L,12,- 6 * L],[6 * L,2 * L ** 2,- 6 * L,4 * L ** 2]])
        MM = Maff / 420 * np.array([[156,22 * L,54,- 13 * L],[22 * L,4 * L ** 2,13 * L,- 3 * L ** 2],[54,13 * L,156,- 22 * L],[- 13 * L,- 3 * L ** 2,- 22 * L,4 * L ** 2]])

        MM = np.delete(MM, [0], axis=1) # removing columns
        MM = np.delete(MM, [0], axis=0) # removing Lines
        KK = np.delete(KK, [0], axis=1) # removing columns
        KK = np.delete(KK, [0], axis=0) # removing Lines

        Mr,Kr,Phi_G,Phi_CB,f_G,f_CB,_,_ = CraigBampton(MM,KK,[1,2], nModesCB=1, fullModesOut=True)
        np.testing.assert_almost_equal(f_G, [0, 5.304894],4)
        np.testing.assert_almost_equal(f_CB, [4.74484],5)





if __name__=='__main__':
    unittest.main()
