# --- Common libraries 
import os
import unittest
import numpy as np
from welib.BEM.highthrust import *

scriptDir=os.path.dirname(__file__)

class Test(unittest.TestCase):
    """ See examples/ for more examples """

    def test_Ct_a(self):
        # Different high thrust corr
        chi = 30*np.pi/180
        F= 1
        a = [0, 0.5, 1]
        Ctgl = Ct_a(a, chi=chi, F=F, method='MomentumGlauertSkew')
        Ctbl = Ct_a(a, chi=chi, F=F, method='BladedCorr'         )
        Cth2 = Ct_a(a, chi=chi, F=F, method='HAWC2'              )

        np.testing.assert_almost_equal(Ctgl , [0      , 1.5275252 , 2.3094011] , 4)
        np.testing.assert_almost_equal(Ctbl , [0      , 1.551637  , 2.839714]  , 4)
        np.testing.assert_almost_equal(Cth2 , [6.4e-6 , 1.546214  , 2.246951]  , 4)

    def test_Ct_a_ConvetionCtU0(self):
        #  Convention
        chi = 30*np.pi/180
        F= 1
        a = [0, 0.5, 1]
        Ctgl = Ct_a(a, chi=chi, F=F, method='MomentumGlauertSkew',convention='CtU0_aU0')
        Cth2 = Ct_a(a, chi=chi, F=F, method='HAWC2'              ,convention='CtU0_aU0')

        np.testing.assert_almost_equal(Ctgl , [0      , 1.239314 , 2.070552] , 4)
        np.testing.assert_almost_equal(Cth2 , [4.8e-6 , 1.259743  , 1.808488]  , 4)

    def test_a_k(self):
        # No Hight thrust
        chi = 30*np.pi/180
        k = [-2, -0.5, 0.5, 2] 
        a, c, vroots = a_k(k, chi=chi, method='MomentumGlauertSkewRoots', outputMore=True) #, phi=phi)
        np.testing.assert_almost_equal(a , [0.551094, 0.280555, 0.280555, 0.551094], 4)
        np.testing.assert_almost_equal(vroots[0,:], [ 0.55109364, -0.91861899, -0.91861899,  0.55109364], 4)
        np.testing.assert_almost_equal(vroots[1,:], [ 0.94141495,  0.28055464,  0.28055464,  0.94141495], 4)
        np.testing.assert_almost_equal(vroots[2,:], [ 0.94141495,  0.98569884,  0.98569884,  0.94141495], 4)
        np.testing.assert_almost_equal(vroots[3,:], [ 2.23274313,  0.98569884,  0.98569884,  2.23274313], 4)

        #a0, k0, k0_lim, tan_chi0 = ak_lim(chi)


    def test_Ct_a_numInverse(self):
        chi= 20 * np.pi/180
        fa_Ct= lambda ct: a_Ct(ct, chi=chi, CT=ct, method='HAWC2')
        fCt_a= lambda a:  Ct_a(a, chi=chi, method='HAWC2') # numInverse 

        a = np.array([0, 0.7, 1])
        Ct = fCt_a(a)
        a2 = [fa_Ct(ct) for ct in Ct]

        np.testing.assert_almost_equal(a, a2, 4)

    def test_Ct_a_numInverse(self):
        fCt_a = lambda a: Ct_a(a, method='Glauert')
        fa_Ct = a_Ct_numInverse(fCt_a)
        
        Ct = np.array([0,0.7,1])
        a  = fa_Ct(Ct)
        Ct2 = fCt_a(a)

        np.testing.assert_almost_equal(Ct, Ct2, 4)


        a = np.array([0,0.7,1])
        Ct = Ct_a(a, method='Glauert')
        a2 = a_Ct(Ct, method='Glauert')  # numInverse
        np.testing.assert_almost_equal(a, a2, 4)

        a = np.array([0,0.7,1])
        Ct = Ct_a(a, method='Spera')
        a2 = a_Ct(Ct, method='Spera')  # numInverse
        np.testing.assert_almost_equal(a, a2, 4)

        a = np.array([0,0.7,1])
        Ct = Ct_a(a, method='Buhl')
        a2 = a_Ct(Ct, method='Buhl')  # numInverse
        np.testing.assert_almost_equal(a, a2, 4)



if __name__ == '__main__':
    unittest.main()
