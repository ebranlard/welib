import unittest
import os
import numpy as np
from welib.yams.sid import FAST2SID


MyDir=os.path.dirname(__file__)

# --------------------------------------------------------------------------------}
# --- TESTS
# --------------------------------------------------------------------------------{
class Test(unittest.TestCase):

    def test_fast2sid_twr(self):
        np.set_printoptions(linewidth=300, precision=9)
        # --- Read data from NREL5MW tower
        EDFile=os.path.join(MyDir,'./../../../data/NREL5MW/onshore/NREL5MW_ED_Onshore.dat')

        old_settings = np.seterr()
        np.seterr(all='ignore')
        sid, _ = FAST2SID(EDFile, Imodes_twr=[(0,1)])
        np.seterr(**old_settings)

        # --- Generalized mass matrix
        np.testing.assert_almost_equal(np.diag(sid.Mtt), [347460.2316]*3, 5)
        np.testing.assert_almost_equal(np.diag(sid.J.M0)/1e8, np.array([7.198598843e8]*2+[3.474602316e5])/1e8, 5)
        np.testing.assert_almost_equal(        sid.Me.M0[1,1], 61094.66490, 5)
#        np.testing.assert_almost_equal(np.diag(sid.Me.M0), [61094.66490]*2, 5) # <<< TODO TODO returns NaN
#        np.testing.assert_almost_equal(freq[0],     0.891449, 5)
#        np.testing.assert_almost_equal(freq[1],     0.891449, 5)
#        np.testing.assert_almost_equal(freq[-1], 5250.756553, 5)
# 
        np.testing.assert_almost_equal(sid.Mrt[0,1], -13265404.838207997, 5) # -m*zCOG
        np.testing.assert_almost_equal(sid.Mgt[0,0],  104625.69072, 5) # -m*zCOG
        np.testing.assert_almost_equal(sid.Mgt[1,1],  104625.69072, 5) # -m*zCOG
        np.testing.assert_almost_equal(sid.Mgr[0,1], 6449889.716099, 5) # -m*zCOG
        np.testing.assert_almost_equal(sid.Mgr[1,0],-6449889.716099, 5) # -m*zCOG
# 
        # --- C3 mass matrix             3  3 12 12 ie
        np.testing.assert_almost_equal(sid.C3[0, 0, 0, 0, 0], 16063.6792 , 5) # -m*zCOG
        np.testing.assert_almost_equal(sid.C3[0, 0, 0, 6, 0], 7901.009   , 5) # -m*zCOG
        np.testing.assert_almost_equal(sid.C3[1, 1, 1, 1, 0], 17921.95635, 5) # -m*zCOG
        np.testing.assert_almost_equal(sid.C3[1, 1, 5, 1, 0], 22014.56673, 5) # -m*zCOG
        np.testing.assert_almost_equal(sid.C3[2, 2, 2, 2, 0], 17921.95635, 5) # -m*zCOG
        np.testing.assert_almost_equal(sid.C3[2, 2,10,10, 0], 34359.12315, 5) # -m*zCOG
        # --- Term for second order Cr (Mgr) terms and Oe
        np.testing.assert_almost_equal(sid.Kr[2,0,1], -61094.66491, 5)
        np.testing.assert_almost_equal(sid.Kr[2,1,0],  61094.66491, 5)
        # --- Terms useful for 0th order of Gr, and 1st order of J
        np.testing.assert_almost_equal(sid.C4[0,2,0], 6449889.7161, 4)
        np.testing.assert_almost_equal(sid.C4[1,2,1], 6449889.7161, 4)
        # --- Omega terms
        np.testing.assert_almost_equal(sid.Kom[0][1,1],   -61094.664906, 5)
        np.testing.assert_almost_equal(sid.Kom[1][0,0],   -61094.664906, 5)
        np.testing.assert_almost_equal(sid.Kom[2][0,0],   -61094.664906, 5)
        np.testing.assert_almost_equal(sid.Kom[3][0,1],    61094.664906, 5)
        np.testing.assert_almost_equal(sid.Kom[4][0,0],    0, 5)
        np.testing.assert_almost_equal(sid.Kom[5][0,0],    0, 5)
        np.testing.assert_almost_equal(sid.GKg['omxx'][0,0],   77201.43393, 5)
        np.testing.assert_almost_equal(sid.GKg['omyy'][0,0],   77201.43393, 5)
        np.testing.assert_almost_equal(sid.GKg['omzz'][0,0],   0, 5)
        np.testing.assert_almost_equal(sid.GKg['omyz'][0,0],   0, 5)

        #print(sid)
        with open('_OUT_SID_TWR_PY.txt','w') as f:
           f.write(str(sid).replace('-0.000000',' 0.000000'))
        try:
            os.remove('_OUT_SID_TWR_PY.txt')
        except:
            pass




    def test_fast2sid_bld(self):
        np.set_printoptions(linewidth=300, precision=9)
        # --- Read data from NREL5MW tower
        EDFile=os.path.join(MyDir,'./../../../data/NREL5MW/onshore/NREL5MW_ED_Onshore.dat')
        _, sid = FAST2SID(EDFile, Imodes_bld=[0,1])
        with open('_OUT_SID_BLD_PY.txt','w') as f:
           f.write(str(sid).replace('-0.000000',' 0.000000'))
        try:
            os.remove('_OUT_SID_BLD_PY.txt')
        except:
            pass

if __name__=='__main__':
    np.seterr(all='raise')
    #Test().test_fast2sid_bld()
    unittest.main()
