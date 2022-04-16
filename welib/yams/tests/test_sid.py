import unittest
import os
import numpy as np
from welib.yams.sid import FASTTower2SID, FASTBlade2SID


MyDir=os.path.dirname(__file__)


def writeSID(sid, filename):

    with open(filename, 'w') as f:
        f.write(str(sid).replace('-0.000000',' 0.000000'))

    Delete=True
    if Delete:
        try:
            os.remove(filename)
        except:
            pass

# --------------------------------------------------------------------------------}
# --- TESTS
# --------------------------------------------------------------------------------{
class Test(unittest.TestCase):

    def test_fast2sid_twr_FEM(self):
        # Use FEM to determine SID
        np.set_printoptions(linewidth=300, precision=9)
        # --- Read data from NREL5MW tower
        EDFile=os.path.join(MyDir,'./../../../data/NREL5MW/onshore/NREL5MW_ED_Onshore.dat')

        old_settings = np.seterr()
        np.seterr(all='ignore')
        sid = FASTTower2SID(EDFile, Imodes_twr=[(0,1)], method='FEM')
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
        writeSID(sid, '_OUT_SID_TWR_PY_FEM.txt')

    def test_fast2sid_twr_SF_OF(self):
        Gravity = 9.80665
        EDFile=os.path.join(MyDir,'./../../../data/NREL5MW/onshore/NREL5MW_ED_Onshore.dat')
        sid = FASTTower2SID(EDFile, method='ShapeFunctions', gravity=Gravity, Imodes_twr=[0,1,2,3])
        writeSID(sid, '_OUT_SID_TWR_PY_SHAPEFUNCTIONS.txt')


    def test_fast2sid_twr_SF_SI(self):
        Gravity = 9.80665
        EDFile=os.path.join(MyDir,'./../../../data/NREL5MW/onshore/NREL5MW_ED_Onshore.dat')
        sid = FASTTower2SID(EDFile, method='ShapeIntegral', Imodes_twr=[0,1,2,3], gravity=Gravity)
        writeSID(sid, '_OUT_SID_TWR_PY_SHAPEINTEGRAL.txt')


    def test_fast2sid_bld_FEM(self):
        np.set_printoptions(linewidth=300, precision=9)
        # --- Read data from NREL5MW tower
        EDFile=os.path.join(MyDir,'./../../../data/NREL5MW/onshore/NREL5MW_ED_Onshore.dat')
        sid = FASTBlade2SID(EDFile, Imodes_bld=[0,1], method='FEM')
        writeSID(sid, '_OUT_SID_BLD_PY_FEM.txt')

    def test_fast2sid_bld_SF_OF(self):
        consistency='OpenFAST'
#         consistency=''
        EDFile=os.path.join(MyDir,'./../../../data/NREL5MW/onshore/NREL5MW_ED_Onshore.dat')
        sid = FASTBlade2SID(EDFile, method='ShapeFunctions', Imodes_bld=[0,1,2], startAtRoot=False, consistency=consistency, AdjBlMs=1)
        writeSID(sid, '_OUT_SID_BLD_PY_SHAPEFUNCTIONS.txt')

        np.testing.assert_almost_equal(sid.J.M0 [0,0], 12319399.637394847)
        np.testing.assert_almost_equal(sid.Oe.M0[0,4],-14620.057176601025)
        np.testing.assert_almost_equal(sid.Ct.M0[0,0],2056.3643491060693)
        np.testing.assert_almost_equal(sid.Gr.M0[0,2],-180861.62393339022)
        np.testing.assert_almost_equal(sid.Ge.M0[0,5],-220.06947247279106)
        np.testing.assert_almost_equal(sid.J.M1 [0,2,0],-90430.81196669511)
        np.testing.assert_almost_equal(sid.Oe.M1[0,0,0],1576.2541285535062)
        np.testing.assert_almost_equal(sid.Ct.M1[0,2,0],32.143137178127645)
        np.testing.assert_almost_equal(sid.Gr.M1[0,0,0],40.94259631698493)




    def test_fast2sid_bld_SF_SI(self):
        consistency='OpenFAST'
#         consistency=''
        EDFile=os.path.join(MyDir,'./../../../data/NREL5MW/onshore/NREL5MW_ED_Onshore.dat')
        sid = FASTBlade2SID(EDFile, method='ShapeIntegral', Imodes_bld=[0,1,2], startAtRoot=False, consistency=consistency, AdjBlMs=1)
        writeSID(sid, '_OUT_SID_BLD_PY_SHAPEINTEGRAL.txt')

        np.testing.assert_almost_equal(sid.J.M0 [0,0], 12319399.637394847)
        np.testing.assert_almost_equal(sid.Oe.M0[0,4],-14620.057176601025)
        np.testing.assert_almost_equal(sid.Ct.M0[0,0],2056.3643491060693)
        np.testing.assert_almost_equal(sid.Gr.M0[0,2],-180861.62393339022)
        np.testing.assert_almost_equal(sid.Ge.M0[0,5],-220.06947247279106)
        np.testing.assert_almost_equal(sid.J.M1 [0,2,0],-90430.81196669511)
        np.testing.assert_almost_equal(sid.Oe.M1[0,0,0],1576.2541285535062)
        np.testing.assert_almost_equal(sid.Ct.M1[0,2,0],32.143137178127645)
        np.testing.assert_almost_equal(sid.Gr.M1[0,0,0],40.94259631698493)


    def test_fast2sid_bld_CompareSID(self):
        consistency=''
        EDFile=os.path.join(MyDir,'./../../../data/NREL5MW/onshore/NREL5MW_ED_Onshore.dat')

        sid_SI = FASTBlade2SID(EDFile, method='ShapeIntegral', Imodes_bld=[0,1,2], startAtRoot=False, consistency=consistency)
        sid_SF = FASTBlade2SID(EDFile, method='ShapeFunctions', Imodes_bld=[0,1,2], startAtRoot=False, consistency=consistency)

        np.testing.assert_almost_equal(sid_SI.J.M0 ,sid_SF.J.M0 )
        np.testing.assert_almost_equal(sid_SI.Oe.M0,sid_SF.Oe.M0)
        np.testing.assert_almost_equal(sid_SI.Ct.M0,sid_SF.Ct.M0)
        np.testing.assert_almost_equal(sid_SI.Cr.M0,sid_SF.Cr.M0)
        np.testing.assert_almost_equal(sid_SI.Gr.M0,sid_SF.Gr.M0)
        np.testing.assert_almost_equal(sid_SI.Ge.M0,sid_SF.Ge.M0)

        np.testing.assert_almost_equal(sid_SI.J.M1 ,sid_SF.J.M1 )
        np.testing.assert_almost_equal(sid_SI.Oe.M1,sid_SF.Oe.M1)
        np.testing.assert_almost_equal(sid_SI.Ct.M1,sid_SF.Ct.M1)
        np.testing.assert_almost_equal(sid_SI.Cr.M1,sid_SF.Cr.M1)
        np.testing.assert_almost_equal(sid_SI.Gr.M1,sid_SF.Gr.M1)


if __name__=='__main__':
    #np.seterr(all='raise')
    #Test().test_fast2sid_bld_OF()
    #Test().test_fast2sid_twr_OF()
    unittest.main()
