import unittest
import os
import scipy
import numpy as np
from numpy.linalg import inv

from welib.FEM.fem_beam import applyBC, generalizedMassMatrix, shapeIntegrals
from welib.FEM.fem_beam import geometricalStiffening
from welib.FEM.fem_beam import orthogonalizeModePair, normalize_to_last
from welib.FEM.fem_beam import cbeam_assembly_frame3dlin, cbeam_frame3dlin_Kg
from welib.FEM.frame3dlin import frame3dlin_Mcross, frame3dlin_Kg
import welib.weio as weio
from welib.system.eva import eig
from welib.yams.sid import FEMBeam2SID


MyDir=os.path.dirname(__file__)

# --------------------------------------------------------------------------------}
# --- TESTS
# --------------------------------------------------------------------------------{
class Test(unittest.TestCase):
    def test_beam_linear_element(self):
        np.set_printoptions(linewidth=300, precision=9)
        # --- Read data from NREL5MW tower
        TowerHt=87.6;
        TowerBs=0;
        TwrFile=os.path.join(MyDir,'./../../../data/NREL5MW/5MW_Baseline/NRELOffshrBsline5MW_Onshore_ElastoDyn_Tower.dat')
        twr = weio.FASTInputFile(TwrFile).toDataFrame()
        z = twr['HtFract_[-]']*(TowerHt-TowerBs)
        m = twr['TMassDen_[kg/m]']  # mu
        EIy = twr['TwFAStif_[Nm^2]'] 
        EIz = twr['TwSSStif_[Nm^2]']  # TODO actually EIx

        # --- Create Beam FEM model
        # Derived parameters
        A  = m*0+100                         # Area
        Kv = m*0+100                         # Saint Venant torsion
        E  = 214e9                           # Young modulus  [N/m^2]
        Iy = EIy/E                           # Area moment [m^4]
        Iz = EIz/E                           # Area moment [m^4]
        nNodes   = len(z)
        nElem    = nNodes-1
        nqe      = 12           # Number of DOF per element
        nqk      = int(nqe/2)   # Number of DOF per nodes
        nDOF_tot = nNodes*nqk   # Total number of DOF without constraint (BC)
        nDOF     = nDOF_tot-nqk # Total number of DOF with constraint (BC) applied

        # Nodes positions
        xNodes = np.zeros((3,nNodes))
        xNodes[2,:]=z

        MM, KK, xNodes, DCM, Elem2Nodes, Nodes2DOF, Elem2DOF = cbeam_assembly_frame3dlin(xNodes, m, Iy, Iz=Iz, A=A, Kv=Kv, E=E)

        np.testing.assert_almost_equal(MM[0,0], 17921.9563543, 5)
        np.testing.assert_almost_equal(MM[-1,-1], 7590.2188  , 5)
        np.testing.assert_almost_equal(MM[11,11], 30565.98330, 5)
        np.testing.assert_almost_equal(MM[20,20], 26585.67290, 5)
        np.testing.assert_almost_equal(KK[7,7]/1e10, 1.91655, 5)
        np.testing.assert_almost_equal(KK[10,10]/1e11, 4.893305, 5)
        np.testing.assert_almost_equal(KK[11,11]/1e12, 1.87917, 5)

        # --- Constraints/ BC
        MMr, KKr, Tr,_,_ = applyBC(MM, KK, Elem2Nodes, Nodes2DOF, BC_root=[0,0,0,0,0,0], BC_tip=[1,1,1,1,1,1])
        iStart= 0; 

        # --- Eigenvalues/vectors
        [Q, freq]= eig(KKr, MMr, freq_out=True)
        np.testing.assert_almost_equal(freq[0],     0.891449, 5)
        np.testing.assert_almost_equal(freq[1],     0.891449, 5)
        np.testing.assert_almost_equal(freq[-1], 5250.756553, 5)

        # --- Orthogonalization/ normalization of modes
        Imodes=[0,1]
        Q0,Q1 = orthogonalizeModePair(Q[:,0],Q[:,1], iStart)
        Q[:,0],Q[:,1] = Q0, Q1
        Q= normalize_to_last(Q, Imodes, iStart);


        # --- Export Modes
        #U1 = np.concatenate(([0],Q[0::6,0] )) # along x
        #V1 = np.concatenate(([0],Q[4::6,0] )) # theta y
        #U2 = np.concatenate(([0],Q[1::6,1] )) # along y
        #V2 = np.concatenate(([0],Q[3::6,1] )) # theta x
        #M=np.column_stack([z,U1,V2,U2,V2])
        # np.savetxt('out.csv',M)

        # --- Generalized mass matrix
        # Selecting modes
        if len(Imodes)>0:
            Se= Tr.dot(Q[:, Imodes]) # nDOF_tot x nShapes
        else:
            Se= Tr # All

        Mtt, J0, Mrt, Mgt, Mgr, Mgg, St, Sr= generalizedMassMatrix(xNodes, MM, Se)

        Ct0_ = (Tr.T).dot(MM).dot(St) # Mode mass matrix for all modes

        np.testing.assert_almost_equal(np.diag(Mtt), [347460.2316]*3, 5)
        np.testing.assert_almost_equal(Mgg[1,1], 61094.66490, 4)
        #try:
        #    np.testing.assert_almost_equal(np.diag(Mgg), [61094.66490]*2, 4)
        #except:

        np.testing.assert_almost_equal(np.diag(J0)/1e8, np.array([7.198598843e8]*2+[3.474602316e5])/1e8, 4)
        np.testing.assert_almost_equal(Mrt[0,1], -13265404.838207997, 4) # -m*zCOG
        np.testing.assert_almost_equal(Mgt[0,0],  104625.69072, 4) # -m*zCOG
        np.testing.assert_almost_equal(Mgt[1,1],  104625.69072, 4) # -m*zCOG
        np.testing.assert_almost_equal(Mgr[0,1], 6449889.716099, 4) # -m*zCOG
        np.testing.assert_almost_equal(Mgr[1,0],-6449889.716099, 4) # -m*zCOG

        # --- Shape integrals
        C3, Kr, C4, KFom_ab, Kom, Kom0, Kom0_ = shapeIntegrals(xNodes, Nodes2DOF, Elem2Nodes, Elem2DOF, DCM, m, Se, Sr, Tr)

        # --- C3 mass matrix             3  3 12 12 ie
        np.testing.assert_almost_equal(C3[0, 0, 0, 0, 0], 16063.6792 , 5) # -m*zCOG
        np.testing.assert_almost_equal(C3[0, 0, 0, 6, 0], 7901.009   , 5) # -m*zCOG
        np.testing.assert_almost_equal(C3[1, 1, 1, 1, 0], 17921.95635, 5) # -m*zCOG
        np.testing.assert_almost_equal(C3[1, 1, 5, 1, 0], 22014.56673, 5) # -m*zCOG
        np.testing.assert_almost_equal(C3[2, 2, 2, 2, 0], 17921.95635, 5) # -m*zCOG
        np.testing.assert_almost_equal(C3[2, 2,10,10, 0], 34359.12315, 5) # -m*zCOG
        # --- Term for second order Cr (Mgr) terms and Oe
        np.testing.assert_almost_equal(Kr[2,0,1], -61094.66491, 5)
        np.testing.assert_almost_equal(Kr[2,1,0],  61094.66491, 5)
        # --- Terms useful for 0th order of Gr, and 1st order of J
        np.testing.assert_almost_equal(C4[0,2,0], 6449889.7161, 4)
        np.testing.assert_almost_equal(C4[1,2,1], 6449889.7161, 4)
        # --- Omega terms
        np.testing.assert_almost_equal(KFom_ab[0,0][1,1],  -17921.956354, 5)
        np.testing.assert_almost_equal(KFom_ab[0,0][3,1],   22014.566733, 5)
        np.testing.assert_almost_equal(KFom_ab[0,0][7,1],   -6095.064086, 5)
        np.testing.assert_almost_equal(KFom_ab[0,0][6,6],   0.0, 5)
        np.testing.assert_almost_equal(KFom_ab[2,2][0,0],  -17921.956354, 5)
        np.testing.assert_almost_equal(KFom_ab[2,2][4,0],  -22014.566733, 5)
        np.testing.assert_almost_equal(KFom_ab[2,2][6,0],   -6095.064086, 5)
        np.testing.assert_almost_equal(Kom[0][1,1],   -61094.664906, 5)
        np.testing.assert_almost_equal(Kom[1][0,0],   -61094.664906, 5)
        np.testing.assert_almost_equal(Kom[2][0,0],   -61094.664906, 5)
        np.testing.assert_almost_equal(Kom[3][0,1],    61094.664906, 5)
        np.testing.assert_almost_equal(Kom[4][0,0],    0, 5)
        np.testing.assert_almost_equal(Kom[5][0,0],    0, 5)

        # --- Stiffening terms
        Kinv= Tr.dot(inv(KKr)).dot(Tr.T);
        GKg= geometricalStiffening(xNodes, Kinv, Tr, Se, Nodes2DOF, Elem2Nodes, Elem2DOF, DCM, E, A, Kom0_, Ct0_)

        np.testing.assert_almost_equal(GKg['omxx'][0,0],   77201.43393, 5)
        np.testing.assert_almost_equal(GKg['omyy'][0,0],   77201.43393, 5)
        np.testing.assert_almost_equal(GKg['omzz'][0,0],   0, 5)
        np.testing.assert_almost_equal(GKg['omyz'][0,0],   0, 5)

        # --- Convert to SID 
        sid= FEMBeam2SID(Mtt, J0, Mrt, Mgt, Mgr, Mgg, KK, xNodes, DCM, Se, Kr, Kom0, Kom, C4, GKg)

        #print(sid)
        with open('_OUT_SID_PY.txt','w') as f:
            f.write(str(sid).replace('-0.000000',' 0.000000'))

        try:
            os.remove('_OUT_SID_PY.txt')
        except:
            pass

    #def test_fast2sid(self):
    #    from welib.yams.sid import FAST2SID
    #    np.set_printoptions(linewidth=300, precision=9)
    #    # --- Read data from NREL5MW tower
    #    EDFile=os.path.join(MyDir,'')
    #    sid = FAST2SID(EDFile, Imodes_twr=[0,1])



if __name__=='__main__':
    unittest.main()
    #Test.test_fast2sid()
