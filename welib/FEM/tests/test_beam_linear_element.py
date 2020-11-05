import unittest
import os
from welib.FEM.fem_beam import *


MyDir=os.path.dirname(__file__)

def test():
    import welib.weio as weio
    from welib.yams.utils import skew 
    from welib.system.eva import eig


    np.set_printoptions(linewidth=300, precision=9)

    # --- Read data from NREL5MW tower
    TowerHt=87.6;
    TowerBs=0;
    TwrFile=os.path.join(MyDir,'./../../../_data/NREL5MW/data/NREL5MW_ED_Tower_Onshore.dat')
    twr = weio.FASTInputFile(TwrFile).toDataFrame()
    z = twr['HtFract_[-]']*(TowerHt-TowerBs)
    m = twr['TMassDen_[kg/m]']  # mu
    EIy = twr['TwFAStif_[Nm^2]'] 
    EIz = twr['TwSSStif_[Nm^2]']  # TODO actually EIx

    # --- Create Beam FEM model
    # Derived parameters
    A  = m*0+100                         # Area
    Kv = m*0+100                         # Saint Venant torsion
    E  = 214e9                           # Young modulus  
    Iy = EIy/E
    Iz = EIz/E
    nNodes   = len(z)
    nElem    = nNodes-1
    nqe      = 12           # Number of DOF per element
    nqk      = int(nqe/2)   # Number of DOF per nodes
    nDOF_tot = nNodes*nqk   # Total number of DOF without constraint (BC)
    nDOF     = nDOF_tot-nqk # Total number of DOF with constraint (BC) applied

    # Nodes positions
    xNodes = np.zeros((3,nNodes))
    xNodes[2,:]=z

    MM, KK, DCM, Elem2Nodes, Nodes2DOF = continuousBeam_linearElements(xNodes, m, Iy, Iz=Iz, A=A, Kv=Kv, E=E)

    #np.testing.assert_almost_equal(Ke[0,0,0]/1e12   , 2.44292, 5)
    #np.testing.assert_almost_equal(Ke[-1,-1,-1]/1e10, 5.58489, 5)
    #np.testing.assert_almost_equal(Me[0,0,0], 16063.67920)
    #np.testing.assert_almost_equal(Me[-1,-1,-1], 16843.6269446)
    np.testing.assert_almost_equal(MM[0,0], 17921.9563543, 5)
    np.testing.assert_almost_equal(MM[-1,-1], 7590.2188  , 5)
    np.testing.assert_almost_equal(MM[11,11], 30565.98330, 5)
    np.testing.assert_almost_equal(MM[20,20], 26585.67290, 5)
    np.testing.assert_almost_equal(KK[7,7]/1e10, 1.91655, 5)
    np.testing.assert_almost_equal(KK[10,10]/1e11, 4.893305, 5)
    np.testing.assert_almost_equal(KK[11,11]/1e12, 1.87917, 5)

    # --- Constraints/ BC
    Mr, Kr, Tr = applyBC(MM, KK, Elem2Nodes, Nodes2DOF, BC_root=[0,0,0,0,0,0], BC_tip=[1,1,1,1,1,1])
    iStart= 0; 

    # --- Eigenvalues/vectors
    [Q, freq]= eig(Kr, Mr, freq_out=True)
    #E[EF, ef_idx]= sort(EF);
    #EV= V(:, ef_idx);

    np.testing.assert_almost_equal(freq[0],     0.891449, 5)
    np.testing.assert_almost_equal(freq[1],     0.891449, 5)
    np.testing.assert_almost_equal(freq[-1], 5250.756553, 5)

    # --- Orthogonalization/ normalization of modes
    Imodes=[0,1]
    Q[:,0],Q[:,1] = orthogonalizeModePair(Q[:,0],Q[:,1], iStart)
    Q= normalize_to_last(Q, Imodes, iStart);

    #import matplotlib.pyplot as plt
    #fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    #fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    #ax.plot( Q[:,0]   , label='')
    #ax.plot( Q[:,1]   , label='')
    #ax.set_xlabel('')
    #ax.set_ylabel('')
    #ax.legend()
    #ax.tick_params(direction='in')
    #plt.show()


    # --- Guyan Modes (rigid body translation and rotation)
    Sr = np.zeros((nDOF_tot, 3))
    St = np.zeros((nDOF_tot, 3))
    for i in np.arange(nNodes):
        R= skew(xNodes[:,i])
        St[ i*nqk   : i*nqk+3 , :]= np.eye(3)
        Sr[ i*nqk   : i*nqk+3 , :]= -R
        Sr[ i*nqk+3 : i*nqk+6 , :]= np.eye(3)
    # --- Craig Bampton modes
    if len(Imodes)>0:
        Se= Tr.dot(Q[:, Imodes]) # nDOF_tot x nShapes
    else:
        Se= Tr # All
    #print(Se.shape)
    #print(Sr.shape)
    #print(St.shape)
    nShapes = Se.shape[1]

    # --- Reduced mass matrix (Guyan = CB)
    Mtt  = (St.T).dot(MM).dot(St) # Mxx, mE
    J0   = (Sr.T).dot(MM).dot(Sr) # Mrr, Mtt, I0
    Mgg  = (Se.T).dot(MM).dot(Se) # < CB mass matrix, Me
    Mtr  = (Sr.T).dot(MM).dot(St) # Mtr, Mxt, mc0
    Mgt  = (Se.T).dot(MM).dot(St) # Mgt, Mgx, Mxg', Ct0
    Mgr  = (Se.T).dot(MM).dot(Sr) # Mgr, Mgt, Mtg', Cr0
    # CB modes mass matrix for all modes
    Ct0_ = (Tr.T).dot(MM).dot(St) # 
    #print('Ct0_\n',Ct0_)

    #print('Mtt\n' ,Mtt  )
    #print('J0\n'  ,J0  )
    #print('Me\n'  ,Mgg )
    #print('Mtr\n' ,Mtr )
    #print('Mgt\n' ,Mgt )
    #print('Mgr\n',Mgr )
    np.testing.assert_almost_equal(np.diag(Mtt), [347460.2316]*3, 5)
    np.testing.assert_almost_equal(np.diag(Mgg), [61094.66490]*2, 5)
    np.testing.assert_almost_equal(np.diag(J0)/1e8, np.array([7.198598843e8]*2+[3.474602316e5])/1e8, 5)
    np.testing.assert_almost_equal(Mtr[0,1], -13265404.838207997, 5) # -m*zCOG
    np.testing.assert_almost_equal(Mgt[0,0],  104625.69072, 5) # -m*zCOG
    np.testing.assert_almost_equal(Mgt[1,1],  104625.69072, 5) # -m*zCOG
    np.testing.assert_almost_equal(Mgr[0,1], 6449889.716099, 5) # -m*zCOG
    np.testing.assert_almost_equal(Mgr[1,0],-6449889.716099, 5) # -m*zCOG

    
    from welib.yams.sid import FEMBeam2SID

    sid= FEMBeam2SID(Mtt, J0, Mgt, Mgr, Mgg, KK, Imodes, freq, xNodes, DCM, Se, C4=None, Kom=None)

    #print(sid)
    #with open('_OUT_SID_PY.txt','w') as f:
    #    f.write(str(sid).replace('-0.000000',' 0.000000'))



    
# --------------------------------------------------------------------------------}
# --- TESTS
# --------------------------------------------------------------------------------{
class Test(unittest.TestCase):
    def test_beam_linear_element(self):
        test()




if __name__=='__main__':
    unittest.main()
