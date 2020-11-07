import unittest
import os
import scipy
from numpy.linalg import inv

from welib.FEM.fem_beam import * 
from welib.FEM.frame3dlin import frame3dlin_Mcross, frame3dlin_Kg
import welib.weio as weio
from welib.yams.utils import skew 
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

        MM, KK, DCM, Elem2Nodes, Nodes2DOF, Elem2DOF = beam_frame3dlin(xNodes, m, Iy, Iz=Iz, A=A, Kv=Kv, E=E)

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
        MMr, KKr, Tr = applyBC(MM, KK, Elem2Nodes, Nodes2DOF, BC_root=[0,0,0,0,0,0], BC_tip=[1,1,1,1,1,1])
        iStart= 0; 

        # --- Eigenvalues/vectors
        [Q, freq]= eig(KKr, MMr, freq_out=True)
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

        nShapes = Se.shape[1]

        # --- C3 Element mass matrices
        C3 = np.zeros((3,3,12,12,nElem))
        for ie in np.arange(nElem):
            dx= (xNodes[:,ie+1]-xNodes[:,ie]).reshape(3,1)
            le = np.linalg.norm(dx) # element length
            iNode1, iNode2 = Elem2Nodes[ie,:]
            me1 = m[iNode1]*le   # m l = rho * A * l
            me2 = m[iNode2]*le
            c3 = frame3dlin_Mcross(le,me1,me2)
            C3[:,:,:,:,ie]=c3
        
        #                                 3  3 12 12 ie
        np.testing.assert_almost_equal(C3[0, 0, 0, 0, 0], 16063.6792 , 5) # -m*zCOG
        np.testing.assert_almost_equal(C3[0, 0, 0, 6, 0], 7901.009   , 5) # -m*zCOG
        np.testing.assert_almost_equal(C3[1, 1, 1, 1, 0], 17921.95635, 5) # -m*zCOG
        np.testing.assert_almost_equal(C3[1, 1, 5, 1, 0], 22014.56673, 5) # -m*zCOG
        np.testing.assert_almost_equal(C3[2, 2, 2, 2, 0], 17921.95635, 5) # -m*zCOG
        np.testing.assert_almost_equal(C3[2, 2,10,10, 0], 34359.12315, 5) # -m*zCOG

        # --- Term for second order Cr (Mgr) terms and Oe
        # (5.252) S. 233, (6.401) S. 338
        KFr= np.zeros((3,nDOF_tot,nDOF_tot))
        Kr  =np.zeros((3,nShapes,nShapes))
        for ia in np.arange(3):
            for ie in np.arange(nElem):
                lmn= [0,1,2]
                for l in np.arange(3):
                    m= lmn[1];
                    n= lmn[2];
                    IDOF  = Elem2DOF[ie,:]
                    R     = DCM[:,:,ie]
                    RR    = scipy.linalg.block_diag(R,R,R,R)
                    Gamma = DCM[:,:,ie]

                    KFr[ia][np.ix_(IDOF,IDOF)] += (RR.T).dot( -C3[m, n,:,:,ie] + C3[n, m,:,:,ie]).dot(RR) * Gamma[l, ia]
                    lmn= np.roll(lmn,-1) #circshift(lmn, [0 -1]);
            # % (6.483) S. 367
            Kr[ia,:,:]= (Se.T).dot(KFr[ia]).dot(Se)

        np.testing.assert_almost_equal(Kr[2,0,1], -61094.66491, 5)
        np.testing.assert_almost_equal(Kr[2,1,0],  61094.66491, 5)

        # --- Terms useful for 0th order of Gr, and 1st order of J
        #(6.490) S. 368; (6.531) S. 379 or (6.515) S. 375
        C4= np.zeros((3, 3, nShapes))
        for l in np.arange(nShapes):
            for ia in np.arange(3):
                for ib in np.arange(3):
                    C4[ia, ib, l]= -(Sr[:, ia].T).dot(KFr[ib]).dot(Se[:, l]);

        np.testing.assert_almost_equal(C4[0,2,0], 6449889.7161, 4)
        np.testing.assert_almost_equal(C4[1,2,1], 6449889.7161, 4)

        # --- 
        # (5.268) S. 237
        KFom_ab= np.zeros((3,3, nDOF_tot, nDOF_tot)) # = C6
        for ia in np.arange(3):
            for ib in np.arange(3):
                for l in np.arange(3):
                    for m in np.arange(3):
                        for ie in np.arange(nElem):
                            IDOF  = Elem2DOF[ie,:]
                            R     = DCM[:,:,ie]
                            RR    = scipy.linalg.block_diag(R,R,R,R)
                            Gamma = DCM[:,:,ie]
                            if l==m:
                                m_= l+1;
                                if m_>2: m_= 0
                                n_= m_+1;
                                if n_>2: n_= 0
                                Xi= -(C3[m_, m_,:,:,ie]+C3[n_, n_,:,:,ie]) # (5.266) S. 236
                            else:
                                Xi= C3[m, l,:,:,ie];

                            Madd = (RR.T).dot(Xi).dot(RR) * Gamma[l, ia]*Gamma[m, ib]
                            KFom_ab[ia, ib][np.ix_(IDOF,IDOF)] += Madd
#                             if ia==0 and ib==0:
#                                 if ie==0 or ie==1:
#                                     print(Madd)



        np.testing.assert_almost_equal(KFom_ab[0,0][1,1],  -17921.956354, 5)
        np.testing.assert_almost_equal(KFom_ab[0,0][3,1],   22014.566733, 5)
        np.testing.assert_almost_equal(KFom_ab[0,0][7,1],   -6095.064086, 5)
        np.testing.assert_almost_equal(KFom_ab[0,0][6,6],   0.0, 5)

        np.testing.assert_almost_equal(KFom_ab[2,2][0,0],  -17921.956354, 5)
        np.testing.assert_almost_equal(KFom_ab[2,2][4,0],  -22014.566733, 5)
        np.testing.assert_almost_equal(KFom_ab[2,2][6,0],   -6095.064086, 5)


        # --- DOF undisplaced values
        ZF0= np.zeros((nDOF_tot,1))
        for iNode in np.arange(nNodes):
            IDOF=Nodes2DOF[iNode][:3] # translational DOF only
            ZF0[IDOF,0]= xNodes[:,iNode];

        # ---  (5.271) S. 237
        KFom = np.zeros((6,nDOF_tot, nDOF_tot))
        Kom  = np.zeros((6,nShapes,nShapes))
        Kom0 = np.zeros((nShapes, 6))
        Kom0_= np.zeros((Tr.shape[1], 6));
        for i in np.arange(6):
            if i<3:
                KFom[i]= KFom_ab[i, i]
            else:
                a= i-3;
                b= a+1;
                if b>2: b= 0
                KFom[i]= KFom_ab[a, b] + KFom_ab[a, b].T
            Kom[i]= (Se.T).dot(KFom[i]).dot(Se);
            Kom0 [:, i]= (Se.T).dot(KFom[i]).dot(ZF0).ravel()
            Kom0_[:, i]= (Tr.T).dot(KFom[i]).dot(ZF0).ravel()

        np.testing.assert_almost_equal(Kom[0][1,1],   -61094.664906, 5)
        np.testing.assert_almost_equal(Kom[1][0,0],   -61094.664906, 5)
        np.testing.assert_almost_equal(Kom[2][0,0],   -61094.664906, 5)
        np.testing.assert_almost_equal(Kom[3][0,1],    61094.664906, 5)
        np.testing.assert_almost_equal(Kom[4][0,0],    0, 5)
        np.testing.assert_almost_equal(Kom[5][0,0],    0, 5)


        # --- Stiffening terms
        Kinv= Tr.dot(inv(KKr)).dot(Tr.T);

        iMaxDim = np.argmax(np.max(np.abs(xNodes),axis=1)-np.min(np.abs(xNodes),axis=1)) 
        # only axial stiffening
        # 6.330 S. 319

        def geo_stiff_wrap(Tload):
            return beam_frame3dlin_Kg(Tload, xNodes, Elem2Nodes, Elem2DOF, DCM, E, A)

        # Stiffness from tip load
        Fend_ax = np.zeros((nDOF_tot, 1))
        iNode=nNodes-1 # Load node
        DOF=Nodes2DOF[iNode,:]
        Fend_ax[DOF[iMaxDim], 0]= 1 # Unit loads at tip
        
        GKg=dict()
        GKg['Fend'] = (Se.T).dot( geo_stiff_wrap( -Kinv.dot(Fend_ax))                 ).dot(Se)
        GKg['t_ax'] = (Se.T).dot( geo_stiff_wrap( -Kinv.dot(Tr.dot(Ct0_[:, iMaxDim])))).dot(Se)
        GKg['omxx'] = (Se.T).dot( geo_stiff_wrap( -Kinv.dot(Tr.dot(Kom0_[:, 0])))     ).dot(Se) 
        GKg['omyy'] = (Se.T).dot( geo_stiff_wrap( -Kinv.dot(Tr.dot(Kom0_[:, 1])))     ).dot(Se) 
        GKg['omzz'] = (Se.T).dot( geo_stiff_wrap( -Kinv.dot(Tr.dot(Kom0_[:, 2])))     ).dot(Se) 
        GKg['omxy'] = (Se.T).dot( geo_stiff_wrap( -Kinv.dot(Tr.dot(Kom0_[:, 3])))     ).dot(Se) 
        GKg['omxz'] = (Se.T).dot( geo_stiff_wrap( -Kinv.dot(Tr.dot(Kom0_[:, 4])))     ).dot(Se) 
        GKg['omyz'] = (Se.T).dot( geo_stiff_wrap( -Kinv.dot(Tr.dot(Kom0_[:, 5])))     ).dot(Se) 

        np.testing.assert_almost_equal(GKg['omxx'][0,0],   77201.43393, 5)
        np.testing.assert_almost_equal(GKg['omyy'][0,0],   77201.43393, 5)
        np.testing.assert_almost_equal(GKg['omzz'][0,0],   0, 5)
        np.testing.assert_almost_equal(GKg['omyz'][0,0],   0, 5)


        sid= FEMBeam2SID(Mtt, J0, Mtr, Mgt, Mgr, Mgg, KK, Imodes, freq, xNodes, DCM, Se, Kr, Kom0, Kom, C4, GKg)

        #print(sid)
        with open('_OUT_SID_PY.txt','w') as f:
            f.write(str(sid).replace('-0.000000',' 0.000000'))
# 
# 
# 

if __name__=='__main__':
    unittest.main()
