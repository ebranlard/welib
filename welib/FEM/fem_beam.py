""" 
Classes and tools to easily set up a FEM model made of beam elements


References:


 [2]  Richard Schwertassek,  Oskar Wallrapp
      "Dynamik Flexibler Mehrkoerpersysteme : Methoden Der Mechanik 



"""
import numpy as np
import scipy
from welib.FEM.utils import skew 

# --------------------------------------------------------------------------------}
# ---  
# --------------------------------------------------------------------------------{
def cbeam_frame3dlin(xNodes, m, Iy, Iz=None, A=None, Kv=None, E=None, G=None, phi=None):
    """
    Assemble a FEM system for a continuous beam using frame3d linear elements
    Elements are assumed to be connected continuously from 1st node to last

    xNodes: (3 x nNodes) position of the nodes.
    m     : (nNodes) linear mass per length

    phi   : rotation of principal axes wrt mean line (tangent) of the beam
    
    """
    from welib.FEM.frame3dlin import frame3dlin_KeMe
    import scipy 

    assert(xNodes.shape[0]==3)

    nNodes   = xNodes.shape[1]
    nElem    = nNodes-1
    nqe      = 12           # Number of DOF per element
    nqk      = int(nqe/2)   # Number of DOF per nodes
    nDOF_tot = nNodes*nqk   # Total number of DOF without constraint (BC)

    # --- Default values
    if Iz is None: Iz=Iy
    if E is None:  E = 211e9       # Young modulus
    if G is None:  G = E/2/(1+0.3) # Young modulus
    if A is None:  A= m*0+100      # Area
    if Kv is None: Kv= m*0+100     # Saint Venant torsion

    # --- Coordinates system / direction cosine of each element
    # Putting "z" along x
    DCM = elementDCMfromBeamNodes(xNodes,phi=phi)

    # --- Distribution of DOFs on nodes and elements
    Nodes2DOF=np.zeros((nNodes,6), dtype=int)
    for i in np.arange(nNodes):
        Nodes2DOF[i,:]=np.arange( i*6, (i+1)*6) 
    Elem2DOF=np.zeros((nElem,12),dtype=int)
    for i in np.arange(nElem):
        Elem2DOF[i,:]=np.concatenate((Nodes2DOF[i,:], Nodes2DOF[i+1,:]))

    Elem2Nodes=np.zeros((nElem,2), dtype=int)
    for i in np.arange(nElem):
        Elem2Nodes[i,:]=(i,i+1)

    # --- Element mass matrices
    Me = np.zeros((12,12,nElem))
    Ke = np.zeros((12,12,nElem))
    for ie in np.arange(nElem):
        dx= (xNodes[:,ie+1]-xNodes[:,ie]).reshape(3,1)
        le = np.linalg.norm(dx) # element length
        iNode1, iNode2 = Elem2Nodes[ie,:]
        me1 = m[iNode1]*le   # m l = rho * A * l
        me2 = m[iNode2]*le
        A1  = A[iNode1]
        A2  = A[iNode2]
        Kv1 = Kv[iNode1]
        Kv2 = Kv[iNode2]
        Iy1 = Iy[iNode1]
        Iy2 = Iy[iNode2]
        Iz1 = Iz[iNode1]
        Iz2 = Iz[iNode2]
        ke,me = frame3dlin_KeMe(E,G,Kv1,Kv2,A1,A2,Iy1,Iy2,Iz1,Iz2,le,me1,me2, R=None)
        #ke,me= frame3dlin_KeMe(me1, me2, le)
        Me[:,:,ie]=me
        Ke[:,:,ie]=ke

    # --- Assembly
    MM = np.zeros((nDOF_tot,nDOF_tot))
    KK = np.zeros((nDOF_tot,nDOF_tot))
    for ie in np.arange(nElem):
        IDOF = Elem2DOF[ie,:]
        R    = DCM[:,:,ie]
        RR   = scipy.linalg.block_diag(R,R,R,R)
        Mg   = (RR.T).dot(Me[:,:,ie]).dot(RR)
        Kg   = (RR.T).dot(Ke[:,:,ie]).dot(RR)
        MM[np.ix_(IDOF,IDOF)] += Mg
        KK[np.ix_(IDOF,IDOF)] += Kg


    return MM, KK, DCM, Elem2Nodes, Nodes2DOF, Elem2DOF

def cbeam_frame3dlin_Kg(Tload, xNodes, Elem2Nodes, Elem2DOF, DCM, E, A, FEMmodel='frame3d_lin'):
    """ 
    Geometric stiffness due a load Tload on all the DOFs
    """
    from welib.FEM.frame3dlin import frame3dlin_Kg # TODO switch between implementation

    nDOF_tot = len(Tload)
    nElem = Elem2Nodes.shape[0]
    Kg= np.zeros((nDOF_tot,nDOF_tot))

    # --- Element mass matrices
    for ie in np.arange(nElem):
        # Going from load in global to load in local
        IDOF = Elem2DOF[ie,:]
        R    = DCM[:,:,ie]
        RR   = scipy.linalg.block_diag(R,R,R,R)
        Te   = RR.dot(Tload[IDOF])
        # Element geometrical stiffness matrix in global 
        dx = (xNodes[:,ie+1]-xNodes[:,ie]).reshape(3,1)
        L  = np.linalg.norm(dx)                         # element length
        iNode1, iNode2 = Elem2Nodes[ie,:]
        A1 = A[iNode1]
        A2 = A[iNode2]
        Kge_gl  = frame3dlin_Kg(E,A1,A2,L,Te[0],Te[6],R=DCM[:,:,ie])
        # Assembly
        Kg[np.ix_(IDOF,IDOF)] += Kge_gl
    return Kg


# --------------------------------------------------------------------------------}
# ---  
# --------------------------------------------------------------------------------{
def applyBC(MM, KK, Elem2Nodes, Nodes2DOF, BC_root=[0,0,0,0,0,0], BC_tip=[1,1,1,1,1,1]):
    """ Apply simply boundary conditions at tip and root

    INPUTS:
        MM, KK: mass matrix, stiffness matrix
        BC: 6-arrays for the BC of each DOF
                "0" = fixed
                "1" = free
            default: cantilever, root clamped and tip free
    OUTPUTS:
        Mr, Kr : (nr x nr) reduced mass and stiffness matrix
        Tr     : (n x nr) reduction matrix such that  Mr = Tr' MM Tr
 
    """
    nDOF_tot = MM.shape[0]
    Tr=np.eye(nDOF_tot)

    # Root BC
    ie=0
    iNode, _ = Elem2Nodes[ie,:]
    IDOF_to_remove = [i for i,iBC in zip(Nodes2DOF[iNode,:], BC_root) if iBC==0]
    Tr = np.delete(Tr, IDOF_to_remove, axis=1) # removing columns
    # Tip BC
    #ie=nElem-1
    _, iNode = Elem2Nodes[-1,:]
    IDOF_to_remove = [i for i,iBC in zip(Nodes2DOF[iNode,:], BC_tip) if iBC==0]
    Tr = np.delete(Tr, IDOF_to_remove, axis=1) # removing columns

    Mr = (Tr.T).dot(MM).dot(Tr)
    Kr = (Tr.T).dot(KK).dot(Tr)
    return Mr, Kr, Tr



# --------------------------------------------------------------------------------}
# ---  
# --------------------------------------------------------------------------------{
def generalizedMassMatrix(xNodes, MM, Se):
    """ 
    Generalized mass matrix from a FEM representation when the structure is undeflected.
    xNodes: Position of the nodes (3 x nNodes)
    MM:   FEM Mass Matrix (nDOF x nDOF)
    Se:   FEM Modes (nDOF x nModes) (elastic modes, e)
    
    """
    dpn=6 # Number of DOF per nodes

    assert(xNodes.shape[0]==3)
    nDOF=MM.shape[0]

    # --- Rigid body modes (t: translation, r:rotation)
    St = np.zeros((nDOF, 3))
    Sr = np.zeros((nDOF, 3))
    for i in np.arange(xNodes.shape[1]):
        R= skew(xNodes[:,i])
        St[i*dpn   : i*dpn+3, :]= np.eye(3)
        Sr[i*dpn   : i*dpn+3, :]= -R
        Sr[i*dpn+3 : i*dpn+6, :]= np.eye(3)
    # Se: Selected modes (e:elastic)

    # --- Generalized mass matrix
    # Rigid body part             # Different Notations:
    Mtt  = (St.T).dot(MM).dot(St) # Mxx, mE
    J0   = (Sr.T).dot(MM).dot(Sr) # Mrr, Mtt, I0
    Mrt  = (Sr.T).dot(MM).dot(St) # Mrt, Mxt, mc0
    # Flexible part
    Mgt  = (Se.T).dot(MM).dot(St) # Mgt, Mgx, Mxg', Ct0
    Mgr  = (Se.T).dot(MM).dot(Sr) # Mgr, Mgt, Mtg', Cr0
    Mgg  = (Se.T).dot(MM).dot(Se) # Mgg, Me
    return Mtt, J0, Mrt, Mgt, Mgr, Mgg, St, Sr


def shapeIntegrals(xNodes, Nodes2DOF, Elem2Nodes, Elem2DOF, DCM, m, Se, Sr, Tr):
    """ 
    Compute main shape integrals from FEM implementation
    (frame3dlin for C3 for now) 

    See [2] for equations and details

    Inspired by a matlab implementation by J. Geilser:
        https://github.com/jgeisler0303/FEMBeam
    
    """
    from welib.FEM.frame3dlin import frame3dlin_Mcross

    # init
    nElem    = Elem2Nodes.shape[0]
    nNodes   = xNodes.shape[1]
    nShapes  = Se.shape[1]
    nDOF_tot = Se.shape[0]

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

    # --- Term for second order Cr (Mgr) terms and Oe
    # (5.252) S. 233, (6.401) S. 338
    KFr= np.zeros((3,nDOF_tot,nDOF_tot))
    Kr  =np.zeros((3,nShapes,nShapes))
    for ia in np.arange(3):
        for ie in np.arange(nElem):
            lmn= [0,1,2]
            for l in np.arange(3):
                m_= lmn[1];
                n_= lmn[2];
                IDOF  = Elem2DOF[ie,:]
                R     = DCM[:,:,ie]
                RR    = scipy.linalg.block_diag(R,R,R,R)
                Gamma = DCM[:,:,ie]

                KFr[ia][np.ix_(IDOF,IDOF)] += (RR.T).dot( -C3[m_, n_,:,:,ie] + C3[n_, m_,:,:,ie]).dot(RR) * Gamma[l, ia]
                lmn= np.roll(lmn,-1) #circshift(lmn, [0 -1]);
        # % (6.483) S. 367
        Kr[ia,:,:]= (Se.T).dot(KFr[ia]).dot(Se)

    # --- Terms useful for 0th order of Gr, and 1st order of J
    #(6.490) S. 368; (6.531) S. 379 or (6.515) S. 375
    C4= np.zeros((3, 3, nShapes))
    for l in np.arange(nShapes):
        for ia in np.arange(3):
            for ib in np.arange(3):
                C4[ia, ib, l]= -(Sr[:, ia].T).dot(KFr[ib]).dot(Se[:, l]);

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

    return C3, Kr, C4, KFom_ab, Kom, Kom0, Kom0_ 


def geometricalStiffening(xNodes, Kinv, Tr, Se, Nodes2DOF, Elem2Nodes, Elem2DOF, DCM, E, A, Kom0_=None, Ct0_=None):
    """ 
    Axial stiffening terms
    See [2] 6.330 S. 319
    
    """
    def geo_stiff_wrap(Tload):
        return cbeam_frame3dlin_Kg(Tload, xNodes, Elem2Nodes, Elem2DOF, DCM, E, A)

    nDOF_tot = Kinv.shape[0]
    nNodes   = Nodes2DOF.shape[0]
    iMaxDim = np.argmax(np.max(np.abs(xNodes),axis=1)-np.min(np.abs(xNodes),axis=1)) 

    # Stiffness from tip load
    Fend_ax = np.zeros((nDOF_tot, 1))
    iNode=nNodes-1 # Load node
    DOF=Nodes2DOF[iNode,:]
    Fend_ax[DOF[iMaxDim], 0]= 1 # Unit loads at tip
    
    # All axial stiffening contributions
    GKg=dict()
    GKg['Fend'] = (Se.T).dot( geo_stiff_wrap( -Kinv.dot(Fend_ax))                 ).dot(Se)
    GKg['t_ax'] = (Se.T).dot( geo_stiff_wrap( -Kinv.dot(Tr.dot(Ct0_[:, iMaxDim])))).dot(Se)
    GKg['omxx'] = (Se.T).dot( geo_stiff_wrap( -Kinv.dot(Tr.dot(Kom0_[:, 0])))     ).dot(Se) 
    GKg['omyy'] = (Se.T).dot( geo_stiff_wrap( -Kinv.dot(Tr.dot(Kom0_[:, 1])))     ).dot(Se) 
    GKg['omzz'] = (Se.T).dot( geo_stiff_wrap( -Kinv.dot(Tr.dot(Kom0_[:, 2])))     ).dot(Se) 
    GKg['omxy'] = (Se.T).dot( geo_stiff_wrap( -Kinv.dot(Tr.dot(Kom0_[:, 3])))     ).dot(Se) 
    GKg['omxz'] = (Se.T).dot( geo_stiff_wrap( -Kinv.dot(Tr.dot(Kom0_[:, 4])))     ).dot(Se) 
    GKg['omyz'] = (Se.T).dot( geo_stiff_wrap( -Kinv.dot(Tr.dot(Kom0_[:, 5])))     ).dot(Se) 

    return GKg


# TODO verify that these are DCM and not the transpose
def elementDCMfromBeamNodes(xNodes, phi=None):
    """ Generate element Direction cosine matricse (DCM) 
    from a set of ordered node coordinates defining a beam mean line

    INPUTS:
        xNodes: 3 x nNodes
        phi (optional): nNodes angles about mean line to rotate the section axes
    OUTPUTS:
        DCM:  3 x 3 x (nNodes-1)
    """
    def null(a, rtol=1e-5):
        u, s, v = np.linalg.svd(a)
        rank = (s > rtol*s[0]).sum()
        return v[rank:].T.copy()

    assert(xNodes.shape[0]==3)
    nElem=xNodes.shape[1]-1
    DCM = np.zeros((3,3,nElem))
    for i in np.arange(nElem):
        dx= (xNodes[:,i+1]-xNodes[:,i]).reshape(3,1)
        le = np.linalg.norm(dx) # element length
        e1 = dx/le # tangent vector
        if i==0:
            e1_last = e1
            e2_last = null(e1.T)[:,0].reshape(3,1) # x,z-> y , y-> -x 
        # normal vector
        de1 = e1 - e1_last
        if np.linalg.norm(de1)<1e-8:
            e2 = e2_last
        else:
            e2 = de1/np.linalg.norm(de1)
        # Rotation about e1
        if phi is not None:
            R  = np.cos(phi[i])*np.eye(3) + np.sin(phi[i])*skew(e1) + (1-np.cos(phi[i]))*e1.dot(e1.T);
            e2 = R.dot(e2)
        # Third vector
        e3=np.cross(e1.ravel(),e2.ravel()).reshape(3,1)
        DCM[:,:,i]= np.column_stack((e1,e2,e3)).T;
        e1_last= e1
        e2_last= e2
    return DCM



# --------------------------------------------------------------------------------}
# --- Mode tools 
# --------------------------------------------------------------------------------{
def maxModeAmplitudes(q, iDOFstart):
    """ Return max magnitude of components of a mode """
    MaxMag=np.array([0.,0.,0.,0.,0.,0.])
    for i in np.arange(6): 
        MaxMag[i] = np.sum(np.abs(q[iDOFstart+i::6]))
    return MaxMag

def normalize_to_last(Q, Imodes, iStart):
    for iimode, imode in enumerate(Imodes):
        mag = maxModeAmplitudes(Q[:,imode], iStart)[:3]
        iMax= np.argmax(mag);
        v_= Q[iStart+iMax::6, imode];
        Q[:, imode]= Q[:, imode]/v_[-1]
    return Q

def orthogonalizeModePair(Q1, Q2, iStart):
    # Find magnitudes to see in which direction the mode is the most
    mag = maxModeAmplitudes(Q1, iStart)[:3]
    idx= np.argsort(mag)[-1::-1]
    k11 = sum(Q1[iStart+idx[0]-1::6]);
    k12 = sum(Q1[iStart+idx[1]-1::6]);
    k21 = sum(Q2[iStart+idx[0]-1::6]);
    k22 = sum(Q2[iStart+idx[1]-1::6]);
    Q1_ = k11*Q1 + k21*Q2
    Q2_ = k12*Q1 + k22*Q2
    return Q1_, Q2_

# --------------------------------------------------------------------------------}
# ---  
# --------------------------------------------------------------------------------{
if __name__=='__main__':
    pass
