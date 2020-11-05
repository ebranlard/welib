""" 
Classes and tools to easily set up a FEM model made of beam elements


"""
import numpy as np

# --------------------------------------------------------------------------------}
# ---  
# --------------------------------------------------------------------------------{
def continuousBeam_linearElements(xNodes, m, Iy, Iz=None, A=None, Kv=None, E=None, G=None, phi=None):
    """ """
    from welib.FEM.utils import elementDCMfromBeamNodes
    from welib.FEM.frame3d import frame3dlin_KeMe
    import scipy 

    E  = 214e9                           # Young modulus  

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

        

    # --- Coordinates system / direction cosine of each element
    # Putting "z" along x
    DCM = elementDCMfromBeamNodes(xNodes)

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

    return MM, KK, DCM, Elem2Nodes, Nodes2DOF



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
