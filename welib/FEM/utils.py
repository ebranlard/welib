import numpy as np

def skew(x):
    x=np.asarray(x).ravel()
    """ Returns the skew symmetric matrix M, such that: cross(x,v) = M v """
    return np.array([[0, -x[2], x[1]],[x[2],0,-x[0]],[-x[1],x[0],0]])


# --------------------------------------------------------------------------------}
# --- Direction cosine matrix 
# --------------------------------------------------------------------------------{
def DCM(P1, P2, main_axis='z'):
    """ Computes directional cosine matrix between two points, where `main_axis` is the direction between the two"""
    return DCM_T(P1,P2, main_axis=main_axis).transpose()

def DCM_T(P1, P2, main_axis='z'):
    """
    Computes transpose of directional cosine matrix, where `main_axis` is the direction between the two points.
    Transforms from element to global coordinates:  xg = DC_T.xe,  Kg = DC_T.Ke.DC_T^t
    NOTE that this is the transpose of what is normally considered the Direction Cosine Matrix  
    """
    Dx = P2[0]-P1[0]
    Dy = P2[1]-P1[1]
    Dz = P2[2]-P1[2]
    Dxy = np.sqrt( Dx**2 + Dy**2 )
    L   = np.sqrt( Dx**2 + Dy**2 + Dz**2)
    if L==0:
        raise Exception('Length is zero')

    R = np.zeros((3,3))

    if main_axis=='z':
        if Dxy == 0.0: # TODO tolerance?
            if Dz < 0:   #x is kept along global x
                R[0, 0] =  1.0
                R[1, 1] = -1.0
                R[2, 2] = -1.0
            else:
                R[0, 0] = 1.0
                R[1, 1] = 1.0
                R[2, 2] = 1.0
        else:
            R[0, 0] =  Dy/Dxy
            R[0, 1] = +Dx*Dz/(L*Dxy)
            R[0, 2] =  Dx/L
            R[1, 0] = -Dx/Dxy
            R[1, 1] = +Dz*Dy/(L*Dxy)
            R[1, 2] =  Dy/L
            R[2, 0] = 0.0
            R[2, 1] = -Dxy/L
            R[2, 2] = +Dz/L
    else:
        raise NotImplementedError()
    return R


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
# --- Rigid transformations 
# --------------------------------------------------------------------------------{
def rigidTransformationMatrix(DOF, refPoint, DOF2Nodes, Nodes) :
    """ 
    Returns a rigid body transformation matrix from nDOF to 6 reference DOF: T_ref (6 x nDOF), such that Uref = T_ref.U_subset
    Typically called to get: 
      - the transformation from the interface points to the TP point
      - the transformation from the bottom nodes to SubDyn origin (0,0,)

    INPUTS:
      DOF(nDOF)   : DOF indices that are used to create the transformation matrix
      RefPoint(3) : Coordinate of the reference point 
      DOF2Nodes   : array(nDOF x 4), where columns are: index, node, number of DOF per node, index of DOF in node
      Nodes       : array(nNodes x 3),  x,y,z coordinates of Nodes
    OUTPUTS:
      T_ref(nDOF,6): matrix that relates the subset of DOFs to the reference point
    """
    T_ref = np.zeros((len(DOF), 6))
    for i,iDOF in enumerate(DOF):
        iNode       = DOF2Nodes[iDOF,1] # node 
        nDOFPerNode = DOF2Nodes[iDOF,2] # number of DOF per node
        iiDOF       = DOF2Nodes[iDOF,3] # dof index for this joint (1-6 for cantilever)
        if iiDOF<1 or iiDOF>6:
            import pdb; pdb.set_trace()
            raise Exception('Node DOF number is not valid. DOF: {} Node: {}'.format(iDOF, iNode))
        if nDOFPerNode!=6:
            raise Exception('Node doesnt have 6 DOFs. DOF: {} Node: {}'.format(iDOF, iNode))
        dx = Nodes[iNode, 0] - refPoint[0]
        dy = Nodes[iNode, 1] - refPoint[1]
        dz = Nodes[iNode, 2] - refPoint[2]
        T_ref[i,: ] = rigidTransformationLine(dx, dy, dz, iiDOF)
    return T_ref



def rigidTransformationLine(dx,dy,dz,iLine):
    """ """
    if   iLine ==1: Line = (1, 0, 0, 0 ,  dz, -dy)
    elif iLine ==2: Line = (0, 1, 0,-dz,  0 ,  dx)
    elif iLine ==3: Line = (0, 0, 1, dy, -dx,  0 )
    elif iLine ==4: Line = (0, 0, 0, 1 ,  0 ,  0 )
    elif iLine ==5: Line = (0, 0, 0, 0 ,  1 ,  0 )
    elif iLine ==6: Line = (0, 0, 0, 0 ,  0 ,  1 )
    return Line

def rigidTransformationTwoPoints(Pj, Pk):
    """ Rigid transformation matrix between DOFs of node j and k where node j is the leader node. """
    T_rigid = np.eye(6)
    T_rigid[ 0, 4 ] =  (Pk[2] - Pj[2])
    T_rigid[ 0, 5 ] = -(Pk[1] - Pj[1])
    T_rigid[ 1, 3 ] = -(Pk[2] - Pj[2])
    T_rigid[ 1, 5 ] =  (Pk[0] - Pj[0])
    T_rigid[ 2, 3 ] =  (Pk[1] - Pj[1])
    T_rigid[ 2, 4 ] = -(Pk[0] - Pj[0])
    return T_rigid


