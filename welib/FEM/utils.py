import numpy as np

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

# !> Computes directional cosine matrix DirCos
# !! Transforms from element to global coordinates:  xg = DC.xe,  Kg = DC.Ke.DC^t
# !! Assumes that the element main direction is along ze.
# !! NOTE that this is the transpose of what is normally considered the Direction Cosine Matrix  
# SUBROUTINE GetDirCos(P1, P2, DirCos, L_out, ErrStat, ErrMsg)
#    REAL(ReKi) ,      INTENT(IN   )  :: P1(3), P2(3)      ! (x,y,z) global positions of two nodes making up an element
#    REAL(FEKi) ,      INTENT(  OUT)  :: DirCos(3, 3)      ! calculated direction cosine matrix
#    REAL(ReKi) ,      INTENT(  OUT)  :: L_out             ! length of element
#    INTEGER(IntKi),   INTENT(  OUT)  :: ErrStat           ! Error status of the operation
#    CHARACTER(*),     INTENT(  OUT)  :: ErrMsg            ! Error message if ErrStat /= ErrID_None
#    REAL(FEKi)                       :: Dx, Dy, Dz, Dxy,L! distances between nodes
#    ErrMsg  = ""
#    ErrStat = ErrID_None
#    
#    Dx=P2(1)-P1(1)
#    Dy=P2(2)-P1(2)
#    Dz=P2(3)-P1(3)
#    Dxy = sqrt( Dx**2 + Dy**2 )
#    L   = sqrt( Dx**2 + Dy**2 + Dz**2)
#    
#    IF ( EqualRealNos(L, 0.0_FEKi) ) THEN
#       ErrMsg = ' Same starting and ending location in the element.'
#       ErrStat = ErrID_Fatal
#       RETURN
#    ENDIF
#    
#    IF ( EqualRealNos(Dxy, 0.0_FEKi) ) THEN 
#       DirCos=0.0_FEKi    ! whole matrix set to 0
#       IF ( Dz < 0) THEN  !x is kept along global x
#          DirCos(1, 1) =  1.0_FEKi
#          DirCos(2, 2) = -1.0_FEKi
#          DirCos(3, 3) = -1.0_FEKi
#       ELSE
#          DirCos(1, 1) = 1.0_ReKi
#          DirCos(2, 2) = 1.0_ReKi
#          DirCos(3, 3) = 1.0_ReKi
#       ENDIF 
#    ELSE
#       DirCos(1, 1) =  Dy/Dxy
#       DirCos(1, 2) = +Dx*Dz/(L*Dxy)
#       DirCos(1, 3) =  Dx/L
#       
#       DirCos(2, 1) = -Dx/Dxy
#       DirCos(2, 2) = +Dz*Dy/(L*Dxy)
#       DirCos(2, 3) =  Dy/L
#      
#       DirCos(3, 1) = 0.0_FEKi
#       DirCos(3, 2) = -Dxy/L
#       DirCos(3, 3) = +Dz/L
#    ENDIF
#    L_out= real(L, ReKi)
# 
# END SUBROUTINE GetDirCos
# !------------------------------------------------------------------------------------------------------
# !> Rigid transformation matrix between DOFs of node j and k where node j is the leader node.
# SUBROUTINE GetRigidTransformation(Pj, Pk, TRigid, ErrStat, ErrMsg)
#    REAL(ReKi),       INTENT(IN   )  :: Pj(3)         ! (x,y,z) positions of leader node
#    REAL(ReKi),       INTENT(IN   )  :: Pk(3)         ! (x,y,z) positions of follower node
#    REAL(ReKi),       INTENT(  OUT)  :: TRigid(6,6)   ! Transformation matrix such that xk = T.xj
#    INTEGER(IntKi),   INTENT(  OUT)  :: ErrStat       ! Error status of the operation
#    CHARACTER(*),     INTENT(  OUT)  :: ErrMsg        ! Error message if ErrStat /= ErrID_None
#    ! Local
#    !REAL(ReKi) :: L             ! length of element
#    !REAL(ReKi) :: DirCos(3, 3)  ! direction cosine matrix
#    !REAL(ReKi) :: R0(3,3) 
#    integer(IntKi) :: I
#    ErrStat = ErrID_None
#    ErrMsg  = ""
# 
#    ! --- Formulation using Delta of Global coordinates
#    Trigid=0; do I = 1,6; Trigid(I,I) = 1; enddo
#    Trigid ( 1, 5 ) =  (Pk(3) - Pj(3))
#    Trigid ( 1, 6 ) = -(Pk(2) - Pj(2))
#    Trigid ( 2, 4 ) = -(Pk(3) - Pj(3))
#    Trigid ( 2, 6 ) =  (Pk(1) - Pj(1))
#    Trigid ( 3, 4 ) =  (Pk(2) - Pj(2))
#    Trigid ( 3, 5 ) = -(Pk(1) - Pj(1))
# 
#    ! --- Formulation bty transforming the "local" matrix into a global one
#    !call GetDirCos(Pj, Pk, R0, L, ErrStat, ErrMsg)
#    !TRigid = 0 ; do I = 1,6; TRigid(I,I) = 1; enddo
#    !TRigid (1, 5) =  L
#    !TRigid (2, 4) = -L
#    !TRigid(1:3,4:6) =  matmul( R0 , matmul(TRigid(1:3,4:6), transpose(R0)) )
# 
#    ! --- Formulation using L and Rotation matrix
#    !TRigid = 0; do I = 1,6; TRigid(I,I) = 1; enddo
#    !TRigid ( 1, 5 ) =  L*R0(3,3)
#    !TRigid ( 1, 6 ) = -L*R0(2,3)
#    !TRigid ( 2, 4 ) = -L*R0(3,3)
#    !TRigid ( 2, 6 ) =  L*R0(1,3)
#    !TRigid ( 3, 4 ) =  L*R0(2,3)
#    !TRigid ( 3, 5 ) = -L*R0(1,3)
# END SUBROUTINE GetRigidTransformation
# SUBROUTINE RigidTransformationLine(dx,dy,dz,iLine,Line)
#    real(ReKi),               INTENT(IN)  :: dx,dy,dz
#    integer(IntKi)     ,      INTENT(IN)  :: iLine 
#    Real(ReKi), dimension(6), INTENT(OUT) :: Line
#    SELECT CASE (iLine)
#       CASE (1); Line = (/1.0_ReKi, 0.0_ReKi, 0.0_ReKi, 0.0_ReKi,       dz,      -dy/)
#       CASE (2); Line = (/0.0_ReKi, 1.0_ReKi, 0.0_ReKi,      -dz, 0.0_ReKi,       dx/)
#       CASE (3); Line = (/0.0_ReKi, 0.0_ReKi, 1.0_ReKi,       dy,      -dx, 0.0_ReKi/)
#       CASE (4); Line = (/0.0_ReKi, 0.0_ReKi, 0.0_ReKi, 1.0_ReKi, 0.0_ReKi, 0.0_ReKi/)
#       CASE (5); Line = (/0.0_ReKi, 0.0_ReKi, 0.0_ReKi, 0.0_ReKi, 1.0_ReKi, 0.0_ReKi/)
#       CASE (6); Line = (/0.0_ReKi, 0.0_ReKi, 0.0_ReKi, 0.0_ReKi, 0.0_ReKi, 1.0_ReKi/)
#       CASE DEFAULT
#          Line=-99999999_ReKi
#          print*,'Error in RigidTransformationLine'
#          STOP
# !          ErrStat = ErrID_Fatal
# !          ErrMsg  = 'Error calculating transformation matrix TI '
# !          return
#    END SELECT
# END SUBROUTINE

