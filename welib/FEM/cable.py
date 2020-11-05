
# SUBROUTINE ElemK_Cable(A, L, E, T0, DirCos, K)
#    REAL(ReKi), INTENT( IN) :: A, L, E
#    REAL(ReKi), INTENT( IN) :: T0 ! Pretension [N]
#    REAL(FEKi), INTENT( IN) :: DirCos(3,3) !< From element to global: xg = DC.xe,  Kg = DC.Ke.DC^t
#    REAL(FEKi), INTENT(OUT) :: K(12, 12) 
#    ! Local variables
#    REAL(FEKi) :: L0, Eps0, EAL0, EE
#    REAL(FEKi) :: DC(12, 12)
# 
#    Eps0 = T0/(E*A)
#    L0   = L/(1+Eps0)  ! "rest length" for which pretension would be 0
#    EAL0 = E*A/L0
#    EE   = EAL0* Eps0/(1+Eps0)
# 
#    K(1:12,1:12)=0.0_FEKi
# 
#    ! Note: only translational DOF involved (1-3, 7-9)
#    K(1,1)= EE
#    K(2,2)= EE
#    K(3,3)= EAL0
# 
#    K(1,7)= -EE
#    K(2,8)= -EE
#    K(3,9)= -EAL0
# 
#    K(7,1)= -EE
#    K(8,2)= -EE
#    K(9,3)= -EAL0
# 
#    K(7,7)= EE
#    K(8,8)= EE
#    K(9,9)= EAL0
# 
# 
#    DC = 0.0_FEKi
#    DC( 1: 3,  1: 3) = DirCos
#    DC( 4: 6,  4: 6) = DirCos
#    DC( 7: 9,  7: 9) = DirCos
#    DC(10:12, 10:12) = DirCos
#    
#    K = MATMUL( MATMUL(DC, K), TRANSPOSE(DC) ) ! TODO: change me if DirCos convention is  transposed
# END SUBROUTINE ElemK_Cable
# !> 
# SUBROUTINE ElemF_Cable(T0, DirCos, F)
#    REAL(ReKi), INTENT( IN ) :: T0          !< Pretension load [N]
#    REAL(FEKi), INTENT( IN)  :: DirCos(3,3) !< From element to global: xg = DC.xe,  Kg = DC.Ke.DC^t
#    REAL(FEKi), INTENT( OUT) :: F(12)       !< returned loads. 1-6 for node 1; 7-12 for node 2.
#    ! Local variables
#    REAL(FEKi) :: DC(12, 12)
# 
#    F(1:12) = 0.0_FEKi  ! init 
#    F(3) = +T0  
#    F(9) = -T0 
# 
#    DC = 0.0_FEKi
#    DC( 1: 3,  1: 3) = DirCos
#    DC( 4: 6,  4: 6) = DirCos
#    DC( 7: 9,  7: 9) = DirCos
#    DC(10:12, 10:12) = DirCos
# 
#    F = MATMUL(DC, F)! TODO: change me if DirCos convention is  transposed
# 
# END SUBROUTINE ElemF_Cable

# !> Element stiffness matrix for pretension cable
# SUBROUTINE ElemM_Cable(A, L, rho, DirCos, M)
#    REAL(ReKi), INTENT( IN) :: A,rho
#    REAL(FEKi), INTENT( IN) :: L
#    REAL(FEKi), INTENT( IN) :: DirCos(3,3) !< From element to global: xg = DC.xe,  Kg = DC.Ke.DC^t
#    REAL(FEKi), INTENT(OUT) :: M(12, 12) 
#    ! Local variables
#    REAL(FEKi) :: DC(12, 12)
#    REAL(FEKi) :: t
# 
#    t = rho*A*L;
# 
#    M(1:12,1:12) = 0.0_FEKi
# 
#    M( 1,  1) = 13._FEKi/35._FEKi * t
#    M( 2,  2) = 13._FEKi/35._FEKi * t
#    M( 3,  3) = t/3.0_FEKi
# 
#    M( 7,  7) = 13._FEKi/35._FEKi * t
#    M( 8,  8) = 13._FEKi/35._FEKi * t
#    M( 9,  9) = t/3.0_FEKi
# 
#    M( 1,  7) =  9._FEKi/70._FEKi * t
#    M( 2,  8) =  9._FEKi/70._FEKi * t
#    M( 3,  9) = t/6.0_FEKi
# 
#    M( 7,  1) =  9._FEKi/70._FEKi * t 
#    M( 8,  2) =  9._FEKi/70._FEKi * t
#    M( 9,  3) = t/6.0_FEKi
#    
#    DC = 0.0_FEKi
#    DC( 1: 3,  1: 3) = DirCos
#    DC( 4: 6,  4: 6) = DirCos
#    DC( 7: 9,  7: 9) = DirCos
#    DC(10:12, 10:12) = DirCos
#    
#    M = MATMUL( MATMUL(DC, M), TRANSPOSE(DC) ) ! TODO: change me if DirCos convention is  transposed
# END SUBROUTINE ElemM_Cable
