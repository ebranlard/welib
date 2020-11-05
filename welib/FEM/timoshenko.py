
# SUBROUTINE ElemM_Beam(A, L, Ixx, Iyy, Jzz, rho, DirCos, M)
#    REAL(ReKi), INTENT( IN) :: A, L, Ixx, Iyy, Jzz, rho
#    REAL(FEKi), INTENT( IN) :: DirCos(3,3) !< From element to global: xg = DC.xe,  Kg = DC.Ke.DC^t
#    REAL(FEKi), INTENT(OUT) :: M(12, 12)
# 
#    REAL(FEKi) :: t, rx, ry, po
#    REAL(FEKi) :: DC(12, 12)
#    
#    t = rho*A*L;
#    rx = rho*Ixx;
#    ry = rho*Iyy;
#    po = rho*Jzz*L;   
# 
#    M(1:12,1:12) = 0.0_FEKi
#       
#    M( 9,  9) = t/3.0_FEKi
#    M( 7,  7) = 13.0_FEKi*t/35.0_FEKi + 6.0_FEKi*ry/(5.0_FEKi*L)
#    M( 8,  8) = 13.0_FEKi*t/35.0_FEKi + 6.0_FEKi*rx/(5.0_FEKi*L)
#    M(12, 12) = po/3.0_FEKi
#    M(10, 10) = t*L*L/105.0_FEKi + 2.0_FEKi*L*rx/15.0_FEKi
#    M(11, 11) = t*L*L/105.0_FEKi + 2.0_FEKi*L*ry/15.0_FEKi
#    M( 2,  4) = -11.0_FEKi*t*L/210.0_FEKi - rx/10.0_FEKi
#    M( 1,  5) =  11.0_FEKi*t*L/210.0_FEKi + ry/10.0_FEKi
#    M( 3,  9) = t/6.0_FEKi
#    M( 5,  7) =  13._FEKi*t*L/420._FEKi - ry/10._FEKi
#    M( 4,  8) = -13._FEKi*t*L/420._FEKi + rx/10._FEKi
#    M( 6, 12) = po/6._FEKi
#    M( 2, 10) =  13._FEKi*t*L/420._FEKi - rx/10._FEKi
#    M( 1, 11) = -13._FEKi*t*L/420._FEKi + ry/10._FEKi
#    M( 8, 10) =  11._FEKi*t*L/210._FEKi + rx/10._FEKi
#    M( 7, 11) = -11._FEKi*t*L/210._FEKi - ry/10._FEKi
#    M( 1,  7) =  9._FEKi*t/70._FEKi - 6._FEKi*ry/(5._FEKi*L)
#    M( 2,  8) =  9._FEKi*t/70._FEKi - 6._FEKi*rx/(5._FEKi*L)
#    M( 4, 10) = -L*L*t/140._FEKi - rx*L/30._FEKi 
#    M( 5, 11) = -L*L*t/140._FEKi - ry*L/30._FEKi
#    
#    M( 3,  3) = M( 9,  9)
#    M( 1,  1) = M( 7,  7)
#    M( 2,  2) = M( 8,  8)
#    M( 6,  6) = M(12, 12)
#    M( 4,  4) = M(10, 10)
#    M( 5,  5) = M(11, 11)
#    M( 4,  2) = M( 2,  4)
#    M( 5,  1) = M( 1,  5)
#    M( 9,  3) = M( 3,  9)
#    M( 7,  5) = M( 5,  7)
#    M( 8,  4) = M( 4,  8)
#    M(12,  6) = M( 6, 12)
#    M(10,  2) = M( 2, 10)
#    M(11,  1) = M( 1, 11)
#    M(10,  8) = M( 8, 10)
#    M(11,  7) = M( 7, 11)
#    M( 7,  1) = M( 1,  7)
#    M( 8,  2) = M( 2,  8)
#    M(10,  4) = M( 4, 10)
#    M(11,  5) = M( 5, 11)
#    
#    DC = 0.0_FEKi
#    DC( 1: 3,  1: 3) = DirCos
#    DC( 4: 6,  4: 6) = DirCos
#    DC( 7: 9,  7: 9) = DirCos
#    DC(10:12, 10:12) = DirCos
#    
#    M = MATMUL( MATMUL(DC, M), TRANSPOSE(DC) ) ! TODO change me if direction cosine is transposed
# 
# END SUBROUTINE ElemM_Beam


# !> Element stiffness matrix for classical beam elements
# !! shear is true  -- non-tapered Timoshenko beam 
# !! shear is false -- non-tapered Euler-Bernoulli beam 
# SUBROUTINE ElemK_Beam(A, L, Ixx, Iyy, Jzz, Shear, kappa, E, G, DirCos, K)
#    REAL(ReKi), INTENT( IN) :: A, L, Ixx, Iyy, Jzz, E, G, kappa
#    REAL(FEKi), INTENT( IN) :: DirCos(3,3) !< From element to global: xg = DC.xe,  Kg = DC.Ke.DC^t
#    LOGICAL   , INTENT( IN) :: Shear
#    REAL(FEKi), INTENT(OUT) :: K(12, 12) 
#    ! Local variables
#    REAL(FEKi)                            :: Ax, Ay, Kx, Ky
#    REAL(FEKi)                            :: DC(12, 12)
#    
#    Ax = kappa*A
#    Ay = kappa*A
#    
#    K(1:12,1:12) = 0.0_FEKi
#    
#    IF (Shear) THEN
#       Kx = 12.0_FEKi*E*Iyy / (G*Ax*L*L)
#       Ky = 12.0_FEKi*E*Ixx / (G*Ay*L*L)
#    ELSE
#       Kx = 0.0_FEKi
#       Ky = 0.0_FEKi
#    ENDIF
#       
#    K( 9,  9) = E*A/L
#    K( 7,  7) = 12.0_FEKi*E*Iyy/( L*L*L*(1.0_FEKi + Kx) )
#    K( 8,  8) = 12.0_FEKi*E*Ixx/( L*L*L*(1.0_FEKi + Ky) )
#    K(12, 12) = G*Jzz/L
#    K(10, 10) = (4.0_FEKi + Ky)*E*Ixx / ( L*(1.0_FEKi+Ky) )  
#    K(11, 11) = (4.0_FEKi + Kx)*E*Iyy / ( L*(1.0_FEKi+Kx) )
#    K( 2,  4) = -6._FEKi*E*Ixx / ( L*L*(1.0_FEKi+Ky) )
#    K( 1,  5) =  6._FEKi*E*Iyy / ( L*L*(1.0_FEKi+Kx) )
#    K( 4, 10) = (2.0_FEKi-Ky)*E*Ixx / ( L*(1.0_FEKi+Ky) )
#    K( 5, 11) = (2.0_FEKi-Kx)*E*Iyy / ( L*(1.0_FEKi+Kx) )
#    
#    K( 3,  3)  = K(9,9)
#    K( 1,  1)  = K(7,7)
#    K( 2,  2)  = K(8,8)
#    K( 6,  6)  = K(12,12)
#    K( 4,  4)  = K(10,10)
#    K(5,5)  = K(11,11)
#    K(4,2)  = K(2,4)
#    K(5,1)  = K(1,5)
#    K(10,4) = K(4,10)
#    K(11,5) = K(5,11)
#    K(12,6)= -K(6,6)
#    K(10,2)=  K(4,2)
#    K(11,1)=  K(5,1)
#    K(9,3) = -K(3,3)
#    K(7,1) = -K(1,1)
#    K(8,2) = -K(2,2)
#    K(6, 12) = -K(6,6)
#    K(2, 10) =  K(4,2)
#    K(1, 11) =  K(5,1)
#    K(3, 9)  = -K(3,3)
#    K(1, 7)  = -K(1,1)
#    K(2, 8)  = -K(2,2)
#    K(11,7) = -K(5,1)
#    K(10,8) = -K(4,2)
#    K(7,11) = -K(5,1)
#    K(8,10) = -K(4,2)
#    K(7,5) = -K(5,1)
#    K(5,7) = -K(5,1)
#    K(8,4) = -K(4,2)
#    K(4,8) = -K(4,2)
#    
#    DC = 0.0_FEKi
#    DC( 1: 3,  1: 3) = DirCos
#    DC( 4: 6,  4: 6) = DirCos
#    DC( 7: 9,  7: 9) = DirCos
#    DC(10:12, 10:12) = DirCos
#    
#    K = MATMUL( MATMUL(DC, K), TRANSPOSE(DC) ) ! TODO: change me if DirCos convention is  transposed
#    
# END SUBROUTINE ElemK_Beam
# 
