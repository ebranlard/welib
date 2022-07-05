import numpy as np
import scipy


def timoshenko_KeMe(L, A, Ixx, Iyy, Jzz, kappa, E, G, rho, R=None, main_axis='z'): 
    Me = timoshenko_Me(L, A, Ixx, Iyy, Jzz, rho, R=R, main_axis=main_axis)
    Ke = timoshenko_Ke(L, A, Ixx, Iyy, Jzz, kappa, E, G, shear=True, R=R, main_axis=main_axis)
    return Ke, Me

def timoshenko_Me(L, A, Ixx, Iyy, Jzz, rho, R=None, main_axis='z'):
    """ Element Mass matrix for Timoshenko beam elements

    NOTE: beam is along z

    INPUTS:
     - E : Young'se (elastic) modulus
     - L : Elemente length
     - A : Cross seection area

    OPTIONAL INPUTeS
        R   : Transformation matrix (3x3) from global coord to element coord: x_e = R.x_g
             if provided, element matrix is provided in global coord

    """
    if main_axis!='z':
        raise NotImplementedError()

    t = rho*A*L
    rx = rho*Ixx
    ry = rho*Iyy
    po = rho*Jzz*L   

    Me = np.zeros((12,12))
       
    Me[ 8,  8] = t/3.
    Me[ 6,  6] = 13.*t/35. + 6.*ry/(5.*L)
    Me[ 7,  7] = 13.*t/35. + 6.*rx/(5.*L)
    Me[11, 11] = po/3.
    Me[ 9,  9] = t*L*L/105. + 2.*L*rx/15.
    Me[10, 10] = t*L*L/105. + 2.*L*ry/15.
    Me[ 1,  3] = -11.*t*L/210. - rx/10.
    Me[ 0,  4] =  11.*t*L/210. + ry/10.
    Me[ 2,  8] = t/6.
    Me[ 4,  6] =  13*t*L/420 - ry/10
    Me[ 3,  7] = -13*t*L/420 + rx/10
    Me[ 5, 11] = po/6
    Me[ 1,  9] =  13*t*L/420 - rx/10
    Me[ 0, 10] = -13*t*L/420 + ry/10
    Me[ 7,  9] =  11*t*L/210 + rx/10
    Me[ 6, 10] = -11*t*L/210 - ry/10
    Me[ 0,  6] =  9*t/70 - 6*ry/(5*L)
    Me[ 1,  7] =  9*t/70 - 6*rx/(5*L)
    Me[ 3,  9] = -L*L*t/140 - rx*L/30 
    Me[ 4, 10] = -L*L*t/140 - ry*L/30
    Me[ 2,  2] = Me[ 8,  8]
    Me[ 0,  0] = Me[ 6,  6]
    Me[ 1,  1] = Me[ 7,  7]
    Me[ 5,  5] = Me[11, 11]
    Me[ 3,  3] = Me[ 9,  9]
    Me[ 4,  4] = Me[10, 10]
    Me[ 3,  1] = Me[ 1,  3]
    Me[ 4,  0] = Me[ 0,  4]
    Me[ 8,  2] = Me[ 2,  8]
    Me[ 6,  4] = Me[ 4,  6]
    Me[ 7,  3] = Me[ 3,  7]
    Me[11,  5] = Me[ 5, 11]
    Me[ 9,  1] = Me[ 1,  9]
    Me[10,  0] = Me[ 0, 10]
    Me[ 9,  7] = Me[ 7,  9]
    Me[10,  6] = Me[ 6, 10]
    Me[ 6,  0] = Me[ 0,  6]
    Me[ 7,  1] = Me[ 1,  7]
    Me[ 9,  3] = Me[ 3,  9]
    Me[10,  4] = Me[ 4, 10]

    if (R is not None):
        RR = scipy.linalg.block_diag(R,R,R,R)
        Me = np.transpose(RR).dot(Me).dot(RR)
        Me = (Me + Me.T)/2 # enforcing symmetry

    return Me


def timoshenko_Ke(L, A, Ixx, Iyy, Jzz, kappa, E, G, shear=True, R=None, main_axis='z'):
    """ Element stiffness matrix for Timoshenko beam elements

    NOTE: Beam directed along z

    INPUTS:
     - E : Young's (elastic) modulus
     - L : Element length
     - A : Cross section area
     - shear: if true, Timoshenko, else Euler-Bernoulli
     - R    : Transformation matrix (3x3) from global coord to element coord: x_e = R.x_g
              if provided, element matrix is provided in global coord
     """
    if main_axis!='z':
        raise NotImplementedError()
    Ke = np.zeros((12,12))
    
    Ax = kappa*A
    Ay = kappa*A
    if (shear):
       Kx = 12*E*Iyy / (G*Ax*L*L)
       Ky = 12*E*Ixx / (G*Ay*L*L)
    else:
       Kx = 0
       Ky = 0
    Ke[ 8 , 8 ] = E*A/L
    Ke[ 6 , 6 ] = 12.0*E*Iyy/( L*L*L*(1.0 + Kx))
    Ke[ 7 , 7 ] = 12.0*E*Ixx/( L*L*L*(1.0 + Ky))
    Ke[11 , 11] = G*Jzz/L
    Ke[9 , 9  ] = (4.0 + Ky                    ) *E*Ixx / ( L*(1.0+Ky )  )
    Ke[10 , 10] = (4.0 + Kx                    ) *E*Iyy / ( L*(1.0+Kx )  )
    Ke[ 1 , 3 ] = -6.*E*Ixx / ( L*L*(1.0+Ky   ))
    Ke[ 0 , 4 ] = 6.*E*Iyy / ( L*L*(1.0+Kx   ))
    Ke[ 3 , 9 ] = (2.0-Ky                      ) *E*Ixx / ( L*(1.0+Ky )  )
    Ke[ 4 , 10] = (2.0-Kx                      ) *E*Iyy / ( L*(1.0+Kx )  )

    Ke[ 2 , 2 ] =  Ke[8 , 8 ]
    Ke[ 0 , 0 ] =  Ke[6 , 6 ]
    Ke[ 1 , 1 ] =  Ke[7 , 7 ]
    Ke[ 5 , 5 ] =  Ke[11, 11]
    Ke[ 3 , 3 ] =  Ke[9, 9]
    Ke[4  , 4 ] =  Ke[10, 10]

    Ke[3  , 1 ] =  Ke[1 , 3 ]
    Ke[4  , 0 ] =  Ke[0 , 4 ]
    Ke[9 , 3  ] =  Ke[3 , 9]
    Ke[10 , 4 ] =  Ke[4 , 10]

    Ke[11 , 5 ] = -Ke[5 , 5 ]
    Ke[9 , 1  ] =  Ke[3 , 1 ]
    Ke[10 , 0 ] =  Ke[4 , 0 ]
    Ke[8  , 2 ] = -Ke[2 , 2 ]
    Ke[6  , 0 ] = -Ke[0 , 0 ]
    Ke[7  , 1 ] = -Ke[1 , 1 ]
    Ke[5  , 11] = -Ke[5 , 5 ]
    Ke[1  , 9 ] =  Ke[3 , 1 ]
    Ke[0  , 10] =  Ke[4 , 0 ]
    Ke[2  , 8 ] = -Ke[2 , 2 ]
    Ke[0  , 6 ] = -Ke[0 , 0 ]
    Ke[1  , 7 ] = -Ke[1, 1  ]
    Ke[10 , 6 ] = -Ke[4, 0  ]
    Ke[9 , 7  ] = -Ke[3 , 1 ]
    Ke[6  , 10] = -Ke[4 , 0 ]
    Ke[7  , 9 ] = -Ke[3 , 1 ]

    Ke[6  , 4 ] = -Ke[4 , 0 ]
    Ke[4  , 6 ] = -Ke[4 , 0 ]
    Ke[7  , 3 ] = -Ke[3 , 1 ]
    Ke[3  , 7 ] = -Ke[3 , 1 ]

    if (R is not None):
        RR = scipy.linalg.block_diag(R,R,R,R)
        #Ke = np.transpose(RR).dot(Ke).dot(RR)
        Ke = np.transpose(RR).dot( ( Ke.dot(RR) ) )
        Ke = (Ke + Ke.T)/2 # enforcing symmetry 
    return Ke


def timoshenko_Fe_g(L, A, rho, g, R=np.eye(3), main_axis='z'):
    """
    Calculates the lumped forces and moments due to gravity on a given element.
    The element has two nodes, with the loads for both elements stored in array F. 
    Indexing of F is:
           (Fx_1, Fy_1, Fz_1, Mx_1, My_1, Mz_1, Fx_2, Fy_2, Fz_2, Mx_2, My_2, Mz_2)
    INPUTS
        L :    Element length
        A :    Cross section area 
        rho:   material density [kg/m^3]
        g:     acceleration of gravity (along z vertical)
        R:     Transformation matrix (3x3) from global coord to element coord: x_e = R.x_g
               if provided, element matrix is provided in global coord
    """
    F = np.zeros(12)
    w = rho*A*g       # weight per unit length
    # lumped forces on both nodes (z component only):
    F[2] = -0.5*L*w 
    F[8] = F[2]
    # lumped moments on node 1 (x and y components only):
    TempCoeff = L*L*w/12
    F[3] = -TempCoeff * R[2,1] # = -L*w*Dy/12
    F[4] =  TempCoeff * R[2,0] # =  L*w*Dx/12
    # lumped moments on node 2: (note the opposite sign of node 1 moment)
    F[9]  = -F[3]
    F[10] = -F[4]
    #F(12) is 0 for g along z alone
    return F

# SUBROUTINE ElemG(A, L, rho, DirCos, F, g)
#    REAL(ReKi), INTENT( IN ) :: A     !< area
#    REAL(ReKi), INTENT( IN ) :: L     !< element length
#    REAL(ReKi), INTENT( IN ) :: rho   !< density
#    REAL(FEKi), INTENT( IN)  :: DirCos(3,3)      !< From element to global: xg = DC.xe,  Kg = DC.Ke.DC^t
#    REAL(ReKi), INTENT( IN ) :: g     !< gravity
#    REAL(FEKi), INTENT( OUT) :: F(12) !< returned loads. positions 1-6 are the loads for node 1 ; 7-12 are loads for node 2.
#    REAL(FEKi) :: TempCoeff
#    REAL(FEKi) :: w            ! weight per unit length
#    
#    F = 0.0_FEKi      ! initialize whole array to zero, then set the non-zero portions
#    w = rho*A*g       ! weight per unit length
#    
#    ! lumped forces on both nodes (z component only):
#    F(3) = -0.5_FEKi*L*w 
#    F(9) = F(3)
#           
#    ! lumped moments on node 1 (x and y components only):
#    ! bjj: note that RRD wants factor of 1/12 because of boundary conditions. Our MeshMapping routines use factor of 1/6 (assuming generic/different boundary  
#    !      conditions), so we may have some inconsistent behavior. JMJ suggests using line2 elements for SubDyn's input/output meshes to improve the situation.
#    TempCoeff = L*L*w/12.0_FEKi ! let's not calculate this twice  
#    F(4) = -TempCoeff * DirCos(2,3) ! = -L*w*Dy/12._FEKi   !bjj: DirCos(2,3) = Dy/L
#    F(5) =  TempCoeff * DirCos(1,3) ! =  L*w*Dx/12._FEKi   !bjj: DirCos(1,3) = Dx/L
# 
#       ! lumped moments on node 2: (note the opposite sign of node 1 moment)
#    F(10) = -F(4)
#    F(11) = -F(5)
#    !F(12) is 0 for g along z alone
#    
# END SUBROUTINE ElemG

