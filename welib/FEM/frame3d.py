import numpy as np
    

def frame3dlin_KeMe(E,G,Kv1,Kv2,A1,A2,Iy1,Iy2,Iz1,Iz2,L,me1,me2,R=None):
    """
    Linear frame element directed along x
    Values are given at both nodes

    INPUTS
        "1" at node 1, "2" at node 2
        E  : Young's (elastic) modulus
        Gs : Shear modulus. For an isotropic material G = E/2(nu+1) with nu the Poission's ratio
        Kv : Saint-Venant's torsion constant, Polar moment of i
        L  : Element length
        A  : Cross section area
        Ix : Polar  second moment of area,local x-axis. Ix=\iint(y^2+z^2) dy dz [m4]
        Iy : Planar second moment of area,local y-axis. Iy=\iint z^2 dy dz [m4]
        Iz : Planar second moment of area,local z-axis. Iz=\iint y^2 dy dz [m4]
    OPTIONAL INPUTS
        R   : Transformation matrix (3x3) from global coord to element coord: x_e = R.x_g
             if provided, element matrix is provided in global coord
    OUTPUTS
        ke: Element stiffness matrix (12x12)
        me: Element mass matrix (12x12)

    """
    # --- Stifness matrix
    ke = np.array([
        [((A2+A1)*E)/(2*L)  , 0                       , 0                       , 0                    , 0                       , 0                       , -((A2+A1)*E)/(2*L) , 0                       , 0                       , 0                    , 0                       , 0]                       , 
        [0                  , ((6*Iz2+6*Iz1)*E)/L**3  , 0                       , 0                    , 0                       , ((2*Iz2+4*Iz1)*E)/L**2  , 0                  , -((6*Iz2+6*Iz1)*E)/L**3 , 0                       , 0                    , 0                       , ((4*Iz2+2*Iz1)*E)/L**2]  , 
        [0                  , 0                       , ((6*Iy2+6*Iy1)*E)/L**3  , 0                    , -((2*Iy2+4*Iy1)*E)/L**2 , 0                       , 0                  , 0                       , -((6*Iy2+6*Iy1)*E)/L**3 , 0                    , -((4*Iy2+2*Iy1)*E)/L**2 , 0]                       , 
        [0                  , 0                       , 0                       , ((Kv2+Kv1)*G)/(2*L)  , 0                       , 0                       , 0                  , 0                       , 0                       , -((Kv2+Kv1)*G)/(2*L) , 0                       , 0]                       , 
        [0                  , 0                       , -((2*Iy2+4*Iy1)*E)/L**2 , 0                    , ((Iy2+3*Iy1)*E)/L       , 0                       , 0                  , 0                       , ((2*Iy2+4*Iy1)*E)/L**2  , 0                    , ((Iy2+Iy1)*E)/L         , 0]                       , 
        [0                  , ((2*Iz2+4*Iz1)*E)/L**2  , 0                       , 0                    , 0                       , ((Iz2+3*Iz1)*E)/L       , 0                  , -((2*Iz2+4*Iz1)*E)/L**2 , 0                       , 0                    , 0                       , ((Iz2+Iz1)*E)/L]         , 
        [-((A2+A1)*E)/(2*L) , 0                       , 0                       , 0                    , 0                       , 0                       , ((A2+A1)*E)/(2*L)  , 0                       , 0                       , 0                    , 0                       , 0]                       , 
        [0                  , -((6*Iz2+6*Iz1)*E)/L**3 , 0                       , 0                    , 0                       , -((2*Iz2+4*Iz1)*E)/L**2 , 0                  , ((6*Iz2+6*Iz1)*E)/L**3  , 0                       , 0                    , 0                       , -((4*Iz2+2*Iz1)*E)/L**2] , 
        [0                  , 0                       , -((6*Iy2+6*Iy1)*E)/L**3 , 0                    , ((2*Iy2+4*Iy1)*E)/L**2  , 0                       , 0                  , 0                       , ((6*Iy2+6*Iy1)*E)/L**3  , 0                    , ((4*Iy2+2*Iy1)*E)/L**2  , 0]                       , 
        [0                  , 0                       , 0                       , -((Kv2+Kv1)*G)/(2*L) , 0                       , 0                       , 0                  , 0                       , 0                       , ((Kv2+Kv1)*G)/(2*L)  , 0                       , 0]                       , 
        [0                  , 0                       , -((4*Iy2+2*Iy1)*E)/L**2 , 0                    , ((Iy2+Iy1)*E)/L         , 0                       , 0                  , 0                       , ((4*Iy2+2*Iy1)*E)/L**2  , 0                    , ((3*Iy2+Iy1)*E)/L       , 0]                       , 
        [0                  , ((4*Iz2+2*Iz1)*E)/L**2  , 0                       , 0                    , 0                       , ((Iz2+Iz1)*E)/L         , 0                  , -((4*Iz2+2*Iz1)*E)/L**2 , 0                       , 0                    , 0                       , ((3*Iz2+Iz1)*E)/L]
        ])
    # --- Mass matrix
    me = np.array([
        [(me2+3*me1)/12 , 0                      , 0                       , 0              , 0                           , 0                           , (me2+me1)/12   , 0                       , 0                      , 0              , 0                           , 0]                           , 
        [0              , (3*me2+10*me1)/35      , 0                       , 0              , 0                           , (7*L*me2+15*L*me1)/420      , 0              , (9*me2+9*me1)/140       , 0                      , 0              , 0                           , -(6*L*me2+7*L*me1)/420]      , 
        [0              , 0                      , (3*me2+10*me1)/35       , 0              , -(7*L*me2+15*L*me1)/420     , 0                           , 0              , 0                       , (9*me2+9*me1)/140      , 0              , (6*L*me2+7*L*me1)/420       , 0]                           , 
        [0              , 0                      , 0                       , (me2+3*me1)/12 , 0                           , 0                           , 0              , 0                       , 0                      , (me2+me1)/12   , 0                           , 0]                           , 
        [0              , 0                      , -(7*L*me2+15*L*me1)/420 , 0              , (3*L**2*me2+5*L**2*me1)/840 , 0                           , 0              , 0                       , -(7*L*me2+6*L*me1)/420 , 0              , -(L**2*me2+L**2*me1)/280    , 0]                           , 
        [0              , (7*L*me2+15*L*me1)/420 , 0                       , 0              , 0                           , (3*L**2*me2+5*L**2*me1)/840 , 0              , (7*L*me2+6*L*me1)/420   , 0                      , 0              , 0                           , -(L**2*me2+L**2*me1)/280]    , 
        [(me2+me1)/12   , 0                      , 0                       , 0              , 0                           , 0                           , (3*me2+me1)/12 , 0                       , 0                      , 0              , 0                           , 0]                           , 
        [0              , (9*me2+9*me1)/140      , 0                       , 0              , 0                           , (7*L*me2+6*L*me1)/420       , 0              , (10*me2+3*me1)/35       , 0                      , 0              , 0                           , -(15*L*me2+7*L*me1)/420]     , 
        [0              , 0                      , (9*me2+9*me1)/140       , 0              , -(7*L*me2+6*L*me1)/420      , 0                           , 0              , 0                       , (10*me2+3*me1)/35      , 0              , (15*L*me2+7*L*me1)/420      , 0]                           , 
        [0              , 0                      , 0                       , (me2+me1)/12   , 0                           , 0                           , 0              , 0                       , 0                      , (3*me2+me1)/12 , 0                           , 0]                           , 
        [0              , 0                      , (6*L*me2+7*L*me1)/420   , 0              , -(L**2*me2+L**2*me1)/280    , 0                           , 0              , 0                       , (15*L*me2+7*L*me1)/420 , 0              , (5*L**2*me2+3*L**2*me1)/840 , 0]                           , 
        [0              , -(6*L*me2+7*L*me1)/420 , 0                       , 0              , 0                           , -(L**2*me2+L**2*me1)/280    , 0              , -(15*L*me2+7*L*me1)/420 , 0                      , 0              , 0                           , (5*L**2*me2+3*L**2*me1)/840]
        ])
    return ke, me



def frame3d_KeMe(E,G,Kv,EA,EIx,EIy,EIz,L,A,Mass,R = None): 
    """ 
    Stiffness and mass matrices for Hermitian beam element with 6DOF per node
    Euler-Bernoulli beam model. The torsion is de-coupled to the rest, not like Timoshenko.
        
    The beam coordinate system is such that the cross section is assumed to be in the y-z plane
        
    (ux uy uz thetax thetay thetaz)
    (ux1 uy1 uz1 tx1 ty1 tz1 ux2 uy2 yz2 tx2 ty2 tz2)
        
    The torsional equation is fully decoupled and as follows:
      Ipx/3A Mass txddot  + G Kv /L tx = Torsional Moment
        
    INPUTS
        E : Young's (elastic) modulus
        Gs: Shear modulus. For an isotropic material G = E/2(nu+1) with nu the Poission's ratio
        Kv: Saint-Venant's torsion constant, Polar moment of i
        L :    Element length
        A :    Cross section area
        Mass :    Element mass = rho * A * L [kg]
        EA  : Young Modulus times Cross section.
        EIx : Young Modulus times Polar  second moment of area,local x-axis. Ix=\iint(y^2+z^2) dy dz [m4]
        EIy : Young Modulus times Planar second moment of area,local y-axis. Iy=\iint z^2 dy dz [m4]
        EIz : Young Modulus times Planar second moment of area,local z-axis. Iz=\iint y^2 dy dz [m4]
    OPTIONAL INPUTS
        R   : Transformation matrix (3x3) from global coord to element coord: x_e = R.x_g
             if provided, element matrix is provided in global coord
    OUTPUTS
        ke: Element stiffness matrix (12x12)
        me: Element mass matrix (12x12)
        
    AUTHOR: E. Branlard
    """
    
    # --- Stiffness matrix
    a = EA / L
    b = 12 * EIz / L ** 3
    c = 6 * EIz / L ** 2
    d = 12 * EIy / L ** 3
    e = 6 * EIy / L ** 2
    f = G * Kv / L
    g = 2 * EIy / L
    h = 2 * EIz / L
    # NOTE: OK with
    #          - Serano beam3e function
    #          - Matlab FEM Book
    #          - frame3d_6j
    #          - Panzer-Hubele
    # NOTE: compatible with Timoshenko with shear offsets
    #    ux1 uy1 uz1 tx1 ty1 tz1   ux2  uy2 yz2 tx2 ty2 tz2
    ke = np.array([
        [a  , 0  , 0  , 0  , 0     , 0     , -a , 0  , 0  , 0  , 0     , 0]     , 
        [0  , b  , 0  , 0  , 0     , c     , 0  , -b , 0  , 0  , 0     , c]     , 
        [0  , 0  , d  , 0  , -e    , 0     , 0  , 0  , -d , 0  , -e    , 0]     , 
        [0  , 0  , 0  , f  , 0     , 0     , 0  , 0  , 0  , -f , 0     , 0]     , 
        [0  , 0  , -e , 0  , 2 * g , 0     , 0  , 0  , e  , 0  , g     , 0]     , 
        [0  , c  , 0  , 0  , 0     , 2 * h , 0  , -c , 0  , 0  , 0     , h]     , 
        [-a , 0  , 0  , 0  , 0     , 0     , a  , 0  , 0  , 0  , 0     , 0]     , 
        [0  , -b , 0  , 0  , 0     , -c    , 0  , b  , 0  , 0  , 0     , -c]    , 
        [0  , 0  , -d , 0  , e     , 0     , 0  , 0  , d  , 0  , e     , 0]     , 
        [0  , 0  , 0  , -f , 0     , 0     , 0  , 0  , 0  , f  , 0     , 0]     , 
        [0  , 0  , -e , 0  , g     , 0     , 0  , 0  , e  , 0  , 2 * g , 0]     , 
        [0  , c  , 0  , 0  , 0     , h     , 0  , -c , 0  , 0  , 0     , 2 * h]
        ])
    # ---
    # SOURCE: frame3d_6j.m -  Structural analysis of 3d frame by Ace_ventura  - https://se.mathworks.com/matlabcentral/fileexchange/49559-structural-analysis-of-3d-frames
    # NOTE: compatible with above
    # GIx=G*Kv;
    # ke2=[EA/L      0             0             0      0            0            -EA/L 0             0             0      0            0
    #      0         12*EIz/(L^3)  0             0      0            6*EIz/(L^2)  0     -12*EIz/(L^3) 0             0      0            6*EIz/(L^2)
    #      0         0             12*EIy/(L^3)  0      -6*EIy/(L^2) 0            0     0             -12*EIy/(L^3) 0      -6*EIy/(L^2) 0
    #      0         0             0             GIx/L  0            0            0     0             0             -GIx/L 0            0
    #      0         0             -6*EIy/(L^2)  0      4*EIy/L      0            0     0             6*EIy/(L^2)   0      2*EIy/L      0
    #      0         6*EIz/(L^2)   0             0      0            4*EIz/L      0     -6*EIz/(L^2)  0             0      0            2*EIz/L
    #      -EA/L     0             0             0      0            0            EA/L  0             0             0      0            0
    #      0         -12*EIz/(L^3) 0             0      0            -6*EIz/(L^2) 0     12*EIz/(L^3)  0             0      0            -6*EIz/(L^2)
    #      0         0             -12*EIy/(L^3) 0      6*EIy/(L^2)  0            0     0             12*EIy/(L^3)  0      6*EIy/(L^2)  0
    #      0         0             0             -GIx/L 0            0            0     0             0             GIx/L  0            0
    #      0         0             -6*EIy/(L^2)  0      2*EIy/L      0            0     0             6*EIy/(L^2)   0      4*EIy/L      0
    #      0         6*EIz/(L^2)   0             0      0            2*EIz/L      0     -6*EIz/(L^2)  0             0      0            4*EIz/L      ]; #formation of element stiffness matrix IN MEMBER AXIS

        # ---
    # SOURCE: Panzer-Hubele - Generating a Parametric Finite Element Model of a 3D Cantilever Timoshenko Beam Using Matlab
    # NOTE: compatible with above
    # Py=0; Pz=0; It=Kv; l=L; Iz=EIz/E; Iy=EIy/E;
    # K11 = zeros(6,6);
    # K11(1,1) = E * A/l ;
    # K11(2,2) = 12 * E * Iz/(l^3 * (1+Py)) ;
    # K11(2,6) = 6 * E * Iz/(l^2 * (1+Py)) ;
    # K11(3,3) = 12 * E * Iy/(l^3 * (1+Pz)) ;
    # K11(3,5) = -6 * E * Iy/(l^2 * (1+Pz)) ;
    # K11(4,4) = G * It/l ;
    # K11(5,5) = (4+Pz) * E * Iy/(l * (1+Pz)) ;
    # K11(5,3) = K11(3,5) ;
    # K11(6,6) = (4+Py) * E * Iz/(l * (1+Py)) ;
    # K11(6,2) = K11(2,6) ;

        # K22 = -K11 + 2 * diag(diag(K11));
    # K21 = K11 - 2 * diag(diag(K11));
    # K21(5,5) = (2-Pz) * E * Iy/(l * (1+Pz)) ;
    # K21(6,6) = (2-Py) * E * Iz/(l * (1+Py)) ;
    # K21(2,6) = -K21(6,2);
    # K21(3,5) = -K21(5,3);
    # ke3 = [K11, K21'; K21, K22];
    # ---  Mass matrix
    # SOURCE: What-When-How-FEM-For-Frames. NOTE: the sign was reveresed in front of 35*r2!!!, to be consistent with Panzer-Hubele with Iy and Iz=0
    a = L / 2
    a2 = a ** 2
    r2 = EIx / E / A
    me = Mass / 2 / 105 * np.array( [
                [70 , 0       , 0       , 0       , 0       , 0       , 35 , 0       , 0       , 0       , 0       , 0]       , 
                [0  , 78      , 0       , 0       , 0       , 22 * a  , 0  , 27      , 0       , 0       , 0       , -13 * a] , 
                [0  , 0       , 78      , 0       , -22 * a , 0       , 0  , 0       , 27      , 0       , 13 * a  , 0]       , 
                [0  , 0       , 0       , 70 * r2 , 0       , 0       , 0  , 0       , 0       , 35 * r2 , 0       , 0]       , 
                [0  , 0       , -22 * a , 0       , 8 * a2  , 0       , 0  , 0       , -13 * a , 0       , -6 * a2 , 0]       , 
                [0  , 22 * a  , 0       , 0       , 0       , 8 * a2  , 0  , 13 * a  , 0       , 0       , 0       , -6 * a2] , 
                [35 , 0       , 0       , 0       , 0       , 0       , 70 , 0       , 0       , 0       , 0       , 0]       , 
                [0  , 27      , 0       , 0       , 0       , 13 * a  , 0  , 78      , 0       , 0       , 0       , -22 * a] , 
                [0  , 0       , 27      , 0       , -13 * a , 0       , 0  , 0       , 78      , 0       , 22 * a  , 0]       , 
                [0  , 0       , 0       , 35 * r2 , 0       , 0       , 0  , 0       , 0       , 70 * r2 , 0       , 0]       , 
                [0  , 0       , 13 * a  , 0       , -6 * a2 , 0       , 0  , 0       , 22 * a  , 0       , 8 * a2  , 0]       , 
                [0  , -13 * a , 0       , 0       , 0       , -6 * a2 , 0  , -22 * a , 0       , 0       , 0       , 8 * a2]
                ])

    ## Element in global coord
    if (R is not None):
        RR = np.blkdiag(R,R,R,R)
        me = np.transpose(RR).dot(me.dot(RR))
        ke = np.transpose(RR).dot(ke.dot(RR))
    return ke, me
    
    # ---
    # SOURCE: Panzer-Hubele - Generating a Parametric Finite Element Model of a 3D Cantilever Timoshenko Beam Using Matlab
    # Iz=EIz/E; Iy=EIy/E; Ip=EIx/E; l=L;
    # M11 = zeros(6,6);
    # M11(1,1) = 1/3;
    # M11(2,2) = 13/35 + 6 * Iz/(5 * A * l^2);
    # M11(3,3) = 13/35 + 6 * Iy/(5 * A * l^2);
    # M11(4,4) = Ip/(3 * A);
    # M11(5,5) = l^2/105 + 2 * Iy/(15 * A);
    # M11(6,6) = l^2/105 + 2 * Iz/(15 * A);
    # M11(6,2) = 11 * l/210 + Iz/(10 * A * l);
    # M11(2,6) = M11(6,2) ;
    # M11(5,3) = -11 * l/210 - Iy/(10 * A * l);
    # M11(3,5) = M11(5,3) ;
    # M22 = -M11 + 2 * diag(diag(M11));
    # M21 = zeros(6,6);
    # M21(1,1) = 1/6;
    # M21(2,2) = 9/70 - 6 * Iz/(5 * A * l^2);
    # M21(3,3) = 9/70 - 6 * Iy/(5 * A * l^2);
    # M21(4,4) = Ip/(6 * A);
    # M21(5,5) = -l^2/140 - Iy/(30 * A);
    # M21(6,6) = -l^2/140 - Iz/(30 * A);
    # M21(6,2) = -13 * l/420 + Iz/(10 * A * l);
    # M21(2,6) = -M21(6,2);
    # M21(5,3) = 13 * l/420 - Iy/(10 * A * l);
    # M21(3,5) = -M21(5,3);
    # me= Mass * [M11, M21'; M21, M22];
    # keyboard
        #   b=[ ex(2)-ex(1); ey(2)-ey(1); ez(2)-ez(1) ];
    #   L=sqrt(b'*b);  n1=b/L;
    #   lc=sqrt(eo*eo'); n3=eo/lc;
    #  # Elemenr load vector
    #      if nargin==5;   eq=[0 0 0 0];  end
    #  eq = [qx qy qz qw];    distributed loads
    #      qx=eq(1); qy=eq(2); qz=eq(3); qw=eq(4);
    #    fle=L/2*[qx qy qz qw -1/6*qz*L 1/6*qy*L qx qy qz qw 1/6*qz*L -1/6*qy*L]';
        #  #
    #     n2(1)=n3(2)*n1(3)-n3(3)*n1(2);
    #     n2(2)=-n1(3)*n3(1)+n1(1)*n3(3);
    #     n2(3)=n3(1)*n1(2)-n1(1)*n3(2);
    # #
    #     An=[n1';
    #         n2;
    #         n3];
    # #
    #     Grot=[  An     zeros(3) zeros(3) zeros(3);
    #        zeros(3)   An     zeros(3) zeros(3);
    #        zeros(3) zeros(3)   An     zeros(3);
    #        zeros(3) zeros(3) zeros(3)   An    ];
    #  #
    #     Ke1=Grot'*Kle*Grot;  fe1=Grot'*fle;


# !------------------------------------------------------------------------------------------------------
# !> calculates the lumped forces and moments due to gravity on a given element:
# !! the element has two nodes, with the loads for both elements stored in array F. Indexing of F is:
# !!    Fx_n1=1,Fy_n1=2,Fz_n1=3,Mx_n1= 4,My_n1= 5,Mz_n1= 6,
# !!    Fx_n2=7,Fy_n2=8,Fz_n2=9,Mx_n2=10,My_n2=11,Mz_n2=12
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
# !------------------------------------------------------------------------------------------------------


if __name__=='__main__':


    pass
