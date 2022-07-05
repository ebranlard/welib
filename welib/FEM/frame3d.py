import numpy as np
import sympy
import scipy
    
# --------------------------------------------------------------------------------}
# --- Shape functions, displacement field, energy
# --------------------------------------------------------------------------------{
def b1(s) :
    return 1-s 
def b4(s) :
    return s 
def b2(s) :
    return 1 -3*s**2 + 2*s**3
def b3(s,L) :
    return L*s*(1-s)**2       
def b5(s) :
    return 3*s**2 - 2*s**3     
def b6(s,L) :
    return  L*s**2*(s-1)       

def frame3d_N(x,L):
    """  Interpolation matrix from nodal DOF to deflections for a 3d frame
     [   u_x   ]  = N(x) * [ux1,uy1,uz1,tx1,...,tz2]^t
     [   u_y   ]
     [   u_z   ]
     [\theta_x ]
    """
    if isinstance(x, sympy.Symbol):
        NN = sympy.zeros(4,12)
    else:
        NN = np.zeros((4,12))
    NN[0,0]  = b1(x/L)
    NN[1,1]  = b2(x/L)
    NN[2,2]  = b2(x/L)
    NN[3,3]  = b1(x/L)
    NN[2,4]  = -b3(x/L,L)
    NN[1,5]  = b3(x/L,L)
    NN[0,6]  = b4(x/L)
    NN[1,7]  = b5(x/L)
    NN[2,8]  = b5(x/L)
    NN[3,9]  = b4(x/L)
    NN[2,10] = -b6(x/L,L)
    NN[1,11] = b6(x/L,L)
    return NN


# ---
def frame3d_u(x,q,L):
    """ Returns deflections for a 3d frame:
    [   u_x   ]  = N(x) * q   with q = [ux1,uy1,uz1,tx1,...,tz2]^t
    [   u_y   ]
    [   u_z   ]
    [\theta_x ]
    """
    s=x/L
    N = frame3d_N(x,L)
    if isinstance(N, sympy.Matrix):
        return N * q
    else:
        return N.dot(q)

def h(x, u_2, u_3, u_5, u_6, L):
    """ Transverse displacement field """
    s=x/L
    return  u_2*b2(s) + u_3*b3(s,L) +u_5*b5(s) + u_6*b6(s,L)


# --------------------------------------------------------------------------------}
# --- Element formulation 
# --------------------------------------------------------------------------------{

def frame3d_KeMe(E,G,Kv,EA,EIx,EIy,EIz,L,A,Mass,T=0,R=None, main_axis='x'): 
    """ 
    Stiffness and mass matrices for Hermitian beam element with 6DOF per node
    Beam directed along x
    See folder derivation for sympy derivations.

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
        T   : Axial loads (along x) to use for the geometric stiffness matrix Kg
        R   : Transformation matrix (3x3) from global coord to element coord: x_e = R.x_g
             if provided, element matrix is provided in global coord
    OUTPUTS
        ke: Element stiffness matrix (12x12)
        me: Element mass matrix (12x12)
        Kg: Element geometrical stiffness matrix (12x12)

        
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
    Ke = np.array([
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
    # ---  Mass matrix
    a = L / 2
    a2 = a ** 2
    r2 = EIx / E / A
    Me = Mass / 2 / 105 * np.array( [
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

    # --- Geometrical stiffness matrix
    Kg= T* np.array([
        [0 , 0         , 0         , 0 , 0      , 0      , 0 , 0         , 0         , 0 , 0      , 0]      , 
        [0 , 6./(5*L)  , 0         , 0 , 0      , 1./10  , 0 , -6./(5*L) , 0         , 0 , 0      , 1./10]  , 
        [0 , 0         , 6./(5*L)  , 0 , -1./10 , 0      , 0 , 0         , -6./(5*L) , 0 , -1./10 , 0]      , 
        [0 , 0         , 0         , 0 , 0      , 0      , 0 , 0         , 0         , 0 , 0      , 0]      , 
        [0 , 0         , -1./10    , 0 , 2*L/15 , 0      , 0 , 0         , 1./10     , 0 , -L/30  , 0]      , 
        [0 , 1./10     , 0         , 0 , 0      , 2*L/15 , 0 , -1./10    , 0         , 0 , 0      , -L/30] , 
        [0 , 0         , 0         , 0 , 0      , 0      , 0 , 0         , 0         , 0 , 0      , 0]      , 
        [0 , -6./(5*L) , 0         , 0 , 0      , -1./10 , 0 , 6./(5*L)  , 0         , 0 , 0      , -1./10] , 
        [0 , 0         , -6./(5*L) , 0 , 1./10  , 0      , 0 , 0         , 6./(5*L)  , 0 , 1./10  , 0]      , 
        [0 , 0         , 0         , 0 , 0      , 0      , 0 , 0         , 0         , 0 , 0      , 0]      , 
        [0 , 0         , -1./10    , 0 , -L/30  , 0      , 0 , 0         , 1./10     , 0 , 2*L/15 , 0]      , 
        [0 , 1./10     , 0         , 0 , 0      , -L/30  , 0 , -1./10    , 0         , 0 , 0      , 2*L/15]
       ])



    ## Element in global coord
    if (R is not None):
        RR = scipy.linalg.block_diag(R,R,R,R)
        Me = np.transpose(RR).dot(Me).dot(RR)
        Ke = np.transpose(RR).dot(Ke).dot(RR)
        Kg = np.transpose(RR).dot(Kg).dot(RR)
    return Ke, Me, Kg
    
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
    #     Grot=[  An     zeros(3) zeros(3) zeros(3);
    #        zeros(3)   An     zeros(3) zeros(3);
    #        zeros(3) zeros(3)   An     zeros(3);
    #        zeros(3) zeros(3) zeros(3)   An    ];
    #     Ke1=Grot'*Kle*Grot;  fe1=Grot'*fle;



if __name__=='__main__':


    pass
