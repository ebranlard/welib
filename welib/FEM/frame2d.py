import numpy as np

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
# ---
def uax(x, u_1, u_4, L):
    """ longitudinal/axial displacement field"""
    s=x/L
    return  u_1*b1(s) + u_4*b4(s)

def h(x, u_2, u_3, u_5, u_6, L):
    """ Transverse displacement field """
    s=x/L
    return  u_2*b2(s) + u_3*b3(s,L) +u_5*b5(s) + u_6*b6(s,L)


# --------------------------------------------------------------------------------}
# --- Element formulation 
# --------------------------------------------------------------------------------{
def frame2d_KeMe(EA,EI,L,Mass,T=0,theta=None, MMFormulation='consistent'):
    r""" 
    Stiffness and mass matrices for Hermitian beam element with 3DOF per node.
      - Elongation (along x)
      - Bending in one transverse direction (y)
    Euler-Bernoulli beam model. 
        
    The beam coordinate system is such that the cross section is assumed to be in the y-z plane
    (along x)
        
    Nodal DOF   : (ux uy theta)
    Element DOFs: (ux1 uy1 t1 ux2 uy2 t2) or (u1,u2,u3,u4,u5,u6)
        
    INPUTS
        EA  : Young Modulus times Cross section.
        EI  : Young Modulus times Planar second moment of area,local y-axis. Iy=\iint z^2 dy dz [m4]
        L   :    Element length
        Mass :    Element mass = rho * A * L [kg]

    OPTIONAL INPUTS
        theta: orientation angle of element in global coordinate system
               if provided, element matrix is provided in global coord
        T   : Axial load to be used for the compuation of the geometrical stiffness

    OUTPUTS
        Ke: Element stiffness matrix (6x6)
        Me: Element mass matrix      (6x6)
        Kg: Element geometrical stiffness matrix (6,6)
        
    AUTHOR: E. Branlard
    """
    # NOTE: matrices determined using sympy, see scripts in current folder 

    # --- Stiffness matrices
    a = EA / L
    c = EI / (L ** 3)
    Ke = np.array([
              [a  , 0      , 0        , -a , 0       , 0]        , 
              [0  , 12*c   , 6*L*c    , 0  , - 12*c  , 6*L*c]    , 
              [0  , 6*L*c  , 4*L**2*c , 0  , - 6*L*c , 2*L**2*c] , 
              [-a , 0      , 0        , a  , 0       , 0]        , 
              [0  , - 12*c , - 6*L*c  , 0  , 12*c    , - 6*L*c]  , 
              [0  , 6*L*c  , 2*L**2*c , 0  , - 6*L*c , 4*L**2*c]
        ])
    #Ke = np.array( [
    #        [EA/L  , 0           , 0          , -EA/L , 0           , 0]          , 
    #        [0     , 12*EI/L**3  , 6*EI/L**2  , 0     , -12*EI/L**3 , 6*EI/L**2]  , 
    #        [0     , 6*EI/L**2   , 4*EI/L     , 0     , -6*EI/L**2  , 2*EI/L]     , 
    #        [-EA/L , 0           , 0          , EA/L  , 0           , 0]          , 
    #        [0     , -12*EI/L**3 , -6*EI/L**2 , 0     , 12*EI/L**3  , -6*EI/L**2] , 
    #        [0     , 6*EI/L**2   , 2*EI/L     , 0     , -6*EI/L**2  , 4*EI/L]
    #        ])
    Kg = np.array([
            [0 , 0          , 0        , 0 , 0          , 0]        , 
            [0 , 6*T/(5*L)  , T/10     , 0 , -6*T/(5*L) , T/10]     , 
            [0 , T/10       , 2*L*T/15 , 0 , -T/10      , -L*T/30]  , 
            [0 , 0          , 0        , 0 , 0          , 0]        , 
            [0 , -6*T/(5*L) , -T/10    , 0 , 6*T/(5*L)  , -T/10]    , 
            [0 , T/10       , -L*T/30  , 0 , -T/10      , 2*L*T/15]
        ])

    # --- Mass Matrix
    if MMFormulation=='consistent':
        # consistent mass matrix
        mm = Mass / 420
        ma = Mass / 6
        Me = np.array([
            [2*ma , 0         , 0           , ma   , 0         , 0]           , 
            [0    , 156*mm    , 22*L*mm     , 0    , 54*mm     , - 13*L*mm]   , 
            [0    , 22*L*mm   , 4*L**2*mm   , 0    , 13*L*mm   , - 3*L**2*mm] , 
            [ma   , 0         , 0           , 2*ma , 0         , 0]           , 
            [0    , 54*mm     , 13*L*mm     , 0    , 156*mm    , - 22*L*mm]   , 
            [0    , - 13*L*mm , - 3*L**2*mm , 0    , - 22*L*mm , 4*L**2*mm]
           ])
         #Me =np.arrya([
         #    [m/3 , 0           , 0           , m/6 , 0           , 0]           , 
         #    [0   , 13*m/35     , 11*L*m/210  , 0   , 9*m/70      , -13*L*m/420] , 
         #    [0   , 11*L*m/210  , L**2*m/105  , 0   , 13*L*m/420  , -L**2*m/140] , 
         #    [m/6 , 0           , 0           , m/3 , 0           , 0]           , 
         #    [0   , 9*m/70      , 13*L*m/420  , 0   , 13*m/35     , -11*L*m/210] , 
         #    [0   , -13*L*m/420 , -L**2*m/140 , 0   , -11*L*m/210 , L**2*m/105]
         #    ])
    elif MMFormulation =='lumped':
        # lumped mass matrix
        Me = Mass * np.diag([1/2, 1/2, 0, 1/2, 1/2, 0])
    elif MMFormulation =='diagonal':
        # diagonal mass matrix
        Me = Mass * np.diag([1/2, 1/2, L**2/78, 1/2, 1/2, L**2/78])
    
    # --- Conversion to global system if requested
    if theta is not None:
        R = np.array([
            [np.cos(theta)   , np.sin(theta) , 0 , 0               , 0             , 0]   , 
            [- np.sin(theta) , np.cos(theta) , 0 , 0               , 0             , 0]   , 
            [0               , 0             , 1 , 0               , 0             , 0]   , 
            [0               , 0             , 0 , np.cos(theta)   , np.sin(theta) , 0]   , 
            [0               , 0             , 0 , - np.sin(theta) , np.cos(theta) , 0]   , 
            [0               , 0             , 0 , 0               , 0             , 1]])

        Me = np.transpose(R) * Me * R
        Ke = np.transpose(R) * Ke * R
        Kg = np.transpose(R) * Kg * R

    return Ke, Me, Kg
