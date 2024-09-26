""" 
Set of utils useful for Structural and Multi-Body dynamics

"""

import numpy as np

# --- Definitions to ease comparison with sympy versions
from numpy import cos ,sin
from welib.yams.rotations import R_x, R_y, R_z

# def Matrix(m):
#     return np.asarray(m)

def skew(x, symb=False):
    """ Returns the skew symmetric matrix M, such that: cross(x,v) = M v 
    [ 0, -z , y]
    [ z,  0 ,-x]
    [-y,  x , 0]
    """
    if not symb:
        x=np.asarray(x).ravel()
        return np.array([[0, -x[2], x[1]],[x[2],0,-x[0]],[-x[1],x[0],0]])
    else:
        from sympy import Matrix
        if hasattr(x,'shape') and len(x.shape)==2:
            if x.shape[0]==3:
                return Matrix(np.array([[0, -x[2,0], x[1,0]],[x[2,0],0,-x[0,0]],[-x[1,0],x[0,0],0]]))
            else:
                raise Exception('fSkew expect a vector of size 3 or matrix of size 3x1, got {}'.format(x.shape))
        else:
            return Matrix(np.array([[0, -x[2], x[1]],[x[2],0,-x[0]],[-x[1],x[0],0]]))

def skew2(x):
    """ Returns the skew(x).skew(x) 
         [ 0 -z  y]   [ 0 -z  y]   [  -y**2 - z**2  ,       xy       ,         xz      ]
    S2 = [ z  0 -x] . [ z  0 -x] = [       yx       , - x**2 - z**2  ,         yz      ]
         [-y  x  0]   [-y  x  0]   [       zx       ,       zy       ,   - x**2 - y**2 ]
    """
    x,y,z=np.asarray(x).flatten()
    return np.array( [
            [ - y**2 - z**2  ,       x*y      ,         x*z      ],
            [       y*x      , - x**2 - z**2  ,         y*z      ],
            [       z*x      ,       z*y      ,   - x**2 - y**2 ]])

def extractVectFromSkew(M):
    """ Return a 3-vector from a skew matrix """
    # [ 0, -z , y]
    # [ z,  0 ,-x]
    # [-y,  x , 0]
    M = np.asarray(M)
    x = 0.5*( M[2,1]-M[1,2])
    y = 0.5*( M[0,2]-M[2,0])
    z = 0.5*(-M[0,1]+M[1,0])
    return np.array([x,y,z])


def buildRigidBodyMassMatrix(Mass, J_P, COG=None, symb=False): 
    """ 
    Simply builds the mass matrix of a rigid body (i.e. mass matrix) Eq.(15) of [1] 
    INPUTS:
      - Mass: (scalar) mass of the body
      - JP: (3x3 matrix) full inertia matrix at Point P
      - COG: (3-vector) x,y,z position of center of mass
    """
    if symb:
        from sympy import zeros, Matrix
        MM = zeros(6,6)
    else:
        MM = np.zeros((6,6));
    S = Mass*skew(COG, symb=symb)
    MM[0,0] = Mass
    MM[1,1] = Mass
    MM[2,2] = Mass
    MM[0:3,3:6] = -S;
    MM[3:6,0:3] = S ; # transpose(S)=-S;
    MM[3:6,3:6] = J_P ;
    return MM

def rigidBodyMassMatrixAtP(m=None, J_G=None, Ref2COG=None, symb=False):
    """ 
    Rigid body mass matrix (6x6) at a given reference point: 
      the center of gravity (if Ref2COG is None) 


    INPUTS:
     - m/tip: (scalar) body mass 
                     default: None, no mass
     - J_G: (3-vector or 3x3 matrix), diagonal coefficients or full inertia matrix
                     with respect to COG of body! 
                     The inertia is transferred to the reference point if Ref2COG is not None
                     default: None 
     - Ref2COG: (3-vector) x,y,z position of center of gravity (COG) with respect to a reference point
                     default: None, at first/last node.
    OUTPUTS:
      - M66 (6x6) : rigid body mass matrix at COG or given point 
    """
    # TODO SYMPY
    # Default values
    if m is None: m=0
    if Ref2COG is None: Ref2COG=(0,0,0)
    if J_G is None: J_G=np.zeros((3,3))
    if len(J_G.flatten()==3): J_G = np.eye(3).dot(J_G)

    M66 = np.zeros((6,6))
    x,y,z = Ref2COG
    Jxx,Jxy,Jxz = J_G[0,:]
    _  ,Jyy,Jyz = J_G[1,:]
    _  ,_  ,Jzz = J_G[2,:]
    M66[0, :] =[   m     ,   0     ,   0     ,   0                 ,  z*m                , -y*m                 ]
    M66[1, :] =[   0     ,   m     ,   0     , -z*m                ,   0                 ,  x*m                 ]
    M66[2, :] =[   0     ,   0     ,   m     ,  y*m                , -x*m                ,   0                  ]
    M66[3, :] =[   0     , -z*m    ,  y*m    , Jxx + m*(y**2+z**2) , Jxy - m*x*y         , Jxz  - m*x*z         ]
    M66[4, :] =[  z*m    ,   0     , -x*m    , Jxy - m*x*y         , Jyy + m*(x**2+z**2) , Jyz  - m*y*z         ]
    M66[5, :] =[ -y*m    , x*m     ,   0     , Jxz - m*x*z         , Jyz - m*y*z         , Jzz  + m*(x**2+y**2) ]
    return M66

def identifyRigidBodyMM(MM):
    """ 
    Based on a 6x6 mass matrix at a reference point:
     - Identify the position of the center of mass
     - Compute the inertia at the center of mass
    """
    mass = MM[0,0]
    # Using average of two coeffs to get estimate of COG
    xCM = 0.5*( MM[1,5]-MM[2,4])/mass
    zCM = 0.5*( MM[0,4]-MM[1,3])/mass
    yCM = 0.5*(-MM[0,5]+MM[2,3])/mass
    # Distance from refopint to COG
    Ref2COG=np.array((xCM,yCM,zCM))
    # Inertia at ref oint
    J_P = MM[3:6,3:6].copy()
    # Inertia at COG
    J_G = translateInertiaMatrixToCOG(J_P, mass, r_PG=Ref2COG ) 
    return mass, J_G, Ref2COG


def translateRigidBodyMassMatrix(M, r_P1P2):
    """ 
    Translate a 6x6 rigid mass matrix from point 1 to point 2
    r_P1P2: vector from point1 to point2 
    """
    # First identify main properties (mass, inertia, location of center of mass from previous ref point)
    mass, J_G, Ref2COG = identifyRigidBodyMM(M)
    # New COG location from new (x,y) ref point
    Ref2COG -= np.asarray(r_P1P2)
    # Compute mass matrix 
    M_new =  rigidBodyMassMatrixAtP(mass, J_G, Ref2COG)
    return M_new


def rotateRigidBodyMassMatrix(M_ss, R_s2d):
    """ 
    Rotate a 6x6 rigid body mass matrix from a source (s) coordinate to destination (d) coordinate system

    INPUTS:
     - M_ss: mass matrix in source coordinate, 6x6 array
     - R_s2d: transformation matrix source two destination
    """
    R66_s2d = np.block(R_s2d)
    M_dd = R66_s2d.dot(M_ss).dot(R66_s2d.T)
    return M_dd



# --------------------------------------------------------------------------------}
# --- Inertia functions 
# --------------------------------------------------------------------------------{
def translateInertiaMatrix(I_A, Mass, r_BG, r_AG = np.array([0,0,0])):
    """
    Transform inertia matrix with respect to point A to the inertia matrix with respect to point B
    NOTE: the vectors and the inertia matrix needs to be expressed in the same coordinate system. 
    NOTE: one of the vector r_BG or r_AG may be empty or 0 instead of [0,0,0];
    NOTE: if r_AG is not provided it is assumed to be 0, i.e. A=G
    To avoid this confusion you can use translateInertiaMatrixFromCOG  and translateInertiaMatrixToCOG
    
    INPUTS:
       I_A  : Inertia matrix 3x3 in the coordinate system A
       Mass : Mass of the body
       r_BG: vector from point B to COG of the body
    
    OPTIONAL INPUTS:
       r_AG: vector from point A to point G
    """
    if len(r_AG) < 3:
        r_AG = np.array([0,0,0])
    if len(r_BG) < 3:
        r_BG = np.array([0,0,0])   
    I_B = I_A - Mass*(skew2(r_BG)-skew2(r_AG))
    #I_G = translateInertiaMatrixToCOG(I_A, Mass, r_AG)
    #I_B = translateInertiaMatrixFromCOG(I_G, Mass, -np.array(r_BG))
    return I_B

def translateInertiaMatrixToCOG(I_P, Mass, r_PG): 
    """ Transform inertia matrix with respect to point P to the inertia matrix with respect to the COG
    NOTE: the vectors and the inertia matrix needs to be expressed in the same coordinate system.
    
    INPUTS:
      I_P  : Inertia matrix 3x3 with respect to point P
      Mass : Mass of the body
      r_PG: vector from P to COG 
    """
    I_G = I_P + Mass * skew2(r_PG)
    return I_G

def translateInertiaMatrixFromCOG(I_G, Mass, r_GP): 
    """
    Transform inertia matrix with respect to COG to the inertia matrix with respect to point P
    NOTE: the vectors and the inertia matrix needs to be expressed in the same coordinate system.
    INPUTS:
       I_G  : Inertia matrix 3x3 with respect to COG
       Mass : Mass of the body
       r_GP: vector from COG of the body to point P
    """
    I_P = I_G - Mass * skew2(r_GP)
    return I_P
    


# --------------------------------------------------------------------------------}
# --- Loads 
# --------------------------------------------------------------------------------{
def transferLoadsZPoint(ls, z, phi_x, phi_y, phi_z, rot_type='default'):
    """ 
    Used to transfer loads from HydroF*i to the Platform ref point

    HydroF*i are translated similarly to the ref point, but due to rotation, an extra laver arm is present

    Undisplaced  ->   Displaced:
           P     ->      P
           |              \
           |               \
           0     ->         0

    HydroF*i are computed at "0" above. 
    They need to be transfered to P

    z: destination to source (z_s - z_d) = (z_0 - z_P)
    """
    # --- Rotation
    s_b=[0,0,z]
    from welib.yams.rotations import BodyXYZ_A, smallRot_OF, smallRot_A
    if rot_type=='default' or rot_type=='bodyXYZ':
        #  BodyXYZ_A.dot(s_b)
        r = (  z*np.sin(phi_y) , -z * np.sin(phi_x) * np.cos(phi_y),  z *np.cos(phi_x)* np.cos(phi_y))

    elif rot_type=='smallRot_OF':
        #  smallRot_OF.T .dot(s_b)
        r = np.zeros((3,len(phi_z)))
        for i,(phi_x1,phi_y1, phi_z1) in enumerate(zip(phi_x,phi_y,phi_z)):
            R_b2g = smallRot_OF(phi_x1, phi_y1, phi_z1).T
            r[:,i] = R_b2g.dot(s_b)

    elif rot_type=='smallRot':
        #  smallRot_A .dot(s_b)
        r = (z*phi_y, -z*phi_x , z)

    else:
        raise Exception()
    ld    = np.zeros(ls.shape)
    # Forces
    ld[0] = ls[0]
    ld[1] = ls[1]
    ld[2] = ls[2]
    # Moments
    ld[3] = ls[3] + r[1] * ls[2] - r[2] * ls[1]
    ld[4] = ls[4] + r[2] * ls[0] - r[0] * ls[2]
    ld[5] = ls[5] + r[0] * ls[1] - r[1] * ls[0]
    return ld
