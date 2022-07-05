"""
Rotation matrices for different rotational coordinates conventions: 
 - EulerP : Euler parameters
 - BodyZXZ : Euler 
 - BodyXYZ : Bryant 
 - BodyZYX :  
"""
import numpy as np
from numpy import cos, sin, tan, arccos, trace



# --------------------------------------------------------------------------------}
# --- Euler parameters 
# --------------------------------------------------------------------------------{
def EulerP_A(e0,e1,e2,e3):
    """ Transformation matrix body to global for Euler parameters
    Shabana (2.33), Nikravesh (6.19)
    """
    A = 2*np.array([
       [e0**2+e1**2-0.5 , e1*e2-e0*e3     , e1*e3+e0*e2    ], 
       [e1*e2+e0*e3     , e0**2+e2**2-0.5 , e2*e3-e0*e1    ], 
       [e1*e3-e0*e2     , e2*e3+e0*e1     , e0**2+e3**2-0.5]])
    return A

def EulerP_G(e0,e1,e2,e3):
    """ Angular velocity matrix such that omega_global = G theta_dot for Euler parameters
        G = 2 E
    Shabana  (2.35, 2.54)
    """
    G = 2*np.array([
        [-e1 , e0, -e3,  e2], 
        [-e2 , e3,  e0, -e1],
        [-e3 ,-e2,  e1,  e0]])
    return G 

def EulerP_Gb(e0,e1,e2,e3):
    """ Angular velocity matrix such that omega_local = Gb theta_dot for Euler parameters
        Gb = 2 Eb
    Shabana  (2.35, 2.54)
    """
    Gb = 2* np.array([
        [-e1,  e0,  e3, -e2],
        [-e2, -e3,  e0,  e1],
        [-e3,  e2, -e1,  e0 ]]);
    return Gb

def EulerP_Gbinv(e0,e1,e2,e3):
    """ Inverse of "G" matrix, such that  p_dot = Gbinv omega_local """
    return  0.5* EulerP_Eb(e0,e1,e2,e3).T 

def EulerP_E(e0,e1,e2,e3):
    """ 
    Shabana (2.35), Nikravesh "G" (6.38)
    """
    E = np.array([
        [-e1 , e0, -e3,  e2], 
        [-e2 , e3,  e0, -e1],
        [-e3 ,-e2,  e1,  e0]])
    return E

def EulerP_Eb(e0,e1,e2,e3):
    """ 
    Shabana (2.35), Nikravesh "L" (6.39)
    """
    Eb = np.array([
        [-e1,  e0,  e3, -e2],
        [-e2, -e3,  e0,  e1],
        [-e3,  e2, -e1,  e0 ]]);
    return Eb
 
def EulerP_toBodyZXZ(e0,e1,e2,e3):
    """ 
    From Euler parameters to Euler angles
    Shabana 2.107
    """
    theta = arccos(2*(e0**2 + e3**2)-1) 
    phi   = arccos(-2*(e2*e3- e0*e1)/sin(theta)) 
    psi   = arccos( 2*(e2*e3+ e0*e1)/sin(theta)) 
    return phi, theta, psi

def EulerP_toBodyXYZ(e0,e1,e2,e3):
    A = EulerP_A(e0,e1,e2,e3)
    raise NotImplementedError()
    return phi_x, phi_y, phi_z

def EulerP_toBodyZYX(e0,e1,e2,e3):
    A = EulerP_A(e0,e1,e2,e3)
    raise NotImplementedError()
    return phi_x, phi_y, phi_z


def EulerP_fromA(A):
    """ 
    Find the 4 Euler parameters from a transformation matrix
    """
    trA = trace(A)
    # Squared coeffs,  1 = e_02 + e_12 + e_22 +e_32
    e_02 = (trA +1)/4
    e_12 = (1+2*A[0,0]-trA)/4
    e_22 = (1+2*A[1,1]-trA)/4
    e_32 = (1+2*A[2,2]-trA)/4
    # Selecting max coeffs to avoid division by 0
    imax = np.argmax([e_02, e_12,e_22, e_32])
    if imax==0:
        e0 = np.sqrt(e_02)
        e1 = (A[2,1]-A[1,2])/(4*e0)
        e2 = (A[0,2]-A[2,0])/(4*e0)
        e3 = (A[1,0]-A[0,1])/(4*e0)
    elif imax==1:
        e1 = np.sqrt(np.abs(e_12))
        e0 = (A[2,1]-A[1,2])/(4*e1)
        e2 = (A[1,0]+A[0,1])/(4*e1)
        e3 = (A[2,0]+A[0,2])/(4*e1)
    elif imax==2:
        e2 = np.sqrt(e_22)
        e0 = (A[0,2]-A[2,0])/(4*e2)
        e1 = (A[0,1]+A[1,0])/(4*e2)
        e3 = (A[1,2]+A[1,2])/(4*e2)
    else:
        e3 = np.sqrt(e_32)
        e0 = (A[0,1]-A[1,0])/(4*e3)
        e1 = (A[0,2]+A[2,0])/(4*e3)
        e2 = (A[1,2]+A[2,1])/(4*e3)

    return e0,e1,e2,e3


# --------------------------------------------------------------------------------}
# --- Euler angles ZXZ 
# --------------------------------------------------------------------------------{
def BodyZXZ_A(phi, theta, psi):
    """ Transformation matrix body to global for rotation of type Body and order ZXZ"""
    A = np.zeros((3,3))
    A[0,:] = [-sin(phi)*sin(psi)*cos(theta)+cos(phi)*cos(psi),-sin(phi)*cos(psi)*cos(theta)-sin(psi)*cos(phi),sin(phi)*sin(theta)]
    A[1,:] = [sin(phi)*cos(psi)+sin(psi)*cos(phi)*cos(theta),-sin(phi)*sin(psi)+cos(phi)*cos(psi)*cos(theta),-sin(theta)*cos(phi)]
    A[2,:] = [sin(psi)*sin(theta),sin(theta)*cos(psi),cos(theta)]
    return A

def BodyZXZ_G(phi, theta, psi):
    """ Angular velocity matrix such that omega_global = G theta_dot for rotation of type Body and order ZXZ"""
    G = np.zeros((3,3))
    G[0,:] = [0,cos(phi),sin(phi)*sin(theta)]
    G[1,:] = [0,sin(phi),-sin(theta)*cos(phi)]
    G[2,:] = [1,0,cos(theta)]
    return G

def BodyZXZ_Gb(phi, theta, psi):
    """ Angular velocity matrix such that omega_local = Gb theta_dot for rotation of type Body and order ZXZ"""
    Gb = np.zeros((3,3))
    Gb[0,:] = [sin(psi)*sin(theta),cos(psi),0]
    Gb[1,:] = [sin(theta)*cos(psi),-sin(psi),0]
    Gb[2,:] = [cos(theta),0,1]
    return Gb

def BodyZXZ_Ginv(phi, theta, psi):
    """ Inverse of "G" matrix, such that  theta_dot = Ginv omega_global, for rotation of type Body and order ZXZ"""
    Ginv = np.zeros((3,3))
    Ginv[0,:] = [-sin(phi)/tan(theta),cos(phi)/tan(theta),1]
    Ginv[1,:] = [cos(phi),sin(phi),0]
    Ginv[2,:] = [sin(phi)/sin(theta),-cos(phi)/sin(theta),0]
    return Ginv

def BodyZXZ_Gbinv(phi, theta, psi):
    """ Inverse of "Gb" matrix, such that  theta_dot = Gbinv omega_local, for rotation of type Body and order ZXZ"""
    Gbinv = np.zeros((3,3))
    Gbinv[0,:] = [sin(psi)/sin(theta),cos(psi)/sin(theta),0]
    Gbinv[1,:] = [cos(psi),-sin(psi),0]
    Gbinv[2,:] = [-sin(psi)/tan(theta),-cos(psi)/tan(theta),1]
    return Gbinv


def BodyZXZ_toEuler(phi, theta, psi):
    """ 
    Shabana 2.106
    """
    #e0,e1,e2,e2 = EulerP_fromA(BodyZXZ_A(phi,theta,psi)) # Alternative, keep me
    e0 = cos(theta/2)* cos((phi+psi)/2)
    e1 = sin(theta/2)* cos((phi-psi)/2)
    e2 = sin(theta/2)* sin((phi-psi)/2)
    e3 = cos(theta/2)* sin((phi+psi)/2)
    return e0,e1,e2,e3


# --------------------------------------------------------------------------------}
# --- Bryant angles XYZ
# --------------------------------------------------------------------------------{
def BodyXYZ_A(phi_x, phi_y, phi_z):
    """ Transformation matrix body to global for rotation of type Body and order XYZ"""
    A = np.zeros((3,3))
    A[0,:] = [cos(phi_y)*cos(phi_z),-sin(phi_z)*cos(phi_y),sin(phi_y)]
    A[1,:] = [sin(phi_x)*sin(phi_y)*cos(phi_z)+sin(phi_z)*cos(phi_x),-sin(phi_x)*sin(phi_y)*sin(phi_z)+cos(phi_x)*cos(phi_z),-sin(phi_x)*cos(phi_y)]
    A[2,:] = [sin(phi_x)*sin(phi_z)-sin(phi_y)*cos(phi_x)*cos(phi_z),sin(phi_x)*cos(phi_z)+sin(phi_y)*sin(phi_z)*cos(phi_x),cos(phi_x)*cos(phi_y)]
    return A

def BodyXYZ_G(phi_x, phi_y, phi_z):
    """ Angular velocity matrix such that omega_global = G theta_dot for rotation of type Body and order XYZ"""
    G = np.zeros((3,3))
    G[0,:] = [1,0,sin(phi_y)]
    G[1,:] = [0,cos(phi_x),-sin(phi_x)*cos(phi_y)]
    G[2,:] = [0,sin(phi_x),cos(phi_x)*cos(phi_y)]
    return G

def BodyXYZ_Gb(phi_x, phi_y, phi_z):
    """ Angular velocity matrix such that omega_local = Gb theta_dot for rotation of type Body and order XYZ"""
    Gb = np.zeros((3,3))
    Gb[0,:] = [cos(phi_y)*cos(phi_z),sin(phi_z),0]
    Gb[1,:] = [-sin(phi_z)*cos(phi_y),cos(phi_z),0]
    Gb[2,:] = [sin(phi_y),0,1]
    return Gb

def BodyXYZ_Ginv(phi_x, phi_y, phi_z):
    """ Inverse of "G" matrix, such that  theta_dot = Ginv omega_global, for rotation of type Body and order XYZ"""
    Ginv = np.zeros((3,3))
    Ginv[0,:] = [1,sin(phi_x)*tan(phi_y),-cos(phi_x)*tan(phi_y)]
    Ginv[1,:] = [0,cos(phi_x),sin(phi_x)]
    Ginv[2,:] = [0,-sin(phi_x)/cos(phi_y),cos(phi_x)/cos(phi_y)]
    return Ginv

def BodyXYZ_Gbinv(phi_x, phi_y, phi_z):
    """ Inverse of "Gb" matrix, such that  theta_dot = Gbinv omega_local, for rotation of type Body and order XYZ"""
    Gbinv = np.zeros((3,3))
    Gbinv[0,:] = [cos(phi_z)/cos(phi_y),-sin(phi_z)/cos(phi_y),0]
    Gbinv[1,:] = [sin(phi_z),cos(phi_z),0]
    Gbinv[2,:] = [-cos(phi_z)*tan(phi_y),sin(phi_z)*tan(phi_y),1]
    return Gbinv

def BodyXYZ_toEuler(phi_x, phi_y, phi_z):
    """ 
    Convert to Bryant angles to Euler parameters
    """
    return EulerP_fromA(BodyXYZ_A(phi_x,phi_y,phi_z))

# --------------------------------------------------------------------------------}
# --- Angles ZYX
# --------------------------------------------------------------------------------{
def BodyZYX_A(phi_x, phi_y, phi_z):
    """ Transformation matrix body to global for rotation of type Body and order ZYX"""
    A = np.zeros((3,3))
    A[0,:] = [cos(phi_y)*cos(phi_z),sin(phi_x)*sin(phi_y)*cos(phi_z)-sin(phi_z)*cos(phi_x),sin(phi_x)*sin(phi_z)+sin(phi_y)*cos(phi_x)*cos(phi_z)]
    A[1,:] = [sin(phi_z)*cos(phi_y),sin(phi_x)*sin(phi_y)*sin(phi_z)+cos(phi_x)*cos(phi_z),-sin(phi_x)*cos(phi_z)+sin(phi_y)*sin(phi_z)*cos(phi_x)]
    A[2,:] = [-sin(phi_y),sin(phi_x)*cos(phi_y),cos(phi_x)*cos(phi_y)]
    return A

def BodyZYX_G(phi_x, phi_y, phi_z):
    """ Angular velocity matrix such that omega_global = G theta_dot for rotation of type Body and order ZYX"""
    G = np.zeros((3,3))
    G[0,:] = [cos(phi_y)*cos(phi_z),-sin(phi_z),0]
    G[1,:] = [sin(phi_z)*cos(phi_y),cos(phi_z),0]
    G[2,:] = [-sin(phi_y),0,1]
    return G

def BodyZYX_Gb(phi_x, phi_y, phi_z):
    """ Angular velocity matrix such that omega_local = Gb theta_dot for rotation of type Body and order ZYX"""
    Gb = np.zeros((3,3))
    Gb[0,:] = [1,0,-sin(phi_y)]
    Gb[1,:] = [0,cos(phi_x),sin(phi_x)*cos(phi_y)]
    Gb[2,:] = [0,-sin(phi_x),cos(phi_x)*cos(phi_y)]
    return Gb

def BodyZYX_Ginv(phi_x, phi_y, phi_z):
    """ Inverse of "G" matrix, such that  theta_dot = Ginv omega_global, for rotation of type Body and order ZYX"""
    Ginv = np.zeros((3,3))
    Ginv[0,:] = [cos(phi_z)/cos(phi_y),sin(phi_z)/cos(phi_y),0]
    Ginv[1,:] = [-sin(phi_z),cos(phi_z),0]
    Ginv[2,:] = [cos(phi_z)*tan(phi_y),sin(phi_z)*tan(phi_y),1]
    return Ginv

def BodyZYX_Gbinv(phi_x, phi_y, phi_z):
    """ Inverse of "Gb" matrix, such that  theta_dot = Gbinv omega_local, for rotation of type Body and order ZYX"""
    Gbinv = np.zeros((3,3))
    Gbinv[0,:] = [1,sin(phi_x)*tan(phi_y),cos(phi_x)*tan(phi_y)]
    Gbinv[1,:] = [0,cos(phi_x),-sin(phi_x)]
    Gbinv[2,:] = [0,sin(phi_x)/cos(phi_y),cos(phi_x)/cos(phi_y)]
    return Gbinv

def BodyZYX_toEuler(phi_x, phi_y, phi_z):
    """ 
    Convert to Bryant angles to Euler parameters
    """
    return EulerP_fromA(BodyXYZ_A(phi_x, phi_y, phi_z))




# --------------------------------------------------------------------------------}
# --- Rodriguez (OpenFAST functions)
# --------------------------------------------------------------------------------{
def Rodriguez_A(a):
    """
    ! calculates rotation matrix R to rotate unit vertical vector to direction of input vector `a`
    a(3): input vector
    R(3,3): rotation matrix from Rodrigues's rotation formula
    """
    factor = a.dot(a)
    if factor==0:
        return np.eye(3) # Return the identity if the vector is zero
    if a[0]==0 and a[1]==0:   # return identity if vertical
        R = np.eye(3)
        if a[2] < 0:
            R = -R
    else:  
        R = np.zeros((3,3))
        vec = a/np.sqrt(factor) # normalize a
        vec[2] += 1
        factor = 2.0/(vec.dot(vec))
        R[0,0] = factor*vec[0]*vec[0] - 1
        R[0,1] = factor*vec[0]*vec[1]
        R[0,2] = factor*vec[0]*vec[2]
        R[1,0] = factor*vec[1]*vec[0]
        R[1,1] = factor*vec[1]*vec[1] - 1
        R[1,2] = factor*vec[1]*vec[2]
        R[2,0] = factor*vec[2]*vec[0]
        R[2,1] = factor*vec[2]*vec[1]
        R[2,2] = factor*vec[2]*vec[2] - 1


# --------------------------------------------------------------------------------}
# --- Small Rotations
# --------------------------------------------------------------------------------{
def smallRot_A(tx, ty, tz):
    """ 
    small angle approximation: I+skew(theta)
    from body to global
    """
    R=np.zeros((3,3))
    R[0,:] =( 1., -tz,  ty )
    R[1,:] =( tz,  1., -tx )
    R[2,:] =(-ty,  tx,  1. )
    return R


def smallRot_OF(theta1, theta2, theta3):
    r"""

    NOTE: this is a DCM from global to body!!!

    !>  This routine computes the 3x3 transformation matrix, \f$TransMat\f$,
    !!   to a coordinate system \f$x\f$ (with orthogonal axes \f$x_1, x_2, x_3\f$)
    !!   resulting from three rotations (\f$\theta_1\f$, \f$\theta_2\f$, \f$\theta_3\f$) about the
    !!   orthogonal axes (\f$X_1, X_2, X_3\f$) of coordinate system \f$X\f$.  All angles
    !!   are assummed to be small, as such, the order of rotations does
    !!   not matter and Euler angles do not need to be used.  This routine
    !!   is used to compute the transformation matrix (\f$TransMat\f$) between
    !!   undeflected (\f$X\f$) and deflected (\f$x\f$) coordinate systems.  In matrix
    !!   form:
    !! \f{equation}{   
    !!   \left\{ \begin{matrix} x_1 \\ x_2 \\ x_3 \end{matrix} \right\} =
    !!   \left[ TransMat(\theta_1, \theta_2, \theta_3) \right]
    !!   \left\{ \begin{matrix} X_1 \\ X_2 \\ X_3 \end{matrix} \right\}   
    !! \f}
    !!
    !! The transformation matrix, \f$TransMat\f$, is the closest orthonormal
    !!   matrix to the nonorthonormal, but skew-symmetric, Bernoulli-Euler
    !!   matrix:
    !! \f{equation}{   A =
    !!   \begin{bmatrix}    1      &  \theta_3 & -\theta_2 \\
    !!                   -\theta_3 &  1        &  \theta_1 \\
    !!                    \theta_2 & -\theta_1 &  1 \end{bmatrix}   
    !! \f}
    !!   In the Frobenius Norm sense, the closest orthornormal matrix is:
    !!      \f$TransMat = U V^T\f$,
    !!   where the columns of \f$U\f$ contain the eigenvectors of\f$ AA^T\f$ and the
    !!   columns of \f$V\f$ contain the eigenvectors of \f$A^TA\f$ (note that \f$^T\f$ = transpose).
    !!   This result comes directly from the Singular Value Decomposition
    !!   (SVD) of \f$A = USV^T\f$ where \f$S\f$ is a diagonal matrix containing the
    !!   singular values of \f$A\f$, which are sqrt( eigenvalues of \f$AA^T\f$ ) =
    !!   sqrt( eigenvalues of \f$A^TA\f$ ).

    INPUTS:
      - Theta1        : \f$\theta_1\f$: the small rotation about \f$X_1\f$, (rad).
      - Theta2        : \f$\theta_2\f$: the small rotation about \f$X_2\f$, (rad).
      - Theta3        : \f$\theta_3\f$: the small rotation about \f$X_3\f$, (rad).
    OUTPUTS:
      - TransMat (3,3): The resulting transformation matrix from \f$X\f$ to \f$x\f$, (-).
    """
    LrgAngle  = 0.4 # Threshold for when a small angle becomes large (about 23deg).  This comes from: COS(SmllAngle) ~ 1/SQRT( 1 + SmllAngle^2 ) and SIN(SmllAngle) ~ SmllAngle/SQRT( 1 + SmllAngle^2 ) results in ~5% error when SmllAngle = 0.4rad.
    if (abs(theta1) > LrgAngle ) or ( abs(theta2) > LrgAngle ) or ( abs(theta3) > LrgAngle ):
        raise Exception('Small angle assumption violated in SUBROUTINE SmllRotTrans() due to a large rotation')
    # Compute some intermediate results:
    theta11      = theta1*theta1
    theta22      = theta2*theta2
    theta33      = theta3*theta3
    SqrdSum      = theta11 + theta22 + theta33
    SQRT1SqrdSum = np.sqrt( 1. + SqrdSum )
    ComDenom     = SqrdSum*SQRT1SqrdSum
    theta12S     = theta1*theta2*(SQRT1SqrdSum - 1.)
    theta13S     = theta1*theta3*(SQRT1SqrdSum - 1.)
    theta23S     = theta2*theta3*(SQRT1SqrdSum - 1.)
    #Define the transformation matrix:
    TransMat=np.zeros((3,3))
    if ( ComDenom == 0.0 ):
        TransMat[0,:] = ( 1., 0., 0.)
        TransMat[1,:] = ( 0., 1., 0.)
        TransMat[2,:] = ( 0., 0., 1.)
    else:
        TransMat[0,0] = ( theta11*SQRT1SqrdSum + theta22              + theta33              )/ComDenom
        TransMat[1,1] = ( theta11              + theta22*SQRT1SqrdSum + theta33              )/ComDenom
        TransMat[2,2] = ( theta11              + theta22              + theta33*SQRT1SqrdSum )/ComDenom
        TransMat[0,1] = (  theta3*SqrdSum + theta12S )/ComDenom
        TransMat[1,0] = ( -theta3*SqrdSum + theta12S )/ComDenom
        TransMat[0,2] = ( -theta2*SqrdSum + theta13S )/ComDenom
        TransMat[2,0] = (  theta2*SqrdSum + theta13S )/ComDenom
        TransMat[1,2] = (  theta1*SqrdSum + theta23S )/ComDenom
        TransMat[2,1] = ( -theta1*SqrdSum + theta23S )/ComDenom
    return TransMat

