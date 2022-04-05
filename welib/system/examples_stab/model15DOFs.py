
import numpy as np


def bendingMode(s, ModeNr = None): 
    """
    Compute clamped-free bending modes for Modenr<5
     - s: spanwise coordinate  (in meters)

     NOTE: not normlized to unity!!!
    """
    temp = np.array([1.8751,4.6941,7.8548,10.9955])
    lambda_ = temp[ModeNr-1]
    F = np.sinh(lambda_) + np.sin(lambda_)
    G = np.cosh(lambda_) + np.cos(lambda_)
    R = s[-1]
    c1 = G / F # 0.7340955 for mode 1
    return np.cosh(s * lambda_ / R) - np.cos(s * lambda_ / R) - c1 * (np.sinh(s * lambda_ / R) - np.sin(s * lambda_ / R))


def torsionalMode(s, ModeNr = 1): 
    """
    Returns analytical torsional mode for a beam with constant properties
     - s: spanwise coordinate  (in meters)
    
    """
    return np.sin((2 * ModeNr - 1) * np.pi * s / (2 * s[-1]))



def diagonalStructuralMatrices(p,Omega,lambda1,c1): 
    """ 
    lambda1,c1 coefficients for analytical mode shape TODO replace by num integration of phi, Only used for Mb and Kb
    """
    # Turbine parameters
    R    = p['R']
    c    = p['c']
    a_cg = p['acg']
    EIx  = p['EIx']
    EIy  = p['EIy']
    GK   = p['GK']
    m    = p['m']
    J    = p['J']
    Ls   = p['Ls']
    Gs   = p['Gs']
    Mn   = p['Mn']
    Ix   = p['Ix']
    Iy   = p['Iy']
    Iz   = p['Iz']
    M    = p['Mtot']
    Kx   = p['Kx']
    Ky   = p['Ky']
    Gx   = p['Gx']
    Gy   = p['Gy']
    Gz   = p['Gz']
    gxy  = p['gxy']
    # Blade mass matrix
#   flap, edge torsion
    Mb = np.zeros((3,3))
    Mb[0,0] = - m * R * (4.0 * np.exp((2 * lambda1)) * c1 + 1.0 - c1 ** 2 * np.exp((4 * lambda1)) + 4.0 * c1 ** 2 * np.cos(lambda1) * np.exp(lambda1) + 4.0 * np.cos(lambda1) * np.exp((3 * lambda1)) + 4.0 * c1 ** 2 * np.sin(lambda1) * np.exp(lambda1) + 4.0 * np.sin(lambda1) * np.exp((3 * lambda1)) + 4.0 * c1 ** 2 * np.exp((2 * lambda1)) * np.cos(lambda1) * np.sin(lambda1) + 2.0 * c1 * np.exp((4 * lambda1)) - np.exp((4 * lambda1)) - 4.0 * c1 ** 2 * np.exp((3 * lambda1)) * np.cos(lambda1) - 8.0 * lambda1 * np.exp((2 * lambda1)) + 4.0 * c1 ** 2 * np.exp((3 * lambda1)) * np.sin(lambda1) - 4.0 * np.cos(lambda1) * np.sin(lambda1) * np.exp((2 * lambda1)) + 8.0 * c1 * np.sin(lambda1) * np.exp(lambda1) - 8.0 * c1 * np.cos(lambda1) ** 2 * np.exp((2 * lambda1)) + c1 ** 2 - 4.0 * np.cos(lambda1) * np.exp(lambda1) + 4.0 * np.sin(lambda1) * np.exp(lambda1) - 8.0 * c1 * np.exp((3 * lambda1)) * np.sin(lambda1) + 2.0 * c1) / lambda1 * np.exp(- (2 * lambda1)) / 8.0
    Mb[0,1] = 0
    Mb[0,2] = 2.0 * lambda1 * R * a_cg * m / (- np.pi ** 4 + (16 * lambda1 ** 4)) * (- 8.0 * np.pi * lambda1 + 2.0 * np.cos(lambda1) * c1 * np.pi ** 2 - np.exp(- lambda1) * c1 * np.pi ** 2 + 2.0 * np.sin(lambda1) * np.pi ** 2 - 4.0 * (lambda1 ** 2) * np.exp(lambda1) + 4.0 * np.exp(- lambda1) * (lambda1 ** 2) * c1 + 8.0 * np.sin(lambda1) * (lambda1 ** 2) + np.pi ** 2 * np.exp(lambda1) - c1 * np.pi ** 2 * np.exp(lambda1) + 8.0 * np.cos(lambda1) * (lambda1 ** 2) * c1 - np.exp(- lambda1) * np.pi ** 2 + 4.0 * (lambda1 ** 2) * c1 * np.exp(lambda1) + 4.0 * np.exp(- lambda1) * (lambda1 ** 2))
    Mb[1,0] = 0
    Mb[1,1] = Mb[0,0]
    Mb[1,2] = 0
    Mb[2,0] = Mb[0,2]
    Mb[2,1] = 0
    Mb[2,2] = R * m * a_cg ** 2 / 2.0 + R * J / 2.0
    # Blade stiffness matrix
    Kb = np.zeros((3,3))
    Kb[0,0] = np.exp(- (2 * lambda1)) * (- (12 * EIx * lambda1 ** 4) - (12 * EIx * lambda1 ** 4 * c1 ** 2) + 96.0 * EIx * (lambda1 ** 5) * np.exp((2 * lambda1)) - (24 * EIx * lambda1 ** 4 * c1) - 48.0 * EIx * (lambda1 ** 4) * np.cos(lambda1) * np.exp(lambda1) + 48.0 * EIx * (lambda1 ** 4) * np.sin(lambda1) * np.exp(lambda1) + 48.0 * EIx * (lambda1 ** 4) * np.cos(lambda1) * np.sin(lambda1) * np.exp((2 * lambda1)) + 96.0 * EIx * (lambda1 ** 4) * c1 * np.cos(lambda1) ** 2 * np.exp((2 * lambda1)) + 96.0 * EIx * (lambda1 ** 4) * c1 * np.sin(lambda1) * np.exp(lambda1) - 48.0 * EIx * (lambda1 ** 4) * (c1 ** 2) * np.exp((2 * lambda1)) * np.cos(lambda1) * np.sin(lambda1) - 24.0 * EIx * (lambda1 ** 4) * c1 * np.exp((4 * lambda1)) - 48.0 * EIx * (lambda1 ** 4) * np.exp((2 * lambda1)) * c1 + 12.0 * EIx * (lambda1 ** 4) * (c1 ** 2) * np.exp((4 * lambda1)) + 48.0 * EIx * (lambda1 ** 4) * (c1 ** 2) * np.sin(lambda1) * np.exp(lambda1) + 12.0 * EIx * (lambda1 ** 4) * np.exp((4 * lambda1)) + 48.0 * EIx * (lambda1 ** 4) * np.cos(lambda1) * np.exp((3 * lambda1)) + 48.0 * EIx * (lambda1 ** 4) * np.sin(lambda1) * np.exp((3 * lambda1)) - 48.0 * EIx * (lambda1 ** 4) * (c1 ** 2) * np.exp((3 * lambda1)) * np.cos(lambda1) + 48.0 * EIx * (lambda1 ** 4) * (c1 ** 2) * np.cos(lambda1) * np.exp(lambda1) - 96.0 * EIx * (lambda1 ** 4) * c1 * np.exp((3 * lambda1)) * np.sin(lambda1) + 48.0 * EIx * (lambda1 ** 4) * (c1 ** 2) * np.exp((3 * lambda1)) * np.sin(lambda1) + (6 * m * Omega ** 2 * R ** 4 * c1 ** 2 * lambda1) + (12 * m * Omega ** 2 * R ** 4 * c1 * lambda1) + (3 * m * Omega ** 2 * R ** 4) + (3 * m * Omega ** 2 * R ** 4 * c1 ** 2) + (6 * m * Omega ** 2 * R ** 4 * c1) + (6 * m * Omega ** 2 * R ** 4 * lambda1) - 3.0 * m * (Omega ** 2) * (R ** 4) * np.exp((4 * lambda1)) + 32.0 * m * (Omega ** 2) * (R ** 4) * (lambda1 ** 3) * (c1 ** 2) * np.exp((2 * lambda1)) + 48.0 * m * (Omega ** 2) * (R ** 4) * c1 * np.sin(lambda1) * np.exp((3 * lambda1)) * lambda1 - 48.0 * m * (Omega ** 2) * (R ** 4) * np.cos(lambda1) * lambda1 * np.exp(lambda1) - 48.0 * m * (Omega ** 2) * (R ** 4) * (lambda1 ** 2) * c1 * np.exp((2 * lambda1)) + 6.0 * m * (Omega ** 2) * (R ** 4) * np.exp((4 * lambda1)) * (c1 ** 2) * lambda1 + 24.0 * m * (Omega ** 2) * (R ** 4) * np.cos(lambda1) ** 2 * lambda1 * np.exp((2 * lambda1)) - 48.0 * m * (Omega ** 2) * (R ** 4) * np.cos(lambda1) * np.exp((3 * lambda1)) * lambda1 - 12.0 * m * (Omega ** 2) * (R ** 4) * np.exp((4 * lambda1)) * c1 * lambda1 + 12.0 * m * (Omega ** 2) * (R ** 4) * (c1 ** 2) * np.exp((2 * lambda1)) * lambda1 - 48.0 * m * (Omega ** 2) * (R ** 4) * c1 * np.cos(lambda1) * np.sin(lambda1) * lambda1 * np.exp((2 * lambda1)) + 48.0 * m * (Omega ** 2) * (R ** 4) * c1 * np.sin(lambda1) * lambda1 * np.exp(lambda1) + 48.0 * m * (Omega ** 2) * (R ** 4) * (c1 ** 2) * np.sin(lambda1) * lambda1 * np.exp(lambda1) - 48.0 * m * (Omega ** 2) * (R ** 4) * c1 * np.cos(lambda1) * lambda1 * np.exp(lambda1) + 6.0 * m * (Omega ** 2) * (R ** 4) * lambda1 * np.exp((4 * lambda1)) + 48.0 * m * (Omega ** 2) * (R ** 4) * np.exp((3 * lambda1)) * c1 * np.cos(lambda1) * lambda1 - 24.0 * m * (Omega ** 2) * (R ** 4) * (c1 ** 2) * np.exp((2 * lambda1)) * lambda1 * np.cos(lambda1) ** 2 - 48.0 * m * (Omega ** 2) * (R ** 4) * (c1 ** 2) * np.sin(lambda1) * np.exp((3 * lambda1)) * lambda1 + 12.0 * m * (Omega ** 2) * (R ** 4) * np.exp((2 * lambda1)) * c1 - 3.0 * m * (Omega ** 2) * (R ** 4) * (c1 ** 2) * np.exp((4 * lambda1)) + 24.0 * m * (Omega ** 2) * (R ** 4) * np.exp(lambda1) * (c1 ** 2) * np.cos(lambda1) + 24.0 * m * (Omega ** 2) * (R ** 4) * np.cos(lambda1) * np.exp((3 * lambda1)) + 24.0 * m * (Omega ** 2) * (R ** 4) * (c1 ** 2) * np.sin(lambda1) * np.exp(lambda1) + 24.0 * m * (Omega ** 2) * (R ** 4) * np.sin(lambda1) * np.exp((3 * lambda1)) + 12.0 * m * (Omega ** 2) * (R ** 4) * (c1 ** 2) * np.exp((2 * lambda1)) * np.cos(lambda1) * np.sin(lambda1) + 6.0 * m * (Omega ** 2) * (R ** 4) * c1 * np.exp((4 * lambda1)) - 24.0 * m * (Omega ** 2) * (R ** 4) * (c1 ** 2) * np.exp((3 * lambda1)) * np.cos(lambda1) - 12.0 * m * (Omega ** 2) * (R ** 4) * lambda1 * np.exp((2 * lambda1)) + 24.0 * m * (Omega ** 2) * (R ** 4) * (c1 ** 2) * np.exp((3 * lambda1)) * np.sin(lambda1) - 12.0 * m * (Omega ** 2) * (R ** 4) * np.cos(lambda1) * np.sin(lambda1) * np.exp((2 * lambda1)) + 48.0 * m * (Omega ** 2) * (R ** 4) * c1 * np.sin(lambda1) * np.exp(lambda1) - 24.0 * m * (Omega ** 2) * (R ** 4) * c1 * np.cos(lambda1) ** 2 * np.exp((2 * lambda1)) - 24.0 * m * (Omega ** 2) * (R ** 4) * np.cos(lambda1) * np.exp(lambda1) + 24.0 * m * (Omega ** 2) * (R ** 4) * np.exp(lambda1) * np.sin(lambda1) - 48.0 * m * (Omega ** 2) * (R ** 4) * c1 * np.exp((3 * lambda1)) * np.sin(lambda1)) / (R ** 3) / lambda1 / 96.0
    Kb[1,1] = np.exp(- (2 * lambda1)) * ((6 * m * Omega ** 2 * R ** 4 * c1 ** 2 * lambda1) + (12 * m * Omega ** 2 * R ** 4 * c1 * lambda1) + (15 * m * Omega ** 2 * R ** 4) + (15 * m * Omega ** 2 * R ** 4 * c1 ** 2) + (30 * m * Omega ** 2 * R ** 4 * c1) + (6 * m * Omega ** 2 * R ** 4 * lambda1) - 15.0 * m * (Omega ** 2) * (R ** 4) * np.exp((4 * lambda1)) - 48.0 * EIy * (lambda1 ** 4) * np.exp((2 * lambda1)) * c1 + 12.0 * EIy * (lambda1 ** 4) * (c1 ** 2) * np.exp((4 * lambda1)) - 48.0 * EIy * (lambda1 ** 4) * np.cos(lambda1) * np.exp(lambda1) + 48.0 * EIy * (lambda1 ** 4) * np.exp(lambda1) * np.sin(lambda1) - 24.0 * EIy * (lambda1 ** 4) * c1 * np.exp((4 * lambda1)) + 48.0 * EIy * (lambda1 ** 4) * np.cos(lambda1) * np.exp((3 * lambda1)) + 48.0 * EIy * (lambda1 ** 4) * np.sin(lambda1) * np.exp((3 * lambda1)) - 12.0 * EIy * (lambda1 ** 4) - 12.0 * EIy * (lambda1 ** 4) * (c1 ** 2) - 24.0 * EIy * (lambda1 ** 4) * c1 + 96.0 * EIy * (lambda1 ** 5) * np.exp((2 * lambda1)) + 12.0 * EIy * (lambda1 ** 4) * np.exp((4 * lambda1)) - 48.0 * EIy * (lambda1 ** 4) * (c1 ** 2) * np.exp((2 * lambda1)) * np.cos(lambda1) * np.sin(lambda1) + 32.0 * m * (Omega ** 2) * (R ** 4) * (lambda1 ** 3) * (c1 ** 2) * np.exp((2 * lambda1)) - 48.0 * EIy * (lambda1 ** 4) * (c1 ** 2) * np.exp((3 * lambda1)) * np.cos(lambda1) + 48.0 * EIy * (lambda1 ** 4) * np.exp(lambda1) * (c1 ** 2) * np.cos(lambda1) - 96.0 * EIy * (lambda1 ** 4) * c1 * np.exp((3 * lambda1)) * np.sin(lambda1) + 48.0 * EIy * (lambda1 ** 4) * (c1 ** 2) * np.exp((3 * lambda1)) * np.sin(lambda1) + 48.0 * m * (Omega ** 2) * (R ** 4) * c1 * np.sin(lambda1) * np.exp((3 * lambda1)) * lambda1 - 48.0 * m * (Omega ** 2) * (R ** 4) * np.cos(lambda1) * lambda1 * np.exp(lambda1) - 48.0 * m * (Omega ** 2) * (R ** 4) * (lambda1 ** 2) * c1 * np.exp((2 * lambda1)) + 6.0 * m * (Omega ** 2) * (R ** 4) * np.exp((4 * lambda1)) * (c1 ** 2) * lambda1 + 24.0 * m * (Omega ** 2) * (R ** 4) * np.cos(lambda1) ** 2 * lambda1 * np.exp((2 * lambda1)) - 48.0 * m * (Omega ** 2) * (R ** 4) * np.cos(lambda1) * np.exp((3 * lambda1)) * lambda1 - 12.0 * m * (Omega ** 2) * (R ** 4) * np.exp((4 * lambda1)) * c1 * lambda1 + 12.0 * m * (Omega ** 2) * (R ** 4) * (c1 ** 2) * np.exp((2 * lambda1)) * lambda1 - 48.0 * m * (Omega ** 2) * (R ** 4) * c1 * np.cos(lambda1) * np.sin(lambda1) * lambda1 * np.exp((2 * lambda1)) + 48.0 * m * (Omega ** 2) * (R ** 4) * c1 * np.sin(lambda1) * lambda1 * np.exp(lambda1) + 48.0 * m * (Omega ** 2) * (R ** 4) * (c1 ** 2) * np.sin(lambda1) * lambda1 * np.exp(lambda1) - 48.0 * m * (Omega ** 2) * (R ** 4) * c1 * np.cos(lambda1) * lambda1 * np.exp(lambda1) + 6.0 * m * (Omega ** 2) * (R ** 4) * lambda1 * np.exp((4 * lambda1)) + 48.0 * m * (Omega ** 2) * (R ** 4) * np.exp((3 * lambda1)) * c1 * np.cos(lambda1) * lambda1 - 24.0 * m * (Omega ** 2) * (R ** 4) * (c1 ** 2) * np.exp((2 * lambda1)) * lambda1 * np.cos(lambda1) ** 2 - 48.0 * m * (Omega ** 2) * (R ** 4) * (c1 ** 2) * np.sin(lambda1) * np.exp((3 * lambda1)) * lambda1 + 60.0 * m * (Omega ** 2) * (R ** 4) * np.exp((2 * lambda1)) * c1 - 15.0 * m * (Omega ** 2) * (R ** 4) * (c1 ** 2) * np.exp((4 * lambda1)) + 72.0 * m * (Omega ** 2) * (R ** 4) * np.exp(lambda1) * (c1 ** 2) * np.cos(lambda1) + 72.0 * m * (Omega ** 2) * (R ** 4) * np.cos(lambda1) * np.exp((3 * lambda1)) + 72.0 * m * (Omega ** 2) * (R ** 4) * (c1 ** 2) * np.sin(lambda1) * np.exp(lambda1) + 72.0 * m * (Omega ** 2) * (R ** 4) * np.sin(lambda1) * np.exp((3 * lambda1)) + 60.0 * m * (Omega ** 2) * (R ** 4) * (c1 ** 2) * np.exp((2 * lambda1)) * np.cos(lambda1) * np.sin(lambda1) + 30.0 * m * (Omega ** 2) * (R ** 4) * c1 * np.exp((4 * lambda1)) - 72.0 * m * (Omega ** 2) * (R ** 4) * (c1 ** 2) * np.exp((3 * lambda1)) * np.cos(lambda1) - 108.0 * m * (Omega ** 2) * (R ** 4) * lambda1 * np.exp((2 * lambda1)) + 72.0 * m * (Omega ** 2) * (R ** 4) * (c1 ** 2) * np.exp((3 * lambda1)) * np.sin(lambda1) - 60.0 * m * (Omega ** 2) * (R ** 4) * np.cos(lambda1) * np.sin(lambda1) * np.exp((2 * lambda1)) + 144.0 * m * (Omega ** 2) * (R ** 4) * c1 * np.sin(lambda1) * np.exp(lambda1) - 120.0 * m * (Omega ** 2) * (R ** 4) * c1 * np.cos(lambda1) ** 2 * np.exp((2 * lambda1)) - 72.0 * m * (Omega ** 2) * (R ** 4) * np.cos(lambda1) * np.exp(lambda1) + 72.0 * m * (Omega ** 2) * (R ** 4) * np.exp(lambda1) * np.sin(lambda1) - 144.0 * m * (Omega ** 2) * (R ** 4) * c1 * np.exp((3 * lambda1)) * np.sin(lambda1) + 48.0 * EIy * (lambda1 ** 4) * np.cos(lambda1) * np.sin(lambda1) * np.exp((2 * lambda1)) + 96.0 * EIy * (lambda1 ** 4) * c1 * np.cos(lambda1) ** 2 * np.exp((2 * lambda1)) + 96.0 * EIy * (lambda1 ** 4) * c1 * np.sin(lambda1) * np.exp(lambda1) + 48.0 * EIy * (lambda1 ** 4) * (c1 ** 2) * np.sin(lambda1) * np.exp(lambda1)) / (R ** 3) / lambda1 / 96.0
    Kb[2,2] = (m * Omega ** 2 * a_cg ** 2 * R) / 2.0 + GK * np.pi ** 2 / R / 8.0
    # Blade gyroscopic matrix
    Gb = np.zeros((3,3))
    # Nacelle/tower top mass matrix
    Mt = np.array([
        [3*m*R+M  , 0       , 0                                       , 3*m*R*Ls                                , 0                       , 0]                     , 
        [0        , 3*m*R+M , 0                                       , 0                                       , 0                       , 0]                     , 
        [0        , 0       , 0.5*m*R*(R**2+3*a_cg**2)+3*m*R*Ls**2+Ix , 0                                       , 0                       , 0]                     , 
        [3*m*R*Ls , 0       , 0                                       , 0.5*m*R*(R**2+3*a_cg**2)+3*m*R*Ls**2+Iz , 0                       , 0]                     , 
        [0        , 0       , 0                                       , 0                                       , m*R**3+3*m*R*a_cg**2+Iy , m*R*(R**2+3*a_cg**2)]  , 
        [0        , 0       , 0                                       , 0                                       , m*R*(R**2+3*a_cg**2)    , m*R*(R**2+3*a_cg**2)]])
    # Nacelle/tower top gyroscopic matrix
    Gt = np.array([
        [0 , 0 , 0                          , 0                           , 0 , 0]   , 
        [0 , 0 , 0                          , 0                           , 0 , 0]   , 
        [0 , 0 , 0                          , -m*R*Omega*(R**2+3*a_cg**2) , 0 , 0]   , 
        [0 , 0 , m*R*Omega*(R**2+3*a_cg**2) , 0                           , 0 , 0]   , 
        [0 , 0 , 0                          , 0                           , 0 , 0]   , 
        [0 , 0 , 0                          , 0                           , 0 , 0]])
    # Nacelle/tower top stiffness matrix
    # Gx: spring stiffness for tilt
    # Gz: spring stiffness for yaw 
    #  utx   uty    tilt     yaw   roll  DTTorsion
    Kt = np.array([
        [Kx  , 0     , 0     , 0  , gxy , 0]    , 
        [0   , Ky    , - gxy , 0  , 0   , 0]    , 
        [0   , - gxy , Gx    , 0  , 0   , 0]    , 
        [0   , 0     , 0     , Gz , 0   , 0]    , 
        [gxy , 0     , 0     , 0  , Gy  , 0]    , 
        [0   , 0     , 0     , 0  , 0   , Gs]])
    return Mb,Kb,Gb,Mt,Gt,Kt


def azimuthDependentStructuralMatrix(p,Omega, psi = 0, N=500): 
    """
    N: Number of radial stations used for integration
    """
    # Shape functions used in the modal expansion
    R = p['R']
    z = np.linspace(0,R,N)
    h = z[2] - z[1]
    phi_f = bendingMode(z,1)
    phi_e = phi_f
    phi_t = torsionalMode(z,1)
    # Approximate integrals by trapez formula
    Mtb = 0
    Gbt = 0
    Gtb = 0
    Ktb = 0
    dMtb_old,dGbt_old,dGtb_old,dKtb_old = dMatrix(z[0],phi_f[0],phi_e[0],phi_t[0],psi,p,Omega)
    for i in np.arange(1,N):
        dMtb,dGbt,dGtb,dKtb = dMatrix(z[i],phi_f[i],phi_e[i],phi_t[i],psi,p,Omega)
        Mtb = Mtb + h / 2 * (dMtb_old + dMtb)
        Gtb = Gtb + h / 2 * (dGtb_old + dGtb)
        Gbt = Gbt + h / 2 * (dGbt_old + dGbt)
        Ktb = Ktb + h / 2 * (dKtb_old + dKtb)
        dMtb_old = dMtb
        dGbt_old = dGbt
        dGtb_old = dGtb
        dKtb_old = dKtb 
    return Mtb,Gbt,Gtb,Ktb   

    
def dMatrix(z = None, phi_f = None, phi_e = None, phi_t = None, psi = None, p = None, Omega = None): 
    """
     Kernel of integrals in Mtb, Gbt, Gtb, and Ktb
    """
    # Turbine parameters
    a_cg = p['acg']
    m    = p['m']
    Ls   = p['Ls']
    # Mass coupling between rotor and tower
    Mtb = m * np.array([
                [0                                       , np.cos(psi)*phi_e    , 0                                          ], 
                [phi_f                                   , 0                    , -a_cg*phi_t                                ], 
                [-phi_f*(a_cg*np.sin(psi)+z*np.cos(psi)) , Ls*np.sin(psi)*phi_e , a_cg*phi_t*(a_cg*np.sin(psi)+z*np.cos(psi))], 
                [-phi_f*(a_cg*np.cos(psi)-z*np.sin(psi)) , Ls*np.cos(psi)*phi_e , a_cg*phi_t*(a_cg*np.cos(psi)-z*np.sin(psi))], 
                [0                                       , phi_e*z              , 0                                          ], 
                [0                                       , phi_e*z              , 0]])
    # Gyroscopic coupling from tower to rotor
    Gbt = m * Omega * np.array([
        [0 , 0 , -2*phi_f*(a_cg*np.cos(psi)-z*np.sin(psi))     , 2*phi_f*(a_cg*np.sin(psi)+z*np.cos(psi))       , 0            , 0           ], 
        [0 , 0 , 0                                             , 0                                              , 2*a_cg*phi_e , 2*a_cg*phi_e], 
        [0 , 0 , 2*a_cg*phi_t*(a_cg*np.cos(psi)-z*np.sin(psi)) , -2*a_cg*phi_t*(a_cg*np.sin(psi)+z*np.cos(psi)) , 0            , 0           ]])
    # Gyroscopic coupling from rotor to tower
    Gtb = m * Omega * np.array([
        [0 , - 2 * np.sin(psi) * phi_e      , 0], 
        [0 , 0                              , 0], 
        [0 , 2 * Ls * np.cos(psi) * phi_e   , 0], 
        [0 , - 2 * Ls * np.sin(psi) * phi_e , 0], 
        [0 , - 2 * a_cg * phi_e             , 0], 
        [0 , - 2 * a_cg * phi_e             , 0]])
    # Deflection proportional inertia coupling from rotor to tower
    Ktb = m * Omega ** 2 * np.array([
        [0                                                , - np.cos(psi) * phi_e      , 0                                                    ] , 
        [0                                                , 0                          , 0                                                    ] , 
        [- phi_f * (a_cg * np.sin(psi) + z * np.cos(psi)) , - Ls * np.sin(psi) * phi_e , a_cg * phi_t * (a_cg * np.sin(psi) + z * np.cos(psi))] , 
        [- phi_f * (a_cg * np.cos(psi) - z * np.sin(psi)) , - Ls * np.cos(psi) * phi_e , a_cg * phi_t * (a_cg * np.cos(psi) - z * np.sin(psi))] , 
        [0                                                , 0                          , 0                                                    ] , 
        [0                                                , 0                          , 0]])
    return Mtb,Gbt,Gtb,Ktb



def MCKmat(p, Omega, psi1 = 0): 
    """
    Return full system matrices

    """
    ndof = 15
    nB = 3
    M = np.zeros((ndof,ndof))
    C = np.zeros((ndof,ndof))
    K = np.zeros((ndof,ndof))
    psi2 = psi1 + 2 * np.pi / 3
    psi3 = psi2 + 2 * np.pi / 3
    # Mass matrix blade tower
    Mtb1,Gbt1,Gtb1,Ktb1 = azimuthDependentStructuralMatrix(p,Omega,psi1)
    Mtb2,Gbt2,Gtb2,Ktb2 = azimuthDependentStructuralMatrix(p,Omega,psi2)
    Mtb3,Gbt3,Gtb3,Ktb3 = azimuthDependentStructuralMatrix(p,Omega,psi3)
    # mass matrix blade blade
    lambda1 = 1.8751
    c1 = 0.7341
    Mb,Kb,Gb,Mt,Gt,Kt = diagonalStructuralMatrices(p, Omega, lambda1, c1)
    
    # Zeros matrices
    Z3  = np.zeros((3,3))
    Z36 = np.zeros((3,6))
    # Filling the matrix blocks
    M = np.block([
        [Mb   , Z3   , Z3   , np.transpose(Mtb1)] , 
        [Z3   , Mb   , Z3   , np.transpose(Mtb2)] , 
        [Z3   , Z3   , Mb   , np.transpose(Mtb3)] , 
        [Mtb1 , Mtb2 , Mtb3 , Mt]])
    # Filling the stiffness matrix
    # NOTE: Gtb/=Gbt.T they are not transpose of each other
    C = np.block([
        [Z3   , Z3   , Z3   , Gbt1] , 
        [Z3   , Z3   , Z3   , Gbt2] , 
        [Z3   , Z3   , Z3   , Gbt3] , 
        [Gtb1 , Gtb2 , Gtb3 , Gt]])
    # Filling the stiffness matrix
    K = np.block([
        [Kb   , Z3   , Z3   , Z36]  , 
        [Z3   , Kb   , Z3   , Z36]  , 
        [Z3   , Z3   , Kb   , Z36]  , 
        [Ktb1 , Ktb2 , Ktb3 , Kt]])
    return M,C,K

def Bmat(psi1,Omega): 
    """ 
     See Hansen 2003 Eq. 14 15, 16
    """
    psi2 = psi1 + (2 * np.pi / 3)
    psi3 = psi2 + (2 * np.pi / 3)
    B = np.block([
        [np.eye(3    , 3)  , np.eye(3 , 3) * np.cos(psi1) , np.eye(3 , 3) * np.sin(psi1) , np.zeros((3 , 6))] , 
        [np.eye(3    , 3)  , np.eye(3 , 3) * np.cos(psi2) , np.eye(3 , 3) * np.sin(psi2) , np.zeros((3 , 6))] , 
        [np.eye(3    , 3)  , np.eye(3 , 3) * np.cos(psi3) , np.eye(3 , 3) * np.sin(psi3) , np.zeros((3 , 6))] , 
        [np.zeros((6 , 9)) , np.eye(6 , 6)]]) # TODO Omega t?
    Bdot = np.block([
        [np.zeros((3 , 3)) , Omega * np.eye(3 , 3) * - np.sin(psi1) , Omega * np.eye(3 , 3) * np.cos(psi1) , np.zeros((3 , 6))] , 
        [np.zeros((3 , 3)) , Omega * np.eye(3 , 3) * - np.sin(psi2) , Omega * np.eye(3 , 3) * np.cos(psi2) , np.zeros((3 , 6))] , 
        [np.zeros((3 , 3)) , Omega * np.eye(3 , 3) * - np.sin(psi3) , Omega * np.eye(3 , 3) * np.cos(psi3) , np.zeros((3 , 6))] , 
        [np.zeros((6 , 9)) , np.zeros((6      , 6))]])
    mu = np.block([
        [(1 / 3) * np.eye(3 , 3)  , np.zeros((3        , 3)) , np.zeros((3        , 3)) , np.zeros((3 , 3)) , np.zeros((3 , 3))]  , 
        [np.zeros((3        , 3)) , (2 / 3) * np.eye(3 , 3)  , np.zeros((3        , 3)) , np.zeros((3 , 3)) , np.zeros((3 , 3))]  , 
        [np.zeros((3        , 3)) , np.zeros((3        , 3)) , (2 / 3) * np.eye(3 , 3)  , np.zeros((3 , 3)) , np.zeros((3 , 3))]  ,  # Should be 2/3
        [np.zeros((3        , 3)) , np.zeros((3        , 3)) , np.zeros((3        , 3)) , np.eye(3    , 3)  , np.zeros((3 , 3))]  , 
        [np.zeros((3        , 3)) , np.zeros((3        , 3)) , np.zeros((3        , 3)) , np.zeros((3 , 3)) , np.eye(3    , 3)]])
    R = np.linalg.inv(B).dot(Bdot)
    # invB=inv(B);
    # transB=B';
    # invbt=inv(transB)
    # mu=invB*invbt
    
    return B,mu,Bdot,R

def MbCbKbmat(B, M, C, K, mu): 
    """
     See Hansen 2003 Eq. 16
    """
    MB = mu.dot(np.transpose(B)).dot(M).dot(B)
    CB = mu.dot(np.transpose(B)).dot(C).dot(B)
    KB = mu.dot(np.transpose(B)).dot(K).dot(B)
    return MB,CB,KB

def MCKTildemat(B, MB, CB, KB, mu, R): 
    """ 
     Compute tilde matrices, ie the matrices of the system in z (Coleman coordinates)
     See Hansen 2003 Eq. 16
    """
    MBT = MB
    CBT = 2 * MB.dot(R) + CB
    KBT = MB.dot(R).dot(R) + CB.dot(R) + KB
    return MBT,CBT,KBT

def Coleman2Comp(a0, a1, b1):
    """ """
    A0    = 0.5 * np.sqrt( np.real(a0) ** 2 + np.imag(a0) ** 2)
    ABW   = 0.5 * np.sqrt((np.real(b1) - np.imag(a1)) ** 2 + (np.imag(b1) + np.real(a1)) ** 2)
    AFW   = 0.5 * np.sqrt((np.real(b1) + np.imag(a1)) ** 2 + (np.real(a1) - np.imag(b1)) ** 2)
    phi0  = np.arctan2( np.imag(a0),  np.real(a0))
    phiBW = np.arctan2(np.imag(a1) - np.real(b1),np.imag(b1) + np.real(a1))
    phiFW = np.arctan2(np.real(b1) + np.imag(a1),np.real(a1) - np.imag(b1))
    return A0,ABW,AFW,phi0,phiBW,phiFW
    
def aerodynMatrix(par,Omega,W): 
    # computation parameters
    N = 500
    z = np.linspace(0,par['R'],N)
    h = z[1] - z[0]
    # bending modes
    phi_e = bendingMode(z,1)  
    phi_f = bendingMode(z,1)  
    phi_t = torsionalMode(z,1)
    for i in np.arange(N):
        dAeroMat =dAerodynMatrix(par,z[i],Omega,W,phi_f[i],phi_e[i],phi_t[i])
        if i==0:
            aeroMat= dAeroMat
        else:
            dAeroMat =dAerodynMatrix(par,z[i],Omega,W,phi_f[i],phi_e[i],phi_t[i])
            for k,v in aeroMat.items():
                aeroMat[k] +=  h / 2 * (dAeroMat_old[k] + dAeroMat[k])
        dAeroMat_old = dAeroMat
    return aeroMat
    
    
def dAerodynMatrix(par, z, Omega, W, phi_f, phi_e, phi_t): 
    # aerodynamic paramters
    c   = par['c']
    Ls  = par['Ls']
    rho = par['rho']
    # computer parameters
    lambda_ = z * Omega / W
    alpha0  = np.arctan(1 / lambda_)
    psi0    = alpha0
    
    U0 = np.sqrt(W ** 2 + (z * Omega) ** 2)
    Cl,Cd,dCl,dCd = liftDrag(alpha0)
    CX0 = Cl * np.sin(psi0) - Cd * np.cos(psi0)
    CY0 = Cl * np.cos(psi0) + Cd * np.sin(psi0)
    CXp0 = dCl * np.sin(psi0) - dCd * np.cos(psi0)
    CYp0 = dCl * np.cos(psi0) + dCd * np.sin(psi0)
    dAeroMat=dict()
    # computer matrix
    dAeroMat['Cab'] = 0.25 * c * rho * W * np.array([
        [4*phi_f**2*CY0+2*phi_f**2*CYp0*lambda_-2*phi_f**2*CX0*lambda_         , 2*phi_f*phi_e*CYp0-2*phi_f*phi_e*CX0-4*phi_f*phi_e*CY0*lambda_ , -c*phi_f*phi_t*CYp0*lambda_] , 
        [4*phi_f*phi_e*CX0+2*phi_f*phi_e*CXp0*lambda_+2*phi_f*phi_e*CY0*lambda_ , 2*phi_e**2*CXp0+2*phi_e**2*CY0-4*phi_e**2*CX0*lambda_          , -c*phi_e*phi_t*CXp0*lambda_] , 
        [0                                                                     , 0                                                              , 0]])
    dAeroMat['Kab'] = 0.5 * c * rho * U0 ** 2 * np.array([
        [0 , 0 , -phi_t*phi_f*CYp0] , 
        [0 , 0 , -phi_e*phi_t*CXp0] , 
        [0 , 0 , 0]])
    dAeroMat['Cat']=0.125*c*rho*W*np.array([
            [6*CXp0+6*CY0-12*CX0*lambda_          , 0                                           , -6*lambda_*z*CY0-6*lambda_*z*CXp0-12*z*CX0                                                                     , 6*Ls*CY0+(6*Ls-3*c*lambda_)*CXp0-12*Ls*CX0*lambda_                                                             , 0                                            , 0]                                             , 
            [0                                    , 24*CY0+12*CYp0*lambda_-12*CX0*lambda_       , 0                                                                                                              , 0                                                                                                              , -24*lambda_*z*CY0+12*z*CYp0-12*z*CX0         , -24*lambda_*z*CY0+12*z*CYp0-12*z*CX0 ]         , 
            [12*lambda_*z*CY0-6*z*CYp0+6*z*CX0    , 0                                           , (6*Ls**2+12*z**2)*CY0+6*z**2*CYp0*lambda_+(6*Ls**2-3*Ls*c*lambda_)*CXp0+(-6*z**2*lambda_-12*Ls**2*lambda_)*CX0 , 18*Ls*lambda_*z*CY0+3*(-2*Ls+c*lambda_)*z*CYp0+6*Ls*lambda_*z*CXp0+18*z*Ls*CX0                                 , 0                                            , 0                                    ]         , 
            [6*Ls*CXp0+6*Ls*CY0-12*Ls*CX0*lambda_ , 0                                           , -18*Ls*lambda_*z*CY0-3*(-2*Ls+c*lambda_)*z*CYp0-6*Ls*lambda_*z*CXp0-18*z*Ls*CX0                                , (6*Ls**2+12*z**2)*CY0+6*z**2*CYp0*lambda_+(6*Ls**2-3*Ls*c*lambda_)*CXp0+(-6*z**2*lambda_-12*Ls**2*lambda_)*CX0 , 0                                            , 0                                    ]         , 
            [0                                    , 12*lambda_*z*CY0+12*lambda_*z*CXp0+24*z*CX0 , 0                                                                                                              , 0                                                                                                              , 12*z**2*CY0+12*z**2*CXp0-24*z**2*CX0*lambda_ , 12*z**2*CY0+12*z**2*CXp0-24*z**2*CX0*lambda_ ] , 
            [0                                    , 12*lambda_*z*CY0+12*lambda_*z*CXp0+24*z*CX0 , 0                                                                                                              , 0                                                                                                              , 12*z**2*CY0+12*z**2*CXp0-24*z**2*CX0*lambda_ , 12*z**2*CY0+12*z**2*CXp0-24*z**2*CX0*lambda_]])
    dAeroMat['Kat']    = 0.25*c*rho*W*W*np.array([
        [0 , 0 , 0                                    , 3*CY0+6*CY0*lambda_**2+6*CX0*lambda_-3*CXp0     , 0 , 0]  , 
        [0 , 0 , 0                                    , 0                                               , 0 , 0]  , 
        [0 , 0 , -3*Ls*CY0+6*Ls*CX0*lambda_-3*Ls*CXp0 , (6*CX0*lambda_**2-6*CY0*lambda_+3*CYp0+3*CX0)*z , 0 , 0]  , 
        [0 , 0 , (-3*CYp0+6*CY0*lambda_+3*CX0)*z      , -3*Ls*CY0+6*Ls*CX0*lambda_-3*Ls*CXp0            , 0 , 0]  , 
        [0 , 0 , 0                                    , 0                                               , 0 , 0]  , 
        [0 , 0 , 0                                    , 0                                               , 0 , 0]])
    dAeroMat['Ka1_bt'] = 0.5*c*rho*W*W*np.array([
        [0 , 0 , 0 , -phi_f*CYp0+2*phi_f*CY0*lambda_+phi_f*CX0 , 0 , 0]   , 
        [0 , 0 , 0 , -phi_e*CY0+2*phi_e*CX0*lambda_-phi_e*CXp0 , 0 , 0]   , 
        [0 , 0 , 0 , 0                                         , 0 , 0]])
    dAeroMat['Kb1_bt'] = 0.5*c*rho*W*W*np.array([
        [0 , 0 , -phi_f*CYp0+2*phi_f*CY0*lambda_+phi_f*CX0 , 0 , 0 , 0]  , 
        [0 , 0 , -phi_e*CY0+2*phi_e*CX0*lambda_-phi_e*CXp0 , 0 , 0 , 0]  , 
        [0 , 0 , 0                                         , 0 , 0 , 0]])
    dAeroMat['Ka0_tb'] = 0.5*c*rho*U0**2*np.array([
        [0 , 0 , 0]             , 
        [0 , 0 , -phi_t*CYp0]   , 
        [0 , 0 , 0]             , 
        [0 , 0 , 0]             , 
        [0 , 0 , -z*phi_t*CXp0] , 
        [0 , 0 , -z*phi_t*CXp0]])
    dAeroMat['Ka1_tb'] = 0.5*c*rho*U0**2*np.array([
        [0         , 0          , -phi_t*CXp0]    , 
        [0         , 0          , 0]              , 
        [0         , 0          , phi_t*z*CYp0]   , 
        [phi_f*CX0 , -phi_e*CY0 , -phi_t*Ls*CXp0] , 
        [0         , 0          , 0]              , 
        [0         , 0          , 0]])
    dAeroMat['Kb1_tb'] = 0.5*c*rho*U0**2*np.array([
        [0,0,0],
        [0,0,0],
        [phi_f*CX0,-phi_e*CY0,-phi_t*Ls*CXp0],
        [0,0,-phi_t*z*CYp0],
        [0,0,0],
        [0,0,0]])
    dAeroMat['Ca0_tb'] = 0.5*c*rho*W*np.array([
        [0                                                , 0                                        , 0]                             , 
        [2*phi_f*CY0+phi_f*CYp0*lambda_-phi_f*CX0*lambda_ , phi_e*CYp0-phi_e*CX0-2*phi_e*CY0*lambda_ , -1/2*c*phi_t*CYp0*lambda_]     , 
        [0                                                , 0                                        , 0]                             , 
        [0                                                , 0                                        , 0]                             , 
        [phi_f*(2*CX0+CXp0*lambda_+CY0*lambda_)*z         , -phi_e*(-CXp0-CY0+2*CX0*lambda_)*z       , -1/2*c*phi_t*CXp0*lambda_*z]   , 
        [phi_f*(2*CX0+CXp0*lambda_+CY0*lambda_)*z         , -phi_e*(-CXp0-CY0+2*CX0*lambda_)*z       , -1/2*c*phi_t*CXp0*lambda_*z]])
    dAeroMat['Ca1_tb'] = 0.5*c*rho*W*np.array([
        [2*phi_f*CX0+phi_f*CXp0*lambda_+phi_f*CY0*lambda_          , phi_e*CXp0+phi_e*CY0-2*phi_e*CX0*lambda_          , -1/2*c*phi_t*CXp0*lambda_]    , 
        [0                                                         , 0                                                 , 0]                            , 
        [-phi_f*(2*CY0+CYp0*lambda_-CX0*lambda_)*z                 , phi_e*(-CYp0+CX0+2*CY0*lambda_)*z                 , 1/2*c*phi_t*CYp0*lambda_*z]   , 
        [2*phi_f*Ls*CX0+phi_f*Ls*CXp0*lambda_+phi_f*Ls*CY0*lambda_ , phi_e*Ls*CXp0+phi_e*Ls*CY0-2*phi_e*Ls*CX0*lambda_ , -1/2*c*phi_t*Ls*CXp0*lambda_] , 
        [ 0                                                        , 0                                                 , 0]                            , 
        [0                                                         , 0                                                 , 0]])
    dAeroMat['Cb1_tb'] = 0.5*c*rho*W*np.array([
        [0                                                         , 0                                                 , 0]                            , 
        [0                                                         , 0                                                 , 0]                            , 
        [2*phi_f*Ls*CX0+phi_f*Ls*CXp0*lambda_+phi_f*Ls*CY0*lambda_ , phi_e*Ls*CXp0+phi_e*Ls*CY0-2*phi_e*Ls*CX0*lambda_ , -1/2*c*phi_t*Ls*CXp0*lambda_] , 
        [phi_f*(2*CY0+CYp0*lambda_-CX0*lambda_)*z                  , -phi_e*(-CYp0+CX0+2*CY0*lambda_)*z                , -1/2*c*phi_t*CYp0*lambda_*z]  , 
        [0                                                         , 0                                                 , 0]                            , 
        [0                                                         , 0                                                 , 0]])
    dAeroMat['Ca0_bt'] = 0.5*c*rho*W*np.array([
        [0 , 2*phi_f*CY0+phi_f*CYp0*lambda_-phi_f*CX0*lambda_ , 0 , 0 , -phi_f*(-CYp0+CX0+2*CY0*lambda_)*z , -phi_f*(-CYp0+CX0+2*CY0*lambda_)*z] , 
        [0 , 2*phi_e*CX0+phi_e*CXp0*lambda_+phi_e*CY0*lambda_ , 0 , 0 , -phi_e*(-CXp0-CY0+2*CX0*lambda_)*z , -phi_e*(-CXp0-CY0+2*CX0*lambda_)*z] ,
        [0 , 0 , 0 , 0 , 0 , 0]])
    dAeroMat['Ca1_bt'] = 0.5*c*rho*W*np.array([
        [phi_f*CYp0-phi_f*CX0-2*phi_f*CY0*lambda_ , 0 , -phi_f*(2*CY0+CYp0*lambda_-CX0*lambda_)*z , -phi_f*Ls*CX0+phi_f*Ls*CYp0-1/2*phi_f*c*CYp0*lambda_-2*phi_f*Ls*CY0*lambda_ , 0 , 0]  , 
        [phi_e*CXp0+phi_e*CY0-2*phi_e*CX0*lambda_ , 0 , -phi_e*(2*CX0+CXp0*lambda_+CY0*lambda_)*z , phi_e*Ls*CXp0+phi_e*Ls*CY0-1/2*phi_e*c*CXp0*lambda_-2*phi_e*Ls*CX0*lambda_  , 0 , 0]  , 
        [0                                        , 0 , 0                                         , 0                                                                           , 0 , 0]])
    dAeroMat['Cb1_bt'] = 0.5*c*rho*W*np.array([
        [0 , 0 , -phi_f*Ls*CX0+phi_f*Ls*CYp0-1/2*phi_f*c*CYp0*lambda_-2*phi_f*Ls*CY0*lambda_ , phi_f*(2*CY0+CYp0*lambda_-CX0*lambda_)*z , 0 , 0]  , 
        [0 , 0 , phi_e*Ls*CXp0+phi_e*Ls*CY0-1/2*phi_e*c*CXp0*lambda_-2*phi_e*Ls*CX0*lambda_  , phi_e*(2*CX0+CXp0*lambda_+CY0*lambda_)*z , 0 , 0]  , 
        [0 , 0 , 0                                                                           , 0                                        , 0 , 0]])

    return dAeroMat

def liftDrag(alpha): 
    pi=np.pi
    # Aerodynamic paramters
    CDfric       = 0.005
    alpha_stall  = 14 * pi / 180
    Dalpha_stall = 3 * pi / 180
    # Normal force coefficient and its derivative
    Cn = 2.25 * 2 * pi * np.sin(alpha) / (4 + pi * np.abs(np.sin(alpha)))
    Cnp = 4.5 * pi * np.cos(alpha) / (4 + pi * np.abs(np.sin(alpha))) - np.multiply(4.5 * pi ** 2 * np.sin(alpha) / (4 + pi * np.abs(np.sin(alpha))) ** 2.0 * np.sign(np.sin(alpha)),np.cos(alpha))
    # Separation function and its derivative
    f = 0.5 + 0.5 * np.tanh((alpha_stall - np.abs(alpha)) / Dalpha_stall)
    fp = np.multiply(- 1 / 2 * (1 - np.tanh((alpha_stall - np.abs(alpha)) / Dalpha_stall) ** 2),np.sign(alpha)) / Dalpha_stall
    # Aerodynamic coefficients
    Cl = np.multiply(2 * pi * alpha,f) + np.multiply(np.multiply(Cn,np.cos(alpha)),(1 - f))
    Cd = CDfric * f + np.multiply(np.multiply(3 / 4 * Cn,np.sin(alpha)),(1 - f))
    # Derivatives
    Clp = 2 * pi * f + np.multiply(2 * pi * alpha,fp) + np.multiply(np.multiply(Cnp,np.cos(alpha)),(1 - f)) - np.multiply(np.multiply(Cn,np.sin(alpha)),(1 - f)) - np.multiply(np.multiply(Cn,np.cos(alpha)),fp)
    Cdp = CDfric * fp + np.multiply(np.multiply(3 / 4 * Cnp,np.sin(alpha)),(1 - f)) + np.multiply(np.multiply(3 / 4 * Cn,np.cos(alpha)),(1 - f)) - np.multiply(np.multiply(3 / 4 * Cn,np.sin(alpha)),fp)
    return Cl,Cd,Clp,Cdp



def parameters(): 
    p=dict()
    # Environment
    p['rho'] = 1.225 #ρ Air density 1.225 kg/m 3
    # Rotor pameters
    p['R']   = 50            # R Rotor radius / blade length 50 m
    p['c']   = 3             # c Blade chord length 3 m
    p['acg'] = 1.2           # a cg Distance from torsional point aft to center of gravity 1.2 m
    p['EIx'] = 9879750000.0  # EI x Flapwise bending stiffness along whole blade length 9.87975 GNm 2
    p['EIy'] = 17564040000.0 # EI y Edgewise bending stiffness along whole blade length 17.56404 GNm 2
    p['GK']  = 176480400.0   # GK Torsional stiffness along whole blade length 0.1764804 GNm 2
    p['m']   = 220           # m Blade mass per unit-length 220 kg/m
    p['J']   = 275.75        # J Cross-sectional moment of inertia about center of gravity 275.75 kgm 2
    # Nacelle pameters
    p['Ls'] = 5           # Distance from tower top/drivetrain intersection to rotor center 5 m
    p['Gs'] = 500000000.0 # Torsional stiffness of drivetrain 
    p['Mn'] = 205000.0    # Nacelle mass 
    p['Ix'] = 4500000.0   # Tilt moment of inertia of nacelle
    p['Iy'] = 1200000.0   # Roll moment of inertia of nacelle
    p['Iz'] = 4500000.0   # Yaw moment of inertia of nacelle
    # Damping in logarithmic decrement
    p['etaf'] = 0.01  # Logarithmic decrement of first flapwise blade bending mode
    p['etae'] = 0.01  # Logarithmic decrement of first edgewise blade bending mode 
    p['etab'] = 0.04  # Logarithmic decrement of first blade torsional mode 
    p['etat'] = 0.005 # Logarithmic decrement of first tower bending modes 
    p['etad'] = 0.3   # Logarithmic decrement of first drivetrain torsional mode 
    # Nacelle/tower top pameters
    D    = 5.0            # Outer diameter of tower 
    d    = 4.92           # Inner diameter of tower 
    Et   = 211000000000.0 # t Young’s modulus of tower steel
    nu   = 0.33           # Possion’s ratio of tower steel 
    rhot = 7850           # Density of tower steel 
    H    = 70             # Tower height
    # Tower stiffness pameters derived from elementary loading cases
    # [Dubbel, Handbook of Mech. Eng.]
    EIt = Et * np.pi * (D ** 4 - d ** 4) / 64
    GKt = Et / (2 * (1 + nu)) * np.pi * (D ** 4 - d ** 4) / 32
    p['Kx']  = 12 * EIt / H ** 3  # Lateral stiffness of nacelle support
    p['Ky']  = p['Kx']            # Longitudinal stiffness of nacelle support
    p['Gx']  = 4 * EIt / H        # Tilt stiffness of nacelle support
    p['Gy']  = p['Gx']            # Roll stiffness of nacelle support
    p['Gz']  = GKt / H            # Yaw stiffness of nacelle support
    p['gxy'] = - 6 * EIt / H ** 2 # Coupling stiffness of nacelle support 
    # Nacelle/tower equivalent mass
    mt  = np.pi * (D ** 2 - d ** 2) / 4 * rhot
    h   = np.linspace(0,H,500)
    phi = bendingMode(h,1)
    p['Mtot'] = p['Mn'] + np.trapz(phi**2 * mt, h) / phi[-1] ** 2 #M Equivalent nacelle and tower mass 290.6 ton

    return p







