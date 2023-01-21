import numpy as np




def MBC3_Bmat(nqb, nf, psi1=0, psi2=None, psi3=None, psi0=None, Omega=None, ns=0, ordering='increasing', symb=False): 
    """ 
    Multiblade transformation matrix for three blades
         q = B z

    Assumes that DOFs are ordered as:
       In   q: [q_b1,    q_b2   , q_b3   , q_nonrot] 
       Out  z: [a0_b1-3, a1_b1-3, b1_b1-3, q_nonrot]

        q_b* is of nqb
        q_b* is of nqb

     See Hansen 2003 Eq. 14 15, 16

     INPUTS:
     - nqb: number of dof per blade (length of q_b1)
     - nf: number of dof that are already in the fixed frame (length of q_nonrot)
     - ns: number of shaft degrees of freedom (0, or 2)
     - psi1, psi2, psi3: azimuthal position of each blade. If ordering is provided, then the 
              azimuths of psi2 and psi3 are determined from psi1
     - psi0: azimuth of the rotor, might be different from psi1 if there is a given blade offset azimuth.
             If None, psi0=psi1
     - Omega : rotor speed used for Bdot, Bddot and R. If None, a rotor speed of "1" is used
     - ordering: if 'increasing', then psi2 = psi1+2pi/3, psi3 = psi2+2pi/3
                 if 'decreasing', then psi2 = psi1-2pi/3, psi3 = psi2-2pi/3
     - symb: if true, use sympy, otherwise use numpy

    # TODO return inverse as well
    """
    if Omega is None:
        # NOTE: we return "unit" matrices for Bdot,  Bddot, and R
        Omega = 1
    if ns>0 and psi0 is None:
        psi0 = psi1
    if not ns in [0,2]:
        raise Exception('ns should be 0 or 2')
    if symb:
        import sympy as sp
        from sympy import pi, cos, sin, eye, Rational
        zeros = lambda tup: sp.zeros(tup[0],tup[1])
        block = lambda x: sp.Matrix(sp.BlockMatrix(x))
        array = sp.Matrix
        r13 = Rational(1,3)
        r23 = Rational(2,3)
    else:
        from numpy import pi, cos, sin, eye, zeros
        block = np.block
        array = np.array
        r13 = 1/3
        r23 = 2/3

    if psi2 is None and psi3 is None:
        if ordering=='increasing':
            psi2 = psi1 + (2 * pi / 3)
            psi3 = psi2 + (2 * pi / 3)
        else:
            psi2 = psi1 - (2 * pi / 3)
            psi3 = psi2 - (2 * pi / 3)
    if psi2 is None or psi3 is None:
        raise Exception('Provide both psi2 and psi3, or use ordering')

    # Definining useful matrices
    Iq = eye(nqb)
    If = eye(nf)
    Is = eye(ns)
    Zfq = zeros((nf,nqb))
    Zff = zeros((nf,nf))
    Zfs = zeros((nf,ns))
    Zss = zeros((ns,ns))
    Zsf = zeros((ns,nf))
    Zsq = zeros((ns,nqb))
    Zqs = zeros((nqb,ns))
    Zqf = zeros((nqb,nf))
    Zqq = zeros((nqb,nqb))

    if ns>0:
        RS = array([
            [cos(psi0),-sin(psi0)],
            [sin(psi0), cos(psi0)]])
        RSR= array([
            [0     , -Omega] , 
            [Omega , 0]])
        RSBd = Omega * array([
            [-sin(psi0),-cos(psi0)],
            [ cos(psi0),-sin(psi0)]])
        RSBdd = Omega**2 * array([
            [-cos(psi0), sin(psi0)],
            [-sin(psi0),-cos(psi0)]])
    else:
        RS  = Zss
        RSR = Zss
        RSBd = Zss
        RSBdd = Zss
    B = block([
        [Iq  , Iq * cos(psi1) , Iq * sin(psi1) , Zqs, Zqf]  , 
        [Iq  , Iq * cos(psi2) , Iq * sin(psi2) , Zqs, Zqf]  , 
        [Iq  , Iq * cos(psi3) , Iq * sin(psi3) , Zqs, Zqf]  , 
        [Zsq , Zsq            , Zsq            , RS, Zsf]  , 
        [Zfq , Zfq            , Zfq            , Zfs, If]])

    mu = block([
        [r13*Iq , Zqq    , Zqq    , Zqs , Zqf] , 
        [Zqq    , r23*Iq , Zqq    , Zqs , Zqf] , 
        [Zqq    , Zqq    , r23*Iq , Zqs , Zqf] , 
        [Zsq    , Zsq    , Zsq    , Is  , Zsf] , 
        [Zfq    , Zfq    , Zfq    , Zfs , If ]])

    # B^{-1} = mu B^T
    Binv = mu @ (B.T)

    # --- Time derivatives
    # Bdot = B * R
    Bdot = block([
        [Zqq , -Omega * Iq * sin(psi1) , Omega * Iq * cos(psi1), Zqs  , Zqf] , 
        [Zqq , -Omega * Iq * sin(psi2) , Omega * Iq * cos(psi2), Zqs  , Zqf] , 
        [Zqq , -Omega * Iq * sin(psi3) , Omega * Iq * cos(psi3), Zqs  , Zqf] , 
        [Zsq , Zsq                     , Zsq                   , RSBd, Zsf] ,
        [Zfq , Zfq                     , Zfq                   , Zfs  , Zff]])
    # R = Binv.dot(Bdot)
    R = block([
        [Zqq , Zqq       , Zqq      , Zqs  , Zqf ]  , 
        [Zqq , Zqq       , Omega*Iq , Zqs  , Zqf ]  , 
        [Zqq , -Omega*Iq , Zqq      , Zqs  , Zqf ]  , 
        [Zsq , Zsq       , Zsq      , RSR  , Zsf ]  , 
        [Zfq , Zfq       , Zfq      , Zfs  , Zff]])
    # Bddot = B * R^2 = Bdot * R
    Bddot = block([
        [Zqq , -Omega**2 * Iq * cos(psi1) , -Omega**2 * Iq * sin(psi1), Zqs , Zqf]   , 
        [Zqq , -Omega**2 * Iq * cos(psi2) , -Omega**2 * Iq * sin(psi2), Zqs , Zqf]   , 
        [Zqq , -Omega**2 * Iq * cos(psi3) , -Omega**2 * Iq * sin(psi3), Zqs , Zqf]   , 
        [Zsq ,    Zsq                     ,   Zsq                     , RSBdd, Zsf] , 
        [Zfq ,    Zfq                     ,   Zfq                     , Zfs , Zff]])

    return B, Binv, Bdot, Bddot, mu, R


def MBC3_MCK(MM, CC, KK, B, Binv, Bdot, Bddot, simplify=False):
    """ 
    Apply MBC3 transformation to a second order system with mass, damping and stiffness matrix.
    "tilde" matrices
    See Hansen 2003 Eq. 16

    INPUTS:
     - MM: mass matrix, containing blade, shaft and fixed states
     - CC: damping matrix
     - KK: stiffness matrix
     - B,*: see MBC_Bmat 
    """
    # Compute transformed matrices
    MB = Binv @ MM @ B
    DB = Binv @ CC @ B
    KB = Binv @ KK @ B
    MNR = MB
    DNR = 2 * Binv @ MM @ Bdot +  DB
    KNR = Binv @ MM @ Bddot  + Binv @ CC @ Bdot  +  KB
    if simplify: # only for sympy
        MNR.simplify()
        DNR.simplify()
        KNR.simplify()
    return MNR, DNR, KNR


def Coleman2Comp(a0, a1, b1):
    """ """
    A0    = np.sqrt( np.real(a0) ** 2 + np.imag(a0) ** 2)
    ABW   = 0.5 * np.sqrt((np.real(b1) - np.imag(a1)) ** 2 + (np.imag(b1) + np.real(a1)) ** 2)
    AFW   = 0.5 * np.sqrt((np.real(b1) + np.imag(a1)) ** 2 + (np.real(a1) - np.imag(b1)) ** 2)
    phi0  = np.arctan2( np.imag(a0),  np.real(a0))
    phiBW = np.arctan2(np.imag(a1) - np.real(b1),np.imag(b1) + np.real(a1))
    phiFW = np.arctan2(np.real(b1) + np.imag(a1),np.real(a1) - np.imag(b1))
    return A0,ABW,AFW,phi0,phiBW,phiFW


def MBC3_Rot2Fixed_TS(Omega, time, q, nqb, nf):
    """
    Perform coordinate transformation from rotating DOF to MBC/Fixed DOFs for a time series
    Assumes that DOFs are ordered as:
        q: [q_B1,   ...   , q_Bnb   , q_f]
        q_B1 = [q1,1, ... q1,nqb]

    INPUTS:
      - Omega: rotational speed (rad/s)
      - time : array (nt) time vector
      - q : array (nt, n) degrees of freedom as function of time, arranged according to:
                q: [q_B1,   ...   , q_Bnb   , q_f]
      - nqb: number of DOFs per blade
      - nf: number of non rotating DOFs
    """
    nt, n = q.shape

    n2 = nqb*3+nf
    if n!=n2:
        raise Exception('Input vector `q` should have {} columns instead of {}'.format(n,n2))
    if nt!=len(time):
        raise Exception('Input vector `time` should have {} lines instead of {}'.format(len(time),nt))

    z   = np.zeros((len(time),n))

    for it,t in enumerate(time):
        psi1 = Omega*t 
        B = MBC3_Bmat(nqb=nqb, nf=nf, psi1=psi1, Omega=Omega)[0]
        Binv= np.linalg.inv(B)  # TODO TODO TODO analytical B
        z[it,:] = Binv.dot(q[it,:])
    return z

def MBC3_Fixed2Rot_TS(Omega, time, z, nqb, nf):
    """ """
    nt, n = z.shape
    n2 = nqb*3+nf
    if n!=n2:
        raise Exception('Input vector `z` should have {} columns instead of {}'.format(n,n2))
    if nt!=len(time):
        raise Exception('Input vector `time` should have {} lines instead of {}'.format(len(time),nt))

    q   = np.zeros((len(time),n))
    for it,t in enumerate(time):
        psi1 = Omega*t 
        B = MBC3_Bmat(nqb=nqb, nf=nf, psi1=psi1, Omega=Omega)[0]
        q[it,:] = B.dot(z[it,:])
    return q
