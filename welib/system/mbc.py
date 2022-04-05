import numpy as np


def MBC3_Bmat(nq, nf, psi1=0, Omega=None, symb=False): 
    """ 
    Multiblade transformation matrix for three blades
     See Hansen 2003 Eq. 14 15, 16

     q = B z

     Assumes that DOFs are ordered as:
        q: [q_b1,    q_b2   , q_b3   , q_nonrot]
        z: [a0_b1-3, a1_b1-3, b1_b1-3, q_nonrot]

    # TODO return inverse as well
    """
    if symb:
        from sympy import pi, cos, sin, eye, zeros, Rational
    else:
        from numpy import pi, cos, sin, eye, zeros

    psi2 = psi1 + (2 * pi / 3)
    psi3 = psi2 + (2 * pi / 3)
    if not symb:
        B = np.block([
            [eye(nq)  , eye(nq) * cos(psi1) , eye(nq) * sin(psi1) , zeros((nq, nf))] , 
            [eye(nq)  , eye(nq) * cos(psi2) , eye(nq) * sin(psi2) , zeros((nq, nf))] , 
            [eye(nq)  , eye(nq) * cos(psi3) , eye(nq) * sin(psi3) , zeros((nq, nf))] , 
            [zeros((nf , 3*nq)) , eye(nf)]]) # TODO Omega t for shaft DOFs?

        mu = np.block([
            [(1 / 3) * eye(nq) , zeros((nq , nq)), zeros((nq, nq)), zeros((nq, nf))], 
            [zeros((nq, nq)), (2 / 3) * eye(nq)  , zeros((nq, nq)), zeros((nq, nf))], 
            [zeros((nq, nq)), zeros((nq, nq)) , (2 / 3) * eye(nq) , zeros((nq, nf))],  # NOTE: should be 2/3
            [zeros((nf, nq)), zeros((nf, nq)) , zeros((nf, nq)), eye(nf)]])

        if Omega is None:
            Bdot = None
        else:
            Bdot = np.block([
                [zeros((nq, nq)) , Omega * eye(nq) * - sin(psi1) , Omega * eye(nq) * cos(psi1) , zeros((nq, nf))] , 
                [zeros((nq, nq)) , Omega * eye(nq) * - sin(psi2) , Omega * eye(nq) * cos(psi2) , zeros((nq, nf))] , 
                [zeros((nq, nq)) , Omega * eye(nq) * - sin(psi3) , Omega * eye(nq) * cos(psi3) , zeros((nq, nf))] , 
                [zeros((nf , 3*nq)), zeros((nf, nf))]])

        Binv = mu.dot(B.T)
        R = Binv.dot(Bdot) # TODO analytical

    else:
        Iq = eye(nq)
        If = eye(nf)
        Z0 = zeros(nq,nq)
        Z1 = zeros(nq,nf)
        Z2 = zeros(nf,3*nq)
        # B matrix
        B0= Iq
        B0= B0.row_join(Iq*cos(psi1))
        B0= B0.row_join(Iq*sin(psi1))
        B0= B0.row_join(Z1)
        B1= Iq
        B1= B1.row_join(Iq*cos(psi2))
        B1= B1.row_join(Iq*sin(psi2))
        B1= B1.row_join(Z1)
        B2= Iq
        B2= B2.row_join(Iq*cos(psi3))
        B2= B2.row_join(Iq*sin(psi3))
        B2= B2.row_join(Z1)
        B3= Z2
        B3= B3.row_join(If)
        B=B0
        B=B.col_join(B1)
        B=B.col_join(B2)
        B=B.col_join(B3)

        # Mu
        mu0= Rational(1,3)*Iq
        mu0= mu0.row_join(Z0)
        mu0= mu0.row_join(Z0)
        mu0= mu0.row_join(Z1)
        mu1= Z0
        mu1= mu1.row_join(Rational(2,3)*Iq)
        mu1= mu1.row_join(Z0)
        mu1= mu1.row_join(Z1)
        mu2= Z0
        mu2= mu2.row_join(Z0)
        mu2= mu2.row_join(Rational(2,3)*Iq)
        mu2= mu2.row_join(Z1)
        mu3= Z2
        mu3= mu3.row_join(If)
        mu=mu0
        mu=mu.col_join(mu1)
        mu=mu.col_join(mu2)
        mu=mu.col_join(mu3)

        Bdot = None
        Binv = mu * (B.T)
        R = None
        #R = Binv.dot(Bdot) # TODO analytical


    return B, Binv, Bdot, mu, R

def MBC3_Rot2Fixed_TS(Omega, time, q, nqb, nb, nf):
    """
    Transfer coordinate transformation from rotating DOF to MBC/Fixed DOFs
    for a time series
    Assumes that DOFs are ordered as:
        q: [q_B1,   ...   , q_Bnb   , q_f]
        q_B1 = [q1,1, ... q1,nqb]

    INPUTS:
        Omega: rotational speed (rad/s)
        time : array (nt) time vector
        q : array (nt, n) degrees of freedom as function of time, arranged according to:
                q: [q_B1,   ...   , q_Bnb   , q_f]
        nqb: number of DOFs per blade
        nb: number of blades
        nf: number of non rotating DOFs
    """
    nt, n = q.shape

    n2 = nqb*nb+nf
    if n!=n2:
        raise Exception('Input vector `q` should have {} columns instead of {}'.format(n,n2))
    if nt!=len(time):
        raise Exception('Input vector `time` should have {} lines instead of {}'.format(len(time),nt))

    z   = np.zeros((len(time),n))

    for it,t in enumerate(time):
        psi1 = Omega*t 
        B = MBC3_Bmat(nq=nqb, nf=nf, psi1=psi1, Omega=Omega)[0]
        Binv= np.linalg.inv(B)  # TODO TODO TODO analytical B
        z[it,:] = Binv.dot(q[it,:])
    return z

def MBC3_Fixed2Rot_TS(Omega, time, z, nqb, nb, nf):
    """ """
    nt, n = z.shape
    n2 = nqb*nb+nf
    if n!=n2:
        raise Exception('Input vector `z` should have {} columns instead of {}'.format(n,n2))
    if nt!=len(time):
        raise Exception('Input vector `time` should have {} lines instead of {}'.format(len(time),nt))

    q   = np.zeros((len(time),n))
    for it,t in enumerate(time):
        psi1 = Omega*t 
        B = MBC3_Bmat(nq=nqb, nf=nf, psi1=psi1, Omega=Omega)[0]
        q[it,:] = B.dot(z[it,:])
    return q
