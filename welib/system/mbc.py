import numpy as np


def MBC3_Bmat(ndofq, ndoff, psi1=0, Omega=None): 
    """ 
    Multiblade transformation matrix for three blades
     See Hansen 2003 Eq. 14 15, 16

     q = B z

     Assumes that DOFs are ordered as:
        q: [q_b1,    q_b2   , q_b3   , q_nonrot]
        z: [a0_b1-3, a1_b1-3, b1_b1-3, q_nonrot]

    # TODO return inverse as well
    """
    psi2 = psi1 + (2 * np.pi / 3)
    psi3 = psi2 + (2 * np.pi / 3)
    B = np.block([
        [np.eye(ndofq)  , np.eye(ndofq) * np.cos(psi1) , np.eye(ndofq) * np.sin(psi1) , np.zeros((ndofq, ndoff))] , 
        [np.eye(ndofq)  , np.eye(ndofq) * np.cos(psi2) , np.eye(ndofq) * np.sin(psi2) , np.zeros((ndofq, ndoff))] , 
        [np.eye(ndofq)  , np.eye(ndofq) * np.cos(psi3) , np.eye(ndofq) * np.sin(psi3) , np.zeros((ndofq, ndoff))] , 
        [np.zeros((ndoff , 3*ndofq)) , np.eye(ndoff)]]) # TODO Omega t for shaft DOFs?
#     Bdot = np.block([
#         [np.zeros((3 , 3)) , Omega * np.eye(3 , 3) * - np.sin(psi1) , Omega * np.eye(3 , 3) * np.cos(psi1) , np.zeros((3 , 6))] , 
#         [np.zeros((3 , 3)) , Omega * np.eye(3 , 3) * - np.sin(psi2) , Omega * np.eye(3 , 3) * np.cos(psi2) , np.zeros((3 , 6))] , 
#         [np.zeros((3 , 3)) , Omega * np.eye(3 , 3) * - np.sin(psi3) , Omega * np.eye(3 , 3) * np.cos(psi3) , np.zeros((3 , 6))] , 
#         [np.zeros((6 , 9)) , np.zeros((6      , 6))]])
#     mu = np.block([
#         [(1 / 3) * np.eye(3 , 3)  , np.zeros((3        , 3)) , np.zeros((3        , 3)) , np.zeros((3 , 3)) , np.zeros((3 , 3))]  , 
#         [np.zeros((3        , 3)) , (2 / 3) * np.eye(3 , 3)  , np.zeros((3        , 3)) , np.zeros((3 , 3)) , np.zeros((3 , 3))]  , 
#         [np.zeros((3        , 3)) , np.zeros((3        , 3)) , (2 / 3) * np.eye(3 , 3)  , np.zeros((3 , 3)) , np.zeros((3 , 3))]  , 
#         [np.zeros((3        , 3)) , np.zeros((3        , 3)) , np.zeros((3        , 3)) , np.eye(3    , 3)  , np.zeros((3 , 3))]  , 
#         [np.zeros((3        , 3)) , np.zeros((3        , 3)) , np.zeros((3        , 3)) , np.zeros((3 , 3)) , np.eye(3    , 3)]])
#     R = np.linalg.inv(B).dot(Bdot)

    return B

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
        q : array (nt, nDOF) degrees of freedom as function of time, arranged according to:
                q: [q_B1,   ...   , q_Bnb   , q_f]
        nqb: number of DOFs per blade
        nb: number of blades
        nf: number of non rotating DOFs
    """
    nt, nDOF = q.shape

    nDOF2 = nqb*nb+nf
    if nDOF!=nDOF2:
        raise Exception('Input vector `q` should have {} columns instead of {}'.format(nDOF,nDOF2))
    if nt!=len(time):
        raise Exception('Input vector `time` should have {} lines instead of {}'.format(len(time),nt))

    z   = np.zeros((len(time),nDOF))

    for it,t in enumerate(time):
        psi1 = Omega*t 
        B = MBC3_Bmat(ndofq=nqb, ndoff=nf, psi1=psi1, Omega=Omega)
        Binv= np.linalg.inv(B)  # TODO TODO TODO analytical B
        z[it,:] = Binv.dot(q[it,:])
    return z

def MBC3_Fixed2Rot_TS(Omega, time, z, nqb, nb, nf):
    """ """
    nt, nDOF = z.shape
    nDOF2 = nqb*nb+nf
    if nDOF!=nDOF2:
        raise Exception('Input vector `z` should have {} columns instead of {}'.format(nDOF,nDOF2))
    if nt!=len(time):
        raise Exception('Input vector `time` should have {} lines instead of {}'.format(len(time),nt))

    q   = np.zeros((len(time),nDOF))
    for it,t in enumerate(time):
        psi1 = Omega*t 
        B = MBC3_Bmat(ndofq=nqb, ndoff=nf, psi1=psi1, Omega=Omega)
        q[it,:] = B.dot(z[it,:])
    return q
