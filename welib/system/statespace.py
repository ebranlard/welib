import numpy as np
from numpy.linalg import inv
# --------------------------------------------------------------------------------}
# --- Functions for state space model integrations
# --------------------------------------------------------------------------------{
def StateMatrix(Minv,C,K):
    n = len(Minv)
    A = np.zeros((2*n,2*n))
    A[:n,n:] =  np.identity(n)
    A[n:,:n] = -np.dot(Minv,K)
    A[n:,n:] = -np.dot(Minv,C)
    return A


def vec_interp(t,vTime,vF):
    """ Interpolate vector known at discrete values (vTime, vF) to a given time `t` """
    F    = np.zeros(vF.shape[0])
    for iDOF,F_DOF in enumerate(vF):
        F[iDOF] = np.interp(t,vTime,F_DOF)
    return F

def B_interp(t,Minv,vTime,vF):
    """ Interpolate B-vector from loads known at discrete values (vTime, vF) at a given time `t` """
    nDOF=len(vF)
    B = np.zeros((2*nDOF,1))
    F = vec_interp(t,vTime,vF)
    B[nDOF:,0] = np.dot(Minv,F)
    return B

def B_reg(t,Minv,F):
    """ Return B vector from loads at time t and mass matrix """
    nDOF=len(F)
    B = np.zeros((2*nDOF,1))
    B[nDOF:,0] = np.dot(Minv,F)
    return B

def dxdt(q, t, A, M, vTime, vF): 
    B = B_interp(t, M, vTime, vF)
    dq = np.dot(A,q)+B
    return dq
 
def odefun(t, dq):
    return dxdt(dq, t, A, vTime, vF)
