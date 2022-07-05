import numpy as np
from numpy.linalg import inv
# --------------------------------------------------------------------------------}
# --- Functions for state space model integrations
# --------------------------------------------------------------------------------{
def StateMatrix(Minv=None,C=None,K=None,M=None):
    if M is not None:
        Minv=inv(M)
    if C is None:
        C=np.zeros(K.shape)
    n = len(Minv)
    A = np.zeros((2*n,2*n))
    A[:n,n:] =  np.identity(n)
    A[n:,:n] = -np.dot(Minv,K)
    A[n:,n:] = -np.dot(Minv,C)
    return A

def vec_interp(t,vTime,vF):
    """ Interpolate vector known at discrete values (vTime, vF) to a given time `t` 
    TODO use and interpolant if possible
      t : scalar!
      vTime : 1d-array of length nt 
      vF : nDOF x nt array
    """
    F    = np.zeros(vF.shape[0])
    for iDOF,F_DOF in enumerate(vF):
        F[iDOF] = np.interp(t,vTime,F_DOF)
    return F

def B_interp(t,Minv,vTime,vF, flat=False):
    """ Interpolate B-vector from loads known at discrete values (vTime, vF) at a given time `t` """
    nDOF=len(vF)
    if flat:
        B = np.zeros(2*nDOF)
        F = vec_interp(t,vTime,vF)
        B[nDOF:] = np.dot(Minv,F)
    else:
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
