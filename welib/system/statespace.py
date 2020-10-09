
# --------------------------------------------------------------------------------}
# --- Functions for state space model integrations
# --------------------------------------------------------------------------------{
def StateMatrix(M,C,K):
    A = np.zeros((2*len(M),2*len(M)))
    A[:len(M),len(M):]  =  np.identity(len(M))
    A[len(M):,:len(M)]  = -np.dot(inv(M),K)
    A[len(M):,len(M):]  = -np.dot(inv(M),C)
    return A


def vec_interp(t,vTime,vF):
    """ Interpolate vector known at discrete values (vTime, vF) to a given time `t` """
    F    = np.zeros(vF.shape[0])
    for iDOF,F_DOF in enumerate(vF):
        F[iDOF] = np.interp(t,vTime,F_DOF)
    return F

def B_interp(t,M,vTime,vF):
    """ Interpolate B-vector from loads known at discrete values (vTime, vF) at a given time `t` """
    nDOF=len(vF)
    B = np.zeros((2*nDOF,1))
    F = vec_interp(t,vTime,vF)
    B[nDOF:,0] = np.dot(inv(M),F)
    return B

def B_reg(t,M,F):
    """ Return B vector from loads at time t and mass matrix """
    nDOF=len(F)
    B = np.zeros((2*nDOF,1))
    B[nDOF:,0] = np.dot(inv(M),F)
    return B

def dxdt(q, t, A, M, vTime, vF): 
    B = B_interp(t, M, vTime, vF)
    dq = np.dot(A,q)+B
    return dq
 
def odefun(t, dq):
    return dxdt(dq, t, A, vTime, vF)
