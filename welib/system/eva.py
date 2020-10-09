""" 
Eigenvalue analyses tools for mechnical system: 
   mass matrix M, stiffness matrix K and possibly damping matrix C
"""
import numpy as np
pi=np.pi
from scipy import linalg
import pandas as pd    

def eig(K,M):
    """ performs eigenvalue analysis and return sam values as matlab """
    D,Q = linalg.eig(K,M)
    # --- rescaling TODO, this can be made smarter
    for j in range(M.shape[1]):
        q_j = Q[:,j]
        modalmass_j = np.dot(q_j.T,M).dot(q_j)
        Q[:,j]= Q[:,j]/np.sqrt(modalmass_j)
    Lambda=np.dot(Q.T,K).dot(Q)
    lambdaDiag=np.diag(Lambda) # Note lambda might have off diganoal values due to numerics
    I = np.argsort(lambdaDiag)
    # Sorting eigen values
    Q=Q[:,I]
    Lambda = np.diag(lambdaDiag[I]) # enforcing purely diagonal
    return Q,Lambda

def eigMCK(M,C,K, method='diag_beta'): 
    """ """
    if method.lower()=='diag_beta':
        ## using K, M and damping assuming diagonal beta matrix (Rayleigh Damping)
        Q, Lambda   = eig(K,M) # provide scaled EV
        freq        = np.sqrt(np.diag(Lambda))/(2*pi)
        betaMat     = np.dot(Q,C).dot(Q.T)
        xi          = (np.diag(betaMat)*pi/(2*pi*freq))
        xi[xi>2*pi] = np.NAN
        zeta        = xi/(2*pi)
        freq_d      = freq*np.sqrt(1-zeta**2)
    #    return Q, Lambda,freq, betaMat,xi,zeta
#         if method.lower()=='full_matrix':
#             ## Method 2 - Damping based on K, M and full D matrix
#             Q,e = polyeig(K,C,M)
#             zeta = - real(e) / np.abs(e)
#             freq_d = imag(e) / 2 / pi
#             # Sorting
#             freq_d,idx = __builtint__.sorted(freq_d)
#             zeta = zeta(idx)
#             # Keeping only positive frequencies
#             bValid = freq_d > 1e-08
#             freq_d = freq_d(bValid)
#             zeta = zeta(bValid)
#             # Undamped frequency and pseudo log dec
#             freq = freq_d / np.sqrt(1 - zeta ** 2)
#             xi = 2 * pi * zeta
            # logdec2 = 2*pi*dampratio_sorted./sqrt(1-dampratio_sorted.^2);
# valid=freq_sorted>0.1; # Treshold value
    return freq_d,zeta,Q,freq,xi


if __name__=='__main__':
    nDOF   = 2
    M      = np.zeros((nDOF,nDOF))
    K      = np.zeros((nDOF,nDOF))
    C      = np.zeros((nDOF,nDOF))
    M[0,0] = 435024.04730258
    M[1,1] = 42864056.19657615
    C[0,0] = 7255.30090655
    K[0,0] = 2751727.25652762
    z = eigMCK(M,C,K)
    print(z[0],z[3])
