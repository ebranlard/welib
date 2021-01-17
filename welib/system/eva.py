""" 
Eigenvalue analyses tools for mechnical system: 
   mass matrix M, stiffness matrix K and possibly damping matrix C
"""
import pandas as pd    

import numpy as np
from scipy import linalg

def polyeig(*A):
    """
    Solve the polynomial eigenvalue problem:
        (A0 + e A1 +...+  e**p Ap)x = 0

    Return the eigenvectors [x_i] and eigenvalues [e_i] that are solutions.

    Usage:
        X,e = polyeig(A0,A1,..,Ap)

    Most common usage, to solve a second order system: (K + C e + M e**2) x =0
        X,e = polyeig(K,C,M)

    """
    if len(A)<=0:
        raise Exception('Provide at least one matrix')
    for Ai in A:
        if Ai.shape[0] != Ai.shape[1]:
            raise Exception('Matrices must be square')
        if Ai.shape != A[0].shape:
            raise Exception('All matrices must have the same shapes');

    n = A[0].shape[0]
    l = len(A)-1 
    # Assemble matrices for generalized problem
    C = np.block([
        [np.zeros((n*(l-1),n)), np.eye(n*(l-1))],
        [-np.column_stack( A[0:-1])]
        ])
    D = np.block([
        [np.eye(n*(l-1)), np.zeros((n*(l-1), n))],
        [np.zeros((n, n*(l-1))), A[-1]          ]
        ]);
    # Solve generalized eigenvalue problem
    e, X = linalg.eig(C, D);
    if np.all(np.isreal(e)):
        e=np.real(e)
    X=X[:n,:]

    # Sort eigen values
    #I = np.argsort(e)
    #X = X[:,I]
    #e = e[I]

    # Scaling each mode by max
    X /= np.tile(np.max(np.abs(X),axis=0), (n,1))

    return X, e


def eig(K,M=None, freq_out=False, sort=True):
    """ performs eigenvalue analysis and return same values as matlab 

    returns:
       Q     : matrix of column eigenvectors
       Lambda: matrix where diagonal values are eigenvalues
              frequency = np.sqrt(np.diag(Lambda))/(2*np.pi)
         or
    frequencies (if freq_out is True)
    """
    if M is not None:
        D,Q = linalg.eig(K,M)
    else:
        D,Q = linalg.eig(K)
    # --- rescaling TODO, this can be made smarter
    if M is not None:
        for j in range(M.shape[1]):
            q_j = Q[:,j]
            modalmass_j = np.dot(q_j.T,M).dot(q_j)
            Q[:,j]= Q[:,j]/np.sqrt(modalmass_j)
        Lambda=np.dot(Q.T,K).dot(Q)
        lambdaDiag=np.diag(Lambda) # Note lambda might have off diganoal values due to numerics
        I = np.argsort(lambdaDiag)
        # Sorting eigen values
        if sort:
            Q          = Q[:,I]
            lambdaDiag = lambdaDiag[I]
        if freq_out:
            Lambda = np.sqrt(lambdaDiag)/(2*np.pi) # frequencies [Hz]
        else:
            Lambda = np.diag(lambdaDiag) # enforcing purely diagonal
    else:
        Lambda = np.diag(D)

    return Q,Lambda


def eigA(A, nq=None, nq1=None, full=False):
    """
    Perform eigenvalue analysis on a "state" matrix A
    where states are assumed to be orderd as {q, q_dot, q1}

    if full is True, the entire eigenvectors are returned, otherwise, 
    only the part associated with q and q1 are returned
    """
    n,m = A.shape

    if m!=n:
        raise Exception('Matrix needs to be squared')
    if nq is None:
        if nq1 is None:
            nq1=0
        nq = int((n-nq1)/2)
    else:
        nq1 = n-2*nq
    if n!=2*nq+nq1 or nq1<0:
        raise Exception('Number of 1st and second order dofs should match the matrix shape (n= 2*nq + nq1')
    Q, Lambda = eig(A, sort=False)
    v = np.diag(Lambda)

    if not full:
        Q=np.delete(Q, slice(nq,2*nq), axis=0)

    # Selecting eigenvalues with positive imaginary part (frequency)
    Ipos = np.imag(v)>0
    Q = Q[:,Ipos]
    v = v[Ipos]

    # Frequencies and damping based on compled eigenvalues
    omega_0 = np.abs(v)              # natural cylic frequency [rad/s]
    freq_d  = np.imag(v)/(2*np.pi)   # damped frequency [Hz]
    zeta    = - np.real(v)/omega_0   # damping ratio
    freq_0  = omega_0/(2*np.pi)      # natural frequency [Hz]

    # Sorting
    I = np.argsort(freq_d)
    freq_d = freq_d[I]
    freq_0 = freq_0[I]
    zeta   = zeta[I]
    Q      = Q[:,I]
    return freq_d, zeta, Q, freq_0 


def eigMCK(M, C, K, method='diag_beta'): 
    """ """
    if method.lower()=='diag_beta':
        ## using K, M and damping assuming diagonal beta matrix (Rayleigh Damping)
        Q, Lambda   = eig(K,M, sort=False) # provide scaled EV, important, no sorting here!
        freq        = np.sqrt(np.diag(Lambda))/(2*np.pi)
        betaMat     = np.dot(Q,C).dot(Q.T)
        xi          = (np.diag(betaMat)*np.pi/(2*np.pi*freq))
        xi[xi>2*np.pi] = np.NAN
        zeta        = xi/(2*np.pi)
        freq_d      = freq*np.sqrt(1-zeta**2)
        # Sorting here
        I = np.argsort(freq_d)
        freq   = freq[I]
        freq_d = freq_d[I]
        zeta   = zeta[I]
        xi     = xi[I]
        Q      = Q[:,I]

    #    return Q, Lambda,freq, betaMat,xi,zeta
    elif method.lower()=='full_matrix':
        ## Method 2 - Damping based on K, M and full D matrix
        Q,e = polyeig(K,C,M)
        #omega0 = np.abs(e)
        zeta = - np.real(e) / np.abs(e)
        freq_d = np.imag(e) / (2*np.pi)
        # Sorting
        I = np.argsort(freq_d)
        freq_d = freq_d[I]
        zeta   = zeta[I]
        Q      = Q[:,I]
        # Keeping only positive frequencies
        bValid = freq_d > 1e-08
        freq_d = freq_d[bValid]
        zeta   = zeta[bValid]
        Q      = Q[:,bValid]
        # Undamped frequency and pseudo log dec
        freq = freq_d / np.sqrt(1 - zeta**2)
        xi = 2 * np.pi * zeta
        # logdec2 = 2*pi*dampratio_sorted./sqrt(1-dampratio_sorted.^2);
    else:
        raise NotImplementedError()
    return freq_d,zeta,Q,freq,xi


if __name__=='__main__':
    np.set_printoptions(linewidth=300, precision=4)
    nDOF   = 2
    M      = np.zeros((nDOF,nDOF))
    K      = np.zeros((nDOF,nDOF))
    C      = np.zeros((nDOF,nDOF))
    M[0,0] = 430000;
    M[1,1] = 42000000;
    C[0,0] = 7255;
    C[1,1] = M[1,1]*0.001;
    K[0,0] = 2700000.;
    K[1,1] = 200000000.;

    freq_d, zeta, Q, freq, xi = eigMCK(M,C,K)
    print(freq_d)
    print(Q)


    #M = diag([3,0,0,0], [0, 1,0,0], [0,0,3,0],[0,0,0, 1])
    M = np.diag([3,1,3,1])
    C = np.array([[0.4 , 0 , -0.3 , 0] , [0 , 0  , 0 , 0] , [-0.3 , 0 , 0.5 , -0.2 ] , [ 0 , 0 , -0.2 , 0.2]])
    K = np.array([[-7  , 2 , 4    , 0] , [2 , -4 , 2 , 0] , [4    , 2 , -9  , 3    ] , [ 0 , 0 , 3    , -3]])

    X,e = polyeig(K,C,M)
    print('X:\n',X)
    print('e:\n',e)
    # Test taht first eigenvector and valur satisfy eigenvalue problem:
    s = e[0];
    x = X[:,0];
    res = (M*s**2 + C*s + K).dot(x) # residuals
    assert(np.all(np.abs(res)<1e-12))

