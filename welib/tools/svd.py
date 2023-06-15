""" 
Set of tools for svd
"""
import numpy as np



def svd(A, econ=True):
    """ svd with similar interface to MATLAB"""
    U, S, VT = np.linalg.svd(X, full_matrices = not econ)
    S = np.diag(S)
    V = VT.T
    return U, S, V

def truncate(U, S, VT, r):
    """ Truncate SVD to order r """
    Ur  = U [:,:r]
    VTr = VT[:r,:]
    if S.ndim == 1:
        Sr  = S [:r]
    else:
        Sr  = S [:r,:r]
    return Ur, Sr, VTr

def evaluate(Ur, Sr, VTr):
    """ evaluate SVD (typically reduced order) """
    if Sr.ndim == 1:
        Sr = np.diag(Sr)
    return Ur @ Sr @ VTr

def truncEvaluate(U, S, VT, r):
    if S.ndim == 1:
        S = np.diag(S)
    return U[:,:r] @ S[:r,:r] @ VT[:r,:]


def svd_regression(b, U, S, VT):
    """ """
    xtilde = VT.T @ np.linalg.inv(np.diag(S)) @ U.T @ b 
    return xtilde


def plotS(S):
    """ 
    S:diagonal matrix of singular values
    """
    import matplotlib.pyplot as plt

    if S.ndim == 1:
        s=S
    else:
        s = np.diag(S)

    fig,axes = plt.subplots(1, 2, sharey=False, figsize=(12.0,4.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.30)
    ax = axes[0]
    ax.plot(s , label='')
    ax.set_yscale('log')
    ax.set_xlabel('Index')
    ax.set_ylabel('Singular value magnitude')

    ax = axes[1]
    ax.plot(np.cumsum(s)/np.sum(s))
    ax.set_xlabel('Index')
    ax.set_ylabel('Contribution to total sum')
    return axes
