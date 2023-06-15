""" 
Tools for linear regression
"""
import numpy as np




def regression(A, b, method='svd'):
    """ 
    return xtilde such that A xtilde =~ b

    A: is array of shape n x m
    b: is array of shape n
    x: is array of shape m

    One can then compare Axtilde and b

    """
    if method=='svd':
        U, S, VT = np.linalg.svd(A, full_matrices = False)
        xtilde = VT.T @ np.linalg.inv(np.diag(S)) @ U.T @ b  # TODO improve without inv


    elif method=='pinv':
        xtilde = np.linalg.pinv(A) @ b

    else:
        raise NotImplementedError(method)

    return xtilde
