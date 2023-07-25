"""
Functions for verifying nalu wind outputs
"""

import numpy as np
from scipy.linalg import eigh

import yaml
import matplotlib.pyplot as plt

def modal_time_series(Mmat, Kmat, Cmat, x0, thist):
    """
    Calculate a time series response with modal analysis for initial displacement
    Only considers unforced.

    Inputs:
      Mmat - Mass matrix
      Kmat - Stiffness matrix
      Cmat - Damping matrix
      x0 - initial displacements
      thist - time history to evaluate at

    Outputs:
      xhist - displacement history

    NOTES:
      1. Assumes only initial displacement and zero initial velocity
      2. Assumes diagonal modal damping. May be incorrect unless Cmat is appropriately constructed.
    """

    sub_inds = [0, Mmat.shape[0]-1]

    eigvals, eigvecs = eigh(Kmat, Mmat, subset_by_index=sub_inds)
    
    # Spectral Matrix
    Lambda = np.diag(eigvals)
    ModalDamping = eigvecs.T @ Cmat @ eigvecs
    
    # frequencies
    omega = np.sqrt(eigvals)
    zeta = np.diag(ModalDamping)/2/omega
    
    # damped frequency
    omega_d = omega * np.sqrt(1 - zeta**2)

    # Initial Modal Displacements
    modal_q0 = eigvecs.T @ Mmat @ x0

    qhist = np.zeros((x0.shape[0], thist.shape[0]))
    qdothist = np.zeros((x0.shape[0], thist.shape[0]))

    # Evaluate sdof response
    for i in range(x0.shape[0]):
        qhist[i, :] += np.exp(-zeta[i]*omega[i]*thist)*modal_q0[i]*np.cos(omega_d[i]*thist)

        qdothist[i, :] += np.exp(-zeta[i]*omega[i]*thist)*modal_q0[i] \
                           *(-np.cos(omega_d[i]*thist)*omega[i]*zeta[i] \
                             -omega_d[i]*np.sin(omega_d[i]*thist) )
    
    xhist = eigvecs @ qhist
    vhist = eigvecs @ qdothist

    return xhist,vhist

def load_gen_time(struct_file, x0, thist):

    iea15mw_sviv2d = yaml.load(open(struct_file),Loader=yaml.UnsafeLoader)
    Mmat = np.array(iea15mw_sviv2d['mass_matrix']).reshape(3,3)
    Kmat = np.array(iea15mw_sviv2d['stiffness_matrix']).reshape(3,3)
    Cmat = np.array(iea15mw_sviv2d['damping_matrix']).reshape(3,3)

    xhist,vhist = modal_time_series(Mmat, Kmat, Cmat, x0, thist)

    return xhist,vhist

def plot_compare(t, xytheta, xhist):
    """
    Plot comparison between two different time series data sets
    """

    print(t.shape)
    print(xytheta.shape)
    print(xhist.shape)

    fig,ax = plt.subplots(1, 1, sharey=False, figsize=(8.4,5.8)) # (6.4,4.8)

    ax.plot(t, xytheta[0, :], label='Numerical x')
    ax.plot(t, xhist[0, :], '--', label='Analytical x')

    ax.plot(t, xytheta[1, :], label='Numerical y')
    ax.plot(t, xhist[1, :], '--', label='Analytical y')
    ax.legend()

    fig,ax = plt.subplots(1, 1, sharey=False, figsize=(8.4,5.8)) # (6.4,4.8)

    ax.plot(t, xytheta[2, :], label='Numerical Theta')
    ax.plot(t, xhist[2, :], '--', label='Analytical Theta')
    ax.legend()
