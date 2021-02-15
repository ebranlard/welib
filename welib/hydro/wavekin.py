import numpy as np
import pandas as pd
from scipy.optimize import fsolve


def wavenumber(f, h, g=9.81):   
    """ solves the dispersion relation, returns the wave number k
    INPUTS:
      omega: wave cyclic frequency [rad/s], scalar or array-like
      h : water depth [m]
      g: gravity [m/s^2]
    OUTPUTS:
      k: wavenumber
    """
    omega = 2*np.pi*f
    if hasattr(omega, '__len__'): 
        k = np.array([fsolve(lambda k: om**2/g - k*np.tanh(k*h), (om**2)/g)[0] for om in omega])
    else:
        func = lambda k: omega**2/g - k*np.tanh(k*h)
        k_guess = (omega**2)/g 
        k = fsolve(func, k_guess)[0]
    return k

# Functions 
def elevation2d(a, f, k, eps, t, x=0):
    """  wave elevation (eta) 
    INPUTS:
      a : amplitudes,  scalar or array-like of dimension nf
      f : frequencies, scalar or array-like of dimension nf 
      k : wavenumbers, scalar or array-like of dimension nf
      t : time, scalar or array-like of dimension nt
      x : longitudinal position, scalar
    """ 
    t = np.asarray(t)
    a = np.asarray(a)
    omega = 2*np.pi * np.asarray(f)
    eta = np.zeros(len(t))
    for ai,oi,ki,ei in zip(a,omega,k,eps):
        eta += ai * np.cos(oi*t - ki*x + ei)
    return eta     

if __name__ == '__main__':
    pass
