""" 
4-step Runge-Kutta-Nystrom time integration of 2nd order differential system

f is a function with signature f(t,x) 

IMPORTANT: x is assumed to be the concatenation of [pos,vel]

"""
import numpy as np
from scipy.optimize import OptimizeResult as OdeResultsClass 

def rk_nystrom_step(f, t, pos, vel, acc, dt, nDOF):
    # --- Estimates at t+dt/2
    A=dt/2*acc;
    # Speed changed to v+A, position changed to x+b
    b=dt/2*(vel+A/2)
    B=dt/2 * f(t+dt/2, np.concatenate((pos+b, vel+A)) )
    B=B[nDOF:]
    # speed changed to V+B
    C=dt/2 * f(t+dt/2, np.concatenate((pos+b, vel+B)))
    C=C[nDOF:]
    # --- Estimates at t+dt
    # speed changed to v+2*C, position changed to x+d
    d=dt*(vel+C);
    D=dt/2 * f(t+dt , np.concatenate((pos+d , vel+2*C)))
    D=D[nDOF:]
    # final estimate
    pos = pos+dt*(vel+1/3*(A+B+C)) 
    vel = vel+1/3*(A+2*B+2*C+D)  
    acc = f(t+dt , np.concatenate((pos, vel)))
    acc = acc[nDOF:]
    return pos, vel, acc

def rk_nystrom_integrate(f, t, x0, acc0=None):
    nt = len(t)
    nx = len(x0)
    y  = np.zeros((nx, nt))
    y[:,0] = x0
    nDOF = int(nx/2)
    pos = x0[:nDOF]
    vel = x0[nDOF:]
    if acc0 is None:
        f0 = f(t[0], x0)
        #acc = np.zeros(nDOF) # TODO chose
        acc = f0[nDOF:]
    else:
        acc = acc0

    for k in range(nt-1):
        dt = t[k+1]-t[k]
        pos, vel, acc = rk_nystrom_step(f, t[k], pos, vel, acc, dt, nDOF)
        y[:nDOF,k+1] = pos
        y[nDOF:,k+1] = vel

    return OdeResultsClass(t=t, y=y) # To mimic result class of solve_ivp
