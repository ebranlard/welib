import numpy as np
import pandas as pd
from .Polar import Polar as Pol
from .DynamicStall import *

def dynstall_mhh_dxdt_smd(t, x, U, p):

    """ Time derivative of states for continous formulation for airfoil attached
    to 3-DOF Spring-Mass-Damper system at pitch axis"""

    U_dot = 0.0 #Seems to be zero for some reason - Not sure why

    # States
    x1=x[0] # Downwash memory term 1
    x2=x[1] # Downwash memory term 2
    x3=x[2] # Clp', Lift coefficient with a time lag to the attached lift coeff
    x4=x[3] # f'' , Final separation point function

    if (p['prescribed_oscillations']):
        x[4] = 0.0
        x[5] = 0.0
        x[6] = np.radians(50.0) + np.radians(10.0) * np.sin(2.0 * np.pi * 0.4 * t)
        x[7] = 0.0
        x[8] = 0.0
        x[9] = 2.0 * np.pi * 0.4 * np.radians(10.0) * np.cos(2.0 * np.pi * 0.4 * t)

    xflap = x[4]
    xedge = x[5]
    theta = x[6]
    xflap_dot = x[7]
    xedge_dot = x[8]
    omega = x[9]

    if (not p['prescribed_oscillations']):
        theta = x[6] + p['reference_aoa']

    alpha_34 = np.arctan2( U - omega * (0.75 - p['x_pitch']) * p['chord'] * np.sin(theta), omega * (0.75 - p['x_pitch']) * p['chord'] * np.cos(theta) )
    alpha_34 = theta
    cldyn, cddyn, cmdyn = dynstall_mhh_outputs_simple(t, x[:4], U, U_dot, omega, alpha_34, p)
    force_SMD = np.dot( p['force_transform'], 0.5 * p['rho'] * U**2 * p['chord'] * np.r_[ cldyn, cddyn, cmdyn * p['chord'] + (p['x_pitch'] - 0.25) * (cldyn * np.cos(theta) + cddyn * np.sin(theta) )  ])

    # Parameters
    alpha0 = p['alpha0']
    Cla    = p['Cla']
    c      = p['chord']
    A1     = p['A1']
    A2     = p['A2']
    b1     = p['b1']
    b2     = p['b2']
    F_st   = p['F_st']
    # Variables derived from inputs
    U  = max(U, 0.01)
    Tu = max(c/(2*U), 1e-4)                                     # Eq. 23
    Tf     = p['Tf0']*Tu  # OLD was twice: Tf = p['Tf0']*c/U
    Tp     = p['Tp0']*Tu  # OLD was twice: Tp = p['Tp0']*c/U
    # Variables derived from states
    if p['alpha0_in_x1x2']:
        alphaE  = alpha_34*(1-A1-A2)+ x1 + x2  # Eq. 12
    else:
        alphaE  = (alpha_34-alpha0)*(1-A1-A2)+ x1 + x2  + alpha0  # Eq. 12

#     alphaE = u['alphaE'](t) # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< HACK HACK TODO TODO TODO TODO TODO

    Clp     = Cla * (alphaE-alpha0) + np.pi * Tu * omega      # Eq. 13
    alphaF  = x3/Cla+alpha0                                   # p. 13
    fs_aF   = F_st(alphaF)                                    # p. 13
    if(fs_aF<0):
        print('Problematic fs:',fs_aF)
    x4 = np.clip(x4, 1e-16, 1.0) # Constraining x4 between 0 and 1 increases numerical stability
    # State equation
    xdot = [0]*10
    if p['alpha0_in_x1x2']:
        xdot[0] = -1/Tu * (b1 + c * U_dot/(2*U**2)) * x1 + b1 * A1 / Tu * alpha_34
        xdot[1] = -1/Tu * (b2 + c * U_dot/(2*U**2)) * x2 + b2 * A2 / Tu * alpha_34
    else:
        xdot[0] = -1/Tu * (b1 + c * U_dot/(2*U**2)) * x1 + b1 * A1 / Tu * (alpha_34-alpha0)
        xdot[1] = -1/Tu * (b2 + c * U_dot/(2*U**2)) * x2 + b2 * A2 / Tu * (alpha_34-alpha0)
    xdot[2] = -1/Tp                             * x3 + 1/Tp * Clp
    xdot[3] = -1/Tf                             * x4 + 1/Tf * fs_aF

    #Now for the derivatives of the SMD
    xdot[4] = xflap_dot
    xdot[5] = xedge_dot
    xdot[6] = omega

    smd_acc = np.dot(p['m_inv'], force_SMD) - np.dot(p['m_inv c'], np.r_[xflap_dot, xedge_dot, omega]) - np.dot(p['m_inv k'], np.r_[xflap, xedge, theta])
    xdot[7], xdot[8], xdot[9] = smd_acc[0], smd_acc[1], smd_acc[2]

    if (p['prescribed_oscillations']):
        xdot[7] = 0.0
        xdot[8] = 0.0
        xdot[9] = 2.0 * np.pi * 0.4 * np.radians(10.0) * np.sin(2.0 * np.pi * 0.4 * t)


    return xdot
