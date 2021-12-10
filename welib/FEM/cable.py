"""
Pretension cable elements

Reference: 
    J. Jonkman, E. Branlard, M. Hall, G. Hayman, A. Platt, A. Robertson
    2020
    Implementation of Substructure Flexibility and Member-Level Load Capabilities for Floating Offshore Wind Turbines in OpenFAST
    NREL-TP-5000-76822
"""

import numpy as np
import scipy



def cable3d_Ke(L, A, E, T0=0, R=None, main_axis='z'):
    """ 
    Element stiffness  matrix for pretension cable

    INPUTS:
        L :    Element length [m]
        A :    Cross section area [m^2]
        T0:    Pretension [N]
        R:     Transformation matrix (3x3) from global coord to element coord: x_e = R.x_g
               if provided, element matrix is provided in global coord

    """
    Eps0 = T0/(E*A)
    L0   = L/(1+Eps0)  # "rest length" for which pretension would be 0
    EAL0 = E*A/L0
    EE   = EAL0* Eps0/(1+Eps0)

    Ke = np.zeros((12,12))

    # Note: only translational DOF involved (1-3, 7-9)
    Ke[0,0]= EE
    Ke[1,1]= EE
    Ke[2,2]= EAL0
    Ke[0,6]= -EE
    Ke[1,7]= -EE
    Ke[2,8]= -EAL0
    Ke[6,0]= -EE
    Ke[7,1]= -EE
    Ke[8,2]= -EAL0
    Ke[6,6]= EE
    Ke[7,7]= EE
    Ke[8,8]= EAL0

    if (R is not None):
        RR = scipy.linalg.block_diag(R,R,R,R)
        Ke = np.transpose(RR).dot(Ke).dot(RR)
        Ke = (Ke + Ke.T)/2 # enforcing symmetry 

    return Ke

def cable3d_Me(L, A, rho, R=None, main_axis='z'):
    """ Element mass matrix for pretension cable
    INPUTS:
        L :    Element length [m]
        A :    Cross section area [m^2]
        rho:   Material density [kg/m^3]
        R:     Transformation matrix (3x3) from global coord to element coord: x_e = R.x_g
               if provided, element matrix is provided in global coord
    """
    t = rho*A*L
    Me = np.zeros((12,12))

    Me[ 0,  0] = 13/35 * t
    Me[ 1,  1] = 13/35 * t
    Me[ 2,  2] = t/3
    Me[ 6,  6] = 13/35 * t
    Me[ 7,  7] = 13/35 * t
    Me[ 8,  8] = t/3
    Me[ 0,  6] =  9/70 * t
    Me[ 1,  7] =  9/70 * t
    Me[ 2,  8] = t/6
    Me[ 6,  0] =  9/70 * t 
    Me[ 7,  1] =  9/70 * t
    Me[ 8,  2] = t/6

    if (R is not None):
        RR = scipy.linalg.block_diag(R,R,R,R)
        Me = np.transpose(RR).dot(Me).dot(RR)
        Me = (Me + Me.T)/2 # enforcing symmetry 

    return Me


def cable3d_Fe_T0(T0, R=None, main_axis='z'):
    """ 
    Element stiffness  matrix for pretension cable

    INPUTS:
        T0:    Pretension [N]
        R:     Transformation matrix (3x3) from global coord to element coord: x_e = R.x_g
               if provided, element matrix is provided in global coord
    """
    Fe = np.zeros(12)
    Fe[2] = +T0  
    Fe[8] = -T0 

    if R is not None:
        Fe[0:3] = R.transpose().dot(Fe[0:3])
        Fe[6:9] = R.transpose().dot(Fe[6:9])

    return Fe

