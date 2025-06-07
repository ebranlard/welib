""" 

See also welib.wt_theory.idealrotors

References:
  [1] Wind Energy Explained, 3rd edition, Chapter 3

"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

def ADwr_CPopt_a_lambdar_eq(a, lambda_r):
    """ 
    Actuator disk with wake rotation, for optimal CP we need a'(1-a) to be max
    leading to:
        lambda_r^2 = (1-a)(4a-1)^2 / (1-3a)
    or, rewritten as:
         16 * a**3 - 24 * a**2 + a * (9 - 3 * lambda_r**2) + lambda_r**2 -1 =0

    NOTE:
       The analytical solution is provided in Branlard 2017

       It should be basically:

        phi   = 2./3. * np.arctan(1./(lambda_r))
        chord = 8*np.pi * r/(B*Cl_design) * (1-np.cos(phi))
        sigma = B*chord / (2*np.pi*r)
        a = 1/(1+4*np.sin(phi)**2 / (sigma * Cl_design * np.cos(phi)) )

    """

    return 16 * a**3 - 24 * a**2 + a * (9 - 3 * lambda_r**2)+  lambda_r**2 -1

def planform_ClCd(r, R,  TSR_design, Cl_design, Cd_design, B=3):
    """ 
    NOTE:  a and aprime will be the same as planfrom_wakerot
           They are not changed by Cd or F becasue they are "actuator disc" inductions
    """

    lambda_r =r/R * TSR_design

    # Axial induction factor from optimal CP maximization of a'(1-a)
    a = np.zeros(len(r))
    for i, tsr in enumerate(lambda_r):
        a[i] = fsolve(lambda aa: ADwr_CPopt_a_lambdar_eq(aa, tsr), 0.3)[0]

    # analytical solution
    phi   = 2./3. * np.arctan(1./(lambda_r))
    a2 = 1/(1+ np.sin(phi)**2/(  (1-np.cos(phi))*np.cos(phi)   ) )
    np.testing.assert_almost_equal(a, a2, 13)

    
    # Tangential induction based on orthogonality of induced velocity for an AD
    aprime  = (1 - 3 * a) / (4 * a - 1)                     # Only true for optimal rotor
    #aprime2= 1/2*(-1+np.sqrt(1+4/lambda_r**2 * a * (1-a))) # Always true 
    #np.testing.assert_almost_equal(aprime, aprime2, 13)

#     phi    = np.arctan((1 - a) / ((1 + aprime) * lambda_r))   # radians




    cn     = Cl_design * np.cos(phi) + Cd_design * np.sin(phi)
    F      = (2/np.pi) * np.arccos(np.exp(-B / 2 * (R - r) / (r * np.sin(phi))))
    #c      = (8 * np.pi * R * F * a * lambda_r * (np.sin(phi))**2) / ((1 - a) * B * TSR_design * Cn)
    chord  = a/(1-a) * (8 * np.pi * F * r * (np.sin(phi))**2) / (B*cn)

    phi *= 180/np.pi            # [deg]
    return chord, phi,  a, aprime 


def planform_nowakerot(r, R, TSR_design, Cl_design, B=3):
    """ 
    Chord and twist distribution for ideal rotor without wake rotation

    INPUTS:
     - r         : radial positions (array-like)
     - R         : rotor radius (scalar)
     - TSR_design: design tip-speed ratio=Omega R/U0 (scalar)
     - Cl_design : design lift coefficient
     - B         : number of blades
    OUTPUTS:
     - chord: chord [m] (array-like)
     - phi:   tiwst [rad] (array-like)
    """
    r = np.asarray(r)
    lambda_r = TSR_design * r/R

    phi   = np.arctan(2./(3.*lambda_r))
    chord = 8*np.pi * r * np.sin(phi) / (3*B*Cl_design*lambda_r)

    a  = 1/3 * np.ones(len(r))
    ap = np.zeros(len(r))

    phi *= 180/np.pi            # [deg]
    return chord, phi, a, ap


def planform_wakerot(r, R, TSR_design, Cl_design, B=3):
    """ 
    Chord and twist distribution for ideal rotor with wake rotation

    INPUTS:
     - r         : radial positions (array-like)
     - R         : rotor radius (scalar)
     - TSR_design: design tip-speed ratio=Omega R/U0 (scalar)
     - Cl_design : design lift coefficient
     - B         : number of blades
    OUTPUTS:
     - chord: chord [m] (array-like)
     - phi:   twist [rad] (array-like)
    """
    lambda_r = TSR_design * r/R

    phi   = 2./3. * np.arctan(1./(lambda_r))
    chord = 8*np.pi * r/(B*Cl_design) * (1-np.cos(phi))

    sigma = B*chord / (2*np.pi*r)

    a = 1/(1+4*np.sin(phi)**2 / (sigma * Cl_design * np.cos(phi)) )
    ap= (1-3*a)/(4*a-1) # Based on orthogonality and optimal a only!

    phi *= 180/np.pi            # [deg]
    return chord, phi, a, ap



if __name__ == '__main__':
    # --------------------------------------------------------------------------------}
    # --- Blade planform,  No Wake Rot  and With Wake rot
    # --------------------------------------------------------------------------------{
    TSR = 7
    B   = 3
    Cl  = 1
    R   = 100
    r_hub  = 0.1*R
    #r   = np.arange(0.1,1,0.05)*R
    r   = np.linspace(r_hub, R, 100)
    r_bar = r/R
    alpha =7

    # No Wake Rot
    chord_nw, phi_nw, a, ap = planform_nowakerot(r, R, TSR, Cl, B)
    twist_nw = phi_nw-phi_nw[-1] # Alternative: remove alpha_d!
    # Wake Rot
    chord_wr, phi_wr = planform_wakerot(r, R, TSR, Cl, B)
    twist_wr = phi_wr-phi_wr[-1] # Alternative: remove alpha_d!

    import matplotlib.pyplot as plt

    #fig,axes = plt.subplots(1, 1, sharey=False, figsize=(10.4,3.9)) # (6.4,4.8)
    #fig.subplots_adjust(left=0.12, right=0.98, top=0.98, bottom=0.15, hspace=0.10, wspace=0.13)
    fig,axes = plt.subplots(1, 1, sharey=False, figsize=(6.4,3.9)) # (6.4,4.8)
    fig.subplots_adjust(left=0.09, right=0.97, top=0.98, bottom=0.15, hspace=0.10, wspace=0.13)
    ax=axes
    ax.plot(r_bar, chord_nw/R ,'k'    , label='Without wake rotation')
    ax.plot(r_bar, chord_wr/R ,'k--'  , label='With wake rotation')
    ax.set_xlabel('Non-dimensionalized blade radius, r/R')
    ax.set_ylabel('Non-dimensionalized blade chord, c/R   ')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_yticks([0,0.1,0.2,0.3])
    ax.set_ylim([0,0.3])
    ax.set_xlim([0,1])
    ax.legend()

    fig,axes = plt.subplots(1, 1, sharey=False, figsize=(6.4,3.9)) # (6.4,4.8)
    fig.subplots_adjust(left=0.09, right=0.97, top=0.98, bottom=0.15, hspace=0.10, wspace=0.13)
    ax=axes
    ax.plot(r_bar, twist_nw, 'k' , label='Without wake rotation')
    ax.plot(r_bar, twist_wr, 'k--' , label='With wake rotation')
    ax.set_xlabel('Non-dimensionalized blade radius, r/R')
    ax.set_ylabel('Blade twist angle, deg')
    ax.set_xlim([0,1])
    ax.set_ylim([0,40])
    ax.legend()
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    fig.save()

    plt.show()




