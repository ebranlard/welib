""" 

References:
  [1] Wind Energy Explained, 3rd edition, Chapter 3

"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

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

    return chord, phi


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
     - phi:   tiwst [rad] (array-like)
    """
    lambda_r = TSR_design * r/R

    phi   = 2./3. * np.arctan(1./(lambda_r))
    chord = 8*np.pi * r/(B*Cl_design) * (1-np.cos(phi))

    # a = 1/(1+4 np.sin(phi)**2 / (sigma * Cl * np.cos(phi))
    # ap= (1-3a)/(4a-1)

    return chord, phi




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
    chord, phi = planform_nowakerot(r, R, TSR, Cl, B)
    phi   = phi*180/np.pi
    twist = phi-phi[-1]
    chord_nw = chord
    twist_nw = twist
    # Wake Rot
    chord, phi = planform_wakerot(r, R, TSR, Cl, B)
    phi   = phi*180/np.pi
    twist = phi-phi[-1]
    chord_wr = chord
    twist_wr = twist


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




