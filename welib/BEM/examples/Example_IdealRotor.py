""" 
Blade planform for "ideal rotor" without wake rotation (a'=0) or with wake rotation

Equations obtained by simplifying the BEM equations (Cd=0, F=1)
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Local 
from welib.BEM.idealrotors import *

TSR = 7      # Design tip speed ratio
Cl  = 1      # Design lift coefficient
alpha = 7    # Design angle of attack
B   = 3      # number of blades
R   = 100    # Rotor radius
r_hub = 0.1*R # hub radius


r     = np.linspace(r_hub, R, 100)
r_bar = r/R

# No Wake Rot
chord, phi = planform_nowakerot(r, R, TSR, Cl, B)
phi   = phi*180/np.pi
twist = phi-phi[-1]
chord_nw = chord
twist_nw = twist
# Wake Rot
chord, phi = planform_wakerot(r, R, TSR, Cl, B)
phi    = phi*180/np.pi
pitch0 = phi[-1] # pitch is set such that twist is 0 at blade tip, it's a choice 
twist  = phi-pitch0
chord_wr = chord
twist_wr = twist



#fig,axes = plt.subplots(1, 1, sharey=False, figsize=(10.4,3.9)) # (6.4,4.8)
#fig.subplots_adjust(left=0.12, right=0.98, top=0.98, bottom=0.15, hspace=0.10, wspace=0.13)
fig,axes = plt.subplots(2, 1, sharex=True, sharey=False, figsize=(6.4,3.9)) # (6.4,4.8)
fig.subplots_adjust(left=0.09, right=0.97, top=0.94, bottom=0.15, hspace=0.17, wspace=0.13)
ax=axes[0]
ax.plot(r_bar, chord_nw/R ,'k'    , label='Without wake rotation')
ax.plot(r_bar, chord_wr/R ,'k--'  , label='With wake rotation')
ax.set_ylabel('Chord, c/R [-]')
ax.set_yticks([0,0.1,0.2,0.3])
ax.set_ylim([0,0.3])
ax.set_xlim([0,1])
ax.legend()
ax.set_title('BEM Theory - Ideal rotor planform')

ax=axes[1]
ax.plot(r_bar, twist_nw, 'k' , label='Without wake rotation')
ax.plot(r_bar, twist_wr, 'k--' , label='With wake rotation')
ax.set_ylabel('Twist angle [deg]')
ax.set_xlabel('Blade radius r/R [-]')
ax.set_xlim([0,1])
ax.set_ylim([0,40])

for ax in axes:
    ax.tick_params(direction='in', top=True, right=True)
    # ax.spines['right'].set_visible(False)
    # ax.spines['top'].set_visible(False)

if __name__ == '__main__':
    plt.show()
if __name__=="__export__":
    from welib.tools.repo import export_figs_callback
    export_figs_callback(__file__)
