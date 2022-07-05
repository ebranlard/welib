""" 
Amplitude response to force vibrations for different damping values
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Welib
from welib.system.singledof import *
import sympy as sp


vzeta=[1,0.5,0.375,0.25,0.15,0.0]
# zeta=[1,0.5,0.3,0.18,0.13,0.0]

frat=np.linspace(0,5,100)
F0_over_k=1

fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)


vsty=['-',':',':','--']

vdashes= [
        [1,0],
        [1,1],
        [3,2],
        [4,4],
        [4,2,1,2],
        [5,3,2,3],
        ]

for i,zeta in enumerate(vzeta):

    H0,phi = forced_vibration_particular_cst(frat, F0_over_k, zeta)

    ax.plot(frat, H0 ,  dashes=vdashes[i] ,c='k'  , label=r'$\zeta = {:.2f}$'.format(zeta))
ax.set_xlabel(r'Frequency ratio $\omega/\omega_n$')
ax.set_ylabel(r'Non-dimensional amplitude, xk/F')
ax.set_xlim([0,5])
ax.set_ylim([0,3.5])
ax.legend()
ax.set_title('System - 2nd order - forced vibrations')




if __name__ == '__main__':
    plt.show()
if __name__ == '__test__':
    pass
if __name__=="__export__":
    from welib.tools.repo import export_figs_callback
    export_figs_callback(__file__)
