""" 
Frequency response of a mass spring damper system, put into state space form, 
given some misc forcing inputs

"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Welib
from welib.system.statespacelinear import *
import sympy as sp

# Parameters defining the system
m = 250.0           # system mass
k = 40.0            # spring constant
b = 60.0            # damping constant

# System matrices
A = [[0, 1.], [-k/m, -b/m]]
B = [[0,0], [1/m, 2/m]]
C = [[1., 0],  # ouput 1 is position
     [0., 1]]  # ouput 2 is velocity
D=B


# --- Initialize a linear time invariant system
sys= LinearStateSpace(A,B,C,D)

# --- Compute frequency response
# G: magnitude, phi: phase
omega  = np.linspace(0.01,10, 1000)
G, phi = sys.frequency_response(omega)
f      = omega/(2*np.pi)

# --- Bode Plot
fig,axes = plt.subplots(2, 1, sharex=True, figsize=(6.4,4.8)) # (6.4,4.8)
fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
ax =axes[0]
ax.plot(f, G[0,0],       label=r'00  $u_1$ -> $x $'     )
ax.plot(f, G[1,0],       label=r'10  $u_1$ -> $\dot{x}$')
ax.plot(f, G[0,1], '--', label=r'01  $u_2$ -> $x $'     )
ax.plot(f, G[1,1], '--', label=r'11  $u_2$ -> $\dot{x}$')
ax.set_ylabel('Amplitude')
ax.set_yscale('log')
ax.set_xscale('log')
ax =axes[1]
ax.plot(f, phi[0,0]*180/np.pi      , label=r'00  $u_1$ -> $x $'     )
ax.plot(f, phi[1,0]*180/np.pi      , label=r'10  $u_1$ -> $\dot{x}$')
ax.plot(f, phi[0,1]*180/np.pi, '--', label=r'01  $u_2$ -> $x $'     )
ax.plot(f, phi[1,1]*180/np.pi, '--', label=r'11  $u_2$ -> $\dot{x}$')
ax.set_xscale('log')
ax.set_xlabel('Frequency [Hz]')
ax.set_ylabel('Phase [deg]')
ax.set_title('System - 2nd order -Frequency response LTI')
for ax in axes:
    ax.tick_params(direction='in')

ax.legend()




if __name__ == '__main__':
    plt.show()
if __name__ == '__test__':
    pass
if __name__=="__export__":
    pass
    #from welib.tools.repo import export_figs_callback
    #export_figs_callback(__file__)
