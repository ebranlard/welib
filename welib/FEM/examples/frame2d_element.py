"""
This examples plots the shape functions for a 2d frame element
and shows an example of deflection for some given nodal inputs.
"""

import numpy as np
from welib.FEM.frame2d import *
import matplotlib.pyplot as plt

def main():
    #fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    #fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    L=100

    vx = np.linspace(0,L,100)
    fig,ax = plt.subplots(1,1)
    ax.plot(vx/L, b2(vx/L),     '-',label=r'$N_2$')
    ax.plot(vx/L, b3(vx/L,L)/L, '--',label=r'$N_3/L$')
    ax.plot(vx/L, b5(vx/L),     '-',label=r'$N_5$')
    ax.plot(vx/L, b6(vx/L,L)/L, '--',label=r'$N_6/L$')
    ax.set_xlabel('x/L [-]')
    ax.set_ylabel('Transverse shape functions')
    ax.tick_params(direction='in')
    ax.legend()

    fig,ax = plt.subplots(1,1)
    ax.plot(vx/L, b1(vx/L), '-',  label=r'$N_1$')
    ax.plot(vx/L, b4(vx/L), '--', label=r'$N_4$')
    ax.set_xlabel('x/L [-]')
    ax.set_ylabel('Longitudinal shape functions')
    ax.tick_params(direction='in')
    ax.legend()

    fig,ax = plt.subplots(1,1)
    ax.plot(vx/L, h(vx, 1, 2/L, 3,-5/L,L),   label=r'$h(x)\ for\ q=(1,2/L,3,-5/L)$')
    ax.set_xlabel('x/L [-]')
    ax.set_ylabel('Beam transverse deflection [m]')
    ax.tick_params(direction='in')
    ax.legend()

    #fig,ax = plt.subplots(1,1)
    #ax.plot(vx/L, uax(vx,1,2,1),   label=r'$u(x)\ for\ q=(1,2)$')
    #ax.set_xlabel('')
    #ax.set_ylabel('')
    #ax.tick_params(direction='in')
    #ax.legend()
if __name__ == '__main__':
    main()
    plt.show()
if __name__ == '__test__':
    main()
