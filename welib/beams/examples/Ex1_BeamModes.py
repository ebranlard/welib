"""
Compute and plots the analytical modes of a uniform cantilever beam with and without a top mass
"""
import numpy as np
from welib.beams.theory import *
import matplotlib.pyplot as plt


def main():
    L  = 100    
    EI = 1.86e+12
    m  = 8.82e+03
    w  = 10
    Mtop=100000

    fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)

    freq,x,U,V,K = UniformBeamBendingModes('unloaded-clamped-free',EI,m,A=1,L=L,nModes=6)
    ax.plot(x, U[0,:], '-', label='First mode')
    ax.plot(x, U[1,:], '--', label='Second mode')

    freq,x,U,V,K = UniformBeamBendingModes('unloaded-topmass-clamped-free',EI,m,A=1,L=L,nModes=6,Mtop=Mtop)
    ax.plot(x, U[0,:], '-', label='First mode, with top-mass')
    ax.plot(x, U[1,:], '--', label='Second mode, with top-mass')

    ax.set_xlabel('Beam span [m]')
    ax.set_ylabel('Deflection [m]')
    ax.legend()
    ax.tick_params(direction='in')
    ax.set_title('Beam - Analytical mode shapes of a beam')



main()

if __name__ == '__main__':
    plt.show()
if __name__ == '__test__':
    try:
        plt.close('all')
    except:
        pass
if __name__=="__export__":
    from welib.tools.repo import export_figs_callback
    export_figs_callback(__file__)

