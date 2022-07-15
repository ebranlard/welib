"""
Compute and plots the analytical modes of a uniform beam for different boundary conditions
"""
import numpy as np
from welib.beams.theory import *
import matplotlib
import matplotlib.pyplot as plt


def main():

    BC=np.array([
        [ 'clamped-clamped', 'clamped-hinged', 'clamped-free', 'clamped-guided'],
        ['hinged-hinged', 'hinged-guided', 'guided-guided',None],
        ['free-free', 'free-hinged','free-guided',None]
        ])

    L  = 100    
    EI = 1.86e+12
    m  = 8.82e+03

    fig,axes = plt.subplots(BC.shape[0], BC.shape[1], sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.02, right=0.98, top=0.94, bottom=0.02, hspace=0.31, wspace=0.19)

    for i in range(BC.shape[0]):
        for j in range(BC.shape[1]):
            ax=axes[i,j]
            bc = BC[i,j]
            print('bc',bc)
            if bc is None:
                ax.set_axis_off()
                continue
            bc1,bc2 = bc.split('-')
            if bc2=='free' and bc1!='free':
                norm='tip'
            else:
                norm='max'
            # Compute frequencies and mode shapes for a given boundary condition
            freq,x,U,V,K = UniformBeamBendingModes('unloaded-'+bc,EI,m,A=1,L=L, nModes=3, norm=norm)
            ax.plot(x/L, U[0,:], '-',  label='First mode')
            ax.plot(x/L, U[1,:], '--', label='Second mode')
            ax.plot(x/L, U[2,:], '-.', label='Third mode')
            ax.set_title(bc, fontsize=10)
            ax.set_xlim([0,1])
            ax.set_ylim([-1,1])
            ax.set_xticks([])
            ax.set_yticks([])
            ax.tick_params(direction='in')
            if i==2 and j==2:
                ax.legend(loc='center left', bbox_to_anchor=(1.1, 0.5))
    fig.suptitle('Beam - Analytical mode shapes different BC')



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

