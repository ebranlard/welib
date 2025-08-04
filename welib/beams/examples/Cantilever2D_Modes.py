""" 
Coupled modes of a cantilever beam with bending in two directions, arbitrary stiffness distribution and structural coupling the directions
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import welib.weio as weio
import welib.beams.cantilever2d as c2d


def main(plot=True):

    # --- Beam Data 
    scriptDir = os.path.dirname(__file__)
    df = weio.read(os.path.join(scriptDir, '../../../data/NREL5MW/5MW_Baseline/NRELOffshrBsline5MW_Blade.dat')).toDataFrame()
    R = 63
    r_hub =1.5
    r    = df['BlFract_[-]']*(R-r_hub)  + r_hub
    beta = df['StrcTwst_[deg]']*np.pi/180
    m    = df['BMassDen_[kg/m]']
    EI1  = df['FlpStff_[Nm^2]']
    EI2  = df['EdgStff_[Nm^2]']

    # --- Compute coupled modes
    n_modes = 4
    freqs, modes = c2d.compute_modes(n_modes, r, EI1, EI2, m, beta)

    # --- Plot coupled modes
    if plot:
        fig,axes = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(6.4,4.8))
        fig.subplots_adjust(left=0.12, right=0.95, top=0.85, bottom=0.11, hspace=0.20, wspace=0.10)
        for i, (f,m) in enumerate(zip(freqs,modes)):
            print(f'Mode {i+1}: Frequency = {freqs[i]:.3f} Hz')
            ii = i // 2  # Row index
            jj = i % 2   # Column index
            axes[ii,jj].plot(r/R, modes[i][0], label=r'$u_y$')
            axes[ii,jj].plot(r/R, modes[i][1], label=r'$u_z$')
            axes[ii,jj].set_title(f'Mode {i+1}')
            axes[ii,jj].tick_params(direction='in', top=True, right=True, labelright=False, labeltop=False, which='both')
        axes[0,0].legend()
        axes[1,0].set_xlabel(r'Beam span $x$ [-]')
        axes[1,1].set_xlabel(r'Beam span $x$ [-]')
        axes[0,0].set_ylabel(r'Beam deflection $u$ [-]')
        axes[1,0].set_ylabel(r'Beam deflection $u$ [-]')

        fig.suptitle('Beam - 2D - NREL5MW coupled blade modes')

    return freqs, modes


if __name__ == '__main__':
    freqs, modes = main(plot=True)
    plt.show()

if __name__ == '__test__':
    freqs, modes = main(plot=False)

if __name__=='__export__':
    freqs, modes = main(plot=True)
    from welib.tools.repo import export_figs_callback
    export_figs_callback(__file__)
