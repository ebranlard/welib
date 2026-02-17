""" 
Deflection of a cantilever beam with bending in two directions, arbitrary stiffness distribution and structural coupling the directions
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import welib.weio as weio
import welib.beams.cantilever2d as c2d
from welib.BEM.steadyBEM import SteadyBEM 

def main(plot=True):
    scriptDir = os.path.dirname(__file__)

    # --- Beam Data 
    df = weio.read(os.path.join(scriptDir, '../../../data/NREL5MW/5MW_Baseline/NRELOffshrBsline5MW_Blade.dat')).toDataFrame()
    R = 63
    r_hub =1.5
    r    = df['BlFract_[-]']*(R-r_hub)  + r_hub
    beta = df['StrcTwst_[deg]']*np.pi/180
    m    = df['BMassDen_[kg/m]']
    EI1  = df['FlpStff_[Nm^2]']
    EI2  = df['EdgStff_[Nm^2]']

    # --- Loading data
    MainFASTFile           = os.path.join(scriptDir,'../../../data/NREL5MW/Main_Onshore.fst')
    BEM = SteadyBEM(filename=MainFASTFile) # Initialize based on OpenFAST parameters, done on AeroDyn grid
    out = BEM.calcOutput(Omega=12.1, pitch=0, V0=10)
    df = out.radialDataFrame(r_new=r) # reinterpolate to structural nodes NOTE: we should remove hub..
    pz = df['fn_[N/m]']
    py = df['ft_[N/m]']

    # --- Compute Deflections
    uy, uz = c2d.deflection(py, pz, r, beta, EI1, EI2)

    # --- Plot loading and deflection
    if plot:
        fig,axes = plt.subplots(2, 1, sharex=True, sharey=False, figsize=(6.4,4.8))
        fig.subplots_adjust(left=0.12, right=0.95, top=0.92, bottom=0.11, hspace=0.20, wspace=0.10)
        axes[1].plot(r/R, pz, label=r'$p_n$')
        axes[1].plot(r/R, py, label=r'$p_t$')
        axes[0].plot(r/R, uz, label=r'$u_n$')
        axes[0].plot(r/R, uy, label=r'$u_t$')
        for ax in axes:
            ax.tick_params(direction='in', top=True, right=True, labelright=False, labeltop=False, which='both')
            ax.legend()
        axes[1].set_xlabel(r'Beam span [-]')
        axes[0].set_ylabel(r'Beam deflection $u$ [-]')
        axes[1].set_ylabel(r'Beam loading $p$ [N/m]')
        fig.suptitle('Beam - 2D - NREL5MW deflections')
    return uy, uz


if __name__ == '__main__':
    main(plot=True)
    plt.show()

if __name__ == '__test__':
    main(plot=False)

if __name__=='__export__':
    main(plot=True)
    from welib.tools.repo import export_figs_callback
    export_figs_callback(__file__)
