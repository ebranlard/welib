"""
Plot the fully separated, inviscid and sepeartion funciton
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Local 
from welib.airfoils.Polar import Polar
from welib.tools.colors import fColrs


scriptDir=os.path.dirname(__file__)


def main_sepf(test=False):
    polarFile_in = os.path.join(scriptDir,'../data/FFA-W3-241-Re12M.dat')
#     polarFile_in = os.path.join(scriptDir,'../data/63-235.csv')

    pol = Polar(polarFile_in, compute_params=True, verbose=False)
    print(pol)
    #  - method: 'max', 'optim', 'leastsquare', 'leastsquare_constraint'
    cl_fs, f_st = pol.cl_fully_separated(method='optim')
    #f_st=(P.cl-cl_fs)/(P.cl_inv-cl_fs);

    # --- Plot
    fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    ax.plot(pol.alpha, pol.cl     ,'-', c=fColrs(1), label= r'$C_{l,{st}}$ (steady)')
    ax.plot(pol.alpha, pol.cl_inv ,'-', c=fColrs(2), label= r'$C_{l,{inv}}$ (inviscid)')
    ax.plot(pol.alpha,      cl_fs ,'-', c=fColrs(3), label= r'$C_{l,{fs}}$ (fully-sep)')
    ax.plot(pol.alpha, f_st       ,'--', c=fColrs(4), label= r'$f_{sep}$ (separation function)')
#     ax.plot([pol.alpha0()], [0], 'ko')
#     ax.axvline(pol.alpha0)
    

    ax.tick_params(direction='in', top=True, right=True)
    ax.set_xlabel(r'Angle of attack, $\alpha$ [deg]')
#     ax.set_ylabel(r'Lift coefficient, $C_l$ [-]')
    ax.set_title(r'Airfoils - Separation function')
    ax.set_xlim([-50,50])
    ax.set_ylim([-1.5,2])
    ax.legend()

    return pol

main_sepf()

if __name__ == '__main__':
    plt.show()
if __name__ == '__test__':
    pass
if __name__=="__export__":
    from welib.tools.repo import export_figs_callback
    export_figs_callback(__file__)
