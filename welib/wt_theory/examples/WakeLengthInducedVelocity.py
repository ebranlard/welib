# --- General
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
# --- Local
from welib.vortilib.elements.VortexCylinder import vc_tang_u, cylinder_tang_u


def induction_vs_length(vL,Zcp=0, gamma_t=-1, yGround=None):
    """ Returns uzL/uz0 induction for a wake of length L compared to infinite (0)"""
    Xcp=0.000001
    Ycp=0
    vL = np.asarray(vL)
    uzL =np.zeros(vL.shape)
    _,uz0     = vc_tang_u       (Xcp,Ycp,Zcp,gamma_t=gamma_t,R=1          ,polar_out=True,epsilon=0)
    for iL, L in enumerate(vL):
        _,uzL[iL] =  cylinder_tang_u(Xcp,Ycp,Zcp,gamma_t=gamma_t,R=1,z1=0,z2=L,polar_out=True,epsilon=0)
    # Adding ground
    if yGround is not None:
        Ycp=Ycp-yGround
        _,uzg     = vc_tang_u       (Xcp,Ycp,Zcp,gamma_t=gamma_t,R=1          ,polar_out=True,epsilon=0)
        uz0+=uzg
        uzLg =np.zeros(vL.shape)
        for iL, L in enumerate(vL):
            _,uzLg[iL] =  cylinder_tang_u(Xcp,Ycp,Zcp,gamma_t=gamma_t,R=1,z1=0,z2=L,polar_out=True,epsilon=0)
        uzL+=uzLg

    return uzL, uz0 ,uzL/uz0*100


def plotInduction(U0=1, yGround=-1.5, rel='Abs'):
    vL_inR = np.linspace(0,25,100) # in R
    uzL,uz0             , r_00 = induction_vs_length(vL_inR, Zcp = 0)
    # uzL_2p5D,uz0_2p5D   , r_25 = induction_vs_length(vL, Zcp = -5)
    # uzLg    ,uz0g       , r_00g= induction_vs_length(vL, Zcp = 0 , yGround=-1.5)
    # uzL_2p5Dg,uz0_2p5Dg , r_25g= induction_vs_length(vL, Zcp = -5, yGround=-1.5)

    fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    xCrits_inR=np.array([3,4,6])*2
    u = uzL/uz0*100
    ax.plot(vL_inR/2, u, 'k',label='')
    for xCrit_inR in xCrits_inR: 
        rCrit = np.interp(xCrit_inR, vL_inR, u)
        ax.plot(np.array([0        ,xCrit_inR])/2      ,[rCrit,rCrit], 'k--', lw=0.5)
        ax.plot(np.array([xCrit_inR,xCrit_inR])/2,[0   , rCrit], 'k--', lw=0.5)
        ax.text(0.5, rCrit+0.04, '{:.1f}%'.format(rCrit))
    # ax.legend()
    ax.set_xlim([0,12])
    ax.set_xlabel('Wake length x/D [-]')
    if rel =='Rel':
        ax.set_ylabel(r'Rel. error in induced velocity $(u_\infty-u_L)/u_\infty$ [%]')
        ax.set_yscale('log')
        ax.set_ylim([0,5])
    else:
        ax.set_ylabel(r'Induced velocity fraction $u_L/u_\infty$ [%]')
        ax.set_ylim([95,100])
    ax.tick_params(direction='in', top=True, right=True)
    ax.set_title('WT Theory - Induced velocity vs Wake length')

if __name__ == '__main__':
    #fig.savefig('WakeLengthRotorInduction{}.pdf'.format(rel))
    plotInduction()
    plt.show()

if __name__=="__export__":
    plotInduction()

    from welib.tools.repo import export_figs_callback
    export_figs_callback(__file__)


