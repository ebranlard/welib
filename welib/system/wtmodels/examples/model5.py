import numpy as np
import matplotlib.pyplot as plt
# Local 
from welib.system.wtmodels.model5CS import *
from welib.tools.eva import eigMCK
from welib.system.mbc import Coleman2Comp



def Campbell5CS(vOmega=None):
    if vOmega is None:
        vOmega = np.linspace(0,1,11)

    # --- System parameters, TODO
    kt  = 11369784    # [Nm/rad]
    om0 = 0.8*2*np.pi # [rad/s] Blade edgewise natural frequency
    kx  = 200000      # [N/m] Lateral spring constant for M
    ky  = 250000      # Vertical spring constant for M
    l   = 30          # [m] Distance between m and M
    m   = 500         # [kg] Mass of each rotating mass
    M   = 50000       # [kg] Nacelle/generalized tower mass
    kb = m*l**2*om0**2

    plane='XYneg'; ordering='decreasing'; method='analytical'
    #plane='XYneg'; ordering='decreasing'; method='numerical'


    # ---- Coleman
    f_NR = np.empty( (len(vOmega),  MNR2.shape[0]) )
    for i, Omega in enumerate(vOmega):
        MNR, CNR, KNR = systemMatricesNR(M, m, l, kb, kx, ky, Omega, plane=plane, ordering=ordering, method=method)
        MNR2 = MNR
        CNR2 = CNR
        KNR2 = KNR

        fd, zeta, Q, f0 = eigMCK(MNR2, CNR2, KNR2, method='state_space')
                                         
        f_NR[i,:] = f0

        Qprint = np.zeros(Q.shape)
        for iMode in range(Q.shape[1]):
            a0,a1,b1 = Q[:3,iMode]
            x        = Q[3,iMode]
            A0,ABW,AFW,phi0,phiBW,phiFW = Coleman2Comp(a0, a1, b1)
            Qprint[0,iMode] = abs(x)
            Qprint[1:4,iMode] = [A0,ABW,AFW]

    # --- Plot Campbell
    fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8))
    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)

    vOmega_RPM = vOmega* 30/np.pi
    vOmega_Hz = vOmega/(2*np.pi)

    COLRS = plt.rcParams['axes.prop_cycle'].by_key()['color']

    ax.plot(vOmega_Hz, f_NR[:,0], '-', c=COLRS[0], label='Tower side-side') 
    ax.plot(vOmega_Hz, f_NR[:,1], '-', c=COLRS[1], label='Tower fore-aft') 
    ax.plot(vOmega_Hz, f_NR[:,2], ':', c=COLRS[2], label='Backward whirling') 
    ax.plot(vOmega_Hz, f_NR[:,3], '-', c=COLRS[2], label='Symmetric') 
    ax.plot(vOmega_Hz, f_NR[:,4], '--', c=COLRS[2], label='Forward whirling') 
    ax.plot(vOmega_Hz,   vOmega/2/np.pi, ':', c='0.2', label='1P')
    ax.plot(vOmega_Hz, 3*vOmega/2/np.pi, '--', c='0.2', label='3P')
    ax.grid('on', ls=':')
    ax.set_xlabel(f'Rotational frequency [Hz]')
    ax.set_ylabel('Frequencies [Hz]')
    fig.tight_layout()


    return vOmega, f_NR, fig



if __name__ == '__main__':
    from scipy.linalg import eig
    vOmega, f_NR, fig = Campbell5CS()
    plt.show()
