import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

MyDir=os.path.dirname(__file__)

from welib.airfoils.Polar import * 

# --------------------------------------------------------------------------------}
# ---  
# --------------------------------------------------------------------------------{
def prescribed_oscillations():
    #FFA-W3-241 airfoil Dyna Stall
    P=Polar.fromfile(os.path.join(MyDir,'../data/FFA-W3-241-Re12M.dat'),compute_params=True)

    # Parameters
    omega       = 12.57
    T           = 2*np.pi/omega
    tau         = 0.08
    valpha_mean = [5,10,15,20,25,30,35]
    t_max       = 1.3*T                  # simulation length
    dt          = 0.01                   # time step
    XLIM        = [0,40]

    # Derived params
    vt       = np.arange(0,t_max,dt)
    Cl       = np.zeros((len(valpha_mean),len(vt)))
    Cl_num   = np.zeros((len(valpha_mean),len(vt)))
    valpha_t = np.zeros((len(valpha_mean),len(vt)))

    # Loop on alpham and time 
    for ia,alpham in enumerate(valpha_mean):
        valpha_t[ia,:]=alpham+2*np.sin(omega*vt)
        fs_prev = P.f_st_interp(alpham)# init with steady value

        def dyna_stall_oye(t,fs):
            """ d(fs)/dt = 1/tau (fs_st - fs) """
            alpha_t = np.interp(t, vt, valpha_t[ia,:])
            f_st    = P.f_st_interp  (alpha_t)
            return 1/tau *( f_st - fs)

        # Integration using solve_ivp
        sol = solve_ivp(dyna_stall_oye, t_span=[0, t_max], y0=[fs_prev], t_eval=vt)
        for it,fs in enumerate(sol.y[0,:]):
            alpha_t = valpha_t[ia,it]
            Clinv   = P.cl_inv_interp(alpha_t)
            Clfs    = P.cl_fs_interp (alpha_t)
            Cl_num[ia,it]=fs*Clinv+(1-fs)*Clfs               

        # Simple discrete integration
        for it,alpha_t in enumerate(valpha_t[ia,:]):
            Cl[ia,it],fs_prev = P.dynaStallOye_DiscreteStep(alpha_t,tau,fs_prev,dt)

    fig=plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(P.alpha  , P.cl  , label='Cl static', LineWidth=2)
    for ia,alpham in enumerate(valpha_mean):
        if ia==0:
            lbl='Cl dynamic'
        else:
            lbl=''
        ax.plot(valpha_t[ia,:], Cl_num[ia,:], 'k-'   , label=lbl)
        #ax.plot(valpha_t[ia,:], Cl[ia,:], 'k--'   , label=lbl)
    ax.set_xlabel('Alpha [deg]')
    ax.set_ylabel('Cl [-]')
    ax.set_xlim(XLIM)
    ax.set_ylim([0,2.2])
    ax.set_title('Oye dynamic stall model')

if __name__ == '__main__':
    prescribed_oscillations()
    plt.show()
if __name__ == '__test__':
    prescribed_oscillations()
if __name__=="__export__":
    prescribed_oscillations()
    from welib.tools.repo import export_figs_callback
    export_figs_callback(__file__)

