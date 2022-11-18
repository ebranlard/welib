import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

MyDir=os.path.dirname(__file__)

import welib.airfoils
from welib.airfoils.Polar import * 
from welib.airfoils.DynamicStall import * 

# --------------------------------------------------------------------------------}
# ---  
# --------------------------------------------------------------------------------{
def prescribed_oscillations_ris_r_1792():
    """
        See Riso-R-1792 report, p23-26 Figure 4.1 4.2
    """
    radians=True
    #FFA-W3-241 airfoil Dyna Stall
    P=Polar(os.path.join(MyDir,'../data/DU21_A17.csv'), compute_params=True, radians=radians)

    if radians:
        deg_scale=np.pi/180
    else:
        deg_scale=1

    #print('alpha0:',P._alpha0/deg_scale)

    K_omega    = 0.05    # reduced frequency k=omega*c/(2U0)
    U0         = 1
    chord      = 0.1591 # gives a omega of 12.57 with k=0.1, U0=10
    DeltaAlpha = 4.5
    # Parameters
    omega       = 2*U0*K_omega/chord
    T           = 2*np.pi/omega 
    tau_oye     = 4 * chord/U0
    valpha_mean = [0,10]
    valpha_mean = [-.5,9.5]
    t_max       = 1.3*T                  # simulation length
    dt          = 0.01                   # time step
    XLIM        = np.array([0,40])


    # Derived params
    vt       = np.arange(0,t_max,dt)
    Cl_mhh   = np.zeros((len(valpha_mean),len(vt)))
    Cd_mhh   = np.zeros((len(valpha_mean),len(vt)))
    Cm_mhh   = np.zeros((len(valpha_mean),len(vt)))
    Cl_oye   = np.zeros((len(valpha_mean),len(vt)))
    valpha_t = np.zeros((len(valpha_mean),len(vt)))

    # Loop on alpham and time 
    for ia,alpham in enumerate(valpha_mean):
        valpha_t[ia,:]   = (alpham+DeltaAlpha*np.sin(omega*vt))*deg_scale
        valpha_dot_t     = (2*omega*np.cos(omega*vt) )*deg_scale 
        fs_prev = P.fs_interp(alpham*deg_scale)# init with steady value

        # Oye's Parameters and Inputs
        p_oye = dynstall_oye_param_from_polar(P, tau=tau_oye)
        u_oye=dict()
        u_oye['alpha']    = lambda t: np.interp(t, vt, valpha_t[ia,:])

        # MHH Parameters and Inputs
        u=dict()
        p = dynstall_mhh_param_from_polar(P, chord, constants='Jones')
        u['U']         = lambda t: U0
        u['U_dot']     = lambda t: 0
        u['alpha']     = lambda t: np.interp(t, vt, valpha_t[ia,:])
        u['omega']     = lambda t: np.interp(t, vt, valpha_dot_t)
        u['alpha_34']  = lambda t: np.interp(t, vt, valpha_t[ia,:]) # using alpha

        y0_mhh=[0,0,0,0]
        y0_oye=[0]
        y0_mhh = dynstall_mhh_steady(0,u,p)
        y0_oye = [fs_prev]



        # Oye - Integration using solve_ivp
        sol_oye = solve_ivp(lambda t,x: dynstall_oye_dxdt(t,x,u_oye,p_oye), t_span=[0, max(vt)], y0=y0_oye, t_eval=vt)
        for it,t in enumerate(vt):
            Cl_oye[ia,it] = dynstall_oye_output(vt[it],sol_oye.y[0,it],u_oye,p_oye)

        # Integration using solve_ivp
        sol_mhh = solve_ivp(lambda t,x: dynstall_mhh_dxdt(t,x,u,p), t_span=[0, t_max], y0=y0_mhh, t_eval=vt)
        for it,t in enumerate(vt):
            Cl_mhh[ia,it],Cd_mhh[ia,it],Cm_mhh[ia,it] = dynstall_mhh_outputs(t,sol_mhh.y[:,it],u,p)

        
    XLIM=[np.array([-5.10,4.10]),np.array([4.90,14.10])]
    # Cl
    YLIM=[np.array([-0.2,1.2]),np.array([1,1.7])]
    fig=plt.figure();exec("try:import pybra.figlib;pybra.figlib.fig_move(fig)\nexcept:pass")
    for ia,alpham in enumerate(valpha_mean):
        ax = fig.add_subplot(3,2,ia+1)
        ax.plot(P.alpha/deg_scale  , P.cl  , label='Cl static', linewidth=2)
        #ax.plot(valpha_t[ia,:]/deg_scale,Cl_oye[ia,:],'k--' ,label='Cl dynamic (Oye)')
        ax.plot(valpha_t[ia,:]/deg_scale,Cl_mhh[ia,:],'k-',label = 'Cl dynamic (MHH)')
        ax.set_ylabel('Cl [-]')
        ax.set_xlim(XLIM[ia])
        ax.set_ylim(YLIM[ia])
        ax.grid()
    # Cd
    YLIM=[np.array([-0.005,0.020]),np.array([-0.02,0.1])]
    for ia,alpham in enumerate(valpha_mean):
        ax = fig.add_subplot(3,2,ia+1+2)
        ax.plot(P.alpha/deg_scale  , P.cd  , label='Cd static', linewidth=2)
        ax.plot(valpha_t[ia,:]/deg_scale,Cd_mhh[ia,:],'k-',label = 'Cd dynamic (MHH)')
        ax.set_ylabel('Cd [-]')
        ax.set_xlim(XLIM[ia])
        ax.set_ylim(YLIM[ia])
        ax.grid()
    # Cm
    YLIM=[np.array([-0.15,-0.115]),np.array([-0.14,-0.08])]
    for ia,alpham in enumerate(valpha_mean):
        ax = fig.add_subplot(3,2,ia+1+4)
        ax.plot(P.alpha/deg_scale  , P.cm  , label='Static', linewidth=2)
        ax.plot(valpha_t[ia,:]/deg_scale,Cm_mhh[ia,:],'k-',label = 'Dynamic (MHH)')
        ax.set_xlabel('Alpha [deg]')
        ax.set_ylabel('Cm [-]')
        ax.set_xlim(XLIM[ia])
        ax.set_ylim(YLIM[ia])
        ax.grid()
    ax.legend()
    return Cl_mhh, Cd_mhh, Cm_mhh


if __name__ == '__main__':
    Cl_mhh, Cd_mhh, Cm_mhh = prescribed_oscillations_ris_r_1792()
    plt.show()
if __name__ == '__test__':
    Cl_mhh, Cd_mhh, Cm_mhh = prescribed_oscillations_ris_r_1792()
    np.testing.assert_almost_equal(Cl_mhh[0,1299], 0.9520, 3)
    np.testing.assert_almost_equal(Cd_mhh[0,1299], 0.0106, 3)
    np.testing.assert_almost_equal(Cm_mhh[0,1299],-0.1397, 3)
