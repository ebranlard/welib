import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

MyDir=os.path.dirname(__file__)

import welib.airfoils
from welib.airfoils.Polar import * 
from welib.airfoils.DynamicStall import * 




def step_change():
    # --- 
    # We use a step from alpha0 to alpha0+2, testing mainly the circulatory response (history)
    # Oye's dynamic stall model will not give a proper response here:
    #  - We are in the linear region, f close to 1, resulting in mainly Cl_inv
    #  - There is no induction build-up (induction history) in the Oye's stall model itself
    #    Fs is continuous and progressive but not Cl since it uses the current alpha. 
    radians=True
    P=Polar.fromfile(os.path.join(MyDir,'../data/FFA-W3-241-Re12M.dat'),compute_params=True,to_radians=radians)
    if radians:
        deg_scale=np.pi/180
    else:
        deg_scale=1
    U0      = 10
    chord   = 0.1591
    alpha1  = P._alpha0 
    alpha2  = alpha1+2*deg_scale
    tau_oye = 3 * chord/U0
    tau_t   = np.linspace(0,30,1000)
    vt      = chord * tau_t / (2*U0)
       
    # Wagner function flatplate
    # Wagner - Garrick approximation (Garrick 1938)
    Cl_wag_Garr=(tau_t+2)/(tau_t+4);
    # Wagner - R.T Jones approximation (Jones 1938)
    #Cl_wag_Jones=1-0.165*npexp(-0.0455*tau_t)-0.335*np.exp(-0.3*tau_t);
    Cl_wag_Jones=wagner(tau_t, constants='Jones')

    # Oye's Parameters and Inputs
    p_oye = dynstall_oye_param_from_polar(P, tau=tau_oye)
    u_oye=dict()
    u_oye['alpha']     = lambda t: alpha1 if t<=0 else alpha2 

    # MHH Parameters and Inputs
    p = dynstall_mhh_param_from_polar(P, chord, constants='Jones')
    u=dict()
    u['U']         = lambda t: U0
    u['U_dot']     = lambda t: 0 
    u['alpha']     = lambda t: alpha1 if t<=0 else alpha2 
    u['omega']     = lambda t: 0
    u['alpha_34']  = u['alpha']

    # Steady values
    Cl_st1 = P.cl_interp(alpha1)
    Cl_st2 = P.cl_interp(alpha2)
    fs1    = P.f_st_interp(alpha1)        # init with steady value
    fs2    = P.f_st_interp(alpha2)        # init with steady value
    y0_oye = [0.] # <<<<<<<<<<<<<<<<<<<<<<<<< do not init to fs1

    Cl_oye  = np.zeros(len(vt))


    # --- MHH
    p['alpha0_in_x1x2']=True
    df_mhh1  = dynstall_mhh_sim(vt, u, p, prefix = '', method = 'continuous')
    df_mhh2  = dynstall_mhh_sim(vt, u, p, prefix = '', method = 'discrete')
    p['alpha0_in_x1x2']=False
    df_mhh3  = dynstall_mhh_sim(vt, u, p, prefix = '', method = 'continuous')
    df_mhh4  = dynstall_mhh_sim(vt, u, p, prefix = '', method = 'discrete')

    # --- Oye
    sol_oye = solve_ivp(lambda t,x: dynstall_oye_dxdt(t,x,u_oye,p_oye), t_span=[0, max(vt)], y0=y0_oye, t_eval=vt)
    for it,t in enumerate(vt):
        Cl_oye[it] = dynstall_oye_output(vt[it],sol_oye.y[0,it],u_oye,p_oye)

    print('Cl steady:',Cl_st1,Cl_st2,df_mhh1['Cl_[-]'][0],Cl_oye[0])
    print('Fs steady:',fs1,fs2)

    fig=plt.figure();exec("try:import pybra.figlib;pybra.figlib.fig_move(fig)\nexcept:pass")
    ax = fig.add_subplot(111)
    ax.plot(tau_t,Cl_wag_Jones,'k'  ,label='Wagner (Jones approx.)')
    ax.plot(tau_t,Cl_oye[:]/Cl_st2,'-' ,label = 'Cl dynamic (Oye)')
    ax.plot(tau_t[:],df_mhh1['Cl_[-]'][:]/Cl_st2 ,'--',label = 'Cl dynamic (MHH)')
    ax.plot(tau_t[:],df_mhh2['Cl_[-]'][:]/Cl_st2 ,'--',label = 'Cl dynamic (MHH, discrete)')
    #ax.plot(tau_t,sol_oye.y[0,:]   ,'-' ,label = 'Fs (Oye)')
    #ax.plot(tau_t,sol_mhh.y[3,:] ,'--',label = 'Fs (MHH)')
    ax.set_xlabel('Dimensionless time [-]')
    ax.set_ylabel('Cl [-]')
    #ax.plot(tau_t,Cl_wag_Garr ,'k--',label='Wagner - Garr')
    plt.ylim([0.3,1.1])
    plt.legend()

    return Cl_wag_Jones, df_mhh1, df_mhh2, df_mhh3, df_mhh4, Cl_st2



if __name__ == '__main__':
    Cl_wag_Jones, df_mhh1, df_mhh2, df_mhh3, df_mhh4, Cl_st2 = step_change()
    np.testing.assert_almost_equal(df_mhh1['Cl_[-]'][1:]/Cl_st2, Cl_wag_Jones[1:], 2)
    np.testing.assert_almost_equal(df_mhh2['Cl_[-]'][1:]/Cl_st2, Cl_wag_Jones[1:], 2)
    np.testing.assert_almost_equal(df_mhh3['Cl_[-]'][1:]/Cl_st2, Cl_wag_Jones[1:], 2)
    np.testing.assert_almost_equal(df_mhh4['Cl_[-]'][1:]/Cl_st2, Cl_wag_Jones[1:], 2)
    plt.show()
if __name__ == '__test__':
    step_change()
