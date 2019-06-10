import unittest
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import pdb

MyDir=os.path.dirname(__file__)
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))

import welib.airfoils
from welib.airfoils.Polar import * 
from welib.airfoils.DynamicStall import * 

# Wagner - R.T Jones approximation (Jones 1938)
#Cl_wag_Jones=1-0.165*npexp(-0.0455*tau_t)-0.335*np.exp(-0.3*tau_t);
A1_Jones=0.165
A2_Jones=0.335
b1_Jones=0.0455
b2_Jones=0.3



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
#     A1_Jones=0.165
#     A2_Jones=0.335
#     b1_Jones=0.0455
#     b2_Jones=0.3

    Cl_wag_Jones=1-A1_Jones*np.exp(-b1_Jones*tau_t)-A2_Jones*np.exp(-b2_Jones*tau_t);

    # Oye's Parameters and Inputs
    p_oye = dynstall_oye_param_from_polar(P, tau=tau_oye)
    u_oye=dict()
    u_oye['alpha']     = lambda t: alpha1 if t<=0 else alpha2 

    # MHH Parameters and Inputs
    p = dynstall_mhh_param_from_polar(P, chord, tau_chord=chord/U0, Jones=True)
    u=dict()
    u['U']         = lambda t: U0
    u['U_dot']     = lambda t: 0 
    u['alpha']     = lambda t: alpha1 if t<=0 else alpha2 
    u['alpha_dot'] = lambda t: 0
    u['alpha_34']  = u['alpha']

    # Steady values
    Cl_st1 = P.cl_interp(alpha1)
    Cl_st2 = P.cl_interp(alpha2)
    fs1    = P.f_st_interp(alpha1)        # init with steady value
    fs2    = P.f_st_interp(alpha2)        # init with steady value
    y0_mhh = dynstall_mhh_steady(0,u,p)
    y0_oye = [0.] # <<<<<<<<<<<<<<<<<<<<<<<<< do not init to fs1

    Cl_mhh  = np.zeros(len(vt))
    Cl_oye  = np.zeros(len(vt))
    # Integration using solve_ivp
    sol_mhh = solve_ivp(lambda t,x: dynstall_mhh_dxdt(t,x,u,p), t_span=[0, max(vt)], y0=y0_mhh, t_eval=vt)
    for it,t in enumerate(vt):
        Cl_mhh[it],_,_ = dynstall_mhh_outputs(t,sol_mhh.y[:,it],u,p)

    # Integration using solve_ivp
    sol_oye = solve_ivp(lambda t,x: dynstall_oye_dxdt(t,x,u_oye,p_oye), t_span=[0, max(vt)], y0=y0_oye, t_eval=vt)
    for it,t in enumerate(vt):
        Cl_oye[it] = dynstall_oye_output(vt[it],sol_oye.y[0,it],u_oye,p_oye)

    print('Cl steady:',Cl_st1,Cl_st2,Cl_mhh[0],Cl_oye[0])
    print('Fs steady:',fs1,fs2)

    fig=plt.figure();exec("try:import pybra.figlib;pybra.figlib.fig_move(fig)\nexcept:pass")
    ax = fig.add_subplot(111)
    ax.plot(tau_t,Cl_wag_Jones,'k'  ,label='Wagner (Jones approx.)')
    ax.plot(tau_t,Cl_oye[:]/Cl_st2,'-' ,label = 'Cl dynamic (Oye)')
    ax.plot(tau_t[1:],Cl_mhh[1:]/Cl_st2 ,'--',label = 'Cl dynamic (MHH)')
    #ax.plot(tau_t,sol_oye.y[0,:]   ,'-' ,label = 'Fs (Oye)')
    #ax.plot(tau_t,sol_mhh.y[3,:] ,'--',label = 'Fs (MHH)')
    ax.set_xlabel('Dimensionless time [-]')
    ax.set_ylabel('Cl [-]')
    #ax.plot(tau_t,Cl_wag_Garr ,'k--',label='Wagner - Garr')
    plt.ylim([0.3,1.1])
    plt.legend()



# --------------------------------------------------------------------------------}
# ---  
# --------------------------------------------------------------------------------{
def prescribed_oscillations_ris_r_1792():
    """
        See Riso-R-1792 report, p23-26 Figure 4.1 4.2
        """
    radians=True
    #FFA-W3-241 airfoil Dyna Stall
    P=Polar.fromfile(os.path.join(MyDir,'../data/DU21_A17.csv'),compute_params=True,to_radians=radians)

    if radians:
        deg_scale=np.pi/180
    else:
        deg_scale=1

    print('alpha0:',P._alpha0/deg_scale)

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
        fs_prev = P.f_st_interp(alpham*deg_scale)# init with steady value

        # Oye's Parameters and Inputs
        p_oye = dynstall_oye_param_from_polar(P, tau=tau_oye)
        u_oye=dict()
        u_oye['alpha']    = lambda t: np.interp(t, vt, valpha_t[ia,:])

        # MHH Parameters and Inputs
        u=dict()
        p = dynstall_mhh_param_from_polar(P, chord, tau_chord=chord/U0, Jones=True)
        u['U']         = lambda t: U0
        u['U_dot']     = lambda t: 0
        u['alpha']     = lambda t: np.interp(t, vt, valpha_t[ia,:])
        u['alpha_dot'] = lambda t: np.interp(t, vt, valpha_dot_t)
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
        ax.plot(P.alpha/deg_scale  , P.cl  , label='Cl static', LineWidth=2)
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
        ax.plot(P.alpha/deg_scale  , P.cd  , label='Cd static', LineWidth=2)
        ax.plot(valpha_t[ia,:]/deg_scale,Cd_mhh[ia,:],'k-',label = 'Cd dynamic (MHH)')
        ax.set_ylabel('Cd [-]')
        ax.set_xlim(XLIM[ia])
        ax.set_ylim(YLIM[ia])
        ax.grid()
    # Cm
    YLIM=[np.array([-0.15,-0.115]),np.array([-0.14,-0.08])]
    for ia,alpham in enumerate(valpha_mean):
        ax = fig.add_subplot(3,2,ia+1+4)
        ax.plot(P.alpha/deg_scale  , P.cm  , label='Static', LineWidth=2)
        ax.plot(valpha_t[ia,:]/deg_scale,Cm_mhh[ia,:],'k-',label = 'Dynamic (MHH)')
        ax.set_xlabel('Alpha [deg]')
        ax.set_ylabel('Cm [-]')
        ax.set_xlim(XLIM[ia])
        ax.set_ylim(YLIM[ia])
        ax.grid()
    ax.legend()


# --------------------------------------------------------------------------------}
# ---  
# --------------------------------------------------------------------------------{
def prescribed_oscillations():
    radians=True
    #FFA-W3-241 airfoil Dyna Stall
#     P=Polar.fromfile(os.path.join(MyDir,'../data/FFA-W3-241-Re12M.dat'),compute_params=True,to_radians=radians)
    P=Polar.fromfile(os.path.join(MyDir,'../data/DU21_A17.csv'),compute_params=True,to_radians=radians)

    if radians:
        deg_scale=np.pi/180
    else:
        deg_scale=1

    K_omega    = 0.1    # reduced frequency k=omega*c/(2U0)
    U0         = 10
    chord      = 0.1591 # gives a omega of 12.57 with k=0.1, U0=10
    DeltaAlpha = 4
    # Parameters
    omega       = 2*U0*K_omega/chord
    T           = 2*np.pi/omega 
    #omega=0
    #T           = 20
    tau         = 0.08
    tau_oye     = 4 * chord/U0
    valpha_mean = [5,10,15,20,25,30,35,50,60,80,100,110,130,150]
    #valpha_mean = [20]
    t_max       = 1.3*T                  # simulation length
    dt          = 0.01                   # time step
    XLIM        = np.array([0,40])
    XLIM        = np.array([0,180])


    # Derived params
    vt       = np.arange(0,t_max,dt)
    Cl_mhh   = np.zeros((len(valpha_mean),len(vt)))
    Cl_oye   = np.zeros((len(valpha_mean),len(vt)))
    valpha_t = np.zeros((len(valpha_mean),len(vt)))

    # Loop on alpham and time 
    for ia,alpham in enumerate(valpha_mean):
        valpha_t[ia,:]   = (alpham+DeltaAlpha*np.sin(omega*vt))*deg_scale
        valpha_dot_t     = (2*omega*np.cos(omega*vt) )*deg_scale
        fs_prev = P.f_st_interp(alpham*deg_scale)# init with steady value

        # Oye's Parameters and Inputs
        p_oye = dynstall_oye_param_from_polar(P, tau=tau_oye)
        u_oye=dict()
        u_oye['alpha']    = lambda t: np.interp(t, vt, valpha_t[ia,:])

        # MHH Parameters and Inputs
        p = dynstall_mhh_param_from_polar(P, chord, tau_chord=chord/U0, FAST=True)
        u=dict()
        u['U']         = lambda t: U0
        u['U_dot']     = lambda t: 0
        u['alpha']     = lambda t: np.interp(t, vt, valpha_t[ia,:])
        u['alpha_dot'] = lambda t: np.interp(t, vt, valpha_dot_t)
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
            Cl_mhh[ia,it],_,_ = dynstall_mhh_outputs(t,sol_mhh.y[:,it],u,p)

        #print('alpham    ', alpham*deg_scale, y0_mhh[0]+y0_mhh[1])
        #print('fst oye   ', y0_oye)
        #print('steady mhh', y0_mhh)
        #print('Cl mhh t0 ', dynstall_mh_outputs(t,y0_mhh,u,p))
        #print('Cl mhh t0 ', Cl_mh[ia,0])
        #print('Cl oye t0 ', Cl_oye[ia,0])
        #print('Cl oye t0 ', P.cl_interp(alpham*deg_scale))


        #fig=plt.figure()
        #ax = fig.add_subplot(111)
#       #  ax.plot(vt,sol_mh.y[0,:],label='x1')
        #ax.plot(vt,sol_mh.y[0,:]+sol_mh.y[1,:],label='alphaE')
        ##  alphaF  = x3/Cla+alpha0                                  # p. 13
        #ax.plot(vt,sol_mh.y[2,:]/P._linear_slope+P._alpha0,label='alphaF')
        #ax.plot(vt,valpha_t[ia,:],label='alpha')
        #ax.plot(vt,sol_mh.y[3,:],label='x4')
        #ax.plot(vt,sol.y[0,:]  ,label='f_st')
        ##ax.plot(vt,Cl_mh[ia,:],label='Cl_mh')
        ##ax.plot(vt,Cl_oye[ia,:],label='Cl_oy')
        #fig.legend()
        #print(sol_mh.y[:,-1])
        #print(dynstall_mh_steady(0,u,p))
        #print(Cl_mh[ia,-1])
        #print(P.cl_interp(alpham*deg_scale))



    fig=plt.figure();exec("try:import pybra.figlib;pybra.figlib.fig_move(fig,'Center')\nexcept:pass")
    ax = fig.add_subplot(111)
    ax.plot(P.alpha/deg_scale  , P.cl  , label='Cl static', LineWidth=2)
    for ia,alpham in enumerate(valpha_mean):
        if ia==0:
            lbl1='Cl dynamic (Oye)'
            lbl2='Cl dynamic (MHH)'
        else:
            lbl1=''
            lbl2=''
        ax.plot(valpha_t[ia,:]/deg_scale, Cl_oye[ia,:], 'k-'   , label=lbl1)
        ax.plot(valpha_t[ia,:]/deg_scale, Cl_mhh[ia,:] , 'k--'   , label=lbl2)
    ax.set_xlabel('Alpha')
    ax.set_ylabel('Cl [-]')
    plt.xlim(XLIM)
    plt.ylim([0,2.2])
    plt.legend()

if __name__ == '__main__':
    #prescribed_oscillations_ris_r_1792()
    prescribed_oscillations()
    #step_change()
    plt.show()
