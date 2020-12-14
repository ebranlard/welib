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
    DeltaAlpha = 10
    # Parameters
    omega       = 2*U0*K_omega/chord
    T           = 2*np.pi/omega 
    #omega=0
    tau_oye     = 4 * chord/U0
    #valpha_mean = [5,10,15,20,25,30,35,50,60,80,100,110,130,150]
    valpha_mean = [-12,-22,10]
    #valpha_mean = [20]
    t_max       = 1.6*T                  # simulation length
    dt          = T/500                   # time step
    #XLIM        = np.array([0,40])
    XLIM        = np.array([-40,40])
    YLIM        = np.array([-2,2])


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
    from welib.tools.colors import fColrs

    fig,ax = plt.subplots(1, 1, sharey=False, figsize=(8.4,5.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    ax.plot(P.alpha/deg_scale  , P.f_st ,     label='f_st')
    ax.plot(P.alpha/deg_scale  , P.cl  , 'k-',label='Cl static', LineWidth=2)
    ax.plot(P.alpha/deg_scale  , P.cl_fs  ,'k--',  label='Cl fully separated')
    ax.plot(P.alpha/deg_scale  , P.cl_inv ,'k-.',  label='Cl inviscid')
    iStart=np.argmin(np.abs(vt-(vt[-1]-T)))-1
    for ia,alpham in enumerate(valpha_mean):
        if ia==0:
            lbl1='Cl dynamic (Oye)'
            lbl2='Cl dynamic (MHH)'
        else:
            lbl1=''
            lbl2=''
        col=fColrs(ia+1)
        ax.plot(valpha_t[ia,iStart:]/deg_scale, Cl_oye[ia,iStart:], ':' , color=col  , label=lbl1)
        ax.plot(valpha_t[ia,iStart:]/deg_scale, Cl_mhh[ia,iStart:] ,'-', color=col  , label=lbl2)
    ax.tick_params(direction='in')
    ax.set_xlabel('Alpha [deg]')
    ax.set_ylabel('Cl [-]')
    plt.xlim(XLIM)
    plt.ylim(YLIM)
    plt.legend()

if __name__ == '__main__':
    prescribed_oscillations()
    plt.show()
