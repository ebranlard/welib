import sys, os, yaml
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy import linalg

MyDir=os.path.dirname(__file__)

import welib.airfoils
from welib.airfoils.Polar import *
from welib.airfoils.DynamicStall import *
from welib.airfoils.DynamicStall_SMD import *
# --------------------------------------------------------------------------------}
# ---
# --------------------------------------------------------------------------------{
def prescribed_oscillations():
    radians=True
    #FFA-W3-211 airfoil Dyna Stall
    P=Polar.fromfile(os.path.join(MyDir,'../data/IEA-15-240-RWT_AeroDyn15_Polar_35.dat'),compute_params=True,to_radians=radians)
    #FFA-W3-241 airfoil Dyna Stall
    #P=Polar.fromfile(os.path.join(MyDir,'../data/FFA-W3-241-Re12M.dat'),compute_params=True,to_radians=radians)
#    P=Polar(os.path.join(MyDir,'../data/DU21_A17.csv'), compute_params=True, radians=radians)

    if radians:
        deg_scale=np.pi/180
    else:
        deg_scale=1

    K_omega    = 2.0 * np.pi * 4.0 * 2.8480588747143942 / (2.0 * 15.0)     # reduced frequency k=omega*c/(2U0)
    U0         = 15.0
    chord      = 2.8480588747143942 # gives a omega of 12.57 with k=0.1, U0=10
    DeltaAlpha = 10
    # Parameters
    omega       = 2*U0*K_omega/chord
    T           = 10.0 #2*np.pi/omega
    #omega=0
    tau_oye     = 4 * chord/U0
    #valpha_mean = [5,10,15,20,25,30,35,50,60,80,100,110,130,150]
    valpha_mean = [50]
    #valpha_mean = [20]
    t_max       = T                  # simulation length
    dt          = T/500                   # time step
    #XLIM        = np.array([0,40])
    XLIM        = np.array([-90,90])
    YLIM        = np.array([-3,3])


    # Derived params
    vt       = np.arange(0,t_max,dt)
    Cl_mhh   = np.zeros((len(valpha_mean),len(vt)))
    Cl_oye   = np.zeros((len(valpha_mean),len(vt)))
    valpha_t = np.zeros((len(valpha_mean),len(vt)))

    # Loop on alpham and time
    for ia,alpham in enumerate(valpha_mean):
        valpha_t[ia,:]   = np.radians(alpham+DeltaAlpha*np.sin(omega*vt))
        valpha_dot_t     = np.radians(DeltaAlpha*omega*np.cos(omega*vt) )
        fs_prev = P.fs_interp(alpham*deg_scale)# init with steady value

        # Oye's Parameters and Inputs
        p_oye = dynstall_oye_param_from_polar(P, tau=tau_oye)
        u_oye=dict()
        u_oye['alpha']    = lambda t: np.interp(t, vt, valpha_t[ia,:])

        # MHH Parameters and Inputs
        p = dynstall_mhh_param_from_polar(P, chord, constants='OpenFAST')
        u=dict()
        u['U']         = lambda t: U0
        u['U_dot']     = lambda t: 0
        u['alpha']     = lambda t: np.interp(t, vt, valpha_t[ia,:])
        u['omega'] = lambda t: np.interp(t, vt, valpha_dot_t)
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
    ax.plot(P.alpha/deg_scale  , P.fs  ,     label='f_s')
    ax.plot(P.alpha/deg_scale  , P.cl  , 'k-',label='Cl static', linewidth=2)
    ax.plot(P.alpha/deg_scale  , P.cl_fs  ,'k--',  label='Cl fully separated')
    ax.plot(P.alpha/deg_scale  , P.cl_inv ,'k-.',  label='Cl inviscid')
    iStart=np.argmin(np.abs(vt-(vt[-1]-T)))-1
    iStart=0
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
    ax.set_xlim(XLIM)
    ax.set_ylim(YLIM)
    ax.legend()
    ax.set_title('Airfoils - MGH dynamic stall model')

def sviv_2d_prescribed_oscillations(u_infty, aoa_0):
    radians=True
    #FFA-W3-211 airfoil Dyna Stall
    P=Polar.fromfile(os.path.join(MyDir,'../data/IEA-15-240-RWT_AeroDyn15_Polar_35.dat'),compute_params=True,to_radians=radians)

    if radians:
        deg_scale=np.pi/180
    else:
        deg_scale=1

    U0         = u_infty
    chord      = 2.8480588747143942
    # Parameters
    T           = 10.0 #10 seconds - might be arbitrary
    t_max       = 10.0                 # simulation length
    dt          = T/5000                   # time step
    XLIM        = np.array([-90,90])
    YLIM        = np.array([-3,3])


    # Derived params
    vt       = np.arange(0,t_max,dt)

    # MHH Parameters and Inputs
    p = dynstall_mhh_param_from_polar(P, chord, constants='OpenFAST')

    #Need to add SMD parameters here
    iea15mw_sviv2d = yaml.load(open('chord_3dof.yaml'),Loader=yaml.UnsafeLoader)
    p['x_pitch'] = -iea15mw_sviv2d['displacement'][0]
    m_matrix = np.array(iea15mw_sviv2d['mass_matrix']).reshape(3,3)
    m_inv = linalg.inv(m_matrix)
    k_matrix = np.array(iea15mw_sviv2d['stiffness_matrix']).reshape(3,3)
    c_matrix = np.array(iea15mw_sviv2d['damping_matrix']).reshape(3,3)
    f_matrix = np.array(iea15mw_sviv2d['force_transform_matrix']).reshape(3,3)
    p['m_inv'] = m_inv
    p['m_inv c'] = np.matmul(m_inv, c_matrix)
    p['m_inv k'] = np.matmul(m_inv, k_matrix)
    p['force_transform'] = f_matrix
    p['rho'] = 1.225
    p['prescribed_oscillations'] = True

    y0_mhh = np.r_[ dynstall_mhh_steady_simple(u_infty, aoa_0, p), 0.0, 0.0, aoa_0, 0.0, 0.0, 0.0]
    # Integration using solve_ivp
    sol_mhh = solve_ivp(lambda t,x: dynstall_mhh_dxdt_smd(t,x,u_infty,p), t_span=[0, t_max], y0=y0_mhh, t_eval=vt)

    Cl_mhh = np.zeros_like(vt)
    Cd_mhh = np.zeros_like(vt)
    Cm_mhh = np.zeros_like(vt)
    alpha34_mhh = np.zeros_like(vt)
    for it,t in enumerate(vt):
        sol_mhh.y[9,it] = 2.0 * np.pi * 0.4 * np.radians(10.0) * np.cos(2.0 * np.pi * 0.4 * t)
        sol_mhh.y[6,it] = np.radians(50.0) + np.radians(10.0) * np.sin(2.0 * np.pi * 0.4 * t)
        omega = sol_mhh.y[9,it]
        theta = sol_mhh.y[6,it]
        #alpha34_mhh[it] = np.arctan2( u_infty - omega * (0.75 - p['x_pitch']) * p['chord'] * np.sin(theta), omega * (0.75 - p['x_pitch']) * p['chord'] * np.cos(theta) )
        alpha34_mhh[it] = theta
        Cl_mhh[it],Cd_mhh[it],Cm_mhh[it] = dynstall_mhh_outputs_simple(t,sol_mhh.y[:4,it], u_infty, 0.0, omega, alpha34_mhh[it], p)

    Cl_mhh_orig = np.zeros_like(vt)
    valpha_t =  np.radians(50.0) + np.radians(10.0) * np.sin(2.0 * np.pi * 0.4 * vt)
    valpha_dot_t =  2.0 * np.pi * 0.4 * np.radians(10.0) * np.cos(2.0 * np.pi * 0.4 * vt)

    u=dict()
    u['U']     = lambda t: U0
    u['U_dot'] = lambda t: 0
    u['alpha'] = lambda t: np.interp(t, vt, valpha_t)
    u['omega'] = lambda t: np.interp(t, vt, valpha_dot_t)
    u['alpha_34']  = lambda t: np.interp(t, vt, valpha_t)
    y0_mhh = dynstall_mhh_steady(0,u,p)
    # Integration using solve_ivp
    sol_mhh_orig = solve_ivp(lambda t,x: dynstall_mhh_dxdt(t,x,u,p), t_span=[0, t_max], y0=y0_mhh, t_eval=vt)
    for it,t in enumerate(vt):
        Cl_mhh_orig[it],_,_ = dynstall_mhh_outputs(t,sol_mhh_orig.y[:,it],u,p)


    fig = plt.figure()
    plt.plot(vt, sol_mhh.y[0,:], label='x1')
    plt.plot(vt, sol_mhh_orig.y[0,:], '--', label='x1_orig')
    plt.plot(vt, sol_mhh.y[1,:], label='x2')
    plt.plot(vt, sol_mhh_orig.y[1,:], '--', label='x2_orig')
    plt.plot(vt, sol_mhh.y[2,:], label='x3')
    plt.plot(vt, sol_mhh_orig.y[2,:], '--', label='x3_orig')
    plt.plot(vt, sol_mhh.y[3,:], label='x4')
    plt.plot(vt, sol_mhh_orig.y[3,:], '--', label='x4_orig')
    plt.legend(loc=0)
    plt.tight_layout()

    # fig = plt.figure()
    # plt.plot(vt, sol_mhh.y[6,:], label='Torsion disp - radians')
    # plt.plot(vt, alpha34_mhh, label='alpha34 - radians')
    # plt.plot(vt, Cl_mhh, label='Cl')
    # plt.legend(loc=0)
    # plt.tight_layout()

    # fig = plt.figure()
    # plt.plot(vt, Cl_mhh, label='Cl')
    # plt.plot(vt, Cd_mhh, label='Cd')
    # plt.plot(vt, Cm_mhh, label='Cm')
    # plt.legend(loc=0)
    # plt.tight_layout()

    from welib.tools.colors import fColrs

    fig,ax = plt.subplots(1, 1, sharey=False, figsize=(8.4,5.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    ax.plot(P.alpha/deg_scale  , P.fs  ,     label='f_s')
    ax.plot(P.alpha/deg_scale  , P.cl  , 'k-',label='Cl static', linewidth=2)
    ax.plot(P.alpha/deg_scale  , P.cl_fs  ,'k--',  label='Cl fully separated')
    ax.plot(P.alpha/deg_scale  , P.cl_inv ,'k-.',  label='Cl inviscid')
    col=fColrs(0)
    ax.plot(np.degrees(u['alpha'](vt)), Cl_mhh, ':' , color=col  , label='SMD')
    ax.plot(np.degrees(u['alpha'](vt)), Cl_mhh_orig, '--' ,color=fColrs(1), label='MHH orig')
    ax.tick_params(direction='in')
    ax.set_xlabel('Alpha [deg]')
    ax.set_ylabel('Cl [-]')
    ax.set_xlim(XLIM)
    ax.set_ylim(YLIM)
    ax.legend()
    ax.set_title('FFA-W3-2011 - SMD - MGH dynamic stall model')

def sviv_2d(u_infty):
    radians=True
    #FFA-W3-211 airfoil Dyna Stall
    P=Polar.fromfile(os.path.join(MyDir,'../data/IEA-15-240-RWT_AeroDyn15_Polar_35.dat'),compute_params=True,to_radians=radians)

    if radians:
        deg_scale=np.pi/180
    else:
        deg_scale=1

    U0         = u_infty
    chord      = 2.8480588747143942
    # Parameters
    T           = 30.0 #10 seconds - might be arbitrary
    t_max       = 30.0                 # simulation length
    dt          = T/5000                   # time step
    #XLIM        = np.array([0,40])
    XLIM        = np.array([-90,90])
    YLIM        = np.array([-3,3])


    # Derived params
    vt       = np.arange(0,t_max,dt)

    # MHH Parameters and Inputs
    p = dynstall_mhh_param_from_polar(P, chord, constants='OpenFAST')

    #Need to add SMD parameters here

    iea15mw_sviv2d = yaml.load(open('chord_3dof.yaml'),Loader=yaml.UnsafeLoader)
    p['x_pitch'] = -iea15mw_sviv2d['displacement'][0]
    m_matrix = np.array(iea15mw_sviv2d['mass_matrix']).reshape(3,3)
    m_inv = linalg.inv(m_matrix)
    k_matrix = np.array(iea15mw_sviv2d['stiffness_matrix']).reshape(3,3)
    c_matrix = np.array(iea15mw_sviv2d['damping_matrix']).reshape(3,3)
    f_matrix = np.array(iea15mw_sviv2d['force_transform_matrix']).reshape(3,3)
    p['m_inv'] = m_inv
    p['m_inv c'] = np.matmul(m_inv, c_matrix)
    p['m_inv k'] = np.matmul(m_inv, k_matrix)
    p['force_transform'] = f_matrix
    p['rho'] = 1.225
    p['reference_aoa'] = np.radians(iea15mw_sviv2d['angle'])
    p['prescribed_oscillations'] = False

    y0_mhh = np.r_[ dynstall_mhh_steady_simple(u_infty, np.radians(50.0), p), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    # Integration using solve_ivp
    sol_mhh = solve_ivp(lambda t,x: dynstall_mhh_dxdt_smd(t,x,u_infty,p), t_span=[0, t_max], y0=y0_mhh, t_eval=vt)

    Cl_mhh = np.zeros_like(vt)
    Cd_mhh = np.zeros_like(vt)
    Cm_mhh = np.zeros_like(vt)
    alpha34_mhh = np.zeros_like(vt)
    for it,t in enumerate(vt):
        omega = sol_mhh.y[9,it]
        theta = -sol_mhh.y[6,it] + p['reference_aoa']
        #alpha34_mhh[it] = np.arctan2( u_infty - omega * (0.75 - p['x_pitch']) * p['chord'] * np.sin(theta), omega * (0.75 - p['x_pitch']) * p['chord'] * np.cos(theta) )
        alpha34_mhh[it] = theta
        Cl_mhh[it],Cd_mhh[it],Cm_mhh[it] = dynstall_mhh_outputs_simple(t,sol_mhh.y[:4,it], u_infty, 0.0, omega, alpha34_mhh[it], p)

    #Could return some combination of sol_mhh and Cl_mhh, Cd_mhh, Cm_mhh here for analysis
    print('alpha34_mhh:')
    print(alpha34_mhh)

    fig = plt.figure()
    plt.plot(vt, sol_mhh.y[0,:], label='x1')
    plt.plot(vt, sol_mhh.y[1,:], label='x2')
    plt.plot(vt, sol_mhh.y[2,:], label='x3')
    plt.plot(vt, sol_mhh.y[3,:], label='x4')
    plt.legend(loc=0)
    plt.tight_layout()

    fig = plt.figure()
    plt.plot(vt, sol_mhh.y[4,:], label='Flap disp')
    plt.plot(vt, sol_mhh.y[5,:], label='Edge disp')
    plt.legend(loc=0)
    plt.tight_layout()

    fig = plt.figure()
    plt.plot(vt, -sol_mhh.y[6,:] + p['reference_aoa'], label='Torsion disp - radians')
    plt.plot(vt, alpha34_mhh, label='alpha34 - radians')
    plt.plot(vt, Cl_mhh, label='Cl')
    plt.legend(loc=0)
    plt.tight_layout()

    fig = plt.figure()
    plt.plot(vt, Cl_mhh, label='Cl')
    plt.plot(vt, Cd_mhh, label='Cd')
    plt.plot(vt, Cm_mhh, label='Cm')
    plt.legend(loc=0)
    plt.tight_layout()

    from welib.tools.colors import fColrs

    fig,ax = plt.subplots(1, 1, sharey=False, figsize=(8.4,5.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    ax.plot(P.alpha/deg_scale  , P.fs  ,     label='f_s')
    ax.plot(P.alpha/deg_scale  , P.cl  , 'k-',label='Cl static', linewidth=2)
    ax.plot(P.alpha/deg_scale  , P.cl_fs  ,'k--',  label='Cl fully separated')
    ax.plot(P.alpha/deg_scale  , P.cl_inv ,'k-.',  label='Cl inviscid')
    col=fColrs(0)
    ax.plot(np.degrees(-sol_mhh.y[6,:] + p['reference_aoa']), Cl_mhh, ':' , color=col  , label='SMD')
    ax.tick_params(direction='in')
    ax.set_xlabel('Alpha [deg]')
    ax.set_ylabel('Cl [-]')
    ax.set_xlim(XLIM)
    ax.set_ylim(YLIM)
    ax.legend()
    ax.set_title('FFA-W3-2011 - SMD - MGH dynamic stall model')


if __name__ == '__main__':
    sviv_2d_prescribed_oscillations(15.0, np.radians(50.0))
    # sviv_2d(15.0)
    # prescribed_oscillations()
    plt.show()
