"""
References:
   [1] Branlard 2017 Wind Turbine Aerodynamics and Vorticity-based method Chap10, p198
   [2] Branlard et al. 2022: Dynamic inflow and unsteady aerodynamics models for modal and stability analyses in OpenFAST

"""
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

def tau1_oye(a_bar, R, U0):
    """ """
    return  1.1/( 1-1.3*min(a_bar,0.5) ) * R/U0

def tau1_dbemt(a_bar, R, Un0_disk):
    # We need to extrapolate the radius and disk velocity to the i+1 timestep
    # We will grab the i+1 version of vind,s and the disk averaged induction by using the
    # the already updated states of the BEMT module.
    #bjj: I believe u(1) is at t, which seems inconsistant with this comment
    U0n_disk = min(U0n_disk,0.1)
    return min( 1.1/( 1-1.3*min(a_disk,0.5) ) * R/Un0_disk , 100) # plateau 100s

def tau2_oye(r_bar, tau1):
    return (0.39-0.26*r_bar**2)*tau1

def dynflow_oye_dxdt(t, x, u, tau1, tau2, k=0.6):
    """ 
    Continuous version of Oye's dynamic inflow model. Two states:
        x[0]: Vred:      Vint = Vred + k*Vqs
        x[1]: Vdyn
    See Reference [2]
    """
    x_qs     = u['Vqs'](t)
    xred = x[0]
    xdyn = x[1]
    xred_dot =               - 1/tau1*xred + (1-k)/tau1*x_qs
    xdyn_dot  = -1/tau2*xdyn + 1/tau2*xred +     k/tau2*x_qs
    return [xred_dot, xdyn_dot]

def dynflow_oye_steady(t, u, p):
    Vqs  = u['Vqs'](t)
    x    = np.zeros(2)
    x[0] = Vqs * (1-p['k'])   # Vred   # Vint=Vred + kVqs
    x[1] = Vqs                # Vdyn
    return x

def dyninflow_oye_update_discr(t, dt, xd_old, u, tau1, tau2, k=0.6):
    """
    Discrete formualtion of Oye's dynamics inflow model
    See Reference [1]

    xd_old: old discrete states
         0: dynamic value
         1: intermediate value
         2: quasi steady value
    u: 
         0: new quasi steady value
    """
    # Oye's dynamic inflow model, discrete time integration
    Vdyn_old = xd_old[0]
    Vint_old = xd_old[1]
    Vqs__old = xd_old[2]
    Vqs__new = u['Vqs'](t)
    H        = Vqs__new + k * tau1 * (Vqs__new - Vqs__old) /dt
    Vint_new = H + (Vint_old - H) * np.exp(-dt/tau1) # intermediate velocity
    Vdyn_new = Vint_new + (Vdyn_old - Vint_new) * np.exp(-dt/tau2)
    # 
    # 
    xd_new = np.zeros(3)
    xd_new[0] = Vdyn_new
    xd_new[1] = Vint_new
    xd_new[2] = Vqs__new
    return xd_new

def dyninflow_oye_discr_steady(t,u):
    """ Return steady state values for Oye's discrete dynamic flow model """
    Vqs = u['Vqs'](t)
    xd = np.zeros(3)
    xd[0] = Vqs
    xd[1] = Vqs
    xd[2] = Vqs
    return xd


def dyninflow_oye_sim(time, u, p, x0=None, prefix='', method='discrete'):
    """ Perform simulation using Oye dynamic inflow model

    INPUTS:
     - time: time vector
     - u: dictionary of input functions of time
         keys 'Vqs' :quasi steady value
     - p: dictionary of parameters
     - x0: initial conditions for the model. 
           If None, the steady steady values are used
     - method: 'continuous' or 'discrete' to chose a formulation
     - prefix: prefix used for channel names of the dataframe. Use 'AB1N001' to match OpenFAST.
    OUTPUTS:
     - df: dataframe with outputs similar to UA module of OpenFAST
    """

    # --- Parameters. For now use constant taus
    tau1 = p['tau1']
    tau2 = p['tau2']
    k    = p['k']

    # --- Prepare outputs
    df=pd.DataFrame()
    df['Time_[s]'] = time

    if method=='discrete':
        # --- Initial conditions for states
        if x0 is None:
            x0 = dyninflow_oye_discr_steady(0,u)
        # --- Time Integration of states
        y  = np.zeros((3,len(time)))
        xd = x0
        y[:,0] = xd
        for it,t in enumerate(time[1:]):
            dt = t - time[it] # Note: time[it] is in fact t-dt
            xd = dyninflow_oye_update_discr(t, dt, xd, u, tau1=tau1, tau2=tau2, k=k)
            y[:,it+1] = xd
        # --- Calc output
        df[prefix + 'Vdyn_[m/s]'] = y[0,:]
        df[prefix + 'Vint_[m/s]'] = y[1,:]
        df[prefix + 'Vqs_[m/s]']  = y[2,:]
    elif method=='continuous':
        # --- Initial conditions for states
        if x0 is None:
            x0 = dynflow_oye_steady(0, u, p)
        # --- Time integration of states
        system = lambda t,x: dynflow_oye_dxdt(t, x, u, tau1=tau1, tau2=tau2, k=k)
        sol = solve_ivp(system , t_span=[0, max(time)], y0=x0, t_eval=time, method='LSODA')
        y=sol.y
        # --- Calc output
        df[prefix + 'Vred_[m/s]']   = y[0,:]
        df[prefix + 'Vdyn_[m/s]'] = y[1,:]
        df[prefix + 'Vqs_[m/s]']  = [u['Vqs'](t) for t in time]
        df[prefix + 'Vint_[m/s]'] = df[prefix + 'Vred_[m/s]'] + k*df[prefix + 'Vqs_[m/s]']
    elif method=='continuous_v2':
        # --- Initial conditions for states
        if x0 is None:
            x0 = dynflow_oye_v2_steady(0, u)
        # --- Time integration of states
        system = lambda t,x: dynflow_oye_v2_dxdt(t, x, u, tau1=tau1, tau2=tau2, k=k)
        sol = solve_ivp(system , t_span=[0, max(time)], y0=x0, t_eval=time)
        y=sol.y
        # --- Calc output
        df[prefix + 'Vdyn_[m/s]']     = y[0,:]
        df[prefix + 'dVdyn/dt_[m/s]'] = y[1,:]
        df[prefix + 'Vqs_[m/s]']  = [u['Vqs'](t) for t in time]
        df[prefix + 'dVqs_dt_[m/s]']  = [u['dVqs_dt'](t) for t in time]
    elif method=='continuous_v4':
        # --- Initial conditions for states
        if x0 is None:
            x0 = dynflow_oye_v4_steady(0, u, p)
        # --- Time integration of states
        system = lambda t,x: dynflow_oye_v4_dxdt(t, x, u, tau1=tau1, tau2=tau2, k=k)
        sol = solve_ivp(system , t_span=[0, max(time)], y0=x0, t_eval=time, method='LSODA')
        y=sol.y
        # --- Calc output
        df[prefix + 'Vdyn_[m/s]'] = y[0,:]
        df[prefix + 'Vint_[m/s]'] = y[1,:]
        df[prefix + 'Vqs_[m/s]']  = [u['Vqs'](t) for t in time]
    else:
        raise NotImplementedError(method)


    return df


# class DBEMT():
#     class Dummy():
#         pass
#     def __init__(self,dt,r,R,mode,tau1_cst=None)
#         self.p=Dummy()
#         self.x=Dummy()
#         self.o=Dummy()
#         self.m=Dummy()
#         self.p.dt                    = dt
#         self.p.numNodes, self.p.numBlades = r.shape
#         self.p.k_0ye                 = 0.6
#         self.p.tau1_const            = tau1_const
#         self.p.mode = mode
#         self.p.R = R
#         
#         if self.p.mode == 'tauConst':
#             self.p.spanRatio=r/R
#             if tau1_const is None:
#                 raise Exception('Please provide tau1_cst when `mode` is constant')
#         self.x.vind   = np.zeros((2,p.numNodes, p.numBlades))# This is the axial and tangential induced velocity at node i on blade j
#         self.x.vind_1 = np.zeros((2,p.numNodes, p.numBlades))# This is the axial and tangential induced velocity at node i on blade j
#         self.o.areStatesInitialized = False
#         self.o.tau1 = 0.0
#         self.m.FirstWarn_tau1 = True
# 
# 
#     def UpdateStates(self,u1,u2):
#         if self.mode=='tauConst':
#              self.o.tau1 = self.p.tau1_const
#         else:
#              self.o.tau1 = tau1_dbemt(u1.a_disk, self.p.R, u1.Un0_disk)
#         #    if ( .not. OtherState%areStatesInitialized(i,j) ) then
#         #       x%vind_1(1,i,j) = u(1)%vind_s(1)
#         #       x%vind_1(2,i,j) = u(1)%vind_s(2)
#         #       x%vind(1,i,j) = u(1)%vind_s(1)
#         #       x%vind(2,i,j) = u(1)%vind_s(2)
#         #       OtherState%areStatesInitialized(i,j) = .true.
#         #       return
#         #    end if
#         #    if ( p%DBEMT_Mod == DBEMT_tauConst ) then
#         #       spanRatio   = p%spanRatio(i,j)
#         #    else
#         #       spanRatio = u(1)%spanRatio
#         #    end if
#         #    do indx=1,2 !Axial and tangential components.  jmj questions if this should be done for tangential induction.
#         #       ! TODO: Deal with initialization so that we avoid spikes???
#         #       B =  ( u(2)%vind_s(indx) - u(1)%vind_s(indx) ) / p%dt                                              ! Eqn. 1.17c  ! bjj: note that interpOrder will affect this numerical derivative
#         #       A = u(1)%vind_s(indx) + B*p%k_0ye*OtherState%tau1                                                  ! Eqn. 1.17b
#         #       
#         #       C0 = x%vind_1(indx,i,j) - A + B*OtherState%tau1                                                    ! Eqn. 1.21b
#         #       
#         #       x%vind_1(indx,i,j) = C0*exp(-p%dt/OtherState%tau1) + A + B*(p%dt-OtherState%tau1)                  ! Eqn. 1.21a, but this is using p%dt instead of t
#         #       
#         #       k_tau = 0.39 - 0.26*spanRatio**2                                                                   ! Eqn. 1.23b
#         #       tau2 = k_tau*OtherState%tau1                                                                       ! Eqn. 1.7 or Eqn 1.23a
#         #       C0_2 = x%vind(indx,i,j) - C0/(1-k_tau) - A + B*(OtherState%tau1 + tau2)                            ! Eqn. 1.24c 
#         #       x%vind(indx,i,j) = C0_2*exp(-p%dt/tau2) + A + B*(p%dt-OtherState%tau1-tau2) + (C0/(1-k_tau))*exp(-p%dt/OtherState%tau1)  ! Eqn. 1.25
# 
#     def CalcOutput(self):
#         pass
#         #if ( .not. OtherState%areStatesInitialized(i,j) ) :
#         #   y_vind = u%vind_s
#         #else:
#         #   y_vind(:) = x%vind(:,i,j)
# 


def run_astep(a1, a2, r_bar, tau1=None, dt=0.1, tstep=10, tmax=100, U0=10, R=50, k=0.6):
    """
    Run step of axial induction, assuming "U0" constant, and tau1 and tau2 constant
    """
    from welib.tools.colors import fColrs
    p=dict()
    if tau1 is None:
        tau1 = tau1_oye(a1, R, U0)
    p['tau1'] = tau1
    p['tau2'] = tau2 = tau2_oye(r_bar,tau1)
    p['k']    = k

    time = np.arange(0,tmax+dt/2, dt)
    # Inputs for dyna stall 
    u=dict()
    #u['Vqs']     = lambda t: a1*U0 if t<tstep else a2*U0
    u['Vqs']     = lambda t: a1 if t<tstep else a2
    u['dVqs_dt'] = lambda t: 0


    df_d = dyninflow_oye_sim(time, u, p, x0=None, prefix='', method='discrete')
    df_c = dyninflow_oye_sim(time, u, p, x0=None, prefix='', method='continuous')
    p['k']=0.9
    df_c9= dyninflow_oye_sim(time, u, p, x0=None, prefix='', method='continuous')
    p['k']=0.1
    df_c1= dyninflow_oye_sim(time, u, p, x0=None, prefix='', method='continuous')

    fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
#     ax.plot(df_d['Time_[s]'], df_d['Vdyn_[m/s]'], '-' , label='dynamic discrete')
#     ax.plot(df_d['Time_[s]'], df_d['Vint_[m/s]'], '-' , label='intermediate discrete')
    ax.plot(df_c ['Time_[s]'], df_c ['Vdyn_[m/s]'], '-' , c=fColrs(1), label='dynamic continuous')
    ax.plot(df_c1['Time_[s]'], df_c1['Vdyn_[m/s]'], '-' , c=fColrs(2), label='dynamic continuous 1')
    ax.plot(df_c9['Time_[s]'], df_c9['Vdyn_[m/s]'], '-' , c=fColrs(3), label='dynamic continuous 9')

    ax.plot(df_c ['Time_[s]'], df_c ['Vint_[m/s]'], ':' , c=fColrs(1), label='intermediate continuous')
    ax.plot(df_c1['Time_[s]'], df_c1['Vint_[m/s]'], ':' , c=fColrs(2), label='intermediate continuous 1')
    ax.plot(df_c9['Time_[s]'], df_c9['Vint_[m/s]'], ':' , c=fColrs(3), label='intermediate continuous 9')

    ax.plot(df_c ['Time_[s]'], df_c ['Vred_[m/s]'],   '--', c=fColrs(1), label='reduced continuous')
    ax.plot(df_c1['Time_[s]'], df_c1['Vred_[m/s]'],   '--', c=fColrs(2), label='reduced continuous 1 ')
    ax.plot(df_c9['Time_[s]'], df_c9['Vred_[m/s]'],   '--', c=fColrs(3), label='reduced continuous 9')

    ax.plot(df_d['Time_[s]'], df_d['Vqs_[m/s]'] , 'k-', label='quasi-steady')
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.legend()
    ax.legend()
    ax.set_title('DynInflow step')
    plt.show()

#         ax.plot((vt-tstep)/scale,[a_qs(t) for t in vt],'k', label='Quasi-steady')
#         ax.plot((vt-tstep)/scale,sol.y[0,:], label='tau1={:4.1f} - U_0={:2.0f}'.format(tau1,U0))









def example_astep_WS():
    a_mean = 0.4
    r_bar=0.5
    tstep=10

    U0     = 1
    vt = np.linspace(0,100,10000)
    vU0 = np.arange(5,25,5)
    vR      = [65]
    vR      = [50, 100]
    R      = 50

    fig,ax=plt.subplots(1,1)
    for i,U0 in enumerate(vU0):
        # Inputs for dyna stall 
        fU0        = lambda t: U0
        a_qs      = lambda t: a_mean if t<tstep else a_mean/2
        da_qs_dt  = lambda t: 0
        a_bar     = lambda t: a_mean if t<tstep else a_mean/2
        dtau1_dt  = lambda t: 0
        # Initial value
        y0=[a_qs(0),da_qs_dt(0)] # a, adot
        y0=[a_mean,0] # a, adot
        #y0=[0,0] # a, adot

        tau1 = tau1_oye(a_mean, R, U0)
        tau2 = tau2_oye(r_bar,tau1)
        print('tau1',tau1, 'tau2',tau2,'R/U0',R/U0)

        # regular system
        system = lambda t,x: dynflow_oye_dxdt(t,x, a_qs(t) , a_bar(t), da_qs_dt(t), dtau1_dt(t), fU0(t), r_bar*R, R)
        sol = solve_ivp(system , t_span=[0, max(vt)], y0=y0, t_eval=vt)

        scale=1
        if i==0:
            ax.plot((vt-tstep)/scale,[a_qs(t) for t in vt],'k', label='Quasi-steady')
        ax.plot((vt-tstep)/scale,sol.y[0,:], label='tau1={:4.1f} - U_0={:2.0f}'.format(tau1,U0))
    ax.legend()
    ax.set_title('DynInflow Influence of WindSpeed')



def example_astep_radius():
    a_mean = 0.5
    tstep  = 20
    U0     = 8
    vt     = np.linspace(0,100,10000)
    R      = 50
    vr_bar = np.arange(0.1,1,0.2)

    fig,ax=plt.subplots(1,1)
    for i,r_bar in enumerate(vr_bar):
        # Inputs for dyna stall 
        fU0        = lambda t: U0
        a_qs      = lambda t: a_mean if t<tstep else a_mean/2
        da_qs_dt  = lambda t: 0
        a_bar     = lambda t: a_mean if t<tstep else a_mean/2
        dtau1_dt  = lambda t: 0
        # Initial value
        y0=[a_qs(0),da_qs_dt(0)] # a, adot
        y0=[a_mean,0] # a, adot
        #y0=[0,0] # a, adot

        tau1 = tau1_oye(a_mean, R, U0)
        tau2 = tau2_oye(r_bar,tau1)
        print('tau1',tau1, 'tau2',tau2,'R/U0',R/U0)

        # regular system
        system = lambda t,x: dynflow_oye_dxdt(t,x, a_qs(t) , a_bar(t), da_qs_dt(t), dtau1_dt(t), fU0(t), r_bar*R, R)
        sol = solve_ivp(system , t_span=[0, max(vt)], y0=y0, t_eval=vt)

        scale=tau1
        if i==0:
            ax.plot((vt-tstep)/scale,[a_qs(t) for t in vt],'k', label='Quasi-steady')
        ax.plot((vt-tstep)/scale,sol.y[0,:], label='r_bar={:3.1f}'.format(r_bar))
    ax.set_xlabel('t/tau1 [-]')
    ax.set_xlim([-1,4])
    ax.legend()
    ax.set_title('DynInflow Influence of Radius')

def dyninflow_case(Case='Sine', T=5, U0=10, R=65,r_bar=0.5,a_mean=0.2,a_ampl=0.15,U0_ampl=3):
    # Inputs for dyna stall 
    if Case=='Sine':
        Case='SineT{:.1f}'.format(T)
        omega = 2*np.pi/T
        a_qs      = lambda t: a_mean+a_ampl*np.sin(omega*t)
        da_qs_dt  = lambda t: a_ampl*omega*np.cos(omega*t)
        a_bar     = lambda t: a_mean
        dtau1_dt  = lambda t: 0
        fU0       = lambda t: U0
        dU0_dt    = lambda t: 0
        # Initial value
        y0=[a_mean/2,0] # a, adot
        tmax=20*T
    elif Case=='WindSine': # Both wind and a have sine
        Case='WindSineT{:.1f}'.format(T)
        omega = 2*np.pi/T
        a_qs      = lambda t: a_mean+a_ampl*np.sin(omega*t)
        da_qs_dt  = lambda t: a_ampl*omega*np.cos(omega*t)
        a_bar     = lambda t: a_mean
        dtau1_dt  = lambda t: 0 # TODO TODO
        fU0       = lambda t: U0- U0_ampl*np.sin(omega*t)
        dU0_dt    = lambda t: -U0_ampl*omega*np.cos(omega*t)
        y0=[a_mean/2,0] # a, adot
        tmax=20*T
    elif Case=='StepDown':
        a_qs      = lambda t: a_mean if t<10 else a_mean/2
        da_qs_dt  = lambda t: 0
        a_bar     = lambda t: a_mean if t<10 else a_mean/2
        dtau1_dt  = lambda t: 0
        fU0       = lambda t: U0
        dU0_dt    = lambda t: 0
        # Initial value
        y0=[a_mean,0] # a, adot
        tmax=50
    elif Case=='StepUp':
        a_qs      = lambda t: a_mean/2 if t<10 else a_mean
        da_qs_dt  = lambda t: 0
        a_bar     = lambda t: a_mean/2 if t<10 else a_mean
        dtau1_dt  = lambda t: 0
        fU0       = lambda t: U0
        dU0_dt    = lambda t: 0
        # Initial value
        y0=[a_mean/2,0] # a, adot
        tmax=50
    vt = np.linspace(0,tmax,1000)
    return Case, a_qs, da_qs_dt, a_bar, dtau1_dt, fU0, dU0_dt, vt, y0

def example_algo_investigation(Case='Sine', T=5, U0=10, R=65,r_bar=0.5,a_mean=0.2,a_ampl=0.15,U0_ampl=3):

    Case, a_qs, da_qs_dt, a_bar, dtau1_dt, fU0, dU0_dt, vt, y0 = dyninflow_case(Case, T, U0, R,r_bar,a_mean,a_ampl,U0_ampl)

    tau1 = tau1_oye(a_mean, R, U0)
    tau2 = tau2_oye(r_bar ,tau1)
    print('tau1',tau1, 'tau2',tau2)

    # regular system
    system = lambda t,x: dynflow_oye_dxdt(t, x, a_qs(t) , a_bar(t), da_qs_dt(t), dtau1_dt(t), fU0(t), r_bar*R, R)
    sol = solve_ivp(system , t_span=[0, max(vt)], y0=y0, t_eval=vt)

    # No derivative
    system = lambda t,x: dynflow_oye_dxdt(t, x, a_qs(t) , a_bar(t), 0          , dtau1_dt(t), fU0(t), r_bar*R, R)
    sol_no_deriv = solve_ivp(system , t_span=[0, max(vt)], y0=y0, t_eval=vt)

    # Implicit a_bar
    system = lambda t,x: dynflow_oye_dxdt(t, x, a_qs(t) , x[0]    , da_qs_dt(t), dtau1_dt(t), fU0(t), r_bar*R, R)
    sol_implicit = solve_ivp(system , t_span=[0, max(vt)], y0=y0, t_eval=vt)

    # Implicit a_bar
    system = lambda t,x: dynflow_oye_dxdt(t, x, a_qs(t) , x[0]    ,     0      , dtau1_dt(t), fU0(t), r_bar*R, R)
    sol_implicit_no_deriv = solve_ivp(system , t_span=[0, max(vt)], y0=y0, t_eval=vt)



    fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)

    ax.plot(vt,sol.y[0,:], label='a')
    ax.plot(vt,sol_no_deriv.y[0,:], '-.',label='a (no deriv)')
    ax.plot(vt,sol_implicit.y[0,:], '--',label='a (implicit)')
    ax.plot(vt,sol_implicit_no_deriv.y[0,:], ':',label='a (implicit no deriv)')
    ax.plot(vt,[a_qs(t) for t in vt], label='a_qs')
#     plt.plot(vt,sol.y[1,:], label='da/dt')
#     plt.plot(vt,da_qs_dt(vt), label='da_qs/dt')
    ax.legend()
    ax.set_title('DynInflow {} AlgoInvest'.format(Case))
    return ax


def example_tau1_abar():
    va_bar = np.linspace(0,0.6,100)
    tau1=[tau1_oye(a_bar, R=1, U0=1) for a_bar in va_bar]

    plt.figure()
    plt.plot(va_bar,tau1)
    plt.xlabel('a_bar [-]')
    plt.ylabel('tau1 /(R/U0) [-]')
    plt.title('DynInflow Tau1 vs a_bar')

def example_tau2_radius():

    vr_bar = np.linspace(0,1,100)
    tau2=tau2_oye(vr_bar,tau1=1)

    plt.figure()
    plt.plot(vr_bar,tau2)
    plt.xlabel('r/R [-]')
    plt.ylabel('tau2 / tau1 [-]')
    plt.title('DynInflow Tau2 vs Tau1')


if __name__=='__main__':
    from welib.tools.clean_exceptions import *
    run_astep(a1=0.4, a2=0.5, r_bar=1, tau1=None, dt=0.1, tstep=10, tmax=100, U0=10, R=50, k=0.6)

    example_astep_WS()
#     example_astep_radius()
#     example_tau2_radius()
#     example_tau1_abar()
#     example_algo_investigation(Case='StepUp'  )
#     example_algo_investigation(Case='StepDown')
#     example_algo_investigation(Case='Sine',T=20)
#     example_algo_investigation(Case='Sine',T=5)
#     example_algo_investigation(Case='WindSine',T=5)

#     export2png()



    plt.show()
