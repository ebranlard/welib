"""
Example of dynamic inflow model, used in isolation for a given radial station 

This example assumes a constant wind speed, and a step of quasi steady (qs) axial induction

"""
import numpy as np
import matplotlib.pyplot as plt
# Local 
from welib.dyninflow.DynamicInflow import *


def main():
    # Parameters
    U0     = 10 # mean wind speed [m/s]
    R      = 65  # rotor radius [m]
    r_bar  = 0.5 # spanwise station investigated [-]
    a_mean = 0.2 # axial induction used for the step

    # Define functions of time (NOTE: these functions would be set within an unsteady BEM algorithm)
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

    tau1 = tau1_oye(a_mean,R/U0)
    tau2 = tau2_oye(r_bar ,tau1)

    # integrate system
    system = lambda t,x: dynflow_oye_dxdt(t,x, a_qs(t) , a_bar(t), da_qs_dt(t), dtau1_dt(t), fU0(t), r_bar*R, R)
    sol = solve_ivp(system , t_span=[0, max(vt)], y0=y0, t_eval=vt)

    # Plot
    fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    ax.plot(vt,sol.y[0,:], label='a')
    ax.plot(vt,[a_qs(t) for t in vt], label='a_qs')
    ax.legend()
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Axial induction [-]')
    ax.set_title('Dynamic Inflow (Oye) - induction step')
    return ax

main()

if __name__ == '__main__':
    plt.show()
if __name__=="__test__":
    pass
if __name__=="__export__":
    from welib.tools.repo import export_figs_callback
    export_figs_callback(__file__)
