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
    a1    = 0.1 # axial induction used before the step
    a2    = 0.2 # axial induction used after the step
    tstep = 10  # time at which the step of axial induction occurs
    tmax  = 50  # simulation time

    # Derived parameters
    p=dict()
    p['tau1'] = tau1_oye(a1, R, U0)
    p['tau2'] = tau2_oye(r_bar ,p['tau1'])
    p['k']    = 0.6

    # Define functions of time (NOTE: these functions would be set within an unsteady BEM algorithm)
    # Look at the function dyninflow_oye_sim for lower level interfaces
    u=dict()
    u['Vqs']  = lambda t: a1*U0 if t<tstep else a2*U0  
    time = np.linspace(0,tmax,1000)

    df_d = dyninflow_oye_sim(time, u, p, x0=None, prefix='', method='discrete')
    df_c = dyninflow_oye_sim(time, u, p, x0=None, prefix='', method='continuous')
    print(df_d)
    # Plot
    fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    ax.plot(time, df_d['Vdyn_[m/s]']/U0, '-' , label=r'$a_{dyn}$ (discrete formulation)')
    ax.plot(time, df_c['Vdyn_[m/s]']/U0, ':' , label=r'$a_{dyn}$ (continuous formulation)')
    ax.plot(time, df_c['Vqs_[m/s]']/U0 , 'k-', label=r'$a_{qs}$')
    ax.legend()
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Axial induction [-]')
    ax.set_ylim([0.09,0.21])
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
