""" 
Time integration of the Lorenz system

solve_ivp and ode_int tend to give different results after a while

"""

import numpy as np
from scipy.integrate import odeint, solve_ivp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from welib.system.chaos import dqdt_lorenz
 
 
# --- Parameters
sigma  = 10.0
beta   = 8.0 / 3.0
rho    = 28.0
tmax   = 40
dt     = 0.01
method = 'LSODA'         # RK23, RK45, BDF, LSODA, Radau, DOP853
y0     = [1.0, 1.0, 1.0] # Initial state of the system

# --- Derived parameters
p = (sigma, beta, rho)  # Parameters of the system
t_eval = np.arange(0.0, tmax, dt)
 
# --- Time integration
res1 = odeint(dqdt_lorenz, y0, t_eval, p, tfirst=True) # NOTE: uses LSODA method
res2 = solve_ivp(dqdt_lorenz, (t_eval[0],t_eval[-1]), y0, args=p, method=method, t_eval=t_eval, atol=1.49e-8, rtol=1.49e-8)
 
# --- Plots
fig = plt.figure()
ax = fig.add_subplot(1, 2, 1, projection='3d')
ax.plot(res1[:, 0], res1[:, 1], res1[:, 2])
ax.set_title('odeint')
 
ax = fig.add_subplot(1, 2, 2, projection='3d')
ax.plot(res2.y[0, :], res2.y[1, :], res2.y[2, :])
ax.set_title('solve_ivp')

fig.suptitle('System - Lorenz attractor')



if __name__ == '__main__':

    fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    ax.plot(t_eval, res1[:,0]    , label='')
    ax.plot(t_eval, res2.y[0,:]    , label='')
    ax.set_xlabel('Time ')
    ax.set_ylabel('x')

    plt.show()

if __name__=="__export__":
    from welib.tools.repo import export_figs_callback
    export_figs_callback(__file__)

