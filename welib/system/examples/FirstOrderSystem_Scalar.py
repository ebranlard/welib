"""
Example showing the responses of a first order system defined by:

   x_dot(t) + 1/tau x = f(t) 

"""
import numpy as np
import matplotlib.pyplot as plt
from welib.system.firstorder import *
from welib.tools.signal_analysis import impulse, step, ramp, hat, correlated_signal
from welib.tools.tictoc import Timer

# --- Script parameters
simpleResponses = True
randomSignal    = False
tau = 0.1       # time constant of first order system
time = np.linspace(-3*tau,10*tau,1001) # time span
x0  = 0         # Initial condition, x[time[0]]
# --- Parameters for simple inputs
t0  = 2*tau     # initial time where input starts (for simple inputs)
A   = 1         # amplitude for all simple inputs
omega = 30      # for sine 
phi   = np.pi/3 # for sine
T     = 2*tau   # for hat

# --- Derived parameters
t_norm = (time-t0)/tau

if simpleResponses:
    # --- Analytical integrations
    x_d = impulse_response(time, tau, t0=t0, A=A, x0=x0)
    x_s = step_response   (time, tau, t0=t0, A=A, x0=x0)
    x_r = ramp_response   (time, tau, t0=t0, A=A, x0=x0)
    x_o = sine_response   (time, tau, t0=t0, A=A, x0=x0, omega=omega, phi=phi)
    x_h = hat_response    (time, tau, t0=t0, A=A, x0=x0, T=T)

    # --- Numerical integrations of simple functions
    u_d = impulse(time, t0=t0, A=A, epsilon=0.08) # NOTE: smooth, and sensitive
    t_d0 =(np.array([np.min(time), t0,t0,t0, np.max(time)])-t0)/tau # fake dirac for plotting
    u_d0 = np.array([0, np.min(u_d),np.max(u_d)*1.1,np.min(u_d), 0])
    u_s = step   (time, t0=t0, A=A)
    u_r = ramp   (time, t0=t0, A=A)
    u_h = hat    (time, t0=t0, A=A, T=T, method='sum')
    u_o = A*np.sin(omega*(time-t0)+phi); u_o[time<t0]=0

    # ----  if method in ['LSODA','RK45','RK23','DOP853','Radau','BDF']:
    #x_D = integrate(time, tau, u_d, x0=x0, method='RK45') # NOTE sensitive
    x_D = integrate(time, tau, u_d, x0=x0, method='DOP853') # NOTE sensitive
    x_S = integrate(time, tau, u_s, x0=x0, method='LSODA')
    x_R = integrate(time, tau, u_r, x0=x0, method='LSODA')
    #x_H = integrate(time, tau, u_h, x0=x0, method='DOP853') # NOTE sensitive, Funny behavior for T=1
    x_H = integrate(time, tau, u_h, x0=x0, method='RK45') # NOTE sensitive, Funny behavior to T=2
    x_O = integrate(time, tau, u_o, x0=x0, method='LSODA')

    # --- Numerical integration using Duhamel/convolution
    x_D2 = integrate(time, tau, u_d, x0=x0, method='convolution')
    x_S2 = integrate(time, tau, u_s, x0=x0, method='convolution')
    x_R2 = integrate(time, tau, u_r, x0=x0, method='convolution')
    x_H2 = integrate(time, tau, u_h, x0=x0, method='convolution')
    x_O2 = integrate(time, tau, u_o, x0=x0, method='convolution')

if randomSignal:
    # --- Numerical integrations of correlated, random signal
    u_m = correlated_signal(coeff=0.9, n=len(time), seed=129) # seed provided for reproducibility
    #with Timer('Num'):
    #x_M  = integrate(time, tau, u_m, x0=x0, method='DOP853')
    x_M  = integrate(time, tau, u_m, x0=x0, method='RK23')
    #with Timer('Conv'):
    x_M2 = integrate(time, tau, u_m, x0=x0, method='convolution')


# --------------------------------------------------------------------------------}
# --- Individual plots 
# --------------------------------------------------------------------------------{
def plot_input_output(t_norm, u, x_th=None, x_num=None, x_conv=None, fig=None, axes=None):
    if fig is None:
        fig,axes = plt.subplots(2, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
        fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    axes=np.atleast_1d(axes)
#     axes[0].plot(t_d0  , u_d0,'-'  , label=r'')
    axes[0].plot(t_norm , u     , 'k-' , label=r'')
    if x_th is not None:
        axes[1].plot(t_norm , x_th  , 'k-'  , label=r'Analytical')
    if x_num is not None:
        axes[1].plot(t_norm , x_num , '--' , label=r'Numerical')
    if x_conv is not None:
        axes[1].plot(t_norm , x_conv, ':'  , label=r'Numerical (convolution)')
    axes[0].set_ylabel('')
    axes[1].legend()
    axes[1].set_xlabel(r'$(t-t_0)/\tau$ [-]')
    return fig, axes

if simpleResponses:
    plot_input_output(t_norm, u_d, x_th=x_d, x_num=x_D, x_conv=x_D2)
    plot_input_output(t_norm, u_s, x_th=x_s, x_num=x_S, x_conv=x_S2)
    plot_input_output(t_norm, u_r, x_th=x_r, x_num=x_R, x_conv=x_R2)
    plot_input_output(t_norm, u_o, x_th=x_o, x_num=x_O, x_conv=x_O2)
    plot_input_output(t_norm, u_h, x_th=x_h, x_num=x_H, x_conv=x_H2)

if randomSignal:
    plot_input_output(t_norm, u_m, x_th=None, x_num=x_M, x_conv=x_M2)

 
# --------------------------------------------------------------------------------}
# ---  
# --------------------------------------------------------------------------------{
if simpleResponses:
    fig,axes = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    axes=np.atleast_2d(axes)
    ax=axes[0,0]
    ax.plot(t_norm, x_d*tau        , label=r'Impulse response ($x \tau$)')
    ax.plot(t_norm, x_s/x_s[-1]    , label=r'Step response ($x/x_\inf$)')
    ax.plot(t_norm, x_r/x_r[-1]    , label=r'Ramp response ($x/x_{10}$)')
    ax.set_xlabel(r'$t/\tau$ [-]')
    ax.set_ylabel('Normalized responses [-]')
    ax.legend()



if __name__ == '__main__':
    plt.show()
if __name__ == '__test__':
    pass
if __name__=="__export__":
    pass
    #from welib.tools.repo import export_figs_callback
    #export_figs_callback(__file__)
