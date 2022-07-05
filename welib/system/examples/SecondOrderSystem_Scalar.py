"""
Example showing the responses of a first order system defined by:

   x_dot(t) + 1/tau x = f(t) 

"""
import numpy as np
import matplotlib.pyplot as plt
#from welib.system.singledof import *
from welib.system.secondorder import *
from welib.tools.signal_analysis import impulse, step, ramp, hat, correlated_signal
from welib.tools.tictoc import Timer

# --- Script parameters
simpleResponses = True
randomSignal    = True

outD=True  # outupt derivatives
m           = 250.0 # system mass
k           = 40.0  # spring constant
c           = 100.0 # damping constant
zeta   = c/(2*np.sqrt(k*m))
omega0 = np.sqrt(k/m)
T0     = 2*np.pi/omega0
print('zeta',zeta, 'omega0',omega0,'T0',T0)
b      = 1/m # NOTE: should be 1/m for a mechanical system

time = np.linspace(-3*T0,10*T0,10001) # time span # Convolution needs enough time steps
q0  = [0,0]         # Initial condition, x[time[0]] # TODO
# --- Parameters for simple inputs
t0  = 0*T0     # initial time where input starts (for simple inputs)
A   = -1         # amplitude for all simple inputs
omega = 30      # for sine 
phi   = np.pi/3 # for sine
T     = 2*T0   # for hat

# --- Derived parameters
t_norm = (time-t0)/T


if simpleResponses:
    # --- Analytical integrations
    x_d = impulse_response(time, omega0, zeta, b=b, t0=t0, A=A, q0=q0, both=outD)
    x_s = step_response   (time, omega0, zeta, b=b, t0=t0, A=A, q0=q0, both=outD)
#     x_r = ramp_response   (time, tau, t0=t0, A=A, q0=q0)
#     x_o = sine_response   (time, tau, t0=t0, A=A, q0=q0, omega=omega, phi=phi)
#     x_h = hat_response    (time, tau, t0=t0, A=A, q0=q0, T=T)

    # --- Numerical integrations of simple functions
    u_d = impulse(time, t0=t0, A=A, epsilon=0.01*T) # NOTE: smooth, and sensitive
    u_s = step   (time, t0=t0, A=A) 
#     u_r = ramp   (time, t0=t0, A=A)
#     u_h = hat    (time, t0=t0, A=A, T=T, method='sum')
#     u_o = A*np.sin(omega*(time-t0)+phi); u_o[time<t0]=0

    # ----  if method in ['LSODA','RK45','RK23','DOP853','Radau','BDF']:
    x_D = integrate(time, omega0, zeta, u_d, b=b, q0=q0, method='RK23' , both=outD) # NOTE sensitive
    x_S = integrate(time, omega0, zeta, u_s, b=b, q0=q0, method='LSODA', both=outD)
#     x_R = integrate(time, tau, u_r, x0=x0, method='LSODA')
    # # x_H = integrate(time, tau, u_h, x0=x0, method='DOP853') # NOTE sensitive, Funny behavior for T=1
#     x_H = integrate(time, tau, u_h, x0=x0, method='RK45') # NOTE sensitive, Funny behavior to T=2
#     x_O = integrate(time, tau, u_o, x0=x0, method='LSODA')

    # --- Numerical integration using Duhamel/convolution
    x_D2 = integrate(time, omega0, zeta, u_d, b=b, q0=q0, method='convolution', both=outD)
    x_S2 = integrate(time, omega0, zeta, u_s, b=b, q0=q0, method='convolution', both=outD) # NOTE: convolution difficult for xdot, need more point
#     x_R2 = integrate(time, tau, u_r, x0=x0, method='convolution')
#     x_H2 = integrate(time, tau, u_h, x0=x0, method='convolution')
#     x_O2 = integrate(time, tau, u_o, x0=x0, method='convolution')

if randomSignal:
    # --- Numerical integrations of correlated, random signal
    u_m = correlated_signal(coeff=0.992, n=len(time), seed=129) # seed provided for reproducibility
    #with Timer('Num'):
    #x_M  = integrate(time, tau, u_m, x0=x0, method='DOP853')
    x_M  = integrate(time, omega0, zeta, u_m, b=b, q0=q0, method='RK23')
    #with Timer('Conv'):
    x_M2 = integrate(time, omega0, zeta, u_m, b=b, q0=q0, method='convolution')


# --------------------------------------------------------------------------------}
# --- Individual plots 
# --------------------------------------------------------------------------------{
def plot_input_output(t_norm, u, x_th=None, x_num=None, x_conv=None, fig=None, axes=None, outD=False):
    if fig is None:
        if outD:
            fig,axes = plt.subplots(3, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
            fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
        else:
            fig,axes = plt.subplots(2, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
            fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    if outD:
        axes=np.atleast_1d(axes)
        axes[0].plot(t_norm , u     , 'k-' , label=r'')
        if x_th is not None:
            axes[1].plot(t_norm , x_th[0], 'k-'  , label=r'Analytical')
            axes[2].plot(t_norm , x_th[1], 'k-'  , label=r'Analytical')
        if x_num is not None:
            axes[1].plot(t_norm , x_num[0] , '--' , label=r'Numerical')
            axes[2].plot(t_norm , x_num[1] , '--' , label=r'Numerical')
        if x_conv is not None:
            axes[1].plot(t_norm , x_conv[0], ':'  , label=r'Numerical (convolution)')
            axes[2].plot(t_norm , x_conv[1], ':'  , label=r'Numerical (convolution)')
        axes[0].set_ylabel('Input')
        axes[1].set_ylabel('Output x')
        axes[2].set_ylabel('Output xdot')
        axes[1].legend()
        axes[1].set_xlabel(r'$(t-t_0)/T$ [-]')
    else:
        axes=np.atleast_1d(axes)
        axes[0].plot(t_norm , u     , 'k-' , label=r'')
        if x_th is not None:
            axes[1].plot(t_norm , x_th  , 'k-'  , label=r'Analytical')
        if x_num is not None:
            axes[1].plot(t_norm , x_num , '--' , label=r'Numerical')
        if x_conv is not None:
            axes[1].plot(t_norm , x_conv, ':'  , label=r'Numerical (convolution)')
        axes[0].set_ylabel('Input')
        axes[1].set_ylabel('Output')
        axes[1].legend()
        axes[1].set_xlabel(r'$(t-t_0)/T$ [-]')
    return fig, axes

if simpleResponses:
    plot_input_output(t_norm, u_d, x_th=x_d, x_num=x_D, x_conv=x_D2, outD=outD)
    plot_input_output(t_norm, u_s, x_th=x_s, x_num=x_S, x_conv=x_S2, outD=outD)
#     plot_input_output(t_norm, u_r, x_th=x_r, x_num=x_R, x_conv=x_R2, outD=outD)
#     plot_input_output(t_norm, u_o, x_th=x_o, x_num=x_O, x_conv=x_O2, outD=outD)
#     plot_input_output(t_norm, u_h, x_th=x_h, x_num=x_H, x_conv=x_H2, outD=outD)

if randomSignal:
    plot_input_output(t_norm, u_m, x_th=None, x_num=x_M, x_conv=x_M2)

 
# --------------------------------------------------------------------------------}
# ---  
# --------------------------------------------------------------------------------{
# if simpleResponses:
#     fig,axes = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
#     fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
#     axes=np.atleast_2d(axes)
#     ax=axes[0,0]
#     ax.plot(t_norm, x_d*tau        , label=r'Impulse response ($x \tau$)')
#     ax.plot(t_norm, x_s/x_s[-1]    , label=r'Step response ($x/x_\inf$)')
# #     ax.plot(t_norm, x_r/x_r[-1]    , label=r'Ramp response ($x/x_{10}$)')
#     ax.set_xlabel(r'$t/\tau$ [-]')
#     ax.set_ylabel('Normalized responses [-]')
#     ax.legend()
# 


if __name__ == '__main__':
    plt.show()
if __name__ == '__test__':
    pass
if __name__=="__export__":
    pass
    #from welib.tools.repo import export_figs_callback
    #export_figs_callback(__file__)
