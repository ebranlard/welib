
"""
Example showing the responses of a first order system defined by:

   x_dot(t) + 1/tau x = f(t) 

"""
import numpy as np
import matplotlib.pyplot as plt
#from welib.system.singledof import *
from welib.system.secondorder import *
from welib.tools.signal import impulse, step, ramp, hat, correlated_signal
from welib.tools.tictoc import Timer

# --- Script parameters
t0=0
A=1
q0=[0,0]
omega0=2*np.pi
zeta=0.25
b=omega0**2
# 
# T=1
# omega=2*np.pi

outD=False  # outupt derivatives

time = np.linspace(0,10,500)

# x_d = impulse_response(time, omega0, zeta, b=b, t0=t0, A=A, q0=q0, both=outD)
x_s = step_response   (time, omega0, zeta, b=b, t0=t0, A=A, q0=q0, both=outD)


fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
ax.plot(time, x_s    , label='')
ax.set_xlabel('')
ax.set_ylabel('')
ax.legend()

if __name__ == '__main__':
    plt.show()
