""" 
Time-integrate a mass spring damper system with arbitrary forcing
using Duhamel's integral method, or a time integration scheme

"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Welib
from welib.system.mech_system import MechSystem
from welib.system.statespacelinear import LinearStateSpace
from welib.tools.signal_analysis import correlated_signal, step
from scipy.interpolate import interp1d


# --- Parameters defining the system and simulation
m           = 250.0 # system mass
k           = 40.0  # spring constant
c           = 60.0  # damping constant
nPeriods    = 5
nPerPeriods = 300
useStep = True

# --- Derived parameters
omega0 = np.sqrt(k/m)
T = 2*np.pi/omega0
time = np.linspace(0,nPeriods*T,nPerPeriods*nPeriods+1)

# --- Initialize a linear time invariant system
if useStep:
    F = step(time)
else:
    F = correlated_signal(coeff=0.9, n=len(time), seed=129) # seed provided for reproducibility
    F=F-F[0]

# --- Setup a Mechanical system object, to easily integrate
sys= MechSystem(M=m, C=c, K=k)
sys.setForceTimeSeries(time,F)
print(sys)
resn = sys.integrate(time, method='LSODA')
resd = sys.integrate(time, method='duhamel')

# --- Do the same with the linear state space class
sysl = LinearStateSpace(A=sys.A, B=sys.B)
sysl.setInputTimeSeries(time,F)
print(sysl)
resln = sysl.integrate(time, method='LSODA')
resld = sysl.integrate(time, method='impulse')



# --- Plot
fig, axes = sys.plot(res  = resn , label = 'MechSys Numerical int.')
fig, axes = sys.plot(res  = resd , label = 'MechSys convolution (Duhamel)' ,  fig=fig, axes=axes, ls = '--')
fig, axes = sysl.plot(res = resln, label = 'Lin sys Numerical int.'        ,  fig=fig, axes=axes, ls = ':')
fig, axes = sysl.plot(res = resld, label = 'Lin sys convolution'           ,  fig=fig, axes=axes, ls = '-.', c = 'k')
fig.subplots_adjust(left=0.14, right=0.99, top=0.98, bottom=0.10, hspace=0.20, wspace=0.20)
axes[0].legend()
axes[0].set_title('System - 2nd order - Duhamel or numerical')

# sys.plot_forcing()
# sysl.plot_inputs()

if __name__ == '__main__':
    plt.show()

if __name__ == '__test__':
    pass
    #try:
    #    plt.close()
    #except:
    #    pass
    #pass
if __name__=="__export__":
    pass
    #from welib.tools.repo import export_figs_callback
    #export_figs_callback(__file__)
