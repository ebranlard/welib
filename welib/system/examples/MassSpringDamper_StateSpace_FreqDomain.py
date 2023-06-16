""" 
Frequency response of a mass spring damper system, put into state space form, 
given two forcing inputs

"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Welib
from welib.system.statespacelinear import *
import sympy as sp


# --- Parameters for Bode plot
omega     = 10**(np.linspace(-2,2, 1000)) # frequencies used to evaluate the transfer function "analytically"
omega_num = 10**(np.linspace(-1,1, 3 ))   # frequencies used to evaluate the transfer function "numerically"
tmin      = 30                            # minimum time before numerical evaluation of frequency response for "numerical" evaluation

# --- Parameters defining the system
m = 250.0           # system mass
k = 40.0            # spring constant
b = 50.0            # damping constant
# System matrices
A = np.array([[0, 1.], [-k/m, -b/m]])
B = np.array([[0,0], [1/m, 100/m]])
# C = np.eye(2)  # ouput 1 is position, output 2 is velocity
C = np.zeros((3,2))  # ouput 1 is position, output 2 is velocity, output 3 is acceleration
D = np.zeros((3,2))
C[0,0], C[1,1] = 1, 1
C[2,:] = A[1,:] 
D[2,:] = B[1,:] 
sY =[r'$y_1=x$'    , r'$y_2=\dot{x}$', r'$y_3=\ddot{x}$']
sU =[r'$u_1=F_1$'  , r'$u_2=F_2$'] 

# --- Derived parameters
nU, nY = B.shape[1], C.shape[0] 
Colrs = np.array(plt.rcParams['axes.prop_cycle'].by_key()['color'])[:nY*nU].reshape(nY,nU) # colors
LS =['-','--'] # line style, per inputs
MK =['o','.']  # markeres, per inputs

# --- Initialize a linear time invariant system
sys= LinearStateSpace(A,B,C,D)
print(sys)


# --- "Analytical" frequency response (better)
# (G: magnitude, phi: phase)
G, phi = sys.frequency_response(omega, deg=True)
freq   = omega/(2*np.pi)

# --- "Numerical" frequency response (just for fun)
G_num , phi_num = sys.frequency_response(omega_num, method='numerical', deg=True)
freq_num   = omega_num/(2*np.pi)

# --- Bode Plot
fig,axes = plt.subplots(2, 1, sharex=True, figsize=(6.4,4.8)) # (6.4,4.8)
fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
ax = axes[0]
for iu in range(nU): # loop on inputs
    for iy in range(nY): # loop on outputs
        # Amplitudes
        axes[0].plot(freq_num, G_num[iy,iu,:]  ,    MK[iu], color=Colrs[iy,iu]  )
        axes[0].plot(freq    , G    [iy,iu,:]  , ls=LS[iu], color=Colrs[iy,iu], label='{} -> {}'.format(sU[iu], sY[iy]))
        # Phases
        axes[1].plot(freq_num, phi_num[iy,iu,:],    MK[iu], color=Colrs[iy,iu])
        axes[1].plot(freq    , phi    [iy,iu,:], ls=LS[iu], color=Colrs[iy,iu])
axes[0].set_ylabel('Amplitude ratio [-]')
axes[1].set_ylabel('Phase [deg]')    
axes[0].set_yscale('log')
axes[0].legend(ncol=2)
for ax in axes:
    ax.tick_params(direction='in')
    ax.set_xscale('log')
    ax.set_xlabel('Frequency [Hz]')
axes[0].set_title('System -  LTI Bode plot - 2nd order mass spring damper')    




if __name__ == '__main__':
    plt.show()
if __name__ == '__test__':
    pass
if __name__=="__export__":
    from welib.tools.repo import export_figs_callback
    export_figs_callback(__file__)
