""" 
Integrate the equation of motions for a pendulum "on a cart" with prescribed motion of the cart

    [Cart]  <--> Prescribed motion
      o
     /
    / pendulum (rod)
   /

"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from welib.system.mech_system import MechSystem
import welib.weio as weio

# --- Parameters
l        = 5      # Rod length
kp_over_m= 0*10   # Spring stifness at card/rod interface (divided by mass)
zeta     = 0.01   # Damping ratio for card/rod interface
g        = 1*9.8  # acceleration of gravity
theta_0  = 00     # intial position of the pendulum [deg]
A_forced = 1.000  # Amplitude of forced displacements
T_forced = 3.0    # Period of forced displacements
fps      = 2      # Kind of frame per seconds for animation
# --- Derived params
mp    = 1                            # Mass of the rod (only relevant to scale stiffness and damping)
kp    = kp_over_m * mp
JzzA  = mp*l**2/3                    # Mass moment of inertia of the rod
Omega = 2*np.pi/T_forced             # Frequency of forced "displacements"
dt    = T_forced/50
k_eff = (kp + mp*l/2*g)
omega0=np.sqrt(k_eff/JzzA)
cp    = zeta*2.0*np.sqrt(k_eff*JzzA) # Damping at card/rod interface
time  = np.arange(0, 50, dt)
# Matrices
M=np.array([[JzzA]])
K=np.array([[kp]])
C=np.array([[cp]])

print('Natural period: T = {:.1f}s'.format(2*np.pi/(omega0)))
print('Forced period:  T = {:.1f}s'.format(T_forced))


# --- Definition of prescribed motion and force
def motion(t): 
    """ returns cart pos vel and acc (harmonic oscillation)"""
    return A_forced*np.sin(Omega*t),A_forced*Omega*np.cos(Omega*t),-A_forced*Omega**2*np.sin(Omega*t)

def Fx(t,x,xdot):
    """ Returns force on pendulum"""
    theta      = x[0,0]
    Force      = np.zeros((1,1))
    _,_,acc    = motion(t)
    Force[0,0] = mp*l/2* (-acc* np.cos(theta) - g * np.sin(theta) )
    return Force

# --- Define a system and perform time integration
sys=MechSystem(M,C,K)
#sys.setForceTimeSeries(time,Fx_lin)
sys.setForceFunction(Fx)
sys.setInitialConditions([theta_0*np.pi/180],[0])
res,_=sys.integrate(time, method='RK45') # **options):

# --- Post process
theta = res.y[0,:]
x_cart= motion(time)[0]
x_end = l*np.sin(theta)  + x_cart 
y_end =-l*np.cos(theta) 

if __name__ == '__main__':
    # --- Plot
    fig,axes = plt.subplots(2, 1, sharey=False, figsize=(6.4,5.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    axes[0].plot([-A_forced, A_forced],[0,0],'k--')
    axes[0].set_xlabel('x [m]')
    axes[0].set_ylabel('y [m]')
    axes[0].tick_params(direction='in')
    axes[1].plot(res.t,x_end,        label='x tip')
    axes[1].plot(res.t,x_cart, 'k-', label='x cart')
    axes[1].plot(res.t,res.y[0,:], label='theta')
    axes[1].set_xlabel('time [s]')
    axes[1].set_ylabel('')
    axes[1].legend()
    axes[1].tick_params(direction='in')
    lnc1, = axes[0].plot([], [], 'ks',lw=2,ms=9)
    lnc2, = axes[1].plot([], [], 'ks',lw=2)
    lnt1, = axes[0].plot([], [], 'o',lw=2)
    lnt2, = axes[1].plot([], [], 'bo',lw=2)
    lnp,  = axes[0].plot([], [], 'k-',lw=2)

    def init():
        axes[0].set_xlim(-l-abs(A_forced), l+abs(A_forced))
        axes[0].set_ylim(-l*1.1, l*0.2)
        axes[0].set_aspect('equal')
        return lnp,lnc1,lnc2,lnt1,lnt2,

    def update(i):
        i=fps*int(i)
        lnp.set_data([x_cart[i],x_end[i]],[0,y_end[i]])
        lnt1.set_data([x_end[i]],[y_end[i]])
        lnc1.set_data([x_cart[i]],[0])
        lnt2.set_data([time[i]],[x_end[i] ])
        lnc2.set_data([time[i]],[x_cart[i]])
        return lnp,lnc1,lnc2,lnt1,lnt2
    ani = FuncAnimation(fig, update, frames=np.arange(0,len(time)/fps), init_func=init, blit=True)
    plt.show()

if __name__ == '__test__':
    pass



