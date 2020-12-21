"""  
 Pendulum modelling a pinned rigid body under gravity
 
 Equation of motions:
  
     JO theta_ddot + Cp theta_dot + Kp theta + M *g *l/2 sin\theta =0

"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Local 
from welib.system.mech_system import MechSystem

# --- Parameters
Mtot = 5.4e6 # [kg]     body mass
L    = 150   # [m]      body length
rho  = 7850  # [kg/m^3] material density
g    = 10  # [m/s-2]  acceleration of gravity
kp   = 0e9   # [N/m]  pin stiffness
cp   = 0e9  # [N/m.s] pin damping
theta0    = 10*np.pi/180 # Initial angle    [rad]
thetadot0 =  0*np.pi/180 # Initial velocity [rad/s]

# --- Derived parameters
D  = np.sqrt(4*Mtot/(rho*np.pi*L))
JO = Mtot * L**2/3   # moment of inertia wrt to O (around x or y axis) [kg m^2]
JG = Mtot * L**2/12  # moment of inertia wrt to G (around x or y axis) [kg m^2]
zG = -L/2 # position of center of mass
omega = np.sqrt(Mtot * g * L/2 / JO) # eigen frequency

# Matrices
M=np.array([[JO]])
K=np.array([[kp]])
C=np.array([[cp]])

print('Natural period: T = {:.1f}s'.format(2*np.pi/(omega)))
time  = np.linspace(0, 5*2*np.pi/(omega), 100)

# --- Definition of RHS
def Fx(t,x,xdot):
    """ Returns force on pendulum"""
    theta      = x[0,0]
    Force      = np.zeros((1,1))
    Force[0,0] = - g* Mtot*L/2*np.sin(theta) 
    return Force

# --- Define a system and perform time integration
sys=MechSystem(M, C, K, Fx, x0=[theta0], xdot0=[thetadot0] )
res=sys.integrate(time, method='RK45') # **options):

sys.plot()


if __name__=="__main__":
    plt.show()
if __name__=="__test__":
    pass

