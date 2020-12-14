""" 
Variables used for the monopile

"""

import numpy as np
from scipy.optimize import fsolve

# --- Parameters
L   = 100         # [m] Total monopile length
h   = 50          # [m] Wather depth
m   = 9000        # [kg/m] Mass per length
rho = 7850        # [kg/m^3]
E   = 2.1e+11 # [N/m2]
EI  = 2.0e12     # N m^2
# 
# 
D = 8.1966401
t = 0.04476779


# --- Solve for D and t
print('m  (calc):',rho * np.pi *(D*t-t**2))
print('EI (calc):',E*np.pi/64 *(D**4 - (D-2*t)**4))

I_wanted=EI/E
D_t_solve = fsolve(
        lambda x : [
            np.pi/64 *(x[0]**4 - (x[0]-2*x[1])**4) - I_wanted, 
            rho*np.pi *(x[0]*x[1]-x[1]**2) - m
        ], [8,0.045])
print('D,t:',D_t_solve)
# 
# I_res=np.pi/64*(D**4 - (D-2*t_solve)**4)
# print('I_res:  ', I_res)
# print('I_want: ', I_wanted)

# 
# 
# HydroCoeffs:
# CD = 1.0
# CM = 2.0 
# CP =???
# 
# Wave kinematics
T=12s
f=1/12 # 0.083333333
rho_water = 1025 
# wave height H=6m
# 
