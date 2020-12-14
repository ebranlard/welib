# See
# J. Peraire, S. Widnall
# 16.07 Dynamics
# Fall 2008
# Version 2.0
# Lecture L30 - 3D Rigid Body Dynamics: Tops and Gyroscopes

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Local 
import sympy
from sympy import symbols
from sympy.physics.mechanics import dynamicsymbols, ReferenceFrame, Point
from sympy.physics.mechanics import inertia, RigidBody
from sympy import init_session
from sympy import init_printing

print=lambda x: sympy.pprint(x, use_unicode=False,wrap_line=False)
# init_session() 
# init_printing() 


theta, phi, psi = dynamicsymbols('theta, phi, psi')

inertial_frame  = ReferenceFrame('O')
gyro_frame   = ReferenceFrame('e')
gyro_frame.orient(inertial_frame, 'Body', (phi,  theta, psi), 'ZXZ')

print('>>>>')
print(gyro_frame.dcm(inertial_frame)) # From Inertial to Euler

print('>>>>')


omega=gyro_frame.ang_vel_in(inertial_frame)  # Angular velocity of a frame in another
print(omega.to_matrix(gyro_frame))


# upper_leg_frame.orient(lower_leg_frame, 'Axis', (theta2, lower_leg_frame.z))
# torso_frame = ReferenceFrame('T')
# torso_frame.orient(upper_leg_frame, 'Axis', (theta3, upper_leg_frame.z))


