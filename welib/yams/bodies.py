"""
Generic bodies classes
These classes will be used for more advanced classes:
    - new and old YAMS body classes for Sympy
    - YAMS body for numerical yams
"""


__all__ = ['Body','InertialBody','RigidBody','FlexibleBody']

# --- For harmony with sympy
import numpy as np
from numpy import eye, cross, cos ,sin
def Matrix(m):
    return np.asarray(m)
def zeros(m,n):
    return np.zeros((m,n))




# --------------------------------------------------------------------------------}
# --- Generic Body 
# --------------------------------------------------------------------------------{
class Body(object):
    def __init__(self, name=''):
        self.r_O  = None       # position of body origin in global coordinates
        self.mass = None       # mass
        self.MM   = zeros(6,6) # Mass Matrix

    def __repr__(B):
        s='<Generic {} object>:\n'.format(type(self).__name__)
        return s

    @property
    def Mass(self):
        print('Warning: `Mass` is an old interface, use `mass` instead')
        return self.mass

# --------------------------------------------------------------------------------}
# --- Ground Body 
# --------------------------------------------------------------------------------{
class InertialBody(Body):
    def __init__(self):
        Body.__init__(self, name='Grd')

# --------------------------------------------------------------------------------}
# --- Rigid Body 
# --------------------------------------------------------------------------------{
class RigidBody(Body):
    def __init__(self, name, mass, J_G, rho_G):
        """
        Creates a rigid body 
        """
        Body.__init__(self, name)
        self.mass = mass
        #self.s_G_inB = rho_G
        #self.J_G_inB = J_G
        #self.J_O_inB = translateInertiaMatrixFromCOG(self.J_G_inB, mass, self.s_G_inB)

# --------------------------------------------------------------------------------}
# --- Flxible Body 
# --------------------------------------------------------------------------------{
class FlexibleBody(Body):
    def __init__(self, name):
        """
        Creates a rigid body 
        """
        Body.__init__(self, name)
