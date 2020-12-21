"""
Generic bodies classes
These classes will be used for more advanced classes:
    - new and old YAMS body classes for Sympy
    - YAMS body for numerical yams
"""
from welib.yams.utils import translateInertiaMatrixToCOG, translateInertiaMatrixFromCOG
from welib.yams.utils import rigidBodyMassMatrix 
# from welib.yams.utils import skew


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
    def __init__(self, name='', r_O=None):
        self.name = name
        self._r_O  = None       # position of body origin in global coordinates
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
    def __init__(self, name, mass, J, s_OG, R_b2g=np.eye(3), s_OP=None, r_O=[0,0,0]):
        """
        Creates a rigid body 

        INPUTS:
         - name: name of body (string)
         - mass: body mass (float)
         - J: inertia tensor (array-like) in body frame, defined either at:
               - center of mass G, located at s_OG from the body origin
               - OR point P, located at s_OP from the body origin
             J may be defined as:
               - a 3x3 matrix
               - a 3-vector (Jxx, Jyy, Jzz) representing the diagonal values
               - a 6-vector (Jxx, Jyy, Jzz, Jxy, Jyz, Jzx) representing the diagonal values
         - s_OG: vector from body origin to body COG in body frame 

         - s_OP: vector from body origin to point where inertia is defined,
                 in body frame
                 (only if inertia is not defined at COG).
         - r_O: vector from global origin to body origin, in global coordinates
         - R_b2g : transformation matrix from body to gobal coordinates

        """
        Body.__init__(self, name)
        self.mass  = mass
        self._r_O   = np.asarray(r_O).ravel()
        self._s_OG = np.asarray(s_OG).ravel()
        self._R_b2g = np.asarray(R_b2g)

        # Ensuring a 3x3 inertia matrix
        J = np.asarray(J)
        Jflat=J.ravel()
        if len(Jflat)==3:
            J = np.diag(Jflat)
        elif len(Jflat)==6:
            J = np.diag(Jflat[:3])
            J[0,1]=J[1,0]=Jflat[3]
            J[1,2]=J[2,1]=Jflat[4]
            J[1,3]=J[3,1]=Jflat[5]
            
        # inertia at COG
        if s_OP is not None:
            s_PG=  self._s_OG - s_OP
            self._J_G = translateInertiaMatrixToCOG(J, mass, s_PG)
        else:
            self._J_G = J

    @property    
    def pos_global(self):
        """ Position of origin in global coordinates """
        return self._r_O

    @pos_global.setter
    def pos_global(self, r_O):
        self._r_O = r_O
        
    @property
    def R_b2g(self):
        """ Transformation matrix from body to global """
        return self._R_b2g

    @R_b2g.setter
    def R_b2g(self, R_b2g):
        self._R_b2g = R_b2g

    @property
    def R_g2b(self):
        """ Transformation matrix from global to body """
        return self._R_b2g.transpose() 

    def pos_local(self, x_gl):
        """ return position vector from origin of body, in body coordinates, of a point in global """
        return self.R_g2b.dot(x_gl - self._r_O)

    # --------------------------------------------------------------------------------
    # --- Inertia
    # --------------------------------------------------------------------------------
    @property    
    def masscenter(self):
        """ Position of mass center in body frame"""
        return self._s_OG

    @property
    def masscenter_pos_global(self):
        """ return masscenter position from inertial frame """
        return self._r_O + self.R_b2g.dot(self._s_OG)

    @property    
    def inertia(self):
        return self.inertia_at([0,0,0])

    @property    
    def masscenter_inertia(self):
        """ Returns inertia matrix at COG in body frame"""
        return self._J_G

    def inertia_at(self, s_OP, R_f2g=None):
        """ returns body inertia at a given point, and given frame (default body frame)
        INPUTS:
         - s_OP: point coordinates from body origin in body coordinates
         - R_f2g: transformation matrix from a given frame when inertia is wanted to global
        """
        # 
        s_GP =   np.asarray(s_OP) - self._s_OG
        J = translateInertiaMatrixFromCOG(self._J_G, self.mass, s_GP)
        if R_f2g is not None:
            R_b2f = np.dot(R_f2g.T, self.R_b2g)
            J = R_b2f.dot(J).dot(R_b2f.T)
        return J

    @property
    def mass_matrix(self):
        """ Body mass matrix at origin"""
        return rigidBodyMassMatrix(self.mass, self._J_G, self._s_OG)

    def __repr__(self):
        s='<{} {} object>:\n'.format(type(self).__name__, self.name)
        s+=' * pos_global:            {} (origin)\n'.format(np.around(self.pos_global,6))
        s+=' * masscenter:            {} (body frame)\n'.format(np.around(self.masscenter,6))
        s+=' * masscenter_pos_global: {} \n'.format(np.around(self.masscenter_pos_global,6))
        s+=' - mass:         {}\n'.format(self.mass)
        s+=' * R_b2g: \n {}\n'.format(self.R_b2g)
        s+=' * masscenter_inertia: \n{}\n'.format(np.around(self.masscenter_inertia,6))
        s+=' * inertia: (at origin)\n{}\n'.format(np.around(self.inertia,6))
        s+='Useful getters: inertia_at, mass_matrix\n'
        return s

    def combine(self, other, name=None, R_b2g=np.eye(3), r_O=None):
        """ Combine two rigid bodies and form a new rigid body
        
        """
        M = self.mass + other.mass
        x_G = (self.mass * self.masscenter_pos_global + other.mass * other.masscenter_pos_global)/M

        if name is None:
            name=self.name + other.name

        # Inertias in new body frame and at new COG
        s_O1_G = self.pos_local(x_G)
        s_O2_G = other.pos_local(x_G)
        J1 = self.inertia_at(s_O1_G, R_b2g)
        J2 = other.inertia_at(s_O2_G, R_b2g)
        #print('s_O1_G ',s_O1_G)
        #print('s_O2_G ',s_O2_G)
        #print('J1\n ',J1)
        #print('J2\n ',J2)
        #print('J12\n ',J1+J2)

        if r_O is None:
            # Putting origin of new body at COG of common body
            r_O  = x_G
            s_OG = [0,0,0]
        else:
            s_OG = (R_b2g.T).dot(x_G-r_O)
        return RigidBody(name, M, J1+J2, s_OG, r_O=r_O, R_b2g=R_b2g)


# --------------------------------------------------------------------------------}
# --- Flexible Body 
# --------------------------------------------------------------------------------{
class FlexibleBody(Body):
    def __init__(self, name):
        """
        Creates a Flexible body 
        """
        Body.__init__(self, name)
