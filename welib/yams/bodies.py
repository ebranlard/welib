"""
Generic bodies classes
These classes will be used for more advanced classes:
    - new and old YAMS body classes for Sympy
    - YAMS body for numerical yams
"""


__all__ = ['Body','InertialBody','RigidBody','FlexibleBody']


class Body(object):
    def __init__(self, name=''):
        self.r_O  = None # position of body origin in global coordinates
        self.mass = None # mass

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
class GroundBody(Body):
    def __init__(B):
        super(GroundBody,B).__init__('Grd')

# --------------------------------------------------------------------------------}
# --- Rigid Body 
# --------------------------------------------------------------------------------{
class RigidBody(Body):
    def __init__(B, Name, Mass, J_G, rho_G):
        """
        Creates a rigid body 
        """
        super(RigidBody,B).__init__()
        B.s_G_inB = rho_G
        B.J_G_inB = J_G
        B.J_O_inB = fTranslateInertiaMatrixFromCOG(B.J_G_inB, Mass, B.s_G_inB)
        B.MM = fGMRigidBody(Mass,B.J_O_inB,B.s_G_inB)
        B.DD = np.zeros((6,6))
        B.KK = np.zeros((6,6))
