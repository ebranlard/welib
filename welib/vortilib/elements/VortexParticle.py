""" 
3D vortex particles! 
For 2D see VortexPoint

"""
import numpy as np
import unittest

#     MINDENOM=0.0
__MINNORM=1e-4

def vp_u_raw(CP, Pv, Alpha, RegFunction, RegParam):
    """
    Induced velocity from one vortex particle on one control point (CP
    """
    DP = CP-Pv
    fourpi_inv=1/(4*np.pi)

    rDeltaP = np.sqrt(DP[0]**2 + DP[1]**2 + DP[2]**2) # norm
    if (rDeltaP<__MINNORM): #--- Exactly on the Singularity 
       return np.zeros(3)
    else: #--- Normal Procedure 
       C = np.zeros(3)
       C[0] = Alpha[1] * DP[2] - Alpha[2] * DP[1]
       C[1] = Alpha[2] * DP[0] - Alpha[0] * DP[2]
       C[2] = Alpha[0] * DP[1] - Alpha[1] * DP[0]
       if RegFunction==0:# No mollification
          r3_inv     = 1./(rDeltaP**3)
          ScalarPart = r3_inv*fourpi_inv
       elif RegFunction==1: # Exponential mollifier
          r3_inv     = 1./(rDeltaP**3)
          E          = np.exp(-rDeltaP**3/RegParam**3)
          ScalarPart = (1.-E)*r3_inv*fourpi_inv
       elif RegFunction==2: # Compact support
          r3_inv     = 1./np.sqrt(RegParam**6+rDeltaP**6)
          ScalarPart = r3_inv*fourpi_inv
       else:
           raise Exception('Wrong regularization function for particles {}'.format(RegFunction))
       return C*ScalarPart


def vp_u(CPs, Pv, Alpha, RegFunction=0, RegParam=0):
    """ Induced velocity from one vortex particle on several control points

    CPs   : n x 3
    Pv    : 3   Position of vortex particle
    Alpha : 3 Intensity of vortex particle

    RegFunction: Regularization function:
                 0: None
                 1: Exponential
                 2: Compact
    OUTPUTS:
        u: (n x 3) velocity, shape of Xcp
    """
    u = np.zeros(CPs.shape)
    # TODO perform numpy vectorization
    for i,CP in enumerate(CPs):
        u[i,:] = vp_u_raw(CP, Pv, Alpha, RegFunction, RegParam)
    return u

# --------------------------------------------------------------------------------}
# --- TESTS
# --------------------------------------------------------------------------------{
class TestVortexParticles(unittest.TestCase):
    def test_VP_RegFunctions(self):
        import warnings
#         warnings.filterwarnings('error')
        # --- One vortex Particle
        # test, 0 on singularity
        Alpha=[0,0,1]
        PPart=np.array([0,0,0])
        U  = vp_u_raw(PPart, PPart, Alpha, RegFunction = 0, RegParam = 0)
        np.testing.assert_equal(U.ravel(), np.zeros(3))
        U  = vp_u_raw(PPart, PPart, Alpha, RegFunction = 1, RegParam = 0)
        np.testing.assert_equal(U.ravel(), np.zeros(3))
        U  = vp_u_raw(PPart, PPart, Alpha, RegFunction = 2, RegParam = 0)
        np.testing.assert_equal(U.ravel(), np.zeros(3))

if __name__ == "__main__":
    unittest.main()


