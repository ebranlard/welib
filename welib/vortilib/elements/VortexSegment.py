"""
Reference:
    [1] E. Branlard - Wind Turbine Aerodynamics and Vorticity Based Method, Springer, 2017
"""
# --- General
import unittest
import numpy as np

# --------------------------------------------------------------------------------{
def vs_u_raw(CP, Pa, Pb, Gamma, RegFunction=0, RegParam=0, nt=None, RegParamW=None):
    """ Induced velocity from a vortex segment on one control point
    See fUi_VortexSegment11_smooth

    CPs : 3 x n

    RegFunction: Regularization function:
                 0: None
                 1: Rankine
                 2: Lamb-Oseen
                 3: Vatistas
                 4: Denominator offset
    """
    DPa = CP.ravel()-Pa.ravel()
    DPb = CP.ravel()-Pb.ravel()
    xa, ya, za = DPa[0], DPa[1], DPa[2]
    xb, yb, zb = DPb[0], DPb[1], DPb[2]

    norm_a      = np.sqrt(xa * xa + ya * ya + za * za)
    norm_b      = np.sqrt(xb * xb + yb * yb + zb * zb)
    denominator = norm_a * norm_b * (norm_a * norm_b + xa * xb + ya * yb + za * zb)

    if (denominator < 1e-17):
        return np.zeros((1,3))
    if (norm_a < 1e-08 or norm_b < 1e-08):
        return np.zeros((1,3))

    crossprod       = np.array([[ya * zb - za * yb, za * xb - xa * zb, xa * yb - ya * xb]])
    # Singular model
    if RegFunction==0:
        Kv = Gamma / (4.0 * np.pi) * (norm_a + norm_b) / denominator
        return Kv * crossprod
    # --- Regularization
    if nt is None:  
        # Regularization models, based on orthogonal distance to segment h2
        norm2_r0        = (xa - xb)**2 + (ya - yb)**2 + (za - zb)**2
        norm2_crossprod = crossprod[0,0]**2 + crossprod[0,1]**2 + crossprod[0,2]**2
        h2              = norm2_crossprod/norm2_r0 # Orthogonal distance (r1 x r2)/r0
        eps2 = h2/RegParam**2
        if RegFunction==1:
            if (eps2 < 1):
                Kv = eps2
            else:
                Kv = 1.0
        elif RegFunction==2:
            Kv = 1.0 - np.exp(-1.25643 * eps2)
        elif RegFunction==3:
            Kv = eps2 / np.sqrt(1 + eps2**2)
        elif RegFunction==4:
            Kv = 1.0
            denominator = denominator + RegParam**2 * norm2_r0
        else:
            raise NotImplementedError('Regularization for segments {}'.format(RegFunction))
    else:
        # --- Using 2D Gaussian
        nt = np.asarray(nt)
        es = np.array([[(xa - xb), (ya - yb), (za - zb)]])
        nw = np.cross(es,nt)
        nw = nw/np.linalg.norm(nw)
        rt = np.dot(DPa.ravel(),nt.ravel()) # distance along nt component
        rw = np.dot(DPa.ravel(),nw.ravel()) # distance along nw component

        eps2 = rt**2/RegParam**2 + rw**2/RegParamW**2
        #eps2 = (rt**2+rw**2)/RegParam**2 

        if RegFunction==2:
            Kv = 1.0 - np.exp(-1.25643 * eps2 )
        elif RegFunction==3:
            Kv = eps2 / np.sqrt(1 + eps2**2)
        else:
            raise NotImplementedError('Regularization (2D Gaussian) for segments {}'.format(RegFunction))

    Kv = Gamma * Kv / (4.0 * np.pi) * (norm_a + norm_b) / denominator





    return Kv * crossprod

def vs_u(Xcp, Ycp, Zcp, Pa, Pb, Gamma, RegFunction=0, RegParam=0, nt=None, RegParamW=None):
    """ Induced velocity from a vortex segment on several control point
    See fUi_VortexSegment11_smooth

    CPs : n x 3

    RegFunction: Regularization function:
                 0: None
                 1: Rankine
                 2: Lamb-Oseen
                 3: Vatistas
                 4: Denominator offset
    OUTPUTS:
        ux, uy, uz: velocity, shape of Xcp
    """
    Xcp = np.asarray(Xcp)
    shape_in = Xcp.shape
    Xcp = Xcp.ravel()
    Ycp = np.asarray(Ycp).ravel()
    Zcp = np.asarray(Zcp).ravel()
    ux  = np.zeros(Xcp.shape)
    uy  = np.zeros(Xcp.shape)
    uz  = np.zeros(Xcp.shape)
    for i,(x,y,z) in enumerate(zip(Xcp,Ycp,Zcp)):
        CP=np.array([[x,y,z]])
        u = vs_u_raw(CP,Pa,Pb, Gamma, RegFunction, RegParam, nt, RegParamW)
        ux[i] = u[0,0]
        uy[i] = u[0,1]
        uz[i] = u[0,2]
        
    ux = ux.reshape(shape_in)
    uy = uy.reshape(shape_in)
    uz = uz.reshape(shape_in)
    return ux,uy,uz



# --------------------------------------------------------------------------------}
# --- TESTS
# --------------------------------------------------------------------------------{
class TestVortexSegment(unittest.TestCase):
    def test_VS_RegFunctions(self):
        import warnings
#         warnings.filterwarnings('error')
        # --- One vortex segment
        Pa = np.array([[ 0, 0, -z0]])
        Pb = np.array([[ 0, 0,  z0]])
        # --- test, 0 on singularity
        U  = vs_u_raw(np.array([0,0,0]), Pa, Pb, Gamma = 1, RegFunction = 0, RegParam = 0)
        np.testing.assert_equal(U, np.zeros((1,3)))
        U  = vs_u_raw(Pa, Pa, Pb, Gamma = 1, RegFunction = 0, RegParam = 0)
        np.testing.assert_equal(U, np.zeros((1,3)))
        U  = vs_u_raw(Pb, Pa, Pb, Gamma = 1, RegFunction = 0, RegParam = 0)
        np.testing.assert_equal(U, np.zeros((1,3)))

if __name__ == "__main__":
    unittest.main()


