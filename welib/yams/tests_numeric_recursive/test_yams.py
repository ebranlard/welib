import unittest

import numpy as np
from welib.yams.yams import *

# --------------------------------------------------------------------------------}
# --- TESTS
# --------------------------------------------------------------------------------{
class TestYAMS(unittest.TestCase):

    def test_BMat(self):
        # --- Example taken from Main_TNSB from YAMS
        alpha_y  = -1.466329361065083e-02
        B_N_ref = np.array([ 0.0000E+00, 0.0000E+00, 1.0000E+00, 0.0000E+00, alpha_y, 0.0000E+00]).reshape((6,1))
        B_S_ref = np.array([
          [ 1.466171724553691e-01 ,  0.000000000000000e+00],
          [ 0.000000000000000e+00 ,  0.000000000000000e+00],
          [ 1.002150044745554e+00 ,  0.000000000000000e+00],
          [ 0.000000000000000e+00 , -1.466276815184685e-02],
          [-1.466329361065083e-02 ,  0.000000000000000e+00],
          [ 0.000000000000000e+00 ,  9.998924958364900e-01 ]])
        r_NS_ref=np.array([1.466276815184686e-01, 0.000000000000000e+00, -9.998924958364899e+00]).reshape((3,1))

        R_TN     = R_y(alpha_y)        ;
        q_psi    = 1
        z_NS     = - 10
        r_NS_inN = np.array([0, 0, z_NS]).reshape((3,1))
        # link E-T
        R_ET     = np.eye(3)
        # ---------------------------------------------
        # Link T-N
        r_TN     = np.zeros((3,1))
        r_TN[0]  = 1.0000E+02
        Bx_TN    = np.zeros((3,1))
        Bt_TN    = np.zeros((3,1))
        Bx_TN[2] = 1
        Bt_TN[1] = alpha_y
        B_TN     = np.vstack((Bx_TN,Bt_TN))
        B_T      = np.array([])
        B_N      = fBMatRecursion(B_T, B_TN[:3,:], B_TN[3:,:], R_ET, r_TN)
        np.testing.assert_equal(B_N,B_N_ref)
        R_TN=R_y(alpha_y);
        R_EN=np.dot(R_ET,R_TN)
        # ---------------------------------------------
        # Link N-S
        R_NS = R_z(q_psi+np.pi) # Adding pi here , blade down
        R_ES = np.dot(R_EN, R_NS   )
        r_NS = np.dot(R_EN, r_NS_inN )
        np.testing.assert_almost_equal(r_NS,r_NS_ref)

        Bx_NS=np.array([0,0,0]).reshape((3,1))
        Bt_NS=np.array([0,0,1]).reshape((3,1))
        B_NS =np.vstack((Bx_NS,Bt_NS))

        B_S = fBMatRecursion(B_N, B_NS[:3,:], B_NS[3:,:], R_EN, r_NS);
        np.testing.assert_almost_equal(B_S,B_S_ref)


if __name__=='__main__':
    unittest.main()
