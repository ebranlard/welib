import unittest

import numpy as np
from welib.yams.utils import *


# --------------------------------------------------------------------------------}
# --- TESTS
# --------------------------------------------------------------------------------{
class TestUtils(unittest.TestCase):
    def test_rot(self):
        # --- Identity matrix for 0 rotation
        np.testing.assert_almost_equal(R_x(0),np.eye(3))

    def test_skew(self):
        # --- Testing  \tilde{u} . v  ==  u x v
        u=np.array([1,2,3])
        v=np.array([-1,0,-3])
        np.testing.assert_equal(np.cross(u,v), np.dot(skew(u),v))

    def test_inertia(self):
        # --- Transferring inertia at G to point A and B and then to each other
        I_G=np.diag([1,2,3]); 
        M=2;
        r_OG = np.array([ 1, 2, 10 ])
        r_OA = r_OG + np.array([5, 8 , 2] )
        r_OB = r_OG + np.array([4, -6, -3])
        r_AG = r_OG-r_OA
        r_BG = r_OG-r_OB
        I_A  = translateInertiaMatrix(I_G,M,r_AG)        # I@ A
        I_B  = translateInertiaMatrix(I_G,M,r_BG)        # I@ B
        I_B2 = translateInertiaMatrix(I_A,M,r_BG,r_AG   )# I@B from A
        I_A2 = translateInertiaMatrix(I_B,M,r_AG,r_BG   )# I@A from B
        np.testing.assert_equal(I_A,I_A2)
        np.testing.assert_equal(I_B,I_B2)

        # --- Transfer of inertia at A to COG then back at A
        I_A = np.eye(3)
        M = 12
        r_GA = np.array([3,4,5])
        I_G  = translateInertiaMatrixToCOG  (I_A,M,r_GA)  # I@G from A
        I_A2 = translateInertiaMatrixFromCOG(I_G,M,-r_GA) # I@A from G
        np.testing.assert_equal(I_A,I_A2)

if __name__=='__main__':
    unittest.main()
