import unittest
import os
import numpy as np    
from welib.yams.rotations import *

MyDir=os.path.dirname(__file__)

class Test(unittest.TestCase):
    def test_EulerP(self):
        # Test misc relationships for Euler parameters

        # See Shabana Example 2.7
        e0 = 0.8776
        e1=e2=e3 = 0.2768
        p=np.array([e0,e1,e2,e3]).reshape(4,1)

        A  = EulerP_A(e0,e1,e2,e3)
        G  = EulerP_G(e0,e1,e2,e3)
        Gb = EulerP_Gb(e0,e1,e2,e3)
        E  = EulerP_E(e0,e1,e2,e3)
        Eb = EulerP_Eb(e0,e1,e2,e3)

        # Transformation matrix
        np.testing.assert_almost_equal(A.dot(A.T), np.eye(3), 3)
            
        # G G.T = 4 I_3
        np.testing.assert_almost_equal(G.dot(G.T), np.eye(3)*4, 3)
        np.testing.assert_almost_equal(Gb.dot(Gb.T), np.eye(3)*4, 3)

        # E E.T = Eb Eb.T = I_3
        np.testing.assert_almost_equal(E.dot(E.T), np.eye(3), 3)
        np.testing.assert_almost_equal(Eb.dot(Eb.T), np.eye(3), 3)

        # E p = Eb p = 0
        np.testing.assert_almost_equal(E.dot(p), np.zeros((3,1)), 3)
        np.testing.assert_almost_equal(Eb.dot(p), np.zeros((3,1)), 3)

        # E.T E = I4 - pp.T
        np.testing.assert_almost_equal(E.T.dot(E), np.eye(4) - p.dot(p.T), 3)
        np.testing.assert_almost_equal(Eb.T.dot(Eb), np.eye(4) - p.dot(p.T), 3)

        # Euler parameters from transformation matrix
        e02,e12,e22,e32 = EulerP_fromA(A)
        np.testing.assert_almost_equal(e02, e0, 3)
        np.testing.assert_almost_equal(e12, e1, 3)
        np.testing.assert_almost_equal(e22, e2, 3)
        np.testing.assert_almost_equal(e32, e3, 3)


    def test_EulerConversions(self):
        # Body ZXZ and Euler P
        phi   = np.pi/4
        theta = np.pi/3
        psi   = np.pi/6
        e0,e1,e2,e3 = BodyZXZ_toEuler(phi,theta,psi)
        phi2, theta2, psi2= EulerP_toBodyZXZ(e0,e1,e2,e3)

        np.testing.assert_almost_equal(phi2, phi)
        np.testing.assert_almost_equal(theta2, theta2)
        np.testing.assert_almost_equal(psi2, psi)

        # Compare with alternative
        E1 = EulerP_fromA(BodyZXZ_A(phi,theta,psi))
        E2 = BodyZXZ_toEuler(phi,theta,psi)
        np.testing.assert_almost_equal(E1, E2)

        # --- XYZ
        A = BodyXYZ_A(0, 0, 0)
        E = EulerP_fromA(A)
        A2= EulerP_A(*E)
        np.testing.assert_almost_equal(A, A2, 5)

        # e0=0
        A = BodyXYZ_A(np.pi, 0, 0)
        E = EulerP_fromA(A)
        A2= EulerP_A(*E)
        np.testing.assert_almost_equal(A, A2, 5)
        np.testing.assert_almost_equal(E[0],0, 5)

        # e0,e1=0
        A = BodyXYZ_A(0, np.pi, 0)
        E = EulerP_fromA(A)
        A2= EulerP_A(*E)
        np.testing.assert_almost_equal(A, A2, 5)
        np.testing.assert_almost_equal(E[0],0, 5)
        np.testing.assert_almost_equal(E[1],0, 5)

        # e0,e1,e2=0
        A = BodyXYZ_A(0, 0, np.pi)
        E = EulerP_fromA(A)
        A2= EulerP_A(*E)
        np.testing.assert_almost_equal(A, A2, 5)
        np.testing.assert_almost_equal(E[0],0, 5)
        np.testing.assert_almost_equal(E[1],0, 5)
        np.testing.assert_almost_equal(E[2],0, 5)

        # Corner case e0,e1~=0
        A = BodyXYZ_A( np.pi/3, 0, np.pi)
        E = EulerP_fromA(A)
        A2= EulerP_A(*E)
        np.testing.assert_almost_equal(A, A2, 3)

        # Misc
        A = BodyXYZ_A(np.pi/5, np.pi/3, np.pi/4)
        E = EulerP_fromA(A)
        A2= EulerP_A(*E)
        np.testing.assert_almost_equal(A, A2, 3)



    def test_BodyZXZ(self):
        phi   = np.pi/4
        theta = np.pi/3
        psi   = np.pi/6

        A     = BodyXYZ_A    (phi, theta, psi)
        G     = BodyXYZ_G    (phi, theta, psi)
        Ginv  = BodyXYZ_Ginv (phi, theta, psi)
        Gb    = BodyXYZ_Gb   (phi, theta, psi)
        Gbinv = BodyXYZ_Gbinv(phi, theta, psi)

        np.testing.assert_almost_equal(A.dot(A.T), np.eye(3), 5)
        np.testing.assert_almost_equal(G.dot(Ginv), np.eye(3), 5)
        np.testing.assert_almost_equal(Gb.dot(Gbinv), np.eye(3), 5)


    def test_BodyXYZ(self):
        phix = np.pi/3
        phiy = np.pi/4
        phiz = np.pi/6

        A     = BodyXYZ_A(phix,phiy,phiz)
        G     = BodyXYZ_G(phix,phiy,phiz)
        Ginv  = BodyXYZ_Ginv(phix,phiy,phiz)
        Gb    = BodyXYZ_Gb(phix,phiy,phiz)
        Gbinv = BodyXYZ_Gbinv(phix,phiy,phiz)

        np.testing.assert_almost_equal(A.dot(A.T), np.eye(3), 5)
        np.testing.assert_almost_equal(G.dot(Ginv), np.eye(3), 5)
        np.testing.assert_almost_equal(Gb.dot(Gbinv), np.eye(3), 5)


    def test_BodyZYX(self):
        phix = np.pi/3
        phiy = np.pi/4
        phiz = np.pi/6

        A     = BodyZYX_A(phix,phiy,phiz)
        G     = BodyZYX_G(phix,phiy,phiz)
        Ginv  = BodyZYX_Ginv(phix,phiy,phiz)
        Gb    = BodyZYX_Gb(phix,phiy,phiz)
        Gbinv = BodyZYX_Gbinv(phix,phiy,phiz)

        np.testing.assert_almost_equal(A.dot(A.T), np.eye(3), 5)
        np.testing.assert_almost_equal(G.dot(Ginv), np.eye(3), 5)
        np.testing.assert_almost_equal(Gb.dot(Gbinv), np.eye(3), 5)


if __name__ == '__main__':
    unittest.main()
    #Test().test_EulerConversions()
