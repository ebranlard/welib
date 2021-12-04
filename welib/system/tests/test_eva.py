""" 
Tests for eigenvalue analyses
"""
import unittest
import numpy as np    
import os
from welib.system.eva import *

MyDir=os.path.dirname(__file__)

class Test(unittest.TestCase):

    def test_polyeig(self):
        M = np.diag([3,1,3,1])
        C = np.array([[0.4 , 0 , -0.3 , 0] , [0 , 0  , 0 , 0] , [-0.3 , 0 , 0.5 , -0.2 ] , [ 0 , 0 , -0.2 , 0.2]])
        K = np.array([[-7  , 2 , 4    , 0] , [2 , -4 , 2 , 0] , [4    , 2 , -9  , 3    ] , [ 0 , 0 , 3    , -3]])
        X,e = polyeig(K,C,M)
        #print('X:\n',X)
        #print('e:\n',e)
        # Test that first eigenvector and value satisfy eigenvalue problem:
        s = e[0];
        x = X[:,0];
        res = (M*s**2 + C*s + K).dot(x) # residuals

        self.assertTrue(all(np.abs(res)<1e-12))

        np.testing.assert_almost_equal(e[0], -2.44985, 4)

    def test_eigMCK(self):
        # --- Simple test
        M = np.array([[100.]])
        K = np.array([[10.]])
        omega = np.sqrt(K[0,0]/M[0,0])
        zeta0 = 0.9
        c = 2*M[0,0]*omega*zeta0
        C = np.array([[c]])
        #  Method 1
        freq_d,zeta,Q,freq = eigMCK(M, C, K, method='diag_beta')


        np.testing.assert_almost_equal(freq[0], omega/(2*np.pi), 4)
        np.testing.assert_almost_equal(zeta[0], zeta0, 4)
        np.testing.assert_almost_equal(freq_d[0], omega/(2*np.pi)*np.sqrt(1-zeta0**2), 2)

        #  Method 2
        freq_d,zeta,Q,freq = eigMCK(M, C, K, method='full_matrix')

        np.testing.assert_almost_equal(freq[0], omega/(2*np.pi), 4)
        np.testing.assert_almost_equal(zeta[0], zeta0, 4)
        np.testing.assert_almost_equal(freq_d[0], omega/(2*np.pi)*np.sqrt(1-zeta0**2), 2)

        # --- 
        M = np.diag([3.,1.,3.,1.])
        C = np.array([[0.4 , 0 , 0 , 0] , [0 , 1., 0 , 0] , [0 , 0 , 0.5 , 0 ] , [ 0 , 0 , 0 , 0.2]])
        K = np.array([[ 7. , 0 , 0 , 0] , [0 , 2., 0 , 0] , [0 , 0 ,  9. , 0 ] , [ 0 , 0 , 0 ,  3.]])
        #M = np.diag([100,100,100,100])
        #K = np.diag([10,10,10,10])
        #C = np.diag([0.1*c,0.2*c,0.5*c,1*c])

        # --- Method 1
        freq_d,zeta,Q,freq = eigMCK(M, C, K, method='diag_beta')
        np.testing.assert_almost_equal(freq[0], 0.2250, 4)
        np.testing.assert_almost_equal(freq[3], 0.2756, 4)
        np.testing.assert_almost_equal(freq_d[0], 0.21054, 4)
        np.testing.assert_almost_equal(zeta[0], 0.35355, 4)

        # --- Method 2
        freq_d,zeta,Q,freq = eigMCK(M, C, K, method='full_matrix')
        np.testing.assert_almost_equal(freq[0], 0.2250, 4)
        np.testing.assert_almost_equal(freq[3], 0.2756, 4)
        np.testing.assert_almost_equal(freq_d[0], 0.21054, 4)
        np.testing.assert_almost_equal(zeta[0], 0.35355, 4)


if __name__ == '__main__':
    unittest.main()
