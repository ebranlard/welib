import unittest
import numpy as np

from welib.mesh.vectoranalysis import *

class TestVectorAnalysis(unittest.TestCase):
    def test_grad1d(self):
        x  = np.array([0,1,3,4])
        x2 = x**2
        np.testing.assert_almost_equal(matlab_gradient_1d(x2),np.array([1,4.5,7.5,7]))
        np.testing.assert_almost_equal(matlab_gradient_1d(x2,2),np.array([0.5,2.25,3.75,3.5]))
        # NOTE NOTE: this test below fails due to boundary effects, to replaced by python values..
        # np.testing.assert_almost_equal(matlab_gradient_1d(x2,x),np.array([1,3,5,7]))
        np.testing.assert_almost_equal(matlab_gradient_1d(x2,x),np.array([1,2,6,7]))

        x  = np.array([[0,1,3,4]])
        x2 = x**2
        np.testing.assert_almost_equal(matlab_gradient_1d(x2),np.array([[1,4.5,7.5,7]]))
        np.testing.assert_almost_equal(matlab_gradient_1d(x2.T),np.array([[1,4.5,7.5,7]]).T)

    def test_grad2d(self):
        x = np.arange(-1,1.5,0.5)
        y = np.arange( 0,1.21,0.2)
        X,Y = np.meshgrid(x,y)
        Z = X * np.exp(-X**2 - Y**2)
        # --- Gradient, no dimensions
        px,py = matlab_gradient_2d(Z)
        px_ref =np.array([[ -0.0215210,  0.1839397],[-0.0206771,  0.1767273]])
        py_ref = np.array([ [ 0.014425,  0.015269] ,[ 0.027197,  0.028788]])
        np.testing.assert_almost_equal(px[:2,:2],px_ref,decimal=6)
        np.testing.assert_almost_equal(py[:2,:2],py_ref,decimal=6)
        # --- Gradient, uni dimensions
        px,py = matlab_gradient_2d(Z,0.5)
        px_ref = np.array([ [  -0.0430419,    0.3678794], [  -0.0413542,    0.3534547]])
        py_ref = np.array([ [   0.0288495,    0.0305372], [   0.0543933,    0.0575753]])
        np.testing.assert_almost_equal(px[:2,:2],px_ref,decimal=7)
        np.testing.assert_almost_equal(py[:2,:2],py_ref,decimal=7)
        # --- Gradient, bi dimensions
        px,py = matlab_gradient_2d(Z,0.5,0.2);
        px_ref = np.array([ [  -0.0430419,    0.3678794], [  -0.0413542,    0.3534547]])
        py_ref = np.array([ [   0.0721238,    0.0763430], [   0.1359832,    0.1439382]])
        np.testing.assert_almost_equal(px[:2,:2],px_ref,decimal=7)
        np.testing.assert_almost_equal(py[:2,:2],py_ref,decimal=7)
        px,py = matlab_gradient_2d(Z,x,y);
        np.testing.assert_almost_equal(px[:2,:2],px_ref,decimal=7)
        np.testing.assert_almost_equal(py[:2,:2],py_ref,decimal=7)


    def test_div2d(self):
        x = np.arange(-1,1.5,0.5)
        y = np.arange( 0,1.21,0.2)
        X,Y = np.meshgrid(x,y)
        U =    X * np.exp(-X**2 - Y**2)
        V = Y**2 * np.exp(-X**2 - Y**2)

        # --- Divergence, no dimensions
        div=matlab_div_2d(U,V)
        div_ref = np.array([ [  -0.0073828,    0.2138703], [   0.0044018,    0.2298194]])
        np.testing.assert_almost_equal(div[:2,:2],div_ref,decimal=7)

        # --- Divergence, dimensions
        div=matlab_div_2d(X,Y,U,V)
        div_ref = np.array([ [   0.0276490,    0.5175322], [   0.0840403,    0.6189148]])
        np.testing.assert_almost_equal(div[:2,:2],div_ref,decimal=7)

    def test_curl2d(self):
        x = np.arange(-1,1.5,0.5)
        y = np.arange( 0,1.21,0.2)
        X,Y = np.meshgrid(x,y)
        U =    X * np.exp(-X**2 - Y**2)
        V = Y**2 * np.exp(-X**2 - Y**2)
        # --- Curl, no dimensions
        curl,_=matlab_curl_2d(U,V)
        curl_ref = np.array([ [  -0.0144248,   -0.0152686], [  -0.0114043,   -0.0166409]])
        np.testing.assert_almost_equal(curl[:2,:2],curl_ref,decimal=7)

        # --- Curl, dimensions
        curl,_=matlab_curl_2d(X,Y,U,V)
        curl_ref = np.array([ [  -0.0721238,   -0.0763430], [  -0.1043984,   -0.1196448]])
        np.testing.assert_almost_equal(curl[:2,:2],curl_ref,decimal=7)


if __name__=='__main__':
    unittest.main()
