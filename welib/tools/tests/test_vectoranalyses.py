import unittest
import numpy as np
from welib.tools.fields import *


class TestVectAnalyses(unittest.TestCase):
    def test_int_grad1d(self):
        # Compute 1D grad of a field, integrate and compare the result with original field
        dx=0.01
        xp = np.arange(0,1+dx/2,dx)
        f  = np.exp(xp)
        fx = np.exp(xp)
        fx  = np.gradient(f, xp)      
        fhat = intgrad(fx,xp, constant = 1)
        np.testing.assert_almost_equal(f, fhat, 10)

    def test_int_grad2d(self):
        # Compute 2D grad of a field, integrate and compare the result with original field
        dx=0.01
        dy=0.02
        xp     = np.arange(0,1+dx/2,dx)
        yp     = np.arange(0,1+dx/2,dy)
        # NOTE: mesh grid has shape: len(y) x len(x), i.e. first dim is y
        x,y    = np.meshgrid(xp,yp)       
        f      = np.exp(x+y) + np.sin((x-2*y)*3)
        fy,fx  = np.gradient(f, yp, xp)      
        fhat = intgrad( (fx,fy), (xp,yp), constant=1.0)
        np.testing.assert_almost_equal(f, fhat, 10)

    def test_int_grad3d(self):
        # Compute 3D grad of a field, integrate and compare the result with original field
        xp = np.linspace(0,1,4)
        yp = np.linspace(0,1,5)
        zp = np.linspace(0,1,3)
        [x,y,z] = np.meshgrid(xp,yp,zp);
        f = np.exp(x+y+z) + np.sin((x-2*y+3*z)*3);
        [fy,fx,fz]=np.gradient(f,yp,xp,zp)
        fhat = intgrad((fx,fy,fz),(xp,yp,zp), constant = 1)
        np.testing.assert_almost_equal(f, fhat, 10)

if __name__ == '__main__':
    unittest.main()
