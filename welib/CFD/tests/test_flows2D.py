import unittest
import numpy as np    
import os as os
from welib.CFD.flows2D import *

scriptDir=os.path.dirname(__file__)

class Test(unittest.TestCase):
    def test_flowrate(self):
        theta = np.linspace(0, 2*np.pi, 250, endpoint=True)
        xc = np.cos(theta)
        yc = np.sin(theta)
        uc = xc  # radial field: u = x
        vc = yc  #               v = y
        Q = flowrate2D(xc, yc, uc, vc, verbose=False, ns=-1)
        np.testing.assert_almost_equal(Q, 2*np.pi, 3)

if __name__ == '__main__':
    unittest.main()
