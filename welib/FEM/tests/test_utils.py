import unittest
import os
import numpy as np

from welib.FEM.utils import *

MyDir=os.path.dirname(__file__)

# --------------------------------------------------------------------------------}
# --- TESTS
# --------------------------------------------------------------------------------{
class Test(unittest.TestCase):

    def test_rigid_transformation(self):
        # Transfer loads known at a source point to another point
        Ps = (2,0,0)
        Pd = (0,0,0)
        T = rigidTransformationTwoPoints_Loads(Ps, Pd)

        fs = (4,1,5,50,0,10)
        fd = T.dot(fs)
        np.testing.assert_almost_equal(fd, (4,1,5,50,-10,12), 5)

        Tds = rigidTransformationTwoPoints_Loads(Pd, Ps)
        fs2 = Tds.dot(fd)
        np.testing.assert_almost_equal(fs2, fs, 5)


if __name__=='__main__':
    unittest.main()
