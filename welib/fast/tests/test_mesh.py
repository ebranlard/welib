# --- Common libraries 
import os
import unittest
import numpy as np
from welib.fast.fast_mesh import *

MyDir=os.path.dirname(__file__)

class TestMesh(unittest.TestCase):
    """ See examples/ for more examples """
    def test_rigid_body_motion(self):
        # Test loads for given displacements/velocities/accelerations
        pm = PointMesh(2)
        pm.Position[0,:] = (1,-1,0)
        pm.Position[1,:] = (1,2,-10)
        pm.Connectivity=np.array([[0,1]])
        # --- 
        time =np.linspace(0,10,100)
        q       = (10  , 20  , -2    , 0.1  , 0.2  , 0.01  )
        qd      = (10  , 20  , -2    , 0.1  , 0.2  , 0.01  )
        qdd     = (10  , 20  , -2    , 0.1  , 0.2  , 0.01  )
        pm.rigidBodyMotion(q=q, qd=qd, qdd=qdd, RefPoint=(0,0,0))

        np.testing.assert_almost_equal(pm.TranslationDisp[0,:] ,[ 9.9807925 , 20.02426565, -2.29323805] )
        np.testing.assert_almost_equal(pm.TranslationVel [0,:] ,[ 9.95110973, 20.03913173, -2.29373194] )
        np.testing.assert_almost_equal(pm.TranslationAcc [0,:] ,[ 9.89197203, 20.06801602, -2.28004071] )
        np.testing.assert_almost_equal(pm.RotationVel    [0,:] ,[0.1,  0.2,  0.01] )
        np.testing.assert_almost_equal(pm.RotationAcc    [0,:] ,[0.1,  0.2,  0.01] )


if __name__ == '__main__':
    unittest.main()
