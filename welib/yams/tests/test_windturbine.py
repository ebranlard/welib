import unittest
import os
import numpy as np
from welib.yams.windturbine import *
from welib.yams.utils import *
import welib.weio as weio

MyDir=os.path.dirname(__file__)

# --------------------------------------------------------------------------------}
# --- TESTS
# --------------------------------------------------------------------------------{
class TestWindturbclassmethod(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # Read FAST structural model
        cls.WT = FASTWindTurbine(os.path.join(MyDir,'./../../../data/Spar/Main_Spar_ED.fst'))

    def test_fnd_ED(self):
        # --- Floater
        fnd=self.WT.fnd
        # Values from SubDyn summary file
        Mfnd_ref =  np.array([[   7.500000E+06,   0.000000E+00,   0.000000E+00,   0.000000E+00,  -6.749999E+08,  0.000000E+00],
                              [   0.000000E+00,   7.500000E+06,   0.000000E+00,   6.749999E+08,   0.000000E+00,   0.000000E+00],
                              [   0.000000E+00,   0.000000E+00,   7.500000E+06,   0.000000E+00,   0.000000E+00,   0.000000E+00],
                              [   0.000000E+00,   6.749999E+08,   0.000000E+00,   6.494974E+10,   0.000000E+00,   0.000000E+00],
                              [  -6.749999E+08,   0.000000E+00,   0.000000E+00,   0.000000E+00,   6.494974E+10,   0.000000E+00],
                              [   0.000000E+00,   0.000000E+00,   0.000000E+00,   0.000000E+00,   0.000000E+00,   1.600015E+08]])

        np.testing.assert_allclose(fnd.mass_matrix_at([0,0,-20]),Mfnd_ref, 1e-5)
        np.testing.assert_allclose(fnd.mass, 7.500000E+06)
        np.testing.assert_allclose(fnd.masscenter_pos_global, [0,0,-90])

    def test_twr(self):
        # --- Tower
        twr=self.WT.twr
        np.testing.assert_allclose(twr.pos_global, [0,0,20])
        np.testing.assert_allclose(twr.mass,        217537.844)
        np.testing.assert_allclose(twr.masscenter, [0,0, 28.9555])
        np.testing.assert_allclose(twr.length,      67.6)
        np.testing.assert_allclose(np.diag(twr.inertia), [2.6183e8,2.6183e8, 0.0], 1e-3)
        # NOTE: this is offshore tower
        MM_ref=[[ 2.17538e+05, 0.00000e+00, 0.00000e+00, 0.00000e+00, 6.29892e+06,-0.00000e+00, 5.66449e+04, 8.82302e+05, 0.00000e+00, 0.00000e+00],
                [ 0.00000e+00, 2.17538e+05, 0.00000e+00,-6.29892e+06, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 5.56622e+04, 1.19676e+06],
                [ 0.00000e+00, 0.00000e+00, 2.17538e+05, 0.00000e+00,-0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00],
                [ 0.00000e+00,-6.29892e+06, 0.00000e+00, 2.61831e+08,-0.00000e+00,-0.00000e+00, 0.00000e+00, 0.00000e+00,-2.70207e+06,-4.60419e+07],
                [ 6.29892e+06, 0.00000e+00,-0.00000e+00,-0.00000e+00, 2.61831e+08,-0.00000e+00, 2.74189e+06, 3.43003e+07, 0.00000e+00, 0.00000e+00],
                [-0.00000e+00, 0.00000e+00, 0.00000e+00,-0.00000e+00,-0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00],
                [ 5.66449e+04, 0.00000e+00, 0.00000e+00, 0.00000e+00, 2.74189e+06, 0.00000e+00, 3.16332e+04, 3.23211e+05, 0.00000e+00, 0.00000e+00],
                [ 8.82302e+05, 0.00000e+00, 0.00000e+00, 0.00000e+00, 3.43003e+07, 0.00000e+00, 3.23211e+05, 5.38910e+06, 0.00000e+00, 0.00000e+00],
                [ 0.00000e+00, 5.56622e+04, 0.00000e+00,-2.70207e+06, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 3.08678e+04, 4.21109e+05],
                [ 0.00000e+00, 1.19676e+06, 0.00000e+00,-4.60419e+07, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 4.21109e+05, 9.89031e+06]]
        np.testing.assert_allclose(twr.mass_matrix,MM_ref, 1e-5)

        # convert to rigid body, and check that properties still match
        twr_rigid=twr.toRigidBody()
        for p in['mass','pos_global','masscenter','masscenter_pos_global', 'inertia']:
            np.testing.assert_allclose(getattr(twr,p), getattr(twr_rigid,p))

    def test_bld(self):
        # --- Blade
        pass

if __name__=='__main__':
    np.set_printoptions(linewidth=300, precision=5)
    unittest.main()
