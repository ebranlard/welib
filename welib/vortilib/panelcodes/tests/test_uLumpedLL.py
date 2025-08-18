
import unittest
import numpy as np

from welib.vortilib.panelcodes.uLumpedLL import flat_plate_lin_acc

class TestFlatPlateAcc(unittest.TestCase):
    def test_out_values(self):
        df, wake_out = flat_plate_lin_acc(nstep=200, alpha_deg=5.0, dt_u0_over_c=0.25, U0=50.0, c=1.0, rho=1.0, dxw_factor=0.3)
        ref_clt   = [1.9696, 0.91240, 0.98932]
        ref_Gammat = [0.393939, 0.497633, 0.98909]
        idxs      = [1, 2, -1]
        np.testing.assert_almost_equal(df['Cl_rel'].to_numpy()[idxs],   ref_clt,   decimal=4)
        np.testing.assert_almost_equal(df['Gamma_rel'].to_numpy()[idxs], ref_Gammat, decimal=4)


if __name__ == "__main__":
    unittest.main()