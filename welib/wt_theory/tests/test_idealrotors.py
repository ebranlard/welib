import unittest
import numpy as np    
from welib.wt_theory.idealrotors import *

class TestIdealRotors(unittest.TestCase):
    def test_CP_a_opt(self):
        # AD with wake rotation
        lambda_ = np.linspace(0.01,10,10)

        a1, ap1 = ADMTO_inductions(lambda_, method='fzero')
        a2, ap2 = ADMTO_inductions(lambda_, method='analytical')
        a3, ap3 = ADMTO_inductions(lambda_, method='analytical2')

        np.testing.assert_almost_equal(a1, a2, 10)
        np.testing.assert_almost_equal(a2, a3, 10)
        np.testing.assert_almost_equal(ap2, ap3, 10)

        CP1, a1 = ADMTO_CP(lambda_, method='num_int')
        CP2, a2 = ADMTO_CP(lambda_, method='analytical')
        np.testing.assert_almost_equal(CP1, CP2, 10)


if __name__ == '__main__':
    unittest.main()
