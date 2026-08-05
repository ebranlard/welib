
import unittest
import numpy as np
import os

from welib.yams.models.MTNSB import FASTmodel2MTNSB
from welib.tools.strings import printVec
from welib.tools.strings import printMat

MyDir=os.path.dirname(__file__)

class TestMTNSB(unittest.TestCase):
    def test_MTNSB_FAST_SubDyn_CM(self):
        # SubDyn with Concentrated inertias
        """
        End-to-end example for SubDyn concentrated inertias in YAMS recursive model.
        This example builds a MTNSB model from a FAST file containing SubDyn concentrated
        masses/inertias, performs an eigenanalysis, then disables the concentrated inertias
        on the already-built substructure beam and recomputes the eigenanalysis.

        """
        fstFile = os.path.abspath(os.path.join(MyDir, '../../../data/Monopile/Main_MT100_EVA_CM.fst'))

        # --- Full model with concentrated inertias
        WT = FASTmodel2MTNSB(
            fstFile,
            shapes_sub=[0, 4],
            shapes_twr=[0, 1],
            shapes_bld=[],
            DEBUG=False,
            bStiffening=True,
            main_axis='z',
            fixedShaft=True,
            algo='OpenFAST',
            FEM_method='full',
        ).WT

        q = np.zeros((len(WT.q0), 1))
        WT.grd.setDOF(q)
        freq_d_cm, zeta_cm, Q_cm, freq_0_cm = WT.grd.eva()
        M_cm = WT.grd.M.copy()

        cm_saved = list(getattr(WT.fnd, 'concentrated_inertias', []))

        # --- Removing concentrated inertias
        WT.fnd.concentrated_inertias = []
        WT.fnd.computeMassMatrix(s_G=WT.fnd.s_G0, inPlace=True)
        WT.grd.setDOF(q)
        freq_d_no, zeta_no, Q_no, freq_0_no = WT.grd.eva()
        M_no = WT.grd.M.copy()

        # --- Comparison
        dM00 = M_cm[0, 0] - M_no[0, 0]
        dF = freq_0_cm[:4] - freq_0_no[:4]

        test=False
        if not test:
            print('FAST file:', fstFile)
            print('SubDyn concentrated inertias count:', len(cm_saved))
            print('Global mass matrix shape:', M_cm.shape)
            print('M00 with CM   :', M_cm[0, 0])
            print('M00 without CM:', M_no[0, 0])
            print('dM00          :', dM00)
            printVec('f0 with CM [Hz]', freq_0_cm[:4])
            printVec('f0 no CM  [Hz]', freq_0_no[:4])
            printVec('df0 [Hz]', dF)

        out = {
            'n_cm': len(cm_saved),
            'dM00': dM00,
            'freq_0_cm': freq_0_cm,
            'freq_0_no': freq_0_no,
            'freq_d_cm': freq_d_cm,
            'freq_d_no': freq_d_no,
            'M_cm': M_cm,
            'M_no': M_no,
        }

        # Regression checks:
        #  - The SubDyn input used here contains one concentrated inertia.
        #  - Removing it must reduce the assembled M00 by approximately 1e5 kg.
        #  - The first two frequencies should decrease when CM are included.
        np.testing.assert_equal(out['n_cm'], 1)
        np.testing.assert_allclose(out['dM00'], 1.0e5, rtol=0, atol=1e-6)
        np.testing.assert_array_less(out['freq_0_cm'][:2], out['freq_0_no'][:2])

        return out


if __name__=='__main__':
    TestMTNSB().test_MTNSB_FAST_SubDyn_CM()
#     unittest.main()
