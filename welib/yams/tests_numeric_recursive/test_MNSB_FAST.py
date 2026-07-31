import unittest
import numpy as np
import os

from welib.yams.models.MNSB_FAST import FASTmodel2MNSB
from welib.tools.strings import printMat

MyDir=os.path.dirname(__file__)

class TestMNSB(unittest.TestCase):
    def test_MNSB_FAST(self):

        nDOF=2
        q = np.zeros((nDOF,1)) # TODO, full account of q not done
        q[[0]]=3
        q[[1]]=0.1

        np.set_printoptions(linewidth=500)
        fstFile = os.path.join(MyDir, '../../../data/Monopile/Main_MT100_JONSWAP.fst')
        shapes_sub=[0,4]

        # --- Auto assembly with z axis
        assembly='auto'
        WTA = FASTmodel2MNSB(fstFile, q=q, shapes_sub=shapes_sub, shapes_bld=[], DEBUG=False, bStiffening=True, main_axis='z', assembly=assembly, fixedShaft=True).WT
        # --- Manual assembly with x axis
        assembly='manual'
        WTM = FASTmodel2MNSB(fstFile, shapes_sub=shapes_sub, shapes_bld=[], DEBUG=False, bStiffening=True, main_axis='z', assembly=assembly, fixedShaft=True).WT
        #WTM = WTA.copy()
        #WTM.MM=0
        #WTM.DD=0
        #WTM.KK=0
        #WTM.manual_assembly(q=q, DEBUG=False, fixedShaft=True)

        # --- "Tower"
        MMref=np.array([
           [334285.748379, -4.713536e+06],
           [-4.713536e+06,  8.571429e+07]])
        KKref=np.array([
           [ 2.395185e+07, -1.199358e+09],
           [-1.199358e+09,  7.998260e+10]])
        DDref=np.array([
           [143711.078766, -7.196145e+06],
           [-7.196145e+06,  4.798956e+08]])
        #printMat(WTA.twr.MM[6:,6:]  , var='MM', digits=6)
        #printMat(WTA.twr.KK[6:,6:]  , var='KK', digits=6)
        #printMat(WTA.twr.DD[6:,6:]  , var='DD', digits=6)
        np.testing.assert_almost_equal(WTA.twr.DD    ,WTM.twr.DD)
        np.testing.assert_almost_equal(WTA.twr.KK    ,WTM.twr.KK)
        np.testing.assert_almost_equal(WTA.twr.MM/1e5,WTM.twr.MM/1e5 , 5)

        np.testing.assert_almost_equal(WTA.twr.DD[6:,6:]/1e6,DDref/1e6, 4)
        np.testing.assert_almost_equal(WTA.twr.KK[6:,6:]/1e8,KKref/1e8, 4)
        np.testing.assert_almost_equal(WTA.twr.MM[6:,6:]/1e6,MMref/1e6, 5)

        np.testing.assert_almost_equal(WTA.twr.pos_global.ravel(), (0,0,-50))
        np.testing.assert_almost_equal(WTA.alpha.ravel(), (0,0.1,0))
        # --- nacelle
        #print(WTM.nac)
        #nac_MMref= np.array([[ 240000.,       0.,       0. ,      0.,  420000.,      -0.],
        #                     [      0.,  240000.,       0. ,-420000.,       0.,  456000.],
        #                     [      0.,       0.,  240000. ,      0., -456000.,       0.],
        #                     [      0., -420000.,       0. ,      0.,       0.,       0.],
        #                     [ 420000.,       0., -456000. ,      0.,       0.,       0.],
        #                     [     -0.,  456000.,       0. ,      0.,       0., 2607890.]])
        np.testing.assert_almost_equal(WTA.nac.mass, 0)
        np.testing.assert_almost_equal(WTA.nac.pos_global.ravel(),(3,0,51.0)) # TODO Influence of rotation should actually be felt due to rigid tower of 1m
        #np.testing.assert_almost_equal(WTA.nac.MM,nac_MMref)


        # --------------------------------------------------------------------------------}
        # --- Full system
        # --------------------------------------------------------------------------------{
        #MMref=np.array([
        #   [334285.748379, -4.713536e+06],
        #   [-4.713536e+06,  8.571429e+07]])
        #KKref=np.array([
        #   [ 2.395185e+07, -1.199358e+09],
        #   [-1.199358e+09,  7.998260e+10]])
        #DDref=np.array([
        #   [143711.078766, -7.196145e+06],
        #   [-7.196145e+06,  4.798956e+08]])
        #printMat(WTA.MM  , var='MM', digits=5)
        #printMat(WTA.KK  , var='KK', digits=5)
        #printMat(WTA.DD  , var='DD', digits=5)

        np.testing.assert_almost_equal(WTA.DD,WTM.DD)
        np.testing.assert_almost_equal(WTA.KK,WTM.KK)
        np.testing.assert_almost_equal(WTA.MM/1e5,WTM.MM/1e5 , 5)

        np.testing.assert_almost_equal(WTA.DD/1e6,DDref/1e6, 4)
        np.testing.assert_almost_equal(WTA.KK/1e8,KKref/1e8, 4)
        np.testing.assert_almost_equal(WTA.MM/1e6,MMref/1e6, 5)

if __name__=='__main__':
    unittest.main()
