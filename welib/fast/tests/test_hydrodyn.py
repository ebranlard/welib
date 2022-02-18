# --- Common libraries 
import os
import unittest
import numpy as np
from welib.fast.hydrodyn import *

MyDir=os.path.dirname(__file__)

class Test(unittest.TestCase):
    """ See examples/ for more examples """

    def test_hydro_volume(self):
        # --- Read a HydroDyn model
        filename = os.path.join(MyDir,'../../../data/Spar/Spar_HD.dat')
        WtrDens=1025
        WtrDpth= 320
        Gravity=9.8
        hd = HydroDyn(filename)
        hd.init(Gravity=Gravity, WtrDens=WtrDens, WtrDpth=WtrDpth)
        #hd.writeSummary(filename.replace('.dat','.HD_python.sum'))
        # --- Structure volume
        VMorison = [e.MorisonData['Vinner']     for e in hd.morison.graph.Elements]
        VNoDiv   = [hd.memberVolumeStructure(e) for e in hd.graph.Elements]
        VtotMorison = hd.VolumeStructure(method='Morison')
        VtotNoDiv   = hd.VolumeStructure(method='NoDiv')
        np.testing.assert_almost_equal(VMorison,VNoDiv)
        np.testing.assert_almost_equal(VtotMorison,VtotNoDiv)
        # --- Submerged volume
        VMorison = [e.MorisonData['Vsubmerged'] for e in hd.morison.graph.Elements]
        VNoDiv   = [hd.memberVolumeSubmerged(e) for e in hd.graph.Elements]
        VtotMorison = hd.VolumeSubmerged(method='Morison')
        VtotNoDiv   = hd.VolumeSubmerged(method='NoDiv')
        np.testing.assert_almost_equal(VMorison,VNoDiv)
        np.testing.assert_almost_equal(VtotMorison,VtotNoDiv)

if __name__ == '__main__':
    unittest.main()
