# --- Common libraries 
import os
import unittest
import numpy as np
from welib.fast.hydrodyn import *

MyDir=os.path.dirname(__file__)


def getHDSpar():
    filename = os.path.join(MyDir,'../../../data/Spar/Spar_HD.dat')
    WtrDens=1025
    WtrDpth= 320
    Gravity=9.8
    hd = HydroDyn(filename)
    u, y = hd.init(Gravity=Gravity, WtrDens=WtrDens, WtrDpth=WtrDpth)
    #hd.writeSummary(filename.replace('.dat','.HD_python.sum'))
    return hd, u, y


class TestSpar(unittest.TestCase):
    """ See examples/ for more examples """

    #@classmethod
    #def setUpClass(cls):
    #    cls.hd, cls.u, cls.y = getHDSpar()
    #    return cls

    def test_hydro_volume(self):
        # Test Volume as computed by advanced "Morison" routines
        #  or simplified routines (provided in HydroDyn class)  
        #  without element division
        # --- Read a HydroDyn model
        hd, u, y = getHDSpar()
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



    def test_hydro_loads(self):
        # Test loads for given displacements/velocities/accelerations
        hd, u, y = getHDSpar()
        umesh = u['Morison']['Mesh']
        ymesh = y['Morison']['Mesh']
        q    = (10  , 20  , -2    , 0.1  , 0.2  , 0.01  )
        qd   = (10  , 20  , -2    , 0.1  , 0.2  , 0.01  )
        qdd  = (10  , 20  , -2    , 0.1  , 0.2  , 0.01  )
        umesh.rigidBodyMotion(q=q, qd=qd, qdd=qdd)
        y=hd.calcOutput(t=0, u=u, y=y)
        Fh, Mh = ymesh.mapLoadsToPoint((q[0],q[1],q[2]))
        #print('F_Hydro {:16.3f}{:16.3f}{:16.3f}'.format(*Fh))
        #print('M_Hydro {:16.3f}{:16.3f}{:16.3f}'.format(*Mh))
        np.testing.assert_almost_equal(Fh/1e6, np.array([42733594.104, -441445028.029,  80283363.828])/1e6 )
        np.testing.assert_almost_equal(Mh/1e6, np.array([-29284711118.606, -5247122545.815,  5351864809.826])/1e6)



if __name__ == '__main__':
    #TestSpar().test_hydro_loads()
    unittest.main()
