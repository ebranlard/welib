# --- Common libraries 
import os
import unittest
import numpy as np
from welib.fast.elastodyn import *

MyDir=os.path.dirname(__file__)


class TestED(unittest.TestCase):
    """ See examples/ for more examples """

    def test_ED_blade_params(self):
        EDfilename=os.path.join(MyDir,'../../../data/NREL5MW/onshore/NREL5MW_ED_Onshore.dat')
        p = bladeParameters(EDfilename)
        # Physical quantities / Inertias
        np.testing.assert_almost_equal(p['BldMass'], 16775.66975907849)
        np.testing.assert_almost_equal(p['FirstMom'], 346419.0832041095)
        np.testing.assert_almost_equal(p['SecondMom'], 11242397.13082460)
        np.testing.assert_almost_equal(p['BldCG'], 20.65008957491178)

        # Shape functions
        np.testing.assert_almost_equal(p['ShapeF1_full'][1:3], np.array([9.537890879514419e-05 ,  1.496964261332097e-03]))
        np.testing.assert_almost_equal(p['ShapeF1_full'][-3:-1], np.array([7.907408392607156e-01 , 9.302315172663902e-01]))
        np.testing.assert_almost_equal(p['ddShapeF1_full'][1:3], np.array([ 1.051074892093353e-04,   2.102754797663044e-04]))
        np.testing.assert_almost_equal(p['ddShapeF1_full'][-3:-1], np.array([ 4.957076482520056e-04,  -9.910176341969626e-05]))

        # Generalized quantities
        np.testing.assert_almost_equal(p['MBF'][0,0], 884.8557898849025)
        np.testing.assert_almost_equal(p['MBF'][1,1], 537.0972116447014)
        np.testing.assert_almost_equal(p['KBF'],   np.array([[   1.669817537559324e+04 , -1.453663035386739e+03], [  -1.453663035386742e+03  , 8.495748847887896e+04]]))
        np.testing.assert_almost_equal(p['MBE'][0,0], 1368.695198467347)
        np.testing.assert_almost_equal(p['KBE'][0,0], 66822.08061735689)
        np.testing.assert_almost_equal(p['KBFCent'][1,1], 3.144387706952399e+03)

        #---  Twisted shape funciton
        # x/y, BF1/BF2/BE, node, deriv
        # - Flap 1
        # Root, all deriv
        np.testing.assert_almost_equal(p['TwistedSF'][0,0,0,0],0)
        np.testing.assert_almost_equal(p['TwistedSF'][0,0,0,1],0)
        np.testing.assert_almost_equal(p['TwistedSF'][0,0,0,2],3.200725883815983e-05)
        # Tip, all deriv
        np.testing.assert_almost_equal(p['TwistedSF'][0,0,-1,0], 9.910328208697139e-01)
        np.testing.assert_almost_equal(p['TwistedSF'][0,0,-1,1], 3.815332774490886e-02)
        np.testing.assert_almost_equal(p['TwistedSF'][0,0,-1,2],-4.940445502016021e-04)
        # Mid, all deriv
        np.testing.assert_almost_equal(p['TwistedSF'][0,0,9,0],1.400789659935833e-01)
        np.testing.assert_almost_equal(p['TwistedSF'][0,0,9,1],1.260777581991000e-02)
        np.testing.assert_almost_equal(p['TwistedSF'][0,0,9,2],8.191781538881797e-04)
        # - Edge 1
        # Root, all deriv
        np.testing.assert_almost_equal(p['TwistedSF'][0,2,0,0],                     0)
        np.testing.assert_almost_equal(p['TwistedSF'][0,2,0,1],                     0)
        np.testing.assert_almost_equal(p['TwistedSF'][0,2,0,2], 4.414743766275156e-05)
        # Tip, all deriv
        np.testing.assert_almost_equal(p['TwistedSF'][0,2,-1,0],   0.167501776715005)
        np.testing.assert_almost_equal(p['TwistedSF'][0,2,-1,1],  0.003946589036490)
        np.testing.assert_almost_equal(p['TwistedSF'][0,2,-1,2],                  0)
        # Mid, all deriv
        np.testing.assert_almost_equal(p['TwistedSF'][0,2,9,0],5.141515067563961e-02)
        np.testing.assert_almost_equal(p['TwistedSF'][0,2,9,1],3.261536229664988e-03)
        np.testing.assert_almost_equal(p['TwistedSF'][0,2,9,2],6.728816969020595e-05)

        # --- AxRedBld
        #3, 3, ED['BldNodes']+2))
        np.testing.assert_almost_equal(p['AxRedBld'][:,0,-1], [2.722873037538838e-02,  4.273914168364583e-02,  8.556785740654360e-04])
        np.testing.assert_almost_equal(p['AxRedBld'][:,0,0],[0,0,0]) 
        np.testing.assert_almost_equal(p['AxRedBld'][:,0,9],[1.100604990701515e-03, -1.686447688419782e-03 , 1.552531491358490e-05])
        np.testing.assert_almost_equal(p['AxRedBld'][:,2,-1],[ 8.556785740654360e-04 , 4.524659962852886e-03 , 2.160759559741185e-02])
        np.testing.assert_almost_equal(p['AxRedBld'][:,2,-2],[ 7.654290347321518e-04 , 3.949211069010193e-03 , 2.018155508430989e-02])
        np.testing.assert_almost_equal(p['AxRedBld'][:,2,9],[1.552531491358490e-05 , 2.331092957799932e-05 , 2.931333434244318e-03])

        # --- Frequencies
        np.testing.assert_almost_equal(p['FreqBF'][0,2],0.723649838650248)
        np.testing.assert_almost_equal(p['FreqBF'][1,2],2.038385218375959)
        np.testing.assert_almost_equal(p['FreqBE'][0,:], np.array([1.112056253877955, 1.112056253877955, 1.128488446907760]))
        np.testing.assert_almost_equal(p['CBE']   [0,0],91.32399596067735)


    def test_ED_tower_params(self):
        EDfilename=os.path.join(MyDir,'../../../data/NREL5MW/onshore/NREL5MW_ED_Onshore.dat')
        RotMass = 107107.00927723547
        Gravity = 9.80665
        p = towerParameters(EDfilename, RotMass=RotMass, gravity=Gravity)

        # Physical quantities / Inertias
        np.testing.assert_almost_equal(p['TwrMass'], 347460.2316000000)
        np.testing.assert_almost_equal(p['TwrTpMass'], 347107.0092772355)

        # Shape functions nModes x nNodes x nDeriv
        np.testing.assert_almost_equal(p['TwrFASF'][0,-3:, 0],  [0.871714585606445  ,0.958488717049806  ,1.000000000000000])
        np.testing.assert_almost_equal(p['TwrSSSF'][1,:4, 1], [0, -0.065390448129370,-0.175841634092733,-0.263419614469624])

        # Generalized quantities
        np.testing.assert_almost_equal(np.diag(p['MTFA'])/1e6,np.array([4.010392552595473e+05, 2.755960009968935e+07])/1e6)
        np.testing.assert_almost_equal(np.diag(p['MTSS'])/1e6,np.array([4.008128638209492e+05, 3.558329856444363e+07])/1e6)
        np.testing.assert_almost_equal(       (p['KTSS'])/1e6, np.array([[1.620074333887758e+06, 4.338323386161570e+06],[2.943153808092614e+07, 1.125981782739319e+10]])/1e6)
        np.testing.assert_almost_equal(       (p['KTFA'])/1e6, np.array([[1.722040410360708e+06,-7.442658427757327e+06],[2.400069959451868e+07, 9.087764067630390e+09]])/1e6)
        np.testing.assert_almost_equal(np.diag(p['KTFAGrav']),[ -1.066694734234754e+03 , -1.169976720444022e+06])
        np.testing.assert_almost_equal(np.diag(p['KTSSGravTT']), [-3.763420370177263e-02,-8.551393808571898e+01])
        # 
        np.testing.assert_almost_equal(p['AxRedTFA'][0,0,-3:],[1.158577163756421e-02,1.332102574930896e-02,1.417295499349627e-02])

        # Frequencies
        np.testing.assert_almost_equal(p['FreqTSS'], np.array([[0.874131532646607, 0.306008233418085],[2.845059363346653 ,2.792516489887383]]))
        np.testing.assert_almost_equal(p['CTFA'], np.array([[ 6.095030992638342e+03,-8.145424597927582e+03],[8.494864986527940e+04, 9.945867823326178e+06]]))




if __name__ == '__main__':
    unittest.main()
