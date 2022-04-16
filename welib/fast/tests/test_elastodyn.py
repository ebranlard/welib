# --- Common libraries 
import os
import unittest
import numpy as np
from welib.fast.elastodyn import *
from welib.yams.utils import skew
from welib.yams.flexibility import GMBeam, GKBeam, GKBeamStiffnening

MyDir=os.path.dirname(__file__)


class TestED(unittest.TestCase):
    """ See examples/ for more examples """

    def test_ED_rot_params(self):
        Gravity = 9.80665
        EDfilename=os.path.join(MyDir,'../../../data/NREL5MW/onshore/NREL5MW_ED_Onshore.dat')
        #p = bladeParameters(EDfilename, AdjBlMs=1)  # <<<<<<
        prot,pbld,phub = rotorParameters(EDfilename)

        # --- ElastoDyn Summary file
        #  Rotor Mass            (kg)       109389.842
        #  Rotor Inertia         (kg-m^2) 38677040.613
        #  Mass                  (kg)        17536.614    17536.614    17536.614
        #  Second Mass Moment    (kg-m^2) 11752352.265 11752352.265 11752352.265
        #  First Mass Moment     (kg-m)     362132.653   362132.653   362132.653
        #  Center of Mass        (m)            20.650       20.650       20.650
        np.testing.assert_almost_equal(prot['RotMass'],        109389.842 , 3)
        np.testing.assert_almost_equal(prot['RotIner'],      38677040.613 , 3)
        np.testing.assert_almost_equal(pbld[0]['BldMass'],      17536.614 , 3)
        np.testing.assert_almost_equal(pbld[0]['SecondMom'], 11752352.265 , 3)
        np.testing.assert_almost_equal(pbld[0]['FirstMom'] ,   362132.653 , 3)
        np.testing.assert_almost_equal(pbld[0]['BldCG'] ,           20.650, 3)




    def test_ED_blade_params_NoBldAdj(self):
        # 
        # NOTE:  These parameters were obtained for AdjBlMs=1!!!
        # Kept to avoid redoing all these tests..
        Gravity = 9.80665
        EDfilename=os.path.join(MyDir,'../../../data/NREL5MW/onshore/NREL5MW_ED_Onshore.dat')
        p = bladeParameters(EDfilename, AdjBlMs=1)  # <<<<<<
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

	
        # --------------------------------------------------------------------------------}
        # --- Mass Matrix using GM or not
        # --------------------------------------------------------------------------------{
#         nq=2
#         nNodes=n+2
#         p['Ut'] = np.zeros((nq, 3, nNodes))
#         p['Vt'] = np.zeros((nq, 3, nNodes))
#         p['Kt'] = np.zeros((nq, 3, nNodes))
#         p['U']  = np.zeros((nq, 3, nNodes))
#         p['V']  = np.zeros((nq, 3, nNodes))
#         p['K']  = np.zeros((nq, 3, nNodes))
#         for j,jj,idir,name in zip(range(0,nq), (0,2), (0,1), ('F1','E1')): # direction is x, x, y
#             p['Ut'][j][0,:] = p['TwistedSF'][0, jj, :, 0]  # x
#             p['Ut'][j][1,:] = p['TwistedSF'][1, jj, :, 0]  # y
#             p['Vt'][j][0,:] = p['TwistedSF'][0, jj, :, 1]  # x
#             p['Vt'][j][1,:] = p['TwistedSF'][1, jj, :, 1]  # y
#             p['Kt'][j][0,:] = p['TwistedSF'][0, jj, :, 2]  # x
#             p['Kt'][j][1,:] = p['TwistedSF'][1, jj, :, 2]  # y
#             p['U'][j][idir,:]  = p['Shape'+name+'_full']  
#             p['V'][j][idir,:]  = p['dShape'+name+'_full'] 
#             p['K'][j][idir,:]  = p['ddShape'+name+'_full']


        # --- Calling GM/GK Beam with OpenFAST method
        inertiaAtBladeRoot=True # TODO for loop around that
        rh = 0 # Hub Radius # TODO make this an option if from blade root or not
        s_G0 = np.zeros((3, len(p['s_span'])))
        s_G0[2,:] = p['s_span'] + rh 
        MM, IT = GMBeam(s_G0, p['s_span'], p['m_full'], p['Ut'], rot_terms=True, method='OpenFAST', main_axis='z', U_untwisted=p['U'], M1=True) 
        Gr, Ge, Oe, Oe6 = IT['Gr'], IT['Ge'], IT['Oe'], IT['Oe6']

        KK = GKBeam(p['s_span'], p['EI'], p['K'], bOrth=False, method='OpenFAST')

        KKg_SW = GKBeamStiffnening(p['s_span'], p['V'], Gravity, p['m_full'], Mtop=0, Omega=0, bSelfWeight=True,  bMtop=False, bRot=False, main_axis='z', method='OpenFAST')
        #KKg_TM = GKBeamStiffnening(p['s_span'], p['V'], Gravity, p['m_full'], Mtop=0, Omega=0, bSelfWeight=False, bMtop=True, bRot=False, main_axis='z')
        KKg_Om = GKBeamStiffnening(p['s_span'], p['V'], Gravity, p['m_full'], Mtop=0, Omega=1, bSelfWeight=False, bMtop=False, bRot=True, main_axis='z', method='OpenFAST')

        # --- Call bladeDerivedParameters for "manual" calculation
        p = bladeDerivedParameters(p, inertiaAtBladeRoot=inertiaAtBladeRoot)
#         print('KK\n',KK[6:,6:])
#         print('KKe\n',p['Ke'])
#         print('KKg_SW\n',KKg_SW[6:,6:])
#         print('KKg_Om\n',KKg_Om[6:,6:])
#         print('KKg_Om\n',p['Kg_Om'])

        # --- TODO TODO TODO TODOs
        # - compute general centrifugal stiffening tersm in GM/GKBeam
        # - Update SID/ *parameters functions to use GM/GK beam

#         print('MM',MM[0,0])
#         print('Ms',p['BldMass'])
#         print('J\n',p['J'])
#         print('J\n',MM[3:6,3:6])
#         print('mdCM_GM\n',MM[3:6,0:3])
#         print('mdCM_OF\n',-skew(p['mdCM']))
#         print('me_GM\n',MM[6:,6:])
#         print('me_OF\n',p['Me'])
#         print('Ct_GM\n',MM[0:3,6:])
#         print('Ct_OF\n',p['Ct'].T)
#         print('Cr_GM\n',MM[3:6,6:])
#         print('Cr_OF\n',p['Cr'].T)
#         print('OeM1_GM\n',IT['Oe6M1'])
#         print('OeM1_OF\n',p['OeM1'])
#         for j in np.arange(nq):
#             print('')
#             print('Oe6 GM{}\n'.format(j), Oe6[j])
#             print('Oe6 OF{}\n'.format(j), p['Oe6'][j])
        #for j in np.arange(nq):
        #    print('')
        #    print('mdCM1_GM {}\n'.format(j),IT['mdCM1'][:,j])
        #    print('mdCM1_OF {}\n'.format(j), p['mdCM1'][:,j])
        # --- Compare both "manual" and GM/GKBeam approach
        np.testing.assert_almost_equal(MM[0,0,]           , p['BldMass'])
        np.testing.assert_almost_equal(MM[3:6,3:6]       ,  p['J'])
        np.testing.assert_almost_equal(MM[3:6,0:3]       ,  skew(p['mdCM']))
        np.testing.assert_almost_equal(np.diag(MM[6:,6:]),  np.diag(p['Me']))
        np.testing.assert_almost_equal(MM[0:3,6:],          p['Ct'].T)
        np.testing.assert_almost_equal(MM[3:6,6:],          p['Cr'].T)
        np.testing.assert_almost_equal(IT['Oe6_M1'] , p['Oe_M1'])
        np.testing.assert_almost_equal(IT['mdCM_M1'], p['mdCM_M1'])
        np.testing.assert_almost_equal(KK[6:,6:]  , p['Ke'])


        # --------------------------------------------------------------------------------}
        # --- Using inertia at rotor center
        # --------------------------------------------------------------------------------{
        # --- Calling GM Beam with OpenFAST method
        inertiaAtBladeRoot=False
        rh = p['HubRad'] # Hub Radius # TODO make this an option if from blade root or not
        s_G0 = np.zeros((3, len(p['s_span'])))
        s_G0[2,:] = p['s_span'] + rh 
        MM, IT = GMBeam(s_G0, p['s_span'], p['m_full'], p['Ut'], rot_terms=True, method='OpenFAST', main_axis='z', U_untwisted=p['U'], M1=True) 
        Gr, Ge, Oe, Oe6 = IT['Gr'], IT['Ge'], IT['Oe'], IT['Oe6']
        # --- Call bladeDerivedParameters for "manual" calculation
        p = bladeDerivedParameters(p, inertiaAtBladeRoot=inertiaAtBladeRoot)
        # --- Compare both "manual" and GMBeam approach
        np.testing.assert_almost_equal(MM[0,0,]           , p['BldMass'])
        np.testing.assert_almost_equal(MM[3:6,3:6]       ,  p['J'])
        np.testing.assert_almost_equal(MM[3:6,0:3]       ,  skew(p['mdCM']))
        np.testing.assert_almost_equal(np.diag(MM[6:,6:]),  np.diag(p['Me']))
        np.testing.assert_almost_equal(MM[0:3,6:],          p['Ct'].T)
        np.testing.assert_almost_equal(MM[3:6,6:],          p['Cr'].T)
        np.testing.assert_almost_equal(IT['Oe6_M1'], p['Oe_M1'])
        np.testing.assert_almost_equal(IT['mdCM_M1'], p['mdCM_M1'])
# 


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
        np.testing.assert_almost_equal(np.diag(p['MTFA'])/1e6,np.array([4.010392552595473e+05, 2.755960009968935e+07])/1e6-p['TwrTpMass']/1e6)
        np.testing.assert_almost_equal(np.diag(p['MTSS'])/1e6,np.array([4.008128638209492e+05, 3.558329856444363e+07])/1e6-p['TwrTpMass']/1e6)
        np.testing.assert_almost_equal(       (p['KTFA'])/1e6, np.array([[1.910902930730452e+06   ,2.908366961340314e+06],[2.908366961340315e+06   ,1.024457413318486e+10]])/1e6)
        np.testing.assert_almost_equal(       (p['KTSS'])/1e6, np.array([[1.838181636362126e+06  , 5.692201925965283e+06],[  5.692201925965278e+06 ,  1.319779713916086e+10]])/1e6)
        np.testing.assert_almost_equal(np.diag(p['KTFAGrav']),[ 1.066694734234754e+03 , 1.169976720444022e+06])
        np.testing.assert_almost_equal(np.diag(p['KTSSGravTT']), [1.491565303636175e-02 , 1.838391204221678e+01])

        np.testing.assert_almost_equal(p['AxRedTFA'][0,0,-3:],[1.158577163756421e-02,1.332102574930896e-02,1.417295499349627e-02])

        # Frequencies
        np.testing.assert_almost_equal(p['FreqTSS'], np.array([[0.931115438458394 , 0.346445103037537], [3.080181973131150 , 3.074014122773748]]))
        np.testing.assert_almost_equal(p['CTFA'], np.array([[6.420569660348690e+03  , 2.997894941172623e+03],[9.772015298550068e+03,1.055993187124818e+07]]))

        # --------------------------------------------------------------------------------}
        # --- Mass Matrix using GM or not
        # --------------------------------------------------------------------------------{
        from welib.yams.flexibility import GMBeam, GKBeam
        # --- Call bladeDerivedParameters for "manual" calculation
        p = towerDerivedParameters(p)

        # --- Calling GM/GK Beam with OpenFAST method
        s_G0 = np.zeros((3, len(p['s_span'])))
        s_G0[2,:] = p['s_span']
        MM, IT = GMBeam(s_G0, p['s_span'], p['m_full'], p['U'], rot_terms=True, method='OpenFAST', main_axis='z', M1=True) 
        Gr, Ge, Oe, Oe6 = IT['Gr'], IT['Ge'], IT['Oe'], IT['Oe6']

        KK = GKBeam(p['s_span'], p['EI'], p['K'], bOrth=False, method='OpenFAST')
        Mtop = p['TwrTpMass']
        KKg_SW = GKBeamStiffnening(p['s_span'], p['V'], Gravity, p['m_full'], Mtop=0, Omega=0, bSelfWeight=True,  bMtop=False, bRot=False, main_axis='z')
        KKg_TM = GKBeamStiffnening(p['s_span'], p['V'], Gravity, p['m_full'], Mtop=1, Omega=0, bSelfWeight=False, bMtop=True, bRot=False, main_axis='z')
        #KKg_Om = GKBeamStiffnening(p['s_span'], p['V'], Gravity, p['m_full'], Mtop=0, Omega=1, bSelfWeight=False, bMtop=False, bRot=True, main_axis='z')

        # --- 
#         print('KKg_SW\n',KKg_SW[6:,6:])
#         print('KKg_SW\n',p['Kg_SW'])
#         print('KKg_TM\n',KKg_TM[6:,6:])
#         print('KKg_TM\n',p['Kg_TM'])
#         print('Ke\n',KK[6:,6:])
#         print('Ke\n',p['Ke'])

#         print('MM',MM[0,0])
#         print('Ms',p['TwrMass'])
#         print('J\n',p['J'])
#         print('J\n',MM[3:6,3:6])
#         print('mdCM_GM\n',MM[3:6,0:3])
#         print('mdCM_OF\n',p['mdCM'])
#         print('me_GM\n',MM[6:,6:])
#         print('me_OF\n',p['Me'])
#         print('Ct_GM\n',MM[0:3,6:])
#         print('Ct_OF\n',p['Ct'].T)
#         print('Cr_GM\n',MM[3:6,6:])
#         print('Cr_OF\n',p['Cr'].T)

        # --- Compare both "manual" and GMBeam approach
        np.testing.assert_almost_equal(MM[0,0,]          , p['TwrMass'])
        np.testing.assert_almost_equal(MM[3:6,3:6]/1e6   , p['J']/1e6)
        np.testing.assert_almost_equal(MM[3:6,0:3]       , skew(p['mdCM']))
        np.testing.assert_almost_equal(MM[6:,6:]         , p['Me'])
        np.testing.assert_almost_equal(MM[0:3,6:],         p['Ct'].T)
        np.testing.assert_almost_equal(MM[3:6,6:],         p['Cr'].T)

        np.testing.assert_almost_equal(KK[6:,6:]/1e6,         p['Ke']/1e6)





if __name__ == '__main__':
    unittest.main()
