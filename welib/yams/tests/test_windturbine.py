import unittest
import os
import numpy as np
from welib.yams.windturbine import *
from welib.yams.utils import *
import welib.weio as weio

MyDir=os.path.dirname(__file__)


# --------------------------------------------------------------------------------}
# --- TESTS Wind Turbine, With OpenFAST Algo. Should give perfect match, but some TODOs
# --------------------------------------------------------------------------------{
class TestWindTurbSparOpenFAST(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # Read FAST structural model
        cls.WT = FASTWindTurbine(os.path.join(MyDir,'./../../../data/Spar/Main_Spar_ED.fst'), algo='OpenFAST') #, bldStartAtRotorCenter=False )

    def test_bld(self):
        # --- Blade
        # NOTE: bldes have "R" as origin
        bld0=self.WT.bld[0]
        #print(bld0)
        #print('start_pos',bld0.start_pos)
        #print('s_span',bld0.s_span)
        #print('MM\n',bld0.MM[:6,:6])
        #print('firsdtmom',bld0.first_moment_inertia)
	# Comparison with ElastoDyn summary file
        #   Mass                  (kg)        17600.781    17600.781    17600.781
        #   Second Mass Moment    (kg-m^2) 11688798.933 11688798.933 11688798.933
        #   First Mass Moment     (kg-m)     361234.779   361234.779   361234.779
        #   Center of Mass        (m)            20.524       20.524       20.524
        np.testing.assert_almost_equal(bld0.mass,                              17600.781, 3)
        np.testing.assert_almost_equal(bld0.inertia_at(bld0.start_pos)[0,0],11688798.933, 3)
        np.testing.assert_almost_equal((bld0.first_moment_inertia_from_start),(0,0,361234.779), 3) 
        np.testing.assert_almost_equal((bld0.masscenter-bld0.start_pos)      ,[0,0,20.524], 3)

        # --- Contribution of blade to rotor inertia along x
        HubRad    = 1.5
        PreCone   = -2.5*np.pi/180
        #zCG       = 20.52379
        #BldMass   = 17600.7808000
        #SecondMom = 11688798.933
        SecondMom = bld0.inertia_at(bld0.start_pos)[0,0]
        BldMass   = bld0.mass
        zCG       = bld0.masscenter[2] -bld0.start_pos[2]
        bldIner = (SecondMom + BldMass*HubRad*(2.0*zCG +HubRad))*(np.cos(PreCone)**2)
        np.testing.assert_almost_equal(bldIner                             , 12787728.064229555 , 3)
        np.testing.assert_almost_equal(bld0.inertia_global_at((0,0,0))[0,0], 12787728.064229555 , 3)

        # --- Coordinates of Blade CG in "Hub" coordinates
        # Case 1 (HubRad straight)
        zCG_R = zCG*np.cos(PreCone) + HubRad
        xCG_R = zCG*np.sin(PreCone)
        # Case 2 (HubRad along blade axis) - (OpenFAST)
        zCG_R = (zCG + HubRad)*np.cos(PreCone)
        xCG_R = (zCG + HubRad)*np.sin(PreCone)
        np.testing.assert_almost_equal(bld0.masscenter_pos_global, (xCG_R, 0, zCG_R),4)

        #from welib.fast.elastodyn import bladeParameters, rotorParameters
        #prot, pbld, phub = rotorParameters(self.WT.ED.filename, identicalBlades=True, pbld1=None)
        #print('bldIner     '  ,bldIner  )
        #print('3 bldIner   ',3*bldIner)
        #print('3 bldIner+H ',3*bldIner+phub['HubIner'])
        #print('RotIner     ',prot['RotIner'])
        #print('bldIner     '  ,bldIner  )
        #Jb_O_b = bld0.inertia_at(bld0.start_pos)
        #C=np.cos(PreCone)
        #S=np.sin(PreCone)
        ## --- Inertia of Blade at Blade root, in Global
        #Jx = Jb_O_b[0,0]
        #Jy = Jb_O_b[1,1] # = Jx
        #Jz = Jb_O_b[2,2] # 0 for OpenFAST
        #Jb_O_gl = np.array([
        #        [C**2 *Jx + S**2 * Jz , 0  , -C*S*Jx + C*S*Jz     ], 
        #        [ 0                   , Jy , 0                    ], 
        #        [ -C*S*Jx + C*S*Jz    , 0  , C**2 * Jz + S**2 * Jx]] )
        ##print('Jb_O_gl\n',Jb_O_gl)
        ##J = bld0.inertia_global_at((0,0,1.5))
        ##print('J\n',J)
        ##J = bld0.inertia_at((0,0,0), np.eye(3))
        ##print('J\n',J)

        ## --- Inertia of Blade at Blade COG in body
        #Jb_G_b = np.array([
        #        [Jx - BldMass*zCG**2, 0                , 0 ],
        #        [ 0                 , Jy-BldMass*zCG**2, 0 ],
        #        [ 0                 , 0                , Jz]])
        ##print('Jb_G_b\n',Jb_G_b)
        ##print('Jb_G_b\n',bld0.masscenter_inertia)
        ## --- Inertia of Blade at Blade COG in global
        #Jx = Jb_G_b[0,0]
        #Jy = Jb_G_b[1,1] # = Jx
        #Jz = Jb_G_b[2,2] # 0 for OpenFAST
        #Jb_G_gl = np.array([
        #        [C**2 *Jx + S**2 * Jz , 0  , -C*S*Jx + C*S*Jz     ], 
        #        [ 0                   , Jy , 0                    ], 
        #        [ -C*S*Jx + C*S*Jz    , 0  , C**2 * Jz + S**2 * Jx]] )
        #print('Jb_G_gl\n',Jb_G_gl)
        #print(bld0.inertia_global_at(bld0.masscenter_pos_global))
        #Jb_R_gl = bld0.inertia_global_at((0,0,0))
        #print('Jb_R_gl\n',Jb_R_gl)
        #from welib.fast.elastodyn import bladeParameters
        #pED = bladeParameters(self.WT.ED.filename)
        #print('BldMass   {:15.3f}{:15.3f}'.format(pED['BldMass']  ,    17600.781))
        #print('SecondMom {:15.3f}{:15.3f}'.format(pED['SecondMom'], 11688798.933))
        #print('FirstMom  {:15.3f}{:15.3f}'.format(pED['FirstMom'] ,   361234.779))
        #print('CG        {:15.3f}{:15.3f}'.format(pED['BldCG']    ,       20.524))

    def test_twr(self):
        # --- Tower
        twr=self.WT.twr
        #print(twr)
	# Comparison with ElastoDyn summary file
        # Tower Mass            (kg)       217511.246
        np.testing.assert_almost_equal(twr.mass,        217511.246, 3)
        np.testing.assert_almost_equal(twr.length,      67.6)
        # np.testing.assert_allclose(np.diag(twr.inertia), [2.6183e8,2.6183e8, 0.0], 1e-3)

    def test_rotor(self):
        # NOTE: rotor is hub + rigid blades, with "R" as origin, but define with "N" as global
        rot = self.WT.rot
        # Rotor Mass            (kg)       109582.342
        # Rotor Inertia         (kg-m^2) 38479110.193
        np.testing.assert_almost_equal(rot.mass,           109582.342, 3)
        np.testing.assert_almost_equal(rot.inertia[0,0], 38479110.193, 3) # TODO TODO TODO WRONG BECAUSE OF BLADE DEFINITION


# --------------------------------------------------------------------------------}
# --- TESTS Wind Turbine, Without OpenFAST Algo
# --------------------------------------------------------------------------------{
class TestWindTurbSpar(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # Read FAST structural model
        cls.WT = FASTWindTurbine(os.path.join(MyDir,'./../../../data/Spar/Main_Spar_ED.fst'), algo='')

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
        # NOTE: bldes have "R" as origin
        bld0=self.WT.bld[0]
        # print(bld0)
        #print(bld0.start_pos)
        #print(bld0.end_pos)
        #print(bld0.inertia_at(bld0.start_pos))

	# Comparison with ElastoDyn summary file
        #Mass                  (kg)        17600.773
        #Second Mass Moment    (kg-m^2) 11688785.000
        #First Mass Moment     (kg-m)     361234.438
        #Center of Mass        (m)            20.524
        np.testing.assert_allclose(bld0.mass,                              17600.773, 1e-3)
        np.testing.assert_allclose(bld0.inertia_at(bld0.start_pos)[0,0],11688785.000, 1e-3)
        np.testing.assert_allclose((bld0.masscenter-bld0.start_pos)    ,[0,0,20.524], 1e-3)
	

    def test_rotor(self):
        # NOTE: rotor is hub + rigid blades, with "R" as origin, but define with "N" as global
        rot = self.WT.rot
	# Rotor Mass            (kg)       109582.320
	# Rotor Inertia         (kg-m^2) 38479064.000
        np.testing.assert_allclose(rot.mass,           109582.320, 1e-3)
        np.testing.assert_allclose(rot.inertia[0,0], 38479064.000, 1e-4)

    def test_RNA(self):
        RNA = self.WT.RNA
	#Tower-top Mass        (kg)       349582.312
        np.testing.assert_allclose(RNA.mass,  349582.312, 1e-3)
        np.testing.assert_allclose(RNA.inertia[0,0],  43972846, 1e-3)
        np.testing.assert_allclose(np.around(RNA.masscenter_pos_global,5), [-0.40774, 0, 1.96643], 1e-2)

if __name__=='__main__':
    np.set_printoptions(linewidth=300, precision=5)
    unittest.main()
