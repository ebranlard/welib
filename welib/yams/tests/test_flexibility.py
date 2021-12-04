import unittest
import numpy as np
import os
from welib.yams.flexibility import *
import welib.weio as weio

MyDir=os.path.dirname(__file__)

# --------------------------------------------------------------------------------}
# --- TESTS
# --------------------------------------------------------------------------------{
class Test(unittest.TestCase):
    def test_GMGKBeam(self):
        #  
        #Given mass, stiffness distribution along the span of a beam and a set of shape functions
        #compute the Generalized mass and stiffness matrix for the flexible beam.

        #In that test: 
        #- as shape functions the "theoretical" shape function for a clamped-free beam
        # (see welib.beams.theory)
        #- constant properties along the beam
        #- Beam along x axis
        try:
            import welib.beams.theory as bt
        except:
            print('[FAIL] Loading beam theory')
            pass

        np.set_printoptions(linewidth=500)
            
        # --- Reference data
        MM_ref=np.array([[ 30000.,      0.,      0.00000,      0.,         0.00000,         0.,      0.00000,      0.00000,      0.00000],
                         [     0.,  30000.,      0.00000,      0.,         0.00000,    900000.,      0.00000,      0.00000,      0.00000],
                         [     0.,      0.,  30000.00000,      0.,   -900000.00000,         0.,  11748.96793,  -6494.82063,   3839.68233],
                         [     0.,      0.,      0.00000,6000000.,         0.00000,         0.,      0.00000,      0.00000,      0.00000],
                         [     0.,      0.,-900000.00000,      0.,  36000000.00000,         0.,-512010.35981,  81016.00516, -30396.91796],
                         [     0., 900000.,      0.00000,      0.,         0.00000,  36000000.,      0.00000,      0.00000,      0.00000],
                         [     0.,      0.,  11748.96793,      0.,   -512010.35981,         0.,   7508.18374,     18.30346,     27.42335],
                         [     0.,      0.,  -6494.82063,      0.,     81016.00516,         0.,     18.30346,   7528.42330,     37.54289],
                         [     0.,      0.,   3839.68233,      0.,    -30396.91796,         0.,     27.42335,     37.54289,   7546.66429]])
        
        KKg_ref=np.array([[ 286478.07306 , -4376.65199 , 18360.80780],[-4376.65199,  11281454.27909 ,  -157525.64695],[18360.80780,-157525.64695  ,88662737.01300]])

        MM2_ref=np.array([[30000.00000,     0.,      0.00000,       0.00000,    11730.33344,       0.00000,   -196.26573,  -52.46587,   134.55304],
                          [    0.00000, 30000.,      0.00000,  -90000.00000,        0.00000,  900000.00000,      0.00000,    0.00000,     0.00000],
                          [    0.00000,     0.,  30000.00000,   45000.00000,  -900000.00000,       0.00000,  11748.96793,-6494.82063,  3839.68233],
                          [    0.00000,-90000.,  45000.00000, 6450267.53864, -1800000.00000,-3600000.00000,  25618.35390,-4032.96435,  1537.68181],
                          [ 11730.33344,     0.,-900000.00000,-1800000.00000, 36360214.03092, -180107.01546,-512010.35981,81016.00516,-30396.91796],
                          [    0.00000,900000.,      0.00000,-3600000.00000,  -180107.01546,36090053.50773,      0.00000,    0.00000,     0.00000],
                          [ -196.26573,     0.,  11748.96793,   25618.35390,  -512010.35981,       0.00000,   7508.18374,   18.30346,    27.42335],
                          [  -52.46587,     0.,  -6494.82063,   -4032.96435,    81016.00516,       0.00000,     18.30346, 7528.42330,    37.54289],
                          [  134.55304,     0.,   3839.68233,    1537.68181,   -30396.91796,       0.00000,     27.42335,   37.54289,  7546.66429]])


        # --- Setting up beam properties and shape functions
        nShapes = 3
        nSpan   = 30
        L       = 60
        EI0     = 2E+10 # [Nm^2]
        m       = 5E+2  # [kg/m]
        GKt     = 7e11  # [Nm2]
        jxx     = 1e5   # [kg.m]
        A=1; rho=A*m;

        x=np.linspace(0,L,nSpan);
        # Shape functions (taken as analytical mode shapes)
        freq,s_span,U,V,K = bt.UniformBeamBendingModes('unloaded-clamped-free',EI0,rho,A,L,x=x)
        PhiU = np.zeros((nShapes,3,nSpan)) # Shape
        PhiV = np.zeros((nShapes,3,nSpan)) # Slope
        PhiK = np.zeros((nShapes,3,nSpan)) # Curvature
        for j in np.arange(nShapes):  
            PhiU[j][2,:] = U[j,:] # Setting modes along z
            PhiV[j][2,:] = V[j,:]
            PhiK[j][2,:] = K[j,:]
        m    = m*np.ones(nSpan)
        jxxG = jxx*np.ones(nSpan)
        EI= np.zeros((3,nSpan))
        EI[1,:] = EI0
        EI[2,:] = EI0
        # if ~isempty(p.GKt)
        #       B.GKt=p.GKt*ones(size(B.s_span));
        
        # --- Testing for straight COG
        s_G      = np.zeros((3,nSpan))
        s_G[0,:] = x
        MM = GMBeam(s_G, s_span, m, PhiU, jxxG=jxxG, bUseIW=True, main_axis='x') # Ref uses IW_xm
        KK = GKBeam(s_span, EI, PhiK)

        #np.testing.assert_equal(np.all(MDiff<1e-3),True)
        np.testing.assert_allclose(MM,MM_ref,rtol=1e-5)
        np.testing.assert_allclose(KK[6:,6:],KKg_ref,rtol=1e-5)

        # --- Testing for curved COG
        s_G      = np.zeros((3,nSpan))
        s_G[0,:] = x
        s_G[1,:] = x/20
        s_G[2,:] = x/10
        V_tot=PhiV[0]
        MM = GMBeam(s_G, s_span, m, PhiU, jxxG=jxxG, bUseIW=True, main_axis='x', V=PhiV, bAxialCorr=True,V_tot=V_tot) # Ref uses IW_xm
        ##np.testing.assert_equal(np.all(MDiff<1e-3),True)
        np.testing.assert_allclose(MM,MM2_ref,rtol=1e-5)
        #np.testing.assert_allclose(KK[6:,6:],KKg_ref,rtol=1e-5)

    # --------------------------------------------------------------------------------}
    # --- NREL5MW Tower 
    # --------------------------------------------------------------------------------{
    def test_TowerBeam_SID(self):
        # 
        #In this test we use:
        #- Beam properties from the Tower of the NREL5MW
        #- Shape functions determined using a FEM beam representation (see welib/FEM/tests/test_beam_linear_element.py)
        #- We test for the gyroscopic and centrifugal matrix Gr, Ge, Oe
        #- Beam along z axis
        np.set_printoptions(linewidth=300, precision=9)
        # --- Read data from NREL5MW tower
        TowerHt=87.6;
        TowerBs=0;
        TwrFile=os.path.join(MyDir,'./../../../data/NREL5MW/5MW_Baseline/NRELOffshrBsline5MW_Onshore_ElastoDyn_Tower.dat')
        twr = weio.FASTInputFile(TwrFile).toDataFrame()
        z   = twr['HtFract_[-]']*(TowerHt-TowerBs)
        m   = twr['TMassDen_[kg/m]']  
        nSpan = len(z)

        # --- Shape function taken from FEM
        shapeFile = os.path.join(MyDir,'../../../data/NREL5MW/NREL5MW_Tower_Onshore_FEM_Modes.csv')
        shapes = weio.read(shapeFile).toDataFrame()
        nShapes=2
        PhiU = np.zeros((nShapes,3,nSpan)) # Shape
        PhiV = np.zeros((nShapes,3,nSpan)) # Slope
        s_G      = np.zeros((3,nSpan))
        main_axis='z'
        if main_axis=='x':
            PhiU[0,1,:]=shapes['U1'] # along y
            PhiU[1,2,:]=shapes['U2'] # along z
            PhiV[0,1,:]=shapes['V1'] # along y 
            PhiV[1,2,:]=shapes['V2'] # along z 
            s_G[0,:] = z
        elif main_axis=='z':
            PhiU[0,0,:]=shapes['U1'] # along x
            PhiU[1,1,:]=shapes['U2'] # along y
            PhiV[0,0,:]=shapes['V1'] # along x (around theta y) # TODO double check convention
            PhiV[1,1,:]=shapes['V2'] # along y (around theta x  # TODO double check convention
            s_G[2,:] = z

        # --- Testing for straight COG
        s_span = z
        jxxG= z*0 + m # NOTE: unknown
        #MM = GMBeam(s_G, s_span, m, PhiU, jxxG=jxxG, bUseIW=True, main_axis='x') # Ref uses IW_xm
        Mxx, Mtt, Mxt, Mtg, Mxg, Mgg, Gr, Ge, Oe, Oe6 = GMBeam(s_G, s_span, m, PhiU, jxxG=jxxG, bUseIW=True, main_axis=main_axis, split_outputs=True, rot_terms=True)
        MM, Gr, Ge, Oe, Oe6 = GMBeam(s_G, s_span, m, PhiU, jxxG=jxxG, bUseIW=True, main_axis=main_axis, split_outputs=False, rot_terms=True)
        MM_ref=np.array( # NOTE: this is onshore tower
        [[ 3.474602316e+05,  0.000000000e+00,  0.000000000e+00,  0.000000000e+00,  1.322633773e+07, -0.000000000e+00,  1.046938838e+05,  0.000000000e+00],
         [ 0.000000000e+00,  3.474602316e+05,  0.000000000e+00, -1.322633773e+07,  0.000000000e+00,  0.000000000e+00,  0.000000000e+00,  1.046938838e+05],
         [ 0.000000000e+00,  0.000000000e+00,  3.474602316e+05,  0.000000000e+00, -0.000000000e+00,  0.000000000e+00,  0.000000000e+00,  0.000000000e+00],
         [ 0.000000000e+00, -1.322633773e+07,  0.000000000e+00,  7.182575651e+08, -0.000000000e+00, -0.000000000e+00,  0.000000000e+00, -6.440479129e+06],
         [ 1.322633773e+07,  0.000000000e+00, -0.000000000e+00, -0.000000000e+00,  7.182575651e+08, -0.000000000e+00,  6.440479129e+06,  0.000000000e+00],
         [-0.000000000e+00,  0.000000000e+00,  0.000000000e+00, -0.000000000e+00, -0.000000000e+00,  3.474602316e+05,  0.000000000e+00,  0.000000000e+00],
         [ 1.046938838e+05,  0.000000000e+00,  0.000000000e+00,  0.000000000e+00,  6.440479129e+06,  0.000000000e+00,  6.145588498e+04,  0.000000000e+00],
         [ 0.000000000e+00,  1.046938838e+05,  0.000000000e+00, -6.440479129e+06,  0.000000000e+00,  0.000000000e+00,  0.000000000e+00,  6.145588498e+04]])

        np.testing.assert_allclose(Mxx, MM_ref[0:3,0:3] ,rtol=1e-4)
        np.testing.assert_allclose(Mxt, MM_ref[0:3,3:6] ,rtol=1e-4)
        np.testing.assert_allclose(Mxg, MM_ref[0:3,6::] ,rtol=1e-4)
        np.testing.assert_allclose(Mtt, MM_ref[3:6,3:6] ,rtol=1e-4)
        np.testing.assert_allclose(Mtg, MM_ref[3:6,6: ] ,rtol=1e-4)
        np.testing.assert_allclose(Mgg, MM_ref[6: ,6: ] ,rtol=1e-4)

        np.testing.assert_allclose(Gr[0][0,:],[0,0,-12945834],rtol=1e-5)
        np.testing.assert_allclose(Gr[1][1,:],[0,0,-12945834],rtol=1e-5)
        np.testing.assert_allclose(Ge[0][1,:],[0,0,122911],rtol=1e-5)
        np.testing.assert_allclose(Ge[1][0,:],[0,0,-122911],rtol=1e-5)
        np.testing.assert_allclose(Oe6[0][:],[0,0,0,0,0,6472917],rtol=1e-5)
        np.testing.assert_allclose(Oe6[1][:],[0,0,0,0,6472917,0],rtol=1e-5)

    # --------------------------------------------------------------------------------}
    # --- NREL5MW Blade 
    # --------------------------------------------------------------------------------{
    def test_BladeBeam_SID(self):
        # 
        #In this test we use:
        #- Beam properties from the Blade of the NREL5MW
        #- Shape functions determined using a FEM beam representation (see welib/yams/tests/test_sid.py)
        #- We test for the gyroscopic and centrifugal matrix Gr, Ge, Oe
        #- Beam along z axis
        #
        np.set_printoptions(linewidth=300, precision=9)

        # --- Read data from NREL5MW Blade
        edFile=os.path.join(MyDir,'./../../../data/NREL5MW/offshore/NREL5MW_ED_Offshore.dat')
        parentDir=os.path.dirname(edFile)
        ed = weio.FASTInputFile(edFile)
        TipRad = ed['TipRad']
        HubRad = ed['HubRad']
        BldLen= TipRad-HubRad;
        BldFile = ed['BldFile(1)'].replace('"','')
        BldFile=os.path.join(parentDir,BldFile)
        bld = weio.FASTInputFile(BldFile).toDataFrame()
        z   = bld['BlFract_[-]']*BldLen + HubRad
        m   = bld['BMassDen_[kg/m]']
        nSpan = len(z)

        # --- Shape function taken from FEM
        shapeFile = os.path.join(MyDir,'../../../data/NREL5MW/NREL5MW_Blade_FEM_Modes.csv')
        shapes = weio.read(shapeFile).toDataFrame()
        nShapes=2
        PhiU = np.zeros((nShapes,3,nSpan)) # Shape
        s_G  = np.zeros((3,nSpan))
        main_axis='z'
        if main_axis=='z':
            # Mode 1
            PhiU[0,0,:]=shapes['U1x']
            PhiU[0,1,:]=shapes['U1y']
            # Mode 2
            PhiU[1,0,:]=shapes['U2x']
            PhiU[1,1,:]=shapes['U2y']
            s_G[2,:] = z

        # --- Testing for straight COG
        s_span = z
        jxxG= z*0 + m # NOTE: unknown
        MM, Gr, Ge, Oe, Oe6 = GMBeam(s_G, s_span, m, PhiU, jxxG=jxxG, bUseIW=True, main_axis=main_axis, split_outputs=False, rot_terms=True)

        MM_ref=np.array(
        [[ 1.684475202e+04,  0.000000000e+00,  0.000000000e+00,  0.000000000e+00,  3.707069074e+05, -0.000000000e+00,  1.976263092e+03, -6.366476591e+00],
         [ 0.000000000e+00,  1.684475202e+04,  0.000000000e+00, -3.707069074e+05,  0.000000000e+00,  0.000000000e+00,  4.082717323e+00,  2.915664880e+03],
         [ 0.000000000e+00,  0.000000000e+00,  1.684475202e+04,  0.000000000e+00, -0.000000000e+00,  0.000000000e+00,  0.000000000e+00,  0.000000000e+00],
         [ 0.000000000e+00, -3.707069074e+05,  0.000000000e+00,  1.225603507e+07, -0.000000000e+00, -0.000000000e+00, -1.817970390e+02, -1.200295386e+05],
         [ 3.707069074e+05,  0.000000000e+00, -0.000000000e+00, -0.000000000e+00,  1.225603507e+07, -0.000000000e+00,  8.737268222e+04, -2.574073558e+02],
         [-0.000000000e+00,  0.000000000e+00,  0.000000000e+00, -0.000000000e+00, -0.000000000e+00,  1.684475202e+04,  0.000000000e+00,  0.000000000e+00],
         [ 1.976263092e+03,  4.082717323e+00,  0.000000000e+00, -1.817970390e+02,  8.737268222e+04,  0.000000000e+00,  8.364646779e+02, -9.586291225e-04],
         [-6.366476591e+00,  2.915664880e+03,  0.000000000e+00, -1.200295386e+05, -2.574073558e+02,  0.000000000e+00, -9.586291225e-04,  1.351321852e+03]])

        np.testing.assert_allclose(MM[0:3,0:3], MM_ref[0:3,0:3] ,rtol=1e-4)
        np.testing.assert_allclose(MM[0:3,3:6], MM_ref[0:3,3:6] ,rtol=1e-4)
        np.testing.assert_allclose(MM[0:3,6: ], MM_ref[0:3,6::] ,rtol=1e-4)
        np.testing.assert_allclose(MM[3:6,3:6], MM_ref[3:6,3:6] ,rtol=1e-4)
        np.testing.assert_allclose(MM[3:6,6: ], MM_ref[3:6,6: ] ,rtol=1e-4)
        np.testing.assert_allclose(MM[6::,6::], MM_ref[6: ,6: ] ,rtol=1e-4)
        np.testing.assert_allclose(MM,         MM.T ,rtol=1e-4)

        np.testing.assert_allclose(Gr[0][0,:],[0,0,-174817]  ,rtol   = 1e-5)
        np.testing.assert_allclose(Gr[0][1,:],[0,0,-363.7477],rtol = 1e-5  )
        np.testing.assert_allclose(Gr[1][0,:],[0,0,514.944]  ,rtol  = 1e-5 )
        np.testing.assert_allclose(Gr[1][1,:],[0,0,-240127]  ,rtol  = 1e-5 )

        np.testing.assert_allclose(Ge[0][1,:],[0,0,2087.767],rtol=1e-5)
        np.testing.assert_allclose(Ge[1][0,:],[0,0,-2087.767],rtol=1e-5)
        np.testing.assert_allclose(Oe6[0][:],[0,0,0,0,181.8739,87408.6],rtol=1e-5)
        np.testing.assert_allclose(Oe6[1][:],[0,0,0,0,120063,-257.472],rtol=1e-5)


if __name__=='__main__':
    #Test().test_TowerBeam_SID()
    #Test().test_BladeBeam_SID()
    unittest.main()
