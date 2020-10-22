
import unittest
import numpy as np
from welib.yams.flexibility import *


# --------------------------------------------------------------------------------}
# --- TESTS
# --------------------------------------------------------------------------------{
class Test(unittest.TestCase):
    def test_rot(self):
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


        # --- Setting up mode shapes
        nShapes=3;
        nSpan=30;
        L   = 60  ; EI0 = 2E+10; m = 5E+2;
        GKt = 7e11# [Nm2]
        jxx = 1e5 # [kg.m]
        A=1; rho=A*m;

        x=np.linspace(0,L,nSpan);
        # Mode shapes
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



if __name__=='__main__':
    unittest.main()
