import unittest

import numpy as np
from welib.yams.utils import *
from welib.yams.flexibility import GMBeam
from welib.essentials import *


# --------------------------------------------------------------------------------}
# --- Helper functions 
# --------------------------------------------------------------------------------{
def dimless(MM, M, L):
    MM /= M
    MM[3:,3:] /=L**2
    return MM

def getUBeam(jxx=1, axis='x'):
    nSpan = 200
    L     = 60
    m     = 5E+2                    # [kg/m]
    jxx   = 1e5*jxx                 # [kg.m]
    M     = m*L
    x     = np.linspace(0,L,nSpan);
    m     = m*np.ones(nSpan)
    jxxG  = jxx*np.ones(nSpan)
    s_G      = np.zeros((3,len(x)))
    s_span   = None
    if axis=='z':
        s_G[2,:] = x
    elif axis=='x':
        s_G[0,:] = x
    else:
        pass
    return x, m, jxxG, L, M, s_G, s_span

# --------------------------------------------------------------------------------}
# --- TESTS
# --------------------------------------------------------------------------------{
class TestUtils(unittest.TestCase):
    def test_rot(self):
        # --- Identity matrix for 0 rotation
        np.testing.assert_almost_equal(R_x(0),np.eye(3))

    def test_skew(self):
        # --- Testing  \tilde{u} . v  ==  u x v
        u=np.array([1,2,3])
        v=np.array([-1,0,-3])
        np.testing.assert_equal(np.cross(u,v), np.dot(skew(u),v))

        # --- Testing  -[u~][u~] = [u~]^T [u~] = |u|^2 I - u u^T
        x=u
        S2_1 = (np.dot(x, x) * np.eye(3) - np.outer(x, x))
        S2_2 =   skew(x).T @ skew(x)
        S2_3 = - skew2(x)
        np.testing.assert_equal(S2_1, S2_2)
        np.testing.assert_equal(S2_1, S2_3)




    def test_inertia(self):
        # --- Transferring inertia at G to point A and B and then to each other
        I_G=np.diag([1,2,3]); 
        M=2;
        r_OG = np.array([ 1, 2, 10 ])
        r_OA = r_OG + np.array([5, 8 , 2] )
        r_OB = r_OG + np.array([4, -6, -3])
        r_AG = r_OG-r_OA
        r_BG = r_OG-r_OB
        I_A  = translateInertiaMatrix(I_G,M,r_AG)        # I@ A
        I_B  = translateInertiaMatrix(I_G,M,r_BG)        # I@ B
        I_B2 = translateInertiaMatrix(I_A,M,r_BG,r_AG   )# I@B from A
        I_A2 = translateInertiaMatrix(I_B,M,r_AG,r_BG   )# I@A from B
        np.testing.assert_equal(I_A,I_A2)
        np.testing.assert_equal(I_B,I_B2)

        # --- Transfer of inertia at A to COG then back at A
        I_A = np.eye(3)
        M = 12
        r_GA = np.array([3,4,5])
        I_G  = translateInertiaMatrixToCOG  (I_A,M,-r_GA)  # I@G from A
        I_A2 = translateInertiaMatrixFromCOG(I_G,M, r_GA) # I@A from G
        np.testing.assert_equal(I_A,I_A2)

    def testBeamMM_UBeamz(self):
        # --- Testing for straight uniform beam along z 
        x, m, jxxG, L, M, s_G, s_span = getUBeam(jxx=0, axis='z')
        MM_refz   = rigidBodyMassMatrixAtP(m=M, J_G=(M*L**2/12, M*L**2/12, 0), Ref2COG=(0,0,L/2))
        MM, _     = GMBeam(s_G, s_span, m, jxxG = jxxG, method = 'trapz', main_axis = 'z')
        MM3, _    = GMBeam(s_G, s_span, m, jxxG = jxxG, method = 'Flex', main_axis = 'z')
        MM2, info = rigidBeamMassMatrix(s_G, m, jxxG = jxxG)

        MM      = dimless(MM     ,M,L)
        MM2     = dimless(MM2    ,M,L)
        MM3     = dimless(MM3    ,M,L)
        MM_refz = dimless(MM_refz,M,L)

        #print   (info     )
        #printMat(MM_refz,'MM_refz', digits=4)
        #printMat(MM     ,'MM'   , digits = 4)
        #printMat(MM2    ,'MM2'  , digits = 4 )
        #printMat(MM2-MM ,'delta', digits = 4)
        #printMat(MM2-MM_refz ,'delta', digits = 4)
        np.testing.assert_array_almost_equal(MM , MM2    , 6)
        np.testing.assert_array_almost_equal(MM , MM3    , 4)
        np.testing.assert_array_almost_equal(MM2, MM_refz, 5)

    def testBeamMM_UBeamx(self):
        # --- Testing for straight uniform beam along x 
        x, m, jxxG, L, M, s_G, s_span = getUBeam(jxx=0, axis='x')
        MM_refx = rigidBodyMassMatrixAtP(m=M, J_G=(0, M*L**2/12, M*L**2/12), Ref2COG=(L/2,0,0))
        MM, _     = GMBeam(s_G, s_span, m, jxxG = jxxG, method = 'trapz', main_axis = 'x')
        MM3, _    = GMBeam(s_G, s_span, m, jxxG = jxxG, method = 'Flex', main_axis = 'x')
        MM2, info = rigidBeamMassMatrix(s_G, m, jxxG = jxxG)

        MM      = dimless(MM     ,M,L)
        MM2     = dimless(MM2    ,M,L)
        MM3     = dimless(MM3    ,M,L)
        MM_refx = dimless(MM_refx,M,L)

        #printMat(MM_refx,'MM_refx', digits=4)
        #printMat(MM     ,'MM'   , digits = 4)
        #printMat(MM2    ,'MM2'  , digits = 4 )
        #printMat(MM3    ,'MM3'  , digits = 4 )
        #printMat(MM2-MM ,'delta', digits = 4)
        #printMat(MM2-MM3,'delta', digits = 4)
        #printMat(MM2-MM_refx ,'delta', digits = 4)
        np.testing.assert_array_almost_equal(MM, MM2    , 6)
        np.testing.assert_array_almost_equal(MM , MM3   , 5)
        np.testing.assert_array_almost_equal(MM, MM_refx, 5)

    def testBeamMM_UBeamStraightRot(self):
        # --- Testing for straight uniform beam rotated
        x, m, jxxG, L, M, s_Gb, s_span = getUBeam(jxx=0, axis='x')
        J_Gb = np.diag(np.asarray((0, M*L**2/12, M*L**2/12)))
        r_Gb = np.array((L/2,0,0))
        #aa = 30*np.pi/180
        for aa in np.linspace(0,np.pi, 5):
            R_ib = R_z(aa)
            J_Gi = R_ib @ J_Gb @ R_ib.T
            r_Gi = R_ib @ r_Gb
            s_G   =R_ib @ s_Gb 
            MMr    = rigidBodyMassMatrixAtP(m=M, J_G=J_Gi, Ref2COG=r_Gi)
            MM2, info = rigidBeamMassMatrix(s_G, m)

            MMr     = dimless(MMr    ,M,L)
            MM2     = dimless(MM2    ,M,L)
            #printMat(MMr    ,'MM_ref',digits = 5, xmin=1e-6)
            #printMat(MM2    ,'MM2'   ,digits = 5, xmin=1e-6)
            #printMat(MM2-MMr,'delta' ,digits = 5, xmin=1e-6)
            np.testing.assert_array_almost_equal(MM2, MMr, 5)


    def testInertiaTranslate(self):
        x, m, jxxG, L, M, s_AG, s_span = getUBeam(jxx=0,axis='z')

        r_AB = np.array([10, 20, 30])
        s_BG = s_AG*0
        for i in range(len(x)):
            s_BG[:,i] = s_AG[:,i] + r_AB

        MMA, infoA = rigidBeamMassMatrix(s_AG, m, point='A')
        MMB, infoB = rigidBeamMassMatrix(s_BG, m, point='B')
        r_AG = infoA['r_AG']
        r_BG = infoB['r_BG']
        J_G = infoA['J_G']
        J_A = MMA[3:,3:]
        J_B = MMB[3:,3:]

#         printDict(infoA)
#         printDict(infoB)
#         printMat(MMA , 'MMA'               , digits=4 , nchar=14 )
#         printMat(MMB , 'MMB'               , digits=4 , nchar=14 )
#         printMat(J_A , 'J_A'               , digits=4 , nchar=14 )
#         printMat(J_B , 'J_B'               , digits=4 , nchar=14 )
        np.testing.assert_almost_equal(infoA['J_G'], infoB['J_G'])


        # --- Translate To COG
        J_G2 = translateInertiaMatrixToCOG(J_A, M, r_AG) 
        J_G3 = translateInertiaMatrixToCOG(J_B, M, r_BG) 
        #printMat(J_G2, 'J_G2'              , digits=4 , nchar=14 )
        #printMat(J_G3, 'J_G3'              , digits=4 , nchar=14 )
        #printMat(J_G , 'J_G'               , digits=4 , nchar=14 )
        np.testing.assert_almost_equal(J_G2, J_G)
        np.testing.assert_almost_equal(J_G3, J_G)

        # --- Translate From COG
        J_A3 =  translateInertiaMatrixFromCOG(J_G, M, r_GP=-r_AG) 
        J_B3 =  translateInertiaMatrixFromCOG(J_G, M, r_GP=-r_BG) 
        np.testing.assert_almost_equal(J_B3, J_B)
        np.testing.assert_almost_equal(J_A3, J_A)

        # --- Translate A to B via COG
        J_A2 = translateInertiaMatrix(J_B, M, r_AG, r_AG = r_BG)
        J_B2 = translateInertiaMatrix(J_A, M, r_BG, r_AG = r_AG)
        np.testing.assert_almost_equal(J_B2, J_B)
        np.testing.assert_almost_equal(J_A2, J_A)




    def testInertiaThreeSym(self):
        # --- Test analytical formula for 3-symmetry: inertia of three identical bodies at 0 120 240
        x, m, jxxG, L, M, s_Gb, s_span = getUBeam(jxx=0,axis='z')
        s_Gb[0,:]=10*(x/L)**2
        s_Gb[1,:]=5*(x/L)**2
        # --- Computed reference inertia for body 1 (could be done in any coordinate system)
        MM0, info0 = rigidBeamMassMatrix(s_Gb, m)
        x_OG = info0['r_OG'][0]
        M = MM0[0,0]
        Jxx=MM0[3,3]
        Jyy=MM0[4,4]
        Jzz=MM0[5,5]
        MMr = np.diag((3*M, 3*M, 3*M, 3*Jxx, 3/2 *(Jyy+Jzz) , 3/2 *(Jyy+Jzz)))
        # IMPORTANT, need to add x COG offset
        MMr[5,1] = 3*M * x_OG
        MMr[4,2] =-3*M * x_OG
        MMr[1,5] = 3*M * x_OG
        MMr[2,4] =-3*M * x_OG

        for psi in [0,30, 60,90,120]:
            psi = 30
            R_ib1 = R_x((psi+  0)*np.pi/180)
            s_OG1   =R_ib1 @ s_Gb 
            MM1, info1 = rigidBeamMassMatrix(s_OG1, m)
            
            R_ib2 = R_x((psi+120)*np.pi/180)
            s_OG2   =R_ib2 @ s_Gb 
            MM2, info2 = rigidBeamMassMatrix(s_OG2, m)
            R_ib3 = R_x((psi+240)*np.pi/180)
            s_OG3   =R_ib3 @ s_Gb 
            MM3, info3 = rigidBeamMassMatrix(s_OG3, m)

            MM = MM1 + MM2 + MM3
            #printMat(MM0, 'MM0', digits=4, xmin=1e-5)
            #printMat(MM1, 'MM1', digits=4, xmin=1e-5)
            #printMat(MM2, 'MM2', digits=4, xmin=1e-5)
            #printMat(MM3, 'MM3', digits=4, xmin=1e-5)
            #printMat(MM,  'MM',  digits=4, xmin=1e-5)
            #printMat(MMr, 'MMr',  digits=4, xmin=1e-5)

            np.testing.assert_almost_equal(MM, MMr, 5)


    def testInertiaNSym(self):
        # --- Test analytical formula for N-symmetry
        x, m, jxxG, L, M, s_Gb, s_span = getUBeam(jxx=0,axis='z')
        s_Gb[0,:]=10*(x/L)**2
        s_Gb[1,:]=5*(x/L)**2
        # --- Computed reference inertia for body 1 (could be done in any coordinate system)
        MM0, info0 = rigidBeamMassMatrix(s_Gb, m)
        x_OG = info0['r_OG'][0]
        M = MM0[0,0]
        Jxx=MM0[3,3]
        Jyy=MM0[4,4]
        Jzz=MM0[5,5]

        for n in [3, 4, 5, 10]:
            MMr = np.diag((n*M, n*M, n*M, n*Jxx, n/2 *(Jyy+Jzz) , n/2 *(Jyy+Jzz)))
            # IMPORTANT, need to add x COG offset
            MMr[5,1] = n*M * x_OG
            MMr[4,2] =-n*M * x_OG
            MMr[1,5] = n*M * x_OG
            MMr[2,4] =-n*M * x_OG

            for psiOff in [0, 35]:

                MM =np.zeros((6,6))
                vpsi = np.arange(0, 360, 360/n) # Blade azimuth
                for psii in vpsi: 
                    R_ibi = R_x((psiOff+psii)*np.pi/180)
                    s_OGi   =R_ibi @ s_Gb 
                    MMi, info1 = rigidBeamMassMatrix(s_OGi, m)
                    MM += MMi
                np.testing.assert_almost_equal(MM, MMr, 5)











if __name__=='__main__':
    unittest.main()
#     TestUtils().testBeamMM_UBeamStraightRot()
#     TestUtils().testInertiaTranslate()
