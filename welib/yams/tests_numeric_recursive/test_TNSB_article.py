""" 
Test that builds the 5 DOF model presented in the article:
     [1]: Branlard, Flexible multibody dynamics using joint coordinates and the Rayleigh-Ritz approximation: the general framework behind and beyond Flex, Wind Energy, 2019

""" 
##
import numpy as np
import copy
import unittest
from welib.yams.bodies import FlexibleBody
from welib.yams.yams import *
from welib.yams.TNSB import manual_assembly

def main(DEBUG=False,main_axis='x',nShapes_twr=1,bInit=1):

    # Main Parameters
    nSpan_twr   = 101
    nSpan_bld   = 61
    nShapes_bld = 1 # 0,1,2

    bCompat =False

    bSftMass = 1
    bHubMass = 1
    bNacMass = 1
    bBldMass = 1
    nB       = 2 # 2 or 3
    main_axis =main_axis

    nDOF = 1 + nShapes_twr + nShapes_bld * nB

    q = np.zeros((nDOF,1))
    if bInit:
        q[:] = 1 # Define some kind of initial conditions

    ## --- Strucural and geometrical Inputs
    L_twr   = 100
    EI_twr  = 2*10**12
    m_twr   = 9*10**3
    L_bld   = 60
    EI_bld  = 2*10**10
    m_bld   = 5*10**2
    GKt_bld = 7*10**11
    jxx_bld = 10*5

    r_ET_inE    = np.array([[0]    ,[0],[0]]  )
    if main_axis=='x':
        r_TN_inT    = np.array([[L_twr],[0],[0]]  )
        r_NGnac_inN = np.array([[0]    ,[0],[2.0]])
        r_NS_inN    = np.array([[0]    ,[0],[-10]])
    elif main_axis=='z':
        r_TN_inT    = np.array([[0],[0],[L_twr]] )
        r_NGnac_inN = np.array([[1.0],[0],[0]])
        r_NS_inN    = np.array([[-10],[0],[0]])

    r_SGhub_inS = np.array([[0]    ,[0],[0]]  )
    r_SR_inS    = np.array([[0]    ,[0],[0]]  )
    r_RGhub_inS = np.array([[0]    ,[0],[0]]  )

    M_hub=10**5
    IR_hub = np.zeros((3,3))
    IR_hub[0,0] = 2*10**5
    IR_hub[1,1] = 2*10**5
    IR_hub[2,2] = 3*10**5 
    IR_hub = IR_hub 

    M_nac   = 4*10**5 
    I0_nac=np.zeros((3,3)) 
    I0_nac[0,0]=7*10**6
    I0_nac[1,1]=3*10**6
    I0_nac[2,2]=1*10**6 
    # Inertias not at COG...
    IG_hub = translateInertiaMatrix(IR_hub, M_hub, np.array([0,0,0]), r_RGhub_inS)
    IG_nac = translateInertiaMatrixToCOG(I0_nac,M_nac, r_NGnac_inN)

    ## Derived parameters
    iPsi = nShapes_twr # Index of DOF corresponding to azimuth
    # --------------------------------------------------------------------------------}
    ## --- Creating bodies
    # --------------------------------------------------------------------------------{
    Yaw=RigidBody('YawBearing',0,(0,0,0),(0,0,0));
    # Bld
    # TODO
    # TODO - THIS HAS SOME INITIAL CONDITION IN IT
    #Bld=UniformBeamBody('Blade', nShapes_bld, nSpan_bld, L_bld, EI_bld , m_bld, Mtop=0, jxxG=jxx_bld, GKt=GKt_bld, bCompatibility=bCompat)
    Blds=[]
    Blds.append(Body('B1'))
    #Blds[0].MM = np.array([
    # [  3.0000E+04,   0.0000E+00,   0.0000E+00,   0.0000E+00,   5.2444E+03,   0.0000E+00,  -2.4905E+02,  -1.1333E+03],
    # [  0.0000E+00,   3.0000E+04,   0.0000E+00,  -5.2401E+03,   0.0000E+00,   9.0000E+05,   0.0000E+00,   0.0000E+00],
    # [  0.0000E+00,   0.0000E+00,   3.0000E+04,   0.0000E+00,  -9.0000E+05,   0.0000E+00,   1.1746E+04,  -6.5057E+03],
    # [  0.0000E+00,  -5.2401E+03,   0.0000E+00,   6.0150E+06,   0.0000E+00,  -4.3043E+05,   0.0000E+00,   0.0000E+00],
    # [  5.2444E+03,   0.0000E+00,  -9.0000E+05,   0.0000E+00,   3.6015E+07,   0.0000E+00,  -5.1196E+05,   8.1533E+04],
    # [  0.0000E+00,   9.0000E+05,   0.0000E+00,  -4.3043E+05,   0.0000E+00,   3.6000E+07,   0.0000E+00,   0.0000E+00],
    # [ -2.4905E+02,   0.0000E+00,   1.1746E+04,   0.0000E+00,  -5.1196E+05,   0.0000E+00,   7.5019E+03,   4.2759E+00],
    # [ -1.1333E+03,   0.0000E+00,  -6.5057E+03,   0.0000E+00,   8.1533E+04,   0.0000E+00,   4.2759E+00,   7.5066E+03]])
    Blds[0].MM = np.array([
                 [ 3.0000E+04, 0.0000E+00, 0.0000E+00, 0.0000E+00, 1.1741E+04, 0.0000E+00,-1.9634E+02],
                 [ 0.0000E+00, 3.0000E+04, 0.0000E+00,-1.1746E+04, 0.0000E+00, 9.0000E+05, 0.0000E+00],
                 [ 0.0000E+00, 0.0000E+00, 3.0000E+04, 0.0000E+00,-9.0000E+05, 0.0000E+00, 1.1746E+04],
                 [ 0.0000E+00,-1.1746E+04, 0.0000E+00, 6.0075E+06, 0.0000E+00,-5.1196E+05, 0.0000E+00],
                 [ 1.1741E+04, 0.0000E+00,-9.0000E+05, 0.0000E+00, 3.6008E+07, 0.0000E+00,-5.1196E+05],
                 [ 0.0000E+00, 9.0000E+05, 0.0000E+00,-5.1196E+05, 0.0000E+00, 3.6000E+07, 0.0000E+00],
                 [-1.9634E+02, 0.0000E+00, 1.1746E+04, 0.0000E+00,-5.1196E+05, 0.0000E+00, 7.5019E+03]])


    Blds[0].PhiU =np.zeros(nShapes_bld) # cannot access nf
    Blds[0].KK = np.zeros((8, 8))
    Blds[0].KK[6,6:]= np.array([ 2.8624E+05, -1.0224E+03])
    Blds[0].KK[7,6:]= np.array([-1.0224E+03,  1.1249E+07])
    Blds[0].MM = Blds[0].MM[:6+nShapes_bld,:6+nShapes_bld]
    Blds[0].KK = Blds[0].KK[:6+nShapes_bld,:6+nShapes_bld]
    Blds[0].MM *=bBldMass
    Blds[0].DD = np.zeros((6+nShapes_bld, 6+nShapes_bld))
    for iB in range(nB-1):
        Blds.append(copy.deepcopy(Blds[0]))

    # Generator only
    Gen=RigidBody('Gen', 0, IG_hub, r_SGhub_inS)
    # ShaftHub Body 
    Sft=RigidBody('ShaftHubGen',M_hub,IG_hub,r_SGhub_inS);
    Sft.MM*=bSftMass
    # Nacelle Body
    Nac=RigidBody('Nacelle',M_nac,IG_nac,r_NGnac_inN);
    Nac.MM*=bNacMass
    # Tower Body
    # TODO
    # TODO - THIS HAS SOME INITIAL CONDITION IN IT
    Mtop=sum([B.mass for B in Blds]) + Sft.mass + Nac.mass;
    Twr=UniformBeamBody('Tower', nShapes_twr, nSpan_twr, L_twr, EI_twr , m_twr, Mtop=Mtop, bAxialCorr=False, bStiffening=False, main_axis=main_axis, gravity=0)
    #  Temporary
    x_0=np.array([[0],[0],[0]])
    R_0b=np.eye(3)
    gz=q[0:nShapes_twr,0]
    v_0   = np.zeros(6+nShapes_twr)
    a_v_0 = np.zeros(6+nShapes_twr)
    # TODO: to fully match matlab code, need "UseShapeIntegral" implemented
    Twr.updateKinematics(x_0,R_0b,gz,v_0,a_v_0)
    Twr.computeMassMatrix()
    # --------------------------------------------------------------------------------}
    # --- Manual assembly 
    # --------------------------------------------------------------------------------{
    Struct = manual_assembly(Twr,Yaw,Nac,Gen,Sft,Blds,q,r_ET_inE,r_TN_inT,r_NS_inN,r_SR_inS,main_axis=main_axis,DEBUG=DEBUG)
    return Struct


class TestTNSB(unittest.TestCase):
    def test_TNSB_article(self):
        Struct=main()
        MM=Struct.MM
        KK=Struct.KK
        np.testing.assert_allclose(MM[0,0],7.86e5 ,rtol  = 1e-3)
        np.testing.assert_allclose(MM[1,1],7.23e7 ,rtol  = 1e-3)
        np.testing.assert_allclose(MM[2,2],7.50e3 ,rtol  = 1e-3)
        np.testing.assert_allclose(MM[3,3],7.50e3 ,rtol  = 1e-3)
        np.testing.assert_allclose(MM[0,2],7.71e3 ,rtol  = 1e-3)
        np.testing.assert_allclose(MM[0,3],1.578e4,rtol = 1e-3 )
        np.testing.assert_allclose(KK[0,0],6.01e6 ,rtol  = 1e-3)
        np.testing.assert_allclose(KK[1,1],0.00e0 ,rtol  = 1e-3)
        np.testing.assert_allclose(KK[2,2],2.86e5 ,rtol  = 1e-3)
        np.testing.assert_allclose(KK[2,2],2.86e5 ,rtol  = 1e-3)

        Twr=Struct.Twr
        Twr.gravity=9.81
        Twr.bStiffening=True
        Twr.computeStiffnessMatrix()
        np.testing.assert_allclose(Twr.Mtop,560000.0 ,rtol  = 1e-10)
        np.testing.assert_allclose(Twr.KKg[6,6],-98815.131096 ,rtol  = 1e-10)
        #print(Twr.Mtop)
        #print(Twr.KKg)


if __name__=='__main__':
    np.set_printoptions(linewidth=500)
    unittest.main()
    #Struct=main(DEBUG=True,main_axis='x',nShapes_twr=1,bInit=1)
