import unittest
import numpy as np
from welib.kalman.kalman import *

# --------------------------------------------------------------------------------}
# --- TESTS
# --------------------------------------------------------------------------------{
class Test(unittest.TestCase):

    def test_discretize_exp(self):
        # Discretize a continuous system using the exponential matrix method
        nX=3
        nU=1
        nY=2
        Xx, Xu, Yx, Yu = EmptyStateMat(nX,nU,nY)
        Xx[0,0]=1
        Xx[0,1]=4
        Xx[0,2]=4
        Xx[1,1]=2
        Xx[1,2]=5
        Xx[2,2]=3
        Xu[0,0]=1/2
        Xu[2,0]=1
        dt=0.5
        Xxd,Xud = KFDiscretize(Xx,Xu,dt,method='exponential')
        #print('Xx\n',Xxd)
        #print('Xu\n',Xud)
        Xxd_ref=np.array([
               [1.64872, 4.27824,12.60440],
               [0.00000, 2.71828, 8.81704],
               [0.00000, 0.00000, 4.48169]])
        Xud_ref = np.array([ [-2.00000], [ 0.83333], [ 0.66667]])
        np.testing.assert_almost_equal(Xxd,Xxd_ref,decimal=5)
        np.testing.assert_almost_equal(Xud,Xud_ref,decimal=5)

    def test_discretize_forward(self):
        # Discretize a continuous system using the forward euler method
        nX=3
        nU=1
        nY=2
        Xx, Xu, Yx, Yu = EmptyStateMat(nX,nU,nY)
        Xx[0,0]=0
        Xx[0,1]=1
        Xx[0,2]=0
        Xx[1,1]=0
        Xx[1,2]=1
        Xx[2,2]=0
        Xu[0,0]=1/2
        Xu[2,0]=1
        dt=0.5
        Xxd_xp,Xud_xp = KFDiscretize(Xx,Xu,dt,method = 'exponential'  )
        Xxd_fe,Xud_fe = KFDiscretize(Xx,Xu,dt,method = 'forward_euler')
        np.testing.assert_almost_equal(Xud_fe,Xud_xp,decimal=5)

        # TODO TODO, what to do of A
        #np.testing.assert_almost_equal(Xxd_fe,Xxd_xp,decimal=5)



    def test_build_shaftonly(self):
        # Build continuous matrices for the "Shaft only" case
        np.set_printoptions(linewidth=500)
        nDOF_2nd = 1           # Mech DOFs     :  q  = [psi]
        nDOF_1st = 2*nDOF_2nd  # Mech state    :  x  = [psi,psi_dot]
        nP       = 2           # Extended loads:  p  = [Ma ,Mg]
        nDOF_ext = nDOF_1st+nP # Ext. state    : x~  = [psi,psi_dot ,Ma,Mg]
        nY       = 2           # Outputs       :  y  = [psi_dot,Mg]
        nU       = 0           # Inputs
        J_LSS = 2.0
        # General init
        M,C,K,Ya,Yv,Yq,Yp,Yu,Fp,Fu,Pp,Pq,Pv = EmptySystemMat (nDOF_2nd, nY, nP, nU)
        # Setting values
        M[0,0]  = J_LSS
        Yv[0,0] = 1
        Fp[0,0] = 1
        Fp[0,1] = -1
        Yp[1,1] = 1 # Direct feed through of force

        Xx,Xu,Yx,Yu = BuildSystem_Linear(M,C,K,Ya,Yv,Yq,Fp=Fp,Pp=Pp,Yp=Yp,Method='augmented_first_order')
        # Reference values for test
        Xx_ref, Xu_ref, Yx_ref, Yu_ref = EmptyStateMat(nDOF_ext,nU,nY)
        Xx_ref[0,1] = 1
        Xx_ref[1,2] = 1/J_LSS
        Xx_ref[1,3] = -1/J_LSS
        Xx_ref[1,3] = -1/J_LSS
        Yx_ref[0,1] = 1
        Yx_ref[1,3] = 1
        np.testing.assert_equal(Xx,Xx_ref)
        np.testing.assert_equal(Xu,Xu_ref)
        np.testing.assert_equal(Yx,Yx_ref)
        np.testing.assert_equal(Yu,Yu_ref)

    def test_build_tower1shaft(self):
        # Build continuous matrices for the "Tower (1 mode) + Shaft " case
        nDOF_2nd = 2           # Mech DOFs     :  q  = [u, psi]
        nDOF_1st = 2*nDOF_2nd  # Mech state    :  x  = [u, psi, u_dot, psi_dot]
        nP       = 3           # Extended loads:  p  = [T, Ma ,Mg]
        nDOF_ext = nDOF_1st+nP # Ext. state    : x~  = [u, psi, udot, psi_dot , T, Qa,Qg]
        nY       = 3           # Outputs       :  y  = [u_ddot,psi_dot,Mg]
        nU       = 0           # Inputs
        M_twr = 3.0
        K_twr = 5.0
        C_twr = 10.0
        J_LSS = 2.0
        # General init of zero matrices
        M,C,K,Ya,Yv,Yq,Yp,Yu,Fp,Fu,Pp,Pq,Pv = EmptySystemMat (nDOF_2nd, nY, nP, nU)
        # Setting values
        M[0,0]  = M_twr
        M[1,1]  = J_LSS
        K[0,0]  = K_twr
        C[0,0]  = C_twr
        Ya[0,0] = 1  # uddot = qddot[0]
        Yv[1,1] = 1  # psidot = qdot[1]
        Yp[2,2] = 1  # Direct feed-through of Mg
        Fp[0,0] = 1  # T = p[0]
        Fp[1,1] = 1  # dQ = p[1] -p[2]
        Fp[1,2] = -1 # dQ = p[1] -p[2]

        Xx,Xu,Yx,Yu = BuildSystem_Linear(M,C,K,Ya,Yv,Yq,Fp=Fp,Pp=Pp,Yp=Yp,Method='augmented_first_order')

        # Reference values for test, NOTE: only true because no couplings assumed here
        Xx_ref, Xu_ref, Yx_ref, Yu_ref = EmptyStateMat(nDOF_ext,nU,nY)
        Yx_ref[0,0] = -K_twr/M_twr
        Yx_ref[0,2] = -C_twr/M_twr
        Yx_ref[0,4] =  1/M_twr
        Yx_ref[1,3] =  1           # psi_dot = x~[3]
        Yx_ref[2,6] =  1           # Mq      = x~[6]

        np.testing.assert_almost_equal(Yx,Yx_ref)

#         print('Xx\n',Xx)
#         print('Xu\n',Xu)
#         print('Yx\n',Yx)
#         print('Yu\n',Yu)
if __name__=='__main__':
    unittest.main()
