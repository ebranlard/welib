import numpy as np
import unittest
from scipy.linalg import expm

def EstimateKFTimeStep(u1,y1,z0,Xxd,Xud,Yx,Yu,P0,Q,R): 
    """ Performs one time step of Kalman filter estimation 

    INPUTS:
      u1: inputs at time n
      y1: measurements at time n
      z0: Kalman state estimate at time n-1
    OUTPUTS:
      z1: States at time n
      P1: Process covariance at time n
      Kk: Kalman gain

     Equations number are compared to the following reference:
     [1] Lourens"""
        
    # estimate next step
    z1m   = Xxd.dot(z0)  + Xud.dot(u1)
    y1hat = Yx.dot(z1m)  + Yu.dot(u1)
    P1m   = (Xxd.dot(P0)).dot(Xxd.T) + Q
    
    # Calculate Kalman gain
    # same as Lk from [1] - And their Rtilde_k is G*P1m*G'+R
    Kk =  np.dot(P1m,Yx.T).dot( np.linalg.inv(((Yx.dot(P1m)).dot(Yx.T) + R))) 
    # update estimate with measurement
    z1 = z1m + Kk.dot(y1 - y1hat)
    
    P1 = (np.eye(Xxd.shape[0]) - Kk.dot(Yx) ).dot(P1m)
    return z1,P1,Kk


def KFDiscretize(Xx,Xu,dt,method='exponential'):
    """ Discretize the continuous states matrices Xx, Xu
    
    "Real" system:
        zdot = Xx.x + Xu.u + wd
        zdot = Ac.z + Bc.u + wd
          y  = Gc.z + Jc.u + wn
    "Discretized" system:
        z_{k+1} = Xxd.z_k + Xud.u_k + wd
        z_{k+1} = Ad.z_k  + Bd.u_k + wd
          y_k   = Gd.z_k  + Jd.u_k + wn 

    INPUTS:
        methods: 'exponential', 'eigenvalues', 'forward_euler'

    OUTPUTS:
       Xxd,Xud, disrete versions of Xx and Xu
    """
    # --- A matrix
    if method=='exponential':
        # Using matrix exponential directly
        Xxd = expm(Xx * dt)
        if np.linalg.det(Xx) == 0:
            print('[WARN] Matrix A is singular, using forward euler to discretize B matrix\n' % ())
            # Forward euler
            Xud = dt * Xu
        else:
            mA_B = np.linalg.solve(Xx,Xu)
            Xud = np.dot( (Xx - np.eye(Xx.shape[0])) ,  mA_B)
    elif method=='eigenvalues':
        raise NotImplementedError()
        # Using eigenvalues
        #Q,Lambda = eig(Xx) # Need my version of eigenvalue, so postponed to keep the library standalone
        #Xxd = real(Q * expm(Lambda * dt) / Q)
        #Xud = Xu * dt
    elif method=='forward_euler':
        # Using forward Euler
        Xxd = np.eye(Xx.shape[0]) + Xx * dt 
        Xud = Xu * dt
    else:
        raise Exception('Unknown discretization method: %s',method)

    return Xxd,Xud

def BuildSystem_Linear(M,C,K,Ya,Yv,Yq,Fp=None,Pp=None,Yp=None,Yu=None,Method='default'):
    """ Takes system matrices of a mechanical system, returns a state matrix.
    The state matrix may be an "augmented matrix", in which case Fp, Pp, should be provided

    - Mechanical equation:
       M qddot + Cqdot + Kq = Fp.p + Fu.u 

    - Output equation:
       y = Ya.qddot + Yv.qdot + Yq.q  + Yp.p  + ~Yu.u
       y = Sa.qddot + Sv.qdot + Sd.q  + Yp.p  + ~Yu.u
    - (Augmented load evolution:
         pdot = Pp.p  + Pq.q + Pv.qdot

    State Equation
        xdot = Xx.x + Xu.u + wd
        zdot = Ac.z + Bc.u + wd
    Measurement Equation
          y  = Yx.x + Yu.u + wn
          y  = Gc.z + Jc.u + wn
    """
    nDOF = M.shape[0]
    nY   = Yq.shape[0]
    if Yu is None:
        nU = 0
        Yu = np.zeros((nY,nU))
    else:
        nU = Yu.shape[1]

    if Method=='default':
        Z=np.zeros((nDOF,nDOF))
        I=np.eye(nDOF)
        Xx = np.block( [ [Z , I ], [ mM_K, mM_C] ])
        Xu = np.zeros((2*nDOF,nU))
        Yx = np.block( [ Yq + np.dot(Ya,mM_K),  Yv + np.dot(Ya, mM_C) ] )
    elif Method == 'augmented_first_order':
        # Needs Fp and Pp to be defined!
        if Fp is None or Pp is None:
            raise Exception('Both Fp and Pp needs to be set with augmented first order method')
        nP = Fp.shape[1]
        if Yp is None:
            Yp=np.zeros((nY,nP))

        Z    = np.zeros((nDOF,nDOF))
        Znnp = np.zeros((nDOF,nP  ))
        Znpn = np.zeros((nP  ,nDOF))
        I    = np.eye(nDOF)
        mM_K = np.linalg.solve(-M,K)
        mM_C = np.linalg.solve(-M,C)
        M_Fp  = np.linalg.solve(M,Fp)
        Xx = np.block( [ [Z, I ,Znnp] , [mM_K, mM_C, M_Fp], [Znpn, Znpn, Pp] ])
        Xu = np.zeros((2*nDOF+nP,nU))
        Yx = np.block( [Yq + np.dot(Ya,mM_K), Yv + np.dot(Ya,mM_C), Yp+np.dot(Ya,M_Fp) ])
#         print('Yq..:\n', Yq + np.dot(Ya,mM_K))
#         print('Yv..:\n', Yv + np.dot(Ya,mM_C))
#         print('Fp..:\n', Yp+np.dot(Ya,M_Fp) )
    else:
        raise Exception('Method %s not implemented')
    
    return Xx,Xu,Yx,Yu



def EmptyStateMat(nX,nU,nY):
    """ Returns state matrices with proper dimensions, filled with 0 """
    Xx = np.zeros((nX,nX)) # Ac 
    Yx = np.zeros((nY,nX)) # Gc 
    Xu = np.zeros((nX,nU)) # Xu 
    Yu = np.zeros((nY,nU)) # Jc 
    return Xx,Xu,Yx,Yu

def EmptySystemMat(nDOF_2nd, nY, nP=None, nU=None):
    """ Returns matrices with proper dimensions, filled with 0
    INPUTS:
       - nDOF_2nd: number of "mechanical" degrees of freedoms, when the equations are written 
                   as a secnod order system, e.g. M qddot = F, then nDOF_2nd is the size of qddot
       - nY: Number of outputs
       - nP: Number of extended loads, if states are to be augmented
       - nU: Number of inputs
    NOTE:
        full augmented state vector has size 2*nDOF_2nd + nP
    """
    M=np.zeros((nDOF_2nd,nDOF_2nd))
    C=np.zeros((nDOF_2nd,nDOF_2nd))
    K=np.zeros((nDOF_2nd,nDOF_2nd))
    Ya = np.zeros((nY,nDOF_2nd)) # Sa
    Yv = np.zeros((nY,nDOF_2nd)) # Sv
    Yq = np.zeros((nY,nDOF_2nd)) # Sd

    if (nU is not None and nP is None):
        nP=0
    elif (nU is None and nP is not None):
        nU=0


    if nP is not None:
        Yp = np.zeros((nY,nP))       # Yp
        Fp = np.zeros((nDOF_2nd,nP)) # Sp
        Pp = np.zeros((nP,nP))       # Rp
        Pv = np.zeros((nDOF_2nd,nP)) 
        Pq = np.zeros((nDOF_2nd,nP)) 

    if nU is not None:
        Yu = np.zeros((nY,nU))       
        Fu = np.zeros((nDOF_2nd,nU)) 
        Pu = np.zeros((nP,nU))       

    if (nU is None) and (nP is None):
        return M,C,K,Ya,Yv,Yq
    else:
        return M,C,K,Ya,Yv,Yq,Yp,Yu,Fp,Fu,Pp,Pq,Pv



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
