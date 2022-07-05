import numpy as np
import pandas as pd
from scipy.linalg import expm

def EstimateKFTimeStep(u1,y1,z0,Xxd,Xud,Yx,Yu,P0,Q,R): 
    """ Performs one time step of Kalman filter estimation 

    INPUTS:
      u1: inputs at time n
      y1: measurements at time n
      z0: Kalman state estimate at time n-1
    OUTPUTS:
      z1: States at time n
      P1: Process covariance at time n  (nx x nx)
      Kk: Kalman gain

     Equations number are compared to the following reference:
     [1] Lourens"""
        
    # estimate next step
    z1m   = Xxd.dot(z0)  + Xud.dot(u1)
    y1hat = Yx.dot(z1m)  + Yu.dot(u1)
    P1m   = (Xxd.dot(P0)).dot(Xxd.T) + Q  # (nx x nx)
    
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
    Xx = np.asarray(Xx)
    Xu = np.asarray(Xu)
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


def BuildSystem_Linear_MechOnly(M, C, K, nP=0, nU=0, nY=0, Fp=None):
    """ 
    Takes mechanical system matrices, returns state matrices with only the mechanical part filled.
    The state matrix may be an "augmented matrix" (nP>0)
    The user will have to fill the entries related to:
      - the augmented states "p"
      - the inputs u
      - the outputs y

    Returns Xu, Yu, Yx as zeros! (to be filled by user)

    - Mechanical equation:
       M qddot + Cqdot + Kq = Fp.p + Fu.u 

    State/Output Equations
        xdot = Xx.x + Xu.u + wd
          y  = Yx.x + Yu.u + wn
    """
    nDOF = M.shape[0]
    Z    = np.zeros((nDOF,nDOF))
    Znnp = np.zeros((nDOF,nP  ))
    Znpn = np.zeros((nP  ,nDOF))
    I    = np.eye(nDOF)
    mM_K = np.linalg.solve(-M,K)
    mM_C = np.linalg.solve(-M,C)
    if Fp is not None:
        M_Fp  = np.linalg.solve(M,Fp)
    else:
        M_Fp = np.zeros((nDOF,nP))  # NOTE: to be filled by user
    Pp   = np.zeros((nP  ,nP))  # NOTE: to be filled by user
    Xx = np.block( [ [Z, I ,Znnp] , [mM_K, mM_C, M_Fp], [Znpn, Znpn, Pp] ])

    Xu = np.zeros((2*nDOF+nP,nU))# NOTE: to be filled by user
    Yx = np.zeros((nY,2*nDOF+nP))  # NOTE: to be filled by user
    Yu = np.zeros((nY,nU))       # NOTE: to be filled by user

    return Xx,Xu,Yx,Yu




def EmptyStateMat(nX,nU,nY):
    """ Returns state matrices with proper dimensions, filled with 0 """
    Xx = np.zeros((nX,nX)) # Ac 
    Yx = np.zeros((nY,nX)) # Gc 
    Xu = np.zeros((nX,nU)) # Xu 
    Yu = np.zeros((nY,nU)) # Jc 
    return Xx,Xu,Yx,Yu

def EmptyStateDF(nX,nU,nY,sX,sU,sY):
    """ Returns state matrices with proper dimensions, filled with 0 """
    Xx,Xu,Yx,Yu = EmptyStateMat(nX,nU,nY)
    Xx = pd.DataFrame(data=Xx, index=sX, columns=sX)
    Xu = pd.DataFrame(data=Xu, index=sX, columns=sU)
    Yx = pd.DataFrame(data=Yx, index=sY, columns=sX)
    Yu = pd.DataFrame(data=Yu, index=sY, columns=sU)
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


if __name__=='__main__':
    pass
