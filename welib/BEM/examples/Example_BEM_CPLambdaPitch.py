import os
import numpy as np
import pandas as pd
from welib.BEM.MiniBEM import MiniBEM, FASTFile2MiniBEM

MyDir=os.path.dirname(__file__)

if __name__=="__main__":
    """ Run a parametric study on lambda and pitch, save CP and CT """

    # --- Turbine data
    nB,cone,r,chord,twist,polars,rho,KinVisc = FASTFile2MiniBEM(os.path.join(MyDir,'../../../_data/NREL5MW/Main_Onshore_OF2.fst'))

    vlambda = np.linspace( 2,15,10)
    vpitch  = np.linspace(-10,40,20)
    print(vlambda)
    print(vpitch)

    V0=10
    R=r[-1]

    tilt = 6; # TODO
    V0=V0*np.cos(tilt*np.pi/180)

    CP=np.zeros((len(vlambda),len(vpitch)))
    CT=np.zeros((len(vlambda),len(vpitch)))

    i=0
    ntot=len(vlambda)*len(vpitch)
    for il,lambd in enumerate(vlambda):
        for ip,pitch in enumerate(vpitch):
            i+=1
            Omega=lambd*V0/R * 60/(2*np.pi)
            xdot=0      #[m/s]
            u_turb=0    #[m/s]
            a0,ap0 = None, None
            BEM=MiniBEM(Omega,pitch,V0,xdot,u_turb,
                        nB,cone,r,chord,twist,polars,
                        rho=rho,bTIDrag=True,bAIDrag=True,
                        a_init =a0, ap_init=ap0
                        )
            a0,ap0 = BEM.a, BEM.aprime
            print('{:5d}/{:d} lambda={:4.2f} pitch={:4.0f} CP={:5.3f}'.format(i,ntot,lambd,pitch,BEM.CP))
            CP[il,ip]=BEM.CP
            CT[il,ip]=BEM.CT
    CP[CP<0]=np.nan
    CT[CT<0]=np.nan

    np.savetxt('_CP.csv'    , CP     , delimiter=',')
    np.savetxt('_CT.csv'    , CT     , delimiter=',')
    np.savetxt('_Lambda.csv', vlambda, delimiter=',')
    np.savetxt('_Pitch.csv' , vpitch , delimiter=',')

