""" 
Performs simple BEM simulations of the NREL 5MW turbine for two operating conditions.
"""
# --- Common libraries 
from welib.BEM.MiniBEM import MiniBEM, FASTFile2MiniBEM
import os

MyDir=os.path.dirname(__file__)

if __name__=="__main__":
    """ See examples/ for more examples """

    # --- Read a FAST model to get the main parameters needed
    nB,cone,r,chord,twist,polars,rho,KinVisc = FASTFile2MiniBEM(os.path.join(MyDir,'../../../_data/NREL5MW/Main_Onshore_OF2.fst'))

    # --- Run BEM on a set of operating points
    WS =[5,10]
    RPM=[7,12]
    a0, ap0  = None,None  # inductions, used to speed up parametric study
    for i,ws in enumerate(WS):
        V0        = WS[i]
        Omega     = RPM[i]
        pitch=2     #[deg]
        xdot=0      #[m/s]
        u_turb=0    #[m/s]
        BEM=MiniBEM(Omega,pitch,V0,xdot,u_turb,
                    nB,cone,r,chord,twist,polars,
                    rho=rho,KinVisc=KinVisc,bTIDrag=False,bAIDrag=True,
                    a_init =a0,
                    ap_init=ap0
                    )
        a0, ap0 = BEM.a, BEM.aprime
        print('WS ',V0, 'Power',BEM.Power,'Thrust',BEM.Thrust)

        # --- Save radial distribution to a csv file
        filename='_BEM_ws{}_radial.csv'.format(V0)
        BEM.WriteRadialFile(filename)

