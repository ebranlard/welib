""" 
Performs simple BEM simulations of the NREL 5MW turbine for two operating conditions.
"""
# --- Common libraries 
from welib.BEM.steadyBEM import SteadyBEM, FASTFile2SteadyBEM
import os
import glob
import matplotlib.pyplot as plt

MyDir=os.path.dirname(__file__)

def main(test=False):
    # --- Read a FAST model to get the main parameters needed
    nB,cone,r,chord,twist,polars,rho,KinVisc = FASTFile2SteadyBEM(os.path.join(MyDir,'../../../data/NREL5MW/Main_Onshore.fst'))

    # --- Run BEM on a set of operating points
    if test:
        WS    = [5,10]
        RPM   = [7,12]
        PITCH = [0,0]
    else:
        WS    = [4   ,6   ,8   ,10  ,12  ,14  ,16  ,18  ,20  ,22  ,24  ,26]
        RPM   = [2.38,7.92,9.13,11.3,12.1,12.1,12.1,12.1,12.1,12.0,12.0,12.1]
        PITCH = [0   ,0   ,0   ,0   ,3.4 ,8.4 ,11.8,14.7,17.3,19.6,21.9,23.9 ]

    a0, ap0  = None,None  # inductions, used to speed up parametric study
    for i,ws in enumerate(WS):
        V0     = WS[i]    # [m/s]
        Omega  = RPM[i]   # [rpm]
        pitch  = PITCH[i] # [deg]
        xdot   = 0        # [m/s]
        u_turb = 0        # [m/s]
        BEM=SteadyBEM(Omega,pitch,V0,xdot,u_turb,
                    nB,cone,r,chord,twist,polars,
                    rho=rho,KinVisc=KinVisc,bTIDrag=False,bAIDrag=True,
                    a_init =a0,
                    ap_init=ap0
                    )
        a0, ap0 = BEM.a, BEM.aprime

        # --- Save radial distribution to a csv file
        filename='_BEM_ws{}_radial.csv'.format(V0)
        if not test:
            print('WS ',V0, 'Power',BEM.Power,'Thrust',BEM.Thrust)
            BEM.WriteRadialFile(filename)

if __name__=="__main__":
    main()
if __name__=="__test__":
    main(test=True)
    [os.remove(f) for f in glob.glob(os.path.join(MyDir,'_*.csv'))]
if __name__=="__export__":
    pass
    #main()
    #from welib.tools.repo import export_figs_callback
    #export_figs_callback(__file__)


