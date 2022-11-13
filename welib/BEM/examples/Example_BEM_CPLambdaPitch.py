import os
import numpy as np
import pandas as pd
from welib.BEM.steadyBEM import calcSteadyBEM, FASTFile2SteadyBEM
from welib.weio.rosco_performance_file import ROSCOPerformanceFile
import matplotlib.pyplot as plt

MyDir=os.path.dirname(__file__)


def main(test=False,extra=False):
    """ Run a parametric study on lambda and pitch, save CP and CT """
    exportFmt='rosco' # Export format, rosco or csv
    # --- Turbine data
    nB,cone,r,chord,twist,polars,rho,KinVisc = FASTFile2SteadyBEM(os.path.join(MyDir,'../../../data/NREL5MW/Main_Onshore.fst'))

    if test:
        vlambda = np.linspace( 2,15,3)
        vpitch  = np.linspace(-10,40,2)
    elif extra:
        vlambda = np.linspace( 2 ,16,41)
        vpitch  = np.linspace(-10,40,40)
    else:
        vlambda = np.linspace( 2 ,16,10)
        vpitch  = np.linspace(-10,40,20)

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
            BEM=calcSteadyBEM(Omega,pitch,V0,xdot,u_turb,
                        nB,cone,r,chord,twist,polars,
                        rho=rho,bTIDrag=True,bAIDrag=True,
                        a_init =a0, ap_init=ap0
                        )
            a0,ap0 = BEM.a, BEM.aprime
            print('{:5d}/{:d} lambda={:4.2f} pitch={:4.0f} CP={:5.3f}'.format(i,ntot,lambd,pitch,BEM.CP))
            CP[il,ip]=BEM.CP
            CT[il,ip]=BEM.CT
    CP[CP<0]=0
    CT[CT<0]=0

    # --- Create a ROSCO performance file
    rs = ROSCOPerformanceFile(pitch=vpitch, tsr=vlambda, CP=CP, CT=CT)

    # --- Export performances to a given file format
    if exportFmt.lower()=='csv':
        np.savetxt('_CP.csv'    , CP     , delimiter=',')
        np.savetxt('_CT.csv'    , CT     , delimiter=',')
        np.savetxt('_Lambda.csv', vlambda, delimiter=',')
        np.savetxt('_Pitch.csv' , vpitch , delimiter=',')
    elif exportFmt.lower()=='rosco':
        rs.write('_CPCTCQ.txt')

    # --- Plotting matrix of CP values
    try:
        fig = rs.plotCP3D()
        ax = fig.gca()
        ax.set_title('BEM Steady - CP-lambda-pitch ')
    except:
        print('Plot 3D failed. TODO')
    #fig.colorbar(surf, shrink=0.5, aspect=15)




if __name__=="__main__":
    main()
    plt.show()
if __name__=="__test__":
    main(True)
    try:
        os.remove('_CPCTCQ.txt')
    except:
        pass
if __name__=="__export__":
    main(extra=True)
    from welib.tools.repo import export_figs_callback
    export_figs_callback(__file__)

