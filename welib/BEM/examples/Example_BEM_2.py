""" 
Performs simple BEM simulations of the NREL 5MW turbine for a set of operating conditions.
"""
# --- Common libraries 
import sys
import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
# --- Local libraries
from welib.BEM.steadyBEM import SteadyBEM # option 1: using class
from welib.BEM.steadyBEM import calcSteadyBEM, FASTFile2SteadyBEM # option 2: low level
import welib.weio as weio

MyDir=os.path.dirname(__file__)

def main(test=False):
    """ 
    Performs BEM simulations for different Pitch, RPM and Wind Speed.
    The wind turbine is initialized using a FAST input file (.fst)
    """
    OutDir                 = os.path.join(MyDir,'./')
    MainFASTFile           = os.path.join(MyDir,'../../../data/NREL5MW/Main_Onshore.fst')
    OperatingConditionFile = os.path.join(MyDir,'../../../data/NREL5MW/NREL5MW_Oper.csv')

    #  --- Reading turbine operating conditions, Pitch, RPM, WS  (From FAST)
    df=weio.read(OperatingConditionFile).toDataFrame()
    Pitch = df['BldPitch_[deg]'].values
    Omega = df['RotSpeed_[rpm]'].values
    WS    = df['WS_[m/s]'].values

    filenameOut    = os.path.join(OutDir,'_BEM_IntegratedValues.csv') if not test else None
    radialBasename = os.path.join(OutDir,'_BEM_') if not test else None

    # --- Option 1: using wrapper class
    BEM = SteadyBEM(filename=MainFASTFile) # Initialize based on OpenFAST parameters
    # Change algo options
    BEM.bSwirl     = True  # swirl flow model enabled / disabled
    BEM.bTipLoss   = True  # enable / disable tip loss model
    BEM.bHubLoss   = False # enable / disable hub loss model
    BEM.bAIDrag    = True  # influence on drag coefficient on normal force coefficient
    BEM.bTIDrag    = False
    BEM.relaxation = 0.4
    BEM.bUseCm     = False  # Use Moment
    #print(BEM)
    # Compute performances for multiple operating conditions, store in dataframe, and export to file
    dfOut,_ = BEM.parametric(Omega, Pitch, WS, outputFilename=filenameOut, radialBasename=radialBasename)
    dfRad = BEM.radialDataFrame()

    # --- Option 2: using low level functions instead of class
    if False:
        # --- Extracting information from a FAST main file and sub files
        nB,cone,r,chord,twist,polars,rho,KinVisc = FASTFile2SteadyBEM(MainFASTFile) 
        # --- Running BEM simulations for each operating conditions
        a0 , ap0 = None,None # Initial guess for inductions, to speed up BEM
        dfOut=None
        # Loop on operating conditions
        for i,(V0,RPM,pitch), in enumerate(zip(WS,Omega,Pitch)):
            BEM=calcSteadyBEM(RPM,pitch,V0,xdot=0,u_turb=0,
                        nB=nB,cone=cone,r=r,chord=chord,twist=twist,polars=polars,
                        rho=rho,KinVisc=KinVisc,bTIDrag=False,bAIDrag=True,bUseCm=False,
                        a_init =a0, ap_init=ap0)
            a0, ap0 = out.a, out.aprime
            # Export radial data to file
            dfRad = BEM.radialDataFrame()
            dfOut = out.StoreIntegratedValues(dfOut)

    # --- Plot
    if not test:
        print(dfOut.keys())
        fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
        fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
        ax.plot(dfOut['WS_[m/s]'].values, dfOut['AeroPower_[kW]'].values, label='Power')
        ax.set_xlabel('Wind speed [m/s]')
        ax.set_ylabel('[-]')
        ax.legend()
        ax.set_title('BEM Steady - Performance curve')
        ax.tick_params(direction='in')



    return dfOut, dfRad


if __name__=="__main__":
    dfOut, dfRad = main(test=False)
    plt.show()
if __name__=="__test__":
    dfOut, dfRad = main(test=True)
    np.testing.assert_almost_equal(dfOut['AeroPower_[kW]'].values[-5], 5710.6344, 3)
    np.testing.assert_almost_equal(dfRad['a_prime_[-]'].values[4], 0.1636598, 4)
    [os.remove(f) for f in glob.glob(os.path.join(MyDir,'_*.csv'))]
if __name__=="__export__":
    main()
    from welib.tools.repo import export_figs_callback
    export_figs_callback(__file__)
