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
from welib.BEM.steadyBEM import SteadyBEM, FASTFile2SteadyBEM
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


    # -- Extracting information from a FAST main file and sub files
    nB,cone,r,chord,twist,polars,rho,KinVisc = FASTFile2SteadyBEM(MainFASTFile) 
    BladeData=np.column_stack((r,chord,twist))

    # --- Running BEM simulations for each operating conditions
    a0 , ap0 = None,None # Initial guess for inductions, to speed up BEM
    dfOut=None
    if not os.path.exists(OutDir):
        os.makedirs(OutDir)
    # Loop on operating conditions
    for i,(V0,RPM,pitch), in enumerate(zip(WS,Omega,Pitch)):
        xdot   = 0        # [m/s]
        u_turb = 0        # [m/s] 
        BEM=SteadyBEM(RPM,pitch,V0,xdot,u_turb,
                    nB,cone,r,chord,twist,polars,
                    rho=rho,KinVisc=KinVisc,bTIDrag=False,bAIDrag=True,
                    a_init =a0, ap_init=ap0)
        a0, ap0 = BEM.a, BEM.aprime
        # Export radial data to file
        if not test:
            filenameRadial = os.path.join(OutDir,'_BEM_ws{:02.0f}_radial.csv'.format(V0))
            BEM.WriteRadialFile(filenameRadial)
            print('>>>',filenameRadial)
        dfOut = BEM.StoreIntegratedValues(dfOut)

    # --- Export integrated values to file
    filenameOut= os.path.join(OutDir,'_BEM_IntegratedValues.csv')
    if not test:
        dfOut.to_csv(filenameOut,sep='\t',index=False)
        print('>>>',filenameOut)
        print(dfOut.keys())
        # --- Plot
        fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
        fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
        ax.plot(dfOut['WS_[m/s]'].values, dfOut['AeroPower_[kW]'].values, label='Power')
        ax.set_xlabel('Wind speed [m/s]')
        ax.set_ylabel('[-]')
        ax.legend()
        ax.set_title('BEM Steady - Performance curve')
        ax.tick_params(direction='in')


    return dfOut


if __name__=="__main__":
    main()
    plt.show()
if __name__=="__test__":
    main(test=True)
    [os.remove(f) for f in glob.glob(os.path.join(MyDir,'_*.csv'))]
if __name__=="__export__":
    main()
    from welib.tools.repo import export_figs_callback
    export_figs_callback(__file__)
