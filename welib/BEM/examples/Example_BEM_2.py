""" 
Performs simple BEM simulations of the NREL 5MW turbine for a set of operating conditions.
"""
# --- Common libraries 
import sys
import os
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
    MainFASTFile           = os.path.join(MyDir,'../../../data/NREL5MW/Main_Onshore_OF2.fst')
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
    return dfOut


if __name__=="__main__":
    main()
if __name__=="__test__":
    main(test=True)
