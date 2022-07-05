""" 
Runs a set of OpenFAST linearizations for different operating points
"""
import numpy as np
import os
import welib.weio as weio # https://github.com/ebranlard/weio
import welib.fast.fastlib as fastlib # latest fastlib is found in https://github.com/ebranlard/welib
from   welib.tools.clean_exceptions import *
# --- Parameters for this script
main_file   = 'Main_Onshore.fst'         # Main file in ref_dir, used as a template
ref_dir     = '../../_data/NREL5MW/'         # Folder where the fast input files are located (will be copied)
FAST_EXE    = '../../_data/NREL5MW/openfast2.2_x64s.exe' # Location of a FAST exe (and dll)
trim_points = '../../_data/NREL5MW_Oper.csv'
work_dir    = '../../_data/NREL5MW_Linearizations/'     # Output folder (will be created)

# --- Defining the parametric study  (list of dictionnaries with keys as FAST parameters)
Cases=weio.read(trim_points).toDataFrame()
Cases=Cases.iloc[:-6:1,:]
# Adding 0 wind
# Cases=Cases.append(Cases.iloc[0,:])
# Cases.iloc[-1,:] = 0
# Cases['WS_[m/s]'].loc[-1] = 0.0001
# Cases['BldPitch_[deg]'].loc[-1] = Cases['BldPitch_[deg]'].loc[0]

print(Cases)

Tmin       = 100
nPerPeriod = 6

PARAMS     = []
for i,c in Cases.iterrows():
    if np.mod(i,5)!=1: # skipping some for this example
        continue
    ws    = c['WS_[m/s]']
    rpm   = c['RotSpeed_[rpm]']
    pitch = c['BldPitch_[deg]']
    tt    = c['TTDspFA_[m]']
    Omega = rpm/60*2*np.pi
    if abs(Omega)<0.001:
        LinTimes = [Tmin]
        Tmax     = Tmin+1
    else:
        T = 2*np.pi/Omega
        LinTimes = np.linspace(Tmin,Tmin+T,nPerPeriod+1)[:-1];
        Tmax       = Tmin+1.01*T
    print('rpm',rpm,'Omega',Omega,'T',T)
        
    p=dict()
    p['__name__']       = 'ws_{:.0f}'.format(ws)
    p['TMax']         = Tmax
    p['TStart']       = 0
    p['Linearize']    = True
    p['NLinTimes']    = len(LinTimes)
    p['LinTimes']     = list(LinTimes)
    p['LinInputs']    = 2                   # <<<<<<<<<<<<<<<<<< To get full linearizations
    p['LinOutputs']   = 1                   # <<<<<<<<<<<<<<<<<< To get full linearizations
    p['OutFmt']       = '"ES15.8E2"'
    if abs(ws)<0.001:
        p['CompAero']    = 0
        p['CompInflow']  = 0
    p['InflowFile|WindType'] = 1
    p['InflowFile|HWindSpeed'] = ws
    p['InflowFile|PLexp']  = 0.0                 # <<<<<<<<<<<<<<<<<
    p['EDFile|BlPitch(1)'] = pitch
    p['EDFile|BlPitch(2)'] = pitch
    p['EDFile|BlPitch(3)'] = pitch
    p['EDFile|RotSpeed']   = rpm
    p['EDFile|TTDspFA']    = tt
    p['EDFile|FlapDOF1'] = False
    p['EDFile|FlapDOF2'] = False
    p['EDFile|EdgeDOF']  = False
    p['EDFile|TeetDOF']  = False
    p['EDFile|DrTrDOF']  = False
    p['EDFile|GenDOF']   = True
    p['EDFile|YawDOF']   = False
    p['EDFile|TwFADOF1'] = True
    p['EDFile|TwFADOF2'] = False
    p['EDFile|TwSSDOF1'] = False
    p['EDFile|TwSSDOF2'] = False
    #
    p['AeroFile|WakeMod'] = 1
    p['AeroFile|AFAeroMod'] = 1
    p['AeroFile|Frozenwake'] = True
    p['ServoFile'] = '"NREL5MW_SD_Simple.dat"'
    PARAMS.append(p)
# # --- Generating all files in a workdir
fastfiles=fastlib.templateReplace(PARAMS,ref_dir,outputDir=work_dir,removeRefSubFiles=True,main_file=main_file)
# 
# # --- Creating a batch script just in case
fastlib.writeBatch(os.path.join(work_dir,'_Linearization_RUN_ALL.bat'),fastfiles,fastExe=FAST_EXE)
# # # --- Running the simulations
fastlib.run_fastfiles(fastfiles,fastExe=FAST_EXE,parallel=True,showOutputs=True,nCores=2)
# 
# 
