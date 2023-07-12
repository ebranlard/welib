""" 

NOTE: This script needs polishing and should be moved to fastlib directly

Open an Output file from openfast, this output file should have been generated using a wind file with wind step. 
Wind steps are convenient to generate a "power curve" using ne simulation only.

"""
import os
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import sys

import welib.weio as io
from welib.weio.fast_output_file import FASTOutputFile
#import pyFAST.input_output as io

def powerCurvePostProWindStep(outfile, ColMap={}, tWindow=20, WS=None, KeepAll=False, sWS='WS_[m/s]', sTime='Time_[s]'):
    def renameCol(x):
        for k,v in ColMap.items():
            if x==v:
                return k
        return x

    # --- Reading wind and output file
    print('Reading input...')
    df = FASTOutputFile(outfile).toDataFrame()
    #print(df.columns.values)
    # --- Applying column map
    if len(ColMap)>0:
        for k,v in ColMap.items():
            if v=='':
                v=k
            if v not in df.columns.values:
                raise Exception('Column `{}` not found in output file, cannot apply column mapping'.format(v))
        df.rename(columns=renameCol,inplace=True)
    #units=['('s.split('[')[1].split(']')[0] for s in df.columns.values
    # --- Our main values
    for v in [sWS, sTime]:
        if v not in df.columns.values:
            raise Exception('Column `{}` need to be present in data. Please provide a column map to define it'.format(v))
    vWS = np.round(df[sWS].values,2)
    t   = df[sTime].values
    # --- Finding unique WS values if not provided
    if WS is None:
        WS=np.unique(np.round(vWS,2))
    # --- Looping on wind values, and averaging data over a window
    tol =0.01
    Iall = np.arange(len(t))
    result=None
    for i,ws in enumerate(WS):
        b = abs(vWS-ws)<tol
        I = Iall[b]
        if len(I)>30:
            iEnd = I[-1]
            tEnd = t[iEnd]
            tStart = t[iEnd]-tWindow
            IWindow=(t>tStart) & (t<tEnd)
            print('WS:{:4.2f} - t: {:6.1f}-{:6.1f} - nValues:{}'.format(ws,tStart,tEnd,len(IWindow)))
            MeanValues = pd.DataFrame(df[IWindow].mean()).transpose()
            if i==0:
                result = MeanValues.copy()
            else:
                result = pd.concat((result, MeanValues), axis=0, ignore_index=True)

    # --- Sorting by WS
    result[sWS]=np.round(result[sWS],2)
    result.sort_values([sWS],inplace=True)
    result.reset_index(drop=True,inplace=True) 
    
    # --- Eliminating unnecessary columns
    if len(ColMap)>0 and (not KeepAll):
        ColNames = list(ColMap.keys())
    else:
        ColNames = result.columns.values
    result = result[ColNames]
    
    # --- Rounding
    result=result.round(4)
    return result 

if __name__=='__main__':
    # --- FAST
    ColMap = {}        
    ColMap['WS_[m/s]']         = 'Wind1VelX_[m/s]'
    ColMap['RotSpeed_[rpm]']   = 'RotSpeed_[rpm]'
    ColMap['BldPitch_[deg]']   = 'BldPitch1_[deg]'
    ColMap['GenSpeed_[rpm]']   = ''
    ColMap['GenPower_[kW]']    = 'GenPwr_[kW]'
    ColMap['GenTorque_[kN-m]'] = 'GenTq_[kN-m]'
    ColMap['RtAeroFxh_[kN]']   = 'RtAeroFxh_[N]'   # <<<<<<<<<<<<<<< 1000
    ColMap['RtAeroMxh_[kN-m]'] = 'RtAeroMxh_[N-m]' # <<<<<<<<<<<<<<< 1000
    ColMap['RotThrust_[kN]']   = ''
    ColMap['RotTorq_[kN-m]']   = ''
    ColMap['TTDspFA_[m]']      = ''
    ColMap['TTDspSS_[m]']      = ''
    ColMap['RtAeroCp_[-]']      = ''
    ColMap['RtAeroCt_[-]']      = ''
    ColMap['RtTSR_[-]']      = ''

    OutFiles =[]
    OutFiles +=['PowerCurve_Onshore_Flexible_ROSCO.outb']
    OutFiles +=['PowerCurve_Onshore_Rigid_ROSCO.outb']
    OutFiles +=['PowerCurve_Onshore_Flexible_Discon.outb']
    OutFiles +=['PowerCurve_Onshore_Rigid_Discon.outb']
    for OutFile in OutFiles:
        result = powerCurvePostProWindStep(OutFile, tWindow=20, ColMap=ColMap)
        result['RtAeroFxh_[kN]']    /=1000
        result['RtAeroMxh_[kN-m]']  /=1000
        operFile = OutFile.replace('.outb','.csv').replace('.out','.csv')
        operFile = operFile.replace('PowerCurve','NREL5MW_Oper')
        result.to_csv(operFile, sep='\t', index=False)


if __name__=='__test__':
    # Need openfast.exe, doing nothing
    pass

if __name__=='__export__':
    pass
    #from welib.tools.repo import export_figs_callback
    #export_figs_callback(__file__)


