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
import welib.fast.fastlib as fastlib
from welib.weio.fast_output_file import FASTOutputFile
from welib.essentials import *
#import pyFAST.input_output as io

def powerCurvePostProWindStep(outfile, ColMap={}, tWindow=20, WS=None, KeepAll=False, sWS='WS_[m/s]', sWSRaw='Wind1VelX_[m/s]', sTime='Time_[s]', dataDict=None):
    def renameCol(x):
        for k,v in ColMap.items():
            if x==v:
                return k
        return x

    # --- Reading wind and output file
    print('Reading input...', outfile)
    df = FASTOutputFile(outfile).toDataFrame()
    #units=['('s.split('[')[1].split(']')[0] for s in df.columns.values
    # --- Our main values
    for v in [sWSRaw, sTime]:
        if v not in df.columns.values:
            raise Exception('Column `{}` need to be present in data. Please provide a column map to define it'.format(v))
    vWS = np.round(df[sWSRaw].values,2)
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
            #print('WS:{:4.2f} - t: {:6.1f}-{:6.1f} - nValues:{}'.format(ws,tStart,tEnd,len(IWindow)))
            MeanValues = pd.DataFrame(df[IWindow].mean()).transpose()
            if i==0:
                result = MeanValues.copy()
            else:
                #result=result.append(MeanValues, ignore_index=True)
                result = pd.concat([result, MeanValues], ignore_index=True)

    # --- Sorting by WS
    result[sWSRaw]=np.round(result[sWSRaw],2)
    result.sort_values([sWSRaw],inplace=True)
    result.reset_index(drop=True,inplace=True) 
    
    # --- Apply column map
    dataDict['U'] = result[sWSRaw].values
    result = fastlib.remap_df(result, ColMap, bColKeepNewOnly=False, dataDict=dataDict)
    
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
    TipRad=63
    PreCone=-2.5
    rho = 1.225
    prebendTip = -3.2815226E-04
    R_proj = TipRad*np.cos(PreCone*np.pi/180) - prebendTip*np.sin(PreCone*np.pi/180)
    A = np.pi*R_proj**2
    data={}
    data['rho'] = rho
    data['A'] = A
    data['R'] = R_proj
    # --- FAST OLD OPER
    ColMap = {}        
    ColMap['WS_[m/s]']         = 'Wind1VelX_[m/s]'
    ColMap['RotSpeed_[rpm]']   = 'RotSpeed_[rpm]'
    ColMap['BldPitch_[deg]']   = 'BldPitch1_[deg]'
    ColMap['GenSpeed_[rpm]']   = ''
    ColMap['GenPower_[kW]']    = 'GenPwr_[kW]'
    ColMap['GenTorque_[kN-m]'] = 'GenTq_[kN-m]'
    ColMap['RtAeroFxh_[kN]']   = '{RtAeroFxh_[N]}/1e3'
    ColMap['RtAeroMxh_[kN-m]'] = '{RtAeroMxh_[N-m]}/1e3'
    ColMap['RotThrust_[kN]']   = ''
    ColMap['RotTorq_[kN-m]']   = ''
    ColMap['TTDspFA_[m]']      = ''
    ColMap['TTDspSS_[m]']      = ''
    ColMap['RtAeroCp_[-]']      = ''
    ColMap['RtAeroCt_[-]']      = ''
    ColMap['RtTSR_[-]']      = ''

    
    # --- ROTORSE-LIKE OPER    
    ColMap = {}        
    # Trying to respect how WEIS Output operating conditions
    ColMap['Wind_[m/s]']        = 'Wind1VelX_[m/s]'
    ColMap['Pitch_[deg]']       = 'BldPitch1_[deg]'

    ColMap['Power_[MW]']        = '{GenPwr_[kW]}/1e3'
    ColMap['Aero_Power_[MW]']   = '{RtAeroPwr_[W]}/1e6'  
    ColMap['CP_[-]']            = '{Power_[MW]}     *1e6/(1/2*rho*A*U**3)'
    ColMap['Aero_CP_[-]']       = '{Aero_Power_[MW]}*1e6/(1/2*rho*A*U**3)'
    ColMap['RtAeroCp_[-]']      = ''

    ColMap['Rotor_Speed_[rpm]'] = 'RotSpeed_[rpm]'
    ColMap['GenSpeed_[rpm]']    = ''
    ColMap['Tip_Speed_[m/s]']   = '{RotSpeed_[rpm]}*2*np.pi/60*R'

    ColMap['Thrust_[MN]']       = '{RotThrust_[kN]}/1e3'
    ColMap['Aero_Thrust_[MN]']  = '{RtAeroFxh_[N]}/1e6'

    ColMap['CT_[-]']            = '{Thrust_[MN]}     *1e6/(1/2*rho*A*U**2)'
    ColMap['Aero_CT_[-]']       = '{Aero_Thrust_[MN]}*1e6/(1/2*rho*A*U**2)'
    ColMap['RtAeroCt_[-]']      = ''

    ColMap['Torque_[MNm]']      = '{RotTorq_[kN-m]}/1e3'
    ColMap['GenTorque_[MNm]']   = '{GenTq_[kN-m]}/1e3  '
    ColMap['Aero_Torque_[MNm]'] = '{RtAeroMxh_[N-m]}/1e6'

    ColMap['CQ_[-]']            = '{Torque_[MNm]}     *1e6/(1/2*rho*A*U**2*R)'
    ColMap['Aero_CQ_[-]']       = '{Aero_Torque_[MNm]}*1e6/(1/2*rho*A*U**2*R)'

    ColMap['Blade_Moment_[MNm]'] = '{RootMyc1_[kN-m]}/1e3'
    ColMap['Blade_Moment_Coefficient_[-]'] = '{RootMyc1_[kN-m]}*1e3/(1/2*rho*A*U**2*R)'

    ColMap['R_ref_[m]']          = '{ones}* R'
    ColMap['TSR_[-]']            = '{RotSpeed_[rpm]}*2*np.pi/60*R/U'
    ColMap['RtTSR_[-]']          = ''

    ColMap['TTDspFA_[m]']        = ''
    ColMap['TTDspSS_[m]']        = ''


    OutFiles =[]
    OutFiles +=['_PowerCurve_Onshore_Flexible_Discon.outb']
    OutFiles +=['_PowerCurve_Onshore_Flexible_ROSCO.outb']
    OutFiles +=['_PowerCurve_Onshore_Rigid_Discon.outb']
    OutFiles +=['_PowerCurve_Onshore_Rigid_ROSCO.outb']
    for OutFile in OutFiles:
        result = powerCurvePostProWindStep(OutFile, tWindow=20, ColMap=ColMap, sWS='Wind_[m/s]', dataDict=data)
        operFile = OutFile.replace('.outb','.csv').replace('.out','.csv')
        operFile = operFile.replace('PowerCurve','NREL5MW_Oper')
        print('Writing: ', operFile)        
        result.to_csv(operFile, sep=',', index=False)


if __name__=='__test__':
    # Need openfast.exe, doing nothing
    pass

if __name__=='__export__':
    pass
    #from welib.tools.repo import export_figs_callback
    #export_figs_callback(__file__)


