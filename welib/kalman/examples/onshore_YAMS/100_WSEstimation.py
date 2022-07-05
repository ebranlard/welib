import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import welib.weio as weio # https://github.com/ebranlard/weio

# Local 
from wtDigiTwin.ws_estimator.tabulated import TabulatedWSEstimator

IStudy=[1,3]


if 1 in IStudy:
    try:
        os.mkdir('../../_data/wsest/')
    except:
        pass
    relaxation=0.3
    base='../../_data/NREL5MW'
    Main='../../_data/NREL5MW_SimpleSimulations/TurbWindStep_AllDOF.fst'; 
    # --- Turbine data and estimator
    wse = TabulatedWSEstimator(fst_file=Main)
    wse.load_files(base=base,suffix='')
    print(wse)
    # --- Loading data
    df = weio.read(Main.replace('.fst','.outb')).toDataFrame()
    time      = df['Time_[s]']
    TTvx      = df['NcIMUTVxs_[m/s]']
    ws_ref    = df['RtVAvgxh_[m/s]'] # Rotor avg
    pitch     = df['BldPitch1_[deg]']
    Qaero_ref = df['RtAeroMxh_[N-m]']
    omega     = df['RotSpeed_[rpm]']*2*np.pi/60 # rad/s
    lambda_ref=omega*wse.R/ws_ref
    Qaero  = wse.Torque(ws_ref,pitch,omega)
    Thrust = wse.Thrust(ws_ref,pitch,omega)
    # ----
    print('Estimating...')
    ws_est=np.zeros(omega.shape)
    print(ws_est.shape)
    WS0=ws_ref[0]*0.9 # estimate previous time step
    for i,(Qa,p,o,ws0) in enumerate(zip(Qaero_ref,pitch,omega,ws_ref)):
        ws_hat=wse.estimate(Qa, p, o, WS0,relaxation=relaxation)
        ws_est[i]=ws_hat
        WS0=ws_hat
        if np.mod(i,1000)==0:
            print(i,len(ws_ref),'{:4.1f} {:4.1f}'.format(ws_ref[i],ws_est[i]))

    Qaero2 = wse.Torque(ws_est,pitch,omega)
    # --- Export
    M=np.column_stack((time,ws_ref,ws_est,Qaero_ref,Qaero,Qaero2,omega,pitch))
    header='time,ws_ref,ws,Qaero_ref,Qaero,Qaero2,omega,pitch'
    np.savetxt('../../_data/wsest/_{:.1f}.csv'.format(relaxation),M,delimiter=',',header=header)

if 3 in IStudy:
    base='../../_data/NREL5MW'
    #Main='../../_data/NREL5MW_SimpleSimulations/TurbWindStep_AllDOF.fst'; 
    Main='../../_data/NREL5MW/Main_Onshore.fst'; 
    wse = TabulatedWSEstimator(fst_file=Main)
    wse.load_files(base=base,suffix='')
    print(wse)
# 
plt.show()

