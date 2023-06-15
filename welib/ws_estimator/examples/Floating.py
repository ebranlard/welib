import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Local 
import weio
from welib.ws_estimator.tabulated_floating import *
from welib.tools.figure import *



if __name__ == '__main__':

    setFigureFont(15)

    pklFilename = 'C:/W0/Work/2018-NREL/DigiTwin-Stiesdal/simulations/06-aero/3-CPLambdaPitchPhi/SWT-3p6-130_CPCT_new.pkl'

    # pkl = weio.read(pklFilename)
    # print(pkl)


    R=63
    rho=1.225

    WSE = TabulatedWSEstimatorFloating(R = R, rho = rho, pickleFile = pklFilename)

    print(WSE)

    P = WSE.Power(10, 8*np.pi/30, 0, 5)
    P = WSE.Power([10,10], [8*np.pi/30, 8*np.pi/30],[0,0], [0,5])
    print(P)

    WS    = WSE.OP['WS_[m/s]'].values
    RPM   = WSE.OP['RotSpeed_[rpm]'].values
    Omega = WSE.OP['RotSpeed_[rpm]'].values*np.pi/30
    Pitch = WSE.OP['Pitch_[deg]'].values
    PhiY  = WSE.OP['PhiY_[deg]'].values
    # P0    = WSE.OP['p_[kw]'].values*1000
    P0    = WSE.OP['Paero_i_[W]'].values
    Q0    = WSE.OP['Qaero_i_[Nm]'].values
    print(WSE.OP)

    P = WSE.Power(WS, Omega, Pitch, PhiY)
    Q = WSE.Torque(WS, Omega, Pitch, PhiY)


    # ---
    WSest = np.zeros(WS.shape)



    # ---- DEBUG - TODO TODO TODO PUT SOME OF IT IN CLASS
    # WS=5
    # pitch          =1.3291208804416312
    # omega          =0.7608053160538188
    # phiy           =0.7191349687043711
    # 
    # WS1, P1 = WSE.PowerAt(omega=omega, pitch=pitch, phiy=phiy)
    # WS2, P2 = WSE.PowerAt(omega=omega, pitch=pitch, phiy=0.94)
    # 
    # 
    # print('WS',WSE.WS)
    # print('Omeg',WSE.omega)
    # print('Pitc',WSE.pitch)
    # print('PhiY',WSE.phiy)
    # 
    # ip=np.argmin(np.abs(WSE.pitch-pitch)); 
    # if WSE.pitch[ip]>pitch:
    #     ip=ip-1
    # 
    # io=np.argmin(np.abs(WSE.omega-omega)); 
    # if WSE.omega[io]>omega:
    #     io=io-1
    # 
    # ipy=np.argmin(np.abs(WSE.phiy-phiy)); 
    # if WSE.phiy[io]>phiy:
    #     ipy=ipy-1
    # 
    # iw=np.argmin(np.abs(WSE.WS-WS)); 
    # if WSE.WS[iw]>WS:
    #     iw=iw-1
    # 
    # print('>>> ', (WSE.WS[iw], WS), (WSE.omega[io], omega), (WSE.pitch[ip], pitch), (WSE.phiy[ipy],phiy))
    # 
    # 
    # fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    # fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    # 
    # ax.plot(WSE.WS[iw], WSE.P[iw, io  , ip  , ipy]  , 'k+', label='raw')
    # ax.plot(WSE.WS[iw], WSE.P[iw, io+1, ip  , ipy]  , 'k+')
    # ax.plot(WSE.WS[iw], WSE.P[iw, io  , ip+1, ipy]  , 'k+')
    # ax.plot(WSE.WS[iw], WSE.P[iw, io  , ip  , ipy+1], 'k+')
    # ax.plot(WSE.WS[iw], WSE.P[iw, io+1, ip+1, ipy+1], 'k+')
    # 
    # 
    # 
    # ax.plot(WS1, P1    , label='Interp')
    # ax.plot(WS2, P2    , label='Interp')
    # ax.plot(WSE.OP['WS_[m/s]'], WSE.OP['Paero_i_[W]'], 'k',label='OP')
    # ax.set_xlabel('')
    # ax.set_ylabel('')
    # ax.legend()
    # plt.show()



    # ---
    method= 'crossing-oper'

    # i=16
    i=50
    omega = Omega[i]*1.05
    Qa    = Q0[i]*1.1
    pitch = Pitch[i]*1.5
    phiy  = PhiY[i]*1.2
    WS0   = WS[i]*1.2
    WSest, info = WSE.estimate(Qa, omega, pitch,  phiy , WS0, relaxation=0, WSavg=None, debug=True, method=method, deltaWSMax=10)
    print(info)
    fig = WSE.debugPlot(info=info, HR=False)
    ax=fig.gca()
    ax.set_xlim([2,16])
    ax.legend(fontsize=12, ncol=2, loc='lower center')




    # 
    # 
    # 
    # 
    # method= 'crossing-oper'
    # # for i in range(3,len(WS)):
    # # for i in [7]:
    # for i in [9]:
    # # for i in [26]:
    #     omega = Omega[i]
    #     Qa    = Q0[i]
    #     pitch = Pitch[i]
    #     phiy  = PhiY[i]
    #     WS0   = WS[i]+0.9
    #     WSest[i], info = WSE.estimate(Qa, omega, pitch,  phiy , WS0, relaxation=0, WSavg=None, debug=True, method=method, deltaWSMax=1)
    #     print(WSest[i], WS[i])
    # 
    # print(info)
    # 
    # fig = WSE.debugPlot(info=info, HR=False)
    # ax=fig.gca()
    # ax.set_xlim([2,16])
    # 
    # fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    # fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    # ax.plot(WS, WS, 'k--',label='')
    # ax.plot(WS, WSest, label='')
    # ax.set_xlabel('')
    # ax.set_ylabel('')
    # ax.legend()
    # plt.show()
    # 
    # 
    # # 
    # # 
    # # # ---
    # # fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    # # fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    # # ax.plot(WS, Q0 , '+' , label='ref')
    # # ax.plot(WS, Q  , label='interp')
    # # ax.set_xlabel('')
    # # ax.set_ylabel('')
    # # ax.legend()
    # # 
    # 

    plt.show()
