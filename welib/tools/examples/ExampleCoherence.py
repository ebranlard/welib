""" 
Example on how to compute coherence and cross correlation for a turbulence box

NOTE: this scripts needs a lot of cleanup!

"""


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.signal as sig
# Local 
import welib.weio as weio
from welib.tools.spectral import fft_wrap
from welib.tools.colors import *

def turbBoxCoherence(btsFile):
    ts = weio.read(btsFile)
    print(ts)

    iy0,iz0 = ts._iMid()
    y       = ts['y']
    z       = ts['z']
    yMid    = ts['y'][iy0]
    zMid    = ts['z'][iz0]
    t       = ts['t']
    dt      = t[1]-t[0]
    fs      = 1/dt

    u = ts['u'][0,:,iy0,iz0]
    v = ts['u'][1,:,iy0,iz0]
    w = ts['u'][2,:,iy0,iz0]
    U0=np.mean(u)
    u= u-np.mean(u)
    v= v-np.mean(v)
    w= w-np.mean(w)

    # TODO put this as functions in TurbSimFile

    # --- Cross correlation along y, mid box
    rho_uu_y=np.zeros(len(y))
    rho_vv_y=np.zeros(len(y))
    rho_ww_y=np.zeros(len(y))
    for iy,_ in enumerate(y):
        ud = ts['u'][0,:,iy,iz0]
        vd = ts['u'][1,:,iy,iz0]
        wd = ts['u'][2,:,iy,iz0]
        ud-=np.mean(ud)
        vd-=np.mean(vd)
        wd-=np.mean(wd)
        rho_uu_y[iy] = np.mean(u*ud)/(np.std(u)*np.std(ud))
        rho_vv_y[iy] = np.mean(v*vd)/(np.std(v)*np.std(vd))
        rho_ww_y[iy] = np.mean(w*wd)/(np.std(w)*np.std(wd))

    # --- Cross correlation along z, mid box
    rho_uu_z = np.zeros(len(z))
    rho_vv_z = np.zeros(len(z))
    rho_ww_z = np.zeros(len(z))
    for iz,_ in enumerate(z):
        ud = ts['u'][0,:,iy0,iz]
        vd = ts['u'][1,:,iy0,iz]
        wd = ts['u'][2,:,iy0,iz]
        ud-=np.mean(ud)
        vd-=np.mean(vd)
        wd-=np.mean(wd)
        rho_uu_z[iz] = np.mean(u*ud)/(np.std(u)*np.std(ud))
        rho_vv_z[iz] = np.mean(v*vd)/(np.std(v)*np.std(vd))
        rho_ww_z[iz] = np.mean(w*wd)/(np.std(w)*np.std(wd))

    # --- Plot cross-correlation
    fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    ax.plot(y , rho_uu_y   , label=r'$\rho_{uu}$')
    ax.plot(y , rho_vv_y   , label=r'$\rho_{vv}$')
    ax.plot(y , rho_ww_y   , label=r'$\rho_{ww}$')
    ax.set_xlabel('y')
    ax.set_ylabel('')
    ax.legend()

    fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    ax.plot(z , rho_uu_z   , label=r'$\rho_{uu}$')
    ax.plot(z , rho_vv_z   , label=r'$\rho_{vv}$')
    ax.plot(z , rho_ww_z   , label=r'$\rho_{ww}$')
    ax.set_xlabel('z')
    ax.set_ylabel('')
    ax.legend()
    # 
    fc, chi_uu = sig.csd(u, u, fs=fs, scaling='density') #nperseg=4096, noverlap=2048, detrend='constant')
    fc, chi_vv = sig.csd(v, v, fs=fs, scaling='density') #nperseg=4096, noverlap=2048, detrend='constant')
    fc, chi_ww = sig.csd(w, w, fs=fs, scaling='density') #nperseg=4096, noverlap=2048, detrend='constant')

    # kc=2*np.pi*fc/U0
    # chi_uu= 2*np.pi/U0*chi_uu
    # chi_vv= 2*np.pi/U0*chi_vv
    # chi_ww= 2*np.pi/U0*chi_ww


    # print(chi_uu)

    # fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    # fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    # ax.plot(t, uMid    , label='')
    # ax.set_xlabel('')
    # ax.set_ylabel('')
    # ax.legend()

    # f1, S1, Info  = fft_wrap(t,u,output_type='psd',averaging='Welch', averaging_window='hann',) # <<< TODO Play with nExp
    # f2, S2, Info  = fft_wrap(t,u,output_type='PSD',averaging='none')
    # print(chi_uu/S1)
    # print(Info.__dict__)

    fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.30, wspace=0.20)
    # ax.plot(f1,f1*S1,'r',label ='Welch')
    # ax[1].plot(f2,S2,'k--',label ='None')
    ax.plot(fc,fc*chi_uu,'-',label =r'$\chi_{uu}$')
    ax.plot(fc,fc*chi_vv,'-',label =r'$\chi_{vv}$')
    ax.plot(fc,fc*chi_ww,'-',label =r'$\chi_{ww}$')
    # ax.set_xlabel('Frequency [Hz]')
    # ax.set_ylabel('FFT Amplitude')
    # ax.set_yscale("log", nonpositive='clip')
    ax.set_xscale("log")
    ax.legend()
    # ax.set_title('Signal - FFT')
    # 


    # ---
    diy=1
    iy = iy0+diy
    dy=y[iy]-y[iy0]
    ud = ts['u'][0,:,iy,iz0]
    vd = ts['u'][1,:,iy,iz0]
    wd = ts['u'][2,:,iy,iz0]
    ud-=np.mean(ud)
    vd-=np.mean(vd)
    wd-=np.mean(wd)
    fc, coh_uu_y1 = sig.coherence(u,ud, fs=fs)
    fc, coh_vv_y1 = sig.coherence(v,vd, fs=fs)
    fc, coh_ww_y1 = sig.coherence(w,wd, fs=fs)
    iy = iy+diy
    ud = ts['u'][0,:,iy,iz0]
    vd = ts['u'][1,:,iy,iz0]
    wd = ts['u'][2,:,iy,iz0]
    ud-=np.mean(ud)
    vd-=np.mean(vd)
    wd-=np.mean(wd)
    fc, coh_uu_y2 = sig.coherence(u,ud, fs=fs)
    fc, coh_vv_y2 = sig.coherence(v,vd, fs=fs)
    fc, coh_ww_y2 = sig.coherence(w,wd, fs=fs)

    fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.30, wspace=0.20)
    ax.plot(fc,coh_uu_y1,'-' ,color=fColrs(1),label =r'$coh_{uu}$, '+r'$\Delta y ={:.0f}$'.format(dy))
    ax.plot(fc,coh_vv_y1,'-' ,color=fColrs(2),label =r'$coh_{vv}$, '+r'$\Delta y ={:.0f}$'.format(dy))
    ax.plot(fc,coh_ww_y1,'-' ,color=fColrs(3),label =r'$coh_{ww}$, '+r'$\Delta y ={:.0f}$'.format(dy))
    ax.plot(fc,coh_uu_y2,'--',color=fColrs(1),label =r'$coh_{uu}$, '+r'$\Delta y ={:.0f}$'.format(2*dy))
    ax.plot(fc,coh_vv_y2,'--',color=fColrs(2),label =r'$coh_{vv}$, '+r'$\Delta y ={:.0f}$'.format(2*dy))
    ax.plot(fc,coh_ww_y2,'--',color=fColrs(3),label =r'$coh_{ww}$, '+r'$\Delta y ={:.0f}$'.format(2*dy))
    ax.set_xlabel('Frequency [Hz]')
    # ax.set_ylabel('FFT Amplitude')
    # ax.set_yscale("log", nonpositive='clip')
    ax.set_xscale("log")
    ax.legend()




if __name__ == '__main__':
    # btsFile = 'C:/Work/IEA29/DanAero/Phase_IV.3.2_Wakes/Inflow/B1_5.bts'
    btsFile = './wind_ws07.bts'
    turbBoxCoherence(btsFile)
    plt.show()
    pass
