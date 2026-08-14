""" 
"""
import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
# Local 
from welib.essentials import *
from welib.tools.colors import MW_Orange, fColrs
import welib.weio as weio
from welib.tools.stats import comparison_stats

import pytest

from welib.tools.compare import compare

# Main functions
from welib.yams.windturbine import YAMSSectionLoadCalculator, FASTWindTurbine
from welib.yams.windturbine import monopileSetupFromOpenFAST

scriptDir = os.path.dirname(__file__)


# 
import matplotlib
COLRS   = [fColrs(1), MW_Orange ,fColrs(2), fColrs(3)]
try:
    cmap = matplotlib.colormaps['viridis']
except:
    cmap = matplotlib.cm.get_cmap('viridis')
CMAP = [(cmap(v)[0],cmap(v)[1],cmap(v)[2],0.9) for v in np.linspace(0,1,5+1)]
COLRS=[CMAP[1], CMAP[4]]
colRef=COLRS[0]
colSim=COLRS[1]
LWRef=2.4
LWSim=1.5

def meanabs(x, **kwargs):
    return np.mean(np.abs(x), **kwargs)

def time_plot(t, ref=None, sim=None, label='', other=None, fig=None, ax=None, tRange=None, refLab='OpenFAST', otherLab='Other'):
    figNotProvided = fig is None
    if fig is None:
        fig=plt.figure()
        fig.subplots_adjust(left=0.18, right=0.94, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
        ax = fig.add_subplot(111)
    if ref is not None:
        ax.plot(t, ref,'-'  , c=colRef, lw=LWRef, label=refLab)
    if sim is not None:
        ax.plot(t, sim,'--' , c=colSim, lw=LWSim, label='YAMS')
    if other is not None:
        ax.plot(t, other,':' , label=otherLab)
    if tRange is not None:
        ax.axvline(x = tRange[0] , color = 'k', ls='--')
        ax.axvline(x = tRange[1], color = 'k', ls='--')
    if figNotProvided:
        ax.legend()
        ax.set_xlabel('Time [s]')
    ax.set_ylabel(label)
    ax.tick_params(direction='in')
    return fig, ax


def sec_plotFM(vTime, zDepth, F_sec, F_secRef=None, stat='mean', tRange=None, label='Section Force', other=None, otherLab='Other'):
    if tRange is None:
        IT = np.arange(0,F_sec.shape[1])
    else:
        IT = np.logical_and(vTime>tRange[0], vTime<tRange[1])
        if len(IT)==0:
            IT = np.arange(0, sec.shape[1])
    fstat = {'mean':np.mean, 'max':np.max, 'std':np.std, 'meanabs':meanabs}[stat]
    if F_sec is not None:
        Fx_sec0 = fstat(F_sec   [0,:,IT],axis=0) # 3, nSpan nT
        My_sec0 = fstat(F_sec   [4,:,IT],axis=0)
    if F_secRef is not None:
        Fx_sec1 = fstat(F_secRef[0,:,IT],axis=0) # 3, nSpan nT
        My_sec1 = fstat(F_secRef[4,:,IT],axis=0)
    if other is not None:
        Fx_sec2 = fstat(other[0,:,IT],axis=0) # 3, nSpan nT
        My_sec2 = fstat(other[4,:,IT],axis=0)

    fig,axes = plt.subplots(1, 2, sharey=True, figsize=(6.4,4.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.14, right=0.95, top=0.95, bottom=0.12, hspace=0.20, wspace=0.20)
    ax=axes[0]
    if F_secRef is not None:
        ax.plot(Fx_sec1/1e6, zDepth,  '-' , c=colRef, lw=LWRef, label='OpenFAST')
    if F_sec is not None:
        ax.plot(Fx_sec0/1e6, zDepth,  '--', c=colSim, lw=LWSim, label='YAMS')
    if other is not None:
        ax.plot(Fx_sec2/1e6, zDepth,  ':' ,label=otherLab)

    ax.tick_params(direction='in')

    ax.set_xlabel(label)
    ax.set_ylabel('Vertical position [m]')

    ax=axes[1]
    if F_secRef is not None:
        ax.plot(My_sec1/1e6, zDepth, '-',  c=colRef, lw=LWRef, label='OpenFAST')
    if F_sec is not None:
        ax.plot(My_sec0/1e6, zDepth, '--', c=colSim, lw=LWSim, label='YAMS')
    if other is not None:
        ax.plot(My_sec2/1e6, zDepth,  ':' ,label=otherLab)
    ax.set_xlabel('Section Moment My [MNm]')
    ax.tick_params(direction='in')
    if F_sec is not None:
        ax.legend()
    fig.suptitle(stat)
    return fig

def sec_plot(vTime, zDepth, sec, ref=None, other=None, stat='mean', tRange=None, label='Section Force', fig=None, refLab='OpenFAST', otherLab='Other'):
    if tRange is None:
        IT = np.arange(0, sec.shape[1])
    else:
        IT = np.logical_and(vTime>tRange[0], vTime<tRange[1])
        if len(IT)==0:
            IT = np.arange(0, sec.shape[1])
    fstat = {'mean':np.mean, 'max':np.max, 'std':np.std, 'meanabs':meanabs}[stat]
    sec0 = fstat(sec[:,IT], axis=1)
    if ref is not None:
        sec1 = fstat(ref[:,IT], axis=1)
    if other is not None:
        sec2 = fstat(other[:,IT], axis=1)

    if fig is None:
        fig, ax = plt.subplots(1, 1, sharey=True, figsize=(6.4,4.8)) # (6.4,4.8)
        fig.subplots_adjust(left=0.14, right=0.95, top=0.95, bottom=0.12, hspace=0.20, wspace=0.20)
    else:
        ax=fig.axes[0]

    if ref is not None:
        ax.plot(sec1, zDepth,  '-' , c=colRef, lw=LWRef, label=refLab)
    ax.plot(sec0,     zDepth,  '--', c=colSim, lw=LWSim, label='YAMS')
    if other is not None:
        ax.plot(sec2, zDepth,  '-' , c=colRef, lw=LWRef, label=refLab)
    ax.tick_params(direction='in')
    ax.legend()

    ax.set_xlabel(label)
    ax.set_ylabel('Vertical position [m]')
    fig.suptitle(stat)
    return fig




def section_loads(fstFile, compFile=None, hydroShapeFile=None, subShapes=None, tMin=1, tMax=100, plot=False, dtSamp=0.1):
    """"
    Wrapper, that should rely heavily on windturbine.py (WindTurbineStructure)
    Should compute section loads along tower and monopile
    Should work for a monopile or a floater
    """
    fstFile = os.path.join(scriptDir, fstFile)
    if compFile is not None:
        compFile = os.path.join(scriptDir, compFile)

    Isec = [0,4,7] # Section indices

    # --- Derived parameters
    fstOut = fstFile.replace('.fst','.outb')
    outFile = fstOut.replace('.out','').replace('.outb','')
    outFileSL = outFile + '_SL_YAMS.outb'


    # --- Reference Time series for verification
    df = weio.read(fstOut).toDataFrame()
    df=df[df['Time_[s]']<=tMax]
    df=df[df['Time_[s]']>=tMin]
    df.reset_index(inplace=True)
    # Re sampling at dtSamp to have less datapoints
    t_new = np.arange(df['Time_[s]'].min(), df['Time_[s]'].max(), dtSamp )
    f = interp1d(df['Time_[s]'], df.drop(columns=['Time_[s]']), axis=0, kind='linear', fill_value='extrapolate')
    dfRef = pd.DataFrame(f(t_new), columns=df.columns.drop('Time_[s]'))
    dfRef.insert(0, 'Time_[s]', t_new)
#     dfRef = df

    # --- Start
    WT = FASTWindTurbine(fstFile, algo='OpenFAST', HD_compFile=compFile, SD_FEM_method='cbeam', subShapes=subShapes, verbose=True).WT
    YSL = YAMSSectionLoadCalculator(WT=WT)

    if WT.pSS is not None:
        NOTE('Setting Components', compFile)
        WT.SS_setComponents(compFile)
        NOTE('Setting Compute Eta')
        WT.SS_computeEta(dfRef['Time_[s]'])
        if hydroShapeFile is not None:
            WT.HD_setShapeFunction(hydroShapeFile)


        # --- LEGACY
#         shapes_sub = [0, 4]
#         pSTm, pSSm, pHDm, Sysm, WTm, refm = monopileSetupFromOpenFAST( fstFile, shapes_sub=shapes_sub, TMIN=tMin, TMAX=tMax, compFile=compFile, tuneM=False)
#         WT.monopileSetup = {'WTm':WTm, 'pST': pSTm, 'ref': refm, 'sys': Sysm}
#         #b1=compare(WT.pSS, pSSm, verbose=False)
#         #b2=compare(WT.pHD, pHDm, verbose=False)
#         #b2=compare(WT.fnd, WTm.fnd, verbose=False)
#         # compare(self.WT, WTm, verbose=False, n1='WT', n2='WTm')
        zBeamRef, F_secRef, r_secRef =  WT.fnd.SD.beamSecOutputs(dfRef, verbose=False)
    else:
        F_secRef = None

    dfOut, sections = YSL.fromDF(dfRef, useTopLoadsFromDF=True, useInterfaceLoadsFromDF=True) # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< TODO

    dfOut.export(outFileSL)

    # --------------------------------------------------------------------------------}
    # --- PLOT Tower sections
    # --------------------------------------------------------------------------------{
    if plot:
        if 'TwHt1MLyt_[kN-m]' not in dfRef:
            WARN('Not plotting tower section loads, no data ref')
        else:

            fig,ax = plt.subplots(1, 1, sharey=False, figsize=(12.4,12.8)) # (6.4,4.8)
            fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
            for iiED,iED in enumerate(Isec):
                sT = 'TwHt{}MLyt_[kN-m]'.format(iED+1)

                t1 = dfRef['Time_[s]'].values
                y1 = dfRef[sT].values
                t2 = dfOut['Time_[s]'].values
                y2 = dfOut[sT].values

                ax.plot(t1, y1, 'k-'                        , label='OpenFAST' if iiED==0 else None)
                ax.plot(t2, y2,  '--' , color=fColrs(iiED)  , label='Sim Ht{}'.format(iED+1))

            ax.set_xlabel('')
            ax.set_ylabel('')
            ax.legend()

            outFileFigTwr = outFile + '_SL_TWR.png'
            fig.savefig(outFileFigTwr)


    # --------------------------------------------------------------------------------}
    # --- Plots Monpile
    # --------------------------------------------------------------------------------{
    zDepth   = None
    if WT.pHD is not None:
        zDepth   = WT.pHD['zDepth']
        if plot:
            F_sec    = sections['monopile']['F_sec']
            vTime    = dfRef['Time_[s]']
            tRange = [tMin, tMax]

            # --- Section loads Plot
            IZ = [int(3*len(zDepth)/6)-5, int(2*len(zDepth)/6), int(1*len(zDepth)/6), 0]
            fig = sec_plotFM(vTime, zDepth, F_sec, F_secRef, stat='meanabs', tRange=tRange)
            fig.savefig(outFile + '_SL_MNP_Depth.png')

            fig, axes = plt.subplots(len(IZ), 1, sharey=True, sharex=True, figsize=(8.4,7.5))
            fig.subplots_adjust(left=0.12, right=0.96, top=0.95, bottom=0.05, hspace=0.07, wspace=0.20)
            for ii, iz in enumerate(IZ):
                time_plot (vTime, F_secRef[0, iz, :]/1e6,  F_sec[0, iz, :]/1e6, f'z={zDepth[iz]:.0f}m', tRange=tRange, ax=axes[ii], fig=fig)
                stats, sStats =  comparison_stats(vTime, F_secRef[0,iz,:]/1e6, vTime, F_sec[0,iz,:]/1e6, stats='sigRatio,eps,R2', method='meanabs')
                print(f'z {zDepth[iz]:5.0f}: ', stats)
            axes[-1].set_xlabel('Time [s]')
            fig.suptitle('FxSec [MN]')
            fig.savefig(outFile + '_SL_MNP.png')


            fig, ax = time_plot (vTime, dfRef['Wave1Elev_[m]'], dfOut['Wave1Elev_[m]'], 'Wave elevation [m]', tRange=tRange)
            ax.plot(WT.pSS['eta_time'], WT.pSS['eta'], ':', label='From SS')
            fig.savefig(outFile + '_SL_MNP_Eta.png')

            fig, ax = time_plot (vTime, dfRef['HydroFxi_[N]'], dfOut['HydroFxi_[N]'], 'Hydro Fx [N]', tRange=tRange)
            fig.savefig(outFile + '_SL_MNP_HydroFxi.png')

            fig, ax = time_plot (vTime, F_secRef[4,0,:]/1e6, F_sec[4,0,:]/1e6, label='Sea bed moment [MNm]', tRange=tRange)
            fig.savefig(outFile + '_SL_MNP_M_SeaBed.png')

            stats, sStats =  comparison_stats(vTime, dfRef['Wave1Elev_[m]'], vTime, dfOut['Wave1Elev_[m]'], stats='sigRatio,eps,R2', method='meanabs')
            print('Eta     :', stats)
            stats, sStats =  comparison_stats(vTime, F_secRef[4,0,:], vTime, F_sec[4,0,:], stats='sigRatio,eps,R2', method='meanabs')
            print('M_bot   :', stats)

#     if plot:
#         plt.show()
    return {'dfRef': dfRef, 'dfOut': dfOut, 'F_secRef':F_secRef, 'sections': sections, 'zDepth':zDepth}



def test_floating_tower_TS(plot=False, test=True):
    return

    # Test for tower loads only since this test case has no "monopile", just a floater
    fstFile='06_Jonswap_TS/Main.fst';
    if test:
        out = section_loads(fstFile, plot=plot, tMin=10, tMax=40, dtSamp=0.2)
    else:
        out = section_loads(fstFile, plot=plot, tMin=0, tMax=700, dtSamp=0.05)
    dfRef = out['dfRef']
    dfOut  = out['dfOut']


    # --- Tests
    IsecTwr = [0,4,7] # Section indices
    for iiED,iED in enumerate(IsecTwr):
        t1 = dfRef['Time_[s]'].values
        t2 = dfOut ['Time_[s]'].values
        y1 = np.asarray(dfRef['TwHt{}MLyt_[kN-m]'.format(iED+1)].values, dtype=float)/1000
        y2 = np.asarray(dfOut ['TwHt{}MLyt_[kN-m]'.format(iED+1)  ].values, dtype=float)/1000
        stats, sStats =  comparison_stats(t1, y1, t2, y2, stats='sigRatio,eps,R2', method='mean')
        #print(stats)
        np.testing.assert_array_less(1-np.abs(stats['sigRatio']), 0.014)
        np.testing.assert_array_less(np.abs(stats['eps']), 1.5)
        np.testing.assert_array_less(1-np.abs(stats['R2']), 0.001)

        signal_range = np.ptp(y2)  # Peak-to-peak range
        atol = 0.007 * signal_range  # % of total signal range
        np.testing.assert_allclose(y1, y2, rtol=1e-6, atol=atol)


def test_monopile_only_MT100(plot=False, test=True):
#     fstFile = '05_RegWave_IEA/OF.fst';
#     fstFile = '05_RegWave_IEA/OF_F2_NoRNA.fst'; 
#     fstFile = '05_RegWave_IEA/OF_F3_NoRNA.fst'; 
    # fstFile = '05_RegWave_IEA/OF_F3_RNA.fst'; 
    # fstFile = '05_RegWave_IEA/OF_F3T0_NoRNA.fst';  compFile=None
#     fstFile  ='06_Jonswap_MT100/OF_Long_Hs=8.1_Tp=12.7_h=50.fst'
#     compFile ='06_Jonswap_MT100/UserDefJonswap_Hs=8.1_Tp=12.7_h=50.csv'
    fstFile  = os.path.join(scriptDir, '../../../data/Monopile/Main_MT100_JONSWAP_UserDef.fst')
    compFile = os.path.join(scriptDir, '../../../data/Monopile/Waves/UserDefJonswap_Hs=8.1_Tp=12.7_h=50.csv')
    if test:
        out = section_loads(fstFile, compFile=compFile, tMin=8, tMax=12, plot=plot, subShapes=[0,4])
    else:
        out = section_loads(fstFile, compFile=compFile, tMin=35, tMax=100, plot=plot, subShapes=[0,4])

    dfRef = out['dfRef']
    dfOut = out['dfOut']
    zDepth   = out['zDepth']
    F_secRef = out['F_secRef']
    F_sec    = out['sections']['monopile']['F_sec']

    plt.show()
    # --- Section loads as function of depth
    Fx_sec0 = meanabs(F_sec   [0,:,:],axis=0) # 3, nSpan, nT
    My_sec0 = meanabs(F_sec   [4,:,:],axis=0)
    Fx_sec1 = meanabs(F_secRef[0,:,:],axis=0) # 3, nSpan, nT
    My_sec1 = meanabs(F_secRef[4,:,:],axis=0)

    np.testing.assert_allclose(Fx_sec0/1e6, Fx_sec1/1e6, atol=0.140)
    np.testing.assert_allclose(My_sec0/1e8, My_sec1/1e8, atol=0.040)

    # --- Section loads as function of time and depth
    vTime = dfRef['Time_[s]']
    IZ = [int(3*len(zDepth)/6)-5, int(2*len(zDepth)/6), int(1*len(zDepth)/6), 0]
    for ii, iz in enumerate(IZ):
        stats, sStats =  comparison_stats(vTime, F_secRef[0,iz,:]/1e6, vTime, F_sec[0,iz,:]/1e6, stats='sigRatio,eps,R2', method='meanabs')
        #print(f'z {zDepth[iz]:5.0f}: ', stats)
        np.testing.assert_array_less(np.abs(stats['eps']), 36)
        np.testing.assert_array_less(1-np.abs(stats['R2']), 0.91)
        np.testing.assert_allclose(F_secRef[0,iz,:]/1e6, F_sec[0,iz,:]/1e6, atol=0.48 )

    # Wave elevation and hydro loads are quite accurate
    np.testing.assert_allclose(dfRef['Wave1Elev_[m]']   , dfOut['Wave1Elev_[m]']   , atol=0.001)
    np.testing.assert_allclose(dfRef['HydroFxi_[N]']/1e7, dfOut['HydroFxi_[N]']/1e7, atol=0.020)

    # Sea bed moment
    stats, sStats =  comparison_stats(vTime, F_secRef[4,0,:], vTime, F_sec[4,0,:], stats='sigRatio,eps,R2', method='meanabs')
    np.testing.assert_allclose(F_secRef[4,0,:]/1e9, F_sec[4,0,:]/1e9, atol=0.015)


def test_monopile_tower_IEA(plot=False, test=True):
    return
    # fstFile='06_Jonswap_IEA/OF_F3T1S1_H1A1_Hs=8.1_Tp=12.7.fst'
    #fstFile='06_Jonswap_IEA/OF_F2T1S0_H1A1_Hs=8.1_Tp=12.7.fst'; 
    fstFile='06_Jonswap_IEA/OF_F2T1S1_H1A1_Hs=8.1_Tp=12.7.fst'; 
    compFile =f'Waves/UserDefJonswap_Hs=8.1_Tp=12.7_h=34.csv'
    tMin=20
    if test:
        tMax=100
    else:
        tMax=150

    out = section_loads(fstFile, plot=plot, compFile=compFile, tMin=tMin, tMax=tMax, dtSamp=0.2)


    dfRef = out['dfRef']
    dfOut = out['dfOut']
    zDepth    = out['zDepth']
    F_secRef  = out['F_secRef']
    F_sec     = out['sections']['monopile']['F_sec']
    
    IsecTwr = [0,4,7] # Section indices
    for iiED,iED in enumerate(IsecTwr):
        t1 = dfRef['Time_[s]'].values
        t2 = dfOut ['Time_[s]'].values
        y1 = np.asarray(dfRef['TwHt{}MLyt_[kN-m]'.format(iED+1)].values, dtype=float)/1000
        y2 = np.asarray(dfOut ['TwHt{}MLyt_[kN-m]'.format(iED+1)  ].values, dtype=float)/1000
        stats, sStats =  comparison_stats(t1, y1, t2, y2, stats='sigRatio,eps,R2', method='mean')
        print(stats)

        # TODO AI: Accuracy is low, figure out why
        if test:
            np.testing.assert_array_less(1-np.abs(stats['sigRatio']), 0.08)
            np.testing.assert_array_less(np.abs(stats['eps']), 18)
            np.testing.assert_array_less(1-np.abs(stats['R2']), 0.08)

            signal_range = np.ptp(y2)  # Peak-to-peak range
            atol = 0.30 * signal_range  # % of total signal range
            np.testing.assert_allclose(y1, y2, rtol=1e-3, atol=atol)

    # --- Monopile loads from the same consolidated wrapper

    # --- Section loads as function of depth
    Fx_sec0 = meanabs(F_sec   [0,:,:],axis=0) # 3, nSpan, nT
    My_sec0 = meanabs(F_sec   [4,:,:],axis=0)
    Fx_sec1 = meanabs(F_secRef[0,:,:],axis=0) # 3, nSpan, nT
    My_sec1 = meanabs(F_secRef[4,:,:],axis=0)
    if test:
        np.testing.assert_allclose(Fx_sec0/1e6, Fx_sec1/1e6, atol=2.500)
        np.testing.assert_allclose(My_sec0/1e8, My_sec1/1e8, atol=2.500)

    # --- Section loads as function of time and depth
    vTime = dfRef['Time_[s]']
    IZ = [int(3*len(zDepth)/6)-5, int(2*len(zDepth)/6), int(1*len(zDepth)/6), 0]
    for ii, iz in enumerate(IZ):
        stats, sStats =  comparison_stats(vTime, F_secRef[0,iz,:]/1e6, vTime, F_sec[0,iz,:]/1e6, stats='sigRatio,eps,R2', method='meanabs')
        #print(f'z {zDepth[iz]:5.0f}: ', stats)
        #np.testing.assert_array_less(np.abs(stats['eps']), 33)
        #np.testing.assert_array_less(1-np.abs(stats['R2']), 0.09)
        if test:
            np.testing.assert_allclose(F_secRef[0,iz,:]/1e6, F_sec[0,iz,:]/1e6, atol=2.60 )

    if test:
        # Wave elevation and hydro loads are quite accurate
        np.testing.assert_allclose(dfRef['Wave1Elev_[m]']   , dfOut['Wave1Elev_[m]']   , atol=0.0001)
        np.testing.assert_allclose(dfRef['HydroFxi_[N]']/1e7, dfOut['HydroFxi_[N]']/1e7, atol=0.018)

        # Sea bed moment
        np.testing.assert_allclose(F_secRef[4,0,:]/1e9, F_sec[4,0,:]/1e9, atol=0.050)



if __name__ == '__main__':
    PLOT = True # KEEP ME FOR EASY DEBUG
    TEST=False
    TEST=True
#     test_monopile_tower_IEA(plot=PLOT, test=TEST)  # Was used to develop monopile only 
    test_monopile_only_MT100(plot=PLOT, test=TEST) # Need updating of main code in windturbine and merging of "debug_*" functions
#     test_floating_tower_TS(plot=PLOT, test=TEST)   # Was used to develop tower only
    plt.show()
