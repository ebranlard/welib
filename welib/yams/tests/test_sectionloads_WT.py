""" 
"""
import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
# Local 
from welib.essentials import *
from welib.tools.colors import fColrs
import welib.weio as weio
from welib.tools.stats import comparison_stats
from welib.tools.testing import compareSigWithStats

import pytest

# Main functions
from welib.yams.windturbine import FASTWindTurbine
from welib.yams.section_loads_WT import YAMSSectionLoadCalculator
from welib.yams.section_loads_WT import sec_plotFM, time_plot

scriptDir = os.path.dirname(__file__)

def meanabs(x, **kwargs):
    return np.mean(np.abs(x), **kwargs)

def section_loads(fstFile, compFile=None, hydroShapeFile=None, subShapes=None, tMin=1, tMax=100, plot=False, dtSamp=0.1, fixedShaft=False):
    """"
    Wrapper, that relies heavily on windturbine.py (WindTurbineStructure) 
    to compute section loads along tower and monopile.
    Should work for a monopile or a floater
    """
    tRange = [tMin, tMax]
    # --- Derived parameters
    fstOut    = fstFile.replace('.fst','.outb')
    outBase   = fstOut.replace('.out','').replace('.outb','')
    outFigDir = os.path.join(os.path.dirname(fstOut),  'figs/')
    outFile   = os.path.join(outFigDir, os.path.basename(outBase))
    os.makedirs(outFigDir, exist_ok=True)

    # --- Start
    WT = FASTWindTurbine(fstFile, algo='OpenFAST', HD_compFile=compFile, SD_FEM_method='cbeam', subShapes=subShapes, verbose=True, fixedShaft=fixedShaft).WT
    YSL = YAMSSectionLoadCalculator(WT=WT)

    if WT.pSS is not None:
        NOTE('Setting Components', compFile)
        WT.SS_setComponents(compFile)
#         NOTE('Setting Compute Eta')
#         WT.SS_computeEta(dfRef['Time_[s]'])
        if hydroShapeFile is not None:
            WT.HD_setShapeFunction(hydroShapeFile)
    else:
        F_secRef = None

    df = weio.read(fstOut).toDataFrame() # Reference Time series for verification
    dfOut, sections = YSL.fromDF(df, useTopLoadsFromDF=True, useInterfaceLoadsFromDF=False, dt_resample=dtSamp, tRange=tRange)

    _, F_secRef, _ =  WT.fnd.SD.beamSecOutputs(YSL.dfRef, verbose=False)
    dfOut.export(outBase + '.YAMS.outb')


    # --- Plot Tower sections
    if plot:
        if 'TwHt1HtMLyt_[kN-m]'  in YSL.dfRef:
            YSL.plot_tower_section_loads(IsecTwr=[0,4,8], component='MLyt', figFilename=outFile + '_SL_TWR_SL.png')
            YSL.plot_tower_accelerations(IsecTwr=[0,4,8], component='x', figFilename=outFile+'_SL_TWR_Acc.png')

    # --- Plot Monpile
    if WT.pHD is not None:
        if plot:
            zDepth = YSL.sec['monopile']['z']
            IZ = [int(3*len(zDepth)/6)-5, int(2*len(zDepth)/6), int(1*len(zDepth)/6), 0]
            YSL.plot_monopile_section_loads_stats(IZ=IZ, component=0, figFilename=outFile+'_SL_MNP_SL_STATS.png', tRange=tRange, stat='meanabs')
            YSL.plot_monopile_section_loads      (IZ=IZ, component=0, figFilename=outFile+'_SL_MNP_SL.png', tRange=tRange)

            fig = YSL.plot_comp(sig='Wave1Elev_[m]', ylabel='Wave elevation [m]', figFilename=outFile+'_SL_MNP_Eta.png', tRange=tRange)
            fig.axes[0].plot(WT.pSS['eta_time'], WT.pSS['eta'], ':', label='From SS')
            fig.savefig(outFile + '_SL_MNP_Eta.png')

            fig = YSL.plot_comp(sig='HydroFxi_[N]', ylabel='Hydro Fx [N]', figFilename=outFile+'_SL_MNP_Eta.png', tRange=tRange)
            fig = YSL.plot_comp(sig='M1N1MKye_[N*m]', ylabel='Sea bed moment [MNm]', figFilename=outFile+'_SL_MNP_M_SeaBed.png', tRange=tRange, scale=1e6)
            #fig, ax = time_plot (vTime, F_secRef[4,0,:]/1e6, F_sec[4,0,:]/1e6, label='Sea bed moment [MNm]', tRange=tRange)

    return YSL

def test_monopile_only_MT100(plot=False, test=True):
    fstFile  = os.path.join(scriptDir, '../../../data/Monopile/Main_MT100_JONSWAP_UserDef.fst')
    compFile = os.path.join(scriptDir, '../../../data/Monopile/Waves/UserDefJonswap_Hs=8.1_Tp=12.7_h=50.csv')
    if test:
        YSL = section_loads(fstFile, compFile=compFile, tMin=8, tMax=12, plot=plot, subShapes=[0,4], fixedShaft=True)
    else:
        YSL = section_loads(fstFile, compFile=compFile, tMin=35, tMax=100, plot=plot, subShapes=[0,4], fixedShaft=True)

    dfRef    = YSL.dfRef
    dfOut    = YSL.dfOut
    zDepth   = YSL.sec['monopile']['z']
    F_secRef = YSL.sec['monopile']['F_secRef']
    F_sec    = YSL.sec['monopile']['F_sec']

    # --- Section loads as function of depth
    Fx_sec0 = meanabs(F_sec   [0,:,:],axis=0) # 3, nSpan, nT
    My_sec0 = meanabs(F_sec   [4,:,:],axis=0)
    Fx_sec1 = meanabs(F_secRef[0,:,:],axis=0) # 3, nSpan, nT
    My_sec1 = meanabs(F_secRef[4,:,:],axis=0)
    compareSigWithStats(Fx_sec0, Fx_sec1, t=zDepth, sig='Fx_sec(z)', epsTolP=0.060, atol=0.14, printStats=True, test=test, factor=1e6)
    compareSigWithStats(My_sec0, My_sec1, t=zDepth, sig='My_sec(z)', epsTolP=0.061, atol=0.04, printStats=True, test=test, factor=1e8)

    # --- Section loads as function of time and depth
    vTime = dfRef['Time_[s]']
    IZ = [int(3*len(zDepth)/6)-5, int(2*len(zDepth)/6), int(1*len(zDepth)/6), 0]
    for ii, iz in enumerate(IZ):
        compareSigWithStats(F_sec[0,iz,:], F_secRef[0,iz,:], t=vTime, sig=f'Fx_sec{iz}', epsTolP=35.5, atol=0.48, printStats=True, test=test, factor=1e6)
        compareSigWithStats(F_sec[4,iz,:], F_secRef[4,iz,:], t=vTime, sig=f'My_sec{iz}', epsTolP=35.0, atol=13.6, printStats=True, test=test, factor=1e6)
    compareSigWithStats(F_sec[4,iz,:], F_secRef[4,iz,:], t=vTime, sig=f'My_sec{iz}', epsTolP=0.06, atol=0.015, printStats=True, test=test, factor=1e9)

    # Wave elevation and hydro loads are quite accurate
    compareSigWithStats(dfOut, dfRef, sig='Wave1Elev_[m]', epsTolP=0.001, atol=0.0010, printStats=True, test=True)
    compareSigWithStats(dfOut, dfRef, sig='HydroFxi_[N]',  epsTolP=0.03,  atol=0.13, printStats=True, test=True, factor=1e6)

def test_monopile_tower_IEA(plot=False, test=True):
    fstFile='_06_Jonswap_IEA/OF_F3T1S1_H1A1_Hs=8.1_Tp=12.7.fst'; 
    fstFile  = os.path.join(scriptDir, '../../../data/IEA-22-280-RWT/Jonswap/OF_F3T1S1_H1A1_Hs=8.1_Tp=12.7.fst')
    compFile = os.path.join(scriptDir, '../../../data/IEA-22-280-RWT/Jonswap/UserDefJonswap_Hs=8.1_Tp=12.7_h=34.csv')
    tMin=20
    if test:
        tMin,tMax=5, 10
    else:
        tMax=150
        tMin,tMax=100, 130

    YSL = section_loads(fstFile, plot=plot, compFile=compFile, tMin=tMin, tMax=tMax, dtSamp=0.05, subShapes=[0,4])
    dfRef    = YSL.dfRef
    dfOut    = YSL.dfOut
    zDepth   = YSL.sec['monopile']['z']
    F_secRef = YSL.sec['monopile']['F_secRef']
    F_sec    = YSL.sec['monopile']['F_sec']
        
    # --- Tower Disp
    compareSigWithStats(dfOut, dfRef, sig='TwHt9TDxt_[m]', epsTolP=1e-5, atol=1e-6, printStats=True, test=test)
    compareSigWithStats(dfOut, dfRef, sig='TwHt5TDxt_[m]', epsTolP=1e-5, atol=1e-6, printStats=True, test=test)
    compareSigWithStats(dfOut, dfRef, sig='TwHt1TDxt_[m]', epsTolP=1e-5, atol=1e-6, printStats=True, test=test)
    # TODO TDz not ready
    compareSigWithStats(dfOut, dfRef, sig='TwHt9TPxi_[m]', epsTolP=1e-3, atol=1e-4, printStats=True, test=test)
    compareSigWithStats(dfOut, dfRef, sig='TwHt5TPxi_[m]', epsTolP=1e-3, atol=1e-4, printStats=True, test=test)
    compareSigWithStats(dfOut, dfRef, sig='TwHt1TPxi_[m]', epsTolP=1e-3, atol=1e-4, printStats=True, test=test)


    # --- Tower and RNA accelerations
    compareSigWithStats(dfOut, dfRef, sig='NcIMUTAxs_[m/s^2]', epsTolP=0.012, atol=0.0012, printStats=True, test=test)
    compareSigWithStats(dfOut, dfRef, sig='TwHt9ALxt_[m/s^2]', epsTolP=0.003, atol=0.0005, printStats=True, test=test)
    compareSigWithStats(dfOut, dfRef, sig='TwHt5ALxt_[m/s^2]', epsTolP=0.003, atol=0.0005, printStats=True, test=test)
    compareSigWithStats(dfOut, dfRef, sig='TwHt1ALxt_[m/s^2]', epsTolP=0.003, atol=0.0005, printStats=True, test=test)

    # --- Tower section loads
    IsecTwr = [0,4,8] # Section indices
    for iiED,iED in enumerate(IsecTwr):
        compareSigWithStats(dfOut, dfRef, sig = f'TwHt{iED+1}MLxt_[kN-m]', epsTolP=0.005, atol=0.085, printStats=True, test=test, factor=1e3)
    for iiED,iED in enumerate(IsecTwr):
        compareSigWithStats(dfOut, dfRef, sig = f'TwHt{iED+1}MLyt_[kN-m]', epsTolP=0.05, atol=0.310, printStats=True, test=test, factor=1e4)
    for iiED,iED in enumerate(IsecTwr):
        compareSigWithStats(dfOut, dfRef, sig = f'TwHt{iED+1}MLzt_[kN-m]', epsTolP=0.015, atol=0.050, printStats=True, test=test, factor=1e3)
    for iiED,iED in enumerate(IsecTwr):
        compareSigWithStats(dfOut, dfRef, sig = f'TwHt{iED+1}FLxt_[kN]', epsTolP=0.02, atol=0.310, printStats=True, test=test, factor=1e2)
    for iiED,iED in enumerate(IsecTwr):
        compareSigWithStats(dfOut, dfRef, sig = f'TwHt{iED+1}FLyt_[kN]', epsTolP=0.03, atol=0.110, printStats=True, test=test, factor=1e1)
    for iiED,iED in enumerate(IsecTwr):
        compareSigWithStats(dfOut, dfRef, sig = f'TwHt{iED+1}FLzt_[kN]', epsTolP=0.05, atol=0.001, printStats=True, test=test, factor=1e4)

    # --- Monopile Section loads as function of depth
    Fx_sec0 = meanabs(F_sec   [0,:,:],axis=0) # 3, nSpan, nT
    My_sec0 = meanabs(F_sec   [4,:,:],axis=0)
    Fx_sec1 = meanabs(F_secRef[0,:,:],axis=0) # 3, nSpan, nT
    My_sec1 = meanabs(F_secRef[4,:,:],axis=0)
    compareSigWithStats(Fx_sec0, Fx_sec1, t=zDepth, sig='Fx_sec(z)', epsTolP=0.060, atol=0.21, printStats=True, test=test, factor=1e6)
    compareSigWithStats(My_sec0, My_sec1, t=zDepth, sig='My_sec(z)', epsTolP=0.030, atol=0.05, printStats=True, test=test, factor=1e8)

    # --- Monopile Section loads as function of time and depth
    vTime = dfRef['Time_[s]']
    IZ = [int(3*len(zDepth)/6)-5, int(2*len(zDepth)/6), int(1*len(zDepth)/6), 0] # Sea bed is zero
    for ii, iz in enumerate(IZ):
        compareSigWithStats(F_sec[0,iz,:], F_secRef[0,iz,:], t=vTime, sig=f'Fx_sec{iz}', epsTolP=35.5, atol=0.33, printStats=True, test=test, factor=1e6)
        compareSigWithStats(F_sec[4,iz,:], F_secRef[4,iz,:], t=vTime, sig=f'My_sec{iz}', epsTolP=35.0, atol=11.0, printStats=True, test=test, factor=1e6)
    # TODO FKz Not ready

    compareSigWithStats(dfOut, dfRef, sig='Wave1Elev_[m]', epsTolP=0.001, atol=0.0003, printStats=True, test=True)
    compareSigWithStats(dfOut, dfRef, sig='HydroFxi_[N]',  epsTolP=0.03,  atol=0.22, printStats=True, test=True, factor=1e6)



if __name__ == '__main__':
    #PLOT = True # KEEP ME FOR EASY DEBUG
    PLOT = False # KEEP ME FOR EASY DEBUG
#     TEST=False
    TEST=True
    test_monopile_tower_IEA(plot=PLOT, test=TEST)  # Was used to develop monopile only 
    test_monopile_only_MT100(plot=PLOT, test=TEST) # Need updating of main code in windturbine and merging of "debug_*" functions
    #test_floating_tower_TS(plot=PLOT, test=TEST)   # Was used to develop tower only
    #plt.show()
