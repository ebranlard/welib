"""Monopile/turbine digital twin using augmented Kalman estimation."""

import argparse
import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Welib
from welib.essentials import *

import welib.weio as weio
from welib.kalman.KF_MTNS import KalmanFilterMTNS
from welib.ws_estimator.tabulated import TabulatedWSEstimator
from welib.yams.section_loads_WT import YAMSSectionLoadCalculator



import pytest

scriptDir = os.path.dirname(__file__)

def main(fst_file, tmin=0, tmax=20, show=False,
    comp_file=None, hydro_shape_file=None, aero_map_file=None,
    oper_file=None, lin_file=None, method='YAMS', hacks=None,
    Tp=None,
    tRangeStats=None):
    
    # --- Main parameters

    # --- Default arguments:
    comp_file        = comp_file               or os.path.join(scriptDir, '_simulations/Waves/UserDefJonswap_Hs=8.1_Tp=12.7_h=34.csv')
    aero_map_file    = aero_map_file                or os.path.join(scriptDir, '_simulations/IEA-22-280-RWT/IEA-22-280-RWT_Cp_Ct_Cq.rpf')
    oper_file        = oper_file                    or os.path.join(scriptDir, '_simulations/IEA-22-280-RWT/IEA-22-280-RWT_OperOpenFAST.csv')
    hydro_shape_file = hydro_shape_file or os.path.join(scriptDir, '_data/IEAMonoPile_HydroShapeFunction_Hs=2.5_Tp=10.csv')

    if tRangeStats is None:
         tRangeStats = [tmin, tmax]

    base = os.path.splitext(fst_file)[0] + '_DigitalTwin'

    # --- Read reference output DataFrame 
    out_file = fst_file.replace('.fst', '.outb')
    if not os.path.exists(out_file):
        raise FileNotFoundError(out_file)
    df_ref = weio.read(out_file).toDataFrame()
    df_ref = df_ref[(df_ref['Time_[s]'] >= tmin) & (df_ref['Time_[s]'] <= tmax)]

    nUnderSamp=10
    df_ref=df_ref.iloc[::nUnderSamp,:]                      # reducing sampling
    df_ref.reset_index(inplace=True)

    # Read output file and ensure all mapped columns exist (e.g. for onshore case)
    missing_cols = ['Wave1Elev_[m]', 'HydroFxi_[N]', 'HydroMyi_[N-m]', '-ReactMYss_[N*m]', '-ReactFXss_[N]']
    for col in missing_cols:
        if col not in df_ref.columns:
            WARN('Missing column '+col)
            df_ref[col] = 0.0

    # --------------------------------------------------------------------------------}
    # --- Kalman filter estimation 
    # --------------------------------------------------------------------------------{
    # --- Wind speed estimator (reads tabulated aerodynamic data)
    wse = TabulatedWSEstimator(fstFile=fst_file, operFile=oper_file, aeroMapFile=aero_map_file)
    KF = KalmanFilterMTNS(WSE=wse, hacks=hacks)
    KF.setup_matrices(fst_file, 
                      comp_file=comp_file, hydro_shape_file=hydro_shape_file, Tp=Tp,
                      dfTime=df_ref['Time_[s]'].values, method=method, lin_file=lin_file)

    # --- Loading "Measurements"
    # - Reference file is opened
    # - Measurements are extracted from it
    # - Other signals are extracted from the file, for comparison with estimates. These are referred as "clean" values
    # - Estimate sigmas from measurements (overriden in next section)
    KF.loadMeasurements(measFile=df_ref, tRange=[tmin,tmax], colMap=KF.colMap, timeCol='Time_[s]', raiseIfAbsent=True)

    KF.X_clean['dq_h'] = np.gradient(KF.X_clean['q_h'], KF.dt)
    # --- Storage for plot
    KF.prepareTimeStepping()
    # --- Process and measurement covariances
    dt_ref = 0.01 # NOTE: Q change with dt
    sigs = {'x':{}, 'y':{}, 'Q':{}}
    #sigs['y']['TTacc'] = np.sqrt(1e-3)
    #sigs['Q']['q_s']   = np.sqrt(KF.dt/dt_ref * 1e-6)
    #sigs['Q']['q_p']   = np.sqrt(KF.dt/dt_ref * 1e-6)
    #sigs['Q']['qd_s']  = np.sqrt(KF.dt/dt_ref * 1e-3)
    #sigs['Q']['qd_p']  = np.sqrt(KF.dt/dt_ref * 1e-6)
    #sigs['Q']['q_h']   = np.sqrt(KF.dt/dt_ref * 1e-6)
    #sigs['Q']['qd_h']  = np.sqrt(KF.dt/dt_ref * KF.Sw)
    KF.setupCovariances(
            sigs=sigs,
            useDt=False, Pidentity=True, verbose=False)
    # TODO use sigs above instead
    KF.R[:] = np.diag([1e-3, 2.7e-7, 1e-5, 1e-2, 1e-2, 1e-3, 1e4])
    KF.Q[:] = np.diag([1e-6, 1e-6, 1e-5, 1e-6, 1e-5, 1e-6, 1e-5, 1e-5,
                       1e-6, 0.0001, 1e10])

    # --- Prepare measurements - Create noisy measurements
    KF.setYFromClean(R=KF.R, NoiseRFactor=0)


    # --------------------------------------------------------------------------------}
    # --- Section loads "ideal", everything prescribed
    # --------------------------------------------------------------------------------{
#     print('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SECTION LOADS USING DF_REF')
#     YSL = YAMSSectionLoadCalculator(fstFile=fst_file, WT=KF.WT)
#     with Timer('Section Loads ref'):
#         # NOTE: we cannot use KF.df as the columns have been renamed
#         # NOTE: this will trigger a calculation of the wave elevation
#         df_sl_ref, _ = YSL.fromDF(df_ref, useTopLoadsFromDF=KF.hacks['SL_cleanFtop'], useInterfaceLoadsFromDF=False, accMissing='warn')
#     df_sl_ref.to_outb(base + '_SectionLoads_ref.outb')
#     print('Export:', base + '_SectionLoads_ref.outb')
    # --------------------------------------------------------------------------------}
    # --- Time Loop 
    # --------------------------------------------------------------------------------{
    with Timer('KF time loop'):
        KF.timeLoop()
    df_sl = KF.dfOut
    file_sl = base + '_SectionLoads_KF_timeloop.outb'
    df_sl.to_outb(file_sl)
    print('Export:', file_sl)

    # --------------------------------------------------------------------------------}
    # --- Section Loads as Postpro --- Method 1 "Kalman Filter"
    # --------------------------------------------------------------------------------{
    print('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SECTION LOADS From KF')
    df_sl = KF.calc_sectionLoads(clean=True)
    file_sl = base + '_SectionLoads_KF.outb'
    df_sl.to_outb(file_sl)
    print('Export:', file_sl)

    # --------------------------------------------------------------------------------}
    # --- Section Loads as Postpro --- Method 2 "YSL"
    # --------------------------------------------------------------------------------{
#     print('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SECTION LOADS USING CLEAN')
#     clean = True
#     # Section loads post-processing using the pre-configured WT directly
#     YSL = YAMSSectionLoadCalculator(fstFile=fst_file, WT=KF.WT)
#     df_in2 = YSL.emptyInputDF(len(KF.time), inputFrame='R_xs', units=True)
#     df_in = pd.DataFrame({'Time_[s]': KF.time})
#     units = {'Sg':   ('[m]', '[m/s]', '[m/s^2]'),
#              'Sw':   ('[m]', '[m/s]', '[m/s^2]'),
#              'Hv':   ('[m]', '[m/s]', '[m/s^2]'),
#              'R':    ('[rad]', '[rad/s]', '[rad/s^2]'),
#              'P':    ('[rad]', '[rad/s]', '[rad/s^2]'),
#              'Y':    ('[rad]', '[rad/s]', '[rad/s^2]'),
#              'TFA1': ('[m]', '[m/s]', '[m/s^2]'),
#              'TSS1': ('[m]', '[m/s]', '[m/s^2]'),
#              'Yaw': ('[rad]', '[rad/s]', '[rad/s^2]')}
#     for dof_OF, suffixes in units.items():
#         for prefix, suffix in zip(['Q_', 'QD_', 'QD2_'], suffixes):
#             df_in[prefix + dof_OF + '_' + suffix] = 0.0
#     for name in ['Fadd_R_ys', 'Fadd_R_zs', 'Madd_R_ys', 'Madd_R_zs']:
#         df_in[name] = 0.0
# 
#     if clean:
#         for dof_OF, state in [('Sg', 'x'), ('P', 'phi_y'), ('TFA1', 'q_FA1')]:
#             suffixes = units[dof_OF]
#             df_in['Q_'   + dof_OF + '_' + suffixes[0]] = KF.X_clean[state]
#             df_in['QD_'  + dof_OF + '_' + suffixes[1]] = KF.XD_clean['d' + state]
#             df_in['QD2_' + dof_OF + '_' + suffixes[2]] = KF.XD_clean['dd'+ state]
# #         df_in['Madd_R_xs'] = KF.X_clean['Qaero']
#         df_in['Fadd_R_xs'] = KF.S_clean['Thrust']
#         df_in['Madd_R_xs'] = 0
#     else:
#         for dof_OF, state in [('Sg', 'x'), ('P', 'phi_y'), ('TFA1', 'q_FA1')]:
#             suffixes = units[dof_OF]
#             df_in['Q_'   + dof_OF + '_' + suffixes[0]] = KF.X_hat[state]
#             df_in['QD_'  + dof_OF + '_' + suffixes[1]] = KF.XD_hat['d' + state]
#             df_in['QD2_' + dof_OF + '_' + suffixes[2]] = KF.XD_hat['dd'+ state]
#         df_in['Madd_R_xs'] = KF.X_hat['Qaero']
#         df_in['Fadd_R_xs'] = KF.S_hat['Thrust']
# 
#     with Timer('Section Loads'):
#         df_sl, _ = YSL.fromDF(df_in, useTopLoadsFromDF=KF.hacks['SL_cleanFtop'], useInterfaceLoadsFromDF=False)
#     if clean:
#         file_sl = base + '_SectionLoads_clean.outb'
#     else:
#         file_sl = base + '_SectionLoads_est.outb'
#     df_sl.to_outb(file_sl)
#     print('Export:', file_sl)
        
    # --------------------------------------------------------------------------------}
    # --- Plot
    # --------------------------------------------------------------------------------{
    # Save all estimation results comparison
#     df_kf_all = KF.saveOutputs(filename=base + '_KF.outb', fmt='outb')
#     print('Export:', base + '_KF.outb')
#     pd.concat([df_kf_all, pd.DataFrame(df_sl)], axis=1).to_csv(base + '.csv', index=False)
#     
    statsDict = {}    
    # Plot results
    try:
        fig = KF.plot_X( printStats=True, tRangeStats=tRangeStats, statsDict=statsDict)
        plt.savefig(base + '_KF_X.png')
        KF.plot_Y()
        plt.savefig(base + '_KF_Y.png')
        fig = KF.plot_S(printStats=True, tRangeStats=tRangeStats, statsDict=statsDict)
        plt.savefig(base + '_KF_S.png')
    except Exception as e:
        FAIL('Plotting using KF plot functions failed:'+str(e))
    # KF.plot_P()
    # KF.plot_K()
    # KF.plot_innovation()
    
    if show:
        plt.show()
    return KF, df_ref, df_sl


def test_monopile_tower(test=True):
    if os.getenv('GITHUB_ACTIONS') == 'true':
        NOTE('test Offshore FTNS not ready yet')
        pytest.skip("Skipping local-only test on GitHub Actions")
    fst_file         = os.path.join(scriptDir, '_simulations/06_Jonswap/OF_F3T1S1_H1A1_Hs=8.1_Tp=12.7.fst')
    lin_file         = os.path.join(scriptDir, '_simulations/00_EVA/OF_F3T1S1_H1A1_OnlyWriteOutputs.1.lin')
    comp_file        = os.path.join(scriptDir, '_simulations/Waves/UserDefJonswap_Hs=8.1_Tp=12.7_h=34.csv')
    aero_map_file    = os.path.join(scriptDir, '_simulations/IEA-22-280-RWT/IEA-22-280-RWT_Cp_Ct_Cq.rpf')
    oper_file        = os.path.join(scriptDir, '_simulations/IEA-22-280-RWT/IEA-22-280-RWT_OperOpenFAST.csv')
    hydro_shape_file = os.path.join(scriptDir, '_data/IEAMonoPile_HydroShapeFunction_Hs=2.5_Tp=10.csv')

    #hacks = {'thrust':'clean', 'WSE':'clean_inputs', 'SL_cleanQ':True, 'SL_cleanFtop':True, 'SL_cleanEtaDot':True, 'SL_cleanP':False} # Super hack
    hacks = {'thrust':'clean', 'WSE':'clean_inputs', 'SL_cleanQ':True, 
             'SL_cleanFtop':True, 'SL_cleanEtaDot':True, 'SL_cleanP':True} # Super hack
    #hacks['WSE'] = 'clean_inputs'
    #hacks['thrust'] = 'clean'
    show=True
    if test:
        tRange = [150, 170]
        show=False
    else:
        tRange = [150, 170]
    main(fst_file=fst_file, lin_file=lin_file, 
         comp_file=comp_file,  hydro_shape_file=hydro_shape_file, Tp=12.7,
         aero_map_file=aero_map_file, oper_file=oper_file,
         hacks=hacks, show=show,
         tmin=tRange[0], tmax=tRange[1]
         )

if __name__ == '__main__':
    if len(sys.argv)==1:
        test_monopile_tower(test=False)
        sys.exit(0)
    else:
        raise

    parser = argparse.ArgumentParser()
    parser.add_argument('fst_file', nargs='?', default='_simulations/06_Jonswap/OF_F3T1S1_H1A1_Hs=8.1_Tp=12.7.fst')
    parser.add_argument('--lin-file', type=str, default='_simulations/00_EVA/OF_F3T1S1_H1A1_OnlyWriteOutputs.1.lin')
    parser.add_argument('--tmin', type=float, default=150)
    parser.add_argument('--tmax', type=float, default=170)
    parser.add_argument('--show', action='store_true')
    parser.add_argument('--method', choices=['YAMS', 'OpenFAST'], default='YAMS')


    hacks = {'thrust':None, 'WSE':None}
    #hacks['WSE'] = 'clean_inputs'
    #hacks['thrust'] = 'clean'

    args = parser.parse_args()

    print('Arguments:',args)
    main(args.fst_file, tmin=args.tmin, tmax=args.tmax, show=args.show, method=args.method, lin_file=args.lin_file, hacks=hacks)



