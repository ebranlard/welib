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


def compute_metrics(KF, df_ref, df_sl):
    print("\n" + "="*70)
    print("           DIGITAL TWIN ESTIMATION PERFORMANCE METRICS")
    print("="*70)
    
    metrics = {}
    
    # 1. Wind Speed (RtVAvgxh_[m/s] vs WS)
    if 'RtVAvgxh_[m/s]' in df_ref.columns:
        metrics['Wind Speed'] = (df_ref['RtVAvgxh_[m/s]'].values, KF.S_hat['WS'].values)
        
    # 2. Qaero (RtFldMxh_[N-m] vs Qaero)
    if 'RtFldMxh_[N-m]' in df_ref.columns:
        metrics['Qaero'] = (df_ref['RtFldMxh_[N-m]'].values, KF.X_hat['Qaero'].values)
        
    # 3. Thrust (RtAeroFxh_[N] vs Thrust)
    if 'RtAeroFxh_[N]' in df_ref.columns and 'Thrust' in KF.S_hat.columns:
        metrics['Thrust'] = (df_ref['RtAeroFxh_[N]'].values, KF.S_hat['Thrust'].values)
        
    # 4. q_FA1 (Q_TFA1_[m] vs q_FA1)
    if 'Q_TFA1_[m]' in df_ref.columns:
        metrics['q_FA1'] = (df_ref['Q_TFA1_[m]'].values, KF.X_hat['q_FA1'].values)
        
    # 5. Mid tower bending moment (TwHt5MLyt_[kN-m] vs TwHt5MLyt_[kN-m])
    if 'TwHt5MLyt_[kN-m]' in df_ref.columns and 'TwHt5MLyt_[kN-m]' in df_sl.columns:
        metrics['Mid tower bending moment'] = (df_ref['TwHt5MLyt_[kN-m]'].values, df_sl['TwHt5MLyt_[kN-m]'].values)
        
    # 6. Interface bending moment (TwHt1MLyt_[kN-m] vs TwHt1MLyt_[kN-m])
    if 'TwHt1MLyt_[kN-m]' in df_ref.columns and 'TwHt1MLyt_[kN-m]' in df_sl.columns:
        metrics['Interface bending moment'] = (df_ref['TwHt1MLyt_[kN-m]'].values, df_sl['TwHt1MLyt_[kN-m]'].values)
        
    # 7. Interface shear force (TwHt1FLxt_[kN] vs TwHt1FLxt_[kN])
    if 'TwHt1FLxt_[kN]' in df_ref.columns and 'TwHt1FLxt_[kN]' in df_sl.columns:
        metrics['Interface shear force'] = (df_ref['TwHt1FLxt_[kN]'].values, df_sl['TwHt1FLxt_[kN]'].values)
        
    # 8. Hydro Fx (HydroFxi_[N] vs Fx_h)
    if 'HydroFxi_[N]' in df_ref.columns and 'Fx_h' in KF.S_hat.columns:
        metrics['Hydro Fx'] = (df_ref['HydroFxi_[N]'].values, KF.S_hat['Fx_h'].values)
        
    # 9. Wave elevation (Wave1Elev_[m] vs eta)
    if 'Wave1Elev_[m]' in df_ref.columns and 'eta' in KF.S_hat.columns:
        metrics['Wave elevation'] = (df_ref['Wave1Elev_[m]'].values, KF.S_hat['eta'].values)
        
    # 10. Sea bed moment (-ReactMYss_[N*m] vs M1N1MKye_[N*m])
    if '-ReactMYss_[N*m]' in df_ref.columns and 'M1N1MKye_[N*m]' in df_sl.columns:
        metrics['Sea bed moment'] = (df_ref['-ReactMYss_[N*m]'].values, df_sl['M1N1MKye_[N*m]'].values)
        
    min_baselines = {
        'Wind Speed': 1.0,
        'Qaero': 1000.0,
        'Thrust': 1000.0,
        'q_FA1': 0.1,
        'Mid tower bending moment': 1000.0,
        'Interface bending moment': 1000.0,
        'Interface shear force': 1000.0,
        'Hydro Fx': 1000.0,
        'Wave elevation': 0.1,
        'Sea bed moment': 1000.0
    }
        
    for name, (ref_v, est_v) in metrics.items():
        n = min(len(ref_v), len(est_v))
        ref_v = ref_v[:n]
        est_v = est_v[:n]
        
        # Calculate mean relative error (with a baseline to avoid div by zero)
        baseline = np.mean(np.abs(ref_v))
        if baseline < 1e-6:
            baseline = np.std(ref_v)
        baseline = max(baseline, min_baselines.get(name, 1.0))
            
        mre = np.mean(np.abs(ref_v - est_v)) / baseline * 100
        
        # Also compute R-squared
        ss_res = np.sum((ref_v - est_v)**2)
        ss_tot = np.sum((ref_v - np.mean(ref_v))**2)
        r2 = 1 - (ss_res / ss_tot) if ss_tot > 1e-8 else 1.0
        
        print(f"{name:30s} | Mean Rel Error: {mre:7.2f} % | R^2: {r2:6.3f}")
        
    print("="*70 + "\n")


def main(fst_file, tmin=0, tmax=20, show=False,
    comp_file=None, hydro_shape_file=None, aero_map_file=None,
    oper_file=None, lin_file=None, method='YAMS', hacks=None,
    tRangeStats=None):
    
    # --- Main parameters

    # --- Default arguments:
    comp_file        = comp_file               or os.path.join(scriptDir, '_simulations/Waves/UserDefJonswap_Hs=8.1_Tp=12.7_h=34.csv')
    aero_map_file    = aero_map_file                or os.path.join(scriptDir, '_simulations/IEA-22-280-RWT/IEA-22-280-RWT_Cp_Ct_Cq.rpf')
    oper_file        = oper_file                    or os.path.join(scriptDir, '_simulations/IEA-22-280-RWT/IEA-22-280-RWT_OperOpenFAST.csv')
    hydro_shape_file = hydro_shape_file or os.path.join(scriptDir, '_data/IEAMonoPile_HydroShapeFunction_Hs=2.5_Tp=10.csv')

    if tRangeStats is None:
         tRangeStats = [tmin, tmax]


    # --- Read reference output DataFrame 
    out_file = fst_file.replace('.fst', '.outb')
    if not os.path.exists(out_file):
        raise FileNotFoundError(out_file)
    df_ref = weio.read(out_file).toDataFrame()
    df_ref = df_ref[(df_ref['Time_[s]'] >= tmin) & (df_ref['Time_[s]'] <= tmax)]

    # Read output file and ensure all mapped columns exist (e.g. for onshore case)
    missing_cols = ['Wave1Elev_[m]', 'HydroFxi_[N]', 'HydroMyi_[N-m]', '-ReactMYss_[N*m]', '-ReactFXss_[N]']
    for col in missing_cols:
        if col not in df_ref.columns:
            df_ref[col] = 0.0

    # --------------------------------------------------------------------------------}
    # --- Kalman filter estimation 
    # --------------------------------------------------------------------------------{
    # --- Wind speed estimator (reads tabulated aerodynamic data)
    wse = TabulatedWSEstimator(fstFile=fst_file, operFile=oper_file, aeroMapFile=aero_map_file)
    KF = KalmanFilterMTNS(WSE=wse, hacks=hacks)
    KF.setup_matrices(fst_file, comp_file=comp_file, hydro_shape_file=hydro_shape_file,
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
            useDt=False, Pidentity=True, verbose=True)
    # TODO use sigs above instead
    KF.R[:] = np.diag([1e-3, 2.7e-7, 1e-5, 1e-2, 1e-2, 1e-3, 1e4])
    KF.Q[:] = np.diag([1e-6, 1e-6, 1e-5, 1e-6, 1e-5, 1e-6, 1e-5, 1e-5,
                       1e-6, 0.0001, 1e10])

    # --- Prepare measurements - Create noisy measurements
    KF.setYFromClean(R=KF.R, NoiseRFactor=0)


    # --------------------------------------------------------------------------------}
    # --- Time Loop 
    # --------------------------------------------------------------------------------{
    with Timer('KF time loop'):
        KF.timeLoop()

    # --------------------------------------------------------------------------------}
    # --- Calc Outputs
    # --------------------------------------------------------------------------------{
    # Section loads post-processing using the pre-configured WT directly
    ysl = YAMSSectionLoadCalculator(fstFile=fst_file, WT=KF.WT)
    df_in = pd.DataFrame({'Time_[s]': KF.time})
    units = {'Sg':   ('[m]', '[m/s]', '[m/s^2]'),
             'Sw':   ('[m]', '[m/s]', '[m/s^2]'),
             'Hv':   ('[m]', '[m/s]', '[m/s^2]'),
             'R':    ('[rad]', '[rad/s]', '[rad/s^2]'),
             'P':    ('[rad]', '[rad/s]', '[rad/s^2]'),
             'Y':    ('[rad]', '[rad/s]', '[rad/s^2]'),
             'TFA1': ('[m]', '[m/s]', '[m/s^2]'),
             'TSS1': ('[m]', '[m/s]', '[m/s^2]'),
             'Yaw': ('[rad]', '[rad/s]', '[rad/s^2]')}
    for dof_OF, suffixes in units.items():
        for prefix, suffix in zip(['Q_', 'QD_', 'QD2_'], suffixes):
            df_in[prefix + dof_OF + '_' + suffix] = 0.0
    for dof_OF, state in [('Sg', 'x'), ('P', 'phi_y'), ('TFA1', 'q_FA1')]:
        suffixes = units[dof_OF]
        df_in['Q_'   + dof_OF + '_' + suffixes[0]] = KF.X_hat[state]
        df_in['QD_'  + dof_OF + '_' + suffixes[1]] = KF.XD_hat['d' + state]
        df_in['QD2_' + dof_OF + '_' + suffixes[2]] = KF.XD_hat['dd'+ state]
    df_in['Madd_R_xs'] = KF.X_hat['Qaero']
    df_in['Fadd_R_xs'] = KF.S_hat['Thrust']
    for name in ['Fadd_R_ys', 'Fadd_R_zs', 'Madd_R_ys', 'Madd_R_zs']:
        df_in[name] = 0.0

    with Timer('Section Loads'):
        df_sl, _ = ysl.fromDF(df_in, useTopLoadsFromDF=KF.hacks['useTopLoadsFromDF'], useInterfaceLoadsFromDF=False)
        


    # --------------------------------------------------------------------------------}
    # --- Export and plot
    # --------------------------------------------------------------------------------{
    base = os.path.splitext(fst_file)[0] + '_DigitalTwin'
    
    # Save section loads as .csv and .outb
    from welib.weio.fast_output_file import writeDataFrame
    pd.DataFrame(df_sl).to_csv(base + '_SectionLoads.csv', index=False)
    writeDataFrame(pd.DataFrame(df_sl), base + '_SectionLoads.outb')

    # Save all estimation results comparison
    df_kf_all = KF.saveOutputs(base + '.outb', fmt='outb')
    pd.concat([df_kf_all, pd.DataFrame(df_sl)], axis=1).to_csv(base + '.csv', index=False)
    
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

    # Compute and display performance metrics
    compute_metrics(KF, df_ref, df_sl)
    
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

    hacks = {'thrust':'clean', 'WSE':'clean_inputs', 'useTopLoadsFromDF':True} # Super hack
    #hacks['WSE'] = 'clean_inputs'
    #hacks['thrust'] = 'clean'
    show=True
    if test:
        tRange = [150, 200]
        show=False
    else:
        tRange = [150, 170]
    main(fst_file=fst_file, lin_file=lin_file, 
         comp_file=comp_file,  hydro_shape_file=hydro_shape_file,
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


    hacks = {'thrust':None, 'WSE':None, 'useTopLoadsFromDF':True}
    #hacks['WSE'] = 'clean_inputs'
    #hacks['thrust'] = 'clean'

    args = parser.parse_args()

    print('Arguments:',args)
    main(args.fst_file, tmin=args.tmin, tmax=args.tmax, show=args.show, method=args.method, lin_file=args.lin_file, hacks=hacks)



