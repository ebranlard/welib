"""Monopile/turbine digital twin using augmented Kalman estimation."""

import argparse
import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Welib
from welib.essentials import *

import pytest

# from welib.kalman.KF_MTNS import MonopileOnly, FullStructure_NoWave_NoMonopileDOFs, FullStructure_WithWave
from welib.kalman.KF_MTNS import mainWrapper

tR_short = [210, 250]

scriptDir = os.path.dirname(__file__)

def check(stats, R2=None, eps=None, tolR2=0.02, tolEps=0.5):
    """ Compare R2 and eps (percent) of the channels with reference values """
    for k, v in (R2 or {}).items():
        np.testing.assert_allclose(stats[k]['R2'], v, atol=tolR2, err_msg='R2 of '+k)
    for k, v in (eps or {}).items():
        np.testing.assert_allclose(stats[k]['eps'], v, atol=tolEps, err_msg='eps of '+k)


def test_MonopileTower_WithWave(method='YAMS', tRange=tR_short, hacks=None, show=False, mode='linear'):
    if os.getenv('GITHUB_ACTIONS') == 'true':
        NOTE('test Offshore FTNS not ready yet')
        pytest.skip("MTNS not ready")
        return

    linFiles=[]
    linFiles += [os.path.join(scriptDir, '_simulations/00_EVA/OF_F3T1S1_H1A1.1.lin')]
    linFiles += [os.path.join(scriptDir, '_simulations/00_EVA/OF_F3T1S1_H1A1_Trim4mps.1.lin')]
    linFiles += [os.path.join(scriptDir, '_simulations/00_EVA/OF_F3T1S1_H1A1_Trim7mps.1.lin')]
    linFiles += [os.path.join(scriptDir, '_simulations/00_EVA/OF_F3T1S1_H1A1_Trim10mps.1.lin')]
    linFiles += [os.path.join(scriptDir, '_simulations/00_EVA/OF_F3T1S1_H1A1_Trim12mps.1.lin')]

    fstFile  = os.path.join(scriptDir, '_simulations/06_Jonswap/OF_F3T1S1_H1A1_Hs=8.1_Tp=12.7.fst'); 

    setup_opts = {
        'hydro_states': True,           # Wave estimation enabled
        'monopileDOFs': True,           # Use monopile DOFs (platform surge, pitch)
        'aero_est': True,               # Aerodynamic force estimation
        'q_FA1': True,                  # Tower first bending mode
        'bias_states': True,            # Bias states: Bx, Bphi for monopile
        'sl_mode': mode,             # Section loads mode
        # Explicit bias state covariances:
        'P0_Bx': 1.0,                   # Covariance for platform surge bias (Bx)
        'P0_Bphi': 1.0,                 # Covariance for platform pitch bias (Bphi)
        'P0_Bq': 0.1,                   # Not used (q_FA1 coupled to platform)
        'hybrid': ['KD', 'static', 'hydro'],  # Hybrid model config (blocks to replace from OpenFAST)
        'kh_scale': [1.0, 1.0],         # Wave hydro scale (x2 by default, no scale)
        'static_bias': [0, 0, 0],       # Static bias offset (platform forces)
    }

    _, stats = mainWrapper(fstFile=fstFile, method=method, hacks=hacks, setup_opts=setup_opts, tRange=tRange, show=show, linFiles=linFiles)
    check(stats, R2=dict(x=0.984, q_FA1=0.912, WS=0.969, PtfmIMUAx=0.977, NcIMUAx=0.920, Fx_sb=0.746, My_sb=0.766, Fx_i=0.726, My_i=0.831), eps=dict(x=1.6, q_FA1=4.28, Fx_sb=5.7, My_sb=5.1))




if __name__ == '__main__':

    # --- Monopile Jonswap Hs=8.1 Tp=12.7, Default
    for mode in ['linear', 'exact']:
        test_MonopileTower_WithWave(show=True, mode=mode)
    plt.show()
