""" Unit tests of the monopile/turbine digital twin (KF_MTNS). Tolerances are loose around values obtained with tRange=[210,250] (and [100,400] for longer window) """
import numpy as np
import os
import pytest
from welib.kalman.KF_MTNS import MonopileOnly, FullStructure_NoWave_NoMonopileDOFs, FullStructure_WithWave
from welib.essentials import *

tR_short = [210, 250]
tR_long  = [100, 400]  # Optimal convergence window: 300s allows full convergence without tail noise

def check(stats, R2=None, eps=None, tolR2=0.02, tolEps=0.5):
    """ Compare R2 and eps (percent) of the channels with reference values """
    for k, v in (R2 or {}).items():
        np.testing.assert_allclose(stats[k]['R2'], v, atol=tolR2, err_msg='R2 of '+k)
    for k, v in (eps or {}).items():
        np.testing.assert_allclose(stats[k]['eps'], v, atol=tolEps, err_msg='eps of '+k)

# --------------------------------------------------------------------------------}
# --- Monopile only 
# --------------------------------------------------------------------------------{
def test_case1_MonopileOnly():
    if os.getenv('GITHUB_ACTIONS') == 'true':
        NOTE('test Offshore FTNS not ready yet')
        pytest.skip("MTNS not ready")
        return
    _, s = MonopileOnly(tRange=tR_short, sl_mode='linear')
    check(s, R2=dict(x=0.954, q_h=0.975, dq_h=0.950, Fx_sb=0.620, Fx_h=0.943),
             eps=dict(x=2.67, q_h=1.69, Fx_sb=8.05, My_sb=11.2))

def test_case1_MonopileOnly_OpenFAST():
    if os.getenv('GITHUB_ACTIONS') == 'true':
        NOTE('test Offshore FTNS not ready yet')
        pytest.skip("MTNS not ready")
        return
    _, s = MonopileOnly(tRange=tR_short, method='OpenFAST', sl_mode='linear')
    check(s, R2=dict(x=0.9995, q_h=0.9885, dq_h=0.9856, Fx_h=0.950), eps=dict(x=0.30, q_h=0.97, Fx_h=2.70))


# --------------------------------------------------------------------------------}
# --- Monopile and Tower 
# --------------------------------------------------------------------------------{
def test_case3_MonopileAndTower():
    if os.getenv('GITHUB_ACTIONS') == 'true':
        NOTE('test Offshore FTNS not ready yet')
        pytest.skip("MTNS not ready")
        return
    _, s = FullStructure_WithWave(tRange=tR_long, sl_mode='linear')  # Use optimal [100,400] window
    check(s, R2=dict(x=0.996, q_FA1=0.980, WS=0.990, PtfmIMUAx=0.973, NcIMUAx=0.933, Fx_sb=0.791, My_sb=0.955, Fx_i=0.957, My_i=0.973),
             eps=dict(x=0.80, q_FA1=2.06, Fx_sb=4.74, My_sb=2.36))

def test_case3_Hybrid():
    if os.getenv('GITHUB_ACTIONS') == 'true':
        NOTE('test Offshore FTNS not ready yet')
        pytest.skip("MTNS not ready")
        return
    _, s = FullStructure_WithWave(tRange=tR_long, method='Hybrid', sl_mode='linear')  # Use optimal [100,400] window
    check(s, R2=dict(x=0.995, q_FA1=0.983, NcIMUAx=0.939, Fx_sb=0.788, My_sb=0.956, My_i=0.977), eps=dict(x=0.85, q_FA1=1.74))


# --------------------------------------------------------------------------------}
# --- Tower Only 
# --------------------------------------------------------------------------------{
def test_case2_TowerOnly():
    if os.getenv('GITHUB_ACTIONS') == 'true':
        NOTE('test Offshore FTNS not ready yet')
        pytest.skip("MTNS not ready")
        return
    _, s = FullStructure_NoWave_NoMonopileDOFs(tRange=tR_short, sl_mode='linear')
    check(s, R2=dict(psi=0.994, dpsi=1.0, WS=0.768, Thrust=0.895, Fx_i=0.467, My_i=0.754),
             eps=dict(q_FA1=9.29, WS=7.53, Fx_sb=8.24, My_i=6.88))


def test_case2_TowerOnly_OpenFAST():
    if os.getenv('GITHUB_ACTIONS') == 'true':
        NOTE('test Offshore FTNS not ready yet')
        pytest.skip("MTNS not ready")
        return
    _, s = FullStructure_NoWave_NoMonopileDOFs(tRange=tR_short, method='OpenFAST', sl_mode='linear')
    check(s, R2=dict(q_FA1=0.888, WS=0.767, Thrust=0.895, My_i=0.742), eps=dict(q_FA1=4.52, My_i=7.07))

if __name__ == '__main__':
#     test_case3_MonopileAndTower()
    test_case1_MonopileOnly()
