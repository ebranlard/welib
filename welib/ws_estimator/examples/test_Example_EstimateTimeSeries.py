""" 
These examples illustrate how to perform wind speed estimation on a full time series based on an OpenFAST input file (expecting the signals 'Wind1VelX' and 'RtAeroMxh' in the out file)

Look at estimateTimeSeriesFromOF if you want to implement a similar estimation on a different time series.

"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Local 
import welib.weio as weio
from welib.tools.colors import python_colors, fColrs, lighten_color
from welib.ws_estimator.tabulated import TabulatedWSEstimator

import pytest

scriptDir = os.path.dirname(__file__)

def estimate_time_series(fstFilename, aeroMapFile, operFile, tRange=None, method='crossing-oper', relaxation=0.9):
    """ """
    wse = TabulatedWSEstimator(fstFile=fstFilename, operFile = operFile, aeroMapFile=aeroMapFile)
    #print(wse)
    # --- Time series estimation
    df = wse.estimateTimeSeriesFromOF(fstFilename.replace('.fst','.outb'), relaxation=relaxation, tRange=tRange, method=method, debug=False)
    fig, stats = wse.plotTimeSeriesEstimation()
    return df, stats


def test_onshore_nrel5mw(test=True):
    aeroMapFile = os.path.join(scriptDir, '../../../data/NREL5MW/NREL5MW_CPCTCQ.txt')
    fstFilename = os.path.join(scriptDir, '../../../data/NREL5MW/onshore/Hat.fst')
    operFile    = os.path.join(scriptDir, '../../../data/NREL5MW/NREL5MW_Oper.csv')
    df, stats = estimate_time_series(fstFilename, aeroMapFile, operFile, tRange=None, method='crossing-oper')
    #print(stats)
    np.testing.assert_array_less(stats['WS']['eps']   ,4)
    np.testing.assert_array_less(stats['Qaero']['eps'],4)


def test_monopile_iea22(test=True):
    fstFilename = os.path.join(scriptDir, '../../../data/IEA-22-280-RWT/Jonswap/OF_F3T1S1_H1A1_Hs=8.1_Tp=12.7.fst')
    aeroMapFile = os.path.join(scriptDir, '../../../data/IEA-22-280-RWT/IEA-22-280-RWT/IEA-22-280-RWT_Cp_Ct_Cq.rpf')
    operFile    = os.path.join(scriptDir, '../../../data/IEA-22-280-RWT/IEA-22-280-RWT/IEA-22-280-RWT_OperOpenFAST.csv')
    df, stats = estimate_time_series(fstFilename, aeroMapFile, operFile, tRange=None, method='crossing-oper')
    #print(stats)
    np.testing.assert_array_less(stats['WS']['eps']   ,35)
    np.testing.assert_array_less(stats['Qaero']['eps'],5)


if __name__ == '__main__':
    test_onshore_nrel5mw(test=False)
    test_monopile_iea22(test=False)
    plt.show()

