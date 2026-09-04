"""
This example illustrates how to debug a single estimate.

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

def test_debug():
    # ---

    aeroMapFile = os.path.join(scriptDir, '../../../data/NREL5MW/NREL5MW_CPCTCQ.txt')
    fstFilename = os.path.join(scriptDir, '../../../data/NREL5MW/onshore/Hat.fst')
    operFile    = os.path.join(scriptDir, '../../../data/NREL5MW/NREL5MW_Oper.csv')
    wse = TabulatedWSEstimator(fstFile=fstFilename, operFile = operFile, aeroMapFile=aeroMapFile)

    method= 'crossing-oper'
    #method= 'crossing'
    omega   = 1       # [rad/s]
    Qa      = 1000000 # [Nm]
    pitch   = 0       # [deg]
    WS_last = 0.1
    WSest, info = wse.estimate(Qa, omega=omega, pitch=pitch,  WS0=WS_last, relaxation=0, WSavg=None, debug=True, method=method, deltaWSMax=10)
    print(info)
    fig = wse.debugPlot(info=info, HR=False)
    ax=fig.gca()
    ax.set_xlim([2,16])
    ax.legend(fontsize=12, ncol=2, loc='lower center')

    np.testing.assert_almost_equal(WSest, 6.5991, 3)



if __name__ == '__main__':
    test_debug()
    plt.show()
