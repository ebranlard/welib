""" 
Example to compute mode shapes and frequencies from a SubDyn model

NOTE: unfinished
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Local 
from welib.fast.subdyn import SubDyn
from welib.tools.tictoc import Timer

MyDir=os.path.dirname(__file__)

if __name__ == '__main__':
    np.set_printoptions(precision=5, linewidth=400)

    filename=os.path.join(MyDir,'../../../data/NREL5MW/5MW_Baseline/NRELOffshrBsline5MW_OC4Jacket_SubDyn.dat'); TP=(0,0,20.15); gravity=9.81
    #filename=os.path.join(MyDir,'../../../data/SubDyn/Twr.dat'); TP=(0,0,102.7930); gravity=9.81;
#     filename=os.path.join(MyDir,'../../../data/SubDyn/Jacket.dat');    TP=(0,0,22); gravity=9.81;
#     filename=os.path.join(MyDir,'../../../data/SubDyn/JacketTwr.dat'); TP=(0,0, 102.7930); gravity=9.81;
    #filename=os.path.join(MyDir,'../../../data/SubDyn/SD_Cable_5Joints.dat'); TP=(0,0,0); gravity=9.81;
    #filename=os.path.join(MyDir,'../../../data/SubDyn/TwrSmall.dat'); TP=(0,0,100); gravity=9.81;
    #filename=os.path.join(MyDir,'../../../data/SubDyn/SD_Rigid/SD_Rigid.dat'); TP=(0,0,30); gravity=9.81;

    sd = SubDyn(filename)
    model = sd.init(TP=TP, gravity=gravity)
    with Timer('YAML'):
        sd.toYAML(filename.replace('.dat','.SD.python.yaml'))
    with Timer('JSON'):
        model.toJSON(filename.replace('.dat','.SD.python.json'))

#         print(model.nModesCB)
#         print(model.freq)
#         print(model.f_CB)
#         print(model.f_G)
#         print(model.M_O)


if __name__ == '__test__':
    # --- Test Execution for tower with Sub divisions, Timoshenko Beams, concentrated masses
    filename=os.path.join(MyDir,'../../../data/SubDyn/TwrSmall.dat'); TP=(0,0,100); gravity=9.81;
    sd  = SubDyn(filename)
    model = sd.init(TP= TP, gravity = gravity)
    # Test export
    filebase = os.path.splitext(os.path.basename(filename))[0]
    outbase  = os.path.join(MyDir, '_' + filebase )
    sd.toYAML (outbase+'.python.yaml')
    model.toJSON(outbase+'.python.json')
    try:
        os.delete(outbase+'.yaml')
        os.delete(outbase+'.json')
    except:
        pass
    # Test data
    np.testing.assert_almost_equal(model.M_O [0,0]/1e5, 3.85591, 4)
    np.testing.assert_almost_equal(model.M_O [3,3]/1e9, 1.29432, 4)
    np.testing.assert_almost_equal(model.M_O [5,5]/1e6, 2.03386, 4)
    np.testing.assert_almost_equal(model.freq[0:5], [0.47074,  0.47074,  2.69047,  2.69128,  7.04054], 4)
    np.testing.assert_almost_equal(model.f_CB,      [ 2.72048, 2.72135, 6.96018, 6.96337,13.14011,13.16074,17.27573,20.79389],4)
    np.testing.assert_almost_equal(model.f_G,       [ 0.47338, 0.47338, 4.35527, 4.35587, 9.6198 ,13.50852],4)

    # --- Test for Pretension Cables
    filename=os.path.join(MyDir,'../../../data/SubDyn/SD_Cable_5Joints.dat'); TP=(0,0,0); gravity=9.81;
    sd = SubDyn(filename)
    model = sd.init(TP=TP, gravity=gravity)

    np.testing.assert_almost_equal(model.M_O [0,0]/1e5, 1.8449396, 4)
    np.testing.assert_almost_equal(model.M_O [3,3]/1e4, 9.555369, 4)
    np.testing.assert_almost_equal(model.M_O [5,5]/1e6, 6.01656196, 4)
    np.testing.assert_almost_equal(model.freq[0:5], [3.84067e-12,1.02254e-01,1.02254e-01,7.66803e-01,7.66803e-01], 3)
    np.testing.assert_almost_equal(model.f_CB,      [],4)
    np.testing.assert_almost_equal(model.f_G,    [ 0.,      0.10225, 0.10225, 0.7668,  0.7668, 14.12842] ,4)
