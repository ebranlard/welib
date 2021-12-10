""" 
Example to compute mode shapes and frequencies from a SubDyn model

NOTE: unfinished
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Local 
import welib.weio as weio
from welib.fast.subdyn import SubDyn
from welib.tools.tictoc import Timer
import os




if __name__ == '__main__':
    np.set_printoptions(precision=5, linewidth=400)
    MyDir=os.path.dirname(__file__)
    #filename=os.path.join(MyDir,'../../../data/NREL5MW/5MW_Baseline/NRELOffshrBsline5MW_OC4Jacket_SubDyn.dat')
    #filename=os.path.join(MyDir,'../../../data/NREL5MW/5MW_Baseline/NRELOffshrBsline5MW_OC4Jacket_SubDyn.dat')
    #filename=os.path.join(MyDir,'../../../data/SubDyn/Twr.dat'); TP=(0,0,102.7930); gravity=9.81;
    #filename=os.path.join(MyDir,'../../../data/SubDyn/Jacket.dat'); TP=(0,0,22); gravity=9.81;
    filename=os.path.join(MyDir,'../../../data/SubDyn/TwrSmall.dat'); TP=(0,0,100); gravity=9.81;
    #filename=os.path.join(MyDir,'../../../data/SubDyn/SD_Cable_5Joints.dat'); TP=(0,0,0); gravity=9.81;

    #sd = weio.read(filename)
    #print(sd.keys())
    sd=SubDyn(filename)
    #FEM = sd.beamFEM()
    FEM = sd.FEM(TP=TP, gravity=gravity)

    with Timer('YAML'):
        FEM.toYAML(filename.replace('.dat','.SD.python.yaml'))
    with Timer('JSON'):
        FEM.toJSON(filename.replace('.dat','.SD.python.json'))

#     for e in FEM.Elements:
#         #n1, n2 = e.nodes[0], e.nodes[1]
#         n1, n2 = e.nodeProps[0], e.nodeProps[1]
#         D1, D2 = n1.data['D'], n2.data['D']
#         print(e.ID, D1, D2)


#     print(FEM.keys())
#     print(np.around(FEM['freq'][0:7],2))
    #except:
    #    print(' ')
    #    print(' ')
    #print(' SUBDYN EXAMPLE NOT READY')
    #    print(' ')
    #    print(' ')
    #pass
