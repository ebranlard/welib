""" 
Example to compute mode shapes and frequencies from a SubDyn model

NOTE: unfinished
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Local 
from welib.fast.subdyn import SubDyn
import weio
import os




if __name__ == '__main__':
    MyDir=os.path.dirname(__file__)
    filename=os.path.join(MyDir,'../../../data/NREL5MW/data/NRELOffshrBsline5MW_OC4Jacket_SubDyn.dat')

    # sd = weio.read(filename)
    # print(sd.keys())
    #sd=SubDyn(filename)
    #FEM = sd.beamFEM()
    #print(FEM.keys())
    #print(np.around(FEM['freq'][0:7],2))
    #except:
    #    print(' ')
    #    print(' ')
    print(' SUBDYN EXAMPLE NOT READY')
    #    print(' ')
    #    print(' ')
    #pass
