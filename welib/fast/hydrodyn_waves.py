import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Local 
from welib.tools.clean_exceptions import *
from welib.weio import FASTInputFile


class Waves:

    def __init__(self, File, WtrDpth, MSL2SWL):
        """ 
        File : File content of hydrodyn input file
        """
        self.File         = File

        # Internal
        self.p={}
        self.m={}
        self.p['WtrDpth'] = WtrDpth
        self.p['MSL2SWL'] = MSL2SWL

    def __repr__(self):
        s='<{} object>:\n'.format(type(self).__name__)
        s+='|properties:\n'
        s+='|- File: (input file data)\n'
        s+='|parameters (p):\n'
        s+='| - WtrDpth: {}\n'.format(p.WtrDpth)
        s+='| - MSL2SWL: {}\n'.format(p.MSL2SWL)
        s+='|methods:\n'
        s+='|- init\n'
        return s

    # --------------------------------------------------------------------------------}
    # --- Functions 
    # --------------------------------------------------------------------------------{
    def init(self, initData=None, Gravity = 9.81, WtrDens=1025):
        """
        Initialize Waves model 
        """
        F = self.File
        if F['WaveMod']==0:
            self.p['WaveTime'] = np.array([0,1,2])
            self.p['NStepWave'] = len(self.p['WaveTime'])
        else:
            TMax = F['WaveTMax']
            DT   = F['WaveDT']
            self.p['WaveTime'] = np.arange(0,TMax+DT/2,DT)
            self.p['NStepWave'] = np.ceil(F['WaveTMax']/F['WaveDT'] ).astype(int)

    def calcOutput(self, t, x=None, xd=None, xo=None):
        """ NOTE: """
        pass

if __name__ == '__main__':
    import sys
    if len(sys.argv)>=1:
        filename=sys.argv[1]
    else:
        filename='_SparNoRNA_HD_RefH.dat'
