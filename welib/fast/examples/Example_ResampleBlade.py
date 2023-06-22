import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Local 
import welib.fast.aerodyn as ad

# Get current directory so this script can be called from any location
scriptDir = os.path.dirname(__file__)

old=os.path.join(scriptDir, '../../../data/NREL5MW/5MW_Baseline/NRELOffshrBsline5MW_AeroDyn_blade.dat')
new=os.path.join(scriptDir, '../../../data/NREL5MW/5MW_Baseline/NRELOffshrBsline5MW_AeroDyn_blade_resampled.dat')

# resample blade
r_new = np.linspace(0, 61.5, 100)
ad.resampleBlade(r_new, old, new)



if __name__ == '__main__':
    pass
if __name__ == '__test__':
    os.remove(new)
