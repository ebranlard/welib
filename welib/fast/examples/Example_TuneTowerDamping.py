""" 
Example to tune the ElastoDyn tower damping values (present in the ElastoDyn tower file) 
to obtain the desired damping ratio of the global modes (of the full structure).

NOTE: the method is approximate.

"""
import os
import numpy as np
from welib.fast.tuning import tuneTowerDamping

scriptDir = os.path.dirname(__file__)

fstFile = os.path.join(scriptDir,'../../../data/NREL5MW/Main_Onshore.fst')

zeta_target = [0.02, 0.02, 0.03, 0.04]  # Damping ratios of global modes (not in %)
modes       = ['FA1','FA2','SS1','SS2'] # Mode names corresponding to zeta_target

# Get tuned tower damping ratio for each modes
zt = tuneTowerDamping(fstFile, zeta_target=zeta_target, modes=modes, verbose=False)

if __name__ == '__main__':
    for k,v in zt.items():
        print('TwrDmp{}: {:.4f}%'.format(k,v*100))

if __name__ == '__test__':
    np.testing.assert_almost_equal(zt['FA1']*100, 5.3798, 4)
    np.testing.assert_almost_equal(zt['SS1']*100, 8.0806, 4)
