import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
# Local 
import welib.weio as weio
from welib.fast.elastodyn import fitShapeFunction

scriptDir = os.path.dirname(__file__)

# --- Fit given shape functions
x  = np.linspace(0,1)
phi1 = 1-np.cos(x*np.pi/2) # dummy example of mode 1
phi2 = -26*x**2 + 27*x**4  # dummy example of mode 2
coeffs1, _, _ = fitShapeFunction(x, phi1, plot=True)
coeffs2, _, _ = fitShapeFunction(x, phi2, plot=True)
coeffs3, _, _ = fitShapeFunction(x, phi1, plot=False)
coeffs4, _, _ = fitShapeFunction(x, phi2, plot=False)

# --- Open an ElastoDyn tower file, replace the coefficients with the fitted ones
edFile = os.path.join(scriptDir,'../../../data/NREL5MW/5MW_Baseline/NRELOffshrBsline5MW_Onshore_ElastoDyn_Tower.dat')
twr = weio.read(edFile)
scoeffs1=['TwFAM1Sh(2)','TwFAM1Sh(3)','TwFAM1Sh(4)','TwFAM1Sh(5)','TwFAM1Sh(6)']
scoeffs2=['TwFAM2Sh(2)','TwFAM2Sh(3)','TwFAM2Sh(4)','TwFAM2Sh(5)','TwFAM2Sh(6)']
scoeffs3=['TwSSM1Sh(2)','TwSSM1Sh(3)','TwSSM1Sh(4)','TwSSM1Sh(5)','TwSSM1Sh(6)']
scoeffs4=['TwSSM2Sh(2)','TwSSM2Sh(3)','TwSSM2Sh(4)','TwSSM2Sh(5)','TwSSM2Sh(6)']
for scoeffs, coeffs in zip( (scoeffs1,scoeffs2, scoeffs3, scoeffs4), (coeffs1,coeffs2,coeffs3,coeffs4)):
    for s,c in zip(scoeffs, coeffs):
        twr[s] = c
# Write the New ElastoDyn file
twr.write('_EDtwrFitted.dat')


# --- Open an ElastoDyn blade file, replace the coefficients with the fitted ones
edFile = os.path.join(scriptDir,'../../../data/NREL5MW/5MW_Baseline/NRELOffshrBsline5MW_Blade.dat')
bld = weio.read(edFile)
scoeffs1=['BldFl1Sh(2)','BldFl1Sh(3)','BldFl1Sh(4)','BldFl1Sh(5)','BldFl1Sh(6)']
scoeffs2=['BldFl2Sh(2)','BldFl2Sh(3)','BldFl2Sh(4)','BldFl2Sh(5)','BldFl2Sh(6)']
scoeffs3=['BldEdgSh(2)','BldEdgSh(3)','BldEdgSh(4)','BldEdgSh(5)','BldEdgSh(6)']
for scoeffs, coeffs in zip( (scoeffs1,scoeffs2,scoeffs3), (coeffs1,coeffs2,coeffs3)):
    for s,c in zip(scoeffs, coeffs):
        bld[s] = c
# Write the New ElastoDyn file
bld.write('_EDbldFitted.dat')


if __name__ == '__main__':
    plt.show()

if __name__ == '__test__':
    # more tests are already present in tests/test_elastodyn.py
    try:
        os.remove('_EDtwrFitted.dat')
        os.remove('_EDbldFitted.dat')
    except:
        pass 
