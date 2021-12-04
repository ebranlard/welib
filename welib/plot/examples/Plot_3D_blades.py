""" 
Examples to plot 3D blade based on AeroDyn inputs

NOTE: might be unfinished
"""
import numpy as np
import matplotlib.pyplot as plt

from welib.plot.surface3d import *
from welib.tools.clean_exceptions import *
from weio.fast_input_deck  import FASTInputDeck

thetaMax=2*np.pi
nTheta=14

z1 = 0
z2 = 1
nz = 3
R1=2
R2=0


# --- Reading AeroDyn file, and extract relevent geometrical info
# NOTE: using FASTInputDeck which attempts to have the same interface as WEIS
# ADFile='../../data/NREL5MW/5MW_Baseline/NRELOffshrBsline5MW_AeroDyn.dat'
ADFile='../../../data/NREL5MW/data/NREL5MW_AD15.05.dat'
#from welib.weio.fast_input_deck  import FASTInputDeck

fst = FASTInputDeck()
fst.readAD(ADFile, readlist=['all'], verbose=True)

coords = fst.fst_vt['ac_data']
# bld = fst.fst_vt['AeroDynBlade']['BldAeroNodes']
bld = fst.fst_vt['AeroDynBlade'].toDataFrame()
# [:,6]
r       = bld['BlSpn_[m]'].values.astype(int)
ID      = bld['BlAFID_[-]'].values.astype(int)
chord   = bld['BlChord_[m]'].values
twist   = bld['BlTwist_[deg]'].values*np.pi/180
prebend = bld['BlCrvAC_[m]'].values
sweep   = bld['BlSwpAC_[m]'].values
# print(bld)

df = coords[4].toDataFrame()
x=df['x/c_[-]'].values
y=df['y/c_[-]'].values


# --- Plotting blade
fig = plt.figure()
ax = fig.add_subplot(111, projection = '3d')
# Loop on radial stations
for i in np.arange(len(r)-1): 
    z1=r[i]
    z2=r[i+1]
    #print('z',z1,z2)
    id = ID[i]-1
    df = coords[id].toDataFrame()
    x=df['x/c_[-]'].values * chord[i]
    y=df['y/c_[-]'].values * chord[i]

    X,Y,Z = arbitrary_cylinder(x, y, z1=z1, z2=z2, nz=2, nTheta=25)
    ax.plot_surface(X, Y, Z, color = 'w', rstride = 1, cstride = 1)

# ax.set_xlim3d(-1, 1)
# ax.set_ylim3d(-1, 1)
# ax.set_zlim3d(-1, 1)
plt.show()
