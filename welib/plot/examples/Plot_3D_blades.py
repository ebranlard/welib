""" 
Examples to plot 3D blade based on AeroDyn inputs

NOTE: might be unfinished
"""
import numpy as np
import os
import matplotlib.pyplot as plt
# Local
from welib.plot.surface3d import *

from welib.tools.clean_exceptions import *
from welib.tools.signal_analysis import resample_interp
from welib.weio.fast_input_deck  import FASTInputDeck
from welib.yams.utils import R_x, R_y, R_z

MyDir=os.path.dirname(__file__)

nTheta=25


# --- Turbine geometrical parameters
nB       = 3
yaw      = 0    * np.pi/180
tilt     = -5   * np.pi/180
cone     = -2.5 * np.pi/180
pitch    = 0   * np.pi/180
TwrHt    = 87.6
Twr2Shft = 1.96256
Overhang = -5.0191
HubRad   = 1.5
# yaw      = 0    * np.pi/180
# tilt     = -0   * np.pi/180
# cone     = -0.0 * np.pi/180
# pitch    = 0   * np.pi/180
# TwrHt    = 0
# Twr2Shft = 0
# Overhang = 0
# HubRad   = 0
ADFile=os.path.join(MyDir, '../../../data/NREL5MW/onshore/NREL5MW_AD.dat')
# ADFile='C:/Work/IEA47/IEA15MW/IEA-15-240-RWT/OpenFAST/Case_V.2.x/AD.dat'


# --- Useful function
def rotTranslate(R_ba, X_a, Y_a, Z_a, r_a0=(0,0,0), r_b0=(0,0,0)):
    X_b = R_ba[0,0]*(X_a-r_a0[0]) + R_ba[0,1]*(Y_a - r_a0[1]) + R_ba[0,2]*(Z_a-r_a0[2]) + r_b0[0] 
    Y_b = R_ba[1,0]*(X_a-r_a0[0]) + R_ba[1,1]*(Y_a - r_a0[1]) + R_ba[1,2]*(Z_a-r_a0[2]) + r_b0[1] 
    Z_b = R_ba[2,0]*(X_a-r_a0[0]) + R_ba[2,1]*(Y_a - r_a0[1]) + R_ba[2,2]*(Z_a-r_a0[2]) + r_b0[2]
    return X_b, Y_b, Z_b


# --- Reading AeroDyn file, and extract relevent geometrical info
# NOTE: using FASTInputDeck which attempts to have the same interface as WEIS
# ADFile='../../data/NREL5MW/5MW_Baseline/NRELOffshrBsline5MW_AeroDyn.dat'
#from welib.weio.fast_input_deck  import FASTInputDeck
fst = FASTInputDeck()
fst.readAD(ADFile, readlist=['all'], verbose=True)

coords = fst.fst_vt['ac_data']
# bld = fst.fst_vt['AeroDynBlade']['BldAeroNodes']
bld = fst.fst_vt['AeroDynBlade'].toDataFrame()


# --- Interpolating to more radial stations
# r_old  = bld['BlSpn_[m]'].values+1.3
# r_new  = np.linspace(r_old[0], r_old[-1], 100)
# _, bld = resample_interp(r_old, r_new, df_old=bld)

# [:,6]

r       = bld['BlSpn_[m]'].values
print(r)
# bld['BlSpn_[m]'] = np.linspace(r[0], r[-1], len(r))
# bld['BlTwist_[deg]']*=0 # TODO TODO HACK
# bld['BlCrvAC_[m]'] = np.linspace(0,-10, len(r))
# bld['BlSwpAC_[m]']*=0 # TODO TODO HACK
# bld['BlChord_[m]']=4 # TODO TODO HACK

dr = np.gradient(bld['BlSpn_[m]'])
dx = np.gradient(bld['BlCrvAC_[m]'])
print('dr',dr)
print('dx',dx)

# bld=bld.iloc[:4,:] # HACK

r       = bld['BlSpn_[m]'].values
ID      = bld['BlAFID_[-]'].values.astype(int)
chord   = bld['BlChord_[m]'].values
twist   = bld['BlTwist_[deg]'].values*np.pi/180
prebend = bld['BlCrvAC_[m]'].values
sweep   = bld['BlSwpAC_[m]'].values
bld['BlCrvAng_[deg]'] = np.degrees(np.arctan2(dx,dr)) # Enforce calculation of slope
curve   = bld['BlCrvAng_[deg]'].values*np.pi/180

# print(bld)


# --- Coordinate transformations
I3 = np.eye(3)
R_ti = np.eye(3)     # from inertial to tower
R_nt = R_z(-yaw)     # from tower to nacelle
R_sn = R_y(tilt)     # from nacelle to non rotating shaft NOTE: sign
R_cB = R_y(-cone)    # from straight blade to coned system
R_bc = R_z(-pitch)   # from coned system to blade coordinate system (pitched)
# Derived
R_ni = R_nt.dot(R_ti) # from inertial to nacelle
R_si = R_sn.dot(R_ni) # from inertial to non rotating shaft

# --- Key points of the structure
r_OT_i = [0       ,0,    0    ] # Point T in inertial coord
r_TN_t = [0       ,0,TwrHt    ] # Tower to nacelle in tower coord
r_NS_n = [0       ,0,Twr2Shft ] # Nacelle to shaft in nacelle coord
r_SH_s = [Overhang,0,0        ] # Shaft to hub in shaft coord 
r_HB_c = [0,0,HubRad          ] # Hub center to blade root in coned coordinates (OpenFAST convention)
# to inertial frame
r_ON_i = r_OT_i + R_ti.T.dot(r_TN_t)
r_OS_i = r_ON_i + R_ni.T.dot(r_NS_n)
r_OH_i = r_OS_i + R_si.T.dot(r_SH_s)
#print('Hub Center: ',r_OH_i)

# --- Coordinates of Aerodynamic center in blade coordinate system
xBAC_b = prebend
yBAC_b = sweep
zBAC_b = r

# --- Setting up figure
fig = plt.figure()
ax3d = fig.add_subplot(111, projection = '3d')
fig.subplots_adjust(left=0.00, right=1.0, top=1.0, bottom=0.0, hspace=0.20, wspace=0.20)
ax3d.w_xaxis.set_pane_color((0.0, 0.0, 0.0, 1.0))
ax3d.w_yaxis.set_pane_color((0.0, 0.0, 0.0, 1.0))
ax3d.w_zaxis.set_pane_color((0.0, 0.0, 0.0, 1.0))
# ax3d.set_xticks([])
# ax3d.set_yticks([])
# ax3d.set_zticks([])
ax3d.set_xticklabels([])
ax3d.set_yticklabels([])
ax3d.set_zticklabels([])
ax3d.view_init(elev=16, azim=-136)

# fig2d,ax2d = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)


# --- Loop on azimuthal positions
for iB, psi in enumerate(np.linspace(0,2*np.pi,nB+1)[:-1]):

    # --- Coordinate systems updates based on azimuth
    R_Bs = R_x(-psi)        # from non rotating shaft to straight blade iB
    R_Bi = R_Bs.dot(R_si)   # from inertial to straight blade for blade iB
    R_ci = R_cB.dot(R_Bi)   # from inertial to coned system  for blade iB
    R_bi = R_bc.dot(R_ci)   # from inertial to blade system for blade iB
    r_OB_i = r_OH_i + R_ci.T.dot(r_HB_c) # Blade iB root 

    #print('Blade root:', r_OB_i)

    # Loop on radial stations
    for i in np.arange(len(r)-1): 
        z1=r[i]
        z2=r[i+1]
        id = ID[i]-1
        df = coords[id].toDataFrame()
#         df = coords[-1].toDataFrame()
        AC_frac = 0.25 # TODO we need the position of the AC (not good for cylinders)

        # --- Coordinate sytem updates based on airfoil section
        R_ab = R_y(-curve[i]).dot(R_z(twist[i]))    # from blade system to airfoil section system # NOTE: sign
        R_ai =  R_ab.dot(R_bi)   # from inertial to airfoil system
        R    = R_ai.T
        if False:
            R=np.eye(3)

        r_BAC_b = [xBAC_b[i], yBAC_b[i], zBAC_b[i] ]

        r_OAC_i = r_OB_i + R_bi.T.dot(r_BAC_b) 
        #print('AC: ', r_OAC_i, r_BAC_b)

        # Coordinates in airfoil coordinate system (x/y swapped)
        x_a = df['y/c_[-]'].values * chord[i]  + 0
        y_a = df['x/c_[-]'].values * chord[i]  - AC_frac*chord[i] 
        z_a = x_a*0                            + 0
        ds = (z2-z1)/np.cos(curve[i]) # curvilinear length
        #dz = (z2-z1)
        X_a,Y_a,Z_a = arbitrary_cylinder(x_a, y_a, z1=0, z2=ds, nz=2, nTheta=nTheta)

        # Coordinate in inertial system
        X_i, Y_i, Z_i = rotTranslate(R_ai.T, X_a, Y_a, Z_a, r_a0=(0,0,0),r_b0=r_OAC_i)
        # Coordinate in blade system
        X_b, Y_b, Z_b = rotTranslate(R_ab.T, X_a, Y_a, Z_a, r_a0=(0,0,0),r_b0=r_BAC_b)

        ax3d.plot_surface(X_i, Y_i, Z_i, color = 'w', rstride = 1, cstride = 1)

#         ax2d.plot(x_a, y_a, 'k-')

ax3d.set_xlabel('x')
ax3d.set_ylabel('y')
ax3d.set_zlabel('z')
ax3d.set_title('Plot - 3D blades')
# ax.set_xlim3d(-1, 1)
# ax.set_ylim3d(-1, 1)
# ax.set_zlim3d(-1, 1)
# ax.set_aspect('auto')
axisEqual3D(ax3d)

# ax2d.set_xlabel('x_g')
# ax2d.set_ylabel('y_g')

if __name__=="__main__":
    plt.show()
if __name__=="__test__":
    pass
if __name__=="__export__":
    from welib.tools.repo import export_figs_callback
    export_figs_callback(__file__)
