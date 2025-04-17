"""
Reproduces figure 3/4 from [1] and Figure [5] from [2]

References:
    [1] E. Branlard, A. Meyer Forsting (2015), Using a cylindrical vortex model to assess the induction zone in front of aligned and yawed rotors
    [2] Branlard, Forsting (2020) Assessing the blockage effect of wind turbines and wind farms
using an analytical vortex model

"""
# --- General
import unittest
import matplotlib.pyplot as plt
import numpy as np
# --- Local
from welib.vortilib.elements.VortexCylinder import vc_tang_u
from welib.tools.colors import *

def get_cmap(minSpeed,maxSpeed):
    DS=0.0001
    # MathematicaDarkRainbow=[(60 /255,86 /255,146/255), (64 /255,87 /255,142/255), (67 /255,107/255,98 /255), (74 /255,121/255,69 /255), (106/255,141/255,61 /255), (159/255,171/255,67 /255), (207/255,195/255,77 /255), (223/255,186/255,83 /255), (206/255,128/255,76 /255), (186/255,61 /255,58 /255)]

    #     ManuDarkOrange  = np.array([198 ,106,1   ])/255.;
    #     ManuLightOrange = np.array([255.,212,96  ])/255.;
    # (1,212/255,96/255),  # Light Orange
    # (159/255,159/255,204/255), # Light Blue
    #     MathematicaLightGreen = np.array([158,204,170 ])/255.;
    # (159/255,159/255,204/255), # Light Blue
    seq=[
    (63/255 ,63/255 ,153/255), # Dark Blue
    (159/255,159/255,204/255), # Light Blue
    (158/255,204/255,170/255), # Light Green
    (1,212/255,96/255),  # Light Orange
    (1,1,1),  # White
    (1,1,1),  # White
    (1,1,1),  # White
    (138/255 ,42/255 ,93/255), # DarkRed
    ]
    valuesOri=np.array([
    minSpeed,  # Dark Blue
    0.90,
    0.95,
    0.98,
    1.00-DS , # White
    1.00    , # White
    1.00+DS , # White
    maxSpeed         # DarkRed
    ])
    values=(valuesOri-min(valuesOri))/(max(valuesOri)-min(valuesOri))

    valuesOri=np.around(valuesOri[np.where(np.diff(valuesOri)>DS)[0]],2)

    cmap= make_colormap(seq,values=values)
    return cmap,valuesOri


def main(test=False):
    if test:
        nZ=25
        nX=25
    else:
        nZ=200
        nX=200

    U0  = 1   ;
    R   = 1   ;
    ZMIN=-5*R; ZMAX=5*R;
    XMIN=-5*R; XMAX=5*R;
    CT0 = 0.95;
    gamma_t=-U0*(1-np.sqrt(1-CT0))*0.85;

#     # --- Along z-Axis
#     fig = plt.figure()
#     ax  = fig.add_subplot(111)
#     for CT0 in [0.95,0.4]:
#         gamma_t = -U0*(1-np.sqrt(1-CT0))*0.85;
#         Zcp=np.linspace(ZMIN,ZMAX,nZ)
#         Xcp=Zcp*0
#         Ycp=Zcp*0
#         ur,uz=vc_tang_u(Xcp,Ycp,Zcp,gamma_t,R)
#         ax.plot(Zcp/R,(uz+U0)/U0,label='CT = {}'.format(CT0))
#     ax.set_xlabel('z/R [-]')
#     ax.set_ylabel('U/U0 [-]')
#     ax.legend()
#     ax.set_title('TwoCTsAxis')
#     ax.set_xlim([-5, 3])
#     ax.set_ylim([0.3, 1])

#     # --- Along r-axis
#     fig = plt.figure()
#     ax  = fig.add_subplot(111)
#     for CT0 in [0.95,0.4]:
#         gamma_t = -U0*(1-np.sqrt(1-CT0))*0.85;
#         Xcp=np.linspace(0,2*R,nZ)
#         Zcp=Xcp*0
#         Ycp=Xcp*0
#         ur,uz=vc_tang_u(Xcp,Ycp,Zcp,gamma_t,R)
#         ax.plot(Xcp/R,(uz+U0)/U0,label='u_z, CT = {}'.format(CT0))
#         ax.plot(Xcp/R,ur        ,label='u_r, CT = {}'.format(CT0))
#     ax.set_xlabel('r/R [-]')
#     ax.set_ylabel('U/U0 [-]')
#     ax.legend()
#     ax.set_title('TwoCTsRotor')
#     ax.set_ylim([0, 1.1])
#     ax.set_xlim([0, XMAX])

    # --- Velocity field, CT=0.95
    CT0     = 0.95                       ;
#     gamma_t = -U0*(1-np.sqrt(1-CT0))*0.85;
    gamma_t = -U0*(1-np.sqrt(1-CT0))
    z=np.linspace(3*R,-5*R,nZ)
    x=np.linspace(-2*R,2*R,nX)
    Z,X=np.meshgrid(z,x)
    Y=Z*0;
    ur,uz = vc_tang_u(X,Y,Z,gamma_t,R)
    print('CT0',CT0,'gamma_t',gamma_t/U0)

    urout,uzout=vc_tang_u([1.1*R],[0],[0.1*R],gamma_t,R)
    print('>>uzout',(uzout+U0)/U0)

    # --- Plot the contours of axial induction
    levels=[0.5,0.6,0.7,0.8,0.9,0.95,0.98,0.99,1.01,1.1]
    try:
        cmap,_ = get_cmap(levels[0],levels[-1])
    except:
        cmap='coolwarm'
    fig=plt.figure()
    ax = fig.add_subplot(111)
    im=ax.pcolormesh (z, x, (uz+U0)/U0, vmin=levels[0], vmax=levels[-1], cmap=cmap)
    fig.colorbar(im)
    cs=plt.contour(Z,X,(uz+U0)/U0,levels=levels,colors='k')
    ax.clabel(cs,levels)
    ax.plot([ZMIN,ZMAX],[0,0],'k:')
    ax.plot([3*R,0,0,3*R],[-R,-R,R,R],'k-',linewidth=3)  # Rotor and Cyl
    ax.set_aspect('equal','box')
    ax.set_xlabel('z/R [-]')
    ax.set_ylabel('r/R [-]')
#     ax.set_xlim([ZMIN,R])
#     ax.set_ylim([XMIN,XMAX])
    ax.set_title('NoYawNoSwirlAxialInductionCT{:03d}'.format(int(CT0*100)))

#     # --- Plot the contours of radial induction
#     if test:
#         levels=[0.005,0.01,0.05,0.1,0.2]
#     else:
#         levels=[0.01,0.05]
#     fig=plt.figure()
#     ax = fig.add_subplot(111)
#     cs=plt.contour(Z,X,ur,levels=levels,colors='k')
#     ax.clabel(cs,levels)
#     ax.plot([ZMIN,ZMAX],[0,0],'k:')
#     ax.plot([3*R,0,0,3*R],[-R,-R,R,R],'k-',linewidth=3)  # Rotor and Cyl
#     ax.set_aspect('equal','box')
#     ax.set_xlabel('z/R [-]')
#     ax.set_ylabel('r/R [-]')
#     ax.set_xlim([ZMIN,1])
#     ax.set_ylim([XMIN,XMAX])
#     ax.set_title('NoYawNoSwirlRadialInductionCT{:03d}'.format(int(CT0*100)))


    # ---
    if not test:
        plt.show()
    else:
        plt.close('all')

class Test(unittest.TestCase):
    def test_Article_Induction_NonYaw(self):
        import sys
        if sys.version_info >= (3, 0):
            main(test=True)
        else:
            print('Test skipped due to travis display error')

if __name__ == "__main__":
    main()


