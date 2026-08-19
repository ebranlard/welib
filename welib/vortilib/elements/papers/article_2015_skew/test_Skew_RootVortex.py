"""
Trailed (longitudinal) Vorticity from a skewed vortex cylinder.
This test reproduces the plots of the section 3.3 of the following reference:
    [1] E. Branlard, M. Gaunaa - Cylindrical vortex wake model: skewed cylinder, application to yawed or tilted rotors - Wind Energy, 2015


Coordinate system:
 rotor plane x-y  (y is vertical, x longi)
 wake is z (if not yawed)
""" 

## Parameters
# --- General
import unittest
import matplotlib.pyplot as plt
import numpy as np
from numpy import pi
# --- Local
from welib.vortilib.elements.VortexCylinderSkewed import fV_Root
from welib.tools.colors import fColrs

def pol2cart(rho,phi):
    x= rho * np.cos(phi)
    y= rho * np.sin(phi)
    return (x,y)

def main(test=False):
    #np.warnings.filterwarnings('error')
    ## Params
    Gamma_tot = 1
    Gamma_r = - 1
    visc_model = 0
    t = 0
    R = 1
    bComputeGrad = 0

    ## --- n lines along a circle - AZIMUTHAL SURVEY - WITH ROOT VORTEX
    chi    = 30 * pi / 180
    e      = np.array([np.sin(chi),0,np.cos(chi)])
    n      = 100
    vtheta = np.linspace(0,2 * pi,n)
    npsi   = 100
    vpsi   = np.linspace(0,2 * pi,npsi)
    vr     = np.array([0.2,0.5,1])
    nr     = len(vr)
    vz0    = 0
    m      = np.tan(chi)
    VR,VPSI,VZ0 = np.meshgrid(vr,vpsi,vz0)
    v_z_root ,v_psi_root  = fV_Root(VR,VPSI,VZ0,m,Gamma_r,nout = 2)
    v_z_root0,v_psi_root0 = fV_Root(VR,VPSI,VZ0,0,Gamma_r,nout = 2)

    fig, ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8))
    for ir,r in enumerate(vr):
        ax.plot(vpsi*180/pi,np.squeeze(v_z_root[:,ir,:])*4*pi*R/Gamma_tot, label='$r/R=%2.1f$'%r)
    ax.plot(vpsi*180/pi, np.squeeze(v_z_root0[:,0,:])*4*pi*R/Gamma_tot,'k--',linewidth=2,label=r'Value for $\chi=0$')
    ax.legend()
    ax.set_ylabel(r'$u_{z,r} 4 \pi R / \Gamma_{tot}$')
    ax.set_xlabel(r'Azimuthal position $\psi$ [deg]')
    ax.set_title(r'RootVortexOnly Axial velocity for various radii Chi%2d'%np.round(chi * 180 / pi))
    ax.set_xlim(np.array([0,360]))

    fig, ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8))
    for ir,r in enumerate(vr):
        ax.plot(vpsi*180/pi,np.squeeze(v_psi_root[:,ir,:])*4*pi*R/Gamma_tot, label='$r/R=%2.1f$'%r)
    ax.plot(vpsi*180/pi, np.squeeze(v_psi_root0[:,0,:])*4*pi*R/Gamma_tot,'k--',linewidth=2,label=r'Value for $\chi=0$')
    ax.legend()
    ax.set_ylabel(r'$u_{\psi,r} 4 \pi R / \Gamma_{tot}$')
    ax.set_xlabel(r'Azimuthal position \psi [deg]')
    ax.set_title(r'RootVortexOnlyTangential velocity for various radii Chi%2d'%np.round(chi*180/pi))
    ax.set_ylim(np.array([-8.8,0]))
    ax.set_xlim(np.array([0,360]))

    ## --- Contour plot of axial velocity
    n_azimuth = 180
    n_radial  = 25
    vpsi      = np.linspace(0,360,n_azimuth) * pi / 180
    vr        = np.linspace(0.0,1.5,n_radial) * R
    vz0       = 0
    chi       = 30 * pi / 180
    m         = np.tan(chi)
    VR,VPSI,VZ0 = np.meshgrid(vr,vpsi,vz0)
    VR   = VR.squeeze()
    VPSI = VPSI.squeeze()
    VZ0  = VZ0.squeeze()
    Xplane,Yplane = pol2cart(VR,VPSI)
    # Ui
    uz = fV_Root(VR,VPSI,VZ0,m,Gamma_r,nout=1)

    # Plot
    fig, ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8))
    ContourVal = np.array([0,0.01,0.05,0.1,0.2])
    ContourVal = np.unique(np.sort(np.concatenate((-ContourVal,ContourVal))/ Gamma_tot))
    CS = ax.contour(Xplane/R, Yplane/R, uz/Gamma_tot,levels=ContourVal)
    manual_locations=[ (0.1667575,0.418641), (0.1492567,- 0.4228), (0.0758036,- 0.2096), (0.0213203,0.203710), (- 0.12992,0.702902), (- 0.13028,- 0.7035), (- 0.66038,0.175041), (- 0.66048,- 0.1747)]
    lbs=plt.clabel(CS,inline=1,fontsize=10,manual=manual_locations)
    rots=[ - 350.3, - 6.172, - 6.064, - 12.97, - 36.47, - 323.4, - 336.0, - 23.83]
    for l,r in zip(lbs,rots):
        l.set_rotation(r)
    ax.set_aspect('equal')
    vtheta = np.linspace(0,2 * pi,100)
    ax.plot(R * np.cos(vtheta),R * np.sin(vtheta),'k')
    #ax.set_aspect('square')
    ax.set_xlim(np.array([  1.3,-1.3])) # NOTE: flipped
    ax.set_ylim(np.array([- 1.3,1.3]))
    ax.set_title('RootVortexOnlyContourAxial_chi%03d'%np.round(np.arctan(m) * 180 / pi))
    ax.set_xlabel('$x/R$')
    ax.set_ylabel('$y/R$')


class Test(unittest.TestCase):
    def test_Article_Skew_RootVortex(self):
        import sys
        if sys.version_info >= (3, 0):
            main(test=True)
            plt.close('all')
        else:
            print('Test skipped due to travis display error')

if __name__ == "__main__":
    main()
    plt.show()
