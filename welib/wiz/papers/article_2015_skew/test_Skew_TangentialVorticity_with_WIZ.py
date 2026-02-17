"""
Tangential Vorticity from a skewed vortex cylinder.

This test reproduces SOME OF the plots of the section 3.1 of the following reference using the "wiz" library:
    [1] E. Branlard, M. Gaunaa - Cylindrical vortex wake model: skewed cylinder, application to yawed or tilted rotors - Wind Energy, 2015

LOOK AT: welib.vortilib.elements.papers.article_2015_skew test_Skew_TangentialVorticity.py for more plots

""" 


# --- General
import unittest
import matplotlib.pyplot as plt
import numpy as np
from numpy import pi
# --- Local
from welib.vortilib.elements.VortexCylinderSkewed import fV_Tangential, fKxit, fKzt
from welib.essentials import *


from welib.wiz.WindTurbine import WindTurbine, transform

def setupWT(gamma_t, chi, R, U0=10):
    """ Setup WindTurbine for one tangential cylinder of intensity gamma_t and a skew angle"""
    WT = WindTurbine(R=R,r_hub=[0,0,0],e_shaft_yaw0=[1,0,0],e_vert=[0,0,1])
    WT.update_loading(gamma_t_Ct=lambda x: np.asarray([gamma_t]), Ct=1, nCyl=1)
    WT.update_yaw_pos(-chi) # NOTE: chi is minus yaw pos of nacelle
    WT.update_wind([U0,0,0]) # Useless if only_ind is True
    WT.T_c2g=np.dot(WT.T_wt2g, WT.T_c2wt) # NOTE, I need to add it to my class..
    return WT

def coordinates_pc2g(T_c2g, VR_pc, VPSI_pc, VZ_pc):
    """ Get from "polar cylinder (pc) coordinates" (the ones from the paper) to global (g) 
    TODO: consider origin..
    """
    # First cylinder polar to cylinder cart See figure 1 from paper
    Xc=VR_pc * np.cos(VPSI_pc)
    Yc=VR_pc * np.sin(VPSI_pc)
    Zc=VZ_pc
    # Then to global (OK, origin is 0,0,0)
    Xg, Yg, Zg = transform(T_c2g, Xc, Yc, Zc)
    return Xg, Yg, Zg

def compute_u_pc2c(WT, VR_pc, VPSI_pc, VZ_pc):
    """ Compute velocity based on inputs in polar cylinder (pc) coordinates 
        Output velocity in cylinder coordinates (c) (Cartesian Cylinder)
    """
    # Get coordinates from cylinder polar to global
    Xg, Yg, Zg = coordinates_pc2g(WT.T_c2g, VR_pc, VPSI_pc, VZ_pc)
    # Compute velocity in global
    ux,uy,uz = WT.compute_u(Xg, Yg, Zg, root=False, longi=False, tang=True, only_ind=True)
    # Back to cylinder
    uxc, uyc, uzc = transform(WT.T_c2g.T, ux, uy, uz)
    return uxc, uyc, uzc


# --------------------------------------------------------------------------------}
# --- Main 
# --------------------------------------------------------------------------------{
def main(test=False):
    R           = 1
    gamma_rings = - 1
    ntheta      = 600
    epsilon     = 10 ** - 3
    chi         = 30 * pi / 180
    m           = np.tan(chi)

    IPlot=[]
    IPlot+=[2] # z component azimuthal survey
    IPlot+=[3] # z component azimuthal survey

    WT = setupWT(gamma_rings, chi, R=R)
    print(WT)

    # --------------------------------------------------------------------------------
    # --- Z COMPONENT
    # --------------------------------------------------------------------------------
    if 2 in IPlot:
        # --- Z component Azimuthal survey
        npsi = 3
        nr   = 200
        vpsi = np.linspace(0,pi / 2,npsi)
        vr   = np.linspace(- 2 * R,2 * R,nr) / R
        vr   = vr[np.abs((np.abs(vr) - 1)) > epsilon]
        vz0  = 0
        VR,VPSI,VZ0 = np.meshgrid(vr,vpsi,vz0)

        # Using Vortex Cylinder
        uiz = fV_Tangential(VR,VPSI,VZ0,m,gamma_rings,ntheta,nout=1)

        # Using Wind turbine class
        uxc, uyc, uzc = compute_u_pc2c(WT, VR, VPSI, VZ0)

        # Plot velocity in cylinder to match paper
        plt.figure()
        for ip in np.arange(len(vpsi)):
            plt.plot(vr,np.squeeze(uiz[ip,:,:]), 'k-')
            plt.plot(vr, uzc[ip,:,:],':', label='$\\psi=%3d$ deg.'%np.round(vpsi[ip] * 180 / pi))
        plt.xlabel(r'Radial position $r/R$')
        plt.ylabel(r'Axial induced velocity $-u_{z,t}/\gamma_t$ ')
        plt.title('CastlesAzimuth_chi{:03d}'.format(int(np.round(np.arctan(m)*180/pi))))
        plt.legend()


    if 3 in IPlot:
        # --- Z component Radial survey
        npsi = 360
        vpsi = np.linspace(- pi,pi,npsi)
        vr = np.array([0.1,0.5,0.9])
        vz0 = 0
        VR,VPSI,VZ0 = np.meshgrid(vr,vpsi,vz0)

        #  Using Vortex Cylinder
        uiz = fV_Tangential(VR,VPSI,VZ0,m,gamma_rings,ntheta,nout=1)

        # Using Wind turbine class
        uxc, uyc, uzc = compute_u_pc2c(WT, VR, VPSI, VZ0)

        # Plot velocity in cylinder to match paper
        plt.figure()
        for ip in np.arange(len(vr)):
            plt.plot(vpsi*180/pi,np.squeeze(uiz[:,ip,:]), 'k-')
            plt.plot(vpsi*180/pi,np.squeeze(uzc[:,ip,:]), ':', label='$r/R=%2.1f$'%vr[ip])
        plt.legend()
        plt.xlabel(r'Azimuthal position $\psi$ [deg]')
        plt.ylabel(r'Axial induced velocity  $-u_{z,t}/\gamma_t$')
        plt.xlim(np.array([- 180,180]))
        plt.title('CastlesRadialchi{:03d}'.format(int(np.round(np.arctan(m)*180/pi))))


class Test(unittest.TestCase):
    def test_Article_Skew_TangentialVorticity(self):
        main(test=True)
        plt.close('all')

if __name__ == "__main__":
    main()
    plt.show()
