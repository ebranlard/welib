"""
References:
    [1] E. Branlard - Wind Turbine Aerodynamics and Vorticity Based Method, Springer, 2017
    [2] E. Branlard, M. Gaunaa - Cylindrical vortex wake model: skewed cylinder, application to yawed or tilted rotors - Wind Energy, 2015

Coordinate systems
   c coordinate system used in see [2], rotor in plane z_c=0
   w wind coordinate system where z_w is the wind direction
   theta_yaw : yaw angle, positive around y, pointing upward

   x_c =  x_w cost + z_w sint
   y_c =  y_w
   z_c = -x_w sint + z_w cost

"""
# --- General
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os
# --- Local
from welib.vortilib.elements.VortexCylinderSkewed import svcs_tang_u
from welib.vortilib.elements.VortexRing import ring_u
from welib.tools.colors import darkrainbow as cmap
from welib.tools.tictoc import Timer


def main():
    # --- Parameters
    bWindCoord =True
    bRootVortex =True
    R         = 1
    r_hub     = 0.05*R
    CLIM=[0.3,1.1]
    LIM       = [-2,2] # 
    nx        = 50    # Number of points for velocity evaluation
    CT        = 0.6
    theta_yaw =-25*np.pi/180   # rad
    U0        = 1

    # Number of cylinders per radius
    n_radial  = 10 # 1: tip vortex, 2: tip&root, n: "continuous"

    # --- Derived params
    # chi= theta_yaw*(1+0.3*(1-np.sqrt(1-CT))) # rad
    chi= theta_yaw
    if CT==0.4:
        gamma_t = -0.21341 # CT=0.4
    elif CT==0.6:
        gamma_t = -0.40 # 
    else:
        gamma_t = -0.60414 # CT=0.95
    ny = nx
    m  = np.tan(chi) 


    # --- Analytical circulation
    tau_tip  = 0.4 # good values: 0.1-0.8,  0.1:sharp tip vortex
    tau_hub = 0.4 # good values: 0.1-0.8,  0.1:sharp tip vortex
    rh=np.linspace(0,0.5,100)
    rt=np.tan(rh*np.pi)
    gt_hub=1-np.exp(-rt/tau_hub)
    gt_tip=1-np.exp(-rt/tau_tip)
    dgt_tip = 1/tau_tip *(1/(np.cos(rh/tau_tip)**2))*np.exp(-rt/tau_tip)
    dgt_hub = 1/tau_hub *(1/(np.cos(rh/tau_hub)**2))*np.exp(-rt/tau_hub)
    r0  = np.concatenate((rh,0.5+rh))
    r0 = r0*(R-r_hub)+r_hub
    g0  = np.concatenate((gt_hub ,gt_tip[-1::-1]))
    dg0 = np.concatenate((dgt_hub,-dgt_tip[-1::-1]))
    gamma_hub0=np.trapz(dgt_hub,rh)
    gamma_tip0=np.trapz(dgt_hub,rh)
    # print(gamma_hub)
    # plt.plot(rh,dgt_hub,'+')
    # plt.plot(r0,g0)
    # plt.plot(r0,dg0)
    # plt.show()


    nCyl = 1 # Number of wind turbines
    Xcyl=np.zeros(nCyl)
    Ycyl=np.zeros(nCyl)
    Zcyl=np.zeros(nCyl)
    vR       = np.zeros((nCyl,n_radial))
    vM       = np.ones((nCyl,n_radial)) * m
    vgamma_t = np.zeros((nCyl,n_radial))
    vR[0,:]  = np.linspace(R,r_hub,n_radial)
    if n_radial==1:
        vgamma_t[0,:]= gamma_t
    elif n_radial==2:
        vgamma_t[0,0]= gamma_t
        vgamma_t[0,1]=-gamma_t
    else:
        Gamma         = -gamma_t * np.interp(vR[0,:],r0,g0) *0.3
        vgamma_t[0,:] = -gamma_t * np.interp(vR[0,:],r0,dg0)*0.1 # TODO figure out the scalling
        #plt.plot(vR[0,:],Gamma        ,label='Gamma')
        #plt.plot(vR[0,:],vgamma_t[0,:],label='dGamma')
        #plt.legend()
        #plt.show()

    def Tw2c(x_w,y_w,z_w):
        if bWindCoord:
            x_c =  x_w * np.cos(theta_yaw) + z_w * np.sin(theta_yaw)
            y_c =  y_w
            z_c = -x_w * np.sin(theta_yaw) + z_w * np.cos(theta_yaw)
        else:
            x_c,y_c,z_c = x_w,y_w,z_w
        return x_c,y_c,z_c
    def Tc2w(x_c,y_c,z_c):
        if bWindCoord:
            x_w =  x_c * np.cos(theta_yaw) - z_c * np.sin(theta_yaw)
            y_w =  y_c
            z_w =  x_c * np.sin(theta_yaw) + z_c * np.cos(theta_yaw)
        else:
            x_w,y_w,z_w = x_c,y_c,z_c
        return x_w, y_w, z_w


    # --- Loop on diameters
    for nD in [0,4]:
        z0_w      = nD*2*R #Plane
        # --- Flow field and speed
        x_w = np.linspace(LIM[0],LIM[1],nx)
        y_w = np.linspace(LIM[0],LIM[1],ny)
        [X_w,Y_w]=np.meshgrid(x_w,y_w)
        Z_w=X_w*0+z0_w
        X_c,Y_c,Z_c = Tw2c(X_w,Y_w,Z_w) 

        ux_c,uy_c,uz_c  = svcs_tang_u(X_c,Y_c,Z_c,vgamma_t,vR,vM,Xcyl,Ycyl,Zcyl)
        uz_c=uz_c+U0*np.cos(theta_yaw) # Adding free wind
        ux_c=ux_c+U0*np.sin(theta_yaw)
        ux,uy,uz = Tc2w(ux_c,uy_c,uz_c)

        def plot(ux,uy,uz,label='',clim=None):
            Speed=np.sqrt(uz**2)
            # Temporary HACK until singularity is removed
            print('Min Max: ',np.min(Speed.ravel()),np.max(Speed.ravel()))
            if clim is not None:
                Speed[Speed>clim[1]] = clim[1]
                Speed[Speed<clim[0]] = clim[0]
            print('Min Max: ',np.min(Speed.ravel()),np.max(Speed.ravel()))

            # rotor projection
            vpsi=np.linspace(0,2*np.pi,50)
            xc_w=R*np.cos(vpsi)*np.cos(theta_yaw)
            yc_w=R*np.sin(vpsi)

            dpi=300
            fig=plt.figure()
            ax=fig.add_subplot(111)
            if clim is not None:
                lev=np.linspace(clim[0],clim[1],30)
            else:
                lev=30
            im=ax.contourf(X_w,Y_w,Speed,levels=lev,cmap=cmap)
            ax.plot(xc_w,yc_w,'k--')
            cb=fig.colorbar(im)
            #if clim is not None:
            #    cb.set_clim(CLIM)
            sp=ax.streamplot(x_w,y_w,ux,uy,color='k',linewidth=0.7,density=2)

            ax.set_xlim(LIM)
            ax.set_ylim(LIM)
            ax.set_xlabel('x/R [-]')
            ax.set_ylabel('y/R [-]')
            ax.set_title('z = {}D{}'.format(int(z0_w/(2*R)),label))
            fig.savefig("VC_superp_yaw{:02d}_CT{:03d}_{:d}D{}.png".format(int(np.round(theta_yaw*180/np.pi)),int(CT*100),int(z0_w/(2*R)),label),dpi=dpi)

        plot(ux  ,uy  ,uz  ,' cylinders',clim=CLIM)


if __name__ == '__main__':
    main()
    plt.show()

