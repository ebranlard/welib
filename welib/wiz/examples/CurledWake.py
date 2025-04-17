""" 
Computes the velocity field in different plane downstream of a yawed wind turbine.

The coordinate system used is:
    z: posivite upward
    x: streamwie, along the wind

"""
# --- General
import unittest
import numpy as np
# --- Local

from welib.wiz.WindTurbine import WindTurbine
from welib.wiz.Solver import Ct_const_cutoff

# --------------------------------------------------------------------------------}
import matplotlib.pyplot as plt

def main():
    # --- Parameters for algorithm
    nCyl  = 1
    root  = False
    longi = False
    tang  = True 
    # --- Parameters for plotting
    vD = [0,10] # Wake diameters 
    ny = 40
    nz = 42
    yLIM =[-2,2] # xlim in Radius
    zLIM =[-2,2] # ylim in Radius
    clim =None
    label=''
    # --- Parameters for turbine and turbine loading
    yaw       = -30*np.pi/180  # [rad]
    CT0       = 0.80
    Lambda    = 6
    U0        = 10
    R         = 65
    r_bar_cut = 0.11
    r_bar_tip = 0.90
    vr_bar    = np.linspace(0,1.0,100)
    Ct_AD     = Ct_const_cutoff(CT0,r_bar_cut,vr_bar,r_bar_tip) # TODO change me

    # --- Creating a wind turbine and settting its loading conditions
    WT=WindTurbine(R=R, e_shaft_yaw0=[1,0,0],e_vert=[0,0,1]) # Shaft along x at yaw=0
    WT.update_yaw_pos(yaw)
    WT.update_wind([U0,0,0])
    WT.update_loading(r=vr_bar*R, Ct=Ct_AD, Lambda=Lambda, nCyl=nCyl)
    #ux,uy,uz = WT.compute_u(np.array([0]),np.array([0]),np.array([0]),root=root,longi=longi,tang=tang)

#     plt.figure,
#     plt.plot(WT.r,WT.gamma_t,label='gamma_t')
#     plt.plot(WT.r,WT.gamma_l,label='gamma_l')
#     plt.plot(WT.r,WT.Ct,label='Ct')
#     plt.legend()
#     plt.show()

    Rotor=WT.rotor_disk_points() # Rotor projection
    # --- Loop on downstream diameters
    for iD in vD:
        x0    = iD*2*R  # Plane
        y     = np.linspace(yLIM[0]*R,yLIM[1]*R,ny)
        z     = np.linspace(zLIM[0]*R,zLIM[1]*R,nz)
        [Y,Z] = np.meshgrid(y,z)
        X     = 0*Y+x0
        ux,uy,uz = WT.compute_u(X,Y,Z,root=root,longi=longi,tang=tang)

        Speed=np.sqrt(ux**2)
        #Speed=np.sqrt(ux**2+uy**2+uz**2)
        # Temporary HACK until singularity is removed
        #print('Min Max: ',np.min(Speed.ravel()),np.max(Speed.ravel()))
        #if clim is not None:
        #    Speed[Speed>clim[1]] = clim[1]
        #    Speed[Speed<clim[0]] = clim[0]
        #print('Min Max: ',np.min(Speed.ravel()),np.max(Speed.ravel()))


        dpi=300
        fig=plt.figure()
        ax=fig.add_subplot(111)
        if clim is not None:
            lev=np.linspace(clim[0],clim[1],30)
        else:
            lev=30
        im=ax.contourf(Y/R,Z/R,Speed,levels=lev)
        ax.plot(Rotor[1,:]/R,Rotor[2,:]/R,'k--') # rotor projected
        cb=fig.colorbar(im)
        if clim is not None:
            cb.set_clim(clim)
        sp=ax.streamplot(y/R,z/R,uy,uz,color='k',linewidth=0.7,density=2)

        ax.set_xlim(yLIM[-1::-1]) # NOTE: flipping axis
        ax.set_ylim(zLIM)
        ax.set_xlabel('y/R [-]')
        ax.set_ylabel('z/R [-]')
        ax.set_title('x = {}D{}'.format(int(x0/(2*R)),label))
        #fig.savefig("CurledWake_yaw{:02d}_CT{:03d}_TWR{:02d}_{:d}D{}.png".format(int(np.round(yaw*180/np.pi)),int(CT0*100),int(Lambda),int(x0/(2*R)),label),dpi=dpi)

if __name__ == "__main__":
    main()
    plt.show()
