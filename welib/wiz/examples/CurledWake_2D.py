""" 
Computes the velocity field in different plane downstream of a yawed wind turbine.

The coordinate system used is:
    z: posivite upwards
    x: streamwise, along the wind

"""
# --- General
import unittest
import numpy as np
import matplotlib.pyplot as plt
# --- Local
from welib.wiz.WindTurbine import WindTurbine
from welib.wiz.Solver import Ct_const_cutoff
from welib.vortilib.elements.VortexPoint import vps_u


def ve_u(Ycp,Zcp,R,theta,theta_y,gamma_x,sigma=None):
    """
    Velocity field induced by a 2D ellipse of vorticity with intensity gamma_x.

    INPUTS:
        theta  : linspace from 0 to 2pi
        gamma_x: vorticity distribution as function of theta, same length as theta
    """
    Ycp=np.asarray(Ycp)
    Zcp=np.asarray(Zcp)

    # Location of the ellipse points
    ye      = R*np.cos(theta)*np.cos(theta_y)
    ze      = R*np.sin(theta)

    gamm=R*gamma_x/(2*np.pi)
    # Smooth Kernel
    if sigma is None:
        K= lambda y,z: 1
    else:
        sigma2=sigma**2
        K= lambda y,z: 1-np.exp(-((z-ze)**2 + (y-ye)**2)/sigma2 )

    uy=np.zeros( np.prod(Ycp.shape))
    uz=np.zeros( np.prod(Ycp.shape))
    # Loop on control points
    for i,(ycp,zcp) in enumerate(zip(Ycp.ravel(),Zcp.ravel())):
        denom   = (ycp-ye)**2 +  (zcp-ze)**2
        Ks = K(ycp,zcp)
        DenInv = gamm/(denom)*Ks
        uy[i]= np.trapz(-(zcp-ze)*DenInv, theta)
        uz[i]= np.trapz( (ycp-ye)*DenInv, theta)
    # Reshape at the end
    uy=uy.reshape(Ycp.shape)
    uz=uz.reshape(Ycp.shape)
    return uy,uz

def vcl_u(Xcp,Ycp,Zcp,R,theta,theta_y,gamma_x,sigma=None):
    """
    Velocity field induced by semi infinite skewed vortex cylinder of longitudinal vorticity

    INPUTS:
        theta  : linspace from 0 to 2pi
        gamma_x: vorticity distribution as function of theta, same length as theta
    """
    Xcp=np.asarray(Xcp)
    Ycp=np.asarray(Ycp)
    Zcp=np.asarray(Zcp)

    # Location of the rotor disk "rim" points, where vorticity starts
    x1 =  - R*np.cos(theta)*np.sin(theta_y)
    y1 =  + R*np.cos(theta)*np.cos(theta_y)
    z1 =  + R*np.sin(theta)

    gamm=R*gamma_x/(4*np.pi)

    # Smooth Kernel
    if sigma is None:
        K= lambda y,z: 1
    else:
        sigma2=sigma**2
        K= lambda y,z: 1-np.exp(-((z-z1)**2 + (y-y1)**2)/sigma2 )

    uy=np.zeros( np.prod(Ycp.shape))
    uz=np.zeros( np.prod(Ycp.shape))
    # Loop on control points
    for i,(x,y,z) in enumerate(zip(Xcp.ravel(),Ycp.ravel(),Zcp.ravel())):
        Den1 = np.sqrt((x-x1)**2 + (y-y1)**2 + (z-z1)**2)   # ||r1||
        Den2 = x-x1                                         #  e.r1
        Ks = K(y,z)
        DenInv = gamm/(Den1*(Den1-Den2))*Ks
        uy[i] = np.trapz((z1-z)*DenInv,theta)
        uz[i] = np.trapz((y-y1)*DenInv,theta)
    # Reshape at the end
    uy=uy.reshape(Ycp.shape)
    uz=uz.reshape(Ycp.shape)
    return uy,uz

def main():
    # --- Parameters for algorithm
    nSpan      = 5    # Number of radial vorticity surfaces
    root       = False
    longi      = False
    tang       = True
    nTheta     = 100  # Number of azimuthal angle used for integration
    bSmooth    = True # Smoothen the vorticity?
    SmoothFact = 10    # Smootheness factor
    # --- Parameters for plotting
    Dprobe = 200      # Distance downstream in diameter where we investigate the plane
    ny     = 30
    nz     = 32
    yLIM   = [-2,2] # xlim in Radius
    zLIM   = [-2,2] # ylim in Radius
    clim   = None
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

    # --- Derived parameters
    theta=np.linspace(0,2*np.pi,nTheta)
    a=R
    b=R*np.cos(yaw)
    # Smooth Param obtaine based on Taylor expansion of distance along the ellipse
    if bSmooth:
        SmoothParam=SmoothFact* np.sqrt(a*np.cos(theta)**2+b*np.sin(theta)**2) * 2*np.pi/nTheta
        #plt.figure()
        #plt.plot(theta,SmoothParam)
        #plt.show()
    else:
        SmoothParam=None
    # --- Control points
    x0    = Dprobe*2*R  # Plane
    y     = np.linspace(yLIM[0]*R,yLIM[1]*R,ny)
    z     = np.linspace(zLIM[0]*R,zLIM[1]*R,nz)
    [Y,Z] = np.meshgrid(y,z)
    X     = 0*Y+x0

    # --- Creating a wind turbine and setting its loading conditions
    WT=WindTurbine(R=R, e_shaft_yaw0=[1,0,0],e_vert=[0,0,1]) # Shaft along x at yaw=0
    WT.update_yaw_pos(yaw)
    WT.update_wind([U0,0,0])
    WT.update_loading(r=vr_bar*R, Ct=Ct_AD, Lambda=Lambda, nCyl=nSpan)
    print('gamma_t     :',WT.gamma_t)
    #plt.figure()
    #plt.plot(WT.r,WT.gamma_t)
    #plt.show()

    print('--- Semi-infinite vortex cylinder')
    ux,uy,uz = WT.compute_u(X,Y,Z,root=root,longi=longi,tang=tang)

    # --- Semi-infinite vortex cylinder longi only
    print('--- Semi-infinite vortex cylinder tuned..')
    uy_vcl = np.zeros(X.shape)
    uz_vcl = np.zeros(X.shape)
    for gamma_t,r in zip(WT.gamma_t,WT.r):
        gamma_l     = gamma_t*np.sin(yaw) * np.sin(theta)
        uy_vcl0,uz_vcl0 = vcl_u(X,Y,Z,r,theta,yaw,gamma_l,SmoothParam)
        uy_vcl+=uy_vcl0
        uz_vcl+=uz_vcl0

    # --- Computing velocity field from derived expression of infinite ellispe
    print('--- From 2D-ellipse integral..')
    uy_an = np.zeros(X.shape)
    uz_an = np.zeros(X.shape)
    for gamma_t,r in zip(WT.gamma_t,WT.r):
        gamma_l  = gamma_t*np.sin(yaw) * np.sin(theta)
        uy_an0,uz_an0 = ve_u(Y,Z,r,theta,yaw,gamma_l,SmoothParam)
        uy_an+=uy_an0
        uz_an+=uz_an0

    # --- Numerical integration uing 2D vortex points of intensity Gamma=gamma_l R dtheta
    print('--- Vortex point integrations..')
    dTheta = 2*np.pi /(nTheta-1)
    CP = np.column_stack((Y.ravel(),Z.ravel()))
    uy_vp = np.zeros(X.shape)
    uz_vp = np.zeros(X.shape)
    if not bSmooth:
        SmoothModel=0
    else:
        SmoothModel=2
    for gamma_t,r in zip(WT.gamma_t,WT.r):
        XV      = np.column_stack((r*np.cos(theta[:-1])*np.cos(yaw), r*np.sin(theta[:-1])))
        gamma_l = gamma_t*np.sin(yaw) * np.sin(theta)
        Gammas  = gamma_l * (r*dTheta)
        ui = vps_u(CP,XV,Gammas,SmoothModel=SmoothModel,KernelOrder=2,SmoothParam=SmoothParam)
        uy_vp+=ui[:,0].reshape(X.shape)
        uz_vp+=ui[:,1].reshape(X.shape)


    # --- Plotting results
    clim=[0,3.0]
    def plot(ax,Y,Z,uy,uz,title='',cb=True):
        Rotor=WT.rotor_disk_points() # Rotor projection
        #Speed=np.sqrt(ux**2)
        Speed=np.sqrt(uy**2+uz**2)
        if clim is not None:
            lev=np.linspace(clim[0],clim[1],30)
        else:
            lev=30
        im=ax.contourf(Y/R,Z/R,Speed,levels=lev)
        ax.plot(Rotor[1,:]/R,Rotor[2,:]/R,'k--') # rotor projected
#         if cb:
#             cb=fig.colorbar(im)
#             if clim is not None:
#                 cb.set_clim(clim)
        sp=ax.streamplot(y/R,z/R,uy,uz,color='k',linewidth=0.7,density=2)
        ax.set_xlim(yLIM[-1::-1]) # NOTE: flipping axis
        ax.set_ylim(zLIM)
        ax.set_xlabel('y/R [-]')
        ax.set_ylabel('z/R [-]')
        ax.set_title(title)

    fig=plt.figure()
    ax1=fig.add_subplot(141)
    ax2=fig.add_subplot(142)
    ax3=fig.add_subplot(143)
    ax4=fig.add_subplot(144)
    plot(ax1,Y,Z,uy    ,uz    ,title = 'Full Cylinder (if needed)',cb = False)
    plot(ax2,Y,Z,uy_vp ,uz_vp ,title = 'Vortex Point (numerical)'  ,cb = False )
    plot(ax3,Y,Z,uy_an ,uz_an ,title = 'Analytical ellipse (paper)',cb = False )
    plot(ax4,Y,Z,uy_vcl,uz_vcl,title = 'Longi cylinder (paper)',cb = True   )

    #dpi=300
    #fig.savefig("CurledWake_yaw{:02d}_CT{:03d}_TWR{:02d}_{:d}D{}.png".format(int(np.round(yaw*180/np.pi)),int(CT0*100),int(Lambda),int(x0/(2*R)),label),dpi=dpi)


if __name__ == "__main__":
    main()
    plt.show()
