"""
Compare the velocity field from a vortex ring to the one of a doublet (far field approximation)

Reference:
    [1] E. Branlard - Wind Turbine Aerodynamics and Vorticity Based Method, Springer, 2017
        Section 35.3 - Page 424
"""
# --- General
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
# --- Local
from welib.vortilib.elements.VortexRing    import ring_u
from welib.vortilib.elements.VortexDoublet import doublet_u
from welib.tools.curves import streamQuiver

def main():
    # --- Parameters
    R     = 1
    Gamma = -10
    XLIM=np.asarray([-5*R,5*R]) # 
    nx = 600    # Number of points for velocity evaluation
    nz = nx 

    # --- Derived parameters
    mz = Gamma * R**2 * np.pi # doublet intensity

    # --- Flow field and speed
    #zs = np.linspace(XLIM[0]*1.1,XLIM[1]*1.1,nx).astype(np.float32)
    #xs = np.linspace(XLIM[0]*1.1,XLIM[1]*1.1,nx).astype(np.float32)
    beg = XLIM[0]*1.1
    end = XLIM[1]*1.1
    zs = np.mgrid[beg:end+((end-beg)/(nx-1)):(end-beg)/(nx-1)]
    xs = np.mgrid[beg:end+((end-beg)/(nx-1)):(end-beg)/(nx-1)]
    [Z,X]=np.meshgrid(zs,xs)
    Y=X*0
    uxr,uyr,uzr = ring_u   (X,Y,Z,Gamma=Gamma,R=R,polar_out=False)
    uxd,uyd,uzd = doublet_u(X,Y,Z,m=[0,0,mz])

    Speed_r=np.log(np.sqrt(uzr**2+uxr**2+uyr**2))
    Speed_d=np.log(np.sqrt(uzd**2+uxd**2+uyd**2))
    supermin=min(np.min(Speed_r),np.min(Speed_d))
    supermax=max(np.max(Speed_r),np.max(Speed_d))/2
    print(supermin)
    print(supermax)

    # --- Plotting
    fig,ax = plt.subplots(1,2, sharex=True, sharey=True)
    levels=np.linspace(supermin,supermax,10)
    im1=ax[0].contourf(Z/R,X/R,Speed_r, levels=levels) # the easiest way to get contourf to comply with cmap is to give levels
    im2=ax[1].contourf(Z/R,X/R,Speed_d, levels=levels) # the easiest way to get contourf to comply with cmap is to give levels
    ax[0].plot([0,0],[-1,1],'k-',lw=2)
    ax[1].plot([0],[0],'ko')

    # Streamlines and quiver
    yseed=np.linspace(XLIM[0]*0.9/R,XLIM[1]*0.9/R,7)
    start=np.array([yseed*0,yseed])
    sp=ax[0].streamplot(zs/R,xs/R,uzr,uxr,color='k',start_points=start.T,linewidth=0.7,density=30,arrowstyle='-')
    qv=streamQuiver(ax[0],sp,n=[5,5,5,5,5,5,5],scale=40,angles='xy')

    yseed=np.linspace(XLIM[0]*0.9/R,XLIM[1]*0.9/R,7);
    yseed=np.concatenate((yseed[0:3],yseed[4:]))
    start=np.array([yseed*0,yseed])
    sp=ax[1].streamplot(zs/R,xs/R,uzd,uxd,color='k',start_points=start.T,linewidth=0.7,density=30,arrowstyle='-')
    qv=streamQuiver(ax[1],sp,n=[5,5,5,5,5,5,5],scale=40,angles='xy')
    ax[0].set_xlabel('z/R [-]')
    ax[1].set_xlabel('z/R [-]')
    ax[0].set_ylabel('r/R [-]')

    # fig.subplots_adjust(right=0.8)
    # cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    # fig.colorbar(im1, cax=cbar_ax)

    ax[0].set_xlim(XLIM/R)
    ax[0].set_ylim(XLIM/R)
    ax[1].set_xlim(XLIM/R)
    ax[1].set_ylim(XLIM/R)


if __name__ == '__main__':
    main()
    plt.show()

